/*  fwdmodel.cc - base class for generic forward models

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel.h"

#include <sstream> 
#include "easylog.h"
 
string FwdModel::ModelVersion() const
{
  return "No version info available.";
  // You should overload this function in your FwdModel class
  
  // Something like this:
  // return " $ I d $ "; // without the spaces
  // CVS will automatically replace this with version information that looks
  // like this: $Id: fwdmodel.cc,v 1.38 2014/10/24 15:25:27 chappell Exp $
  // It should probably go in your .cc file, not the header.
}

int FwdModel::Gradient(const ColumnVector& params, Matrix& grad) const
{
  // by default return false -> no gradient is supplied by this model
  return false;
}

int FwdModel::NumOutputs() const
{
    ColumnVector params, result;
    params.ReSize(NumParams()); params = 1; // probably safer than 0
    Evaluate(params, result);
    return result.Nrows();
}

void FwdModel::DumpParameters(const ColumnVector& params, 
                              const string& indent) const
{
    LOG << indent << "Parameters:" << endl;
    vector<string> names;
    NameParams(names);
    
    for (int i = 1; i <= NumParams(); i++)
        LOG << indent << "  " << names[i-1] << " == " << params(i) << endl;
        
    LOG << indent << "Total of " << NumParams() << " parameters." << endl;
}

void FwdModel::pass_in_coords( const ColumnVector& coords )
{
  coord_x = coords(1);
  coord_y = coords(2);
  coord_z = coords(3);
}

#include "fwdmodel_simple.h"
#include "fwdmodel_quipss2.h"
#include "fwdmodel_q2tips.h"
#include "fwdmodel_pcASL.h"
#include "fwdmodel_flobs.h"
#include "fwdmodel_custom.h"

#include "fwdmodel_cest.h"
#ifdef __OXASL
// models associated with Oxford_asl components
#include "fwdmodel_asl_grase.h"
//#include "fwdmodel_asl_buxton.h"
#include "fwdmodel_asl_pvc.h"
#include "fwdmodel_asl_satrecov.h"
#include "fwdmodel_asl_quasar.h"
#include "fwdmodel_asl_rest.h"
#include "fwdmodel_asl_multiphase.h"
#endif /*__OXASL */

#ifdef __DEVEL
// development models (not to be complied in FSL)
#include "fwdmodel_dsc.h"
//#include "fwdmodel_asl_devel.h"
//#include "fwdmodel_asl_2cpt.h"
//#include "fwdmodel_asl_multite.h"
//#include "fwdmodel_cest_devel.h"
//  #include "fwdmodel_asl_dynangio.h" - this is currently broken since changing asl_models, use asl_models_*_archive
//#include "fwdmodel_fasl.h"
//#include "fwdmodel_oef.h"
#include "fwdmodel_biexp.h"
//#include "fwdmodel_deconv.h"
//#include "fwdmodel_flex.h"
//#include "fwdmodel_asl_samira.h"
//#include "fwdmodel_asl_precap.h"

#include "fwdmodel_dce.h"

#endif /* __DEVEL */

// Add your models here

FwdModel* FwdModel::NewFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
    if (name == "simple")
    {
        LOG_ERR("WARNING: the 'simple' forward model has several hard-coded filenames and probably won't work outside of FMRIB\n");
        return new SimpleFwdModel(args);
    }
    else if (name == "quipss2")
    {
        return new Quipss2FwdModel(args);
    } 
    else if (name == "q2tips-dualecho") // NOT the interleaved-echo Q2TIPS described in Woolrich06
    {
        return new Q2tipsFwdModel(args);
    } 
    else if (name == "pcasl-dualecho")
      {
	return new pcASLFwdModel(args);
      }
    else if (name == "linear")
      {
	return new LinearFwdModel(args);
      }
    else if (name == "flobs7")
      {
	// flobs7 = polar coordinate parameterization
	// b1 cos(b2) X1 + b1 sin(b2) X2
	// prior on b2 should probably be modified slightly (unless it's very small).
        return new FlobsFwdModel(args, true);
      }
    else if (name == "flobs5")
      {
	// flobs5 = b1X1 + b1b2X2 parameterization
        return new FlobsFwdModel(args, false);
      }

    else if (name == "cest")
      {
	return new CESTFwdModel(args);
      }
	
#ifdef __OXASL
	//models associated with Oxford_asl components
	else if (name == "buxton")
	{
		return new GraseFwdModel(args);
	}
    //   else if (name == "buxton_simple")
    //	{
    //		return new BuxtonFwdModel(args);
    //	}
    else if (name == "aslpvc")
	{
		return new ASL_PVC_FwdModel(args);
	}
    else if (name == "satrecov")
	{
		return new SatrecovFwdModel(args);
	}
	else if (name == "quasar")
      {
	return new QuasarFwdModel(args);
      }
    else if (name == "aslrest")
      {
	return new ASLFwdModel(args);
      }
    else if (name == "asl_multiphase")
      {
	return new MultiPhaseASLFwdModel(args);
      }

#endif /* __OXASL */
	
#ifdef __DEVEL
// development models (not to be complied in FSL)
    else if (name == "dsc")
      {
	return new DSCFwdModel(args);
      }
    //  else if (name == "devel")
    //   {
    //	return new DevelFwdModel(args);
    //   }

    // else if (name == "cestdevel")
    //  {
    //	return new CESTDevelFwdModel(args);
    //  }
    //else if (name == "dynangio")
    //  {
    //	return new DynAngioFwdModel(args);
    //  }
    // else if (name == "fasl")
    //   {
    //	return new FASLFwdModel(args);
    //   }
    //else if (name == "oef")
    //  {
    //	return new OEFFwdModel(args);
    //  }
    else if (name == "biexp")
      {
	return new BiExpFwdModel(args);
      }
    //else if (name == "deconv")
    //  {
    //	return new DeconvModel(args);
    //  }
    //else if (name == "flex")
    //  {
    //	return new FLEXFwdModel(args);
    //  }
    //else if (name == "samira")
    //  {
    //	return new SamiraFwdModel(args);
    //  }
    //else if (name == "aslprecap")
    //  {
    //	return new PreCapFwdModel(args);
    //  }
    //    else if (name == "twocpt")
    //  {
    //	return new TwoCptFwdModel(args);
    //  }
    //else if (name == "multite")
    //  {
    //	return new multiTEFwdModel(args);
    //  }
    else if (name == "dce")
      {
       return new DCEFwdModel(args);
      }



#endif /* __DEVEL */
    
    else if (name == "custom")
      {
	return new CustomFwdModel(args);
      }
    // Your models go here!
    else
    {
        throw Invalid_option("Unrecognized forward model '" + name + "'\n");
    }  
}

// If you want usage information from --help --model=yourmodel, add it below.
void FwdModel::ModelUsageFromName(const string& name, ArgsType& args)
{
    // Update this to add your own models to the code
    
//    if (name == "simple")
//    {
//        return new SimpleFwdModel(args);
//    }
//    else 
    if (name == "quipss2")
    {
        Quipss2FwdModel::ModelUsage();
    } 
    if (name == "q2tips")
    {
        cout << "Note: --model=q2tips uses exactly the same options as --model=quipss2:\n";
	Quipss2FwdModel::ModelUsage();
    } 
    if (name == "pcasl-dualecho")
      {
	pcASLFwdModel::ModelUsage();
      }
    else if (name == "linear")
      {
	LinearFwdModel::ModelUsage();
      }
    else if (name == "custom")
      {
	CustomFwdModel::ModelUsage();
      }

#ifdef __OXASL
	else if (name == "buxton")
	{
		GraseFwdModel::ModelUsage();
	}
    //   else if (name == "buxton-simple")
    //	{
    //		BuxtonFwdModel::ModelUsage();
    //	}
    else if (name == "asl_pvc")
	{
		ASL_PVC_FwdModel::ModelUsage();
	}
    else if (name == "aslrest")
      {
	ASLFwdModel::ModelUsage();
      }
    else if (name == "asl_multiphase")
      {
	MultiPhaseASLFwdModel::ModelUsage();
      }
#endif /* __OXASL */
	
#ifdef __DEVEL

#endif /* __DEVEL */

    // Your models go here!
    else
    {
        cout << "No model usage information available for --model=" 
            << name << endl;
    }   
}

