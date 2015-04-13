/*  fwdmodel_fasl.cc -  multi-TI functional ASL model

    Michael Chappell, IBME & FMRIB Image Analysis Group

    Copyright (C) 2010 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_fasl.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string FASLFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_fasl.cc,v 1.2 2012/01/13 12:00:59 adriang Exp $";
}

void FASLFwdModel::HardcodedInitialDists(MVNDist& prior, 
    MVNDist& posterior) const
{
    Tracer_Plus tr("FASLFwdModel::HardcodedInitialDists");
    assert(prior.means.Nrows() == NumParams());
    
    // Set priors
    prior.means = 0;
    SymmetricMatrix precisions = IdentityMatrix(NumParams()) * 1e-12;

    // if we have TC difference data then dont infer static mag
    if (dataisdiff) {
      precisions(M0index(),M0index()) = 1e99;
      for (int i=1; i<=Mbasis.Ncols(); i++) {
	precisions(M0index()+i,M0index()+i) = 1e99;
      }
    }

    // We are more informative about arrival time
    prior.means(D0index()) = 0.7;
    precisions(D0index(),D0index()) = 10;
    for (int i=1; i<=Dbasis.Ncols(); i++) {
      prior.means(D0index()+i) = 0;
      precisions(D0index()+i,D0index()+i) = 1e6;//100; //only small changes in deltaT are likely
      }
    
    if (stdevT1b>0)
      {
	    prior.means(T1bIndex()) = fixedT1b;
	    precisions(T1bIndex(),T1bIndex()) = 1/(stdevT1b*stdevT1b);
      }
    if (stdevT1>0)
      {
	    prior.means(T1Index()) = fixedT1;
	    precisions(T1Index(),T1Index()) = 1/(stdevT1*stdevT1);
      }
    if (stdevInvEff>0)
      {
	    prior.means(InvEffIndex()) = fixedInvEff;
	    precisions(InvEffIndex(),InvEffIndex()) = 1/(stdevInvEff*stdevInvEff);
      }

    prior.means(AIndex()) = 1;
    precisions(AIndex(),AIndex()) = 10;

    prior.SetPrecisions(precisions);
    
    // Set informative initial posterior
    posterior = prior;

    posterior.means(M0index()) = 1.5e4; //roughly what we expect the magntiude to be
    posterior.means(Q0index()) = 1e-3*posterior.means(M0index()); //initial value based on static mag
    precisions(Q0index(),Q0index()) = 1; 
    // NB arrival time is sufficicently informative to not need specific intialization    

    //posterior.means(R0index()) = 25;

    posterior.SetPrecisions(precisions);

}    

void FASLFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    Tracer_Plus tr("FASLFwdModel::Evaluate");
    
    // Parameterization used in most recent results:
    // Absolute M and Q change (same units as M0 or Q0):
    ColumnVector StatMag = params(M0index()) - Mbasis * MnOf(params);
    ColumnVector CBF = params(Q0index()) + Qbasis * QnOf(params);
    double D0cpy = params(D0index()); if (D0cpy<0) D0cpy=0; //D0 cannot be negative
    ColumnVector DelT = D0cpy + Dbasis * 100*DnOf(params); //Note scaling of arrival time changes here!

    /* THIS SECTION FOR MULTI ECHO
    // Fractional change in BOLD effect (at TE_2), rather than using % R2* change 
    ColumnVector R2s = -1/echoTime(2) * log( 
        Rbasis * RnOf(params) + exp(-echoTime(2)*params(R0index())));
    */

    double T1b = (stdevT1b>0 ? params(T1bIndex()) : fixedT1b);
    double T1 = (stdevT1>0 ? params(T1Index()) : fixedT1);
    double invEfficiency = (stdevInvEff>0 ? params(InvEffIndex()) : fixedInvEff);

    //kinetic curve modelling
    int npoints = CBF.Nrows();
    ColumnVector Sb(npoints);
    Sb=0.0;
    float lambda=0.9;
    float T_1app = 1/( 1/T1 + 0.01/lambda);
    for (int i=1; i<=npoints; i++) {
      if (dataisdiff) {
	Sb(i) = CBF(i)*2*invEfficiency*kctissue_nodisp(tivec(i),DelT(i),tau,T1b,T_1app);
      }
      else {
	if (rho(i) < 0) {
	  //Sb(i) = 0.0;
	  Sb(i) = -CBF(i)*2*invEfficiency*kctissue_nodisp(tivec(i),DelT(i),tau,T1b,T_1app);
	}
      }
    }

    if (dataisdiff) {
      StatMag=0;
    }
    else {
    //saturation recovery of static magnetization
    double A = params(AIndex());
    for (int i=1; i<=npoints; i++) {
      StatMag(i) = StatMag(i)*(1.0 - A*exp(-tivec(i)/T1));
	}
    }

    ColumnVector S = StatMag + Sb;
    result = S;
      
    //int Ntimes = S.Nrows();
  //int Nte = echoTime.Nrows();
  //if (result.Nrows() != Nte*Ntimes)
  //  result.ReSize(Nte*Ntimes);
  
  /* ASSEMBLE multi TE result  
for (int te = 1; te <= Nte; te++)
    {
      ColumnVector nuisance = Nbasis * NnOf(te, params);
      // Will be all-zero if there are no nuisance regressors
        
      for (int i = 1; i <= Ntimes; i++)
        result( Nte*(i-1) + te ) = 
            S(i) * exp(-echoTime(te) * R2s(i)) + nuisance(i);
      // Fill order: te1 te2 te1 te2 te1 te2 te1 te2 ...
    }
  */

  return;
}

void FASLFwdModel::ModelUsage()
{
    cout << "\nUsage info for --model=quipss2:\n"
      << "Required options:\n"
      << "--bold-basis=<bold_design_file>\n"
      << "--cbf-basis=<cbf_design_file>\n"
      << "--statmag-basis=<statmag_design_file>\n\n"
      << "Optional options:\n"
      << "--nuisance-basis=<nuisance_regressors_design_file> (default: null)\n"
      << "--ti1=<ti1_in_sec>, "
      << "--ti2=<ti2_in_sec> (default: 0.6, 1.5)\n"
      //      << "--te1=<te1_in_millisec>, "
      //<< "--te2=<te2_in_millisec> (default: 9.1, 30)\n"
//      << "--te3=<next_echo_time>, etc.\n"
      << "--tag-pattern=<string_of_Ts_and_Cs> (default: TC)\n"      
      << "--t1b=<T1_of_blood> (default: 1.66), --t1b-stdev=<stdev> (to add it as a parameter)\n"
      //<< "--dt=<bolus_arrival_time>, --dt-stdev (default: --dt=0.5 --dt-stdev=0.25)\n"
      << "--inv-eff=<inversion_efficiency>, --inv-eff-stdev=<stdev> (to add it as a parameter)\n\n"      
      ;
}

FASLFwdModel::FASLFwdModel(ArgsType& args)
{
  Tracer_Plus tr("FASLFwdModel::FASLFwdModel");

    string scanParams = args.ReadWithDefault("scan-params","cmdline");
    string tagPattern;
    ColumnVector tis(1); //will add extra values onto end as needed
    int ndiscarded;
    
    if (scanParams == "cmdline")
    {
	    stdevT1b = convertTo<double>(args.ReadWithDefault("t1b-stdev", "0"));
        fixedT1b = convertTo<double>(args.ReadWithDefault("t1b","1.66"));
	    stdevT1 = convertTo<double>(args.ReadWithDefault("t1-stdev", "0"));
        fixedT1 = convertTo<double>(args.ReadWithDefault("t1","1.3"));
	    stdevInvEff = convertTo<double>(args.ReadWithDefault("inv-eff-stdev","0"));
        fixedInvEff = convertTo<double>(args.ReadWithDefault("inv-eff","1"));
        stdevDt = convertTo<double>(args.ReadWithDefault("dt-stdev","0.25"));
        fixedDt = convertTo<double>(args.ReadWithDefault("dt","0.5"));

	tau = convertTo<double>(args.Read("tau"));
        
        if (stdevInvEff < 0 || stdevDt < 0 || stdevT1b < 0)
            throw Invalid_option("standard deviations must not be negative!");
    
        tagPattern = args.ReadWithDefault("tag-pattern","TC");
        if (tagPattern.find_first_not_of("TCtc") != tagPattern.npos)
            throw Invalid_option("tagpattern string must contain only Ts and Cs!");
        
        /* MULTI ECHO IS CURRENTLY DISABLED - COULD PROBABLY BE ADDED BACK IN LATER
echoTime.ReSize(2);
        echoTime(1) = convertTo<double>(args.ReadWithDefault("te1","9.1"))/1000.0;
        echoTime(2) = convertTo<double>(args.ReadWithDefault("te2","30"))/1000.0;

	while (true) 
	  {
	    int N = echoTime.Nrows()+1;
	    string teString = args.ReadWithDefault("te"+stringify(N), "stop!");
	    if (teString == "stop!") break;

            // This just isn't tested enough (at all)... remove if you dare
            throw Invalid_option(
              "Using more than two echo times is implemented but completely untested... modify the code if you really want to try it.");

	    // Append this TE to the list of TEs
	    ColumnVector tmp(1); 
	    tmp = atof(teString.c_str())/1000.0;
	    echoTime &= tmp; // vertcat

	    // Sanity checks:
	    if (echoTime(N) <= 0.001)
	      throw Invalid_option(
		"Was expecting TE > 1 ms (don't use seconds!)");
	    if (echoTime(N) > 0.500)
	      throw Invalid_option("Was expecting TE < 500 ms");
	  }
	*/

	//TIs - these are TIs and not PLDs (TI = PLD + bolus duration)
      tis(1) = atof(args.Read("ti1","0").c_str());
      
      while (true) //get the rest of the tis
	{
	  int N = tis.Nrows()+1;
	  string tiString = args.ReadWithDefault("ti"+stringify(N), "stop!");
	  if (tiString == "stop!") break; //we have run out of tis
	 
	  // append the new ti onto the end of the list
	  ColumnVector tmp(1);
	  tmp = convertTo<double>(tiString);
	  tis &= tmp; //vertical concatenation

	}
      //timax = tis.Maximum(); //dtermine the final TI

      ndiscarded = convertTo<int>(args.ReadWithDefault("discarded","0")); //number of datapoints discarded from start of data
      dataisdiff = args.ReadBool("dataisdiff");

      //LOGGING HERE!
        LOG << " --t1b=" << fixedT1b << "--t1b-stdev=" << stdevT1b
            << " --inv-eff=" << fixedInvEff << " --inv-eff-stdev=" << stdevInvEff
            << " --dt=" << fixedDt << "--dt-stdev=" << stdevDt
            << " --tag-pattern=" << tagPattern
            ;
	/*	for (int i = 1; i <= echoTime.Nrows(); i++)
	  LOG << " --te" << i << "=" << echoTime(i)*1000.0;
	  LOG << endl;*/
    }

    else
        throw Invalid_option("Only --scan-params=cmdline is accepted at the moment");    
    
    //string rb = args.Read("bold-basis"); // bold basis set
    string qb = args.Read("cbf-basis");  // CBF basis set
    string db = args.ReadWithDefault("delt-basis",qb); //arr time basis set (default to CBF basis set)
    string mb = args.Read("statmag-basis");
    //string nb = args.ReadWithDefault("nuisance-basis", "null");
    
    // LOG_ERR( "    Reading BOLD basis functions: " << rb << endl );
    //if (rb != "null")
    //    Rbasis = read_vest_fabber(rb);
    //else
    //    throw Invalid_option("Currently --bold-basis=null isn't allowed...");
    
    LOG_ERR( "    Reading CBF basis functions: " << qb << endl );
        Qbasis = read_vest_fabber(qb);

    LOG_ERR( "    Reading arrival time basis functions: " << qb << endl );
        Dbasis = read_vest_fabber(db);
        
    LOG_ERR( "    Reading Stat. Mag. basis functions: " << mb << endl );
        Mbasis = read_vest_fabber(mb);

	const int numTR = Qbasis.Nrows();
	/*
    LOG_ERR( "    Reading Nuisance basis functions: " << nb << endl );
    if (nb != "null")
        Nbasis = read_vest_fabber(nb);
    else
        Nbasis.ReSize(numTR, 0);
	*/
    
    // Now we can parse the TagPattern string.
        
    rho.ReSize(numTR);
    
    for (unsigned i = 1; i <= tagPattern.length(); i++)
        rho(i) = (toupper(tagPattern[i-1]) == 'T') ? -1 : 1;
        
    for (int i = tagPattern.length()+1; i<=numTR; i++)
        rho(i) = rho(i-tagPattern.length());

    LOG << "Full tag-control pattern used (" << rho.Nrows() << " TRs): ";
    for (int i = 1; i <= rho.Nrows(); i++)
	LOG << (rho(i)>0 ? "C" : "T");
    LOG << endl;

    // expand the TIs for the full duration

    ColumnVector tivectemp(numTR+ndiscarded);
    tivectemp.Rows(1,tis.Nrows()) = tis; // load in the first block of TIs
    for (int i=tis.Nrows()+1; i<=numTR+ndiscarded; i++) {
      tivectemp(i) = tivectemp(i-tis.Nrows());
    }

    tivec = tivectemp.Rows(ndiscarded+1,numTR+ndiscarded);

    LOG << "Full TIs used: " << tivec.t() << endl;
}

void FASLFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{

}

void FASLFwdModel::NameParams(vector<string>& names) const
{
    names.clear();
    
    names.push_back("Q0");
    for (int i = 1; i <= Qbasis.Ncols(); i++)
      names.push_back("Q_abschg_" + stringify(i));
    
    names.push_back("M0");
    for (int i = 1; i <= Mbasis.Ncols(); i++)
      names.push_back("M_abschg_" + stringify(i));

    names.push_back("D0");
    for (int i = 1; i <= Dbasis.Ncols(); i++)
      names.push_back("D_abschg_" + stringify(i));
    
    //   names.push_back("R0");
    //  for (int i = 1; i <= Rbasis.Ncols(); i++)
    //   names.push_back("BOLD_abschg_" + stringify(i));

    //    for (int i = 1; i <= Nbasis.Ncols(); i++)
    //   names.push_back("Nuisance_signal_" + stringify(i));

    if (stdevInvEff>0)
      names.push_back("InvEff");

    if (stdevT1b>0)
      names.push_back("T1b");

    if (stdevT1>0)
      names.push_back("T1");

    names.push_back("A");

    //   if (stdevDt>0)
    //  names.push_back("dt");
  
    assert(names.size() == unsigned(NumParams())); 
}


float FASLFwdModel::kctissue_nodisp(const float ti, const float delttiss, const float tau, const float T_1b, const float T_1app) const {
  Tracer_Plus tr("FASLFwdModel::kctissue_nodisp");
float kctissue;
 kctissue=0.0;

  // Tissue kinetic curve no dispersion (cASL)
  // Buxton (1998) model

  //float R = 1/T_1app - 1/T_1b; 
  float F = 2 * exp(-ti/T_1app);

      if(ti < delttiss)
	{ kctissue = 0;}
      else if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue = F * T_1app * exp(-delttiss/T_1b) * (1 - exp(-(ti-delttiss)/T_1app));
	}
      else //(ti > delttiss + tau)
	{
	  kctissue = F * T_1app * exp(-delttiss/T_1b) * exp(-(ti-tau-delttiss)/T_1app) * (1 - exp(-tau/T_1app));
	}

  return kctissue;
}
