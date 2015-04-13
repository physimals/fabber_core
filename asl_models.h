/*   asl_model.h Kinetic curve models for ASL

      Michael Chappell - IBME & FMRIB Image Analysis Group

      Copyright (C) 2010-2011 University of Oxford */

/*   CCOPYRIGHT   */

#if !defined(asl_models_h)
#define asl_models_h

#include "miscmaths/miscmaths.h"
#include "miscmaths/miscprob.h"

using namespace MISCMATHS;

#include "utils/tracer_plus.h"

using namespace Utilities;

namespace OXASL {

  // generic AIF model class
  class AIFModel {
  public:
    //evaluate the model
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const = 0;
    //report the number of dispersion parameters 
    virtual int NumDisp() const = 0;
    // return default priors for the parameters
    virtual ColumnVector Priors() const { return priors; }
    virtual string Name() const = 0;
    virtual void SetPriorMean(int paramn, double value) { priors(paramn) = value; }

  protected:
    ColumnVector priors; //list of prior means and precisions - all means first then precisions

  };

  //Specific AIF models
  class AIFModel_nodisp : public AIFModel {
    //AIFModel_nodisp() {}
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl, const ColumnVector dispparam) const;
  virtual int NumDisp() const {return 0;}
  virtual string Name() const { return "None"; }
  };

  class AIFModel_gammadisp : public AIFModel {
  public:
    AIFModel_gammadisp() { priors.ReSize(4); priors << 2 << -0.3 << 1 << 1; } //old prec.s were 10 and 10
    // corresponds to mean s~0.7 and p~0.1
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 2;}
    virtual string Name() const { return "Gamma dispersion kernel"; }
  };

  class AIFModel_gvf : public AIFModel {
  public:
    AIFModel_gvf() { priors.ReSize(4); priors << 2 << 1.6 << 1 << 1; }    
    // corresponds to mean s~0.7 and p~0.7
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 2;}
    virtual string Name() const { return "GVF"; }
  };

  class AIFModel_gaussdisp : public AIFModel {
  public:
    AIFModel_gaussdisp() { priors.ReSize(2); priors << -1.6 << 1; }    
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 1;}
    virtual string Name() const { return "Gauss dispersion kernel"; }
  };

  class AIFModel_spatialgaussdisp : public AIFModel {
  public:
    AIFModel_spatialgaussdisp() { priors.ReSize(2); priors << -1.4 << 1; }    
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 1;}
    virtual string Name() const { return "Spatial gauss dispersion kernel"; }
  };

  class AIFModel_spatialgaussdisp_alternate : public AIFModel {
  public:
    AIFModel_spatialgaussdisp_alternate() { priors.ReSize(2); priors << -1.4 << 1; }    
    virtual double kcblood(const double ti, const double deltblood, const double taub, const double T_1b, bool casl,const ColumnVector dispparam) const;
    virtual int NumDisp() const {return 1;}
    virtual string Name() const { return "Spatial gauss dispersion kernel (alternate)"; }
  };

  //  double kcblood_gallichan(const double ti, const double deltblood, const double taub, const double T_1b, const double xdivVm, bool casl);

  // -------------
  //generic residue function model class
  class ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const = 0;
    //report the number of residue function parameters
    virtual int NumResid() const = 0;
    // return the default priors for the parameters
    virtual ColumnVector Priors() const {return residpriors;}
    virtual string Name() const = 0;
    virtual void SetPriorMean(int paramn, double value) { residpriors(paramn) = value; }

  protected:
    ColumnVector residpriors;
  };

  //specific residue function models
  class ResidModel_wellmix : public ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "Well mixed"; }
  };

  class ResidModel_simple : public ResidModel {
  public:
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;
    
    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "Simple"; }
  };

  class ResidModel_imperm : public ResidModel {
  public:
    ResidModel_imperm() { residpriors.ReSize(2); residpriors << 0.5 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "Impermeable"; }
  };

  class ResidModel_twocpt : public ResidModel {
  public:
    ResidModel_twocpt() { residpriors.ReSize(2); residpriors << 0.8 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "Two comparment (no backflow, no venous output)"; }
  };

  class ResidModel_spa : public ResidModel {
  public:
    ResidModel_spa() { residpriors.ReSize(6); residpriors << 0.02 << 0.03 << 2 << 1e-3 << 1e12 << 10; }
    virtual double resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const;

    virtual int NumResid() const {return 3;}
    virtual string Name() const { return "Single Pass Approximation (2 compartment, no backflow)"; }
  };

  // ------------
  //generic tissue model class
  class TissueModel {
  public:
    //evalute the model
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const = 0;
    // report the number of dipersion parameters
    virtual int NumDisp() const = 0;
    // report the number of residue function parameters (beyond the normal ones)
    virtual int NumResid() const = 0;
    // return default priors for the parameters
    virtual ColumnVector DispPriors() const { return disppriors; }
    virtual ColumnVector ResidPriors() const {return residpriors; }
    virtual string Name() const = 0;
    virtual void SetDispPriorMean(int paramn, double value) { disppriors(paramn) = value; }
    virtual void SetResidPriorMean(int paramn, double value) { residpriors(paramn) = value; }

  protected:
    ColumnVector disppriors; //list of prior means and precisions - all means first then precisions
    ColumnVector residpriors;
  };

  //specific tissue models

  class TissueModel_nodisp_simple : public TissueModel {
    virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return 0;}
    virtual int NumResid() const {return 0;}
    virtual string Name() const { return "No dispersion | Simple"; }
  };

  class TissueModel_nodisp_wellmix : public TissueModel {
  public:
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 0;}
 virtual string Name() const { return "No dispersion | Well mixed"; }
  };

  class TissueModel_nodisp_imperm : public TissueModel {
  public:
    TissueModel_nodisp_imperm() { residpriors.ReSize(2); residpriors << 0.5 << 10; }
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 1;}
 virtual string Name() const { return "No dispersion | Impermeable"; }
  };

  class TissueModel_nodisp_2cpt : public TissueModel {
  public:
    TissueModel_nodisp_2cpt() { residpriors.ReSize(2); residpriors << 0.8 << 10; }
    virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return 0;}
    virtual int NumResid() const {return 1;}
    virtual string Name() const { return "No dispersion | Two compartmentr (no backflow, no venous output)"; }
  };

  class TissueModel_nodisp_spa : public TissueModel {
  public:
    TissueModel_nodisp_spa() { residpriors.ReSize(6); residpriors << 0.02 << 0.03 << 2 << 1e-3 << 1e12 << 10; }
virtual double kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 0;}
 virtual int NumResid() const {return 3;}
 virtual string Name() const { return "No dispersion | Single Pass Approximation (2 compartment no backflow)"; }

  private:
 double Q(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const;
 double R(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const;
  };

  class TissueModel_gammadisp_wellmix : public TissueModel {
  public:
    TissueModel_gammadisp_wellmix() { disppriors.ReSize(4); disppriors << 2 << -0.3 << 10 << 10; } 
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
 virtual int NumDisp() const {return 2;}
 virtual int NumResid() const {return 0;}
  virtual string Name() const { return "Gamma kernel dispersion | Well mixed"; }
  };


  //double kctissue_gvf(const double ti, const double delttiss,const double tau, const double T_1b, const double T_1app, const double s, const double p);
  //double kctissue_gaussdisp(const double ti, const double delttiss, const double tau, const double T_1b, const double T_1app, const double sig1, const double sig2);

  //  a general tissue model that does numerical convolution
  class TissueModel_aif_residue : public TissueModel {
  public:
    TissueModel_aif_residue(AIFModel* paifmodel, ResidModel* presidmodel)
      { aifmodel=paifmodel; residmodel = presidmodel;
	disppriors << aifmodel->Priors();
	residpriors << residmodel->Priors();
      }
    virtual double kctissue(const double ti, const double fcalib, const double delttiss, const double tau, const double T_1b, const double T_1app, const double lambda, const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const;
    virtual int NumDisp() const {return aifmodel->NumDisp();}
    virtual int NumResid() const {return residmodel->NumResid();}
    virtual string Name() const { string name; name = "NUMERICAL CONVOLUTION - " + aifmodel->Name() + residmodel->Name(); return name; }

  protected:
    AIFModel* aifmodel;
    ResidModel* residmodel;
  };

  // useful functions
  double icgf(const double a, const double x);
  double gvf(const double t, const double s, const double p);
  double numerical_integration(ColumnVector integrand, double del, double finalval, double finaldel, string method);


}

#endif
