/* asl_models.cc Kinetic curve models for ASL

Michael Chappell - IBME & FMRIB Analysis Group

Copyright (C) 2010-2011 University of Oxford */

/* CCOPYRIGHT  */

#include "asl_models.h"

namespace OXASL {

// --- Kinetic curve functions ---
//Arterial

  double AIFModel_nodisp::kcblood(const double ti, const double deltblood, const double taub, const double T_1b, const bool casl,const ColumnVector dispparam) const {
  // Non dispersed arterial curve
  Tracer_Plus tr("OXASL:kcblood_nodisp");
  double kcblood = 0.0;


      if(ti < deltblood)
	{ 
	  kcblood = 2 * exp(-deltblood/T_1b) * (0.98 * exp( (ti-deltblood)/0.05 ) + 0.02 * ti/deltblood );
	  // use a arti§fical lead in period for arterial bolus to improve model fitting
	}
      else if(ti >= deltblood && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else      kcblood = 2 * exp(-ti/T_1b); 
	}
      else //(ti > deltblood + tau)
	{
	  // artifical lead out period for taub model fitting
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else      kcblood = 2 * exp(-(deltblood+taub)/T_1b);

	  kcblood *= (0.98 * exp( -(ti - deltblood - taub)/0.05) + 0.02 * (1-(ti - deltblood - taub)/5));
	 
	  if (kcblood<0) kcblood=0; //negative values are possible with the lead out period equation
	}


  return kcblood;
}

  //NOTE: for cASL the version here is an over simplificaiton (just changing the decay term and leaving the rest alone) since it ignores the fact that some blood will be more delayed than the rest due to dispersion
  double AIFModel_gammadisp::kcblood(const double ti,const double deltblood,const double taub,const double T_1b,const bool casl,const ColumnVector dispparam) const {
    // Gamma dispersed arterial curve (pASL)
    Tracer_Plus tr("OXASL:kcblood_gammadisp");
    double kcblood = 0.0;

    //extract dispersion parameters
    double s; double p;
    s = (dispparam.Row(1)).AsScalar();
    s = exp(s);
    double sp = (dispparam.Row(2)).AsScalar();
    sp = exp(sp);
    if (sp>10) sp=10;
    p = sp/s;
        
    double k=1+p*s;

    if(ti < deltblood)
      { 
	kcblood = 0.0;
	}
    else if(ti >= deltblood && ti <= (deltblood + taub))
      { 
	if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	else	  kcblood = 2 * exp(-ti/T_1b);
	
	kcblood *= ( 1 - igamc(k,s*(ti-deltblood))  ); 
      }
      else //(ti > deltblood + taub)
	{
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else	    kcblood = 2 * exp(-ti/T_1b);
	  kcblood *= ( igamc(k,s*(ti-deltblood-taub)) - igamc(k,s*(ti-deltblood)) ) ; 
	  
	}
    
    return kcblood;
  }

  double AIFModel_gvf::kcblood(const double ti,const double deltblood,const double taub,const double T_1b,const bool casl,const ColumnVector dispparam) const {
    //GVF AIF shape
    Tracer_Plus tr("OXASL:kcblood_gvf");
    double kcblood = 0.0;

    //extract dispersion parameters
    double s; double p;
    s = (dispparam.Row(1)).AsScalar();
    s = exp(s);
    double sp = (dispparam.Row(2)).AsScalar();
    sp = exp(sp);
    if (sp>10) sp=10;
    p = sp/s;

    // gamma variate arterial curve
    // NOTE:    taub does not directly affect the shape, jsut scale of this KC - see below
    //          NOT a good idea to use when inferring bolus duration.

      if(ti < deltblood)
	{ 
	  kcblood = 0.0;
	}
      else //if(ti >= deltblood) && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else	    kcblood = 2 * exp(-ti/T_1b);

	  kcblood *= gvf(ti-deltblood,s,p); 
	}
      // we do not have bolus duration within the GVF AIF - the duration is 'built' into the function shape

      kcblood /= taub; //the 'original' bolus duration scales the magtiude of the KC becuase the area under the KC is preserved under dispersion (apart from T1 decay)
  return kcblood;
}
 
  //NOTE: for cASL the version here is an over simplificaiton (just changing the decay term and leaving the rest alone) since it ignores the fact that some blood will be more delayed than the rest due to dispersion
  double AIFModel_gaussdisp::kcblood(const double ti,const double deltblood,const double taub,const double T_1b,const bool casl,const ColumnVector dispparam) const {
  // Gaussian dispersion arterial curve
  // after Hrabe & Lewis, MRM, 2004 for pASL
    // my derrivation based on the spatial gauss solution for cASL
  Tracer_Plus tr("OXASL:kcblood_gaussdisp");
  double kcblood = 0.0;
  double sqrt2 = sqrt(2);

  double sig1;
  sig1 = exp( (dispparam.Row(1)).AsScalar() );
  double erf1;
  double erf2;
  
  if (casl) {
    double a = 1/(sig1*sig1);
    double b = 1/T_1b;

    double Q = (2*a*(deltblood-ti)-b)/(2*a);
    double S = (a*(deltblood-ti)*(deltblood-ti) + b*ti)/a;

    kcblood = 2 * exp(-deltblood/T_1b)*exp(Q*Q-S);

    erf1 = sqrt(a)*(taub+Q);
    erf2 = sqrt(a)*Q;
  }
  else {
    // we will only have diserpsion SD (leading edge)
    // assume trailing edge is related to leading edge as Hrabe did
    double sig2;
   if (deltblood > 0) {
     sig2  = sig1*sqrt( (deltblood+taub)/deltblood );
   }
   else sig2 = sig1;
   
   kcblood = 2 * exp(-ti/T_1b);
   
   erf1 =  (ti-deltblood)/(sqrt2*sig1);
   erf2 = (ti-deltblood-taub)/(sqrt2*sig2);
 }
 
 if (erf1>5) erf1=5;
 if (erf2>5) erf2=5;
 if (erf1<-5) erf1=-5;
 if (erf2<-5) erf2 = -5;
 
 kcblood *= 0.5*( MISCMATHS::erf(erf1)  - MISCMATHS::erf(erf2)  );


  return kcblood;
}

  double AIFModel_spatialgaussdisp_alternate::kcblood(const double ti,const double deltblood,const double taub,const double T_1b,const bool casl,const ColumnVector dispparam) const {
  // Gaussian dispersion arterial curve - in spatial rather than temporal domain
  // after Ozyurt ISMRM 2010 (p4065) for pASL
  // using derrivation from Thijs van Osch for cASL, but developed into closed form solution
    // This is the (orignal, now) alternate version that assumes that dispersion happens with TI
    // NB pASl case is same for both versions
  
  Tracer_Plus tr("OXASL:kcblood_spatialgaussdisp");
  double kcblood = 0.0;

  double k;
  k = exp( (dispparam.Row(1)).AsScalar() );
  double erf1;
  double erf2;

  if (casl) {
    double a = 1/(k*k*ti);
    double b = 1/T_1b;
    
    double Q = (2*a*(deltblood-ti)-b)/(2*a);
    double S = (a*(deltblood-ti)*(deltblood-ti) + b*ti)/a;
    
    kcblood = 2 * exp(-deltblood/T_1b)*exp(Q*Q-S);
    
    erf1 = sqrt(a)*(taub+Q);
    erf2 = sqrt(a)*Q;
  }
  else {
    kcblood = 2 * exp(-ti/T_1b);
    
    erf1 =  (ti-deltblood)/(k*sqrt(ti));
    erf2 = (ti-deltblood-taub)/(k*sqrt(ti));
  }
    
    if (erf1>5) erf1=5;
    if (erf2>5) erf2=5;
    if (erf1<-5) erf1=-5;
    if (erf2<-5) erf2 = -5;
    
    kcblood *= 0.5*( MISCMATHS::erf(erf1)  - MISCMATHS::erf(erf2)  );
    
  return kcblood;
}

double AIFModel_spatialgaussdisp::kcblood(const double ti,const double deltblood,const double taub,const double T_1b,const bool casl,const ColumnVector dispparam) const {
  // Gaussian dispersion arterial curve - in spatial rather than temporal domain
  // after Ozyurt ISMRM 2010 (p4065) for pASL
  // using derrivation from Thijs van Osch for cASL, but developed into closed form solution
    // This is the (orignal, now) alternate version that assumes that dispersion happens with TI
  
  Tracer_Plus tr("OXASL:kcblood_spatialgaussdisp");
  double kcblood = 0.0;

  double k;
  k = exp( (dispparam.Row(1)).AsScalar() );
  double erf1;
  double erf2;

  if (casl) {    
    double dt = 0.01;
    int ndels = floor( ( ti-max(0.0,ti-taub) )/dt ) + 1;
    ColumnVector integrand(ndels);
    double lambda;
      for (int i=1; i<=ndels; i++) {
	lambda = ti-(i-1)*dt;
	if (lambda < 1e-12) {
	  integrand = 0;
	}
	  else {
	    integrand(i) = 1/sqrt(lambda) * exp(-1/(k*k) * ( (deltblood - lambda)*(deltblood - lambda)/lambda ) )*exp(-lambda/T_1b);
	  }
      }
      // final bit
      lambda = max(1e-12,ti-taub);
      double finalintegrand;
      finalintegrand = 1/sqrt(lambda) * exp(-1/(k*k) * ( (deltblood - lambda)*(deltblood - lambda)/lambda ) )*exp(-lambda/T_1b);
      double finaldel = (ti-(ndels-1)*dt) - lambda;
      
      double integral;
      integral = numerical_integration(integrand,dt,finalintegrand,finaldel,"trapezium");
    
      kcblood = integral*dt* 1/sqrt(M_PI) * 1/k;
    
  }
  else {
    kcblood = 2 * exp(-ti/T_1b);
    
    erf1 =  (ti-deltblood)/(k*sqrt(ti));
    erf2 = (ti-deltblood-taub)/(k*sqrt(ti));
    
    if (erf1>5) erf1=5;
    if (erf2>5) erf2=5;
    if (erf1<-5) erf1=-5;
    if (erf2<-5) erf2 = -5;
    
    kcblood *= 0.5*( MISCMATHS::erf(erf1)  - MISCMATHS::erf(erf2)  );


  }

    
  return kcblood;
}

/*
double kcblood_gallichan(const double ti,const double deltblood,const double taub,const double T_1b,const double xdivVm,const bool casl=false) {
  // Model of dispersion based on a geometrical argument from Gallichan MRM 2008
  // Taking equation [6] (so not including QUIPSSII style saturation)
  // including an 'extra' arrival time term as per the paper
  // bolus duration (taub) takes the place of X/V_m and we let it be a variable
    //Tracer_Plus tr("OXASL:kcblood_normdisp");
  double kcblood = 0.0;

  assert(casl==false); // this model is pASL only

      // NOTE: the +xdivVm correction applied to the ti to shift the curve so that the 
      // delay associated with the dispersion parameter has been removed, thus BAT is independent
      // of the dispersion.

      if(ti < deltblood)
	{ 
	  kcblood = 0.0;
	}
      else if((ti >= deltblood) && ti <= (deltblood + taub))
	{ 
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else	    kcblood = 2 * exp(-ti/T_1b);

	  kcblood *= 1 - xdivVm/(ti + xdivVm - deltblood); 
	}
      else //(ti > deltblood + taub)
      	{
	  if (casl) kcblood = 2 * exp(-deltblood/T_1b);
	  else	    kcblood = 2 * exp(-ti/T_1b);

      	  kcblood *= taub/(ti + xdivVm - deltblood) ; 
	}
  return kcblood;
}
  */
  //-----------------------------------------
  //Residue functions (these are primiarly specified for use with the numerical tissue model)
  double ResidModel_wellmix::resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const {
    // Well mixed single compartment
    // Buxton (1998) model
    Tracer_Plus tr("OXASL::residmodel_wellmix");

    double T_1app = 1/( 1/T_1 + fcalib/lambda );
    return exp(-ti/T_1app);
  }

double ResidModel_simple::resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const {
    // Simple impermeable comparment
  // decays with T1b
    Tracer_Plus tr("OXASL::residmodel_simple");

    return exp(-ti/T_1b);;
  }

  double ResidModel_imperm::resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const {
    // impermeable compartment with transit time
    //decays with T1b
    Tracer_Plus tr("OXASL::residmodel_imperm");

    double transit = (residparam.Row(1)).AsScalar();
    double resid = exp(-ti/T_1b);
    if (ti>transit) resid=0.0;
    return resid;
  }

  double ResidModel_twocpt::resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const {
    // Two compartment model - the simplest form of the two cpt model
    // No backflow from tissue to blood
    // no venous outflow
    // From Parkes & Tofts and also St. Lawrence 2000 - both models are the same under these assumptions
    Tracer_Plus tr("OXASL::residmodel_twocpt");

    // extract residue function parameters
    double kw; //exchange rate = PS/vb
    kw = (residparam.Row(1)).AsScalar();
    //double PS; double vb;
    //PS = (residparam.Row(1)).AsScalar();
    //vb = (residparam.Row(2)).AsScalar();

    // calculate the residue function
    double a = kw + 1/T_1b;
    double b = (kw*T_1*T_1b) / (kw*T_1*T_1b + (T_1-T_1b) );
    return b*exp(-ti/T_1) + (1-b)*exp(-a*ti);
    }

  double ResidModel_spa::resid(const double ti, const double fcalib, const double T_1, const double T_1b, const double lambda, const ColumnVector residparam) const {
    // Two compartment model - Single Pass Approximation from St. Lawrence (2000)
    // No backflow from tissue to blood
    // label starts to leave the cappilliary after a capilliary transit time
    Tracer_Plus tr("OXASL::residmodel_spa");

// extract residue function parameters
    double PS; double vb; double tauc;
    PS = (residparam.Row(1)).AsScalar();
    vb = (residparam.Row(2)).AsScalar();
    tauc = (residparam.Row(3)).AsScalar();

    // calcualte residue function
    double a = PS/vb + 1/T_1b;
    double b = (PS*T_1*T_1b) / (PS*T_1*T_1b + (T_1-T_1b)*vb );
    double ER = 1 - exp(-PS/fcalib - (1/T_1b - 1/T_1)*tauc);

    if (ti < tauc) {
      return b*exp(-ti/T_1) + (1-b)*exp(-a*ti);
    }
    else {
      return b*ER*exp(-ti/T_1);
    }

  }

  //----------------------------------
  //Tissue Model
double TissueModel_nodisp_simple::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
  // Tissue kinetic curve - well mixed, but no outflow and decay with T1 blood only
  // (This is just the impermeable model with infinite residence time)
  Tracer_Plus tr("OXASL::kctissue_nodisp_simple");
  double kctissue = 0.0;


  if(ti < delttiss)
    { kctissue = 0;}
  else if(ti >= delttiss && ti <= (delttiss + tau))
    {
      if (casl) kctissue = exp(-delttiss/T_1b) - exp(-ti/T_1b);
      else      kctissue = ti - delttiss;
    }
  else //(ti > delttiss + tau)
    {
      if (casl) kctissue = exp(-ti/T_1b)*(exp(tau/T_1b)-1);
      else      kctissue = tau;
    }

  if (casl)  kctissue *= T_1b;
  else       kctissue *= exp(-ti/T_1b);
  kctissue *= 2;

  return kctissue;
}

  double TissueModel_nodisp_wellmix::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
  // Tissue kinetic curve no dispersion
  // Buxton (1998) model
  Tracer_Plus tr("OXASL::kctissue_nodisp_wellmix");
  double kctissue = 0.0;

  double T_1app = 1/( 1/T_1 + fcalib/lambda );
  double R = 1/T_1app - 1/T_1b; 
  double F = 2 * exp(-ti/T_1app);

      if(ti < delttiss)
	{ kctissue = 0;}
      else if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  if (casl)  kctissue = 2 * T_1app * exp(-delttiss/T_1b) * (1 - exp(-(ti-delttiss)/T_1app));
	  else       kctissue = F/R * ( (exp(R*ti) - exp(R*delttiss)) ) ;
	}
      else //(ti > delttiss + tau)
	{
	  if (casl)  kctissue = 2 * T_1app * exp(-delttiss/T_1b) * exp(-(ti-tau-delttiss)/T_1app) * (1 - exp(-tau/T_1app));
	  else       kctissue = F/R * ( (exp(R*(delttiss+tau)) - exp(R*delttiss))  );
	}
  return kctissue;
}

  double TissueModel_nodisp_imperm::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
  // Tissue kinetic curve no dispersion impermeable vessel
  Tracer_Plus tr("OXASL::kctissue_nodisp_imperm");
  double kctissue = 0.0;

  //extract the pre-cap residence time
  double taup = residparam.AsScalar();

  if (ti>delttiss && ti<delttiss+taup+tau) {
    if (ti<delttiss+tau && ti < delttiss+taup) {
      if (casl) kctissue = exp(-delttiss/T_1b) - exp(-ti/T_1b);
      else      kctissue = ti - delttiss;
    }
    else if (ti<delttiss+tau && ti >= delttiss+taup) {
      if (casl) kctissue = exp(-delttiss/T_1b)*(1-exp(-taup/T_1b));
      else      kctissue = taup;
    }
    else if (ti >= delttiss+tau && ti >= delttiss+taup) {
      if (casl) kctissue = exp(-(ti-tau)/T_1b) - exp(-(delttiss+taup)/T_1b);
      else      kctissue = delttiss+tau+taup-ti;
    }
    else if (ti >= delttiss+tau && ti < delttiss+taup) {
      if (casl) kctissue = exp(-ti/T_1b)*(exp(tau/T_1b)-1);
      else      kctissue = tau;
    }
  }

  if (casl)  kctissue *= T_1b;
  else       kctissue *= exp(-ti/T_1b);
  kctissue *= 2;

  return kctissue;
}

double TissueModel_nodisp_2cpt::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
    // Two compartment model - the simplest form of the two cpt model
    // No backflow from tissue to blood
    // no venous outflow
    // From Parkes & Tofts and also St. Lawrence 2000 - both models are the same under these assumptions
    Tracer_Plus tr("OXASL::Tissuemodel_nodisp_twocpt");

    // extract residue function parameters
    double kw; //exchange rate = PS/vb
    kw = (residparam.Row(1)).AsScalar();
    //double PS; double vb;
    //PS = (residparam.Row(1)).AsScalar();
    //vb = (residparam.Row(2)).AsScalar();

    // calculate the residue function
    double a = kw + 1/T_1b; //alpha
    double b = (kw*T_1*T_1b) / (kw*T_1*T_1b + (T_1-T_1b) ); //beta
    double S;
    if (casl) S = 1/T_1;
    else      S = 1/T_1-1/T_1b;
    double T;
    if (casl) T = a;
    else      T = a - 1/T_1b;
    
    assert(!casl); //we only have pASL at the moment

    double kctissue = 0.0;
      if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue = 2*( b/S*exp(-ti/T_1)*(exp(S*ti)-exp(S*delttiss)) + (1-b)/T*exp(-a*ti)*(exp(T*ti)-exp(T*delttiss)) );
	}
      else if (ti > delttiss + tau)
	{
	  kctissue = 2*( b/S*exp(-ti/T_1)*(exp(S*(delttiss+tau))-exp(S*delttiss)) + (1-b)/T*exp(-a*ti)*(exp(T*(delttiss+tau))-exp(T*delttiss)) );
	}

      if (casl) kctissue *= exp(-delttiss/T_1b);

  return kctissue;
}

double TissueModel_nodisp_spa::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
    // Two compartment model
    // No backflow from tissue to blood
    // venous outflow (alhtough we dont model a venous component to the signal here)
    // St. Lawrence 2000
    Tracer_Plus tr("OXASL::Tissuemodel_nodisp_spa");

    assert(!casl);

    // extract residue function parameters
    double PS; double vb; double tauc;
    PS = (residparam.Row(1)).AsScalar();
    vb = (residparam.Row(2)).AsScalar(); //NB for SPA on the whole vb and PS appear together (as kw), but PS is on its own in ER
    tauc = (residparam.Row(3)).AsScalar();

    // calcualte residue function
   
    double kctissue = 0.0;
    if (ti>delttiss) {
      if (ti < delttiss+tauc && ti <delttiss + tau) {
	kctissue = Q(delttiss,ti,ti,PS,vb,tauc,fcalib,T_1,T_1b);
      }
      else if (ti > delttiss+tauc && ti < delttiss+tau) {
	kctissue = Q(delttiss,delttiss+tauc,ti,PS,vb,tauc,fcalib,T_1,T_1b) + R(delttiss+tauc,ti,ti,PS,vb,tauc,fcalib,T_1,T_1b);
      }
      else if (ti < delttiss + tauc && ti > delttiss+tau) {
	kctissue = Q(delttiss,delttiss+tau,ti,PS,vb,tauc,fcalib,T_1,T_1b);
      }
      else if (ti >= delttiss + tauc && ti >= delttiss+tau && ti < delttiss+tau+tauc){
	if (tauc<=tau) {
	  kctissue = Q(delttiss,delttiss+tauc,ti,PS,vb,tauc,fcalib,T_1,T_1b) + R(delttiss+tauc,delttiss+tau,ti,PS,vb,tauc,fcalib,T_1,T_1b);
	}
	else if (tau<tauc) {
	  kctissue = Q(delttiss,delttiss+tau,ti,PS,vb,tauc,fcalib,T_1,T_1b);
	}
      }
      else if (ti>=delttiss+tau+tauc) {
	kctissue = R(delttiss,delttiss+tau,ti,PS,vb,tauc,fcalib,T_1,T_1b);
      }
    }
    return kctissue;
}

  double TissueModel_nodisp_spa::Q(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const {
    Tracer_Plus tr("OXASL::TissueModel_nodisp_spa::Q");
    double a = PS/vb + 1/T_1b;
    double b = (PS*T_1*T_1b) / (PS*T_1*T_1b + (T_1-T_1b)*vb );
    double S = 1/T_1-1/T_1b;
    double T = PS/vb;
    return b/S*exp(-t3/T_1)*(exp(S*t2)-exp(S*t1)) + (1-b)/T*exp(-a*t3)*(exp(T*t2)-exp(T*t1));
  }

  double TissueModel_nodisp_spa::R(const double t1, const double t2, const double t3,const double PS, const double vb, const double tauc, const double fcalib, const double T_1, const double T_1b) const {
    Tracer_Plus tr("OXASL::TissueModel_nodisp_spa::R");
    double b = (PS*T_1*T_1b) / (PS*T_1*T_1b + (T_1-T_1b)*vb );
    double ER = 1 - exp(-PS/fcalib - (1/T_1b - 1/T_1)*tauc);
    double S = 1/T_1-1/T_1b;
    return b*ER/S*exp(-t3/T_1)*(exp(S*t2)-exp(S*t1));
  }

  double TissueModel_gammadisp_wellmix::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
  Tracer_Plus tr("OXASL::kctissue_gammadisp_wellmix");
  double kctissue = 0.0;

  assert(!casl); //only pASL at the moment!

  //extract dispersion parameters
  double s; double p;
  s = (dispparam.Row(1)).AsScalar();
  s = exp(s);
  double sp = (dispparam.Row(2)).AsScalar();
  sp = exp(sp);
  if (sp>10) sp=10;
  p = sp/s;

  double k=1+p*s;
  double T_1app = 1/( 1/T_1 + fcalib/lambda );
  double A = T_1app - T_1b;
  double B = A + s*T_1app*T_1b;
  if (B<1e-12) B=1e-12;  //really shouldn't happen, but combination of parameters may arise in artefactual voxels?
  double C = pow(s-1/T_1app+1/T_1b,p*s);
  if (s-1/T_1app+1/T_1b<=0) C=1e-12; //really shouldn't happen, but combination of parameters may arise in artefactual voxels?

  //cout << T_1app << " " << A << " " << B << " "<< C << " " << endl ;

      if(ti < delttiss)
	{ kctissue = 0;}
      else if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue = 2* 1/A * exp( -(T_1app*delttiss + (T_1app+T_1b)*ti)/(T_1app*T_1b) )*T_1app*T_1b*pow(B,-k)*
	    (  exp(delttiss/T_1app + ti/T_1b) * pow(s*T_1app*T_1b,k) * ( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
	       exp(delttiss/T_1b + ti/T_1app) * pow(B,k) * ( -1 + igamc(k,s*(ti-delttiss)) )  );  
	  
	}
      else //(ti > delttiss + tau)
	{
	  kctissue = 2* 1/(A*B) *
	    (  exp(-A/(T_1app*T_1b)*(delttiss+tau) - ti/T_1app)*T_1app*T_1b/C*
	       (  pow(s,k)*T_1app*T_1b* ( -1 + exp( (-1/T_1app+1/T_1b)*tau )*( 1 - igamc(k,B/(T_1app*T_1b)*(ti-delttiss)) ) +
					  igamc(k,B/(T_1app*T_1b)*(ti-delttiss-tau)) ) - 
		  exp( -A/(T_1app*T_1b)*(ti-delttiss-tau) )*C*B * ( igamc(k,s*(ti-delttiss-tau)) - igamc(k,s*(ti-delttiss)) ) )  );
	  
      //if (isnan(kctissue(it))) { kctissue(it)=0.0; cout << "Warning NaN in tissue KC"; }
	}
  //cout << kctissue.t() << endl;
  return kctissue;
}

/*
  double kctissue_gvf(const double ti,const double delttiss,const double tau, const double T_1b,const double T_1app,const double s,const double p) {
    // Tissue KC with a GVF AIF
    // NOTE: tau onyl scales the magnitude and does not affacet the overall shape in this model (see also kcblood_gvf)
  //Tracer_Plus tr("OXASL::kctissue_gvf");
  double kctissue = 0.0;

  double k=1+p*s;
  double A = T_1app - T_1b;
  double B = A + s*T_1app*T_1b;
  double C = pow(s-1/T_1app+1/T_1b,p*s);
  double sps = pow(s,k);

  if(ti < delttiss)
	{ kctissue = 0.0;}
      else //if(ti >= delttiss && ti <= (delttiss + tau))
	{
	  kctissue = 2* 1/(B*C) * exp(-(ti-delttiss)/T_1app)*sps*T_1app*T_1b * (1 - igamc(k,(s-1/T_1app-1/T_1b)*(ti-delttiss)));
	}
      // bolus duraiton is specified by the CVF AIF shape and is not an explicit parameter
      //else //(ti > delttiss + tau)
      //	{
      //	  kctissue(it) = exp(-(ti-delttiss-tau)/T_1app) * 2* 1/(B*C) * exp(-(delttiss+tau)/T_1app)*sps*T_1app*T_1b * (1 - igamc(k,(s-1/T_1app-1/T_1b)*(delttiss+tau)));
      //	}

  kctissue /= tau; // the aif is scaled by the duration of the bolus
  return kctissue;
}

  double kctissue_gaussdisp(const double ti,const double delttiss,const double tau,const double T_1b,const double T_1app,const double sig1,const double sig2) {
    // Tissue kinetic curve gaussian dispersion (pASL)
    // Hrabe & Lewis, MRM, 2004
    //Tracer_Plus tr("OXASL::kctissue_gaussdisp");
    double kctissue = 0.0;

    double R = 1/T_1app - 1/T_1b; 
    double sqrt2 = sqrt(2);
    double F = 2 * exp(-ti/T_1app);
    double u1 = (ti-delttiss)/(sqrt2*sig1);
    double u2 = (ti - delttiss - tau)/(sqrt2*sig2);
    
    kctissue = F/(2*R) * (  (erf(u1) - erf(u2))*exp(R*ti) 
				  - (1 + erf(u1 - (R*sig1)/sqrt2))*exp(R*(delttiss+(R*sig1*sig1)/2)) 
				  + (1 + erf(u2 - (R*sig2)/sqrt2))*exp(R*(delttiss+tau+(R*sig2*sig2)/2)) 
				  );
    
  return kctissue;
}
  */

 double TissueModel_aif_residue::kctissue(const double ti,const double fcalib, const double delttiss,const double tau,const double T_1b,const double T_1, const double lambda,const bool casl, const ColumnVector dispparam, const ColumnVector residparam) const {
  Tracer_Plus tr("OXASL::kctissue_aif_residue");
  double kctissue = 0.0;

  // calculate the appropraite time series for the aif and residue
  // assume that the aif is zero before TI=0 (true of most aif models under sensible parameter values)

  // OLD fixed number of discrete time points
  // int nint = 10; //number of intervals for simpsons rule;
  // int niti = 2*nint+1; //the number of time points to use in the integration - for simpsons rule follow the patttern (2*x)+1
  // double delti = ti/(niti-1);

  // fixed time intervals for discretization (these discrete time points are termed iti)
  double delti=0.1; // time interval for the iti
  double dti =  fmod(ti,delti); // the time not covered by the iti
  //double maxiti = ti - dti; //largest iti
  int niti = floor(ti/delti)+1; // number of iti (include the itit at 0)
 
  ColumnVector aifts(niti);
  ColumnVector rests(niti);
  double iti;
  //calacute aif and residue function for all the iti
  for (int i=0; i<niti; i++) { //start from iTI = 0 and go up to iTI=niti*delti)
    iti = i*delti;
    aifts(i+1) = aifmodel->kcblood(iti,delttiss,tau,T_1b,casl,dispparam);
    rests(i+1) = residmodel->resid(iti+dti,fcalib,T_1,T_1b,lambda,residparam); //NB offset on this ready for when we reverse it
  }
  //cacualte the aif and residue at the TI
  double aifti = aifmodel->kcblood(ti,delttiss,tau,T_1b,casl,dispparam);
  //double resti = residmodel->resid(ti,fcalib,T_1,lambda,residparam);

  // do convolution to get the value at this TI
  ColumnVector prod(niti);
  prod = SP(aifts,rests.Reverse());
  double prodti = aifti; //becasue resti at time zero is 1

  string integrate = "trapezoid";
  //if (integrate == "rect") {}
  // rect intergration
    //nothing to do here - this is essentially done as the default option
  if (integrate == "trapezoid") {
    ColumnVector trap(niti);
    trap = 1;
    trap(1) = 0.5;
    trap(niti) = 0.5;
    prod = SP(prod,trap);
  }
  else if (integrate == "simpson") {
  // simpson
    //work out how many time points we can apply simpson to
    int nsimp = floor( (niti-1)/2 )*2 +1;

    ColumnVector simpson(niti);
    simpson(1) = 1/3;
    for (int h=2; h<nsimp; h += 2) {
      simpson(h) = 4/3;
      simpson(h+1) = 2/3;
    }
    simpson(nsimp) = 1/3;

    if (nsimp<niti) {
      //we still have another time point to incorporate do this using trapezium rule
      assert(niti-nsimp==1);
      simpson(nsimp) = simpson(nsimp)+0.5;
      simpson(niti) = 0.5;
    }
    prod = SP(prod,simpson);
  }

    kctissue = prod.Sum()*delti; //NB must be multiplied by the timespacing
    kctissue += (0.5*prodti+0.5*prod(niti))*dti; // plus the last bit (which we will do with trapezium rule)

  return kctissue;
 }

// --- useful general functions ---
  double icgf(const double a,const double x) {
    Tracer_Plus tr("OXASL::icgf");
    
    //incomplete gamma function with a=k, based on the incomplete gamma integral
    
    return MISCMATHS::gamma(a)*igamc(a,x);
  }
  
  double gvf(const double t,const double s,const double p) {
    Tracer_Plus tr("OXASL::gvf");
    
    //The Gamma Variate Function (correctly normalised for area under curve) 
    // Form of Rausch 2000
    // NB this is basically a gamma pdf
    
    if (t<0)    return 0.0;
    else        return pow(s,1+s*p) / MISCMATHS::gamma(1+s*p) * pow(t,s*p) * exp(-s*t);
  }

  double numerical_integration(ColumnVector integrand, double del, double finalval, double finaldel, string method) {
    Tracer_Plus tr("OXASL::numerical_integration");

    int ndel;
    ndel = integrand.Nrows();
    ColumnVector prod(ndel);
    prod = integrand;

    // do numerical intergration on a supplied equispaced vector - with possible extra point at end with spacing: finaldel < del
    if (method == "rect") {}
    // rect intergration
    //nothing to do here - this is essentially done as the default option
    else if (method == "trapezoid") {
      ColumnVector trap(ndel);
      trap = 1;
      trap(1) = 0.5;
      trap(del) = 0.5;
      prod = SP(prod,trap);
    }

    double result;
    result = prod.Sum()*del;
    if (finaldel>1e-12) {
      result += (0.5*finalval+0.5*prod(ndel))*finaldel;
    }

    return result;

  }
} //end namespace
