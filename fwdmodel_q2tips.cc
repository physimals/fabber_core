/*  fwdmodel_q2tips.cc - Implements a Q2TIPS dual-echo ASL model

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_q2tips.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

FactoryRegistration<FwdModelFactory,  Q2tipsFwdModel> 
  Q2tipsFwdModel::registration("q2tips-dualecho");

FwdModel* Q2tipsFwdModel::NewInstance()
{
  return new Q2tipsFwdModel();
}

string Q2tipsFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_q2tips.cc,v 1.2 2007/08/02 15:14:11 adriang Exp $ and "
	+ Quipss2FwdModel::ModelVersion();
}

void Q2tipsFwdModel::Initialize(ArgsType& args)
{
  Quipss2FwdModel::Initialize(args);
}

void Q2tipsFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
    Tracer_Plus tr("Q2tipsFwdModel::Evaluate");
    // Adapted from original_fwdmodel.m
    
    // Parameterization used in most recent results:
    // Absolute M and Q change (same units as M0 or Q0):
    ColumnVector StatMag = params(M0index()) - Mbasis * MnOf(params);
    ColumnVector CBF = params(Q0index()) + Qbasis * QnOf(params);
    // Fractional change in BOLD effect (at TE_2), rather than using % R2* change 
    ColumnVector R2s = -1/echoTime(2) * log( 
        Rbasis * RnOf(params) + exp(-echoTime(2)*params(R0index())));

    // The following are relative magnetizations    
    double pretag = 1; // untagged
    double T1b = (stdevT1b>0 ? params(T1bIndex()) : fixedT1b);
    double invEfficiency = (stdevInvEff>0 ? params(InvEffIndex()) : fixedInvEff);
    ColumnVector bolus = 1 - (1-rho)*invEfficiency*exp(-TI2/T1b); // tag or control
    //    double posttag = 1 - exp(-(TI2-TI1)/T1b); // saturated
    double dt = (stdevDt>0? params(dtIndex()) : fixedDt);
    //    ColumnVector Sb = SP( CBF,  // SP(a,b) means a.*b 
    //        pretag*dt + bolus*TI1 + posttag*(TI2-TI1-dt) );
    double posttagQ2 = ((TI2-TI1-dt) 
			+ T1b*exp(-(TI2-TI1)/T1b)
			- T1b*exp(-dt/T1b));
    ColumnVector Sb = SP( CBF,  // SP(a,b) means a.*b 
        pretag*dt + bolus*TI1 + posttagQ2 );
    ColumnVector S = StatMag + Sb;
      
  int Ntimes = R2s.Nrows();
  int Nte = echoTime.Nrows();
  if (result.Nrows() != Nte*Ntimes)
    result.ReSize(Nte*Ntimes);
    
//  result = 0.0/0.0; // pre-fill with nans to check all overwritten
  
  for (int te = 1; te <= Nte; te++)
    {
      ColumnVector nuisance = Nbasis * NnOf(te, params);
      // Will be all-zero if there are no nuisance regressors
        
      for (int i = 1; i <= Ntimes; i++)
        result( Nte*(i-1) + te ) = 
            S(i) * exp(-echoTime(te) * R2s(i)) + nuisance(i);
      // Fill order: te1 te2 te1 te2 te1 te2 te1 te2 ...
    }

  return; // answer is in the "result" vector
}

vector<string> Q2tipsFwdModel::GetUsage() const {
  // Uses exactly the same options as its parent.
  return Quipss2FwdModel::GetUsage();
}
