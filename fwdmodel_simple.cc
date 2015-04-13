/*  fwdmodel_simple.cc - Implements the simplified ASL model

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include "fwdmodel_simple.h"

#include <iostream>
#include <newmatio.h>
#include <stdexcept>
#include "newimage/newimageall.h"
using namespace NEWIMAGE;
#include "easylog.h"

string SimpleFwdModel::ModelVersion() const
{
  return "$Id: fwdmodel_simple.cc,v 1.19 2012/01/13 12:00:59 adriang Exp $";
}

void SimpleFwdModel::Evaluate(const ColumnVector& params, ColumnVector& result) const
{
//  ColumnVector Qn = id.Qn(params);
  ColumnVector tcModulated = id.Q0(params) 
	* ( 1 + basis * id.Qn(params) / 100.0);
  ColumnVector tcUnmodulated = id.U0(params) 
	* ( 1 + basis * id.Un(params) / 100.0);
  ColumnVector R2s = id.R0(params) * (1 + basis * id.Rn(params) / 100.0);
  ColumnVector S = SP(rho, tcModulated) + tcUnmodulated; // SP means .*
  
  int Ntimes = R2s.Nrows();
  if (result.Nrows() != 2*Ntimes)
    result.ReSize(2*Ntimes);
    
//  result = 0.0/0.0; // pre-fill with nans to check all overwritten
  
  for (int te = 1; te <= echoTime.Nrows(); te++)
    {
      for (int i = 1; i <= Ntimes; i++)
        result( 2*i - 2 + te ) = S(i) * exp(-echoTime(te) * R2s(i));
      // Fill order: te1 te2 te1 te2 te1 te2 te1 te2 ...
    }

/*    
  LOG << "Fwdmodel input:\n" << params;
  LOG << "tcModulated:" << endl;
  LOG << tcModulated;
  LOG << "tcUnmodulated:" << endl;
  LOG << tcUnmodulated;
  LOG << "R2s:" << endl;
  LOG << R2s;
  LOG << "S" << endl << S;
  LOG << "echoTime" << echoTime;
  LOG << "Output:\n" << result;
*/   
//    LOG << "BASIS\n" << basis.t()*basis;

  return;
}

SimpleFwdModel::SimpleFwdModel(ArgsType& args)
{
    
  string scanParams = args.ReadWithDefault("scan-params","hardcoded");
  if (scanParams == "hardcoded")
    {
      LOG << "  Loading hardcoded forward-model priors" << endl;
      echoTime.ReSize(2);
      echoTime << 9.1 << 30;
      echoTime *= 0.001;

      //      LOG << "SimpleFwdModel::SimpleFwdModel isn't implemented yet!" << endl;

      // basis.. not even going to bother trying.  Load from a file.
//      string basisFile = "/home/fs0/adriang/proj/response_fir/bolddesign.mat";
      string basisFile = "/usr/people/woolrich/scratch/tldata/analysis_protocols/response_fromroi/cbvdesign.mat";
      LOG << "    Reading basis functions from file: " << basisFile << endl;
      basis = read_vest_fabber(basisFile);
      //      LOG << "basis == \n" << basis << endl;
      // Nrows = 136, Ncols = 15
      // LOG << "BASIS: " << basis.Nrows() << basis.Ncols() << endl;

      LOG << "      Read " << basis.Ncols() << " basis functions" << endl;
      id.Define(basis.Ncols());

      rho.ReSize(68*2);
      rho(1) = -1;
      for (int i = 2; i <= 68*2; i++)
        rho(i) = rho(i-1) * -1;

      //      LOG << "echoTime:\n" << echoTime << endl;
      //      LOG << "rho:\n" << rho << endl;

    } 
  else
    throw Invalid_option("Only --scan-params=hardcoded is accepted at the moment");
}

void SimpleFwdModel::DumpParameters(const ColumnVector& vec,
                                    const string& indent) const
{
    LOG << indent << "Baseline parameters:" << endl;
    LOG << indent << "  U0 == " << id.U0(vec) << " (baseline unmodulated mag.)\n";
    LOG << indent << "  Q0 == " << id.Q0(vec) << " (baseline modulated mag.)\n";
    LOG << indent << "  R0 == " << id.R0(vec) << " (baseline T2*)\n";
    LOG << indent << "Percent change parameters:" << endl;
    LOG << indent << "  Un == " << id.Un(vec).t();// << "]\n";
    LOG << indent << "  Qn == " << id.Qn(vec).t();// << "]\n";
    LOG << indent << "  Rn == " << id.Rn(vec).t();// << "]\n";
}    

void SimpleFwdModelIdStruct::NameParams(vector<string>& names) const
{
    names.clear();
    
    for (int p = 1; p <= 3; p++)
    {
        string letter;
        switch(p) { 
            case 1: letter="Q"; break;
            case 2: letter="U"; break;
            case 3: letter="R"; break;
        } 
        names.push_back(letter + "0");
        for (int i = 1; i <= Nbasis; i++)
        {
            names.push_back(letter + "_percchg_" + stringify(i));
        }
    }
    assert(names.size() == unsigned(3+3*Nbasis)); 
}
