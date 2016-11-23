/*  convergence.cc - Convergence detectors for FABBER

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "convergence.h"

#include "easylog.h"
#include "dataset.h"

#include "assert.h"
#include <iostream>

using namespace std;

ConvergenceDetector* ConvergenceDetector::NewFromName(const string& name)
{
	ConvergenceDetectorFactory* factory = ConvergenceDetectorFactory::GetInstance();
	ConvergenceDetector* conv = factory->Create(name);
	if (conv == NULL)
	{
		throw Invalid_option("Unrecognized convergence detector: " + name);
	}
	return conv;
}

void CountingConvergenceDetector::Initialize(FabberRunData &params)
{
	m_max_its = convertTo<int> (params.GetStringDefault("max-iterations", "10"));
	LOG << "CountingConvergenceDetector::Max iterations=" << m_max_its << endl;
	if (m_max_its <= 0)
		throw Invalid_option("CountingConvergenceDetector::Max iterations must be positive");
	Reset();
}

bool CountingConvergenceDetector::Test(double)
{
	++m_its;
	if (m_its >= m_max_its)
	{
		m_reason = "Max iterations reached";
		return true;
	}
	else
	{
		return false;
	}
}

void CountingConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
	out << indent << "Starting iteration " << m_its + 1 << " of " << m_max_its << endl << endl;
}

void CountingConvergenceDetector::Reset(double F)
{
	m_its = 0;
	m_reason = "";
}

void FchangeConvergenceDetector::Initialize(FabberRunData &params)
{
	CountingConvergenceDetector::Initialize(params);

	m_min_fchange = convertTo<double> (params.GetStringDefault("min-fchange", "0.01"));
	LOG << "FchangeConvergenceDetector::Minimum F change=" << m_min_fchange << endl;
	if (m_min_fchange <= 0)
		throw Invalid_option("FchangeConvergenceDetector::Minimum F change must be positive");
	Reset();
}

void FchangeConvergenceDetector::Reset(double F)
{
	CountingConvergenceDetector::Reset();
	m_prev_f = F;
	m_save = false;
	m_revert = false;
}

bool FchangeConvergenceDetector::Test(double F)
{
	//    if (F == 1234.5678)
	//        throw logic_error("FchangeConvergenceDetector needs F, but it seems it isn't being calculated!  Internal bug... should have needF = true");
	// F could actually be 1234.5678, but what are the chances of that?
	double diff = F - m_prev_f;
	m_prev_f = F;

	diff = diff > 0 ? diff : -diff;
	if (diff < m_min_fchange)
	{
		m_reason = "Absolute difference less than minimum";
		return true;
	}
	else
		return CountingConvergenceDetector::Test(F);
}

void FchangeConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
	out << indent << "Iteration " << m_its << " of at most " << m_max_its << endl;
	out << indent << "Previous Free Energy == " << m_prev_f << endl;
}

void FreduceConvergenceDetector::Initialize(FabberRunData &params)
{
	FchangeConvergenceDetector::Initialize(params);
	Reset();
}

bool FreduceConvergenceDetector::Test(double F)
{
	double diff = F - m_prev_f;

	if (diff < 0)
	{
		m_reason = "F reduced";
		m_revert = true;
		return true;
	}
	else
	{
		return FchangeConvergenceDetector::Test(F);
	}
}

void FreduceConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
	out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << m_reason << endl;
	out << indent << "Previous Free Energy == " << m_prev_f << endl;
}

void TrialModeConvergenceDetector::Initialize(FabberRunData &params)
{
	FchangeConvergenceDetector::Initialize(params);

	// FIXME for consistency with previous versions but not really correct!
	m_max_its += 1;

	m_max_trials = convertTo<int> (params.GetStringDefault("max-trials", "10"));
	LOG << "TrialModeConvergenceDetector::Max trials=" << m_max_trials << endl;
	if (m_min_fchange <= 0)
		throw Invalid_option("TrialModeConvergenceDetector::Max trials must be positive");
	Reset();
}

void TrialModeConvergenceDetector::Reset(double F)
{
	FchangeConvergenceDetector::Reset();
	m_trials = 0;
	m_save = true; // FIXME why?
	m_trialmode = false;
}

bool TrialModeConvergenceDetector::Test(double F)
{
	double diff = F - m_prev_f;

	double absdiff = diff; //look at magnitude of diff in absdiff
	if (diff < 0)
	{
		absdiff = -diff;
	}

	//if we are not in trial mode
	if (!m_trialmode)
	{
		// if F has reduced then we will enter trial mode
		if (diff < 0)
		{
			++m_trials;
			m_trialmode = true;
			m_revert = true;
			m_save = false;
		}
		return FchangeConvergenceDetector::Test(F);
	}
	// if we are in trial mode
	else
	{
		++m_trials;

		// if F has improved over our previous best then resume iterations
		// FIXME not necessarily improved over previous best, just over
		// previous decrease. Should save=true?
		if (diff > 0)
		{
			m_trialmode = false;
			m_save = true;
			m_revert = false;
			m_trials = 0; //reset number of trials for future use
			return FchangeConvergenceDetector::Test(F);
		}
		//if we have exceeded max trials then stop and output previous best result
		else if (m_trials >= m_max_trials)
		{
			m_reason = "Reached max trials";
			return true;
		}
		// otherwise continue in trial mode for time being. FIXME do we want
		// to stop if difference is < 0 but very small?
		else
		{
			return FchangeConvergenceDetector::Test(F);
		}
	}
}

void TrialModeConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
	out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << m_reason << endl;
	out << indent << "Previous Free Energy == " << m_prev_f << endl;
}

void LMConvergenceDetector::Initialize(FabberRunData &params)
{
	m_max_its = convertTo<int> (params.GetStringDefault("max-iterations", "10"));
	LOG << "LMConvergenceDetector::Max iterations=" << m_max_its << endl;
	if (m_max_its <= 0)
		throw Invalid_option("LMConvergenceDetector::Max iterations must be positive");

	m_max_fchange = convertTo<double> (params.GetStringDefault("max-fchange", "0.01"));
	LOG << "LMConvergenceDetector::Max fchange=" << m_max_fchange << endl;
	if (m_max_fchange <= 0)
		throw Invalid_option("LMConvergenceDetector::Max fchange must be positive");
	Reset();
}

void LMConvergenceDetector::Reset(double F)
{
	m_its = 0;
	prev = F;
	m_save = true;
	m_revert = false;
	alphastart = 1e-6;
	alpha = 0.0;
	alphamax = 1e6;
	LM = false;
}

bool LMConvergenceDetector::Test(double F)
{
	double diff = F - prev;

	// cout << "----" << endl;
	// cout << m_its << endl;
	// cout << "F:" << F << endl;
	// cout << "prev:" << prev << endl;
	// cout << diff << endl;

	double absdiff = diff; //look at magnitude of diff in absdiff
	if (diff < 0)
	{
		absdiff = -diff;
	}

	//if we are not in LM mode
	if (!LM)
	{
		// if F has reduced then we go into LM mode
		if (diff < 0)
		{
			LM = true;
			m_revert = true; //revert to the previous solution and try again with LM adjustment
			alpha = alphastart;
			//cout << "Entering LM" << endl;
			return false;
		}
		//otherwise if we have converged stop
		else if (absdiff < m_max_fchange)
		{
			reason = "F converged";
			m_revert = false;
			return true;
		}
		//otherwise if we have reached max iterations stop
		else if (m_its >= m_max_its)
		{
			reason = "Max iterations reached";
			m_revert = false;
			return true;
		}
		// otherwise carry on to next iteration
		else
		{
			prev = F;
			++m_its;
			return false;
			m_revert = false;
		}
	}
	// if we are in LM mode (NB we dont increase iterations if we are increasing the LM alpha, onyl when we make a sucessful step)
	else
	{
		// if F has improved over our previous best then reduce the alpha and continue from current estimate

		if (diff > 0)
		{
			if (alpha == alphastart)
			{ //leave LM mode if alpha returns to inital value
				LM = false;
			}
			else
			{
				alpha /= 10;
				LM = true;
			}
			m_revert = false;
			prev = F;
			cout << "Reducing LM" << endl;
			++m_its; // if F has improved then we take this as the new 'step' and move onto the next iteration
			return false;
		}
		//if alpha gets too large then we cannot achieve any gain here so stop and revert to previous best solution
		else if (alpha >= alphamax)
		{
			reason = "Reached max alpha";
			m_revert = true;
			//cout << "LM maxed out" << endl;
			return true;
		}
		//otherwise if we have reached max iterations stop
		else if (m_its >= m_max_its)
		{
			reason = "Max iterations reached";
			m_revert = false;
			return true;
		}
		// otherwise continue in LM mode for time being, try increasing alpha
		else
		{
			alpha *= 10;
			m_revert = true; //revert to the previous solution and try again with new LM adjustment
			//cout << "Increasing LM" << endl;
			return false;
		}
	}
}

void LMConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
	out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << reason << endl;
	out << indent << "Previous Free Energy == " << prev << endl;
}
