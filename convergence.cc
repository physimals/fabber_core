/*  convergence.cc - Convergence detectors for FABBER

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "convergence.h"

#include "easylog.h"
#include "rundata.h"

#include <iostream>
#include <string>

using namespace std;

ConvergenceDetector *ConvergenceDetector::NewFromName(const string &name)
{
    ConvergenceDetectorFactory *factory = ConvergenceDetectorFactory::GetInstance();
    ConvergenceDetector *conv = factory->Create(name);
    if (conv == NULL)
    {
        throw InvalidOptionValue("convergence", name, "Unrecognized convergence detector");
    }
    return conv;
}

void ConvergenceDetector::Initialize(FabberRunData &params)
{
    m_log = params.GetLogger();
}
void CountingConvergenceDetector::Initialize(FabberRunData &params)
{
    ConvergenceDetector::Initialize(params);
    m_max_its = convertTo<int>(params.GetStringDefault("max-iterations", "10"));
    if (m_max_its <= 0)
        throw InvalidOptionValue("max_iterations", stringify(m_max_its), "Must be positive");
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

void CountingConvergenceDetector::Dump(ostream &out, const string &indent) const
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

    m_min_fchange = convertTo<double>(params.GetStringDefault("min-fchange", "0.01"));
    if (m_min_fchange <= 0)
        throw InvalidOptionValue("min-fchange", stringify(m_min_fchange), "Must be positive");
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
    //        throw logic_error("FchangeConvergenceDetector needs F, but it seems it isn't being
    //        calculated!  Internal bug... should have needF = true");
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

void FchangeConvergenceDetector::Dump(ostream &out, const string &indent) const
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

void FreduceConvergenceDetector::Dump(ostream &out, const string &indent) const
{
    out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << m_reason
        << endl;
    out << indent << "Previous Free Energy == " << m_prev_f << endl;
}

void TrialModeConvergenceDetector::Initialize(FabberRunData &params)
{
    FchangeConvergenceDetector::Initialize(params);

    // FIXME for consistency with previous versions but not really correct!
    m_max_its += 1;

    m_max_trials = convertTo<int>(params.GetStringDefault("max-trials", "10"));
    if (m_max_trials <= 0)
        throw InvalidOptionValue("max-trials", stringify(m_max_trials), "Must be positive");

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

    // if we are not in trial mode
    if (!m_trialmode)
    {
        if (diff < 0)
        {
            // if F has reduced then we will enter trial mode. Don't overwrite 
            // our best F so far
            ++m_trials;
            m_trialmode = true;
            m_revert = true;
            m_save = false;
            return false;
        }
        else {
            double absdiff = diff > 0 ? diff : -diff;
            if (absdiff < m_min_fchange) {
                // F has increased by less than tolerance - stop here and don't revert
                m_reason = "F increased by less than tolerance";
                m_revert = false;
                m_save = false;
                return true;
            }
            else {
                // F has increased by more than tolerance - continue, and save as best so far 
                m_save = true;
                m_revert = false;
                m_prev_f = F;
                return false;
            }
        }
    }
    // if we are in trial mode
    else
    {
        ++m_trials;

        if (diff > 0)
        {
            double absdiff = diff > 0 ? diff : -diff;
            if (absdiff < m_min_fchange) {
                // if F has improved over our previous best and by less than tolerance
                // - stop here and don't revert
                m_reason = "F increased by less than tolerance during trial mode";
                m_revert = false;
                m_save = false;
                return true;
            }
            else {
                // F has increased over previous best by more than tolerance - continue iterations
                // and save as best so far 
                m_trialmode = false;
                m_trials = 0;
                m_save = true;
                m_revert = false;
                m_prev_f = F;
                return false;
            }
        }
        else if (m_trials >= m_max_trials)
        {
            // if we have exceeded max trials then stop and revert to previous best result
            m_reason = "Reached max trials";
            m_save = false;
            m_revert = true;
            return true;
        }
        else
        {
            // otherwise continue in trial mode for time being.
            m_save = false;
            m_revert = false;
            return false;
        }
    }
}

void TrialModeConvergenceDetector::Dump(ostream &out, const string &indent) const
{
    out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << m_reason
        << endl;
    out << indent << "Previous Free Energy == " << m_prev_f << endl;
}

void LMConvergenceDetector::Initialize(FabberRunData &params)
{
    ConvergenceDetector::Initialize(params);

    m_max_its = convertTo<int>(params.GetStringDefault("max-iterations", "10"));
    if (m_max_its <= 0)
        throw InvalidOptionValue("max-iterations", stringify(m_max_its), "Must be positive");

    m_max_fchange = convertTo<double>(params.GetStringDefault("max-fchange", "0.01"));
    if (m_max_fchange <= 0)
        throw InvalidOptionValue("max-fchange", stringify(m_max_fchange), "Must be positive");
    Reset();
}

void LMConvergenceDetector::Reset(double F)
{
    m_its = 0;
    m_prev = F;
    m_save = true;
    m_revert = false;
    m_alphastart = 1e-6;
    m_alpha = 0.0;
    m_alphamax = 1e6;
    m_LM = false;
}

bool LMConvergenceDetector::Test(double F)
{
    double diff = F - m_prev;

    // cout << "----" << endl;
    // cout << m_its << endl;
    // cout << "F:" << F << endl;
    // cout << "prev:" << prev << endl;
    // cout << diff << endl;

    double absdiff = diff; // look at magnitude of diff in absdiff
    if (diff < 0)
    {
        absdiff = -diff;
    }

    // if we are not in LM mode
    if (!m_LM)
    {
        // if F has reduced then we go into LM mode
        if (diff < 0)
        {
            m_LM = true;
            m_revert = true; // revert to the previous solution and try again with LM adjustment
            m_alpha = m_alphastart;
            // cout << "Entering LM" << endl;
            return false;
        }
        // otherwise if we have converged stop
        else if (absdiff < m_max_fchange)
        {
            m_reason = "F converged";
            m_revert = false;
            return true;
        }
        // otherwise if we have reached max iterations stop
        else if (m_its >= m_max_its)
        {
            m_reason = "Max iterations reached";
            m_revert = false;
            return true;
        }
        // otherwise carry on to next iteration
        else
        {
            m_prev = F;
            ++m_its;
            return false;
        }
    }
    // if we are in LM mode (NB we dont increase iterations if we are increasing the LM m_alpha,
    // onyl when we make a sucessful step)
    else
    {
        // if F has improved over our previous best then reduce the m_alpha and continue from
        // current estimate

        if (diff > 0)
        {
            if (m_alpha == m_alphastart)
            { // leave LM mode if m_alpha returns to inital value
                m_LM = false;
            }
            else
            {
                m_alpha /= 10;
                m_LM = true;
            }
            m_revert = false;
            m_prev = F;
            cout << "Reducing LM" << endl;
            ++m_its; // if F has improved then we take this as the new 'step' and move onto the next
                     // iteration
            return false;
        }
        // if m_alpha gets too large then we cannot achieve any gain here so stop and revert to
        // previous best solution
        else if (m_alpha >= m_alphamax)
        {
            m_reason = "Reached max m_alpha";
            m_revert = true;
            // cout << "LM maxed out" << endl;
            return true;
        }
        // otherwise if we have reached max iterations stop
        else if (m_its >= m_max_its)
        {
            m_reason = "Max iterations reached";
            m_revert = false;
            return true;
        }
        // otherwise continue in LM mode for time being, try increasing m_alpha
        else
        {
            m_alpha *= 10;
            m_revert = true; // revert to the previous solution and try again with new LM adjustment
            // cout << "Increasing LM" << endl;
            return false;
        }
    }
}

void LMConvergenceDetector::Dump(ostream &out, const string &indent) const
{
    out << indent << "Iteration " << m_its << " of at most " << m_max_its << " : " << m_reason
        << endl;
    out << indent << "Previous Free Energy == " << m_prev << endl;
}
