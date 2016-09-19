/*  convergence.h - Convergence detectors for FABBER

 Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include <ostream>

#include "utils.h"
#include "dataset.h"

/**
 * Abstract base class for method of testing whether the free energy maximisation algorithm has converged.
 */
class ConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance();

	static ConvergenceDetector* NewFromName(const string& name);

	/**
	 * Initialize from run parameters
	 */
	virtual void Initialize(FabberRunData &params)=0;

	/**
	 * The key method. Called before iteration with the current free energy
	 *
	 * Returns true if converged, false otherwise
	 */
	virtual bool Test(double F) = 0;
	virtual ~ConvergenceDetector()
	{
	}

	/**
	 * Send information on current progress to output stream
	 */
	virtual void
	DumpTo(std::ostream& out, const std::string indent = "") const = 0;

	/**
	 * Reset as if algorithm was starting from scratch
	 */
	virtual void Reset(double F = -99e99) = 0;

	virtual bool UseF() const = 0;
	virtual bool NeedSave() = 0;
	virtual bool NeedRevert() = 0;
	virtual float LMalpha() = 0;

	/**
	 * Reason convergence reached
	 *
	 * If Test returns true, this should
	 * contain a human readable string giving the reason
	 */
	virtual std::string GetReason() = 0;
};

/**
 * Simple implementation which just carries out a fixed number
 * of iterations
 */
class CountingConvergenceDetector: public ConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance()
	{
		return new CountingConvergenceDetector();
	}

	virtual bool Test(double);

	virtual void Initialize(FabberRunData &params);

	virtual void DumpTo(std::ostream& out, const std::string indent = "") const;

	/**
	 * Set the number of iterations back to zero
	 */
	virtual void Reset(double F = -99e99);

	/**
	 * Does not use the free energy
	 *
	 * @return false
	 */
	virtual bool UseF() const
	{
		return false;
	}

	/**
	 * Do we need to save the last set of parameters?
	 *
	 * @return true if we do, i.e. if the last value
	 *              of F tested was the best so far.
	 */
	virtual bool NeedSave()
	{
		return false;
	}

	/**
	 * Do we need to revert to the previously saved set of parameters?
	 */
	virtual bool NeedRevert()
	{
		return false;
	}

	/**
	 * Used by the LM detector - all others return zero
	 */
	virtual float LMalpha()
	{
		return 0.0;
	}

	/**
	 * Reason convergence reached
	 *
	 * If Test returns true, this should
	 * contain a human readable string giving the reason
	 */
	std::string GetReason()
	{
		return m_reason;
	}

protected:
	int m_its;
	int m_max_its;
	std::string m_reason;
};

/**
 * Converges when the absolute difference between F and the previous
 * value is sufficiently small
 */
class FchangeConvergenceDetector: public CountingConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance()
	{
		return new FchangeConvergenceDetector();
	}
	/**
	 * @return true if F differs from previous value by less than
	 * the configured value
	 */
	virtual bool Test(double F);

	/**
	 * @param maxIts Maximum number of iterations
	 * @param Fchange Convergence when difference between subsequent Fs is less than this
	 */
	virtual void Initialize(FabberRunData &params);

	virtual void DumpTo(std::ostream& out, const std::string indent = "") const;

	/**
	 * Sets number of iterations to zero and the previous free energy to
	 * an unrealistically large number
	 */
	virtual void Reset(double F = -99e99);

	/**
	 * Uses the free energy
	 * @return true
	 */
	virtual bool UseF() const
	{
		return true;
	}

	virtual bool NeedSave()
	{
		return m_save;
	}

	virtual bool NeedRevert()
	{
		return m_revert;
	}

protected:
	double m_prev_f;
	double m_min_fchange;
	bool m_revert;// determines whether we should revert or not if asked
	bool m_save;
};

/**
 * Converges when the absolute difference between F and the previous
 * value is sufficiently small or if F has reduced since the
 * previous iteration
 */
class FreduceConvergenceDetector: public FchangeConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance()
	{
		return new FreduceConvergenceDetector();
	}
	/**
	 * @return true if difference to previous is less than
	 * configured value, or F is less than previous value
	 */
	virtual bool Test(double F);

	/**
	 * @param maxIts Maximum number of iterations
	 * @param Fchange Change if F smaller than this amount means convergence
	 */
	virtual void Initialize(FabberRunData &params);

	virtual void DumpTo(std::ostream& out, const std::string indent = "") const;

protected:
};

/**
 * Convergence detector which gives F a chance to increase
 * after first decrease found
 *
 * Like FreduceConvergenceDetector, however if F is found
 * to have reduced, continues for a maximum number of additional
 * iterations. If F does not go above the previous highest in this
 * time, then convergence is reached. If it does, carries on
 * until either F reduces again, or the absolute difference from
 * the previous is sufficiently small
 */
class TrialModeConvergenceDetector: public FchangeConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance()
	{
		return new TrialModeConvergenceDetector();
	}
	virtual bool Test(double F);
	virtual void Initialize(FabberRunData &params);
	virtual void DumpTo(std::ostream& out, const std::string indent = "") const;
	virtual void Reset(double F = -99e99);

	virtual bool NeedSave()
	{
		return m_save;
	}

protected:
	int m_trials;
	int m_max_trials;
	bool m_trialmode;
};

/**
 * Convergence detector which gives a variable amount of time
 * for F to increase after a decrease.
 *
 * FIXME should subclass FchangeConvergenceDetector as
 * well but not done yet as do not have tests for this
 * class
 */
class LMConvergenceDetector: public ConvergenceDetector
{
public:
	static ConvergenceDetector *NewInstance()
	{
		return new LMConvergenceDetector();
	}
	/**
	 * Convergence is reached if maximum number
	 * of iterations is reached, if the change in F to
	 * the previous step is smaller than Fchange or if
	 * a series of decreases are found.
	 *
	 * When a decrease in F is found, the detector goes
	 * into LM mode, and each decrease is counted.
	 * If an increase is found, one decrease is wiped out.
	 * If the number of decreases returns to zero, the
	 * detector leaves LM mode and continues. If the
	 * number of decreases reaches a maximum, the
	 * detector converges.
	 *
	 * @return true if convergence reached
	 */
	virtual bool Test(double F);
	virtual void Initialize(FabberRunData &params);
	virtual void DumpTo(std::ostream& out, const std::string indent = "") const;
	virtual void Reset(double F = -99e99);
	virtual bool UseF() const
	{
		return true;
	}
	bool NeedSave()
	{
		return m_save;
	}
	bool NeedRevert()
	{
		return m_revert;
	}
	float LMalpha()
	{
		return alpha;
	}
	std::string GetReason()
	{
		return reason;
	}

private:
	int m_its;
	int m_max_its;
	double prev;
	double m_max_fchange;
	std::string reason;
	bool m_save;
	bool m_revert;
	bool LM;

	double alpha;
	double alphastart;
	double alphamax;
};

/**
 * \ref SingletonFactory that returns pointers to \ref ConvergenceDetector.
 */
typedef SingletonFactory<ConvergenceDetector> ConvergenceDetectorFactory;
