/*  tools.h - miscellaneous useful functions & classes

    Adrian Groves, FMRIB Image Analysis Group

    Copyright (C) 2007 University of Oxford  */

/*  CCOPYRIGHT */

#include <newmatio.h>
using namespace NEWMAT;
#include <math.h>
#include <assert.h>

class GenericFunction1D
{
public:
    virtual double Calculate(double x) const = 0;
    double operator() (double x) const { return Calculate(x); }
    // Perhaps add more functions.. e.g. initialGuess, domain of validity, etc.

    // This is useful if your function is very slow to calculate, but you have
    // some cached partial calculations available.  If you have a suitable
    // cached value, store it into guess and return true.  Otherwise return
    // false (and leave guess unchanged).
    virtual bool PickFasterGuess(double* guess, double lower, double upper, bool allowEndpoints = false) const
    { return false; }

    virtual ~GenericFunction1D() { return; }
private:
    // Function's constant data should go here
};

#define REALMAX (1.7976931348623158e+308)

class Guesstimator
{
 public:
  virtual double GetGuess(double lower, double upper, 
			  double atLower, double atUpper) = 0;
  virtual ~Guesstimator() { return; }
};

class BisectionGuesstimator : public Guesstimator
{
 public:
  virtual double GetGuess(double lower, double upper, double, double)
  { assert(lower<upper); return (lower+upper)/2; }
};

class LogBisectionGuesstimator : public Guesstimator
{
 public:
  virtual double GetGuess(double lower, double upper, double, double)
  { assert(lower>0 && upper>lower); 
    double guess = sqrt(lower*upper); 
    if (lower >= guess || guess >= upper)
      {
	cout << "Uh-oh... lower = " << lower
	     << ", guess = " << guess
	     << ", upper = " << upper << endl;
      }
    return guess;
  }
};

class RiddlersGuesstimator : public Guesstimator
{
 public: 
  virtual double GetGuess(double lower, double upper, double atLower, double atUpper);
  RiddlersGuesstimator() { halfDone = false; }
 private:
  bool halfDone; // waiting for f(x3) result?
  double x1, x2, fx1, fx2; // save from phase 1 for phase 2; only valid when halfDone is true.
};

class LogRiddlersGuesstimator : public RiddlersGuesstimator
{
 public:
  virtual double GetGuess(double lower, double upper, double atLower, double atUpper)
    { 
      // Just a log transformation
      return exp(RiddlersGuesstimator::GetGuess(log(lower),log(upper),atLower,atUpper)); 
    }
};

class InterpGuesstimator : public Guesstimator
{
 public:
  virtual double GetGuess(double lower, double upper, 
			  double atLower, double atUpper)
  { 
    double guess = upper - atUpper*(upper-lower)/(atUpper-atLower);

    if (lower >= guess || guess >= upper)
      {
	cout << "Uh-oh... lower = " << lower
	     << ", guess = " << guess
	     << ", upper = " << upper
	     << ", atLower = " << atLower
	     << ", atUpper = " << atUpper << endl;
      }
    return guess;
  }
};

// Note that you have to specify one or more tolerances to get any sensible
// results.  Also note that the ratio tolerances current assume that the
// X or Y value always positive -- otherwise it'll stop too early!
class ZeroFinder
{
public:
    ZeroFinder(const GenericFunction1D& f)
      : fcn(f), searchMin(-REALMAX), searchMax(REALMAX), 
      searchGuess(0), searchScale(REALMAX), searchScaleGrowth(2),
      maxEvaluations(1000000), tolX(REALMAX), tolY(REALMAX), 
      ratioTolX(REALMAX), ratioTolY(REALMAX), 
      guesstimator(new BisectionGuesstimator()), verbosity(2)
      { return; }
    virtual double FindZero() const = 0;
    operator double () const { return FindZero(); }
    virtual ~ZeroFinder() { return; }

    // Set optional parameters through these
    ZeroFinder& InitialGuess(double guess) 
        { searchGuess = guess; return *this; }
    ZeroFinder& SearchMin(double min) 
        { searchMin = min; return *this; }
    ZeroFinder& SearchMax(double max) 
        { searchMax = max; return *this; }
    ZeroFinder& InitialScale(double scale) 
        { searchScale = scale; return *this; }
    ZeroFinder& ScaleGrowth(double growth)
        { assert(growth>1); searchScaleGrowth = growth; return *this; }
    ZeroFinder& MaxEvaluations(int evals)
        { assert(evals>1); maxEvaluations = evals; return *this; }
    ZeroFinder& TolX(double tx)
        { assert(tx>0); tolX = tx; return *this; }
    ZeroFinder& TolY(double ty)
        { assert(ty>0); tolY = ty; return *this; }
    //    ZeroFinder& RatioTolY(double rty) // utterly pointless -- looking for a sign change!
    //      { assert(rty>1); ratioTolY = rty; return *this; }
    ZeroFinder& RatioTolX(double rtx)
      { assert(rtx>1); ratioTolX = rtx; return *this; }
    ZeroFinder& SetGuesstimator(Guesstimator* g)
      { delete guesstimator; guesstimator = g; return *this; }
    ZeroFinder& Verbosity(int v)
      { verbosity = v; return *this; }

protected:
    const GenericFunction1D& fcn;

    // Optional parameters:
    double searchMin;
    double searchMax;
    double searchGuess;
    double searchScale;
    double searchScaleGrowth;
    int maxEvaluations;
    double tolX;
    double tolY;
    double ratioTolX;
    double ratioTolY;
    Guesstimator* guesstimator;
    int verbosity;
};

class DescendingZeroFinder : public ZeroFinder
{
public:
    DescendingZeroFinder(const GenericFunction1D& f) : ZeroFinder(f) { return; }
    virtual double FindZero() const;
};
