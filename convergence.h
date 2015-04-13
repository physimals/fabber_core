/*  convergence.h - Convergence detectors for FABBER

    Adrian Groves and Michael Chappell, FMRIB Image Analysis Group

    Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

class ConvergenceDetector {
 public:
  virtual bool Test(double F) = 0;  // to be called BEFORE each iteration
  virtual ~ConvergenceDetector() { return; }
  virtual void DumpTo(ostream& out, const string indent = "") const = 0;
  virtual void Reset(double F=-99e99) = 0; 
  virtual bool UseF() const = 0;
  virtual bool NeedSave() = 0;
  virtual bool NeedRevert() = 0;
  virtual float LMalpha() = 0;
};

class CountingConvergenceDetector : public ConvergenceDetector {
  public:
    virtual bool Test(double) 
        { return ++its >= max; }
    CountingConvergenceDetector(int maxIts) 
        : max(maxIts) { assert(max>0); Reset(); }
      virtual void DumpTo(ostream& out, const string indent = "") const
        { out << indent << "Starting iteration " << its+1 << " of " << max << endl << endl; }
    virtual void Reset(double F=-99e99) 
        { its = 0; }
    virtual bool UseF() const { return false; }        
    virtual bool NeedSave() { return false; }
    virtual bool NeedRevert() {return false; }
  virtual float LMalpha() {return 0.0;}
  
  private:
    int its;
    const int max;
};

class FchangeConvergenceDetector : public ConvergenceDetector {
  public:
    virtual bool Test(double F);
    FchangeConvergenceDetector(int maxIts, double Fchange)
        : max(maxIts), chg(Fchange)  { assert(max>0); assert(chg>0); Reset(); }
      virtual void DumpTo(ostream& out, const string indent = "") const;
    virtual void Reset(double F=-99e99)
        { its = 0; prev = F; }
    virtual bool UseF() const { return true; }      
    virtual bool NeedSave() { return false; }
    virtual bool NeedRevert() {return false; }
  virtual float LMalpha() {return 0.0;}

  private:
    int its;
    const int max;
    double prev;
    const double chg;
};

inline bool FchangeConvergenceDetector::Test(double F)
{   
//    if (F == 1234.5678)
//        throw logic_error("FchangeConvergenceDetector needs F, but it seems it isn't being calculated!  Internal bug... should have needF = true");
        // F could actually be 1234.5678, but what are the chances of that?        
    double diff = F - prev;
    prev = F; 
    ++its;
    
    diff = diff>0 ? diff : -diff;
    
    return (its >= max) || (diff < chg);
}

inline void FchangeConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
    out << indent << "Iteration " << its << " of at most " << max << endl;
    out << indent << "Previous Free Energy == " << prev << endl;
} 

class FreduceConvergenceDetector : public ConvergenceDetector {
 public:
  virtual bool Test(double F);
  FreduceConvergenceDetector(int maxIts, double Fchange)
    : max(maxIts), chg(Fchange) {assert(max>0); assert(chg>0); Reset(); }
  virtual void DumpTo(ostream& out, const string indent = "") const;
  virtual void Reset(double F=-99e99)
  { its = 0; prev = F; revert = false; }
  virtual bool UseF() const {return true;}
  virtual bool NeedSave() { return true; } //this simple convergence detector always saves every interation
  virtual bool NeedRevert();
  //virtual bool DoNoiseUpdate() {return true;}
  virtual float LMalpha() {return 0.0;}

 private:
  int its;
  const int max;
  double prev;
  const double chg;
  string reason;
  bool revert;// determines whether we should revert or not if asked
};

inline bool FreduceConvergenceDetector::Test(double F)
{
  double diff = F - prev;
  prev = F;
  ++its;
  reason = "blank";
  
  double absdiff = diff; //look at magnitude of diff in absdiff
  if(diff<0) {absdiff = -diff;}

  if (its >= max) {reason = "Max. Iterations Reached";}
  else if (diff < 0) {
    reason = "F reduced - halt here"; 
    revert=true;
  }
  else if (absdiff < chg) {reason = "F Settled";}
  
  return (its >= max) || (absdiff < chg) || (diff < 0);
}

inline bool FreduceConvergenceDetector::NeedRevert()
{
  return revert;
}

inline void FreduceConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
  out << indent << "Iteration " << its << " of at most " << max << " : " << reason << endl;
  out << indent << "Previous Free Energy == " << prev << endl;
}


class TrialModeConvergenceDetector : public ConvergenceDetector {
 public:
  virtual bool Test(double F);
 TrialModeConvergenceDetector(int maxIts, int maxTrials, double Fchange)
   : max(maxIts), maxT(maxTrials), chg(Fchange) {assert(max>0); assert(maxT>0); assert(chg>0); Reset(); }
  virtual void DumpTo(ostream& out, const string indent = "") const;
  virtual void Reset(double F=-99e99)
  { its = 0; prev = F; trials = 0; save = true; revert = false; }
  virtual bool UseF() const {return true;}
  virtual bool NeedSave();
  virtual bool NeedRevert();
  //virtual bool DoNoiseUpdate() {return true;} 
  virtual float LMalpha() {return 0.0;}
 private:
  int its;
  int trials;
  const int max;
  const int maxT;
  double prev;
  const double chg;
  string reason;
  bool save;
  bool revert;
  bool trialmode;
};

inline bool TrialModeConvergenceDetector::NeedSave() { return save; }
inline bool TrialModeConvergenceDetector::NeedRevert() { return revert; }

inline bool TrialModeConvergenceDetector::Test(double F)
{
  double diff = F - prev;

  double absdiff = diff; //look at magnitude of diff in absdiff
  if(diff<0) {absdiff = -diff;}
  
  //if we are not in trial mode
  if (!trialmode) {
    // if F has reduced then we will enter trial mode
    if (diff < 0) {
      ++trials;
      trialmode = true;
      revert = true; save = false;
      return false;
    }
    //otherwise if we have converged stop
    else if (absdiff < chg) {
      reason = "F converged";
      return true;
    }
      //otherwise if we have reached max iterations stop
    else if (its >= max) {
      reason = "Max iterations reached";
      return true;
    }
    // otherwise carry on to next iteration
    else {
      prev = F;
      ++its;
      return false;
    }
  }
    // if we are in trial mode
  else {
    // if F has improved over our previous best then resume iterations
    if ( diff > 0 ) {
      trialmode = false;
      save = true; revert = false;
      ++its;
      trials = 0; //reset number of trials for future use
      prev = F;
      return false;
    }
    //if we have exceeded max trials then stop and output previous best result
    else if (trials >= maxT ) {
      reason = "Reached max trials";
      return true;
    }
    // otherwise continue in trial mode for time being
    else {
      ++trials;
      return false;
    }
  }
}

inline void TrialModeConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
  out << indent << "Iteration " << its << " of at most " << max << " : " << reason << endl;
  out << indent << "Previous Free Energy == " << prev << endl;
}

class LMConvergenceDetector : public ConvergenceDetector {
 public:
  virtual bool Test(double F);
 LMConvergenceDetector(int maxIts, double Fchange)
   : max(maxIts), chg(Fchange) {assert(max>0); assert(chg>0); Reset(); }
  virtual void DumpTo(ostream& out, const string indent = "") const;
  virtual void Reset(double F=-99e99)
  { its = 0; prev = F; save = true; revert = false; alphastart=1e-6; alpha=0.0; alphamax=1e6; LM=false; }
  virtual bool UseF() const {return true;}
  virtual bool NeedSave();
  virtual bool NeedRevert();
  //virtual bool DoNoiseUpdate();
  virtual float LMalpha();
 private:
  int its;
  const int max;
  double prev;
  const double chg;
  string reason;
  bool save;
  bool revert;
  bool LM;

  double alpha;
  double alphastart;
  double alphamax;
};

inline bool LMConvergenceDetector::NeedSave() { return save; }
inline bool LMConvergenceDetector::NeedRevert() { return revert; }
//inline bool LMConvergenceDetector::DoNoiseUpdate() {return !LM; } //dont update the noise when we are in LM mode
inline float LMConvergenceDetector::LMalpha() { return alpha; }

inline bool LMConvergenceDetector::Test(double F)
{
  double diff = F - prev;

  // cout << "----" << endl;
  // cout << its << endl;
  // cout << "F:" << F << endl;
  // cout << "prev:" << prev << endl;
  // cout << diff << endl;

  double absdiff = diff; //look at magnitude of diff in absdiff
  if(diff<0) {absdiff = -diff;}
  
  //if we are not in LM mode
  if (!LM) {
    // if F has reduced then we go into LM mode
    if (diff < 0) {
      LM = true;
      revert=true; //revert to the previous solution and try again with LM adjustment
      alpha=alphastart;
      //cout << "Entering LM" << endl;
      return false;
    }
    //otherwise if we have converged stop
    else if (absdiff < chg) {
      reason = "F converged";
      revert=false;
      return true;
    }
      //otherwise if we have reached max iterations stop
    else if (its >= max) {
      reason = "Max iterations reached";
      revert=false;
      return true;
    }
    // otherwise carry on to next iteration
    else {
      prev = F;
      ++its;
      return false;
      revert=false;
    }
  }
    // if we are in LM mode (NB we dont increase iterations if we are increasing the LM alpha, onyl when we make a sucessful step)
  else {
    // if F has improved over our previous best then reduce the alpha and continue from current estimate

    if ( diff > 0 ) {
      if (alpha == alphastart) { //leave LM mode if alpha returns to inital value
	LM = false;
      }
      else {
	alpha /= 10;
	LM = true;
      }
      revert=false;
      prev = F;
      cout << "Reducing LM" << endl;
      ++its; // if F has improved then we take this as the new 'step' and move onto the next iteration
      return false;
    }
    //if alpha gets too large then we cannot achieve any gain here so stop and revert to previous best solution
    else if (alpha >= alphamax ) {
      reason = "Reached max alpha";
      revert=true;
      //cout << "LM maxed out" << endl;
      return true;
    }
    //otherwise if we have reached max iterations stop
    else if (its >= max) {
      reason = "Max iterations reached";
      revert=false;
      return true;
    }
    // otherwise continue in LM mode for time being, try increasing alpha
    else {
      alpha *= 10;
      revert=true; //revert to the previous solution and try again with new LM adjustment
      //cout << "Increasing LM" << endl;
      return false;
    }
  }
}

inline void LMConvergenceDetector::DumpTo(ostream& out, const string indent) const
{
  out << indent << "Iteration " << its << " of at most " << max << " : " << reason << endl;
  out << indent << "Previous Free Energy == " << prev << endl;
}
