/*  dist_gamma.h - Gamma distribution class/structure

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include "easylog.h"

#include <math.h>
#include <ostream>

class GammaDist : public Loggable {
public:
    explicit GammaDist(EasyLog* log = 0);

    double CalcMean() const;
    double CalcVariance() const;
    //  double CalcLogMoment() { return digamma(b) + log(c); }  // where can I get digamma from?
    void SetMeanVariance(double m, double v);
    void Dump(std::ostream& os) const;

    double b;
    double c;
};

inline std::ostream& operator<<(std::ostream& out, const GammaDist& dist)
{
    dist.Dump(out);
    return out;
}
