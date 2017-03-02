/*  dist_gamma.cc - Gamma distribution class

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#include "dist_gamma.h"

#include <math.h>
#include <ostream>

GammaDist::GammaDist(EasyLog *log)
    : Loggable(log)
    , b(0)
    , c(0)
{
}

double GammaDist::CalcMean() const
{
    return b * c;
}

double GammaDist::CalcVariance() const
{
    return b * b * c;
}

void GammaDist::SetMeanVariance(double m, double v)
{
    b = v / m;
    c = m / b;
}

void GammaDist::Dump(std::ostream &os) const
{
    os << "Noise stdev == " << 1.0 / sqrt(b * c) << " (b==" << b << ", c==" << c << ")" << std::endl;
}
