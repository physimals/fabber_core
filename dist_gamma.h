/*  dist_gamma.h - Gamma distribution class/structure

 Adrian Groves, FMRIB Image Analysis Group

 Copyright (C) 2007-2008 University of Oxford  */

/*  CCOPYRIGHT */

#pragma once
#include <math.h>
#include "easylog.h"

class GammaDist
{
public:
	GammaDist()
	{
		b = c = 0.0;
	}
	double b;
	double c;
	double CalcMean() const
	{
		return b * c;
	}
	double CalcVariance() const
	{
		return b * b * c;
	}
	//  double CalcLogMoment() { return digamma(b) + log(c); }  // where can I get digamma from?
	void SetMeanVariance(double m, double v)
	{
		b = v / m;
		c = m / b;
	}
	void Dump(const string indent = "") const
	{
		LOG << indent << "Noise stdev == " << 1.0 / sqrt(b * c) << " (b==" << b << ", c==" << c << ")" << endl;
	}
};

