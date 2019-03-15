#!/bin/env python

import sys
import traceback

from fabber import self_test, FabberException

save = "--save" in sys.argv
try:
    rundata= {
        "model" : "exp",
        "num-exps" : 2,
        "dt" : 0.02,
        "max-iterations" : 50,
    }
    params = {
        "amp1" : [1, 0.5],    # Amplitude first exponential
        "amp2" : [0.5],      # Amplitude second exponential
        "r1" : [1.0, 0.8],           # Decay rate of first exponential
        "r2" : [6.0],           # Decay rate of second exponential
    }
    test_config = {
        "nt" : 100,           # Number of time points
        "noise" : 0.1,        # Amplitude of Gaussian noise to add to simulated data
        "patchsize" : 20,     # Each patch is 20 voxels along each dimension
    }
    result, log = self_test("exp", rundata, params, save_input=save, save_output=save, invert=True, **test_config)
except FabberException, e:
    print e.log
    traceback.print_exc()
except:
    traceback.print_exc()

