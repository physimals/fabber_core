#!/bin/env python
"""
Simple script to evaluate a model for specified parameter values
"""

usage = """
Usage: 

evaluate.py --param-<param_name>=<param_value> 
            --param-<param_name>=<param_value>  ... 
            --nt=<number of time values to evaluate> 
            [--model=<model name> ... any additional Fabber model options]

A value must be specified for all model parameters that are required for the 
specified model options

Example: Evaluate the simple polynomial model

evaluate.py --model=poly --degree=2 --param-c0=0.5 --param-c1=0.2 --param-c2=0.1 --nt=10
"""

import os, sys
import traceback

sys.path.insert(0, os.environ["FSLDIR"] + "/lib/python")
from fabber import FabberLib, FabberRunData, FabberException

try:
    rundata = FabberRunData()
    params = {}
    model_libs = []
    nt = 10
    for arg in sys.argv[1:]:
        if arg.startswith("--loadmodels"):
            model_libs.append(arg.split("=")[1])
        elif arg.startswith("--nt"):
            nt = int(arg.split("=")[1])
        elif arg.startswith("--param-"):
            key, value = arg[8:].split("=")
            params[key] = float(value)
        elif arg.startswith("--"):
            kv = arg[2:].split("=")
            if len(kv) == 1: rundata[kv[0]] = ""
            else: rundata[kv[0]] = kv[1]
        else:
            raise RuntimeError("Invalid argument: %s" % arg)

    fab = FabberLib(model_libs=model_libs)

    ret = fab.model_evaluate(rundata, params, nt)
    for v in ret: print v
except Exception, e:
    print("Error: %s" % str(e))
    print(usage)
    sys.exit(1)
