#!/bin/env python
"""
Simple script to generate a test image for a model for specified parameter values
"""

usage = """
Usage: 

generate_test_image.py --param-<param_name>=<param_value>,<param_value>,...
                       --param-<param_name>=<param_value>,<param_value>,...  
                       ... more parameters 
                       --nt=<number of time values to evaluate> 
                       [--outfile=output file (default test_image.nii)]
                       [--model=<model name> ... any additional Fabber model options]

A value must be specified for all model parameters that are required for the 
specified model options. No more than 3 parameters can have multiple values
in order to create a 3D Nifti volume

Example: Create a test dataset for the simple polynomial model

generate_test_image.py --model=poly --degree=2 --param-c0=10,20,30 --param-c1=0.5,1,2 --param-c2=5 --nt=10
"""

import os, sys
import traceback

import numpy as np
import nibabel as nib

sys.path.insert(0, os.environ["FSLDIR"] + "/lib/python")
from fabber import FabberRunData, generate_test_data

try:
    rundata = FabberRunData()
    params = {}
    model_libs = []
    outfile="test_image.nii"
    nt = 10
    for arg in sys.argv[1:]:
        if arg.startswith("--loadmodels"):
            model_libs.append(arg.split("=")[1])
        if arg.startswith("--outfile"):
            outfile=arg.split("=")[1]
        elif arg.startswith("--nt"):
            nt = int(arg.split("=")[1])
        elif arg.startswith("--param-"):
            key, value = arg[8:].split("=")
            params[key] = [float(v) for v in value.split(",")]
        elif arg.startswith("--"):
            kv = arg[2:].split("=")
            if len(kv) == 1: rundata[kv[0]] = ""
            else: rundata[kv[0]] = kv[1]
        else:
            raise RuntimeError("Invalid argument: %s" % arg)
except Exception, e:
    print("Error: %s" % str(e))
    print(usage)
    sys.exit(1)

try:
    data = generate_test_data(rundata, params, nt=nt, model_libs=model_libs)
    nii = nib.Nifti1Image(data, np.identity(4))
    nii.to_filename(outfile)
except:
    traceback.print_exc(limit=0)