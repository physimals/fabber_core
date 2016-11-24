#!/bin/python

import sys, os
import nibabel as nib
import numpy as np
import math


def check_file(f, fp1, fp2):
	f1 = nib.load(fp1)
	f2 = nib.load(fp2)
	da1 = f1.get_data()
	da2 = f2.get_data()
	if len(da1.shape) != len(da2.shape):
		print ("%s: Inconsistent number of dimensions: %i!=%i" % (f, len(da1.shape), len(da2.shape)))
		return 1
			
	for n in range(len(da1.shape)):
		if da1.shape[n] != da2.shape[n]:
			print("%s: Inconsistent dimension %i: %i!=%i" % (f, n, da1.shape[n], da2.shape[n]))		
	 	 	return 1

	v1 = da1.flatten()
	v2 = da2.flatten()
	for v in range(len(v1)):
		if abs(v1[v]-v2[v]) > e:
			print("%s: Inconsistent value %i: %f!=%f (within %f)" % (f, v, v1[v], v2[v], e))		
			return 1

	if not quiet: print("%s: OK" % f)
	return 0

def check_dir(d1, d2, recurse):
	ret = 0
	if not quiet: print("Checking %s vs %s\n" % (d1, d2))

	ls1 = set(os.listdir(d1))
	ls2 = set(os.listdir(d2))

	ls = ls1.union(ls2)

	for f in ls:
		diff = 0
		fp1 = os.path.join(d1, f)
		fp2 = os.path.join(d2, f)

		if f not in ls1:
			if not quiet: print("%s: IGNORED (in dir 2 but not 1)" % f)
		elif f not in ls2:
			if not quiet: print("%s: IGNORED (in dir 1 but not 2)" % f)
		elif os.path.isdir(fp1) and os.path.isdir(fp2) and recurse:
			diff = check_dir(fp1, fp2, recurse)
		elif os.path.isfile(fp1) and os.path.isfile(fp2) and f.endswith(".nii.gz"):
			diff = check_file(f, fp1, fp2)
		elif not quiet:
			print("%s: IGNORED (not NIFTI)" % f)
		if diff: ret = 1
	if not quiet: print("\nDIR OK\n")
	return ret
	
e = 0.0001
quiet = 0
if len(sys.argv) in (3, 4):
	d1, d2 = sys.argv[1:3]
	if len(sys.argv) > 3:
		e = float(sys.argv[3])
else:
	print("Usage: diffout.py <dir1> <dir2> [<epsilon>]")
	sys.exit(1)

diff = check_dir(d1, d2, True)
if not quiet:
	if diff: print ("DIFFERENCES FOUND")
	else: print ("NO DIFFERENCES FOUND")
sys.exit(diff)











