import os
import subprocess as sub
import time
import sys
import traceback

import nibabel as nib
import numpy as np
from matplotlib.backends.windowing import FocusManager
		
def check_dims(data, shape):
	""" 
	Check that the dimensions of the given data are consistent with the
	given shape. Only compares up to the largest common dimension, and
	ignores dimensions of 0 (these have not been initialized yet).
	
	Returns a new shape in which the number of dimensions may be expanded
	to account for the data. All future data must be consistent with this
	expanded set of dimensions 
	"""
	for d in range(min(len(shape), len(data.shape))):
		if shape[d] == 0:
			shape[d] = data.shape[d]
		elif data.shape[d] != shape[d]:
			raise Exception("Inconsistent Dimensions: %i (%i != %i)" % (d, data.shape[d], shape[d]))
		
	return shape

class Data:
	"""
	Represents a 3d or 4d volume, either
	input or output data
	"""
	def __init__(self, role, filename):
		self.filename = filename
		self._load()
		self.shape = self.data.shape 
		self.role = role
		self.ndims = len(self.shape)
		self.alpha = 1.0
		if role == "mask": 
			self.visible = False
			self.cm = "autumn"
		elif role == "data":
			self.visible = True
			self.cm = "gray"
		elif role.lower() == "model prediction":
			self.visible = True
			self.cm = "spectral"
		else:
			self.visible = False
			self.cm = "spectral"
		self.min = self.data.min()
		self.max = self.data.max()
		
	def _load(self):
		print("Loading data: " + self.filename)
		try:
			f = nib.load(self.filename)
			self.data = np.ma.masked_values(f.get_data(), 0)
			self.affine = f.affine
		except:
			print("Failed to load data...")
			traceback.print_exc()
			return
		
	def get_slice(self, focus, d):
		x,y,z,t = focus
		if self.ndims == 4:
			if d == "x": return self.data[x,:,:,t]
			if d == "y": return self.data[:,y,:,t]
			if d == "z": return self.data[:,:,z,t]
		else:
			if d == "x": return self.data[x,:,:]
			if d == "y": return self.data[:,y,:]
			if d == "z": return self.data[:,:,z]
	
	def get_value(self, focus):
		x,y,z,t = focus
		if self.ndims == 4:
			return self.data[x,y,z,t]
		else:
			return self.data[x,y,z]
		
	def get_timeseries(self, focus):
		x,y,z,t = focus
		if self.ndims == 4:
			return self.data[x,y,z,:]
		else:
			return []

class Run:
	"""
	A particular run of fabber, with its output
	directory, logfile and output data
	"""
	def __init__(self, dir, shape=[]):
		self.dir = dir
		self.data = {}
		self.shape = shape
		self.logfile = os.path.join(dir, "logfile")
		self.timestamp = self._get_timestamp()
		self.isquick = self._is_quick_run()
		
	def get_log(self):
		f = open(self.logfile)
		log = "".join(f.readlines())
		f.close()
		return log
				
	def scan_output(self):
		if len(self.data) > 0: return
		try:
			outdir_files = [f for f in os.listdir(self.dir) 
					if os.path.isfile(os.path.join(self.dir, f)) and f.endswith(".nii.gz")]
		except:
			print("Could not read output directory: ", self.dir)
			traceback.print_exc()
			return
		
		for fname in outdir_files:
			try:
				name = fname.split(".")[0]
				role = self._get_role_from_name(name)
				item = Data(role, os.path.join(self.dir, fname))
				self.shape = check_dims(item, self.shape)
				self.data[role] = item
			except:
				print("Failed to load output file: ", fname)
	
	def _get_timestamp(self):
		f = open(self.logfile)
		lines = f.readlines()
		f.close()
		for line in lines:
			l = line.strip()
			if l.lower().startswith("start time:"):
				return l[12:]
		raise Exeception("Could not find timestamp in log file")
		
	def _is_quick_run(self):
		return os.path.isfile(os.path.join(self.dir, "QUICKRUN.txt"))
	
	def _get_role_from_name(self, name):
		""" Get file role (human readable unique description of
		    what it contains) from its filesname"""
		if name.startswith("mean_"):
			return name[5:] + " Mean value"
		elif name.startswith("std_"):
			return name[4:] + " Std. dev."
		elif name.startswith("zstat_"):
			return name[6:] + " Zstat"
		elif name.startswith("modelfit"):
			return "Model prediction"
		elif name.startswith("residuals"):
			return "Model residuals"
		else: return name

class RunFile:
	"""
	A fabber options file
	"""
	def __init__(self, filename=None):
		# Fully qualified path to directory containing the .fab file		
		self.fabdir = None
		
		# File name only of the fab file
		self.fabfile = None
		
		# Fully qualified output dir from the options. If specified as a
		# relative path in the .fab file, will be taken as relative to 
		# the directory containing the fab file (self.fabdir)
		self.outdir = None
		
		# Input data, keyed by role
		self.data = {}
		
		# Previous runs detected based on output dir
		self.runs = []
		
		# Mapping of options as key/value pairs
		self.options = {}
		
		# Lines in the file in order, either a key or a comment string
		# Enables file to be recreated with comments in right place
		self.filelines = []
		
		# Input data shape - must be consistent among all data files
		self.shape = [0,0,0,0]
		
		# Used to record what has changed before updating views
		self.changes = set()

		# List of known views of the model
		self.views = []

		# Contextual - really part of view
		self.focus = [0,0,0,0]
		self.current_data = None
		self.current_run = None

		if filename: 
			self.set_filename(filename)
			self._parse()
		else: self._init_default()

	def changed(self, *stuff):
		""" 
		Find out if a particular item has changed, e.g. if fab.changed("focus")...
		"""
		for s in stuff:
			if s in self.changes:
				return True
	
	def get_option(self, key):
		return self.options[key]
		
	def set_option(self, key, value):
		key, value = key.strip(), value.strip()
		if key == "": return
		if key not in self.options:
			self.filelines.append(key)
		self.options[key] = value
		if key[:4]=="data" or key == "mask":
			self._set_data(key, value)
		elif key == "output":
			self._set_outdir()
		self._change("options")		
		self._update_views()

	def set_filename(self, fname):
		self.fabdir, self.fabfile = os.path.split(fname)
		self.fabdir = os.path.abspath(self.fabdir)
		self.filename = os.path.abspath(fname)
		self._set_outdir()
		self._update_views()
		
	def clear_option(self, key):
		key = key.strip()
		if key == "": return
		if key in self.options:
			del self.options[key]
			self.filelines.remove(key)
		self._change("options")		
		self._update_views()
				
	def get_methods(self):
		""" Get known methods. Uses the fabber executable to find out"""
		stdout = self._run_help("--listmethods")
		methods = [line.strip() for line in stdout.split()]
		return methods

	def get_models(self):
		""" Get known inference methods. Uses the fabber executable to find out"""
		stdout = self._run_help("--listmodels")
		models = [line.strip() for line in stdout.split()]
		return models

	def get_options(self, method=None, model=None):
		"""
		Get general options, or options for a method/model
		"""
		args = []
		if method: args.append("--method=%s" % method)
		if model: args.append("--model=%s" % model)
		stdout = self._run_help("--help", *args)
		lines = stdout.split("\n")
		opts = []
		desc = ""
		descnext = False
		for line in lines:
			if line.startswith("--"):
				line2 = line[2:].replace("[", ",").replace("]", "")
				opt = [o.strip() for o in line2.split(",")]
				if opt[3].startswith("DEFAULT="): opt[3] = opt[3].split("=", 1)[1]
				else: opt[3] = ""
				if opt[0] in ["data", "data<n>", "mask", "help", "model", "method", "listmodels", "listmethods"]:
					continue
				if len(opt) == 4:
					opts.append(opt)
					descnext = True
				else: 
					print("Couldn't _parse option: " + line)
			elif descnext:
				opts[-1].append(line.strip())
				descnext = False
			elif not line.startswith("Usage") and not line.startswith("Options") and not line.startswith("Description"):
				desc += line.strip()

		return opts, desc	
			
	def set_exec(self, ex):
		""" 
		Set the path to the fabber executable 
		"""
		self._change("exec")
		self.set_option("fabber", ex)

	def set_visibility(self, item, alpha=None, visible=None, cm=None):
		""" 
		Set visibility of a data item 
		"""
		if alpha is not None: item.alpha = alpha
		if visible is not None: item.visible = visible
		if cm is not None: item.cm = cm
		self._change("data")
		self._update_views()

	def select(self, key):
		""" 
		Select a data item 
		"""
		if key in self.data:
			self.current_data = self.data[key]
		elif key in self.current_run.data:
			self.current_data = self.current_run.data[key]
		else:
			raise RuntimeWarning("Data not found in rundata: %s" % key)
		self._change("selection")		
		self._update_views()

	def update_focus(self, xp=None, yp=None, zp=None, tp=None):
		""" 
		Update the point of focus 
		"""
		if xp: self.focus[0] = xp
		if yp: self.focus[1] = yp
		if zp: self.focus[2] = zp
		if tp: self.focus[3] = tp
		self._change("focus")		
		self._update_views()

	def add_view(self, view):
		"""
		Add a new view and update it
		"""
		self.views.append(view)
		self._change()
		view.update(self)
		self.changes.clear()
		
	def get_runs(self):
		"""
		Get a list of known previous runs. The first in the list is the most recent
		"""
		self.runs = []
		if not self.outdir: return
		dir = self.outdir
		while True:
			try:
				run = Run(dir, shape=self.shape)
			except:
				print("Could not read logfile in: " + dir)
				traceback.print_exc()
				break
			self.runs.append(run)
			dir += "+"
		self.runs.reverse()
		
	def select_run(self, dir=None):
		"""
		Set a particular run as the current one
		"""
		for run in self.runs:
			if run.dir == dir:
				self.current_run = run
				run.scan_output()
				break
		self._change("selection")
		self._update_views()

	def save(self, quick=False):
		"""
		Save .fab file. If quick, save a special version for
		a quick 1-voxel run
		"""
		if not self.filename:
			raise "No file specified for save"
		if quick:
			fname = self.filename + "_temp"
		else:
			fname = self.filename
		print(fname)
		fab = open(fname, "w")
		self._dump(fab, quick)
		self._dump(sys.stdout, quick)
		fab.close()
		return fname

	def run(self, quick=False):
		"""
		Run fabber. If quick do a quick 1-voxel run at the current
		focus
		
		FIXME identifying the actual output directory is not good
		"""
		fname = self.save(quick)
		cmd = [self.options["fabber"], "-f", fname]
		print(cmd)
		print(self.fabdir)
		proc = sub.Popen(cmd, stdout=sub.PIPE, cwd=self.fabdir)
		while 1:
			for line in proc.stdout.readlines():
				print(line)
			if proc.poll() is not None: break
			time.sleep(1)

		self.get_runs()
		if self.runs: 
			self.select_run(self.runs[0].dir)
			self.log = self.runs[0].get_log()
			
		# Was this a temporary quick run?
		if fname != self.filename:
			os.remove(fname)
			self.runs[0].isquick = True
			
		self._change("data")		
		self._update_views()

	def _change(self, *stuff):
		""" 
		Record a _change, e.g. self._change("data", "options"). Views will
		use this info to decide what they need to update
		"""		
		if len(stuff) == 0:
			self.changes = {"exec", "data", "options", "selection", "focus"}
		else:
			for s in stuff:
				self.changes.add(s)

	def _update_views(self):
		"""
		Update all views
		"""
		for view in self.views:
			try:
				view.update(self)
			except:
				print("Failed to update view")
				traceback.print_exc()
		self.changes.clear()


	def _parse(self):		
		fab = open(self.filename, "r")		
		for line in fab.readlines():
			line = line.strip()
			if len(line) == 0 or line[0] == "#":
				self.filelines.append(line)
			else:
				kv = line.split("=", 1)
				key = kv[0].strip()
				if len(kv) > 1: value = kv[1].strip()
				else: value = ""
				self.options[key] = value
				self.filelines.append(key)
				if key[:4]=="data" or key == "mask":
					self._set_data(key, value)
				
		fab.close()
		self._find_fabber()

		self._set_outdir()
		self.get_runs()
		if self.runs: self.select_run(self.runs[0].dir)
		if self.data.has_key("data"): self.select("data")
		elif self.data.has_key("data1"): self.select("data1")
		
		self._change()
		self._update_views()
	
	def _set_data(self, key, value):
		if not os.path.isabs(value):
			value = os.path.join(self.fabdir, value)
		item = Data(key, value)
		check_dims(item, self.shape)
		self.data[key] = item
		for idx, f in enumerate(self.focus):
			if f == 0: 
				self.focus[idx] = self.shape[idx] / 2
				self._change("focus")
			
	def _set_outdir(self):
		if self.options.has_key("output"):
			self.outdir = self.options["output"]
			if not os.path.isabs(self.outdir):
				self.outdir = os.path.join(self.fabdir, self.outdir)
		else:
			print("WARNING: no output dir")

	def _write_temp_mask(self):
		fname = os.path.join(self.fabdir, "fabber_mask_temp.nii.gz")
		data = np.zeros(self.shape[:3])
		data[self.focus[0], self.focus[1], self.focus[2]] = 1
		affine = self.data["data"].affine # FIXME
		img = nib.Nifti1Image(data, affine)
		nib.save(img, fname)
		return fname

	def _dump(self, dev, quick=False):
		for line in self.filelines:
			if len(line) == 0 or line[0] == "#":
				dev.write(line)
			else:
				if self.options[line] == "":
					dev.write("%s" % line)
				elif quick and line=="mask":
					pass
				else:
					dev.write("%s=%s" % (line, self.options[line]))
			dev.write("\n")
		if quick:
			temp_mask = self._write_temp_mask()
			dev.write("mask=%s" % (temp_mask))
			
	def _write_quickrun_file(self, dir):
		"""
		Write a little file which identifies a particular
		run as a 1-voxel test run
		"""
		f = open(os.path.join(dir, "QUICKRUN.txt"), "wc")
		f.write("This data is from a 1-voxel test run\n")
		f.close()
		
	def _run_help(self, *opts):
		"""
		Run Fabber synchronously, this is presumed to be a quick
		run in order to get help options.

		Does not throw on failure because user may need to set
		the location of the executable first
		"""
		cmd = [self.options["fabber"]] + list(opts)
		print cmd
		try:
			p = sub.Popen(cmd, stdout=sub.PIPE, cwd=self.fabdir)
			(stdout, stderr) = p.communicate()
			status = p.wait()
			if status == 0:
				return stdout
			else:
				print("Failed to run fabber: " + str(cmd))
				return ""
		except:
			print("Failed to run fabber: " + str(cmd))
			traceback.print_exc()
			return ""
	
	def _find_fabber(self):
		""" 
		Find a fabber executable
		
		FIXME need to do better than this!
		"""
		if not self.options.has_key("fabber"):
			self.options["fabber"] = "/home/martinc/dev/fabber_core/fabber"
		
	def _init_default(self):
		self.set_filename("newfile.fab")
		
		self.options = {"method" : "vb", "model" : "poly", "degree" : "2", "output" : "out",
					    "overwrite" : "", "save-model-fit" : "", "save-residuals" : "", }
		self.filelines = ["# Generated by fabber_gui.py", "", 
						  "output", "overwrite", "save-model-fit", "save-residuals", ""
						  "method", "model", "degree"]
		self._find_fabber()
		self._set_outdir()
		self._change()
		self._update_views()

		
