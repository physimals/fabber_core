import os
import sys
import warnings
import datetime, time
import collections
import glob
import subprocess as sub
import traceback
import distutils.spawn

from ctypes import *
import numpy as np
import numpy.ctypeslib as npct

import nibabel as nib

if sys.platform.startswith("win"):
    _lib_format = "bin\\%s.dll"
    _bin_format = "bin\\%s.exe"
elif sys.platform.startswith("darwin"):
    _lib_format = "lib/lib%s.dylib"
    _bin_format = "bin/%s"
else:
    _lib_format = "lib/lib%s.so"
    _bin_format = "bin/%s"

def percent_progress(voxel, nvoxels):
    complete = 100*voxel/nvoxels
    sys.stdout.write("\b\b\b\b%3i%%" % complete)
    sys.stdout.flush()
    
def self_test(model, rundata, param_testvalues, save_input=False, save_output=False, invert=True, outfile_format="test_data_%s", **kwargs):
    print("Running self test for model %s" % model)
    rundata["model"] = model
    data, roidata = generate_test_data(rundata, param_testvalues, param_rois=True, auto_load_models=True, **kwargs)
    
    if save_input:
        outfile=outfile_format % model
        print("Saving test data to Nifti files: %s" % outfile)
        data_nii = nib.Nifti1Image(data, np.identity(4))
        data_nii.to_filename(outfile)
        for param in param_testvalues:
            roi_nii = nib.Nifti1Image(roidata[param], np.identity(4))
            roi_nii.to_filename(outfile + "_roi_%s" % param)
    
    log = None
    if invert:
        sys.stdout.write("Inverting test data - running Fabber:  0%%")
        sys.stdout.flush()
        fab = FabberLib(auto_load_models=True)
        if "method" not in rundata: rundata["method"] = "vb"
        if "noise" not in rundata: rundata["noise"] = "white"
        rundata["save-mean"] = ""
        rundata["save-model-fit"] = ""
        rundata["allow-bad-voxels"] = ""
        run = fab.run_with_data(rundata, {"data" : data}, progress_cb=percent_progress)
        print("\n")
        log = run.log
        if save_output:
            data_nii = nib.Nifti1Image(run.data["modelfit"], np.identity(4))
            data_nii.to_filename(outfile + "_modelfit")
        for param, values in param_testvalues.items():
            mean = run.data["mean_%s" % param]
            if save_output:
                data_nii = nib.Nifti1Image(mean, np.identity(4))
                data_nii.to_filename(outfile + "_mean_%s" % param)
            roi = roidata[param]
            print("Parameter: %s" % param)
            for idx, val in enumerate(values):
                out = np.mean(mean[roi==idx+1])
                print("Input %f -> %f Output" % (val, out))
        sys.stdout.flush()
    return log

def generate_test_data(rundata, param_testvalues, nt=10, patchsize=10, 
                       noise=None, patch_rois=False, param_rois=False, **kwargs):
    """ 
    Generate a test Nifti image based on model evaluations

    Returns the image itself - this can be saved to a file using to_filename
    """
    dim_params = []
    dim_values = []
    dim_sizes = []
    fixed_params = {}
    for param, values in param_testvalues.items():
        try:
            val = float(values)
            fixed_params[param] = val
        except:
            if len(values) == 1: fixed_params[param] = values[0]
            else:
                dim_params.append(param)
                dim_values.append(values)
                dim_sizes.append(len(values))

    if len(dim_sizes) > 3: 
        raise RuntimeError("Test image can only have up to 3 dimensions, you supplied %i varying parameters" % len(dim_sizes))
    else:
        for d in range(len(dim_sizes), 3):
            dim_params.append(None)
            dim_values.append([])
            dim_sizes.append(1)

    shape = [d * patchsize for d in dim_sizes]
    data = np.zeros(shape + [nt,])
    if patch_rois: patch_roi_data = np.zeros(shape)
    if param_rois:
        param_roi_data = {}
        for param in dim_params:
            if param is not None: 
                param_roi_data[param] = np.zeros(shape)
    fab = FabberLib(rundata=rundata, **kwargs)

    # I bet there's a neater way to do this!
    patch_label = 1
    for x in range(dim_sizes[0]):
        for y in range(dim_sizes[1]):
            for z in range(dim_sizes[2]):
                pos = [x, y, z]
                for idx, param in enumerate(dim_params):
                    if param is not None:
                        param_value = dim_values[idx][pos[idx]]
                        fixed_params[param] = param_value
                        if param_rois:
                            param_roi_data[param][x*patchsize:(x+1)*patchsize, y*patchsize:(y+1)*patchsize, z*patchsize:(z+1)*patchsize] = pos[idx]+1
                model_curve = fab.model_evaluate(rundata, fixed_params, nt)
                
                data[x*patchsize:(x+1)*patchsize, y*patchsize:(y+1)*patchsize, z*patchsize:(z+1)*patchsize,:] = model_curve
                #if noise is not None:
                #    # Add Gaussian noise
                #    signal_mean = np.mean(model_curve)
                #    noise_data = np.random.normal(0, signal_mean*noise, [patchsize, patchsize, patchsize, nt])
                #    data[x*patchsize:(x+1)*patchsize, y*patchsize:(y+1)*patchsize, z*patchsize:(z+1)*patchsize,:] += noise_data
                if patch_rois: 
                    patch_roi_data[x*patchsize:(x+1)*patchsize, y*patchsize:(y+1)*patchsize, z*patchsize:(z+1)*patchsize] = patch_label
                    patch_label += 1

    if noise is not None:
        # Add Gaussian noise
        #mean_signal = np.mean(data)
        noise = np.random.normal(0, noise, shape + [nt,])
        data += noise

    if patch_rois or param_rois: 
        ret = [data,] 
        if patch_rois: ret.append(patch_roi_data)
        if param_rois: 
            ret.append(param_roi_data)
        return tuple(ret)
    else: return data

def _find_file(f, envdir, newf):
    if f is not None:
        return f
    elif envdir in os.environ:
        newfpath = os.path.join(os.environ[envdir], newf)
        if os.path.isfile(newfpath):
            return newfpath
        else:
            return f
    else:
        return None

def find_fabber():
    """
    Find the Fabber executable, core library and model libraries, or return None if not found

    Returns a tuple of exec, core lib, list of models
    """
    ex, lib, models = None, None, []
    for envdir in ("FABBERDIR", "FSLDIR"):
        ex = _find_file(ex, envdir, _bin_format % "fabber")
        lib = _find_file(lib, envdir, _lib_format % "fabbercore_shared")
        models += glob.glob(os.path.join(os.environ.get(envdir, ""), _lib_format % "fabber_models_*"))

    return ex, lib, models

class FabberException(RuntimeError):
    """
    Thrown if there is an error using the Fabber executable or library
    """

    def __init__(self, msg, errcode=None, log=None):
        self.errcode = errcode
        self.log = log
        if errcode is not None:
            RuntimeError.__init__(self, "%i: %s" % (errcode, msg))
        else:
            RuntimeError.__init__(self, msg)

class RunNotFound(RuntimeWarning):
    """
    A directory looked like a Fabber output directory, but no logfile was found in it
    """

    def __init__(self, dir):
        RuntimeWarning.__init__(self, "Not a Fabber run directory: %s" % dir)
        self.dir = dir

class Model:
    """
    Object which can have views
    """ 
    def __init__(self, name):
        # Used to record what has changed before updating views
        self.changes = set()
        self.name = name
        self.CH_ALL = "all"
        
        # List of known views of the model
        self.views = [] 
     
    def changed(self, *stuff):
        """ 
        Find out if a particular item has changed, e.g. if fab.changed(CH_FOCUS)...
        
        Returns true if it has
        """
        if len(stuff) == 0 or self.CH_ALL in self.changes:
            return True
        else:
            for s in stuff:
                if s in self.changes:
                    return True   
    
    def add_view(self, view):
        """
        Add a new view and update it
        """
        self.views.append(view)
        self._change()
        view.update(self)
        self.changes.clear()
        
    def _change(self, *stuff):
        """ 
        Record a _change, e.g. self._change(CH_INDATA, CH_OPTIONS). Views will
        use this info to decide what they need to update
        """        
        for s in stuff:
             self.changes.add(s)
        if len(stuff) == 0:
             self.changes.add(self.CH_ALL)
             
    def _update_views(self):
        """
        Update all views
        """
        #print "Updating, changes=" + str(self.changes)
        for view in self.views:
            try:
                view.update(self)
            except:
                traceback.print_exc()
        self.changes.clear()
       
class View:
    """ 
    Object which views a model
    """
    def __init__(self, changes=[], *widgets, **kwidgets):
        self.changes = set(changes)
        self.widgets = [w for w in widgets if self._iswidget(w)]
        for name, w in kwidgets.items():
            setattr(self, name, w)
            if self._iswidget(w):
                self.widgets.append(w)
        self.update(None)

    def update(self, obj):
        try:
            for widget in self.widgets:
                widget.blockSignals(True)
            if obj is None: self.set_enabled(False)
            elif not hasattr(self, obj.name) or obj.changed(*self.changes):
                self.set_enabled(True)
                setattr(self, obj.name, obj)
                self.do_update()
        finally:
            for widget in self.widgets:
                widget.blockSignals(False)

    def set_enabled(self, enabled=True, widgets=None):
        if widgets is None: widgets = self.widgets
        for widget in widgets:
                widget.setEnabled(enabled)

    def _iswidget(self, w):
        """
        Crude check to see if this is a widget! Don't want to use
        isinstance and add explicit dependency on QT
        """
        return hasattr(w, "blockSignals") and hasattr(w, "setEnabled")

    def do_update(self, obj_name):
        pass


class FabberRunData(Model, collections.MutableMapping):
    """
    Options for a Fabber run

    Behaves like a dict but supports save to/load from file,
    default values and keeps list of file lines separate
    so comments can be preserved
    """

    def __init__(self, filename=None):
        Model.__init__(self, "rundata")

        # Mapping of options as key/value pairs
        self.options = {}

        # Lines in the file in order, either a key or a comment string
        # Enables file to be recreated with comments in right place
        self.filelines = []

        if filename:
            self.load(filename)
        else:
            self._init_default()

            # self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getstate__(self):
        return (self.options, self.filelines, self.filepath)

    def __setstate__(self, state):
        self.options, self.filelines, self.filepath = state

    def __getitem__(self, key):
        return self.options[key.strip()]

    def __setitem__(self, key, value):
        key = key.strip()
        if key == "": return
        if key not in self.options:
            self.filelines.append(key)
        self.options[key] = value.strip()
        self._change(key)
        self._update_views()

    def __delitem__(self, key):
        key = key.strip()
        if key == "": return
        if key in self.options:
            del self.options[key]
            self.filelines.remove(key)
        self._change(key)
        self._update_views()

    def __iter__(self):
        return iter(self.options)

    def __len__(self):
        return len(self.options)

    def add_comment(self, comment, option=None):
        if option is None:
            idx = 0
        else:
            idx = self.filelines.index(option)

        self.filelines.insert(idx, "# %s" % str(comment))

    def set_file(self, f):
        """ Set the file. Do not parse, this is just used
            as part of 'save as' """
        self.filepath = os.path.abspath(f)
        # FIXME update views?

    def get_filedir(self):
        dir, junk = os.path.split(self.filepath)
        return dir

    def get_filename(self):
        junk, fname = os.path.split(self.filepath)
        return fname

    def load(self, fname):
        self.set_file(fname)
        fabfile = open(self.filepath, "r")
        for line in fabfile.readlines():
            line = line.strip()
            if len(line) == 0 or line[0] == "#":
                self.filelines.append(line)
            else:
                kv = line.split("=", 1)
                key = kv[0].strip()
                if len(kv) > 1:
                    value = kv[1].strip()
                else:
                    value = ""
                self.options[key] = value
                self.filelines.append(key)

        fabfile.close()
        self._change()
        self._update_views()

    def save(self, fpath=None, mask=None):
        """
        Save .fab file. If mask specified use this mask file
        in preference to one in options (used for quick runs)
        """
        if fpath is None:
            fpath = self.filepath
        if not fpath:
            raise RuntimeError("No file specified for save")

        fab = open(fpath, "w")
        self.dump(fab, mask=mask)
        fab.close()

    def dump(self, stream, mask=None):
        """
        Dump to an output stream

        :param stream: Output stream
        :param mask: If specified, replace the mask option with this value. Used for temporary quick runs
        """
        for line in self.filelines:
            if len(line) == 0 or line[0] == "#":
                stream.write(line)
            else:
                if self.options[line] == "":
                    stream.write("%s" % line)
                elif line == "mask" and mask is not None:
                    stream.write("mask=%s" % mask)
                else:
                    stream.write("%s=%s" % (line, self.options[line]))
            stream.write("\n")

    def _init_default(self):
        self.set_file("newfile.fab")

        self.options = {"method": "vb", "noise": "white",
                        "model": "poly", "degree": "2",
                        "output": "out", "overwrite": "",
                        "save-model-fit": ""}
        self.filelines = ["# Generated by fabber_gui.py", "",
                          "output", "overwrite", "save-model-fit", ""
                                                                   "method", "noise", "model", "degree"]
        self._change()
        self._update_views()

class FabberRun:
    """
    Base class for a completed Fabber run, either from the executable or from the library
    """

    def get_log_timestamp(self, log):
        prefixes = ["start time:", "fabberrundata::start time:"]
        timestamp_str = ""
        for line in log.splitlines():
            l = line.strip()
            for prefix in prefixes:
                if l.lower().startswith(prefix):
                    timestamp_str = l[len(prefix):].strip()
                    try:
                        timestamp = time.strptime(timestamp_str)
                        return timestamp, timestamp_str
                    except:
                        warnings.warn("Failed to parse timestamp: '%s'" % timestamp_str)

        warnings.warn("Could not find timestamp in log")
        return datetime.datetime.now(), timestamp_str

class LibRun(FabberRun):
    """
    A fabber library run, with output data and log
    """

    def __init__(self, data, log):
        self.data = data
        self.log = log
        self.timestamp, self.timestamp_str = self.get_log_timestamp(self.log)


class DirectoryRun(FabberRun):
    """
    A run of the fabber executable, with its output directory, logfile and output data.

    The data is not loaded by default, but load_data and load_all_data can be used to do this.
    """

    def __init__(self, dir, load_data=False):
        self.dir = dir
        if not self._is_fabber_dir():
            raise RunNotFound(dir)

        self.files = {}
        self.data = {}
        self.logfile, self.log = self._get_log()
        self.timestamp, self.timestamp_str = self.get_log_timestamp(self.log)
        self.isquick = self._is_quick_run()
        self.params = self._get_params()
        self._scan_output()
        if load_data: self.load_all_data()

    def load_data(self, filename):
        if name not in self.data:
            try:
                self.data[name] = nib.load(self.files[name]).get_data()

            except:
                warnings.warn("Could not load data file: ", f)
        return self.data[name]

    def load_all_data(self):
        for name in self.files:
            self.load_data(name)

    def _is_fabber_dir(self):
        return os.path.isfile(os.path.join(self.dir, "logfile"))

    def _get_log(self):
        logfile = os.path.join(self.dir, "logfile")
        f = open(logfile)
        log = "".join(f.readlines())
        f.close()
        return logfile, log

    def _scan_output(self):
        try:
            outdir_files = [f for f in os.listdir(self.dir)
                            if os.path.isfile(os.path.join(self.dir, f)) and f.endswith(".nii.gz")]
        except:
            warnings.warn("Could not read output directory: ", self.dir)
            traceback.print_exc()
            return

        for fname in outdir_files:
            try:
                name = fname.split(".")[0]
                f = os.path.join(self.dir, fname)
                self.files[name] = f
            except:
                warnings.warn("Not a valid output data file: ", fname)

    def _is_quick_run(self):
        return os.path.isfile(os.path.join(self.dir, "QUICKRUN.txt"))

    def _get_params(self):
        try:
            f = open(os.path.join(self.dir, "paramnames.txt"), "r")
            params = set([p.strip() for p in f.readlines()])
            f.close()
            return params
        except:
            warnings.warn("Failed to get parameters")
            return set()
   
class Fabber:
    def __init__(self, ex=None, lib=None, model_libs=[], rundata=None, auto_load_models=False):
        def_ex, def_lib, models = find_fabber()

        if ex is not None:
            self.ex = ex
            
        elif rundata is not None and "fabber" in rundata:
            self.ex = rundata["fabber"]
        else:
            self.ex = def_ex
        if not os.path.isfile(self.ex):
            raise FabberException("Invalid executable - file not found: %s" % self.ex)

        if lib is not None:
            self.lib = lib
        elif rundata is not None and "fabber_lib" in rundata:
            self.lib = rundata["fabber_lib"]
        else:
            self.lib = def_lib
        if not os.path.isfile(self.lib):
            raise FabberException("Invalid core library - file not found: %s" % self.lib)

        self.model_libs = set(model_libs)
        if rundata is not None and "loadmodels" in rundata:
            self.model_libs.add(rundata["loadmodels"])
        for lib in self.model_libs:
            if not os.path.isfile(lib):
                raise FabberException("Invalid models library - file not found: %s" % lib)
 
        if auto_load_models:
            for model in models:
                self.model_libs.add(model)

class FabberExec(Fabber):
    """
    Encapsulates a Fabber executable

    Provides methods to query models and options and also run a file
    """

    def __init__(self, ex=None, model_libs=[], rundata=None, auto_load_models=False):
        Fabber.__init__(self, ex=ex, model_libs=model_libs, rundata=rundata, auto_load_models=auto_load_models)
            
    def get_methods(self):
        """ Get known inference methods """
        stdout = self._run_help("--listmethods")
        return [line.strip() for line in stdout.split()]

    def get_models(self):
        """ Get known models """
        stdout = self._run_help("--listmodels")
        return [line.strip() for line in stdout.split()]

    def get_options(self, method=None, model=None):
        """
        Get general options, or options for a method/model. Returns a list of dictionaries
        with keys: name, type, optional, default, description
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
                l = [o.strip() for o in line2.split(",")]
                if len(l) < 4:
                    continue
                opt = {}
                opt["name"] = l[0]
                opt["type"] = l[1]
                opt["optional"] = l[2] != "REQUIRED"
                if l[3].startswith("DEFAULT="):
                    opt["default"] = l[3].split("=", 1)[1]
                else:
                    opt["default"] = ""
                if opt["name"] in ["data", "data<n>", "mask", "help", "model", "method", "listmodels", "listmethods"]:
                    continue
                if len(l) == 4:
                    opts.append(opt)
                    descnext = True
                else:
                    descnext = False
                    raise FabberException("Couldn't _parse option: " + line)
            elif descnext:
                opt["description"] = line.strip()
                descnext = False
            elif not line.startswith("Usage") and not line.startswith("Options") and not line.startswith("Description"):
                desc += line.strip()

        return opts, desc

    def get_model_params(self, rundata):
        raise FabberException("get_params not implemented for executable interface")

    def run(self, rundata):
        """
        Run Fabber on the run data specified
        """
        rundata.save()
        workdir = rundata.get_filedir()
        cmd = [self.ex, "-f", rundata.filepath]
        for lib in self.model_libs:
            cmd += " --loadmodels=%s" % lib
        err = ""
        p = sub.Popen(cmd, stdout=sub.PIPE, stderr=sub.PIPE, cwd=workdir)
        while 1:
            (stdout, stderr) = p.communicate()
            status = p.poll()
            if stderr: err += stderr
            if status is not None: break

            if status != 0:
                raise Exception(err)

        # Hack to get the last run
        outdir = rundata["output"]
        if not os.path.isabs(outdir):
            outdir = os.path.join(workdir, outdir)

        return self.get_previous_runs(outdir)[0]

    def get_previous_runs(self, outdir):
        """
        Get a list of known previous runs. The first in the list is the most recent

        This starts with the specified output directory and tries to load a run.
        If successful, it adds + to the filename (as Fabber does) and tries again.
        It continues until it cannot find any more valid run directories.
        """
        runs = []
        while os.path.exists(outdir):
            try:
                run = DirectoryRun(outdir)
                runs.append(run)
            except RunNotFound:
                warnings.warn("Could not find logfile in: %s - ignoring" % outdir)
            outdir += "+"
        runs.sort(key=lambda run: run.timestamp)
        return runs

    def _write_temp_mask(self):
        pass

    """ FIXME
        fname = os.path.join(self.fabdir, "fabber_mask_temp.nii.gz")
        data = np.zeros(self.shape[:3])
        data[self.focus[0], self.focus[1], self.focus[2]] = 1
        affine = self.data["data"].affine # FIXME
        img = nib.Nifti1Image(data, affine)
        nib.save(img, fname)
        return fname
    """

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
        cmd = [self.ex] + list(opts)
        for lib in self.model_libs:
            cmd += " --loadmodels=%s" % lib

        try:
            p = sub.Popen(cmd, stdout=sub.PIPE)
            (stdout, stderr) = p.communicate()
            status = p.wait()
            if status == 0:
                return stdout
            else:
                raise FabberException("Failed to run fabber: " + str(cmd), status)
        except Exception, e:
            raise FabberException("Failed to run fabber: " + str(e))

class FabberLib(Fabber):
    """
    Interface to Fabber in library mode using simplified C-API
    """

    def __init__(self, lib=None, model_libs=[], rundata=None, auto_load_models=False):
        Fabber.__init__(self, lib=lib, model_libs=model_libs, rundata=rundata, auto_load_models=auto_load_models)

        self.errbuf = create_string_buffer(255)
        self.outbuf = create_string_buffer(1000000)
        self.progress_cb_type = CFUNCTYPE(None, c_int, c_int)
        self._init_clib()

    def get_methods(self):
        """ Get known inference methods"""
        self._trycall(self.clib.fabber_get_methods, self.handle, len(self.outbuf), self.outbuf, self.errbuf)
        return self.outbuf.value.splitlines()

    def get_models(self):
        """ Get known models"""
        self._trycall(self.clib.fabber_get_models, self.handle, len(self.outbuf), self.outbuf, self.errbuf)
        return self.outbuf.value.splitlines()

    def get_options(self, method=None, model=None):
        """
        Get known Fabber options

        :param method: If specified, return options for this method
        :param model: If specified, return options for this model

        Only one of method and model should be specified. If neither are specified, generic
        Fabber options are returned.

        :return: Tuple of options, description. Options is a list of options, each in the form of a dictionary.
        Description is a simple text description of the method or model
        """
        if method:
            key = "method"
            value = method
        elif model:
            key = "model"
            value = model
        else:
            key = None
            value = None
        self._trycall(self.clib.fabber_get_options, self.handle, key, value, len(self.outbuf), self.outbuf, self.errbuf)
        opt_keys = ["name", "description", "type", "optional", "default"]
        opts = []
        lines = self.outbuf.value.split("\n")
        for opt in lines[1:]:
            if len(opt) > 0:
                opt = dict(zip(opt_keys, opt.split("\t")))
                opt["optional"] = opt["optional"] == "1"
                opts.append(opt)
        return opts, lines[0]

    def get_model_params(self, rundata):
        """ Get the model parameters, given the specified options"""
        for key, value in rundata.items():
            self._trycall(self.clib.fabber_set_opt, self.handle, key, value, self.errbuf)

        self._trycall(self.clib.fabber_get_model_params, self.handle, len(self.outbuf), self.outbuf, self.errbuf)

        # Reset context because we have set options and don't want them affecting a later call to run()
        self._init_clib()
        return self.outbuf.value.splitlines()

    def model_evaluate(self, rundata, params, nt, indata=None):
        """ """
        for key, value in rundata.items():
            self._trycall(self.clib.fabber_set_opt, self.handle, str(key), str(value), self.errbuf)

        # Get model parameter names
        self._trycall(self.clib.fabber_get_model_params, self.handle, len(self.outbuf), self.outbuf, self.errbuf)
        model_params = self.outbuf.value.splitlines()
        if len(params) != len(model_params):
            raise FabberException("Incorrect number of parameters specified: expected %i (%s)" % (len(model_params), ",".join(model_params)))
        plist = []
        for p in model_params:
            if p not in params:
                raise FabberException("Model parameter %s not specified" % p)
            else:
                plist.append(params[p])
        ret = np.zeros([nt,], dtype=np.float32)
        if indata is None: indata = np.zeros([nt,], dtype=np.float32)
        self._trycall(self.clib.fabber_model_evaluate, self.handle, len(plist), np.array(plist, dtype=np.float32), nt, indata, ret, self.errbuf)

        # Reset context because we have set options and don't want them affecting a later call to run()
        self._init_clib()
        return ret

    def run(self, rundata, progress_cb=None):
        """
        Run fabber on the provided rundata options

        :param rundata: FabberRunData instance
        :param progress_cb: Callable which will be called periodically during processing
        :return: On success, a LibRun instance
        """
        mask = None
        data = {}
        for key, value in rundata.items():
            try:
                f = nib.load(value)
                d = f.get_data()
                if key == "mask":
                    mask = d
                else:
                    data[value] = d
            except:
                if key == "data":
                    raise
                # Otherwise ignore, most options will not be data files
                pass

        return self.run_with_data(rundata, data, mask, progress_cb)

    def run_with_data(self, rundata, data, mask=None, progress_cb=None):
        """
        Run fabber

        :param data: Dictionary of data: string key, Numpy array value
        :param mask: Mask as Numpy array, or None if no mask
        :param progress_cb: Callable which will be called periodically during processing
        :return: On success, a LibRun instance
        """
        if not data.has_key("data"):
            raise Exception("Main voxel data not provided")
        s = data["data"].shape
        nv = s[0] * s[1] * s[2]

        if mask is None: mask = np.ones(nv)
        # Make suitable for passing to int* c function
        mask = np.ascontiguousarray(mask.flatten(), dtype=np.int32)

        for key, value in rundata.items():
            self._trycall(self.clib.fabber_set_opt, self.handle, str(key), str(value), self.errbuf)
        self._trycall(self.clib.fabber_get_model_params, self.handle, len(self.outbuf), self.outbuf, self.errbuf)
        params = self.outbuf.value.splitlines()

        output_items = []
        if "save-mean" in rundata:
            output_items += ["mean_" + p for p in params]
        if "save-std" in rundata:
            output_items += ["std_" + p for p in params]
        if "save-zstat" in rundata:
            output_items += ["zstat_" + p for p in params]
        if "save-noise-mean" in rundata:
            output_items.append("noise_mean")
        if "save-noise-std" in rundata:
            output_items.append("noise_stdevs")
        if "save-free-energy" in rundata:
            output_items.append("freeEnergy")
        if "save-model-fit" in rundata:
            output_items.append("modelfit")
        if "save-residuals" in rundata:
            output_items.append("residuals")
        if "save-mvn" in rundata:
            output_items.append("finalMVN")

        retdata, log = {}, ""
        self._trycall(self.clib.fabber_set_extent, self.handle, s[0], s[1], s[2], mask, self.errbuf)
        for key, item in data.items():
            if len(item.shape) == 3:
                size = 1
            else:
                size = item.shape[3]
            item = np.ascontiguousarray(item.flatten(), dtype=np.float32)
            self._trycall(self.clib.fabber_set_data, self.handle, key, size, item, self.errbuf)

        progress_cb_func = self.progress_cb_type(0)
        if progress_cb is not None:
            progress_cb_func = self.progress_cb_type(progress_cb)

        self._trycall(self.clib.fabber_dorun, self.handle, len(self.outbuf), self.outbuf, self.errbuf, progress_cb_func)
        log = self.outbuf.value
        for key in output_items:
            size = self._trycall(self.clib.fabber_get_data_size, self.handle, key, self.errbuf)
            arr = np.ascontiguousarray(np.empty(nv * size, dtype=np.float32))
            self._trycall(self.clib.fabber_get_data, self.handle, key, arr, self.errbuf)
            if size > 1:
                arr = arr.reshape([s[0], s[1], s[2], size])
            else:
                arr = arr.reshape([s[0], s[1], s[2]])
            retdata[key] = arr

        return LibRun(retdata, log)

    def __del__(self):
        self._destroy_handle()

    def _destroy_handle(self):
        if hasattr(self, "handle"):
            handle = getattr(self, "handle")
            if handle is not None:
                self.clib.fabber_destroy(handle)
                self.handle = None

    def _init_clib(self):
        """
        This is required because currently there is no CAPI function to clear the rundata
        of options
        """
        try:
            self._destroy_handle()
            self.clib = CDLL(str(self.lib))

            # Signatures of the C functions
            c_int_arr = npct.ndpointer(dtype=np.int32, ndim=1, flags='CONTIGUOUS')
            c_float_arr = npct.ndpointer(dtype=np.float32, ndim=1, flags='CONTIGUOUS')

            self.clib.fabber_new.argtypes = [c_char_p]
            self.clib.fabber_new.restype = c_void_p
            self.clib.fabber_load_models.argtypes = [c_void_p, c_char_p, c_char_p]
            self.clib.fabber_set_extent.argtypes = [c_void_p, c_int, c_int, c_int, c_int_arr, c_char_p]
            self.clib.fabber_set_opt.argtypes = [c_void_p, c_char_p, c_char_p, c_char_p]
            self.clib.fabber_set_data.argtypes = [c_void_p, c_char_p, c_int, c_float_arr, c_char_p]
            self.clib.fabber_get_data_size.argtypes = [c_void_p, c_char_p, c_char_p]
            self.clib.fabber_get_data.argtypes = [c_void_p, c_char_p, c_float_arr, c_char_p]
            self.clib.fabber_dorun.argtypes = [c_void_p, c_int, c_char_p, c_char_p, self.progress_cb_type]
            self.clib.fabber_destroy.argtypes = [c_void_p]

            self.clib.fabber_get_options.argtypes = [c_void_p, c_char_p, c_char_p, c_int, c_char_p, c_char_p]
            self.clib.fabber_get_models.argtypes = [c_void_p, c_int, c_char_p, c_char_p]
            self.clib.fabber_get_methods.argtypes = [c_void_p, c_int, c_char_p, c_char_p]

            self.clib.fabber_get_model_params.argtypes = [c_void_p, c_int, c_char_p, c_char_p]
            self.clib.fabber_model_evaluate.argtypes = [c_void_p, c_int, c_float_arr, c_int, c_float_arr, c_float_arr, c_char_p]
        except Exception, e:
            raise FabberException("Error initializing Fabber library: %s" % str(e))

        self.handle = self.clib.fabber_new(self.errbuf)
        if self.handle is None:
            raise FabberException("Error creating fabber context (%s)" % self.errbuf.value)

        for lib in self.model_libs:
            self._trycall(self.clib.fabber_load_models, self.handle, lib, self.errbuf)

    def _trycall(self, call, *args):
        ret = call(*args)
        if ret < 0:
            raise FabberException(self.errbuf.value, ret, self.outbuf.value)
        else:
            return ret
