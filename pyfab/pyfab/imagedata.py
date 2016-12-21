import os
import subprocess as sub
import time
import sys
import traceback
import distutils.spawn
import datetime

import nibabel as nib
import numpy as np

import pylab as plt

from PySide import QtCore, QtGui

from PySide.QtGui import QMessageBox, QLabel, QHBoxLayout, QLineEdit, QVBoxLayout, QFileDialog, QTableWidgetItem, QPlainTextEdit, QSpinBox, QCheckBox, QPushButton

from mvc import Model, View
from model import CH_RUN

"""Change to the  data"""
CH_DATA = "data"

"""Change to the focus point"""
CH_FOCUS = "focus"

"""Change to the selected data item"""
CH_SEL = "sel"

class DataItem:
    """
    Represents a 3d or 4d volume, either
    input or output data
    """
    def __init__(self, filename):
        self.filename = filename
        self._load()
        self.shape = self.data.shape 
        self.ndims = len(self.shape)
        self.alpha = 1.0
        self.min = self.data.min()
        self.max = self.data.max()
        
    def _load(self):
        print("Loading data: " + self.filename)
        f = nib.load(self.filename)
        self.data = np.ma.masked_values(f.get_data(), 0)
        self.affine = f.affine
        
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
         
class FabberImageData(Model, View):
    """
    A collection of image data items, consisting of input data
    and optionally output data as well
    """
    def __init__(self):
        Model.__init__(self, "fab")
        View.__init__(self, [])
        
        # Data items, keyed by role
        self.data = {}
        
        # Data extent. All data in a DataSet must have compatible extent
        self.shape = [0, 0, 0, 0]
        
        # Focus point. Although this relates to visual presentation it
        # is in the model because all views should respect the current 
        # focus
        self.focus = [0, 0, 0, 0]
        
        # Data item which is the current focus. 
        self.focus_data = None
        
        # Current focus run
        self.current_run = None
    
    def do_update(self):
        """
        Run data has changed
        """
        if self.fab.changed("data"):
            self._update_input_data("data")
            self._change(CH_DATA)
        elif self.fab.changed("mask"):
            self._update_input_data("mask")
            self._change(CH_DATA)
        if self.fab.changed(CH_RUN, "output"):
            self._get_last_run()
            self._load_run()
            self._change(CH_DATA)
        self._update_views()
        
    def set_rundir(self, outdir):
        runs = self.fab.get_runs()
        for run in runs:
            if run.outdir == outdir:
                self.run = run
                break
        self._load_run()
        self._change(CH_DATA)
        self._update_views()
        
    def _get_last_run(self):
        runs = self.fab.get_runs()
        if len(runs) > 0:
            self.run = runs[0]
        else: self.run = None

    def _load_run(self):
        """ 
        Load output data from the current run
        """
        if self.run:
            for role, filename in self.run.files.items():
                self._load_data(role, filename)

    def set_visibility(self, item, alpha=None, visible=None, cm=None):
        """ 
        Set visibility of a data item 
        """
        if alpha is not None: item.alpha = alpha
        if visible is not None: item.visible = visible
        if cm is not None: item.cm = cm
        self._change(CH_DATA)
        self._update_views()

    def select(self, key):
        """ 
        Select a data item 
        """
        if key in self.data:
            self.focus_data = self.data[key]
            print("Selected: ", self.focus_data)
        else:
            raise RuntimeWarning("Data not found in ImageDataSet: %s" % key)
        self._change(CH_SEL)        
        self._update_views()

    def update_focus(self, xp=None, yp=None, zp=None, tp=None):
        """ 
        Update the point of focus 
        """
        if xp: self.focus[0] = xp
        if yp: self.focus[1] = yp
        if zp: self.focus[2] = zp
        if tp: self.focus[3] = tp
        self._change(CH_FOCUS)        
        self._update_views()
   
    def _update_input_data(self, key):
        """
        Load a data file
        
        If there is a problem loading the data file, or the
        dimensions are inconsistent, the data will not be
        added and an error is displayed to stderr
        """
        if not self.fab.options.has_key(key):
            # Data does not exist
            if self.data.has_key(key): del self.data[key]
        elif self.data.has_key(key) and self.data[key].filename == self.fab.options[key]:
            # Data file is unchanged
            return
        else:
            self._load_data(key, self.fab.options[key])
        
    def _load_data(self, key, filename):
        print("loading %s from %s" % (key, filename))
        try:
            item = DataItem(filename)
            self._set_data_defaults(item, key)
            self._check_dims(item)
            self.data[key] = item
            for idx, f in enumerate(self.focus):
                if f == 0: 
                    self.focus[idx] = self.shape[idx] / 2
                    self._change(CH_FOCUS)
        except:
            print("Failed to load data: %s, %s" % (key, filename))
            traceback.print_exc()
  
    def _set_data_defaults(self, item, role):
        if role == "mask": 
            item.visible = True
            item.alpha = 0.3
            item.cm = "autumn"
        elif role == "data":
            item.visible = True
            item.cm = "gray"
        elif role.lower() == "model prediction":
            item.visible = True
            item.cm = "spectral"
        else:
            item.visible = False
            item.cm = "spectral"
            
    def _check_dims(self, item):
        """ 
        Check that the dimensions of the given data are consistent with the
        given shape. Only compares up to the largest common dimension, and
        ignores dimensions of 0 (these have not been initialized yet).
        
        Updates the existing shape with the number of dimensions expanded
        to account for the data. All future data must be consistent with this
        expanded set of dimensions 
        """
        for d in range(min(len(self.shape), len(item.shape))):
            if self.shape[d] == 0:
                self.shape[d] = item.shape[d]
            elif item.shape[d] != self.shape[d]:
                raise Exception("Inconsistent Dimensions: %i (%i != %i)" % (d, item.shape[d], self.shape[d]))

class FocusView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [CH_FOCUS,], **kwargs)
        self.slider.valueChanged.connect(self.t_changed)
        self.sbX.valueChanged.connect(self.x_changed)
        self.sbY.valueChanged.connect(self.y_changed)
        self.sbZ.valueChanged.connect(self.z_changed)
        self.sbT.valueChanged.connect(self.t_changed)
        
    def x_changed(self, value):
        self.fab.update_focus(xp=value)
    def y_changed(self, value):
        self.fab.update_focus(yp=value)
    def z_changed(self, value):
        self.fab.update_focus(zp=value)
    def t_changed(self, value):
        self.fab.update_focus(tp=value)
        
    def update_widget(self, w, idx):
        w.setMinimum(0)
        w.setMaximum(self.fab.shape[idx]-1)
        w.setValue(self.fab.focus[idx])
        
    def do_update(self):
        self.update_widget(self.slider, 3)
        self.update_widget(self.sbX, 0)
        self.update_widget(self.sbY, 1)
        self.update_widget(self.sbZ, 2)
        self.update_widget(self.sbT, 3)
        
class ParamValuesView(View):
    """
    Table of parameter values at the focus point
    """
    def __init__(self, **kwargs):
        View.__init__(self, [CH_FOCUS, CH_DATA, CH_RUN], **kwargs)
        
    def populate(self, run):
        self.table.setRowCount(len(run.params))
        print run.params
        for idx, param in enumerate(run.params):
            mean = self.fab.data[param + " Mean value"].get_value(self.fab.focus)
            std = self.fab.data[param + " Std. dev."].get_value(self.fab.focus)
            self.table.setItem(idx, 0, QTableWidgetItem(param))
            self.table.setItem(idx, 1, QTableWidgetItem(str(mean)))
            self.table.setItem(idx, 2, QTableWidgetItem(str(std)))
            
    def do_update(self):
        self.table.setRowCount(0)
        if self.fab.run:
            self.populate(self.fab.run)
         
class DataView(View):
    def __init__(self, logDialog, runDialog, **kwargs):
        View.__init__(self, [CH_DATA,], **kwargs)
        self.logDialog = logDialog
        self.runDialog = runDialog
        self.table.setHorizontalHeaderItem(0, QTableWidgetItem("Role"))    
        self.table.setHorizontalHeaderItem(1, QTableWidgetItem("Filename"))    
        self.chooseRunBtn.clicked.connect(self.choose_run)
        self.logBtn.clicked.connect(self.view_log_pressed)
        self.table.currentItemChanged.connect(self.data_select)
        self.dataBtn.clicked.connect(self.data_clicked)
        self.maskBtn.clicked.connect(self.mask_clicked)
        
    def data_clicked(self):
        filename = QFileDialog.getOpenFileName()[0]
        if filename: self.fab.fab.set_option("data", filename)
        
    def mask_clicked(self):
        filename = QFileDialog.getOpenFileName()[0]
        if filename: self.fab.fab.set_option("mask", filename)
        
    def choose_run(self):
        runs = self.fab.get_runs()
        self.runDialog.table.setRowCount(len(runs))
        for idx, run in enumerate(runs):
            self.runDialog.table.setItem(idx, 0, QTableWidgetItem(run.dir))
            self.runDialog.table.setItem(idx, 1, QTableWidgetItem(run.timestamp))
        self.runDialog.show()
        if self.runDialog.exec_():
            run = runs[self.runDialog.table.currentRow()]
            self.fab.select_run(run.dir)
 
    def view_log_pressed(self):
        self.logDialog.show()
        self.logDialog.raise_()
    
    def data_select(self, cur, prev):
        self.fab.select(self.table.item(cur.row(), 0).text())
        
    def relist(self):
        self.table.setRowCount(0)
        n = 0
        self.table.setRowCount(0)
        for key, data in self.fab.data.items():
            self.table.setRowCount(n+1)
            #if data.visible:
            #    item.setIcon(self.table.style().standardIcon(QtGui.QStyle.SP_FileIcon))        
            self.table.setItem(n, 0, QTableWidgetItem(key))
            self.table.setItem(n, 1, QTableWidgetItem(data.filename))
            n += 1
        if self.fab.run:
            self.logDialog.textBrowser.setText(self.fab.run.get_log())
            self.currentRunEdit.setText(self.fab.run.timestamp)
            if self.fab.run.isquick:
                msgBox = QMessageBox()
                msgBox.setText("This was a quick 1-voxel run, not a full run")
                msgBox.exec_()

    def do_update(self):
        self.relist() 
  
class CurrentDataView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [CH_SEL,], **kwargs)
        self.cb.stateChanged.connect(self.visible_changed)
        self.cmCombo.currentIndexChanged.connect(self.cm_changed)
        self.slider.valueChanged.connect(self.alpha_changed)
        self.maps = None
            
    def do_update(self):
        if not self.maps:
            self.maps=[m for m in plt.cm.datad if not m.endswith("_r")]
            self.maps.sort()
            for map in self.maps:
                self.cmCombo.addItem(map)
            
        self.set_enabled(self.fab.focus_data is not None)
        
        if self.fab.focus_data is not None:
            self.currentEdit.setText(self.fab.focus_data.filename)
            self.slider.setValue(100*self.fab.focus_data.alpha)
            if self.fab.focus_data.visible:
                self.cb.setCheckState(QtCore.Qt.CheckState.Checked)
            else:
                self.cb.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.cmCombo.setCurrentIndex(self.maps.index(self.fab.focus_data.cm))
            self.valueEdit.setText(str(self.fab.focus_data.get_value(self.fab.focus)))
        
    def cm_changed(self, idx):
        text = self.maps[idx]
        self.fab.set_visibility(self.fab.focus_data, cm=text)

    def alpha_changed(self, value):
        self.fab.set_visibility(self.fab.focus_data, alpha=float(value)/100)

    def visible_changed(self):
        self.fab.set_visibility(self.fab.focus_data, visible=self.cb.checkState() == QtCore.Qt.CheckState.Checked)
          
