import traceback
import collections

import nibabel as nib
import numpy as np

import pylab as plt
from PySide import QtCore, QtGui

from mvc import Model, View

"""Change to the  data"""
CH_DATA = "data"

"""Change to the focus point"""
CH_FOCUS = "focus"

"""Change to the selected data item"""
CH_SEL = "sel"

class ImageView:
    """
    Represents a view of a 3d or 4d volume,
    input or output data
    """
    def __init__(self, data=None, filename=None):
        if data is None and filename is None:
            raise RuntimeError("Can't create ImageView without data or filename")
        elif data is not None and filename is not None:
                raise RuntimeError("Can't create ImageView with both data and filename")
        elif data is None:
            self.filename = filename
            data = nib.load(filename).get_data()
        else:
            self.filename = None

        self.data = np.ma.masked_values(data, 0)
        self.shape = self.data.shape 
        self.ndims = len(self.shape)
        self.alpha = 1.0
        self.visible = True
        self.cm = "gray"
        self.min = self.data.min()
        self.max = self.data.max()
        
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

class FabberImageData(Model, View, collections.MutableMapping):
    """
    Input / output image data for Fabber

    Behaves like a dict of ImageView items
    """
    def __init__(self):
        Model.__init__(self, "imdata")
        View.__init__(self, [])

        # Mapping of options as key/value pairs
        self.imdata = {}

        # Shape of main data
        self.shape = None

        # Point of focus
        self.focus = [0, 0, 0, 0]

        # Selected data item which is the current focus.
        self.selected_data = None

        # Current run
        self.run = None

        # self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.imdata[key.strip()]

    def __setitem__(self, key, item):
        self.imdata[key.strip()] = item

    def __delitem__(self, key):
        key = key.strip()
        if key in self.imdata:
            del self.imdata[key]
        self._change(CH_DATA)
        self._update_views()

    def __iter__(self):
        return iter(self.imdata)

    def __len__(self):
        return len(self.imdata)

    def do_update(self):
        """
        Run data has changed
        """
        print("imdata: rundata changed", self.rundata)
        if self.rundata.changed("data"):
            self._update_input_data("data")
            if "data" in self:
                self.shape = self["data"].shape
                self._set_default_focus()
        elif self.rundata.changed("mask"):
            self._update_input_data("mask")
        self._update_views()

    def set_run(self, run):
        """ 
        Load output data from a run
        """
        if self.shape is None:
            raise RuntimeError("Can't load a run without main data")
        self.run = run
        self._clear_output_data()
        for name, filename in self.run.files.items():
            item = ImageView(filename=filename)
            self._check_dims(item.data)
            self._set_data_defaults(item, name)
            self[name] = item

    def set_visibility(self, item, alpha=None, visible=None, cm=None):
        """ 
        Set visibility of a data item 
        """
        if alpha is not None: item.alpha = alpha
        if visible is not None: item.visible = visible
        if cm is not None: item.cm = cm
        self._change(CH_DATA)
        self._update_views()

    def select_data(self, key):
        if key in self:
            self.selected_data = self[key]
        else:
            raise RuntimeWarning("Data not found: %s" % key)
        self._change(CH_SEL)        
        self._update_views()

    def update_focus(self, xp=None, yp=None, zp=None, tp=None):
        if xp: self.focus[0] = xp
        if yp: self.focus[1] = yp
        if zp: self.focus[2] = zp
        if tp: self.focus[3] = tp
        self._change(CH_FOCUS)
        self._update_views()
   
    def _update_input_data(self, key):
        print("Updating input", self.rundata)
        if key in self and key not in self.rundata:
            del self[key]
        elif key in self.rundata and (key not in self or self.rundata[key] != self[key].filename):
            self[key] = ImageView(filename=self.rundata[key])

    def _set_default_focus(self):
       self.focus = [d/2 for d in self.shape]
       self._change(CH_FOCUS)

    def _set_data_defaults(self, item, role):
        if role == "mask":
            item.visible = True
            item.alpha = 0.3
            item.cm = "autumn"
        elif role == "data":
            item.visible = True
            item.cm = "gray"
        else:
            item.visible = False
            item.cm = "spectral"

    def _clear_output_data(self):
        for key in self.keys():
            if key not in ("data", "mask"):
                del self[key]

    def _check_dims(self, item):
        """ 
        Check that the dimensions of the given data are consistent with the
        given shape. Only compares up to the largest common dimension, and
        ignores dimensions of 0 (these have not been initialized yet).
        
        Updates the existing shape with the number of dimensions expanded
        to account for the data. All future data must be consistent with this
        expanded set of dimensions 
        """
        if self.data is None: return

        for d in range(min(len(self.data.shape), len(item.shape))):
            if item.shape[d] != self.shape[d]:
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
        self.imdata.update_focus(xp=value)
    def y_changed(self, value):
        self.imdata.update_focus(yp=value)
    def z_changed(self, value):
        self.imdata.update_focus(zp=value)
    def t_changed(self, value):
        self.imdata.update_focus(tp=value)
        
    def update_widget(self, w, idx):
        w.setMinimum(0)
        if self.imdata.shape is not None:
            w.setMaximum(self.imdata.shape[idx]-1)
        else:
            w.setMaximum(1)
        w.setValue(self.imdata.focus[idx])
        
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
        View.__init__(self, [CH_FOCUS, CH_DATA], **kwargs)
        
    def populate(self, run):
        self.table.setRowCount(len(run.params))
        print run.params
        for idx, param in enumerate(run.params):
            mean = self.imdata[param + " Mean value"].get_value(self.imdata.focus)
            std = self.imdata[param + " Std. dev."].get_value(self.imdata.focus)
            self.table.setItem(idx, 0, QtGui.QTableWidgetItem(param))
            self.table.setItem(idx, 1, QtGui.QTableWidgetItem(str(mean)))
            self.table.setItem(idx, 2, QtGui.QTableWidgetItem(str(std)))
            
    def do_update(self):
        self.table.setRowCount(0)
        if self.imdata.run:
            self.populate(self.imdata.run)
         
class DataView(View):
    def __init__(self, logDialog, runDialog, **kwargs):
        View.__init__(self, [CH_DATA,], **kwargs)
        self.logDialog = logDialog
        self.runDialog = runDialog
        self.table.setHorizontalHeaderItem(0, QtGui.QTableWidgetItem("Role"))
        self.table.setHorizontalHeaderItem(1, QtGui.QTableWidgetItem("Filename"))
        self.chooseRunBtn.clicked.connect(self.choose_run)
        self.logBtn.clicked.connect(self.view_log_pressed)
        self.table.currentItemChanged.connect(self.data_select)
        self.dataBtn.clicked.connect(self.data_clicked)
        self.maskBtn.clicked.connect(self.mask_clicked)
        
    def data_clicked(self):
        filename = QtGui.QFileDialog.getOpenFileName()[0]
        if filename: self.imdata.rundata["data"] = filename
        
    def mask_clicked(self):
        filename = QtGui.QFileDialog.getOpenFileName()[0]
        if filename: self.imdata.rundata["mask"] = filename

    def choose_run(self):
        runs = self.imdata.get_runs()
        self.runDialog.table.setRowCount(len(runs))
        for idx, run in enumerate(runs):
            self.runDialog.table.setItem(idx, 0, QtGui.QTableWidgetItem(run.dir))
            self.runDialog.table.setItem(idx, 1, QtGui.QTableWidgetItem(run.timestamp))
        self.runDialog.show()
        if self.runDialog.exec_():
            run = runs[self.runDialog.table.currentRow()]
            self.imdata.select_run(run.dir)
 
    def view_log_pressed(self):
        self.logDialog.show()
        self.logDialog.raise_()
    
    def data_select(self, cur, prev):
        self.imdata.select(self.table.item(cur.row(), 0).text())
        
    def relist(self):
        self.table.setRowCount(0)
        n = 0
        self.table.setRowCount(0)
        for key, data in self.imdata.items():
            self.table.setRowCount(n+1)
            #if data.visible:
            #    item.setIcon(self.table.style().standardIcon(QtGui.QStyle.SP_FileIcon))        
            self.table.setItem(n, 0, QtGui.QTableWidgetItem(key))
            self.table.setItem(n, 1, QtGui.QTableWidgetItem(data.filename))
            n += 1
        if self.imdata.run:
            self.logDialog.textBrowser.setText(self.imdata.run.get_log())
            self.currentRunEdit.setText(self.imdata.run.timestamp)
            if self.imdata.run.isquick:
                msgBox = QtGui.QMessageBox()
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
            
        self.set_enabled(self.imdata.selected_data is not None)
        
        if self.imdata.selected_data is not None:
            self.currentEdit.setText(self.imdata.selected_data.filename)
            self.slider.setValue(100*self.imdata.selected_data.alpha)
            if self.imdata.selected_data.visible:
                self.cb.setCheckState(QtCore.Qt.CheckState.Checked)
            else:
                self.cb.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.cmCombo.setCurrentIndex(self.maps.index(self.imdata.selected_data.cm))
            self.valueEdit.setText(str(self.imdata.selected_data.get_value(self.imdata.focus)))
        
    def cm_changed(self, idx):
        text = self.maps[idx]
        self.imdata.set_visibility(self.imdata.selected_data, cm=text)

    def alpha_changed(self, value):
        self.imdata.set_visibility(self.imdata.selected_data, alpha=float(value)/100)

    def visible_changed(self):
        self.imdata.set_visibility(self.imdata.selected_data, visible=self.cb.checkState() == QtCore.Qt.CheckState.Checked)
          
