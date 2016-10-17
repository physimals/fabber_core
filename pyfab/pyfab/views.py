import sys
import traceback

from PySide import QtCore, QtGui

from PySide.QtGui import QMessageBox, QLabel, QHBoxLayout, QLineEdit, QVBoxLayout, QFileDialog, QTableWidgetItem, QPlainTextEdit, QSpinBox, QCheckBox, QPushButton

import pylab as plt

class ModelView:
    def __init__(self, *widgets, **kwidgets):
        self.widgets = list(widgets)
        for name, w in kwidgets.items():
            self.widgets.append(w)
            setattr(self, name, w)
            
        self.fab = None
        self.update(None)

    def update(self, fab):
        self.fab = fab
        try:
            self.set_enabled(self.fab is not None)            
            for widget in self.widgets:
                widget.blockSignals(True)
            if self.fab: self.do_update()    
        finally:
            for widget in self.widgets:
                widget.blockSignals(False)

    def set_enabled(self, enabled=True, widgets=None):
        if widgets is None: widgets = self.widgets
        for widget in widgets:
                widget.setEnabled(enabled)
            
    def do_update(self):
        pass

class ParamValuesView(ModelView):
    """
    Table of parameter values at the focus point
    """
    
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        
    def populate(self, run):
        self.table.setRowCount(len(run.params))
        print run.params
        for idx, param in enumerate(run.params):
            mean = run.data[param + " Mean value"].get_value(self.fab.focus)
            std = run.data[param + " Std. dev."].get_value(self.fab.focus)
            self.table.setItem(idx, 0, QTableWidgetItem(param))
            self.table.setItem(idx, 1, QTableWidgetItem(str(mean)))
            self.table.setItem(idx, 2, QTableWidgetItem(str(std)))
            
    def do_update(self):
        if self.fab.changed("focus", "data"):
            self.table.setRowCount(0)
            if self.fab.current_run:
                self.populate(self.fab.current_run)
            
class ExtraOptionsView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        self.models = None
        self.methods = None
        self.modelCombo.currentIndexChanged.connect(self.model_changed)
        self.methodCombo.currentIndexChanged.connect(self.method_changed)

    def model_changed(self):
        self.fab.set_option("model", self.modelCombo.currentText())
        
    def method_changed(self):
        self.fab.set_option("method", self.methodCombo.currentText())
        
    def do_update(self):
        if self.fab.changed("exec"): 
            self.modelCombo.clear()
            self.models = self.fab.get_models()
            for model in self.models:
                self.modelCombo.addItem(model)
        
            self.methodCombo.clear()
            self.methods = self.fab.get_methods()
            for method in self.methods:
                self.methodCombo.addItem(method)
                
        if self.fab.options.has_key("model"):
            cmodel = self.fab.options["model"]
            for idx, model in enumerate(self.models):
                if cmodel == model:
                    self.modelCombo.setCurrentIndex(idx)

        if self.fab.options.has_key("method"):
            cmethod = self.fab.options["method"]
            for idx, method in enumerate(self.methods):
                if cmethod == method:
                    self.methodCombo.setCurrentIndex(idx)

        
def get_label(text, size=None, bold=False, italic=False):
    label = QLabel(text)
    font = label.font()
    font.setBold(bold)
    font.setItalic(italic)
    if size: font.setPointSize(size)
    label.setFont(font)
    return label

class OptionView(ModelView):
    def __init__(self, opt, rescan=False):
        ModelView.__init__(self)
        self.key, self.dtype, self.req, self.default, self.desc = opt
        self.req = (self.req == "REQUIRED")
        self.rescan = rescan

        if self.req:
            self.label = get_label(opt[0], size=10)
        else:
            self.label = QCheckBox(opt[0])
            self.label.stateChanged.connect(self.state_changed)

        self.widgets.append(self.label)
        
        self.desclabel = get_label(opt[4], size=10)
        self.desclabel.resize(400, self.desclabel.height())
        #desclabel.setWordWrap(True)
        self.widgets.append(self.desclabel)

    def state_changed(self):
        self.set_enabled(self.label.checkState() == QtCore.Qt.CheckState.Checked, self.widgets[1:])
        if self.rescan: self.fab._change("exec") # FIXME a hack
        if self.label.checkState() == QtCore.Qt.CheckState.Checked:
            self.changed()
        else:
            self.fab.clear_option(self.key)            
    
    def changed(self):
        self.fab.set_option(self.key, "")

    def do_update(self):
        if not self.req:
            if not self.fab.options.has_key(self.key):
                self.label.setCheckState(QtCore.Qt.CheckState.Unchecked)
            else:
                self.label.setCheckState(QtCore.Qt.CheckState.Checked)
            self.set_enabled(self.label.checkState() == QtCore.Qt.CheckState.Checked, self.widgets[1:])

    def add(self, grid, row):
        grid.addWidget(self.label, row, 0)
        grid.addWidget(self.desclabel, row, 2)

    def is_numbered(self, optname):
        return "<n>" in optname

def get_option_view(opt, **kwargs):
    if opt[1] == "INTEGER":
        return IntegerOptionView(opt, **kwargs)
    elif opt[1] == "BOOL":
        return OptionView(opt, **kwargs)
    elif opt[1] == "FILENAME" :
        return FileOptionView(opt, **kwargs)
    elif opt[1] == "MATRIXFILE":
        return MatrixFileOptionView(opt, **kwargs)
    else:
        return StringOptionView(opt, **kwargs)
        
class NumberedOptionView(OptionView):
    def __init__(self, opt):
        OptionView.__init__(self, opt, rescan=True)

class IntegerOptionView(OptionView):
    def __init__(self, opt, **kwargs):
        OptionView.__init__(self, opt, **kwargs)
        self.sb = QSpinBox()
        self.sb.valueChanged.connect(self.changed)
        self.widgets.append(self.sb)
    
    def changed(self):
        val = str(self.sb.value())
        self.fab.set_option(self.key, val)
        
    def do_update(self):
        OptionView.do_update(self)
        if self.fab.options.has_key(self.key):
            self.sb.setValue(int(self.fab.options[self.key]))

    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addWidget(self.sb, row, 1)

class StringOptionView(OptionView):
    def __init__(self, opt, **kwargs):
        OptionView.__init__(self, opt, **kwargs)
        self.edit = QLineEdit()
        self.edit.editingFinished.connect(self.changed)
        self.widgets.append(self.edit)
        
    def changed(self):
        text = self.edit.text()
        self.fab.set_option(self.key, text)

    def do_update(self):
        OptionView.do_update(self)
        if self.fab.options.has_key(self.key):
            self.edit.setText(self.fab.options[self.key])
        
    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addWidget(self.edit, row, 1)
        
class FileOptionView(StringOptionView):
    def __init__(self, opt, **kwargs):
        StringOptionView.__init__(self, opt, **kwargs)
        self.hbox = QHBoxLayout()
        self.hbox.addWidget(self.edit)
        self.btn = QPushButton("Choose")
        self.hbox.addWidget(self.btn)
        self.widgets.append(self.btn)
        self.btn.clicked.connect(self.choose_file)
    
    def choose_file(self):
        fname = QFileDialog.getOpenFileName()[0]
        if fname: 
            self.edit.setText(fname)
            self.changed()

    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addLayout(self.hbox, row, 1)

class MatrixFileOptionView(FileOptionView):
    def __init__(self, opt, **kwargs):
        FileOptionView.__init__(self, opt, **kwargs)
        self.editBtn = QPushButton("Edit")
        self.hbox.addWidget(self.editBtn)
        self.widgets.append(self.editBtn)
        self.editBtn.clicked.connect(self.edit_file)
    
    def read_vest(self, fname):
        f = None
        try:
            f = open(fname)
            lines = f.readlines()
            nx, ny = 0, 0
            in_matrix = False
            mat = []
            for line in lines:
                if in_matrix:
                    nums = [float(num) for num in line.split()]
                    if len(nums) != nx: raise Exception ("Incorrect number of x values")
                    mat.append(nums)
                elif line.startswith("/Matrix"):
                  if nx == 0 or ny == 0: raise Exception("Missing /NumWaves or /NumPoints")
                  in_matrix = True
                elif line.startswith("/NumWaves"):
                  parts = line.split()
                  if len(parts) == 1: raise Exception("No number following /NumWaves")
                  nx = int(parts[1])
                elif line.startswith("/NumPoints") or line.startswith("/NumContrasts"):
                  parts = line.split()
                  if len(parts) == 1: raise Exception("No number following /NumPoints")
                  ny = int(parts[1])
            if len(mat) != ny:
                raise Exception("Incorrect number of y values")      
        finally:
            if f is not None: f.close()
        return mat

    def edit_file(self):
        print(self.edit.text())
        print(self.read_vest(self.edit.text()))
        
class OptionsView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        self.views = {}
        self.rescan = False
        self.btn.clicked.connect(self.show)
        
    def show(self):
        self.dialog.show()
        self.dialog._raise()

    def del_layout(self, layout):
        while True:
            w = layout.takeAt(0)
            if w is None:
                break
            elif w.widget() is None:
                self.del_layout(w)    
            else:
                w.widget().deleteLater()
                
    def clear(self):
        self.del_layout(self.dialog.grid)
        self.views = {}
        
    def get_concrete_opts(self, base, suffix):
        vals = [0]
        concrete = []
        for opt in self.fab.options.keys():
            if opt.startswith(base) and opt.endswith(suffix):
                try:
                    vals.append(int(opt[:len(opt) - len(suffix)][len(base):]))
                    concrete.append(opt)
                except:
                    traceback.print_exc()
                    print("WARNING: failed to parse %s" % opt)
        for n in range(max(vals)+1):
            if n+1 not in vals:
                break
        concrete.append(base + str(n+1) + suffix)
        return n+1, concrete

    def add_opts(self, opts, startrow):
        row = 0
        for opt in opts:
            if opt[0].find("<n>") >= 0:
                opt_base=opt[0][:opt[0].find("<n>")]
                opt_suffix = opt[0][opt[0].find("<n>") + 3:]
                next, actual = self.get_concrete_opts(opt_base, opt_suffix)
                for key in actual:
                    newopt = [key,] + opt[1:]
                    view = get_option_view(newopt, rescan=True)
                    view.add(self.dialog.grid, row+startrow)
                    self.views[key] = view
                    row += 1
                
            else:
                view = get_option_view(opt)
                view.add(self.dialog.grid, row+startrow)
                self.views[opt[0]] = view
                row += 1
        return row

    def do_update(self):
        if self.fab.changed("exec") or self.rescan:
            self.clear()
            self.opts, d = self.fab.get_options()
            if len(self.opts) == 0:
                msgBox = QMessageBox()
                msgBox.setText("Could not get options from Fabber - check the path to the executable")
                msgBox.exec_()
            self.title = "Fabber General Options"
            self.desc = "These options are not specific to a particular model or inference method"
            self.create_views()
            self.rescan = False

        for view in self.views.values():
            view.update(self.fab)

    def create_views(self):
        req = [opt for opt in self.opts if opt[2] == "REQUIRED"]
        nonreq = [opt for opt in self.opts if opt[2] != "REQUIRED"]
        
        self.dialog.modelLabel.setText(self.title)
        self.dialog.descLabel.setText(self.desc)
        
        if req:
            label = get_label("Required options", size=12, bold=True)
            self.dialog.grid.addWidget(label, 0, 0)
            self.add_opts(req, 1)
        if nonreq:
            label = get_label("Other options", size=12, bold=True)
            self.dialog.grid.addWidget(label, len(req)+1, 0)
            self.add_opts(nonreq, len(req)+2)

        self.dialog.adjustSize()
        
class ComponentOptionsView(OptionsView):
    def __init__(self, otype, text, **kwargs):
        OptionsView.__init__(self, **kwargs)
        self.type = otype
        self.text = text
        self.value = ""
        
    def do_update(self):
        value = self.fab.options.get(self.type,"")
        if self.fab.changed("exec") or self.value != value or self.rescan:
            self.value = value
            self.clear()
            if self.value != "":
                args = {self.type : self.value}
                self.opts, self.desc = self.fab.get_options(**args)
                self.title = "%s: %s" % (self.text, self.value)
                self.create_views()
            self.rescan = False

        for view in self.views.values():
            view.update(self.fab)
    
class InputDataView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        self.table.setHorizontalHeaderItem(0, QTableWidgetItem("Input Data"))    
        self.table.currentItemChanged.connect(self.data_select)
        self.dataBtn.clicked.connect(self.data_clicked)
        self.maskBtn.clicked.connect(self.mask_clicked)
        
    def data_clicked(self):
        filename = QFileDialog.getOpenFileName()[0]
        # FIXME data mode
        if filename: self.fab.set_option("data", filename)
        
    def mask_clicked(self):
        filename = QFileDialog.getOpenFileName()[0]
        if filename: self.fab.set_option("mask", filename)
        
    def data_select(self, cur, prev):
        self.fab.select(cur.text())

    def do_update(self):
        self.table.setRowCount(0)
        n = 0
        for key, data in self.fab.data.items():
            self.table.setRowCount(n+1)
            item = QTableWidgetItem(key)
            if data.visible:
                item.setIcon(self.table.style().standardIcon(QtGui.QStyle.SP_FileIcon))        
            self.table.setItem(n, 0, item)
            n += 1
            
class OutputDataView(ModelView):
    def __init__(self, logDialog, runDialog, **kwargs):
        ModelView.__init__(self, **kwargs)
        self.logDialog = logDialog
        self.runDialog = runDialog
        self.table.setHorizontalHeaderItem(0, QTableWidgetItem("Output Data"))    
        self.chooseRunBtn.clicked.connect(self.choose_run)
        self.logBtn.clicked.connect(self.view_log_pressed)
        self.table.currentItemChanged.connect(self.data_select)
        self.cur_timestamp = ""

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
        self.fab.select(cur.text())
        
    def relist(self):
        self.table.setRowCount(0)
        if self.fab.current_run:
            data = self.fab.current_run.data
            n = 0
            self.table.setRowCount(0)
            for key, data in data.items():
                self.table.setRowCount(n+1)
                item = QTableWidgetItem(key)
                if data.visible:
                    item.setIcon(self.table.style().standardIcon(QtGui.QStyle.SP_FileIcon))        
                self.table.setItem(n, 0, item)
                n += 1
            self.logDialog.textBrowser.setText(self.fab.current_run.get_log())
            if self.fab.current_run.isquick:
                msgBox = QMessageBox()
                msgBox.setText("This was a quick 1-voxel run, not a full run")
                msgBox.exec_()

    def do_update(self):
        if self.fab.current_run and self.cur_timestamp != self.fab.current_run.timestamp:
            self.cur_timestamp = self.fab.current_run.timestamp
            self.dateLabel.setText(self.fab.current_run.timestamp)
            self.relist()
            
    
class CurrentDataView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
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
            
        self.set_enabled(self.fab.current_data is not None)
        
        if self.fab.current_data:
            self.currentEdit.setText(self.fab.current_data.filename)
            self.slider.setValue(100*self.fab.current_data.alpha)
            if self.fab.current_data.visible:
                self.cb.setCheckState(QtCore.Qt.CheckState.Checked)
            else:
                self.cb.setCheckState(QtCore.Qt.CheckState.Unchecked)
            self.cmCombo.setCurrentIndex(self.maps.index(self.fab.current_data.cm))
            self.valueEdit.setText(str(self.fab.current_data.get_value(self.fab.focus)))
        
    def cm_changed(self, idx):
        text = self.maps[idx]
        self.fab.set_visibility(self.fab.current_data, cm=text)

    def alpha_changed(self, value):
        self.fab.set_visibility(self.fab.current_data, alpha=float(value)/100)

    def visible_changed(self):
        self.fab.set_visibility(self.fab.current_data, visible=self.cb.checkState() == QtCore.Qt.CheckState.Checked)

class FocusView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
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
        
class ExecView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        self.execChangeBtn.clicked.connect(self.exec_choose_file)
        
    def do_update(self):
        if self.fab.changed("exec"): 
            self.execEdit.setText(self.fab.options["fabber"])
    
    def exec_choose_file(self):
        fname = QFileDialog.getOpenFileName()[0]
        if fname:
            self.execEdit.setText(fname)
            self.fab.set_exec(self.execEdit.text())
        
class FileView(ModelView):
    def __init__(self, **kwargs):
        ModelView.__init__(self, **kwargs)
        
        self.runBtn.clicked.connect(self.run)
        self.runQuickBtn.clicked.connect(self.run_quick)
        self.saveBtn.clicked.connect(self.save)
        self.saveAsBtn.clicked.connect(self.save_as)
        self.changed = None
        
    def save(self):
        self.fab.save()
    
    def save_as(self):
        # fixme choose file name
        # fixme overwrite
        # fixme clone data
        fname = QFileDialog.getSaveFileName()[0]
        self.fab.set_filename(fname)
        self.fab.save()
        
    def run(self, quick=False):
        try:
            self.fab.run(quick)
        except:
            print sys.exc_info()
            QtGui.QMessageBox.warning(None, "Fabber error", str(sys.exc_info()[1]))
            
    def run_quick(self):
        self.run(quick=True)
    
    def do_update(self):
        self.edit.setText(self.fab.fabfile)
        if self.changed is None:
            self.changed = False
        elif self.fab.changed("options", "data", "exec"):
            self.changed = True
        
        self.saveBtn.setEnabled(self.changed)
        self.runBtn.setEnabled(self.fab.data.has_key("data"))
        self.runQuickBtn.setEnabled(self.fab.data.has_key("data"))
        
        
