import sys, os
import traceback
import re

from PySide import QtCore, QtGui

from .mvc import View
from fabber import FabberLib, find_fabber

NUMBERED_OPTIONS_MAX=20

class ModelMethodView(View):
    def __init__(self, **kwargs):
        View.__init__(self, ["fabber", "model", "method", "loadmodels"], **kwargs)
        self.models = None
        self.methods = None
        self.modelCombo.currentIndexChanged.connect(self.model_changed)
        self.methodCombo.currentIndexChanged.connect(self.method_changed)
        self.auto_load_models = kwargs.get("auto_load_models", False)
        
    def model_changed(self):
        self.rundata["model"] = self.modelCombo.currentText()
        
    def method_changed(self):
        self.rundata["method"] = self.methodCombo.currentText()

    def do_update(self):
        if self.rundata.changed("fabber", "loadmodels"):
            self.modelCombo.clear()
            self.models = FabberLib(rundata=self.rundata, auto_load_models=self.auto_load_models).get_models()
            for model in self.models:
                self.modelCombo.addItem(model)
        
            self.methodCombo.clear()
            self.methods = FabberLib(rundata=self.rundata).get_methods()
            for method in self.methods:
                self.methodCombo.addItem(method)
                
        if "model" in self.rundata:
            cmodel = self.rundata["model"]
            for idx, model in enumerate(self.models):
                if cmodel == model:
                    self.modelCombo.setCurrentIndex(idx)

        if "method" in self.rundata:
            cmethod = self.rundata["method"]
            for idx, method in enumerate(self.methods):
                if cmethod == method:
                    self.methodCombo.setCurrentIndex(idx)

        
def get_label(text, size=None, bold=False, italic=False):
    label = QtGui.QLabel(text)
    font = label.font()
    font.setBold(bold)
    font.setItalic(italic)
    if size: font.setPointSize(size)
    label.setFont(font)
    return label

class OptionView(View):
    def __init__(self, opt, **kwargs):
        View.__init__(self, [opt["name"],], **kwargs)
        self.key = opt["name"]
        self.dtype = opt["type"]
        self.req = not opt["optional"]
        self.default = opt["default"]
        self.desc = opt["description"]
        self.dependents = []
        self.desc_first = kwargs.get("desc_first", False)
        self.label = get_label(opt["name"], size=10)
        self.desclabel = get_label(opt["description"], size=10)
        self.desclabel.setToolTip("--%s" % self.key)
        self.desclabel.resize(400, self.desclabel.height())
        self.desclabel.setWordWrap(True)
        if self.req:
            self.enable_cb = None
        else:
            self.enable_cb = QtGui.QCheckBox()
            self.enable_cb.stateChanged.connect(self.state_changed)
            self.widgets.append(self.enable_cb)
        self.widgets.append(self.label)
        self.widgets.append(self.desclabel)

    def add_dependent(self, dep):
        if not self.enable_cb: return
        self.dependents.append(dep)
        checked = self.enable_cb.checkState() == QtCore.Qt.CheckState.Checked
        dep.set_visible(checked)
        if not checked: dep.enable_cb.setChecked(False)

    def set_visible(self, visible=True, widgets=None):
        if widgets is None: widgets = self.widgets
        for widget in widgets:
                widget.setVisible(visible)

    def state_changed(self):
        # This function is only called if we have a checkbox
        checked = self.enable_cb.checkState() == QtCore.Qt.CheckState.Checked
        self.set_enabled(checked)
        self.enable_cb.setEnabled(True)
        
        for dep in self.dependents:
            dep.set_visible(checked)
            if not checked: dep.enable_cb.setChecked(False)

        if checked:
            self.changed()
        else:
            del self.rundata[self.key]
        
    def changed(self):
        self.rundata[self.key] = ""
        
    def do_update(self):
        if not self.req:
            if not self.key in self.rundata:
                self.enable_cb.setCheckState(QtCore.Qt.CheckState.Unchecked)
            else:
                self.enable_cb.setCheckState(QtCore.Qt.CheckState.Checked)

            self.set_enabled(self.enable_cb.checkState() == QtCore.Qt.CheckState.Checked)
            self.enable_cb.setEnabled(True)
       
    def add(self, grid, row):
        if self.desc_first:
            label = self.desclabel
        else:
            label = self.label
            grid.addWidget(self.desclabel, row, 2)

        if self.req:
            grid.addWidget(label, row, 0)
        else:
            hbox = QtGui.QHBoxLayout()
            hbox.addWidget(self.enable_cb)
            hbox.addWidget(label, row)
            grid.addLayout(hbox, row, 0)

class IntegerOptionView(OptionView):
    def __init__(self, opt, **kwargs):
        OptionView.__init__(self, opt, **kwargs)
        self.sb = QtGui.QSpinBox()
        self.sb.valueChanged.connect(self.changed)
        self.widgets.append(self.sb)
    
    def changed(self):
        val = str(self.sb.value())
        self.rundata[self.key] = val
        
    def do_update(self):
        OptionView.do_update(self)
        if self.key in self.rundata:
            self.sb.setValue(int(self.rundata[self.key]))

    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addWidget(self.sb, row, 1)

class StringOptionView(OptionView):
    def __init__(self, opt, **kwargs):
        OptionView.__init__(self, opt, **kwargs)
        self.edit = QtGui.QLineEdit()
        self.edit.editingFinished.connect(self.changed)
        self.widgets.append(self.edit)
        
    def changed(self):
        # Note that this signal is triggered when the widget
        # is enabled/disabled!
        if self.edit.isEnabled():
            self.rundata[self.key] = self.edit.text()

    def do_update(self):
        OptionView.do_update(self)
        if self.key in self.rundata:
            self.text = self.rundata[self.key]
            self.edit.setText(self.text)
        
    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addWidget(self.edit, row, 1)
        
class FileOptionView(StringOptionView):
    def __init__(self, opt, **kwargs):
        StringOptionView.__init__(self, opt, **kwargs)
        self.hbox = QtGui.QHBoxLayout()
        self.hbox.addWidget(self.edit)
        self.btn = QtGui.QPushButton("Choose")
        self.hbox.addWidget(self.btn)
        self.widgets.append(self.btn)
        self.btn.clicked.connect(self.choose_file)
    
    def choose_file(self):
        dialog = QtGui.QFileDialog()
        dialog.setFileMode(QtGui.QFileDialog.AnyFile)
        if dialog.exec_():
            fname = dialog.selectedFiles()[0]
            self.edit.setText(fname)
            self.changed()

    def add(self, grid, row):
        OptionView.add(self, grid, row)
        grid.addLayout(self.hbox, row, 1)

class MatrixFileOptionView(FileOptionView):
    def __init__(self, opt, **kwargs):
        FileOptionView.__init__(self, opt, **kwargs)
        self.editBtn = QtGui.QPushButton("Edit")
        self.hbox.addWidget(self.editBtn)
        self.widgets.append(self.editBtn)
        self.editBtn.clicked.connect(self.edit_file)
    
    def read_vest(self, fname):
        f = None
        in_matrix = False
        mat = []
        try:
            f = open(fname, "r")
            lines = f.readlines()
            nx, ny = 0, 0
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

        return mat, ""

    def read_ascii(self, fname):
        f = None
        in_matrix = False
        mat = []
        desc = ""
        try:
            f = open(fname, "r")
            lines = f.readlines()
            nx = 0
            for line in lines:
                if line.strip().startswith("#"):
                    desc += line.lstrip("#")
                else:
                    row = [float(n) for n in line.split()]
                    if not in_matrix:
                        nx = len(row)
                        in_matrix = True
                    elif len(row) != nx:
                        raise Exception("Incorrect number of x values: %s" % line)
                    mat.append(row)
        finally:
            if f is not None: f.close()

        return mat, desc

    def write_vest(self, fname, m, desc=""):
        f = None
        try:
            f = open(fname, "w")
            nx, ny = len(m), len(m[0])
            f.write("/NumWaves %i\n" % nx)
            f.write("/NumPoints %i\n" % ny)
            f.write("/Matrix\n")
            for row in m:
                for item in row:
                    f.write("%f " % item)
                f.write("\n")
        finally:
            if f is not None: f.close()

    def write_ascii(self, fname, m, desc=""):
        f = None
        try:
            f = open(fname, "w")
            for line in desc.splitlines():
                f.write("#%s\n" % line)
            for row in m:
                f.write(" ".join([str(v) for v in row]))
                f.write("\n")
        finally:
            if f is not None: f.close()

    def edit_file(self):
        fname = self.edit.text()
        print("editing: ", fname)
        if fname.strip() == "":
            msgBox = QtGui.QMessageBox()
            msgBox.setText("Enter a filename")
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok)
            msgBox.exec_()
            return
        elif not os.path.exists(fname):
            msgBox = QtGui.QMessageBox()
            msgBox.setText("File does not exist - create?")
            msgBox.setStandardButtons(QtGui.QMessageBox.Ok | QtGui.QMessageBox.Cancel)
            msgBox.setDefaultButton(QtGui.QMessageBox.Ok)
            ret = msgBox.exec_()
            if ret != QtGui.QMessageBox.Ok:
                return
            open(fname, "a").close()

        try:
            try:
                mat, desc = self.read_vest(fname)
                ascii = False
            except:
                mat, desc = self.read_ascii(fname)
                ascii = True
            self.mat_dialog.set_matrix(mat, desc)
            if self.mat_dialog.exec_():
                mat, desc = self.mat_dialog.get_matrix()
                print(mat, desc)
                if ascii:
                    self.write_ascii(fname, mat, desc)
                else:
                    self.write_vest(fname, mat, desc)
        except:
            traceback.print_exc()

OPT_VIEW = {
    "INT" : IntegerOptionView,
    "BOOL": OptionView,
    "FILE": FileOptionView,
    "IMAGE": FileOptionView,
    "MVN": FileOptionView,
    "MATRIX": MatrixFileOptionView,
    "": StringOptionView,
}

def get_option_view(opt, **kwargs):
    if opt["type"] in OPT_VIEW:
        return OPT_VIEW[opt["type"]](opt, **kwargs)
    else:
        return OPT_VIEW[""](opt, **kwargs)

class OptionsView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [], **kwargs)
        self.views = {}
        self.ignore_opts = set()
        self.btn.clicked.connect(self.show)
        self.desc_first = kwargs.get("desc_first", False)

    def show(self):
        self.dialog.show()
        self.dialog.raise_()

    def ignore(self, *opts):
        for opt in opts:
            self.ignore_opts.add(opt)

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

    def add_opts(self, opts, startrow):
        row = 0
        for opt in opts:
            if opt["name"].find("<n>") >= 0:
                # This is a numbered option. Create one for each
                opt_base=opt["name"][:opt["name"].find("<n>")]
                opt_suffix = opt["name"][opt["name"].find("<n>") + 3:]
                for n in range(1, NUMBERED_OPTIONS_MAX+1):
                    newopt = dict(opt)
                    newopt["name"] = "%s%i%s" % (opt_base, n, opt_suffix)
                    view = get_option_view(newopt, desc_first=self.desc_first)
                    view.mat_dialog = self.mat_dialog
                    if n > 1:
                        prev.add_dependent(view)

                    view.add(self.dialog.grid, row + startrow)
                    self.views[newopt["name"]] = view
                    prev = view
                    row += 1
            else:
                view = get_option_view(opt, desc_first=self.desc_first)
                view.mat_dialog = self.mat_dialog
                view.add(self.dialog.grid, row+startrow)
                self.views[opt["name"]] = view
                row += 1
        return row

    def do_update(self):
        if self.rundata.changed("fabber"):
            self.clear()
            self.opts, self.desc = FabberLib(rundata=self.rundata).get_options()
            self.opts = [opt for opt in self.opts if opt["name"] not in self.ignore_opts]
            if len(self.opts) == 0:
                msgBox = QtGui.QMessageBox()
                msgBox.setText("Could not get options from Fabber - check the path to the executable")
                msgBox.exec_()
            self.title = "Fabber General Options"
            self.desc = "These options are not specific to a particular model or inference method"
            self.create_views()

        for view in self.views.values():
            view.update(self.rundata)

    def create_views(self):
        req = [opt for opt in self.opts if not opt["optional"]]
        nonreq = [opt for opt in self.opts if opt["optional"]]
        
        self.dialog.modelLabel.setText(self.title)
        self.dialog.descLabel.setText(self.desc)
        
        if req:
            label = get_label("Mandatory options", size=12, bold=True)
            self.dialog.grid.addWidget(label, 0, 0)
            self.add_opts(req, 1)
        if nonreq:
            label = get_label("Non-mandatory options", size=12, bold=True)
            self.dialog.grid.addWidget(label, len(req)+1, 0)
            self.add_opts(nonreq, len(req)+2)

        self.dialog.grid.setAlignment(QtCore.Qt.AlignTop)
        self.dialog.adjustSize()
        
class ComponentOptionsView(OptionsView):
    """
    Options dialog for model or method
    """
    def __init__(self, otype, text, **kwargs):
        OptionsView.__init__(self, **kwargs)
        self.type = otype
        self.text = text
        self.value = ""
        self.auto_load_models = kwargs.get("auto_load_models", False)

    def do_update(self):
        value = self.rundata.get(self.type,"")
        if self.rundata.changed("fabber") or self.value != value:
            self.value = value
            self.clear()
            if self.value != "":
                args = {self.type : self.value}
                self.opts, self.desc = FabberLib(rundata=self.rundata, auto_load_models=self.auto_load_models).get_options(**args)
                self.opts = [opt for opt in self.opts if opt["name"] not in self.ignore_opts]
                self.title = "%s: %s" % (self.text, self.value)
                self.create_views()

        for view in self.views.values():
            view.update(self.rundata)
   
class ModelOptionsView(ComponentOptionsView):
    """
    Options dialog for model
    """
    def __init__(self, **kwargs):
        ComponentOptionsView.__init__(self, "model", "Forward Model", **kwargs)

    def do_update(self):
        ComponentOptionsView.do_update(self)
        self.btn.setText("%s model options" % self.value.upper())

class MethodOptionsView(ComponentOptionsView):
    """
    Options dialog for inference method
    """
    def __init__(self, **kwargs):
        ComponentOptionsView.__init__(self, "method", "Inference Method", **kwargs)

    def do_update(self):
        ComponentOptionsView.do_update(self)
        self.btn.setText("%s method options" % self.value.upper())
   
class ChooseFileView(View):
    def __init__(self, opt, **kwargs):
        self.defaultDir = ""
        self.dialogTitle = "Choose a file"
        self.opt = opt
        View.__init__(self, [opt,], **kwargs)
        self.changeBtn.clicked.connect(self.choose_file)
        
    def do_update(self):
        if self.opt in self.rundata:
            self.edit.setText(self.rundata[self.opt])

    def choose_file(self):
        fname = QtGui.QFileDialog.getOpenFileName(None, self.dialogTitle, self.defaultDir)[0]
        if fname:
            self.edit.setText(fname)
            self.rundata[self.opt] = fname
            self.defaultDir = os.path.dirname(fname)

class ChooseModelLib(View):
    def __init__(self, **kwargs):
        self.opt = "loadmodels"
        View.__init__(self, [self.opt, ], **kwargs)
        ex, lib, model_libs = find_fabber()
        self.defaultDir = os.path.dirname(lib)
        for model_lib in model_libs:
            self.combo.addItem(self.model_name(model_lib), model_lib)
        self.combo.addItem("<No additional models>", "")
        self.changeBtn.clicked.connect(self.choose_file)
        self.combo.currentIndexChanged.connect(self.lib_changed)
        self.combo.setCurrentIndex(-1)

    def model_name(self, model_lib):
        match = re.match(".*fabber_models_(.+)\..+", model_lib, re.I)
        if match:
            return match.group(1);
        else:
            return model_lib

    def select_lib(self, lib):
        if self.combo.findData(lib) < 0:
            self.combo.addItem(self.model_name(lib), lib)
        if self.combo.itemData(self.combo.currentIndex()) != lib:
            self.combo.setCurrentIndex(self.combo.findData(lib))

    def choose_file(self):
        fname = QtGui.QFileDialog.getOpenFileName(None, "Choose model library", self.defaultDir)[0]
        if fname:
            self.defaultDir = os.path.dirname(fname)
            self.select_lib(fname)

    def lib_changed(self, idx):
        if idx >= 0:
            lib = self.combo.itemData(self.combo.currentIndex())
            if lib != "":
                self.rundata[self.opt] = lib
            else:
                del self.rundata[self.opt]
            self.combo.setToolTip(lib)

    def do_update(self):
        self.select_lib(self.rundata.get(self.opt, ""))

class FileView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [], **kwargs)
        
        if self.runBtn: self.runBtn.clicked.connect(self.run)
        if self.runQuickBtn: self.runQuickBtn.clicked.connect(self.run_quick)
        self.saveBtn.clicked.connect(self.save)
        self.saveAsBtn.clicked.connect(self.save_as)
        self.changed = None
        
    def save(self):
        self.rundata.save()
    
    def save_as(self):
        # fixme choose file name
        # fixme overwrite
        # fixme clone data
        fname = QtGui.QFileDialog.getSaveFileName()[0]
        self.rundata.set_file(fname)
        self.rundata.save()
        
    def run(self, focus=None):
        try:
            self.rundata.run(focus=focus)
        except:
            print sys.exc_info()
            QtGui.QMessageBox.warning(None, "Fabber error", str(sys.exc_info()[1]))
            
    def run_quick(self):
        self.run()
    
    def do_update(self):
        self.edit.setText(self.rundata.get_filename())
        if self.changed is None:
            self.changed = False
        else:
            self.changed = True
        
        self.saveBtn.setEnabled(self.changed)
        self.runBtn.setEnabled("data" in self.rundata)
        self.runQuickBtn.setEnabled("data" in self.rundata)
        
        
