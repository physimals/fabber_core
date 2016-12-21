import sys, os
import traceback

from PySide import QtCore, QtGui

from PySide.QtGui import QMessageBox, QLabel, QHBoxLayout, QLineEdit, QVBoxLayout, QFileDialog, QTableWidgetItem, QPlainTextEdit, QSpinBox, QCheckBox, QPushButton

from mvc import View
from model import FabberRunData
from fabber import FabberLib, FabberExec

class ModelMethodView(View):
    def __init__(self, **kwargs):
        View.__init__(self, ["fabber", "model", "method", "loadmodels"], **kwargs)
        self.models = None
        self.methods = None
        self.modelCombo.currentIndexChanged.connect(self.model_changed)
        self.methodCombo.currentIndexChanged.connect(self.method_changed)

    def model_changed(self):
        self.fab.set_option("model", self.modelCombo.currentText())
        
    def method_changed(self):
        self.fab.set_option("method", self.methodCombo.currentText())
        
    def do_update(self):
        if self.fab.changed("fabber", "loadmodels"):
            self.modelCombo.clear()
            self.models = FabberLib(self.fab).get_models()
            for model in self.models:
                self.modelCombo.addItem(model)
        
            self.methodCombo.clear()
            self.methods = FabberLib(self.fab).get_methods()
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

class OptionView(View):
    def __init__(self, opt, rescan=False, **kwargs):
        View.__init__(self, [opt["name"],], **kwargs)
        self.key = opt["name"]
        self.dtype = opt["type"]
        self.req = not opt["optional"]
        self.default = opt["default"]
        self.desc = opt["description"]
        self.rescan = rescan

        if self.req:
            self.label = get_label(opt["name"], size=10)
        else:
            self.label = QCheckBox(opt["name"])
            self.label.stateChanged.connect(self.state_changed)

        self.widgets.append(self.label)
        
        self.desclabel = get_label(opt["description"], size=10)
        self.desclabel.resize(400, self.desclabel.height())
        #desclabel.setWordWrap(True)
        self.widgets.append(self.desclabel)

    def state_changed(self):
        # Only called if label really is a checkbox!
        self.set_enabled(self.label.checkState() == QtCore.Qt.CheckState.Checked)
        self.label.setEnabled(True)
        
        if self.rescan: self.fab._change() # FIXME a hack
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

            self.set_enabled(self.label.checkState() == QtCore.Qt.CheckState.Checked)
            self.label.setEnabled(True)
       
    def add(self, grid, row):
        grid.addWidget(self.label, row, 0)
        grid.addWidget(self.desclabel, row, 2)

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
        # Note that this signal is triggered when the widget
        # is enabled/disabled!
        if self.edit.isEnabled():
            self.fab.set_option(self.key, self.edit.text())

    def do_update(self):
        OptionView.do_update(self)
        if self.fab.options.has_key(self.key):
            self.text = self.fab.options[self.key]
            self.edit.setText(self.text)
        
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
        dialog = QFileDialog(self)
        dialog.setFileMode(QFileDialog.AnyFile)
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
        self.editBtn = QPushButton("Edit")
        self.hbox.addWidget(self.editBtn)
        self.widgets.append(self.editBtn)
        self.editBtn.clicked.connect(self.edit_file)
    
    def read_vest(self, fname):
        f = None
        try:
            f = open(fname, "r")
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

    def write_vest(self, fname, m):
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

    def edit_file(self):
        fname = self.edit.text()
        if fname.strip() == "":
            msgBox = QMessageBox()
            msgBox.setText("Enter a filename")
            msgBox.setStandardButtons(QMessageBox.Ok)
            msgBox.exec_()
            return
        elif not os.path.exists(fname):
            msgBox = QMessageBox()
            msgBox.setText("File does not exist - create?")
            msgBox.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
            msgBox.setDefaultButton(QMessageBox.Ok)
            ret = msgBox.exec_()
            if ret != QMessageBox.Ok:
                return
            open(fname, "a").close()

        try:
            mat = self.read_vest(fname)
            self.mat_dialog.set_matrix(mat)
            if self.mat_dialog.exec_():
                print(self.mat_dialog.get_matrix())
                self.write_vest(fname, self.mat_dialog.get_matrix())
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
    if OPT_VIEW.has_key(opt["type"]):
        return OPT_VIEW[opt["type"]](opt, **kwargs)
    else:
        return OPT_VIEW[""](opt, **kwargs)

class OptionsView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [], **kwargs)
        self.views = {}
        self.ignore_opts = set()
        self.btn.clicked.connect(self.show)

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
        
    def get_concrete_opts(self, base, suffix):
        vals = [0]
        concrete = []
        for opt in self.fab.options.keys():
            if opt.startswith(base) and opt.endswith(suffix):
                try:
                    vals.append(int(opt[:len(opt) - len(suffix)][len(base):]))
                    concrete.append(opt)
                except:
                    # Can happen with empty suffixes, e.g. PSP_byname<n> and PSP_byname<n>_image
                    pass
        for n in range(max(vals)+1):
            if n+1 not in vals:
                break
        concrete.append(base + str(n+1) + suffix)
        return n+1, concrete

    def add_opts(self, opts, startrow):
        row = 0
        for opt in opts:
            if opt["name"].find("<n>") >= 0:
                # This is a numbered option
                opt_base=opt["name"][:opt["name"].find("<n>")]
                opt_suffix = opt["name"][opt["name"].find("<n>") + 3:]
                next, actual = self.get_concrete_opts(opt_base, opt_suffix)
                for key in actual:
                    newopt = dict(opt)
                    newopt["name"] = key
                    view = get_option_view(newopt, rescan=True, mat_dialog=self.mat_dialog)
                    view.add(self.dialog.grid, row+startrow)
                    self.views[key] = view
                    row += 1
                
            else:
                view = get_option_view(opt, mat_dialog=self.mat_dialog)
                view.add(self.dialog.grid, row+startrow)
                self.views[opt["name"]] = view
                row += 1
        return row

    def do_update(self):
        if self.fab.changed("fabber"):
            self.clear()
            self.opts, d = FabberLib(self.fab).get_options()
            self.opts = [opt for opt in self.opts if opt["name"] not in self.ignore_opts]
            if len(self.opts) == 0:
                msgBox = QMessageBox()
                msgBox.setText("Could not get options from Fabber - check the path to the executable")
                msgBox.exec_()
            self.title = "Fabber General Options"
            self.desc = "These options are not specific to a particular model or inference method"
            self.create_views()

        for view in self.views.values():
            view.update(self.fab)

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
        
    def do_update(self):
        value = self.fab.options.get(self.type,"")
        if self.fab.changed("fabber") or self.value != value:
            self.value = value
            self.clear()
            if self.value != "":
                args = {self.type : self.value}
                self.opts, self.desc = FabberLib(self.fab).get_options(**args)
                self.opts = [opt for opt in self.opts if opt["name"] not in self.ignore_opts]
                self.title = "%s: %s" % (self.text, self.value)
                self.create_views()

        for view in self.views.values():
            view.update(self.fab)
   
class ChooseFileView(View):
    def __init__(self, opt, **kwargs):
        View.__init__(self, [opt,], **kwargs)
        self.opt = opt
        self.changeBtn.clicked.connect(self.choose_file)
        
    def do_update(self):
        if self.fab.options.has_key(self.opt):
            self.edit.setText(self.fab.options[self.opt])

    def choose_file(self):
        fname = QFileDialog.getOpenFileName()[0]
        if fname:
            self.edit.setText(fname)
            self.fab.set_option(self.opt, fname)

class FileView(View):
    def __init__(self, **kwargs):
        View.__init__(self, [], **kwargs)
        
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
        self.fab.set_file(fname)
        self.fab.save()
        
    def run(self, focus=None):
        try:
            self.fab.run(focus=focus)
        except:
            print sys.exc_info()
            QtGui.QMessageBox.warning(None, "Fabber error", str(sys.exc_info()[1]))
            
    def run_quick(self):
        self.run()
    
    def do_update(self):
        self.edit.setText(self.fab.get_filename())
        if self.changed is None:
            self.changed = False
        else:
            self.changed = True
        
        self.saveBtn.setEnabled(self.changed)
        self.runBtn.setEnabled(self.fab.options.has_key("data"))
        self.runQuickBtn.setEnabled(self.fab.options.has_key("data"))
        
        
