from PySide import QtCore, QtGui
from PySide.QtGui import QMainWindow, QDialog

from ui_qtd import Ui_LogViewer, Ui_MainWindow, Ui_ModelOptionsDialog, Ui_ChooseRunDialog, Ui_VestDialog

from pyfab.views import *
from pyfab.imagedata import *
from pyfab.plotviews import OrthoView, FitView
from pyfab.model import FabberRunData

class LogViewerDialog(QDialog, Ui_LogViewer):
	def __init__(self, log="", parent=None):
		QDialog.__init__(self, parent)
		self.setupUi(self)
		self.textBrowser.setText(log)

class ModelOptionsDialog(QDialog, Ui_ModelOptionsDialog):
    def __init__(self, parent=None):
        super(ModelOptionsDialog, self).__init__(parent)
        self.setupUi(self)

class ChooseRunDialog(QDialog, Ui_ChooseRunDialog):
    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.setupUi(self)

class MatrixEditDialog(QDialog, Ui_VestDialog):
    def __init__(self, parent=None):
        QDialog.__init__(self, parent)
        self.setupUi(self)
        self.table.setRowCount(1)
        self.table.setColumnCount(1)
        self.AddColBtn.clicked.connect(self.add_col)
        self.DelColBtn.clicked.connect(self.del_col)
        self.AddRowBtn.clicked.connect(self.add_row)
        self.DelRowBtn.clicked.connect(self.del_row)

    def add_col(self):
        self.table.insertColumn(self.table.currentColumn())

    def del_col(self):
        if self.table.columnCount() > 1:
            self.table.removeColumn(self.table.currentColumn())

    def add_row(self):
        self.table.insertRow(self.table.currentRow())

    def del_row(self):
        if self.table.rowCount() > 1:
            self.table.removeRow(self.table.currentRow())

    def set_matrix(self, m):
        self.table.setRowCount(max(1, len(m)))
        if len(m) > 0: self.table.setColumnCount(max(1, len(m[0])))
        for x, row in enumerate(m):
            for y, d in enumerate(row):
                self.table.setItem(x, y, QTableWidgetItem(str(d)))

    def get_matrix(self):
        m = []
        for x in range(self.table.columnCount()):
            row = []
            for y in range(self.table.rowCount()):
                try:
                    data = self.table.item(x, y).text()
                    row.append(float(data))
                except:
                    row.append(0)
                    print("WARNING: non-numeric data, converting to zero")
            m.append(row)
        return m                 
            
class MainWindow(QMainWindow, Ui_MainWindow):
	def __init__(self, parent=None):
		QMainWindow.__init__(self, parent)
		self.setupUi(self)
		
		self.imdata = FabberImageData()
		
		self.imviews = [
		  OrthoView(self.orthoFrame),
		  FitView(self.fitFrame),
		  FocusView(slider=self.timeSlider, sbX=self.sbX, sbY=self.sbY, sbZ=self.sbZ, sbT=self.sbT),
                  ParamValuesView(table=self.paramTable),
		  DataView(LogViewerDialog(), ChooseRunDialog(), table=self.filesTable, dateLabel=self.dateLabel, logBtn=self.showLogBtn, chooseRunBtn=self.chooseRunBtn, dataBtn=self.dataBtn, maskBtn=self.maskBtn),
		  CurrentDataView(slider=self.alphaSlider, cb=self.visibleCB, cmCombo=self.cmCombo, valueEdit=self.dataValueEdit, currentEdit=self.currentEdit),
		]

                self.views = [
		  self.imdata,
                  ModelMethodView(modelCombo=self.modelCombo, methodCombo=self.methodCombo),
		  ComponentOptionsView("model", "Forward model", dialog=ModelOptionsDialog(), btn=self.modelOptionsBtn, mat_dialog=MatrixEditDialog()),
		  ComponentOptionsView("method", "Inference method", dialog=ModelOptionsDialog(), btn=self.methodOptionsBtn, mat_dialog=MatrixEditDialog()),
		  OptionsView(dialog=ModelOptionsDialog(), btn=self.generalOptionsBtn, mat_dialog=MatrixEditDialog()),	
		  FileView(runBtn=self.runBtn, runQuickBtn=self.runQuickBtn, saveBtn=self.saveBtn, saveAsBtn=self.saveAsBtn, edit=self.filenameEdit),
		  ChooseFileView("fabber", changeBtn=self.execChangeBtn, edit=self.execEdit),
		]
		
		self.openBtn.clicked.connect(self.open)
		self.newBtn.clicked.connect(self.new)

	def add_views(self):	
		for view in self.views: self.fab.add_view(view)
		for view in self.imviews: self.imdata.add_view(view)
		
	def refresh(self, fab):
		self.fab = fab
		self.add_views()
	
	def new(self):
		# fixme unsaved changes
		self.refresh(FabberRunData())
		
	def open(self):
		filename = QFileDialog.getOpenFileName()[0]
		if filename: self.refresh(FabberRunData(filename))

