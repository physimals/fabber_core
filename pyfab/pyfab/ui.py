from PySide import QtCore, QtGui
from PySide.QtGui import QMainWindow, QDialog

from ui_qtd import Ui_LogViewer, Ui_MainWindow, Ui_ModelOptionsDialog, Ui_ChooseRunDialog

from pyfab.views import *
from pyfab.plotviews import OrthoView, FitView
from pyfab.model import RunFile

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

class MainWindow(QMainWindow, Ui_MainWindow):
	def __init__(self, parent=None):
		QMainWindow.__init__(self, parent)
		self.setupUi(self)
		
		self.views = [
		  OrthoView(self.orthoFrame),
		  FitView(self.fitFrame),
		  ParamValuesView(table=self.paramTable),
          ExtraOptionsView(modelCombo=self.modelCombo, methodCombo=self.methodCombo),
		  InputDataView(table=self.inputFilesTable, combo=self.dataModeCombo, dataBtn=self.dataBtn, maskBtn=self.maskBtn),
		  OutputDataView(LogViewerDialog(), ChooseRunDialog(), table=self.outputFilesTable, dateLabel=self.dateLabel, logBtn=self.showLogBtn, chooseRunBtn=self.chooseRunBtn),
		  CurrentDataView(slider=self.alphaSlider, cb=self.visibleCB, cmCombo=self.cmCombo, valueEdit=self.dataValueEdit, currentEdit=self.currentEdit),
		  FocusView(slider=self.timeSlider, sbX=self.sbX, sbY=self.sbY, sbZ=self.sbZ, sbT=self.sbT),
		  ComponentOptionsView("model", "Forward model", dialog=ModelOptionsDialog(), btn=self.modelOptionsBtn),
		  ComponentOptionsView("method", "Inference method", dialog=ModelOptionsDialog(), btn=self.methodOptionsBtn),
		  OptionsView(dialog=ModelOptionsDialog(), btn=self.generalOptionsBtn),	
		  FileView(runBtn=self.runBtn, runQuickBtn=self.runQuickBtn, saveBtn=self.saveBtn, saveAsBtn=self.saveAsBtn, edit=self.filenameEdit),
		  ExecView(execChangeBtn=self.execChangeBtn, execEdit=self.execEdit),
		]
		
		self.openBtn.clicked.connect(self.open)
		self.newBtn.clicked.connect(self.new)

	def add_views(self):	
		for view in self.views: self.fab.add_view(view)
		
	def refresh(self, fab):
		self.fab = fab
		self.add_views()
	
	def new(self):
		# fixme unsaved changes
		self.refresh(RunFile())
		
	def open(self):
		filename = QFileDialog.getOpenFileName()[0]
		if filename: self.refresh(RunFile(filename))

