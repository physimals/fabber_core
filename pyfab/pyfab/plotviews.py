import matplotlib
matplotlib.use('Qt4Agg')
import pylab

from matplotlib.backends.backend_qt4agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from PySide.QtGui import QHBoxLayout

from views import ModelView

class OrthoView(ModelView):
	def __init__(self, frame):
		ModelView.__init__(self)
		self.fig = Figure(figsize=(400, 400), facecolor='black')
		self.fig.subplots_adjust(0.01, 0.015, 0.99, 0.99, 0.01, 0.01)
		self.canvas = FigureCanvas(self.fig)
		layout = QHBoxLayout()
		layout.addWidget(self.canvas)
		frame.setLayout(layout)
		self.canvas.mpl_connect("button_press_event", self.onpick)
		
		self.axes = []
		self.axes.append(self.fig.add_subplot(221, axisbg='black'))
		self.axes.append(self.fig.add_subplot(222, axisbg='black'))
		self.axes.append(self.fig.add_subplot(223, axisbg='black'))
		self.tbox = self.fig.add_subplot(224, axisbg='black')
		self.ftext = self.tbox.text(0.2, 0.5, "", color='white')
		self.xhairs, self.yhairs = [0,0,0], [0,0,0]
		for i in range(3):
			self.yhairs[i] = self.axes[i].axvline(0.5, color='white', linestyle='-')
			self.xhairs[i] = self.axes[i].axhline(0.5, color='white', linestyle='-')
			
		self.maxdim = 1
		self.config_axes()

		self.imgs = [{}, {}, {}]

	def onpick(self, event):
		if self.fab is None: return
		x, y = int(event.xdata), int(event.ydata)
		if event.inaxes == self.axes[1]:
			self.fab.update_focus(yp=x, zp=y)
		if event.inaxes == self.axes[0]:
			self.fab.update_focus(xp=x, zp=y)
		if event.inaxes == self.axes[2]:
			self.fab.update_focus(xp=x, yp=y)
		
	def config_axes(self):
		for ax in self.axes:
			ax.get_xaxis().set_visible(False)
			ax.get_yaxis().set_visible(False)
			ax.set_xlim(0, self.maxdim)
			ax.set_ylim(0, self.maxdim)
		self.tbox.get_xaxis().set_visible(False)
		self.tbox.get_yaxis().set_visible(False)
		self.tbox.set_xlim(0, 1)
		self.tbox.set_ylim(0, 1)
		
	def set_extents(self):
		self.maxdim = max(self.fab.shape[:3])
		if self.maxdim > 0: self.config_axes()
		
	def add_slices(self, key, item, focus):
		slices = []
		slices.append(item.get_slice(focus, "y").T)
		slices.append(item.get_slice(focus, "x").T)
		slices.append(item.get_slice(focus, "z").T)
		for i, ax in enumerate(self.axes):
			if item.visible:
				if not self.imgs[i].has_key(key):
					ox = (self.maxdim - slices[i].shape[0])/2
					ex = ox + slices[i].shape[0]
					oy = (self.maxdim - slices[i].shape[1])/2
					ey = oy + slices[i].shape[1]
					self.imgs[i][key] = ax.imshow(slices[i], interpolation="nearest", vmin=item.min, vmax=item.max)
				else:
					self.imgs[i][key].set_data(slices[i])
					
				self.imgs[i][key].set_alpha(item.alpha)
				self.imgs[i][key].set_cmap(item.cm)
				
			elif self.imgs[i].has_key(key):
				# No longer visible
				self.imgs[i][key].remove()
				del self.imgs[i][key]

	def do_update(self):
		self.set_extents()
		focus = self.fab.focus
		if self.fab.current_run:
			outdata = self.fab.current_run.data
		else:
			outdata = {}

		if self.fab.changed("data","focus"):
			# Get rid of data no longer present
			for img in self.imgs:
				for key in img.keys():
					if not self.fab.data.has_key(key) and not outdata.has_key(key):
						img[key].remove()
						del img[key]
						
			for key, item in self.fab.data.items():
				self.add_slices(key, item, focus)
			for key, item in outdata.items():
				self.add_slices(key, item, focus)
						
		if self.fab.changed("focus") and max(self.fab.shape) > 0:
			for i in range(3):
				if i == 0:
					self.yhairs[i].set_xdata(focus[0])
					self.xhairs[i].set_ydata(focus[2])
				elif i == 1:
					self.yhairs[i].set_xdata(focus[1])
					self.xhairs[i].set_ydata(focus[2])
				else:
					self.yhairs[i].set_xdata(focus[0])
					self.xhairs[i].set_ydata(focus[1])

			self.ftext.set_text("X: %i\nY: %i\nZ: %i\nt: %i" % tuple(focus))
		self.canvas.draw_idle()

class FitView(ModelView):
	def __init__(self, frame):
		ModelView.__init__(self)
		self.fig = Figure(figsize=(400, 400))
		self.canvas = FigureCanvas(self.fig)
		layout = QHBoxLayout()
		layout.addWidget(self.canvas)
		frame.setLayout(layout)

		self.ax1 = self.fig.add_subplot(111)
		self.lines = {}

	def plot(self, key, item):
		if item.ndims == 4:
			if item.visible:
				ts = item.get_timeseries(self.fab.focus)
				if not self.lines.has_key(key):
					self.lines[key], = self.ax1.plot(ts)
				else:
					self.lines[key].set_ydata(ts)
					self.ax1.relim()
					self.ax1.autoscale_view()
			elif self.lines.has_key(key):
				self.lines[key].remove()
				del self.lines[key]

	def do_update(self):
		if self.fab.current_run:
			outdata = self.fab.current_run.data
		else:
			outdata = {}
		
		# Delete any data which no longer exists	
		for key in self.lines.keys():
			if not self.fab.data.has_key(key) and not outdata.has_key(key):
				self.lines[key].remove()
				del self.lines[key]
		
		if max(self.fab.shape) > 0:		
			for key, item in self.fab.data.items():
				self.plot(key, item)
			for key, item in outdata.items():
				self.plot(key, item)
			
		self.canvas.draw_idle()

