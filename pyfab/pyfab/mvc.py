import traceback

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
    def __init__(self, changes, *widgets, **kwidgets):
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
            if not obj: self.set_enabled(False)
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

