# region imports
from X2Q2_SP24 import doPlot, simulate
from PyQt5 import QtWidgets as qtw
from PyQt5 import QtGui as qtg
from PyQt5 import QtCore as qtc
from abc import ABC, abstractmethod

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT
from matplotlib.figure import Figure
# endregion

# region RLC circuit classes (MVC)
class circuitModel():
    def __init__(self):
        """
        Model container for the circuit problem.
        """
        self.nodes = []
        self.resistors = []
        self.capacitors = []
        self.inductors = []
        self.voltageSources = []
        self.wires = []


class circuitView():
    def __init__(self, dw=None):
        """
        View class for circuit GUI.
        :param dw: display widgets tuple
        """
        if dw is not None:
            self.setDisplayWidgets(dw)
            self.setupImageLabel()
            self.setupPlot()

    def setDisplayWidgets(self, dw=None):
        """
        Unpack display widgets.
        :param dw: (layout_VertMain, layout_VertInput, form)
        :return: nothing
        """
        if dw is not None:
            # ChatGPT helped write this function.
            self.layout_VertMain, self.layout_VertInput, self.form = dw

    def setupImageLabel(self):
        """
        Displays picture of circuit from Circuit1.png in a label widget
        :return:
        """
        # ChatGPT helped write this function.
        self.pixMap = qtg.QPixmap()
        self.pixMap.load("Circuit1.png")
        self.image_label = qtw.QLabel()
        self.image_label.setPixmap(self.pixMap)
        self.image_label.setAlignment(qtc.Qt.AlignCenter)
        self.layout_VertInput.addWidget(self.image_label)

    def setupPlot(self):
        """
        Create the figure, canvas, axes and toolbar objects and place them on GUI
        :return:
        """
        self.figure = Figure(figsize=(8, 8), tight_layout=True, frameon=True, facecolor='none')
        self.canvas = FigureCanvasQTAgg(self.figure)
        self.ax = self.figure.add_subplot()
        self.toolbar = NavigationToolbar2QT(self.canvas, self.form)
        self.layout_VertMain.addWidget(self.toolbar)
        self.layout_VertMain.addWidget(self.canvas)

    def doPlot(self, args):
        """
        Plot simulation results to GUI canvas.
        :param args: (R, time, solution object)
        :return: nothing
        """
        self.figure.clear()
        self.ax = self.figure.add_subplot()
        doPlot(args, ax=self.ax)
        self.canvas.draw()


class circuitController():
    def __init__(self, args):
        """
        This is the class for a circuitController.  It has a model and view for the circuit.
        :param args: a tuple with input widgets and display widgets
        """
        self.inputWidgets, self.displayWidgets = args

        # ChatGPT helped write this function.
        # unpack input widgets
        (self.le_Inductance,
         self.le_Resistance,
         self.le_Capacitence,
         self.le_Amplitude,
         self.le_Freq,
         self.le_Phase,
         self.le_simTime,
         self.le_simPts) = self.inputWidgets

        self.Model = circuitModel()
        self.View = circuitView(dw=self.displayWidgets)

    def calculate(self):
        """
        Simulates the circuit by calling from X2Q2_SP24 functions.
        Step 1:  read inputs from GUI.
        Step 2:  call simulate.
        Step 3:  call doPlot.
        :return:
        """
        # ChatGPT helped write this function.
        L = float(self.le_Inductance.text())
        R = float(self.le_Resistance.text())
        C = float(self.le_Capacitence.text())
        A = float(self.le_Amplitude.text())
        f = float(self.le_Freq.text())
        p = float(self.le_Phase.text())
        t = float(self.le_simTime.text())
        pts = int(float(self.le_simPts.text()))

        I = simulate(L=L, R=R, C=C, A=A, f=f, p=p, t=t, pts=pts)
        self.View.doPlot((R, I.t, I))

# endregion