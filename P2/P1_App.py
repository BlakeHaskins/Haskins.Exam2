#region imports
from PyQt5.QtWidgets import QApplication, QWidget
from P1_GUI import Ui_MainForm
import sys
from Circuit_Classes import circuitController
#endregion

#region class definitions
class main_window(Ui_MainForm, QWidget):
    def __init__(self):
        """
        Constructor for circuit simulator.
        """
        super().__init__()
        self.setupUi(self)
        # ChatGPT helped write this function.
        self.setWindowTitle("Exam 2 Problem 2 - RLC Circuit Simulator")

        self.inputWidgets = (
            self.le_Inductance,
            self.le_Resistance,
            self.le_Capacitence,
            self.le_Amplitude,
            self.le_Freq,
            self.le_Phase,
            self.le_simTime,
            self.le_simPts
        )

        self.displayWidgets = (
            self.layout_VertMain,
            self.layout_VertInput,
            self
        )

        self.controller = circuitController((self.inputWidgets, self.displayWidgets))
        self.setupSignalsAndSlots()
        self.show()

    def setupSignalsAndSlots(self):
        """
        Connect the push button to the calculate function.
        :return:
        """
        self.pb_Calculate.clicked.connect(self.calculate)

    def calculate(self):
        self.controller.calculate()

#endregion

if __name__ == "__main__":
    app = QApplication.instance()
    if not app:
        app = QApplication(sys.argv)
    app.aboutToQuit.connect(app.deleteLater)
    main_win = main_window()
    sys.exit(app.exec_())