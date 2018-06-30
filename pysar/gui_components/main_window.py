import sys
from PyQt5 import QtWidgets as qw
from PyQt5.QtCore import Qt
from options_tab_widget import OptionsTabWidget
from geoview_widget import GeoViewWidget


class MainWindow(qw.QWidget):

    def __init__(self, parent=None):

        super(MainWindow, self).__init__(parent)

        self.setFixedWidth(1200)

        self.controls_widget = OptionsTabWidget()
        self.plot_widget = GeoViewWidget(iargs=['/Users/joshua/Desktop/pysar/test_data/new_data/velocity.h5'])

        self.layout = qw.QHBoxLayout(self)

        controls_layout = qw.QVBoxLayout()
        controls_layout.addWidget(self.controls_widget)
        plot_button = qw.QPushButton("Plot")
        plot_button.clicked.connect(self.plot_new)
        controls_layout.addWidget(plot_button)


        self.layout.addLayout(controls_layout)
        self.layout.addWidget(self.plot_widget)

    def plot_new(self):
        option_values = list(self.controls_widget.get_options())
        print(option_values)

        options = [option_values[0], "-m", str(option_values[1]), "-M", str(option_values[2]), "-c", option_values[3], "--projection", option_values[4]]
        if option_values[5] is True:
            options.append("--flip-lr")
        if option_values[6] is True:
            options.append("--flip-ud")
        if option_values[7] is True:
            options.append("--wrap")
        if option_values[8] is True:
            options.append("--opposite")

        self.plot_widget.deleteLater()
        self.plot_widget = GeoViewWidget(iargs=options)
        self.layout.addWidget(self.plot_widget)

def main():
    app = qw.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()