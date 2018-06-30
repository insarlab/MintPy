import sys
from PyQt5 import QtWidgets as qw
from PyQt5.QtCore import Qt
from test import tabdemo
from geoview_widget import GeoViewWidget


class MainWindow(qw.QWidget):

    def __init__(self, parent=None):

        super(MainWindow, self).__init__(parent)

        self.setFixedWidth(1200)

        self.controls_widget = tabdemo()
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
        options = self.controls_widget.get_options()
        print(options)
        self.plot_widget.deleteLater()
        self.plot_widget = GeoViewWidget(iargs=[options])
        self.layout.addWidget(self.plot_widget)

def main():
    app = qw.QApplication(sys.argv)
    ex = MainWindow()
    ex.show()
    sys.exit(app.exec_())

if __name__ == '__main__':
    main()