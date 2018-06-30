import sys
from PyQt5 import QtWidgets as qw
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class OptionsTabWidget(qw.QTabWidget):
    def __init__(self, parent=None):
        super(OptionsTabWidget, self).__init__(parent)

        self.file = None

        self.input_tab = qw.QWidget()
        self.output_tab = qw.QWidget()
        self.display_tab = qw.QWidget()
        self.dem_tab = qw.QWidget()
        self.subset_tab = qw.QWidget()
        self.figure_tab = qw.QWidget()
        self.map_tab = qw.QWidget()

        self.input_tab.setFixedHeight(300)
        self.display_tab.setFixedHeight(350)

        self.setFixedWidth(450)
        self.setFixedHeight(650)

        self.addTab(self.input_tab, "Input")
        self.addTab(self.output_tab, "Output")
        self.addTab(self.display_tab, "Display")
        self.addTab(self.dem_tab, "DEM")
        self.addTab(self.subset_tab, "Subset")
        self.addTab(self.figure_tab, "Figure")
        self.addTab(self.map_tab, "Map")

        self.input_tabUI()
        self.output_tabUI()
        self.display_tabUI()
        self.setWindowTitle("tab demo")

    def input_tabUI(self):
        layout = qw.QVBoxLayout()

        select_file_layout = qw.QHBoxLayout()
        select_file_button = qw.QPushButton("Select File")
        select_file_button.clicked.connect(self.select_file)
        file_name_label = qw.QLabel("No File Selected")

        select_file_layout.addWidget(select_file_button)
        select_file_layout.addWidget(file_name_label)

        dataset_layout = qw.QVBoxLayout()

        dset_titles_layout = qw.QHBoxLayout()
        dataset_label = qw.QLabel("Datasets")
        excludes_label = qw.QLabel("Excludes")

        dset_titles_layout.addWidget(dataset_label)
        dset_titles_layout.addWidget(excludes_label)

        dataset_widget_layout = qw.QHBoxLayout()

        dataset_list = qw.QListWidget()
        dataset_list.addItems(['1', '2', '3', '4'])

        exclude_list = qw.QListWidget()
        exclude_list.addItems(['1', '2', '3', '4'])

        dataset_widget_layout.addWidget(dataset_list)
        dataset_widget_layout.addWidget(exclude_list)

        dataset_layout.addLayout(dset_titles_layout)
        dataset_layout.addLayout(dataset_widget_layout)

        mask_file_layout = qw.QHBoxLayout()
        mask_file_button = qw.QPushButton("Select File")
        mask_name_label = qw.QLabel("No File Selected")

        mask_file_layout.addWidget(mask_file_button)
        mask_file_layout.addWidget(mask_name_label)

        layout.addLayout(select_file_layout)
        layout.addLayout(dataset_layout)
        layout.addLayout(mask_file_layout)
        layout.addWidget(qw.QCheckBox("Zero Mask"))

        self.setTabText(0, "Input")
        self.input_tab.setLayout(layout)

    def output_tabUI(self):
        layout = qw.QFormLayout()
        sex = qw.QHBoxLayout()
        sex.addWidget(qw.QRadioButton("Male"))
        sex.addWidget(qw.QRadioButton("Female"))
        layout.addRow(qw.QLabel("Sex"), sex)
        layout.addRow("Date of Birth",qw. QLineEdit())
        self.setTabText(1, "Output")
        self.output_tab.setLayout(layout)

    def display_tabUI(self):
        layout = qw.QVBoxLayout()

        minimum_layout = qw.QHBoxLayout()

        minimum_slider = qw.QSlider(Qt.Horizontal)
        min_value = qw.QLineEdit()
        min_value.setFixedWidth(75)

        minimum_layout.addWidget(qw.QLabel("Minimum"))
        minimum_layout.addWidget(minimum_slider)
        minimum_layout.addWidget(min_value)



        maximum_layout = qw.QHBoxLayout()

        maximum_slider = qw.QSlider(Qt.Horizontal)
        max_value = qw.QLineEdit()
        max_value.setFixedWidth(75)

        maximum_layout.addWidget(qw.QLabel("Maximum"))
        maximum_layout.addWidget(maximum_slider)
        maximum_layout.addWidget(max_value)


        unit_scale_layout = qw.QHBoxLayout()

        unit_label = qw.QLabel("Unit:")
        unit_drop = qw.QComboBox()
        unit_drop.addItem("M")
        unit_drop.addItem("CM")
        unit_drop.addItem("KM")
        unit_drop.addItem("MM")

        unit_scale_layout.addWidget(unit_label)
        unit_scale_layout.addWidget(unit_drop)

        scale_label = qw.QLabel("Scale:")
        scale_drop = qw.QComboBox()
        scale_drop.addItem("1.0")
        scale_drop.addItem("10.0")
        scale_drop.addItem("0.10")
        scale_drop.addItem("100.0")

        unit_scale_layout.addWidget(scale_label)
        unit_scale_layout.addWidget(scale_drop)



        colormap_layout = qw.QHBoxLayout()

        cmap_drop = qw.QComboBox()
        cmap_drop.addItem("1.0")
        cmap_drop.addItem("10.0")
        cmap_drop.addItem("0.10")
        cmap_drop.addItem("100.0")

        colormap_layout.addWidget(qw.QLabel("Colormap"))
        colormap_layout.addWidget(cmap_drop)



        projection_layout = qw.QHBoxLayout()

        proj_drop = qw.QComboBox()
        proj_drop.addItem("1.0")
        proj_drop.addItem("10.0")
        proj_drop.addItem("0.10")
        proj_drop.addItem("100.0")

        projection_layout.addWidget(qw.QLabel("Projection"))
        projection_layout.addWidget(proj_drop)



        orientation_layout = qw.QHBoxLayout()

        lr_flip_check = qw.QCheckBox("Flip LR")
        ud_flip_check = qw.QCheckBox("Flip UD")
        wrap_check = qw.QCheckBox("Wrap")
        oppo_check = qw.QCheckBox("Opposite")

        orientation_layout.addWidget(lr_flip_check)
        orientation_layout.addWidget(ud_flip_check)
        orientation_layout.addWidget(wrap_check)
        orientation_layout.addWidget(oppo_check)




        alpha_layout = qw.QHBoxLayout()

        alpha_slider = qw.QSlider(Qt.Horizontal)
        alpha_value = qw.QLineEdit()
        alpha_value.setFixedWidth(75)

        alpha_layout.addWidget(qw.QLabel("Maximum"))
        alpha_layout.addWidget(alpha_slider)
        alpha_layout.addWidget(alpha_value)



        layout.addLayout(minimum_layout)
        layout.addLayout(maximum_layout)
        layout.addLayout(unit_scale_layout)
        layout.addLayout(colormap_layout)
        layout.addLayout(projection_layout)
        layout.addLayout(orientation_layout)
        layout.addLayout(alpha_layout)

        self.setTabText(2, "Display")
        self.display_tab.setLayout(layout)


    def select_file(self):
        self.file = qw.QFileDialog.getOpenFileName(self, 'Open file', '/', "Data files (*.h5 *.he5)")[0]
        print(self.file)


    def get_options(self):
        return self.file

def main():
    app = qw.QApplication(sys.argv)
    ex = OptionsTabWidget()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
