import sys
from PyQt5 import QtWidgets as qw
from PyQt5.QtCore import Qt, QRegExp
from PyQt5 import QtGui as qg

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure


class OptionsTabWidget(qw.QTabWidget):
    def __init__(self, parent=None):
        super(OptionsTabWidget, self).__init__(parent)

        self.file = None
        self.min = None
        self.max = None
        self.alpha = None
        self.dem_file = None

        self.select_file_button = None
        self.minimum_slider = None
        self.maximum_slider = None
        self.alpha_slider = None
        self.min_value = None
        self.max_value = None
        self.alpha_value = None
        self.cmap_drop = None
        self.proj_drop = None
        self.lr_flip_check = None
        self.ud_flip_check = None
        self.wrap_check = None
        self.oppo_check = None
        self.unit_drop = None
        self.dem_file_button = None
        self.dem_noshade_check = None
        self.dem_nocontour_check = None
        self.dem_smoothing_text = None
        self.dem_contour_text = None

        self.input_tab = qw.QWidget()
        self.output_tab = qw.QWidget()
        self.display_tab = qw.QWidget()
        self.dem_tab = qw.QWidget()
        self.subset_tab = qw.QWidget()
        self.figure_tab = qw.QWidget()
        self.map_tab = qw.QWidget()

        self.input_tab.setFixedHeight(300)
        self.display_tab.setFixedHeight(350)
        self.dem_tab.setFixedHeight(200)

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
        self.dem_tabUI()

        self.setWindowTitle("tab demo")

    def input_tabUI(self):
        layout = qw.QVBoxLayout()

        select_file_layout = qw.QHBoxLayout()
        self.select_file_button = qw.QPushButton("Select File")
        self.select_file_button.clicked.connect(lambda: self.select_file(self.select_file_button, file_name_label))
        file_name_label = qw.QLabel("No File Selected")

        select_file_layout.addWidget(self.select_file_button)
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

        self.minimum_slider = qw.QSlider(Qt.Horizontal)
        self.minimum_slider.setMinimum(-1000)
        self.minimum_slider.setMaximum(10000)
        self.minimum_slider.setValue(1000)

        self.minimum_slider.valueChanged.connect(lambda: self.slider_changed(self.minimum_slider))

        self.min_value = qw.QLineEdit()
        self.min_value.setFixedWidth(75)
        self.min_value.textChanged.connect(self.min_text_changed)

        numbers_validator = qg.QRegExpValidator(QRegExp("-?[0-9]*\.*[0-9]*"), self.min_value)
        self.min_value.setValidator(numbers_validator)

        minimum_layout.addWidget(qw.QLabel("Minimum"))
        minimum_layout.addWidget(self.minimum_slider)
        minimum_layout.addWidget(self.min_value)



        maximum_layout = qw.QHBoxLayout()

        self.maximum_slider = qw.QSlider(Qt.Horizontal)
        self.maximum_slider.setMinimum(-1000)
        self.maximum_slider.setMaximum(10000)
        self.maximum_slider.setValue(1000)

        self.maximum_slider.valueChanged.connect(lambda: self.slider_changed(self.maximum_slider))

        self.max_value = qw.QLineEdit()
        self.max_value.setFixedWidth(75)
        self.max_value.textChanged.connect(self.max_text_changed)
        self.max_value.setValidator(numbers_validator)

        maximum_layout.addWidget(qw.QLabel("Maximum"))
        maximum_layout.addWidget(self.maximum_slider)
        maximum_layout.addWidget(self.max_value)


        unit_scale_layout = qw.QHBoxLayout()

        unit_label = qw.QLabel("Unit:")
        self.unit_drop = qw.QComboBox()
        self.unit_drop.addItem("m")
        self.unit_drop.addItem("cm")
        self.unit_drop.addItem("mm")
        self.unit_drop.addItem("km")

        unit_scale_layout.addWidget(unit_label)
        unit_scale_layout.addWidget(self.unit_drop)

        scale_label = qw.QLabel("Scale:")
        scale_drop = qw.QComboBox()
        scale_drop.addItem("1.0")
        scale_drop.addItem("10.0")
        scale_drop.addItem("0.10")
        scale_drop.addItem("100.0")

        unit_scale_layout.addWidget(scale_label)
        unit_scale_layout.addWidget(scale_drop)





        colormap_layout = qw.QHBoxLayout()

        self.cmap_drop = qw.QComboBox()

        colormaps = ['Accent', 'Accent_r', 'Blues', 'Blues_r', 'BrBG', 'BrBG_r', 'BuGn', 'BuGn_r', 'BuPu', 'BuPu_r', 'CMRmap',
                     'CMRmap_r', 'Dark2', 'Dark2_r', 'GnBu', 'GnBu_r', 'Greens', 'Greens_r', 'Greys', 'Greys_r', 'OrRd', 'OrRd_r',
                     'Oranges', 'Oranges_r', 'PRGn', 'PRGn_r', 'Paired', 'Paired_r', 'Pastel1', 'Pastel1_r', 'Pastel2',
                     'Pastel2_r', 'PiYG', 'PiYG_r', 'PuBu', 'PuBuGn', 'PuBuGn_r', 'PuBu_r', 'PuOr', 'PuOr_r', 'PuRd', 'PuRd_r',
                     'Purples', 'Purples_r', 'RdBu', 'RdBu_r', 'RdGy', 'RdGy_r', 'RdPu', 'RdPu_r', 'RdYlBu', 'RdYlBu_r', 'RdYlGn',
                     'RdYlGn_r', 'Reds', 'Reds_r', 'Set1', 'Set1_r', 'Set2', 'Set2_r', 'Set3', 'Set3_r', 'Spectral',
                     'Spectral_r', 'Vega10', 'Vega10_r', 'Vega20', 'Vega20_r', 'Vega20b', 'Vega20b_r', 'Vega20c', 'Vega20c_r',
                     'Wistia', 'Wistia_r', 'YlGn', 'YlGnBu', 'YlGnBu_r', 'YlGn_r', 'YlOrBr', 'YlOrBr_r', 'YlOrRd', 'YlOrRd_r',
                     'afmhot', 'afmhot_r', 'autumn', 'autumn_r', 'binary', 'binary_r', 'bone', 'bone_r', 'brg', 'brg_r', 'bwr',
                     'bwr_r', 'cool', 'cool_r', 'coolwarm', 'coolwarm_r', 'copper', 'copper_r', 'cubehelix', 'cubehelix_r',
                     'flag', 'flag_r', 'gist_earth', 'gist_earth_r', 'gist_gray', 'gist_gray_r', 'gist_heat', 'gist_heat_r', 'gist_ncar',
                     'gist_ncar_r', 'gist_rainbow', 'gist_rainbow_r', 'gist_stern', 'gist_stern_r', 'gist_yarg',
                     'gist_yarg_r', 'gnuplot', 'gnuplot2', 'gnuplot2_r', 'gnuplot_r', 'gray', 'gray_r', 'hot', 'hot_r', 'hsv', 'hsv_r', 'inferno',
                     'inferno_r', 'jet', 'jet_r', 'magma', 'magma_r', 'nipy_spectral', 'nipy_spectral_r', 'ocean', 'ocean_r', 'pink',
                     'pink_r', 'plasma', 'plasma_r', 'prism', 'prism_r', 'rainbow', 'rainbow_r', 'seismic', 'seismic_r',
                     'spectral', 'spectral_r', 'spring', 'spring_r', 'summer', 'summer_r', 'tab10', 'tab10_r', 'tab20', 'tab20_r',
                     'tab20b', 'tab20b_r', 'tab20c', 'tab20c_r', 'terrain', 'terrain_r', 'viridis', 'viridis_r', 'winter', 'winter_r']

        for cmap in colormaps:
            self.cmap_drop.addItem(cmap)

        colormap_layout.addWidget(qw.QLabel("Colormap"))
        colormap_layout.addWidget(self.cmap_drop)




        projection_layout = qw.QHBoxLayout()

        self.proj_drop = qw.QComboBox()

        projections = ["cea", "mbtfpq", "aeqd", "sinu", "poly", "moerc", "gnom", "moll", "lcc", "tmerc", "nplaea",
                       "gall", "npaeqd", "mill", "merc", "stere", "eqdc", "rotpole", "cyl", "npstere", "spstere", "hammer",
                       "geos", "nsper", "eck4", "aea", "kav7", "spaeqd", "ortho", "cass", "vandg", "laea", "splaea", "robin"]

        for proj in projections:
            self.proj_drop.addItem(proj)

        projection_layout.addWidget(qw.QLabel("Projection"))
        projection_layout.addWidget(self.proj_drop)



        orientation_layout = qw.QHBoxLayout()

        self.lr_flip_check = qw.QCheckBox("Flip LR")
        self.ud_flip_check = qw.QCheckBox("Flip UD")
        self.wrap_check = qw.QCheckBox("Wrap")
        self.oppo_check = qw.QCheckBox("Opposite")

        orientation_layout.addWidget(self.lr_flip_check)
        orientation_layout.addWidget(self.ud_flip_check)
        orientation_layout.addWidget(self.wrap_check)
        orientation_layout.addWidget(self.oppo_check)




        alpha_layout = qw.QHBoxLayout()

        self.alpha_slider = qw.QSlider(Qt.Horizontal)
        self.alpha_slider.setMinimum(0)
        self.alpha_slider.setMaximum(100)
        self.alpha_slider.setValue(100)

        self.alpha_slider.valueChanged.connect(lambda: self.slider_changed(self.alpha_slider))

        self.alpha_value = qw.QLineEdit()
        self.alpha_value.setFixedWidth(75)
        self.alpha_value.textChanged.connect(self.alpha_text_changed)
        self.alpha_value.setValidator(numbers_validator)

        alpha_layout.addWidget(qw.QLabel("Maximum"))
        alpha_layout.addWidget(self.alpha_slider)
        alpha_layout.addWidget(self.alpha_value)



        layout.addLayout(minimum_layout)
        layout.addLayout(maximum_layout)
        layout.addLayout(unit_scale_layout)
        layout.addLayout(colormap_layout)
        layout.addLayout(projection_layout)
        layout.addLayout(orientation_layout)
        layout.addLayout(alpha_layout)

        self.setTabText(2, "Display")
        self.display_tab.setLayout(layout)

    def dem_tabUI(self):

        layout = qw.QVBoxLayout()

        select_file_layout = qw.QHBoxLayout()
        self.dem_file_button = qw.QPushButton("Select File")
        self.dem_file_button.clicked.connect(lambda: self.select_file(self.dem_file_button, dem_file_name_label))
        dem_file_name_label = qw.QLabel("No File Selected")

        select_file_layout.addWidget(self.dem_file_button)
        select_file_layout.addWidget(dem_file_name_label)

        dem_options_layout = qw.QHBoxLayout()

        self.dem_noshade_check = qw.QCheckBox("No Shade")
        self.dem_nocontour_check = qw.QCheckBox("No Contour")

        dem_options_layout.addWidget(self.dem_noshade_check)
        dem_options_layout.addWidget(self.dem_nocontour_check)


        dem_contour_options_layout = qw.QVBoxLayout()

        contour_options_titles_layout = qw.QHBoxLayout()
        smoothing_label = qw.QLabel("Contour Smoothing")
        step_label = qw.QLabel("Contour Step")

        contour_options_titles_layout.addWidget(smoothing_label)
        contour_options_titles_layout.addWidget(step_label)

        contour_options_widget_layout = qw.QHBoxLayout()

        self.dem_smoothing_text = qw.QLineEdit()
        self.dem_contour_text = qw.QLineEdit()

        self.dem_smoothing_text.setText("3.0")
        self.dem_contour_text.setText("200")

        contour_options_widget_layout.addWidget(self.dem_smoothing_text)
        contour_options_widget_layout.addWidget(self.dem_contour_text)

        dem_contour_options_layout.addLayout(contour_options_titles_layout)
        dem_contour_options_layout.addLayout(contour_options_widget_layout)

        layout.addLayout(select_file_layout)
        layout.addLayout(dem_options_layout)
        layout.addLayout(dem_contour_options_layout)

        self.setTabText(3, "DEM")
        self.dem_tab.setLayout(layout)

    def select_file(self, file_button, file_label):
        print("select file")
        if file_button is self.select_file_button:
            self.file = qw.QFileDialog.getOpenFileName(self, 'Open file', '/', "Data files (*.h5 *.he5)")[0]
            file_label.setText(self.file.split("/")[-1])
        elif file_button is self.dem_file_button:
            self.dem_file = qw.QFileDialog.getOpenFileName(self, 'Open file', '/', "Data files (*.h5 *.he5)")[0]
            file_label.setText(self.dem_file.split("/")[-1])

        file_button.setText("Cancel")


    def slider_changed(self, slider):

        if slider is self.minimum_slider:
            self.min = self.minimum_slider.value()
            self.min_value.setText(str(self.min))
        elif slider is self.maximum_slider:
            self.max = self.maximum_slider.value()
            self.max_value.setText(str(self.max))
        else:
            self.alpha = self.alpha_slider.value()/100
            self.alpha_value.setText(str(self.alpha))

    def max_text_changed(self, text):
        self.max = int(text)
        self.maximum_slider.setValue(self.max)

    def min_text_changed(self, text):
        self.min = int(text)
        self.minimum_slider.setValue(self.min)

    def alpha_text_changed(self, text):
        self.alpha = float(text)
        self.alpha_slider.setValue(self.alpha*100)

    def get_options(self):
        return self.file, self.min, self.max, self.cmap_drop.currentText(), self.proj_drop.currentText(), \
               self.lr_flip_check.isChecked(), self.ud_flip_check.isChecked(), self.wrap_check.isChecked(), \
               self.oppo_check.isChecked(), self.unit_drop.currentText(), self.dem_file, \
               self.dem_noshade_check.isChecked(), self.dem_nocontour_check.isChecked(), \
               self.dem_contour_text.text(), self.dem_smoothing_text.text()

def main():
    app = qw.QApplication(sys.argv)
    ex = OptionsTabWidget()
    ex.show()
    sys.exit(app.exec_())


if __name__ == '__main__':
    main()
