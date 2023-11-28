from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astroquery.jplhorizons import Horizons
import twirl
import time
import random
from astropy.visualization import simple_norm

from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, Qt, QDateTime
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QMenu, QFileDialog, QVBoxLayout, QPushButton, QLineEdit, QDialog, QLabel
import threading
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas

import sys
import utils
from solver import PlateSolve

Ui_MainWindow, QtBaseClass = uic.loadUiType("UI.ui")

class Solver(QtWidgets.QMainWindow, Ui_MainWindow):
    def __init__(self):
        QtWidgets.QMainWindow.__init__(self)
        Ui_MainWindow.__init__(self)
        self.setupUi(self)

        self.setTabOrder(self.txtFITS, self.txtRA)
        self.setTabOrder(self.txtRA, self.txtDEC)
        self.setTabOrder(self.txtDEC, self.spinFOV)

        self.plate_solve = PlateSolve()

        self.btnGo.clicked.connect(self.solve_field)
        self.btnFile.clicked.connect(self.get_file)

        self.btnSimbad.clicked.connect(self.choose_object)

        self._abort = threading.Event()

        self.timer = QTimer()
        self.timer.timeout.connect(self.update)
        self.timer_status()

        self.startDate.setDisplayFormat("yyyy-MM-dd hh:mm:ss")
        self.endDate.setDisplayFormat("yyyy-MM-dd hh:mm:ss")

        current_datetime = QDateTime.currentDateTime()
        self.startDate.setDateTime(current_datetime)
        self.endDate.setDateTime(current_datetime)

        self.image_label.setContextMenuPolicy(Qt.CustomContextMenu)
        self.image_label.customContextMenuRequested.connect(self.show_context_menu)

        self.plot_ready = False
    
    def save_image(self):
        image_path, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "PNG Files (*.png);;JPEG Files (*.jpg *.jpeg)")

        if image_path:
            pixmap = self.image_label.pixmap()
            if pixmap:
                pixmap.save(image_path)
    
    def show_context_menu(self, position):
        menu = QMenu(self)
        menu.setStyleSheet("QMenu { background-color: #C0C0C0; }")

        save_action = menu.addAction("Save Image")
        save_action.triggered.connect(self.save_image)

        menu.exec_(self.image_label.mapToGlobal(position))
    
    def choose_object(self):
        dialog = QDialog(self)
        dialog.setWindowTitle("Choose Object")
        dialog.setFixedWidth(280)

        layout = QVBoxLayout(dialog)

        label = QLabel("Insert SIMBAD ID:", dialog)
        label.setStyleSheet("color: rgb(255, 255, 255);")
        layout.addWidget(label)

        input_text = QLineEdit(dialog)
        input_text.setStyleSheet("font-size: 15px; color: white;")
        input_text.setFixedHeight(30)

        layout.addWidget(input_text)

        search_button = QPushButton('Search', dialog)
        search_button.setStyleSheet("color: rgb(255, 255, 255);")
        cancel_button = QPushButton('Cancel', dialog)
        cancel_button.setStyleSheet("color: rgb(255, 255, 255);")

        search_button.clicked.connect(dialog.accept)
        cancel_button.clicked.connect(dialog.reject)

        layout.addWidget(search_button)
        layout.addWidget(cancel_button)

        result = dialog.exec_()

        if result == QDialog.Accepted:
            entered_text = input_text.text()            
            try:
                ra, dec = utils.simbad_radec(entered_text)
                self.txtRA.setText(ra)
                self.txtDEC.setText(dec)
            except Exception as e:
                print("Error Simbad"+str(e))            
    
    def timer_status(self):
        self.timer.start(1000)

    def solve_field(self):  
        self.plot_ready = False 
        if self.plate_solve.is_alive():
            self.plate_solve.stop()
            self.plate_solve.join()

        self.plate_solve = PlateSolve()
        if self.txtFITS.text().endswith(".fits"):
            self.plate_solve.fits = self.txtFITS.text()
            self.plate_solve.folder = None 
            hdu = fits.open(self.txtFITS.text())[0]
            data = hdu.data
            image_shape = data.shape
            if len(image_shape) == 2:
                h, w = image_shape
            if len(image_shape) == 3:
                x, h, w = image_shape
            self.plate_solve.pixel_x = w/2
            self.plate_solve.pixel_y = h/2

            self.plate_solve.start()
        else:
            self.labelPlateSolve.setText("Invalid FITS.")
    
    def update(self):
        if round(time.time(),0) % 5 == 0:
            # PLATE SOLVING
            if self.plate_solve.state and not self.plate_solve.result:
                self.labelPlateSolve.setText("Solving... May take a while, please Wait.")
                self.progressBar.setValue(self.progressBar.value()+random.randint(10, 20))
            elif self.plate_solve.result:
                if isinstance(self.plate_solve.result, tuple):
                    ra_center = self.plate_solve.result[0]
                    dec_center = self.plate_solve.result[1]
                    wcs = self.plate_solve.result[2]
                    path_img = self.txtFITS.text()
                    try:
                        #filename = get_pkg_data_filename(path_img)
                        hdu = fits.open(self.txtFITS.text())[0]
                        self.txtRASolved.setText(utils.hours_to_hms(ra_center))
                        self.txtDECSolved.setText(utils.degrees_to_dms(dec_center))
                        if not self.plot_ready: 
                            self.plot_fits(hdu, ra_center, dec_center, wcs)
                            self.plot_ready = True 
                        self.labelPlateSolve.setText("Completed")
                    except Exception as e:
                        print(e)
                        # self.plate_solve.stop()
                        self.labelPlateSolve.setText("Error"+str(e))
                    self.progressBar.setValue(0)
                elif isinstance(self.plate_solve.result, str):
                    self.txtRASolved.setText('')
                    self.txtDECSolved.setText('')
                    self.labelPlateSolve.setText(self.plate_solve.result)
                    self.progressBar.setValue(0)
                else:
                    self.labelPlateSolve.setText("Error...")
                    self.progressBar.setValue(0)
    
    def get_file(self):
        self.txtFITS.setText(QFileDialog.getOpenFileName(None, "Open File", "C:/", "fits (*.fits)")[0])
    
    def calc_ephem(self):
        opd = {'lon': -45.5825,
                'lat': -22.5344,
                'elevation': 1.864}
        
        try:
            epochs = {'start':self.startDate.dateTime().toString("yyyy-MM-dd hh:mm:ss"), 
                                'stop':self.endDate.dateTime().toString("yyyy-MM-dd hh:mm:ss"),
                                'step':f'{self.spinStep.value()}{self.boxStep.currentText()}'}

            obj = Horizons(id=self.txtOBJ.text(), location=opd,
                        epochs=epochs)

            eph = obj.ephemerides()
            return (eph["targetname"], eph["RA"], eph["DEC"], eph["datetime_str"])
        except:
            return None
        
    def plot_fits(self, hdu, ra_center, dec_center, w):
        header = hdu.header 
        data = hdu.data

        ra, dec = ra_center*15, dec_center
        center = SkyCoord(ra, dec, unit=["deg", "deg"])
        center = [center.ra.value, center.dec.value]

        ra_target, dec_target = None, None
        try:
            ra_target, dec_target = self.txtRA.text(), self.txtDEC.text()
            ra_target, dec_target = utils.hms_to_hours(ra_target)*15, utils.dms_to_degrees(dec_target) 
        except:
            pass

        # if not ra or not dec:
        #     return
        
        data = np.squeeze(data)

        # Create WCS
        wcs = w

        # Create a figure and WCSAxes
        fig = plt.figure(figsize=(8, 8), facecolor='dimgrey')
        
        ax = fig.add_subplot(1, 1, 1, projection=wcs)

        # Display the image
        norm = simple_norm(data, 'linear', percent=99.5)
        ax.imshow(data, norm=norm, origin='lower', cmap="gray")

        # Plot stars and asteroids
        gaias = twirl.gaia_radecs(center, self.spinFOV.value()/60, limit=self.spinStars.value())
        gaias_pixel = np.array(SkyCoord(gaias, unit="deg").to_pixel(wcs)).T

        self.labelPlateSolve.setText("Getting Object Info...")
        small_bodies = self.calc_ephem()
        self.labelPlateSolve.setText("Done! Now Plotting Image.")
        if small_bodies:
            for i in range(len(small_bodies[0])):
                ra = (small_bodies[1][i])
                dec = (small_bodies[2][i])
                name = (small_bodies[0][i])
                date_time = (small_bodies[3][i])
                asteroid = np.array(SkyCoord(ra, dec, unit=["deg", "deg"]).to_pixel(wcs)).T
                if i == 0:
                    color = 'green'
                    ax.annotate(name, asteroid, color='white', xytext=(asteroid[0]+50, asteroid[1]+40),
                        bbox=dict(boxstyle="round", alpha=0.4, color=color), fontsize=10)
                else:
                    color = 'red'
                    ax.annotate(date_time, asteroid, color='white', xytext=(asteroid[0]+50, asteroid[1]+40),
                        bbox=dict(boxstyle="round", alpha=0.4, color=color), fontsize=10)
                
                ax.plot(*asteroid.T, "o", fillstyle="none", ms=18, color=color)
                # Annotate the asteroids
                # ax.annotate(name, asteroid, color='white', xytext=(asteroid[0]+50, asteroid[1]+40),
                #         bbox=dict(boxstyle="round", alpha=0.4, color=color), fontsize=10)
        
        if ra_target and dec_target:
            target = np.array(SkyCoord(ra_target, dec_target, unit=["deg", "deg"]).to_pixel(wcs)).T
            ax.plot(*target.T, "s", fillstyle="none", ms=18, color="gold")
            # Annotate the targets
            ax.annotate("Target", target, color='white', xytext=(target[0]+50, target[1]+40),
                    bbox=dict(boxstyle="round", alpha=0.4, color="gold"), fontsize=10)

        ax.plot(*gaias_pixel.T, "o", fillstyle="none", c="C1", ms=18)
       
        # Set axis labels
        ax.set_xlabel('Right Ascension (J2000)')
        ax.set_ylabel('Declination (J2000)')

        ax.xaxis.label.set_color('white')
        ax.tick_params(axis='x', colors='white')
        ax.yaxis.label.set_color('white')
        ax.tick_params(axis='y', colors='white')

        fig.tight_layout(pad=5)

        canvas = FigureCanvas(fig)
        canvas.draw()

        width, height = int(fig.get_size_inches()[0] * fig.get_dpi()), int(fig.get_size_inches()[1] * fig.get_dpi())
        image = QImage(width, height, QImage.Format_RGB32)
        image.fill(Qt.white)

        painter = QPainter(image)
        canvas.render(painter)
        painter.end()

        # Display the QImage in the QLabel, scaled to fit
        pixmap = QPixmap.fromImage(image)
        pixmap = pixmap.scaled(self.image_label.size(), aspectRatioMode=Qt.KeepAspectRatio)
        self.image_label.setPixmap(pixmap)

        # Display the QImage in the QLabel
        self.image_label.setPixmap(QPixmap.fromImage(image))

    def showDialog(self, msgError):
        """
        Display a message box showing a warning or information message.

        Parameters:
            msgError (str): The message to be displayed in the message box.

        """
        msgBox = QtWidgets.QMessageBox()
        msgBox.setIcon(QtWidgets.QMessageBox.Information)
        msgBox.setText(msgError)
        msgBox.setWindowTitle("Warning")
        msgBox.setStandardButtons(QtWidgets.QMessageBox.Ok)
        msgBox.exec()

if __name__ == "__main__":
    main_app = QtWidgets.QApplication(sys.argv)
    window = Solver()

    window.show()
    sys.exit(main_app.exec_())