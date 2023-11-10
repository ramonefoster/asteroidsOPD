from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import twirl
import time
import random
from astropy.visualization import simple_norm

from PyQt5 import QtWidgets, uic
from PyQt5.QtCore import QTimer, Qt, pyqtSlot
from PyQt5.QtGui import *
from PyQt5.QtWidgets import QFileDialog
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

        self._abort = threading.Event()

        self.timer = QTimer()
        self.timer.timeout.connect(self.update)
        self.timer_status()

        self.plot_ready = False
    
    def timer_status(self):
        self.timer.start(1000)

    def solve_field(self):  
        self.plot_ready = False 
        if self.plate_solve.state:
            self.plate_solve.stop()
            self.plate_solve.join()
            self.plate_solve = PlateSolve()
        else:
            if self.txtFITS.text().endswith(".fits"):
                self.plate_solve.fits = self.txtFITS.text()
                self.plate_solve.folder = None
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
        
    def plot_fits(self, hdu, ra_center, dec_center, w):
        header = hdu.header 
        data = hdu.data

        ra, dec = ra_center*15, dec_center
        center = SkyCoord(ra, dec, unit=["deg", "deg"])
        center = [center.ra.value, center.dec.value]

        RA_chiron, DEC_chiron = self.txtRA.text(), self.txtDEC.text() 
        ra, dec = utils.hms_to_hours(RA_chiron)*15, utils.dms_to_degrees(DEC_chiron) 

        if not ra or not dec:
            return
        
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
        asteroid = np.array(SkyCoord(ra, dec, unit=["deg", "deg"]).to_pixel(wcs)).T

        ax.plot(*gaias_pixel.T, "o", fillstyle="none", c="C1", ms=18)
        ax.plot(*asteroid.T, "o", fillstyle="none", ms=18, color="green")

        # Annotate the asteroids
        ax.annotate(header["OBJECT"], asteroid, color='white', xytext=(asteroid[0]+50, asteroid[1]+40),
                    bbox=dict(boxstyle="round", alpha=0.4, color='green'), fontsize=10)

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