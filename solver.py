from astroquery.astrometry_net import AstrometryNet
from astropy.wcs import WCS

import utils

import os
import threading

class PlateSolve(threading.Thread):
    def __init__(self):
        super().__init__()
        self.folder = None
        self.fits = None
        self.result = None
        self.state = False
        self.stop_flag = threading.Event()
        self.ast = AstrometryNet()

        self.pixel_x = 1024
        self.pixel_y = 1024
        self.ast.api_key = 'epepiumlsedbruam'  # Get an API key at https://nova.astrometry.net/
    
    def update_fields(self, dir, fits):
        self.folder = dir
        self.fits = fits
    
    def get_latest_file(self, folder_path):
        files = [os.path.join(folder_path, file) for file in os.listdir(folder_path)]
        if not files:
            return None
        latest_file = max(files, key=os.path.getctime)
        return latest_file
    
    def stop(self):
        self.state = False 
        self.stop_flag.set()

    def run(self):  
        self.state = True 
        self.result = None     
        try_again = True
        wcs_header = None
        submission_id = None
        while try_again and not self.stop_flag.is_set():
            try:
                if self.folder:
                    fits_path = self.get_latest_file(self.folder)
                elif self.fits:
                    fits_path = self.fits
                else:
                    self.result = "Invalid File or Folder."
                if not submission_id:
                    wcs_header = self.ast.solve_from_image(fits_path,
                                                    submission_id=submission_id, force_image_upload=True)
                else:
                    wcs_header = self.ast.monitor_submission(submission_id,
                                                        solve_timeout=60)
            except TimeoutError as e:
                submission_id = e.args[1]
                self.stop()
            except Exception as e:
                self.result = "Stopped"
                self.stop()
            else:
                try_again = False
        if not self.stop_flag.is_set():
            if wcs_header:      
                w = WCS(wcs_header)     
                sky = w.pixel_to_world(int(self.pixel_x), int(self.pixel_y))
                ra_hours = sky.ra.hour  
                dec_degrees = sky.dec.degree                                
                self.result = (ra_hours, dec_degrees, w)
            else:
                self.result = "No headers."
        else:
            self.result = "Stopped"
    
    

    



