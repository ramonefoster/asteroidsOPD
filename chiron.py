from astropy.coordinates import SkyCoord
import matplotlib.pyplot as plt
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
from astroquery.astrometry_net import AstrometryNet
import twirl
from astropy import units as u
import utils
import astropy.visualization.wcsaxes as axes
from astropy.visualization import simple_norm

def solve(ast, path):  
    try_again = True
    wcs_header = None
    submission_id = None
    while try_again:
        try:            
            if not submission_id:
                wcs_header = ast.solve_from_image(path,
                            submission_id=submission_id, force_image_upload=True)
            else:
                wcs_header = ast.monitor_submission(submission_id,
                            solve_timeout=90)
        except TimeoutError as e:
            submission_id = e.args[1]
        else:
            try_again = False
    
    if wcs_header:
        w = WCS(wcs_header)
        sky = w.pixel_to_world(1024, 1024)
        ra = sky.ra.hour  # Right Ascension 
        dec = sky.dec.degree  # Declination 
        return ra, dec, w
    else:
        return None
    
def plot_fits(hdu, ra_center, dec_center, w):

    header = hdu.header 
    data = hdu.data

    ra, dec = ra_center*15, dec_center
    center = SkyCoord(ra, dec, unit=["deg", "deg"])
    center = [center.ra.value, center.dec.value]

    RA_chiron_agora, DEC_chiron_agora = "0 56 58.44", "+7 41 16.0" 
    RA_chiron, DEC_chiron = "00 57 04.59", "+07 41 58.3" 
    ra, dec = utils.hms_to_hours(RA_chiron)*15, utils.dms_to_degrees(DEC_chiron)  
    ra_next, dec_next = utils.hms_to_hours(RA_chiron_agora)*15, utils.dms_to_degrees(DEC_chiron_agora) 
    data = np.squeeze(data)

    # Create WCS
    wcs = w

    # Create a figure and WCSAxes
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(1, 1, 1, projection=wcs)

    # Display the image
    norm = simple_norm(data, 'linear', percent=99.5)
    ax.imshow(data, norm=norm, origin='lower', cmap="gray")

    # Plot stars and asteroids
    gaias = twirl.gaia_radecs(center, 0.18333, limit=8)
    gaias_pixel = np.array(SkyCoord(gaias, unit="deg").to_pixel(wcs)).T
    asteroid = np.array(SkyCoord(ra, dec, unit=["deg", "deg"]).to_pixel(wcs)).T
    asteroid_next = np.array(SkyCoord(ra_next, dec_next, unit=["deg", "deg"]).to_pixel(wcs)).T

    ax.plot(*gaias_pixel.T, "o", fillstyle="none", c="C1", ms=18)
    ax.plot(*asteroid.T, "o", fillstyle="none", ms=18, color="green")
    ax.plot(*asteroid_next.T, "o", fillstyle="none", ms=18, color="red")

    # Annotate the asteroids
    ax.annotate("Chiron", asteroid, color='white', xytext=(asteroid[0]+50, asteroid[1]+40),
                bbox=dict(boxstyle="round", alpha=0.4, color='green'), fontsize=10)
    ax.annotate("Chiron Next", asteroid_next, color='white', xytext=(asteroid_next[0]+50, asteroid_next[1]+40),
                bbox=dict(boxstyle="round", alpha=0.4, color='red'), fontsize=10)

    # Set axis labels
    ax.set_xlabel('Right Ascension (J2000)')
    ax.set_ylabel('Declination (J2000)')

    plt.show()

ast = AstrometryNet()
ast.api_key = 'epepiumlsedbruam' 

path_img = r'Chiron_6_Clear.fits'
filename = get_pkg_data_filename(path_img)
hdu = fits.open(filename)[0]
ra_center, dec_center, w = solve(ast, path_img)
plot_fits(hdu, ra_center, dec_center, w)
