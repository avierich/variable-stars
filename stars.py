from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import CircularAperture
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates import ICRS, Galactic, FK4, FK5  # Low-level frames
from astropy.coordinates import Angle, Latitude, Longitude  # Angles
from astropy.utils.data import get_pkg_data_filename
from photutils import aperture_photometry
from astropy.wcs.utils import skycoord_to_pixel
import os
import numpy as np

def getStarData(filename, skyCoord) :
        hdu = fits.open(filename)[0]
        wcs = WCS(hdu.header)

        jd = hdu.header['JD']
        
        starX, starY = skycoord_to_pixel(skyCoord, wcs)
        ap_pos = [(starX, starY)]
        apertures = CircularAperture(ap_pos, r =10.)
        apertures.plot()
        phot_table = aperture_photometry(hdu.data, apertures)
        amplitude = phot_table[0][3]

        return jd, amplitude

def plotDir(directory, skyCoord) :
    vardata = []
    vartim = []
    for filename in os.listdir(directory):
        if filename.endswith(".fit"): 
            print('Reading : ' + filename)
            time, amplitude = getStarData(directory + filename, skyCoord)
            vardata.append(amplitude)
            vartim.append(time)
            continue
        else:
            continue

    return vartim, vardata

dy_peg = SkyCoord(90.7142 * u.degree,-39.154 * u.degree, Galactic)
vtime, vamp = plotDir('dy-peg/20150910-pt/V/', dy_peg)
btime, bamp =  plotDir('dy-peg/20150910-pt/B/', dy_peg)
vtime = np.array(vtime)
btime = np.array(btime)
vtime -= min(vtime)
btime -= min(btime)

fig = plt.figure()
plt.plot(vtime, vamp, 'x', label = 'V')
plt.plot(btime,bamp, 'x', label = 'B')
plt.legend()
plt.ylabel('Amplitude')
plt.xlabel('Time (Days)')
plt.show()
