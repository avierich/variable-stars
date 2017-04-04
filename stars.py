from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import CircularAperture
from photutils import CircularAnnulus
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

def buildAperature(wcs, skyCoord) :
    starX, starY = skycoord_to_pixel(skyCoord, wcs)
    ap_pos = [(starX, starY)]
    return [CircularAperture(ap_pos, r =10.), CircularAnnulus(ap_pos, r_in = 15, r_out = 25)]

def plotSky(filename, skyCoords) :
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)
    # Plot the star
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8], projection = wcs)
    tr_fk5 = ax.get_transform("fk5")
    lon = ax.coords[0]
    lat = ax.coords[1]
    lat.set_major_formatter('dd:mm:ss.s')
    lon.set_major_formatter('hh:mm:ss.s')
    plt.imshow(hdu.data, origin='lower', cmap='gray', norm = LogNorm())
    plt.xlabel('RA')
    plt.ylabel('Dec')
   
    apertures = []
    for skyCoord in skyCoords :
        apertures.extend(buildAperature(wcs, skyCoord))

    for app in apertures :
        app.plot()

def getStarData(filename, skyCoords) :
    hdu = fits.open(filename)[0]
    wcs = WCS(hdu.header)

    jd = hdu.header['JD']
  
    instMag = []
    for skyCoord in skyCoords :
        aperatures = buildAperature(wcs, skyCoord)
        phot_table = aperture_photometry(hdu.data, aperatures)
        # Calculate the instrumental magnitude
        instMag.append(-2.5 * np.log(phot_table[0][3]/aperatures[0].area() - phot_table[0][4]/aperatures[1].area()))

    return jd, instMag

def plotDir(directory, skyCoord, plotFirst = False) :
    vardata = []
    vartim = []
    for filename in os.listdir(directory):
        if filename.endswith(".fit"): 
#            print('Reading : ' + filename)
            if plotFirst :
                plotFirst = False
                plotSky(directory + filename, skyCoord)
            time, amplitude = getStarData(directory + filename, skyCoord)
            vardata.append(amplitude)
            vartim.append(time)
            continue
        else:
            continue

    return vartim, vardata

dy_peg = SkyCoord(90.7142 * u.degree,-39.154 * u.degree, Galactic)
GSC_01712_00542 = SkyCoord(90.6931 * u.degree,-39.1854 * u.degree, Galactic)
stars = [dy_peg, GSC_01712_00542]
vtime, vamp = plotDir('dy-peg/20150910-pt/V/', stars, plotFirst = True)
btime, bamp =  plotDir('dy-peg/20150910-pt/B/', stars)
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
