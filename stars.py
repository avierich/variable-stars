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
    return [CircularAperture(ap_pos, r =3.), CircularAnnulus(ap_pos, r_in = 6, r_out = 8)]

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
        phot_table = aperture_photometry(hdu.data, aperatures, error=15*2.2)

        # Calculate the instrumental magnitude
        S = phot_table[0][3]
        errS = phot_table[0][4] 
        B = phot_table[0][5]
        errB = phot_table[0][6]
        areas = aperatures[0].area()/aperatures[1].area()
        instMag.append([-2.5*np.log10(S - areas*B),
                        (2.5*errS/(areas*B - S))**2 + (2.5*areas*errB/(S - areas*B))**2])

    return jd, instMag

def plotDir(directory, skyCoord, plotFirst = False) :
    vardata = []
    vartim = []
    for filename in os.listdir(directory):
        if filename.endswith(".fit"): 
            print('Reading : ' + filename)
            if plotFirst :
                plotFirst = False
                plotSky(directory + filename, skyCoord)
            time, amplitude = getStarData(directory + filename, skyCoord)
            vardata.append(amplitude)
            vartim.append(time)
            continue
        else:
            continue

    return np.array(vartim), np.array(vardata)

def getOffset(starData, knownMag) :
    offsets = []
    for frame in starData :
        starOffset = []
        for i in range(len(knownMag)) :
            starOffset.append(knownMag[i] - frame[i])
        offsets.append(starOffset)

    offsets = np.array(offsets)
    avgOffset = []
    for frame in offsets :
        avgOffset.append([frame.mean(),frame.std()])

    return np.array(avgOffset)



dy_peg = SkyCoord(90.7142 * u.degree,-39.154 * u.degree, Galactic)
GSC_01712_00542 = SkyCoord(90.6931 * u.degree,-39.1854 * u.degree, Galactic)
HD_218587 = SkyCoord(90.7271 * u.degree,-39.2479 * u.degree, Galactic)
GSC_01712_01246 = SkyCoord(90.6167 * u.degree,-39.2009 * u.degree, Galactic)
stars = [dy_peg, GSC_01712_00542, HD_218587, GSC_01712_01246]
vtime, vamp = plotDir('dy-peg/20150910-pt/V/', stars, plotFirst = True)
btime, bamp =  plotDir('dy-peg/20150910-pt/B/', stars)
vtime -= min(vtime)
btime -= min(btime)


offsetB = getOffset(bamp[:,1:,0], [12.5,10.42,13])
offsetV = getOffset(vamp[:,1:,0], [11.7,9.91,11.1])

errV = offsetV[:,1]
errB = offsetB[:,1]

fig = plt.figure()
plt.errorbar(vtime, vamp[:,0,0] + offsetV[:,0],marker = 'o', xerr = None, yerr = errV,  label = 'V',fmt = 'none')
plt.errorbar(btime,bamp[:,0,0] + offsetB[:,0],marker = 'o', xerr = None, yerr =  errB,  label = 'B',fmt = 'none')
plt.legend()
plt.ylabel('Amplitude')
plt.xlabel('Time (Days)')
plt.show()
