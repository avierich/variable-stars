from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import CircularAperture
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_filename
from photutils import aperture_photometry
import os


star_path = 'dy-peg/20150910-PT/B/aligned_Calibrated-T16-jbotte-DY Peg-20150910-020218-B-BIN1-W-030-001.fit'
filename = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
#star_path = filename

# Get WCS corrds
hdu = fits.open(star_path)[0]
wcs = WCS(hdu.header)

firstday = hdu.header['JD']

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

# DY peg from pixels
coordDeg = wcs.all_pix2world([[2012,1368]],1)[0]
print(coordDeg)
c = SkyCoord(ra=coordDeg[0]*u.degree, dec=coordDeg[1]*u.degree, frame='fk5')
print(c.ra.hms)
print(c.dec.dms)


ap_pos = [(2012,1368)]
apertures = CircularAperture(ap_pos, r =10.)
phot_table = aperture_photometry(hdu.data, apertures)

def plotDir(directory):
    vardata = []
    vartim = []
    for filename in os.listdir(directory):
        if filename.endswith(".fit"): 
            # Get WCS corrds
            hdu = fits.open(directory+filename)[0]
            wcs = WCS(hdu.header)

            vartim.append(hdu.header['JD'] - firstday)
            

            # DY peg from pixels
            coordDeg = wcs.all_pix2world([[2012,1368]],1)[0]
            c = SkyCoord(ra=coordDeg[0]*u.degree, dec=coordDeg[1]*u.degree, frame='fk5')

            ap_pos = [(2012,1368)]
            apertures = CircularAperture(ap_pos, r =10.)
            apertures.plot()
            phot_table = aperture_photometry(hdu.data, apertures)
            vardata.append(phot_table[0][3])
            # print(os.path.join(directory, filename))
            continue
        else:
            continue

    return vartim, vardata

vtime, vamp = plotDir('dy-peg/20150910-pt/V/')
btime, bamp =  plotDir('dy-peg/20150910-pt/B/')

fig = plt.figure()
plt.plot(vtime, vamp, label = 'V')
plt.plot(btime,bamp, label = 'B')
plt.legend()
plt.ylabel('Amplitude')
plt.xlabel('Time (Days)')
plt.show()
