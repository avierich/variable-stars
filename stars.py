from astropy.io import fits
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from photutils import CircularAperture
from astropy.wcs import WCS
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.utils.data import get_pkg_data_filename


star_path = 'dy-peg/20150910-PT/B/aligned_Calibrated-T16-jbotte-DY Peg-20150910-020218-B-BIN1-W-030-001.fit'
filename = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')
#star_path = filename

# Get WCS corrds
hdu = fits.open(star_path)[0]
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

# DY peg from pixels
coordDeg = wcs.all_pix2world([[2012,1368]],1)[0]
print(coordDeg)
c = SkyCoord(ra=coordDeg[0]*u.degree, dec=coordDeg[1]*u.degree, frame='fk5')
print(c.ra.hms)
print(c.dec.dms)


ap_pos = [(2012,1368)]
aperatures = CircularAperture(ap_pos, r =10.)
aperatures.plot()

plt.show()
