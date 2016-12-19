import sys
import os
import numpy as np
from astropy.io import fits

try: 
    infile = sys.argv[1]
except:
    sys.exit('Usage: {} FITSFILE'.format(sys.argv[0]))


basename = os.path.basename(infile)
baseroot, _ = os.path.splitext(basename)

hdu = fits.open(infile)[0]
if hdu.data is None:
    hdu = fits.open(infile)[1]
hdr = hdu.header

ny, nx = hdu.data.shape

# Size of chunks
mx, my = 290, 290
xchunks, ychunks = nx//mx, ny//my

# Initialize 3-d array to hold stack of tiles
ntiles = xchunks*ychunks
tilestack = np.zeros((ntiles, my, mx))

# Little nudges to get the peaks to line up for the different y chunks
jshifts = {
    0: -4,
    1: +2,
    2: +3,
    3: +2,
    4: -4,
}

# Put each tile on to the stack
ktile = 0
for jchunk in range(ychunks):
    yslice = slice(jchunk*my, jchunk*my + my)
    jshift = jshifts.get(jchunk, 0)
    for ichunk in range(xchunks):
        xslice = slice(ichunk*mx, ichunk*mx + mx)
        # Roll the data according to the jshift
        tile = np.roll(hdu.data[yslice, xslice], jshift, axis=0)
        # This has side effect of making tile be a copy rather than a view

        # Normalize each tile and dd to the stack
        tilestack[ktile, :, :] = tile / np.nanmedian(tile)
        ktile += 1

# Take the median down the stack to get the pattern
pattern_tile = np.nanmean(tilestack, axis=0)

# Now use copies of pattern_tile to make the pattern noise map
patmap = np.ones_like(hdu.data)
for jchunk in range(ychunks):
    yslice = slice(jchunk*my, jchunk*my + my)
    jshift = jshifts.get(jchunk, 0)
    for ichunk in range(xchunks):
        xslice = slice(ichunk*mx, ichunk*mx + mx)
        # Roll the data back again according to the jshift
        patmap[yslice, xslice] = np.roll(pattern_tile, -jshift, axis=0)

fits.PrimaryHDU(data=patmap,
                header=hdr).writeto(infile.replace('.fits',
                                                   '-pattern.fits'), clobber=True)

fits.PrimaryHDU(data=hdu.data/patmap,
                header=hdr).writeto(infile.replace('.fits',
                                                   '-patfix.fits'), clobber=True)
