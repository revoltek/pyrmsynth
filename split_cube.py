#!/usr/bin/env python3

from pathlib import Path
from astropy.io import fits
import numpy as np

# === CONFIGURATION ===
input_fits = "image_full_low_QU.cube.dirty.fits.fz"          # your main cube
output_dir = Path("split_images")    # output folder
output_dir.mkdir(exist_ok=True)

# === OPEN INPUT CUBE ===
with fits.open(input_fits) as hdul:
    hdul.info()              # shows structure and compressed HDUs
    data = hdul[1].data      # automatically decompressed
    header = hdul[1].header

    # Expect shape (Nfreq, Nstokes, Ndec, Nra)
    print("Data shape:", data.shape)
    n_freq, n_stokes, n_dec, n_ra = data.shape

    for i_f in range(n_freq):
        for i_s in range(n_stokes):
            # Extract one (RA,DEC) plane, but keep freq and stokes as degenerate
            img = data[i_f, i_s, :, :]
            img = img[np.newaxis, np.newaxis, :, :]   # shape = (1,1,Ndec,Nra)

            # Copy header
            hdu_header = header.copy()
            hdu_header['NAXIS'] = 4
            hdu_header['NAXIS1'] = n_ra
            hdu_header['NAXIS2'] = n_dec
            hdu_header['NAXIS3'] = 1
            hdu_header['NAXIS4'] = 1

            # Adjust reference pixels so that this slice corresponds to its freq/stokes
            # (keeping WCS consistent)
            hdu_header['CRPIX3'] = 1.0
            hdu_header['CRPIX4'] = 1.0
            hdu_header['CRVAL3'] = header['CRVAL3'] + (i_s) * header.get('CDELT3', 1)
            hdu_header['CRVAL4'] = header['CRVAL4'] + (i_f) * header.get('CDELT4', 1)

            # Create HDU
            hdu = fits.PrimaryHDU(data=img, header=hdu_header)

            # Write file
            outfile = output_dir / f"image_stk{i_s:02d}_freq{i_f:03d}.fits"
            hdu.writeto(outfile, overwrite=True)
            print(f"Saved {outfile}")
