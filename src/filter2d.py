import os

from astropy.io import fits
import numpy as np


def filter2d(loc, taskid, beam, cube, mosaic=False, overwrite=False):

    # Define some file names:
    # loc = taskid + "/"
    cube_name = 'HI_B0' + str(beam).zfill(2) + '_cube' + str(cube) + '_image'
    if mosaic:
        cube_name = taskid + '_HIcube' + str(cube) + '_image'

    filter2d_name = loc + cube_name + '_filtered-2d.fits'
    mask2d_name = loc + cube_name + '_sofiaFS_mask-2d.fits'
    filter3d_name = loc + cube_name + '_filtered.fits'

    # If the filtered-2d file doesn't exist and we have a template, make the filtered-2d:
    if os.path.isfile(mask2d_name) & (not os.path.isfile(filter2d_name)):

        print("[FILTER2D] Making {}".format(filter2d_name))
        # Get header as template:
        header = fits.getheader(mask2d_name)
        data = fits.getdata(filter3d_name)[1, :, :]
        if np.all(np.isnan(data)):
            data = fits.getdata(filter3d_name)[-10, :, :]
            if np.all(np.isnan(data)):
                data = fits.getdata(filter3d_name)[-500, :, :]

        # Assign a positive value to the filter and set everything else to nan:
        filter2d = np.full(data.shape, np.nan)
        filter2d[np.isnan(data)] = 1.

        hdu_new = fits.PrimaryHDU(data=filter2d, header=header)
        hdu_new.writeto(filter2d_name, overwrite=overwrite)

    elif os.path.isfile(filter2d_name):
        print("\t[FILTER2D] already exists".format(filter2d_name))

    return


# Task yet to be written (if it's actually useful)
def filter2dto3d(taskid, beam, cube, overwrite=False):
    loc = taskid + "/"
    cube_name = 'HI_B0' + str(beam).zfill(2) + '_cube' + str(cube) + '_image'

    filter2d_name = loc + cube_name + '_filtered-2d.fits'
    mask2d_name = loc + cube_name + '_sofiaFS_mask-2d.fits'
    filter3d_name = loc + cube_name + '_filtered.fits'

    return