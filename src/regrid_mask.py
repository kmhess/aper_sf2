from argparse import ArgumentParser, RawTextHelpFormatter
import time as testtime

from astropy.io import fits
import numpy as np
from reproject import reproject_interp


###################################################################

def parse_args():
    parser = ArgumentParser(description="Regrid mosaic mask to an apertif beam. Assumes spectral dimensions are same.",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-t', '--taskid', default='190915041',
                        help='Specify the input taskid (default: %(default)s).')

    parser.add_argument('-b', '--beams', default='0-39',
                        help='Specify a range (0-39) or list (3,5,7,11) of beams on which to do source finding'
                             ' (default: %(default)s).')

    parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                        help='Specify the cubes on which to do source finding (default: %(default)s).')

    # Parse the arguments above
    args = parser.parse_args()
    return args

###################################################################


def main(taskid, beams, cubes):

    loc = 'mos_' + taskid + '/'

    for c in cubes:
        mask_file = loc + taskid + '_HIcube' + str(c) + '_image_sofiaFS_mask_bin.fits'
        mask = fits.open(mask_file)

        for b in beams:
            print("[REGRID_MASK] Regridding {} mosaic mask to cube {}, beam {}".format(taskid, c, b))
            image_file = loc[4:] + 'HI_B0' + str(b).zfill(2) + '_cube' + str(c) + '_image.fits'
            template = fits.getheader(image_file)

            tic1 = testtime.perf_counter()
            mask_reproj, footprint = reproject_interp(mask, template)
            toc1 = testtime.perf_counter()
            print(f"Do mask reprojection: {toc1 - tic1:0.4f} seconds")

            new_mask = fits.PrimaryHDU(data=mask_reproj, header=template)
            new_mask.writeto(image_file[:-5] + '_mask_bin.fits')

        mask.close()

    return


if __name__ == '__main__':
    arguments = parse_args()
    taskid = arguments.taskid
    # Range of beams to work on:
    if '-' in arguments.beams:
        b_range = arguments.beams.split('-')
        beams = np.array(range(int(b_range[1]) - int(b_range[0]) + 1)) + int(b_range[0])
    else:
        beams = [int(b) for b in arguments.beams.split(',')]
    cubes = [int(c) for c in arguments.cubes.split(',')]

    main(taskid, beams, cubes)
