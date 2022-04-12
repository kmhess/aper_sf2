from argparse import ArgumentParser, RawTextHelpFormatter
import time as testtime

from astropy.io import ascii, fits
import numpy as np


###################################################################

def parse_args():
    parser = ArgumentParser(description="Regrid mosaic mask to an apertif beam. Assumes spectral dimensions are same.",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-l', '--loc', default='output/',
                        help='Specify the input working directory relative to where you\'re running.'
                             ' Important for mosaic. (default: %(default)s).')

    parser.add_argument('-t', '--taskid', default='190915041',
                        help='Specify the input taskid (default: %(default)s).')

    parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                        help='Specify the cubes on which to do source finding (default: %(default)s).')

    parser.add_argument('-s', '--sources', default='all',
                        help='Specify the sources included in the binary mask.'
                             ' (default: %(default)s).')

    # Parse the arguments above
    args = parser.parse_args()
    return args

###################################################################


def binary_mask(taskid, cubes, sources):

    loc = 'mos_' + taskid + '/'

    for c in cubes:
        cube_name = taskid + '_HIcube' + str(c) + '_image'
        catalog_file = loc + cube_name + '_sofiaFS_cat.txt'
        catalog = ascii.read(catalog_file, header_start=18)

        if sources == 'all':
            sources = [int(s + 1) for s in range(len(catalog))]
        elif '-' in sources:
            mask_range = sources.split('-')
            sources = [int(s + int(mask_range[0])) for s in range(int(mask_range[1]) - int(mask_range[0]) + 1)]
        else:
            sources = [int(s) for s in sources.split(',')]

        mask_file = loc + cube_name + '_sofiaFS_mask.fits'
        bin_mask_file = loc + cube_name + '_sofiaFS_mask_bin.fits'
        print("[BINARY_MASK] Making {} binary mask including requested sources: {}".format(taskid, sources))
        tic1 = testtime.perf_counter()
        mask = fits.open(mask_file)
        mask_lin = mask[0].data.ravel()
        bin_mask = np.array([1 if m in sources else 0 for m in mask_lin])
        bin_mask = binary.reshape(mask[0].data.shape)
        bin_hdu = fits.PrimaryHDU(data=bin_mask, header=mask[0].header)
        bin_hdu.writeto(bin_mask_file)
        mask.close()
        toc1 = testtime.perf_counter()
        print(f"Do binary mask: {toc1 - tic1:0.4f} seconds")

    return


if __name__ == '__main__':
    arguments = parse_args()
    taskid = arguments.taskid
    cubes = arguments.cubes

    loc = 'mos_' + taskid + '/'                               # Enable while using snakemake
    binary_mask(taskid, cubes, sources=arguments.sources)
