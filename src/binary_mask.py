from argparse import ArgumentParser, RawTextHelpFormatter
import time as testtime

from astropy.io import ascii, fits
import numpy as np


###################################################################

# def parse_args():
parser = ArgumentParser(description="Make a 2d binary mask including just the objects to be cleaned."
                                    " Useful to detect which cubes should be regridded & cleaned.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-l', '--loc', default='output/',
                    help='Specify the input working directory relative to where you\'re running.'
                         ' Important for mosaic. (default: %(default)s).')

parser.add_argument('-t', '--taskid', default='190915041',
                    help='Specify the input taskid or field (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-s', '--sources', default='all',
                    help='Specify the sources included in the binary mask.'
                         ' (default: %(default)s).')

# Parse the arguments above
args = parser.parse_args()
    # return args

###################################################################


# To do: Instead of supplying sources, allow option to read in the edited catalog

# def binary_mask(taskid, cubes, sources, njobs):

taskid = args.taskid
cubes = args.cubes
sources = args.sources

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
    mask_file2d = loc + cube_name + '_sofiaFS_mask-2d.fits'
    bin_mask2d_file = loc + cube_name + '_sofiaFS_mask-2d_bin.fits'
    print("[BINARY_MASK] Making {} binary mask including requested sources: {}".format(taskid, sources))

    # Define what lines of sight need to be kept based on input sources numbers
    mask2d_hdu = fits.open(mask_file2d)
    mask2d = mask2d_hdu[0].data
    xx, yy = range(mask2d.shape[1]), range(mask2d.shape[0])
    x, y = np.meshgrid(xx, yy)
    x, y = x.ravel(), y.ravel()
    mask2d_lin = mask2d.ravel()
    # Generalize to allow for overlapping sources by keeping all lines of sight with any object
    # src2d = np.array([True if m in sources else False for m in mask2d_lin])
    src2d = np.array([True if m > 0 else False for m in mask2d_lin])
    x, y, mask2d_lin = x[src2d], y[src2d], mask2d_lin[src2d]

    tic1 = testtime.perf_counter()
    mask_hdu = fits.open(mask_file)
    mask = mask_hdu[0].data
    bin_mask2d = np.zeros(mask2d.shape)

    for s in range(len(x)):
        mask_lin = mask[:, y[s], x[s]]
        bin_mask_lin = np.array([1 if m in sources else 0 for m in mask_lin])
        bin_mask2d[y[s], x[s]] = np.nanmax(bin_mask_lin)

    toc1 = testtime.perf_counter()
    print(f"Do binary mask: {toc1 - tic1:0.4f} seconds")

    bin_hdu = fits.PrimaryHDU(data=bin_mask2d, header=mask2d_hdu[0].header)
    bin_hdu.writeto(bin_mask2d_file, overwrite=True)

    mask2d_hdu.close()
    mask_hdu.close()

# if __name__ == '__main__':
#     arguments = parse_args()
#     taskid = arguments.taskid
#     cubes = arguments.cubes
#     njobs = arguments.njobs
#
#     loc = 'mos_' + taskid + '/'                               # Enable while using snakemake
#     binary_mask(taskid, cubes, sources=arguments.sources, njobs=njobs)
