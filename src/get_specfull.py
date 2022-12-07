import os
import random
import string

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.io import fits, ascii
from astropy.table import Table
import astropy.units as u
from astropy.wcs import WCS
import numpy as np
from matplotlib import pyplot as plt


def chan2freq(channels, header):
    frequencies = (header['CDELT3'] * (channels - (header['CRPIX3'] - 1)) + header['CRVAL3']) * u.Hz
    return frequencies

###################################################################

parser = ArgumentParser(description="Create new moment maps for (cleaned!) line cubes for a given taskid, beam, cubes",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='190915041',
                    help='Specify the input taskid (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3',
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-s', '--sources', default='all',
                    help='Specify sources to flag if necessary.  Can specify range or list. (default: %(default)s).')

parser.add_argument('-p', '--prefix', default=None,
                    help='Specify sources to flag if necessary.  Can specify range or list. (default: %(default)s).')

###################################################################

# Parse the arguments above
args = parser.parse_args()


# Range of cubes/sources to work on:
taskid = args.taskid
cubes = [int(c) for c in args.cubes.split(',')]
# if args.prefix:
#     prefix = args.prefix + '_'
# else:
#     prefix = ''

mos_loc = 'mos_' + taskid + '/'
for c in cubes:
    if not args.prefix:
        filename = taskid + '_HIcube' + str(c) + '_clean_image'
    else:
        filename = taskid + '_HIcube' + str(c) + '_clean_smooth_image'
    print("[GET_SPECFULL] Reading in mosaicked {} field.".format(taskid))
    mosaic = fits.open(mos_loc + filename + '.fits')
    wcs_mos = WCS(mosaic[0].header).celestial
    channels = np.asarray(range(mosaic[0].data.shape[0]))
    frequency = chan2freq(channels, mosaic[0].header)
    catalog = Table.read(mos_loc + filename + '_cat.xml')
    if not os.path.isdir(mos_loc + filename + '_figures'):
        os.system(f'mkdir {mos_loc}{filename}_figures')

    if args.sources == 'all':
        sources = [str(s + 1) for s in range(len(catalog))]
    elif '-' in args.sources:
        src_range = args.sources.split('-')
        sources = [str(s + int(src_range[0])) for s in range(int(src_range[1]) - int(src_range[0]) + 1)]
    else:
        sources = [str(s) for s in args.sources.split(',')]

    for s in sources:
        # outfile = mos_loc + prefix + filename + '_figures/' + filename + '_' + str(s) + '_specfull.txt'
        outfile = mos_loc + filename + '_figures/' + filename + '_' + str(s) + '_specfull.txt'
        if not os.path.isfile(outfile):
            src_hdu = fits.open(mos_loc + filename + '_cubelets/' + filename + '_' + str(s) + '_mask.fits')
            mask2d = np.sum(src_hdu[0].data, axis=0)
            wcs_src = WCS(src_hdu[0].header).celestial
            src_hdu.close()

            radec1 = wcs_src.pixel_to_world(0, 0)
            radec2 = wcs_src.pixel_to_world(mask2d.shape[1]-1, mask2d.shape[0]-1)
            x1, y1 = wcs_mos.world_to_pixel(radec1)
            x2, y2 = wcs_mos.world_to_pixel(radec2)

            subcube = mosaic[0].data[:, int(np.rint(y1)):int(np.rint(y2))+1, int(np.rint(x1)):int(np.rint(x2))+1]
            spectrum = np.nansum(subcube[:, mask2d != 0], axis=1)
            n_pix = 0 * channels + np.sum(mask2d != 0)

            print("[GET_SPECFULL] Writing *_specfull.txt for source {}.".format(s))
            code = ''.join(random.choices(string.ascii_letters + string.digits, k=6))
            with open(f'temp{code}.txt', 'w') as f:
                f.write("# Integrated source spectrum with noise\n")
                f.write("# Creator: SoFiA-image-pipeline.py\n")  # %s\n#\n" % sofia_version_full)
                f.write("# \n")
                f.write("# The source spectrum, with noise, is calculated by integrating over\n")
                f.write("# the 2D mask of the source in every channel.  This means every row \n")
                f.write("# has the same number of contributing pixels and the noise is the same\n")
                f.write("# for every point.\n")
                f.write("# \n")

                ascii.write([channels, frequency, spectrum, n_pix], f'temp2{code}.txt', format='fixed_width_two_line',
                            names=['chan', 'freq', 'f_sum', 'n_pix'])

            os.system(f"cat temp{code}.txt temp2{code}.txt > {outfile}")
            os.system(f"rm temp{code}.txt temp2{code}.txt")
        else:
            print(f"\t{outfile} already exists.  Will not overwrite.")

    mosaic.close()
