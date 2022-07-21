from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.io import fits
from astropy.table import Table
import numpy as np
from PIL import Image, ImageFont, ImageDraw
from scipy.ndimage import generate_binary_structure, generic_filter
import matplotlib.pyplot as plt


def test_mask(value):
    test = np.any(value > 0)
    return test


foot = np.array(generate_binary_structure(2, 1), dtype=int)

###################################################################

parser = ArgumentParser(description="Create new moment maps for (cleaned!) line cubes for a given taskid, beam, cubes",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='190915041',
                    help='Specify the input taskid (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3',
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-s', '--sources', default='all',
                    help='Specify sources to flag if necessary.  Can specify range or list. (default: %(default)s).')

###################################################################

# Parse the arguments above
args = parser.parse_args()

# Range of cubes/sources to work on:
taskid = args.taskid
cubes = [int(c) for c in args.cubes.split(',')]
if args.sources == 'all':
    sources = [str(s + 1) for s in range(len(catalog))]
elif '-' in args.sources:
    src_range = args.sources.split('-')
    sources = [str(s + int(src_range[0])) for s in range(int(src_range[1]) - int(src_range[0]) + 1)]
else:
    sources = [str(s) for s in args.sources.split(',')]

mos_loc = 'mos_' + taskid + '/'
for c in cubes:
    old_filename = taskid + '_HIcube' + str(c) + '_image'
    filename = taskid + '_HIcube' + str(c) + '_clean_image'
    catalog = Table.read(mos_loc + filename + '_cat.txt', format='ascii', header_start=18)
    hdu_mask2d = fits.open(mos_loc + filename + '_mask-2d.fits')
    hdu_filter2d = fits.open(mos_loc + old_filename + '_filtered-2d.fits')

    for s in sources:
        cat = catalog[catalog['id'] == s]
        Xmin = cat["x_min"]
        Ymin = cat["y_min"]
        Xmax = cat["x_max"]
        Ymax = cat["y_max"]

        # Create subimage of the continuum filtering & 2d mask to raise flag if it affects source
        filter2d = hdu_filter2d[0].data[int(Xmin):int(Xmax) + 1, int(Ymin):int(Ymax) + 1]
        mask2d = hdu_mask2d[0].data[int(Xmin):int(Xmax) + 1, int(Ymin):int(Ymax) + 1]

        fig1 = plt.figure()
        ax1 = fig1.add_subplot(111)
        ax1.imshow(filter2d)
        fig1.savefig('blah.png')

        fig2 = plt.figure()
        ax2 = fig2.add_subplot(111)
        ax2.imshow(mask2d)
        fig2.savefig('blah2.png')

        result = generic_filter(filter2d, test_mask, footprint=foot)
        if np.sum(result * mask2d) > 0:
            print("\tSpatial filtering flag for source {}".format(cat['id']))

            combo_im_name = mos_loc + filename + "_figures/" + filename + '_' + s + '_combo.png'
            print(combo_im_name)
            combo_im = Image.open(combo_im_name)

            I1 = ImageDraw.Draw(combo_im)

            I1.text((170, 42), "!", font=ImageFont.truetype(font='Arial Unicode.ttf', size=48), fill=(255, 0, 0))
            I1.text((170+772, 42), "!", font=ImageFont.truetype(font='Arial Unicode.ttf', size=48), fill=(255, 0, 0))
            I1.text((170+772+870, 42), "!", font=ImageFont.truetype(font='Arial Unicode.ttf', size=48), fill=(255, 0, 0))
            I1.text((170+772+870+845, 42), "!", font=ImageFont.truetype(font='Arial Unicode.ttf', size=48), fill=(255, 0, 0))

            combo_im.save(combo_im_name)

    hdu_mask2d.close()
    hdu_filter2d.close()
