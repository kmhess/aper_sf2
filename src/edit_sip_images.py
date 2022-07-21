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

mos_loc = 'mos_' + taskid + '/'
for c in cubes:
    old_filename = taskid + '_HIcube' + str(c) + '_image'
    filename = taskid + '_HIcube' + str(c) + '_clean_image'
    catalog = Table.read(mos_loc + filename + '_cat.txt', format='ascii', header_start=18)
    hdu_mask2d = fits.open(mos_loc + filename + '_mask-2d.fits')
    hdu_filter2d = fits.open(mos_loc + old_filename + '_filtered-2d.fits')
    cubeDim = hdu_filter2d[0].data.shape

    if args.sources == 'all':
        sources = [str(s + 1) for s in range(len(catalog))]
    elif '-' in args.sources:
        src_range = args.sources.split('-')
        sources = [str(s + int(src_range[0])) for s in range(int(src_range[1]) - int(src_range[0]) + 1)]
    else:
        sources = [str(s) for s in args.sources.split(',')]

    for s in sources:
        cat = catalog[catalog['id'] == int(s)]
        Xc = cat["x"]
        Yc = cat["y"]
        Xmin = cat["x_min"]
        Ymin = cat["y_min"]
        Xmax = cat["x_max"]
        Ymax = cat["y_max"]
        cPixXNew = int(Xc)
        cPixYNew = int(Yc)
        maxX = 2 * max(abs(cPixXNew - Xmin), abs(cPixXNew - Xmax))
        maxY = 2 * max(abs(cPixYNew - Ymin), abs(cPixYNew - Ymax))
        XminNew = cPixXNew - maxX
        if XminNew < 0: XminNew = 0
        YminNew = cPixYNew - maxY
        if YminNew < 0: YminNew = 0
        XmaxNew = cPixXNew + maxX
        if XmaxNew > cubeDim[1] - 1: XmaxNew = cubeDim[1] - 1
        YmaxNew = cPixYNew + maxY
        if YmaxNew > cubeDim[0] - 1: YmaxNew = cubeDim[0] - 1

        # Create subimage of the continuum filtering & 2d mask to raise flag if it affects source
        filter2d = hdu_filter2d[0].data[int(YminNew):int(YmaxNew) + 1, int(XminNew):int(XmaxNew) + 1]
        mask2d = hdu_mask2d[0].data[int(YminNew):int(YmaxNew) + 1, int(XminNew):int(XmaxNew) + 1]

        result = generic_filter(filter2d, test_mask, footprint=foot)

        if np.nansum(result * mask2d) > 0.:
            print("\tSpatial filtering flag for source {}".format(s))

            combo_im_name = mos_loc + filename + "_figures/" + filename + '_' + s + '_combo.png'
            combo_im = Image.open(combo_im_name)

            I1 = ImageDraw.Draw(combo_im)

            # if mac: font='Arial Unicode.ttf'
            # else:
            font = 'DejaVuSans.ttf'
            I1.text((0.051*combo_im.size[0], 0.037*combo_im.size[1]), "!", font=ImageFont.truetype(font=font, size=48),
                    fill=(255, 0, 0))
            I1.text((0.280*combo_im.size[0], 55*combo_im.size[1]), "!", font=ImageFont.truetype(font=font, size=48),
                    fill=(255, 0, 0))
            I1.text((0.539*combo_im.size[0], 55*combo_im.size[1]), "!", font=ImageFont.truetype(font=font, size=48),
                    fill=(255, 0, 0))
            I1.text((0.791*combo_im.size[0], 55*combo_im.size[1]), "!", font=ImageFont.truetype(font=font, size=48),
                    fill=(255, 0, 0))

            combo_im.save(combo_im_name)

    hdu_mask2d.close()
    hdu_filter2d.close()
