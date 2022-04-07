from glob import glob
import os

from modules.functions import chan2freq

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.io import ascii
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS
import matplotlib.pyplot as plt
import numpy as np

from src.filter2d import filter2d


###################################################################

def parse_args():
    parser = ArgumentParser(description="Make plots of sourcefinding.py results.",
                            formatter_class=RawTextHelpFormatter)

    parser.add_argument('-l', '--loc', default='output/',
                        help='Specify the input working directory relative to where you\'re running.'
                             ' Important for mosaic. (default: %(default)s).')

    parser.add_argument('-t', '--taskid', default='190915041',
                        help='Specify the input taskid (default: %(default)s).')

    parser.add_argument('-b', '--beams', default='0-39',
                        help='Specify a range (0-39) or list (3,5,7,11) of beams on which to do source finding'
                             ' (default: %(default)s).')

    parser.add_argument('-n', "--nospline",
                        help="Only controls output name of png!  For book keeping purposes.",
                        action='store_true')

    parser.add_argument('-m', "--mosaic",
                        help="Flag to indicate if operating on a mosaic instead of individual beams.",
                        action='store_true')

    # Parse the arguments above
    args = parser.parse_args()
    return args


# If running by itself use: python -m src/checkmasks -t 191004041 -b 14
def main(loc, taskid, beam=[40], nospline=False, mosaic=False):

    cubes = [1, 2, 3]  # Most sources in 2; nearest galaxies in 3.
    max_cat_len = 35

    HI_restfreq = 1420405751.77 * u.Hz
    optical_HI = u.doppler_optical(HI_restfreq)

    colors = ['purple', 'blue', 'black']
    plt.close('all')
    for b in beam:
        # Define some file names and work space:
        source_per_beam = 0
        results = glob(loc + 'HI_B0{:02}_cube*sofia_mask-2d.fits'.format(b))
        if mosaic:
            results = glob(loc + taskid + '_HIcube*_image_sofia_mask-2d.fits')
        # Output where the program is looking exactly
        print(loc)

        if len(results) > 0:
            print('[CHECKMASKS] Making summary figures/spectra for beam {:02}'.format(b))
            header = fits.getheader(results[0])
            wcs = WCS(header)
            if mosaic:
                fig_im, ax_im = plt.subplots(1, 1, figsize=(15, 12), subplot_kw={'projection': wcs}, squeeze=False)
            else:
                fig_im, ax_im = plt.subplots(1, 3, figsize=(15, 5), subplot_kw={'projection': wcs})

            for c in cubes:
                cube_name = 'HI_image_cube' + str(c)
                if mosaic:
                    cube_name = taskid + '_HIcube' + str(c) + '_image'
                if os.path.isfile(loc + cube_name + '_sofia_cat.txt'):
                    cat = ascii.read(loc + cube_name + '_sofia_cat.txt')
                    if len(cat) > max_cat_len:
                        cat = cat[:max_cat_len]
                    source_per_beam += len(cat)
            fig_spec, ax_spec = plt.subplots(source_per_beam, 1, figsize=(15, 3*source_per_beam), squeeze=False)

            previous = 0
            for c in cubes:
                cube_name = 'HI_B0' + str(b).zfill(2) + '_cube' + str(c) + '_image'
                if mosaic:
                    cube_name = taskid + '_HIcube' + str(c) + '_image'
                # Make robust for saved data:
                print(loc + cube_name + '_filtered-2d.fits')
                if not os.path.isfile(loc + cube_name + '_filtered-2d.fits'):
                    print("\tTry making 2D plot of filtered pixels.")
                    filter2d(loc, taskid, b, c, mosaic=mosaic, overwrite=False)

                if os.path.isfile(loc + cube_name + '_filtered-2d.fits'):
                    print("\tPlotting from 2d filtered file.")
                    hdu_filter = fits.open(loc + cube_name + '_filtered-2d.fits')
                    filter2d_im = hdu_filter[0].data[:, :]
                    filter2d_im[filter2d_im == 1] = 9.
                    hdu_filter.close()
                elif os.path.isfile(loc + cube_name + '_filtered.fits'):
                    print("\tPlotting from 3d filtered file.")
                    hdu_filter = fits.open(loc + cube_name + '_filtered.fits')
                    # SVC data has the 0th channel blanked; so use the first channel instead.
                    filter2d_im = hdu_filter[0].data[1, :, :]
                    filter2d_im[np.isnan(filter2d_im)] = 9.
                    filter2d_im[filter2d_im < 9] = np.nan
                    hdu_filter.close()
                # elif os.path.isfile(loc + cube_name + '_filtered_spline.fits'):
                #     print("\tShould not be here!")
                #     hdu_filter = fits.open(loc + cube_name + '_filtered_spline.fits')
                #     # SVC data has the 0th channel blanked; so use the first channel instead.
                #     filter2d_im = hdu_filter[0].data[1, :, :]
                #     filter2d_im[np.isnan(filter2d_im)] = 9.
                #     filter2d_im[filter2d_im < 9] = np.nan
                #     hdu_filter.close()
                else:
                    print("\tNo continuum filtered file for Beam {:02} Cube {}. Check sourcefinding/ALTA?".format(b, c))
                    # As a precaution, if there is no filtered cube, don't plot source masks either.
                    continue

                if os.path.isfile(loc + cube_name + '_sofia_cat.txt'):
                    cat = ascii.read(loc + cube_name + '_sofia_cat.txt')
                    print("\tFound {} sources in Beam {:02} Cube {}".format(len(cat), b, c))
                    if len(cat) > max_cat_len:
                        print("\tMore than {} candidates: seems this cube is crap".format(max_cat_len))
                        cat = cat[:max_cat_len]
                    hdu_mask = fits.open(loc + cube_name + '_sofia_mask-2d.fits')
                    mask2d = hdu_mask[0].data[:, :]

                    mask2d = np.asfarray(mask2d)
                    mask2d[mask2d < 1] = np.nan

                    # if os.path.isfile(loc + cube_name + '_all_spline.fits'):
                    #     hdu_spline = fits.open(loc + cube_name + '_all_spline.fits')
                    # elif os.path.isfile(loc + cube_name + '_filtered_spline.fits'):
                    #     hdu_spline = fits.open(loc + cube_name + '_filtered_spline.fits')
                    hdu_spline = fits.open(loc + cube_name + '_spline.fits')
                    cube_frequencies = chan2freq(np.array(range(hdu_spline[0].data.shape[0])), hdu_header=hdu_spline[0].header)
                    optical_velocity = cube_frequencies.to(u.km/u.s, equivalencies=optical_HI)

                    a = c
                    if mosaic:
                        a = 1
                    ax_im[0, a-1].imshow(filter2d_im, cmap='Greys_r', vmax=10, vmin=8, origin='lower')
                    ax_im[0, a-1].imshow(mask2d, cmap='gist_rainbow', origin='lower')
                    ax_im[0, a-1].set_title("Beam {:02} Cube {}".format(b, c))
                    ra = ax_im[0, a - 1].coords[0]  # Don't understand why this doesn't work: python 2.7 vs 3?
                    ra.set_format_unit(u.hour)
                    for s in range(len(cat)):
                        ax_im[0, a-1].text(cat['col3'][s] + np.random.uniform(-40, 40),
                                           cat['col4'][s] + np.random.uniform(-40, 40),
                                           cat['col2'][s], color='black', fontsize=16)
                        # Do spectrum sums faster on subcubes:
                        subcube = hdu_spline[0].data[:, int(cat['col8'][s]):int(cat['col9'][s])+1,
                                                        int(cat['col6'][s]):int(cat['col7'][s])+1]
                        submask = mask2d[int(cat['col8'][s]):int(cat['col9'][s])+1,
                                         int(cat['col6'][s]):int(cat['col7'][s])+1]
                        # spectrum = np.nansum(hdu_spline[0].data[:, mask2d == cat['col2'][s]], axis=1)
                        spectrum = np.nansum(subcube[:, submask == cat['col2'][s]], axis=1)
                        maskmin = chan2freq(cat['col10'][s],
                                            hdu_header=hdu_spline[0].header).to(u.km/u.s, equivalencies=optical_HI).value
                        maskmax = chan2freq(cat['col11'][s],
                                            hdu_header=hdu_spline[0].header).to(u.km/u.s, equivalencies=optical_HI).value
                        ax_spec[previous + s, 0].plot([optical_velocity[-1].value, optical_velocity[0].value],
                                                      [0, 0], '--', color='gray')
                        ax_spec[previous + s, 0].plot(optical_velocity, spectrum, c=colors[c-1])
                        ax_spec[previous + s, 0].plot([maskmin, maskmin], [np.nanmin(spectrum), np.nanmax(spectrum)],
                                                      ':', color='gray')
                        ax_spec[previous + s, 0].plot([maskmax, maskmax], [np.nanmin(spectrum), np.nanmax(spectrum)],
                                                      ':', color='gray')
                        ax_spec[previous + s, 0].set_title("Beam {:02}, Cube {}, Source {}".format(b, c,
                                                                                                   cat['col2'][s]))
                        ax_spec[previous + s, 0].set_xlim(optical_velocity[-1].value, optical_velocity[0].value)
                        if (np.nanmax(spectrum) > 2.) | (np.nanmin(spectrum) < -2.):
                            ax_spec[previous + s, 0].set_ylim(np.nanmax(spectrum[cat['col10'][s]:cat['col11'][s]]) * -2,
                                                              np.nanmax(spectrum[cat['col10'][s]:cat['col11'][s]]) * 2)
                        ax_spec[previous + s, 0].set_ylabel("Integrated Flux")
                        if previous + s == source_per_beam - 1:
                            ax_spec[previous + s, 0].set_xlabel("Optical Velocity [km/s]")
                        if s+1 >= max_cat_len:
                            ax_spec[previous + s, 0].text(0.5, 0.05, "Too many spectra to plot, check by hand.",
                                                          ha='center', transform=ax_spec[previous + s, 0].transAxes)
                    previous += len(cat)

                    hdu_mask.close()
                    hdu_spline.close()

                else:
                    print("\tNO sources in Beam {:02} Cube {}".format(b, c))
                    ax_im[c - 1].imshow(filter2d_im, cmap='Greys_r', vmax=10, vmin=8, origin='lower')
                    ax_im[c - 1].set_title("Beam {:02} Cube {}".format(b, c))
                    ra = ax_im[c - 1].coords[0]  # Don't understand why this doesn't work: python 2.7 vs 3?
                    ra.set_format_unit(u.hour)
                    ax_im[c - 1].axis('off')

            if nospline:
                fig_im.savefig(loc + 'HI_image_sofia_summary.png', bbox_inches='tight')
                fig_spec.savefig(loc + 'HI_image_sofia_summary_spec.png', bbox_inches='tight')
            else:
                fig_im.savefig(loc + 'HI_image_sofia_summary_filtspline.png', bbox_inches='tight')
                fig_spec.savefig(loc + 'HI_image_sofia_summary_spec_filtspline.png', bbox_inches='tight')
            plt.close(fig_im)
            plt.close(fig_spec)
        else:
            print("\tNO sources found in any cube for Beam {:02}".format(b))


if __name__ == '__main__':
    arguments = parse_args()
    # Range of beams to work on:
    taskid = arguments.taskid
    mosaic = arguments.mosaic
    loc = arguments.loc
    if '-' in arguments.beams:
        b_range = arguments.beams.split('-')
        beams = np.array(range(int(b_range[1])-int(b_range[0])+1)) + int(b_range[0])
    else:
        beams = [int(b) for b in arguments.beams.split(',')]
    if loc[-1] != '/':
        loc = loc + '/'

    main(loc, taskid, beams, mosaic=mosaic)
