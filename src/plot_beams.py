import os.path
from argparse import ArgumentParser, RawTextHelpFormatter

from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np


###################################################################

# def parse_args():
parser = ArgumentParser(description="Make a 2d binary mask including just the objects to be cleaned."
                                    " Useful to detect which cubes should be regridded & cleaned.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='190915041',
                    help='Specify the input taskid or field (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-w', '--smooth', action='store_true',
                    help='If option is included apply edits to the smoothed SIP images.')

parser.add_argument('-d', '--directory', default='', required=False,
                    help='Specify the directory where taskid/field folders live containing the data (default: %(default)s).')

# Parse the arguments above
args = parser.parse_args()

###################################################################

bmaj = []
bmaj_med = []
bmin = []
bmin_med = []
for b in range(40):
    filename = args.directory + '/' + args.taskid + '/HI_B0'+str(b).zfill(2)+'_cube' + args.cubes + '_spline_clean_image.fits'
    if args.smooth:
        filename = args.directory + '/' + args.taskid + '/HI_B0'+str(b).zfill(2)+'_cube' + args.cubes + '_spline_clean_smooth_image.fits'
    if os.path.isfile(filename):
        file = fits.open(filename)
        if file[0].header.get('BMAJ'):
            bmaj = np.concatenate((bmaj, file[1].data['BMAJ']))
            bmaj_med = np.concatenate((bmaj_med, [file[0].header['BMAJ']]))
            bmin = np.concatenate((bmin, file[1].data['BMIN']))
            bmin_med = np.concatenate((bmin_med, [file[0].header['BMIN']]))
        file.close()
bmaj = bmaj * 3600
bmaj_med = bmaj_med * 3600
bmin = bmin * 3600
bmin_med = bmin_med * 3600

fig = plt.figure(figsize=(5, 7))
ax1 = fig.add_subplot(211)
ax1.hist(bmaj, bins=range(int(np.floor(np.min(bmaj))), int(np.ceil(np.max(bmaj)))+1),
         label='{:.1f}'.format(np.median(bmaj)))
ax1b = ax1.twinx()
ax1b.hist(bmaj_med, bins=range(int(np.floor(np.min(bmaj_med))), int(np.ceil(np.max(bmaj_med)))+1),
          histtype='step', color='orange', label='{:.1f}'.format(np.median(bmaj_med)))
yl, yh = ax1.get_ylim()
ax1.plot([np.median(bmaj), np.median(bmaj)], [yh*0.95, yh*1.05])
ax1.plot([np.median(bmaj_med), np.median(bmaj_med)], [yh*0.95, yh*1.05])
ax1.set_xlabel("Beam major axis [arcsec]")
ax1.set_ylabel("All cleaned channels")
ax1b.set_ylabel("Median per beam")
plt.title("Clean beam distribution for {}".format(args.taskid))
ax1.legend(loc=2)
ax1b.legend(loc=1)

ax2 = fig.add_subplot(212)
ax2.hist(bmin, bins=range(int(np.floor(np.min(bmin))), int(np.ceil(np.max(bmin)))+1),
         label='{:.1f}'.format(np.median(bmin)))
ax2b = ax2.twinx()
ax2b.hist(bmin_med, bins=range(int(np.floor(np.min(bmin_med))), int(np.ceil(np.max(bmin_med)))+1),
          histtype='step', color='orange', label='{:.1f}'.format(np.median(bmin_med)))
ym, yh = ax2.get_ylim()
ax2.plot([np.median(bmin), np.median(bmin)], [yh*0.95, yh+1.05])
ax2.plot([np.median(bmin_med), np.median(bmin_med)], [yh*0.95, yh*1.05])
ax2.set_xlabel("Beam minor axis [arcsec]")
ax2.set_ylabel("All cleaned channels")
ax2b.set_ylabel("Median per beam")
ax2.legend(loc=2)
ax2b.legend(loc=1)

outfile = args.directory + '/mos_' + args.taskid + '/' + args.taskid + '_beams_cube' + args.cubes[0] + '.png'
if args.smooth:
    outfile = args.directory + '/mos_' + args.taskid + '/' + args.taskid + '_beams_cube' + args.cubes[0] + '_smooth.png'    

fig.savefig(outfile, bbox_inches='tight')
