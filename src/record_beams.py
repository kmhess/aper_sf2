import os.path
from argparse import ArgumentParser, RawTextHelpFormatter

from astropy.io import ascii,fits
from astropy.table import Table

###################################################################

# def parse_args():
parser = ArgumentParser(description="Save clean beam information into a text file per beam.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='M1315+3356', required=True,
                    help='Specify the input taskid or field (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

# Parse the arguments above
args = parser.parse_args()
# return args

###################################################################

for b in range(40)[:1]:
    filename = args.taskid + '/HI_B0'+str(b).zfill(2)+'_cube' + args.cubes + '_spline_clean_image.fits'
    targ_filename = args.taskid + '/HI_B0'+str(b).zfill(2)+'_cube' + args.cubes + '.txt'
    if os.path.isfile(filename):
        file = fits.open(filename)
        data = Table()
        data['chan'] = file[1].data['CHAN']
        data['bmaj'] = file[1].data['BMAJ']
        data['bmin'] = file[1].data['BMIN']
        data['bpa'] = file[1].data['BPA']
        file.close()
        if os.path.isfile(targ_filename):
            print(targ_filename)
            ascii.write(data, targ_filename, overwrite=True)
        # else:
        #     os.system('touch {}'.format(targ_filename))
        #     ascii.write(data, targ_filename, overwrite=True)

