import os

from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.io import fits
from astropy.table import Table

###################################################################

# def parse_args():
parser = ArgumentParser(description="Make a 2d binary mask including just the objects to be cleaned."
                                    " Useful to detect which cubes should be regridded & cleaned.",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-f', '--mosaic_file', default='190915041', required=True,
                    help='Specify the mosaic file to which beam info is saved in extensions (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='2', required=True,
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-d', '--directory', default='', required=False,
                    help='Specify the directory where taskid/field folders live containing the data (default: %(default)s).')

###################################################################

# Parse the arguments above
args = parser.parse_args()

mosaic = args.mosaic_file
mosaic_hdu = fits.open(mosaic, mode='update')

field = mosaic.split('/')[-1].split('_')[0] 

print("[BEAMS2MOSAIC] Adding restoring beam properties to mosaic file {}".format(mosaic))
for b in range(40):

    beam_name = str(b).zfill(2)
    
    #trying .txt files first
    end_name = '.txt'
    if 'smooth' in mosaic:
        end_name = '_smooth.txt'
    filename = args.directory + '/' + field + '/HI_B0' + beam_name + '_cube' + args.cubes + end_name
    print('\t',filename)
    if os.path.isfile(filename):
        file = Table.read(filename, format='ascii')
        try:
            if mosaic_hdu['BEAM{}'.format(beam_name)]:
                print("[BEAMS2MOSAIC] BEAM"+str(b).zfill(2)+" extension table already exists.")
                pass
        except KeyError:
            print("[BEAMS2MOSAIC] Adding channel clean beam properties to BEAM"+str(b).zfill(2)+" extension table.")
            col1 = fits.Column(name='BMAJ', format='1E', unit='deg', array=file['bmaj'])
            col2 = fits.Column(name='BMIN', format='1E', unit='deg', array=file['bmin'])
            col3 = fits.Column(name='BPA', format='1E', unit='deg', array=file['bpa'])
            col4 = fits.Column(name='CHAN', format='1J', array=file['chan'])
            beam_hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
            beam_hdu.name = 'BEAM{}'.format(beam_name)
            beam_hdu.header.comments['NAXIS2'] = 'number of channels'
            mosaic_hdu.append(beam_hdu)
    else:
        end_name = 'clean_image.fits'
        if 'smooth' in mosaic:
            end_name = 'clean_smooth_image.fits'
        filename = args.directory + '/' + field + '/HI_B0' + beam_name + '_cube' + args.cubes + '_spline_' + end_name
        print('\t',filename)
        if os.path.isfile(filename):
            file = fits.open(filename)
            try:
                if mosaic_hdu['BEAM{}'.format(beam_name)]:
                    print("[BEAMS2MOSAIC] BEAM"+str(b).zfill(2)+" extension table already exists.")
                    pass
            except KeyError:
                print("[BEAMS2MOSAIC] Adding channel clean beam properties to BEAM"+str(b).zfill(2)+" extension table.")
                col1 = fits.Column(name='BMAJ', format='1E', unit='deg', array=file[1].data['BMAJ'])
                col2 = fits.Column(name='BMIN', format='1E', unit='deg', array=file[1].data['BMIN'])
                col3 = fits.Column(name='BPA', format='1E', unit='deg', array=file[1].data['BPA'])
                col4 = fits.Column(name='CHAN', format='1J', array=file[1].data['CHAN'])
                beam_hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
                beam_hdu.name = 'BEAM{}'.format(beam_name)
                beam_hdu.header.comments['NAXIS2'] = 'number of channels'
                mosaic_hdu.append(beam_hdu)
            file.close()
        else:
            print("[BEAMS2MOSAIC] Beam"+str(b).zfill(2)+" image does not exist.")

mosaic_hdu.flush()
mosaic_hdu.close()
