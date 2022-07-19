from __future__ import print_function
import logging
import os

from modules.natural_cubic_spline import fspline
from src import trim_beam

from argparse import ArgumentParser, RawTextHelpFormatter
import astropy.io.fits as pyfits
import numpy as np

from multiprocessing import Queue, Process, cpu_count
from tqdm.auto import trange

from apercal.libs import lib
from apercal.subs import managefiles
import apercal


def worker(inQueue, outQueue):
    """
    Defines the worker process of the parallelisation with multiprocessing.Queue
    and multiprocessing.Process.
    """
    for i in iter(inQueue.get, 'STOP'):
        status = run(i)
        outQueue.put(status)


def run(i):
    global new_splinecube_data

    try:
        # Do the spline fitting on the z-axis to masked cube
        fit = fspline(np.linspace(1, orig_data.shape[0], orig_data.shape[0]),
                      np.nan_to_num(new_splinecube_data[:, x[i], y[i]]), k=5)
        new_splinecube_data[:, x[i], y[i]] = orig_data[:, x[i], y[i]] - fit
        return 'OK'

    except Exception:
        print("[ERROR] Something went wrong with the Spline fitting [" + str(i) + "]")
        return np.nan


def worker2(inQueue, outQueue):
    """
    Defines the worker process of the parallelisation with multiprocessing.Queue
    and multiprocessing.Process.
    """
    for i in iter(inQueue.get, 'STOP'):
        status = run2(i)
        outQueue.put(status)


def run2(i):
    global new_cleancube_data, new_modelcube_data, new_residualcube_data

    try:
        for name in ['map_{:02}_{:02}'.format(b, minc), 'beam_{:02}_{:02}'.format(b, minc), 'mask_{:02}_{:02}'.format(b, minc)]:
            imsub = lib.miriad('imsub')
            imsub.in_ = name
            imsub.out = name + "_" + str(chan[i]).zfill(4)
            imsub.region = '"images({})"'.format(chan[i] + 1)
            imsub.go()
    except RuntimeError:
        print("\tProblems reading in map/beam/mask data for channel {}.".format(chan[i]))
        return 'OK'

    try:
        # print("[CLEAN2] Cleaning HI emission using SoFiA mask for Sources {}.".format(args.sources))
        clean.map = 'map_{:02}_{:02}_{:04}'.format(b, minc, chan[i])
        clean.beam = 'beam_{:02}_{:02}_{:04}'.format(b, minc, chan[i])
        clean.out = 'model_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])
        clean.cutoff = lineimagestats[2] * 0.5
        clean.niters = 10000
        clean.region = '"mask(mask_{:02}_{:02}_{:04}/)"'.format(b, minc, chan[i])
        clean.go()
    except RuntimeError:
        print("\tProblems CLEANing data for channel {}.".format(chan[i]))
        return

    try:
        # print("[CLEAN2] Restoring line cube.")
        restor.model = 'model_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])
        restor.beam = 'beam_{:02}_{:02}_{:04}'.format(b, minc, chan[i])
        restor.map = 'map_{:02}_{:02}_{:04}'.format(b, minc, chan[i])
        restor.out = 'image_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])
        restor.mode = 'clean'
        restor.go()
    except RuntimeError:
        print("\tProblems RESTORing data for channel {}.".format(chan[i]))

    try:
        # print("[CLEAN2] Making residual cube.")
        out_array = ['image_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])]
        if args.all:
            restor.mode = 'residual'  # Create the residual image
            restor.out = 'residual_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])
            restor.go()
            out_array = ['model_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i]),
                         'image_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i]),
                         'residual_{:02}_{:02}_{:04}'.format(b, minc + 1, chan[i])]
    except RuntimeError:
        print("\tProblems RESTORing data with RESIDUALS for channel {}.".format(chan[i]))

    try:
        for name in out_array:
            fits.op = 'xyout'
            fits.in_ = name
            fits.out = name + '.fits'
            fits.go()
    except RuntimeError:
        print("Problems writing data to FITS for channel {}. Check miriad files instead.".format(chan[i]))

    try:
        new_cleancube_data[chan[i], :, :] = pyfits.getdata('image_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
        if args.all:
            new_modelcube_data[chan[i], :, :] = pyfits.getdata('model_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
            new_residualcube_data[chan[i], :, :] = pyfits.getdata('residual_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
    except RuntimeError:
        print("Couldn't add some channel to new cube because of above errors. Continue to next.".format(chan[i]))

    return 'OK'


###################################################################

parser = ArgumentParser(description="Do source finding in the HI spectral line cubes for a given taskid,"
                                    " beam, cubes", formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='190915041',
                    help='Specify the input taskid or field (default: %(default)s).')

parser.add_argument('-b', '--beams', default='0-39',
                    help='Specify a range (0-39) or list (3,5,7,11) of beams on which to do source finding'
                         ' (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3',
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-s', '--sources', default='all',
                    help='Specify sources to clean.  Can specify range or list. DEPRECATED? (default: %(default)s).')

parser.add_argument('-n', "--nospline",
                    help="Don't do spline fitting; so source finding on only continuum filtered cube.",
                    action='store_true')

parser.add_argument('-o', "--overwrite",
                    help="If option is included, overwrite old clean, model, and residual FITS files.",
                    action='store_true')

parser.add_argument('-j', "--njobs",
                    help="Number of jobs to run in parallel (default: %(default)s) tested on happili-05.",
                    default=18)

parser.add_argument('-a', "--all",
                    help="Write residual and model cubes as well as the cleaned cubes.",
                    action='store_true')

# Parse the arguments above
args = parser.parse_args()
njobs = int(args.njobs)

###################################################################

# Range of cubes/beams to work on:
taskid = args.taskid
cubes = [int(c) for c in args.cubes.split(',')]

if '-' in args.beams:
    b_range = args.beams.split('-')
    beams = np.array(range(int(b_range[1])-int(b_range[0])+1)) + int(b_range[0])
else:
    beams = [int(b) for b in args.beams.split(',')]

overwrite = args.overwrite

cube_name = 'HI_image_cube'
beam_name = 'HI_beam_cube'
alta_dir = '/altaZone/archive/apertif_main/visibilities_default/'

header = ['name', 'id', 'x', 'y', 'z', 'x_min', 'x_max', 'y_min', 'y_max', 'z_min', 'z_max', 'n_pix',
          'f_min', 'f_max', 'f_sum', 'rel', 'flag', 'rms', 'w20', 'w50', 'ell_maj', 'ell_min', 'ell_pa',
          'ell3s_maj', 'ell3s_min', 'ell3s_pa', 'kin_pa', "err_x", "err_y", "err_z", "err_f_sum", 'taskid', 'beam',
          'cube']

catParNames = ("name", "id", "x", "y", "z", "x_min", "x_max", "y_min", "y_max", "z_min", "z_max", "n_pix",
               "f_min", "f_max", "f_sum", "rel", "flag", "rms", "w20", "w50", "ell_maj", "ell_min", "ell_pa",
               "ell3s_maj", "ell3s_min", "ell3s_pa", "kin_pa", "err_x", "err_y", "err_z", "err_f_sum", "taskid", "beam",
               "cube")
catParUnits = ("-", "-", "pix", "pix", "chan", "pix", "pix", "pix", "pix", "chan", "chan", "-",
               "Jy/beam", "Jy/beam", "Jy/beam", "-", "-", "Jy/beam", "chan", "chan", "pix", "pix", "pix",
               "pix", "pix", "deg", "deg", "pix", "pix", "pix", "Jy/beam", "-", "-", "-")
catParFormt = ("%12s", "%7i", "%10.3f", "%10.3f", "%10.3f", "%7i", "%7i", "%7i", "%7i", "%7i", "%7i", "%8i",
               "%10.7f", "%10.7f", "%12.6f", "%8.6f", "%7i", "%12.6f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f",
               "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%10.3f", "%12.6f", "%10s", "%7i", "%7i")  #changed taskid from %10i to %10s

prepare = apercal.prepare()
managefiles.director(prepare, 'ch', taskid)

for b in beams:
    # loc = taskid + '/B0' + str(b).zfill(2) + '/'
    # loc = taskid + '/'
    mask_loc = 'mos_' + taskid + '/'
    # print("\t{}".format(loc))

    # managefiles.director(prepare, 'ch', loc)

    for c in cubes:
        cube_name = 'HI_B0' + str(b).zfill(2) + '_cube' + str(c) + '_image'
        beam_name = 'HI_B0' + str(b).zfill(2) + '_cube' + str(c) + '_psf'

        line_cube = cube_name + '.fits'
        beam_cube = beam_name + '_full.fits'   # Update to the expanded beam (deleted later to save space)
        maskfits = '/mnt/data/' + mask_loc + taskid + '_HIcube2_image_sofiaFS_mask_bin' + str(b).zfill(2) + '_regrid.fits'
        splinefits = cube_name[:-6] + '_spline.fits'

        if os.path.isfile(maskfits):

            mask_expr = '"(<mask_' + str(b).zfill(2) + '_sofia>.eq.-1).or.(<mask_' + str(b).zfill(2) + '_sofia>.ge.0.01)"'

            if (not args.nospline) & ((not os.path.isfile(splinefits)) | args.overwrite):

                print("[CLEAN2] Splinefit while avoiding sources for Beam {:02}, Cube {}.".format(b, c))
                os.system('cp {} {}'.format(line_cube, splinefits))
                print("\t{}".format(splinefits))
                new_splinecube = pyfits.open(splinefits, mode='update')
                maskcube = pyfits.open(maskfits)
                orig = pyfits.open(line_cube)
                orig_data = orig[0].data
                new_splinecube_data = new_splinecube[0].data
                new_splinecube_data[maskcube[0].data > 0] = 0.0
                maskcube.close()

                ################################################
                # Parallelization of spline
                xx, yy = range(new_splinecube_data.shape[1]), range(new_splinecube_data.shape[2])
                x, y = np.meshgrid(xx, yy)
                x, y = x.ravel(), y.ravel()
                ncases = len(x)
                print(" - " + str(ncases) + " cases found")

                if njobs > 1:
                    print(" - Running in parallel mode (" + str(njobs) + " jobs simultaneously)")
                elif njobs == 1:
                    print(" - Running in serial mode")
                else:
                    print("[ERROR] invalid number of NJOBS. Please use a positive number.")
                    exit()

                # Managing the work PARALLEL or SERIAL accordingly
                if njobs > cpu_count():
                    print(
                        "  [WARNING] The chosen number of NJOBS seems to be larger than the number of CPUs in the system!")

                # Create Queues
                print("    - Creating Queues")
                inQueue = Queue()
                outQueue = Queue()

                # Create worker processes
                print("    - Creating worker processes")
                ps = [Process(target=worker, args=(inQueue, outQueue)) for _ in range(njobs)]

                # Start worker processes
                print("    - Starting worker processes")
                for p in ps: p.start()

                # Fill the queue
                print("    - Filling up the queue")
                for i in trange(ncases):
                    inQueue.put((i))

                # Now running the processes
                print("    - Running the processes")
                output = [outQueue.get() for _ in trange(ncases)]

                # Send stop signal to stop iteration
                for _ in range(njobs): inQueue.put('STOP')

                # Stop processes
                print("    - Stopping processes")
                for p in ps: p.join()

                # Updating the Splinecube file with the new data
                print(" - Updating the Splinecube file")
                new_splinecube.data = new_splinecube_data
                new_splinecube.flush()

                # Closing files
                print(" - Closing files")
                new_splinecube.close()
                orig.close()

                ################################################

            if args.nospline:
                f = pyfits.open(filteredfits)
                print("[CLEAN2] Determining the statistics from the filtered Beam {:02}, Cube {}.".format(b, c))
            elif os.path.isfile(splinefits):
                f = pyfits.open(splinefits)
                print("[CLEAN2] Determining the statistics from the spline fitted Beam {:02},"
                      " Cube {}.".format(b, c))
            else:
                print("\tNo spline fitted cube found.  Exiting!")
                exit()
            nchan = f[0].data.shape[0]
            mask = np.ones(nchan, dtype=bool)
            if c == 3:
                mask[376:662] = False
            lineimagestats = [np.nanmin(f[0].data[mask]), np.nanmax(f[0].data[mask]), np.nanstd(f[0].data[mask])]
            f.close()
            print("\tImage min, max, std: {}".format(lineimagestats[:]))

            # Output what exactly is being used to clean the data
            print("\t{}".format(maskfits))
            m = pyfits.open(maskfits)
            # Get channels to clean:
            chan_mask = np.sum(m[0].data, axis=(1, 2))
            chan = np.where(chan_mask > 0)[0]
            m.close()

            # Delete any pre-existing Miriad files.
            os.system('rm -rf model_'+str(b).zfill(2)+'* beam_'+str(b).zfill(2)+'* map_'+str(b).zfill(2) +
                      '* image_'+str(b).zfill(2)+'* mask_'+str(b).zfill(2)+'* residual_'+str(b).zfill(2)+'*')

            print("[CLEAN2] Reading in FITS files, making Miriad mask.")

            fits = lib.miriad('fits')
            fits.op = 'xyin'
            if args.nospline:
                fits.in_ = line_cube
            else:
                fits.in_ = splinefits
            fits.out = 'map_' + str(b).zfill(2) + '_00'
            fits.go()

            if not os.path.isfile(beam_cube):
                print("[CLEAN2] Expanding synthesized beam.")
                trim_beam.main(beam_name+'.fits', beam_cube, 1)
            fits.in_ = beam_cube
            fits.out = 'beam_'+str(b).zfill(2)+'_00'
            fits.go()

            # Work with mask_sofia in current directory...otherwise worried about miriad character length for mask_expr
            fits.in_ = maskfits
            fits.out = 'mask_'+str(b).zfill(2)+'_sofia'
            fits.go()

            maths = lib.miriad('maths')
            maths.out = 'mask_'+str(b).zfill(2)+'_00'
            maths.exp = '"<mask_'+str(b).zfill(2)+'_sofia>"'
            maths.mask = mask_expr
            maths.go()

            if args.all:
                print("[CLEAN2] Initialize clean, model, and residual cubes")
            else:
                print("[CLEAN2] Initialize clean cube")
            if args.nospline:
                dirty_cube = line_cube
                outcube = line_cube[:-5]
            else:
                dirty_cube = splinefits
                outcube = splinefits[:-5]  # + '_rep'

            new_cleanfits = outcube + '_clean_image.fits'
            os.system('cp {} {}'.format(dirty_cube, new_cleanfits))
            print("\t{}".format(new_cleanfits))
            new_cleancube = pyfits.open(new_cleanfits, mode='update')
            new_cleancube_data = new_cleancube[0].data

            if args.all:
                new_modelfits = outcube + '_model.fits'
                os.system('cp {} {}'.format(dirty_cube, new_modelfits))
                pre_model = pyfits.open(new_modelfits, mode='update')
                pre_model[0].data = pre_model[0].data * np.nan
                pre_model.flush()
                pre_model.close()
                print("\t{}".format(new_modelfits))
                new_modelcube = pyfits.open(new_modelfits, mode='update')
                new_modelcube_data = new_modelcube[0].data

                new_residualfits = outcube + '_residual.fits'
                os.system('cp {} {}'.format(dirty_cube, new_residualfits))
                print("\t{}".format(new_residualfits))
                new_residualcube = pyfits.open(new_residualfits, mode='update')
                new_residualcube_data = new_residualcube[0].data

            # Initialize arrays based on cleaned channels
            bmaj_arr = np.zeros(len(chan))
            bmin_arr = np.zeros(len(chan))
            bpa_arr = np.zeros(len(chan))

            nminiter = 1
            clean = lib.miriad('clean')
            restor = lib.miriad('restor')  # Create the restored image
            for minc in range(nminiter):
                print("[CLEAN2] Clean & restor HI emission using SoFiA mask for Sources {}.".format(args.sources))

                ################################################
                # Parallelization of cleaning/restoring
                ncases = len(chan)
                print(" - " + str(ncases) + " cases found")

                if njobs > 1:
                    print(" - Running in parallel mode (" + str(njobs) + " jobs simultaneously)")
                elif njobs == 1:
                    print(" - Running in serial mode")
                else:
                    print("[ERROR] invalid number of NJOBS. Please use a positive number.")
                    exit()

                # Managing the work PARALLEL or SERIAL accordingly
                if njobs > cpu_count():
                    print(
                        "\t[WARNING] The chosen number of NJOBS seems to be larger than the number of CPUs"
                        " in the system!")

                # Create Queues
                print("    - Creating Queues")
                inQueue = Queue()
                outQueue = Queue()

                # Create worker2 processes
                print("    - Creating worker2 processes")
                ps = [Process(target=worker2, args=(inQueue, outQueue)) for _ in range(njobs)]

                # Start worker2 processes
                print("    - Starting worker2 processes")
                for p in ps:
                    p.start()

                # Fill the queue
                print("    - Filling up the queue")
                for i in trange(ncases):
                    inQueue.put(i)

                # Now running the processes
                print("    - Running the processes")
                output = [outQueue.get() for _ in trange(ncases)]

                # Send stop signal to stop iteration
                for _ in range(njobs): inQueue.put('STOP')

                # Stop processes
                print("    - Stopping processes")
                for p in ps: p.join()

                ################################################

            print(" - Updating the clean file")
            new_cleancube.data = new_cleancube_data

            if args.all:
                print(" - Updating the model file")
                new_modelcube.data = new_modelcube_data
                print(" - Updating the residual file")
                new_residualcube.data = new_residualcube_data

                print("[CLEAN2] Updating history of reassembled model, residual, clean cubes")
                for i in range(len(chan)):
                    if os.path.isfile('residual_{:02}_{:04}.fits'.format(minc + 1, chan[i])):
                        residual_chan_hdr = pyfits.getheader('residual_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
                        model_chan_hdr = pyfits.getheader('model_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
                        for hist in residual_chan_hdr[-35:]['HISTORY']:
                            new_residualcube[0].header['HISTORY'] = hist
                        for hist in model_chan_hdr[-26:]['HISTORY']:
                            new_modelcube[0].header['HISTORY'] = hist
                new_residualcube[0].header['HISTORY'] = \
                    'Individual images reassembled using aper_sf2/clean2.py by KMHess'
                new_modelcube[0].header['HISTORY'] = \
                    'Individual images reassembled using aper_sf2/clean2.py by KMHess'

                print(" - Saving the new model file")
                new_modelcube.flush()
                new_modelcube.close()

                print(" - Saving the new residual file")
                new_residualcube.flush()
                new_residualcube.close()

            else:
                print("[CLEAN2] Updating history of reassembled clean cube")

            for i in range(len(chan)):
                if os.path.isfile('image_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i])):
                    clean_chan_hdr = pyfits.getheader('image_{:02}_{:02}_{:04}.fits'.format(b, minc + 1, chan[i]))
                    for hist in clean_chan_hdr[-35:]['HISTORY']:  # Determined through trial and error
                        new_cleancube[0].header['HISTORY'] = hist
                    bmaj_arr[i] = clean_chan_hdr['BMAJ']
                    bmin_arr[i] = clean_chan_hdr['BMIN']
                    bpa_arr[i] = clean_chan_hdr['BPA']
            new_cleancube[0].header['HISTORY'] = 'Individual images reassembled using sourcefinding/clean2.py by KMHess'

            print("[CLEAN2] Adding median beam properties to primary header")
            med_bmaj, med_bmin, med_bpa = np.median(bmaj_arr), np.median(bmin_arr), np.median(bpa_arr)
            new_cleancube[0].header.set('BMAJ', med_bmaj, 'median clean beam bmaj')
            new_cleancube[0].header.set('BMIN', med_bmin, 'median clean beam bmin')
            new_cleancube[0].header.set('BPA', med_bpa, 'median clean beam pa')

            print("[CLEAN2] Adding channel clean beam properties to BEAMS extension table")
            col1 = pyfits.Column(name='BMAJ', format='1E', unit='deg', array=bmaj_arr)
            col2 = pyfits.Column(name='BMIN', format='1E', unit='deg', array=bmin_arr)
            col3 = pyfits.Column(name='BPA', format='1E', unit='deg', array=bpa_arr)
            col4 = pyfits.Column(name='CHAN', format='1J', array=chan)
            beam_hdu = pyfits.BinTableHDU.from_columns([col1, col2, col3, col4])
            beam_hdu.name = 'BEAMS'
            beam_hdu.header.comments['NAXIS2'] = 'number of channels'
            new_cleancube.append(beam_hdu)

            print(" - Saving the new clean file")
            new_cleancube.flush()
            new_cleancube.close()

            # Don't need to output sofia only mask, because same info is in regridded mask.
            # print("[CLEAN2] Writing mask with only cleaned sources")
            # os.system('rm -rf {}_clean_mask.fits'.format(outcube))
            # fits.op = 'xyout'
            # fits.in_ = 'mask_' + str(minc).zfill(2)
            # if not args.nospline:
            #     fits.out = outcube + '_clean_mask.fits'
            #     fits.go()

            # Clean up extra Miriad files
            os.system('rm -rf model_'+str(b).zfill(2)+'* beam_'+str(b).zfill(2)+'* map_'+str(b).zfill(2) +
                      '* image_'+str(b).zfill(2)+'* mask_'+str(b).zfill(2)+'* residual_'+str(b).zfill(2)+'*')
            os.system('rm -rf ' + beam_cube)
        else:
            print("no mask?")

print("[CLEAN2] Done.")
