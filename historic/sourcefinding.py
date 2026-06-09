import os
from argparse import ArgumentParser, RawTextHelpFormatter

from modules.natural_cubic_spline import fspline
from src import checkmasks
from src.get_gal_em_range import find_rms_range

from astropy.io import fits
import dask.array as da
import numpy as np
from multiprocessing import Queue, Process, cpu_count
from tqdm.auto import trange
import time as testtime


dir_name = os.path.dirname(__file__)


def make_param_file(loc_dir=None, cube_name=None, cube=None, mosaic=False, gal_em_range=None):
    param_template = dir_name + '/sofia_parameter_template.par'
    new_paramfile = loc_dir + 'sofia_parameter.par'
    if mosaic == True:
        param_template = dir_name + '/mosaic_parameter_template.par'
        new_paramfile = loc_dir + 'mosaic{}_parameter.par'.format(cube)
    outlog = loc_dir + 'sourcefinding{}.out'.format(cube)

    # Edit parameter file (remove lines that need editing)
    os.system('grep -vwE "(input.data)" ' + param_template + ' > ' + new_paramfile)
    # if mosaic:                          # FILTERED SPLINE ALREADY HAS NOISE CORRECTION APPLIED!
    #     os.system('grep -vwE "(input.noise)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)  #
    os.system('grep -vwE "(output.filename)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)
    if cube == 3:
        os.system('grep -vwE "(linker.radiusZ)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)
        os.system('grep -vwE "(linker.maxSizeXY)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)
        os.system('grep -vwE "(linker.maxSizeZ)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)

    # Add back the parameters needed
    if not args.nospline:
        os.system('echo "input.data                 =  ' + splinefits + '" >> ' + new_paramfile)
        # outroot = cube_name + '_sofia'          #
        outroot = cube_name + '_sofiaFS'
    else:
        os.system('echo "input.data                 =  ' + filteredfits + '" >> ' + new_paramfile)
        outroot = cube_name + '_sofiaB'

    os.system('echo "output.filename            =  ' + outroot + '" >> ' + new_paramfile)
    # if mosaic:                          # FILTERED SPLINE ALREADY HAS NOISE CORRECTION APPLIED!
    #     os.system('echo "input.noise                =  ' + noisefits + '" >> ' + new_paramfile)  #
    if cube == 3:
        # Figures out which channels to exclude
        gal_em_range_str = ',{min_chan},{max_chan}'.format(min_chan=str(gal_em_range[0]), max_chan=str(gal_em_range[1]))
        if (gal_em_range[0]!=0) and (gal_em_range[1]!=0):
            print("\tExcluding Galactic emission between channels {}".format(gal_em_range_str))
            os.system('grep -vwE "(flag.region)" ' + new_paramfile + ' > temp && mv temp ' + new_paramfile)
            os.system('echo "flag.region                =  0,2600,0,2000{}" >> '.format(gal_em_range_str) + new_paramfile)
        os.system('echo "linker.radiusZ             =  6" >> ' + new_paramfile)
        os.system('echo "linker.maxSizeXY           =  280" >> ' + new_paramfile)
        os.system('echo "linker.maxSizeZ            =  385" >> ' + new_paramfile)

    return new_paramfile, outlog


def worker(inQueue, outQueue):

    """
    Defines the worker process of the parallelisation with multiprocessing.Queue
    and multiprocessing.Process.
    """
    for i in iter(inQueue.get, 'STOP'):

        status = run(i)

        outQueue.put((status))


def run(i):
    global splinecube_data

    try:
        # Do the spline fitting on the z-axis to masked cube
        fit = fspline(np.linspace(1, orig_data.shape[0], orig_data.shape[0]),
                      np.nan_to_num(splinecube_data[:, x[i], y[i]]), k=5)
        splinecube_data[:, x[i], y[i]] = orig_data[:, x[i], y[i]] - fit
        return 'OK'

    except Exception:
        print("[ERROR] Something went wrong with the Spline fitting [" + str(i) + "]")
        return np.nan


###################################################################

parser = ArgumentParser(description="Do source finding in the HI spectral line cubes for a given taskid, beam, cubes",
                        formatter_class=RawTextHelpFormatter)

parser.add_argument('-t', '--taskid', default='190915041', required=True,
                    help='Specify the input taskid or field in case of mosaic (default: %(default)s).')

parser.add_argument('-b', '--beams', default='0-39', required=False,
                    help='Specify a range (0-39) or list (3,5,7,11) of beams on which to do source finding'
                         ' (default: %(default)s).')

parser.add_argument('-c', '--cubes', default='1,2,3', required=True,
                    help='Specify the cubes on which to do source finding (default: %(default)s).')

parser.add_argument('-r', "--chanrange", default=None, required=False,
                    help='Specify a channel range (e.g., 480-500) to exclude during source finding. Specify 0-0 if no'
                         ' channels should be excluded.'
                         ' (default: calculated from Galactic HI exclusion algorithm).')

parser.add_argument('-o', "--overwrite", required=False,
                    help="If option is included, overwrite old continuum filtered and/or spline fitted file if either"
                         " exists.",
                    action='store_true')

parser.add_argument('-n', "--nospline", required=False,
                    help="Don't do spline fitting; so source finding on only continuum filtered cube. UNTESTED HERE.",
                    action='store_true')

parser.add_argument('-j', "--njobs", required=False,
                    help="Number of jobs to run in parallel (default: %(default)s) tested on happili-05.",
                    default=18)

parser.add_argument('-m', "--mosaic", required=False,
                    help="If option is included, operate on a mosaic field.",
                    action='store_true')

parser.add_argument('-d', '--directory', default='', required=False,
                    help='Specify the directory where taskid/field folders live containing the data (default: %(default)s).')

###################################################################

# Parse the arguments above
args = parser.parse_args()
njobs = int(args.njobs)

# Range of cubes/beams to work on:
taskid = args.taskid
cubes = [int(c) for c in args.cubes.split(',')]
if '-' in args.beams:
    b_range = args.beams.split('-')
    beams = np.array(range(int(b_range[1]) - int(b_range[0]) + 1)) + int(b_range[0])
else:
    beams = [int(b) for b in args.beams.split(',')]
overwrite = args.overwrite
d = args.directory

# If operating on a mosaic, give a dummy beam.
if args.mosaic:
    beams = [40]

# Main source finding code for all cubes/beams
for b in beams:
    # Define some file names and work space:
    # loc = taskid + '/B0' + str(b).zfill(2) + '/'
    loc = d + '/' + taskid + '/'
    # Snakemake required input! (for now)
    if args.mosaic:                                             # Enable while using snakemake
        loc = d + '/mos_' + taskid + '/'                               # Enable while using snakemake
    for c in cubes:
        # cube_name = 'HI_image_cube' + str(c)
        cube_name = 'HI_B0' + str(b).zfill(2) + '_cube' + str(c) + '_image'
        noisefits = None
        if args.mosaic:
            cube_name = taskid + '_HIcube' + str(c) + '_image'
            noisefits = loc + taskid + '_HIcube' + str(c) + '_noise.fits'
        if b == 40:
            print("[SOURCEFINDING] Working on full mosaic field {} Cube {}".format(taskid, c))
        else:
            print("[SOURCEFINDING] Working on Beam {:02} Cube {}".format(b, c))

        sourcefits = loc + cube_name + '.fits'
        filteredfits = loc + cube_name + '_filtered.fits'
        splinefits = loc + cube_name + '_filtered_spline.fits'              #DELETE THIS WHEN NECESSARY *******************************
        # splinefits = loc + cube_name + '_spline.fits'                     #
        # Output exactly where sourcefinding is starting
        print('\t' + sourcefits)
        new_filter_par = 'filtering{}.par'.format(c)

        # Check to see if the continuum filtered file exists.  If not, make it  with SoFiA-2
        if (not overwrite) & (os.path.isfile(filteredfits) | os.path.isfile(splinefits)):
            print("[SOURCEFINDING] Continuum filtered file exists and will not be overwritten.")
        elif os.path.isfile(sourcefits):
            print("[SOURCEFINDING] Making continuum filtered file.")
            os.system('grep -vwE "(input.data)" ' + dir_name + '/filtering.par > ' + loc + new_filter_par)
            if noisefits:
                os.system('grep -vwE "(input.noise)" ' + loc + new_filter_par + ' > temp && mv temp ' +
                          loc + new_filter_par)
            os.system('echo "input.data                 =  ' + sourcefits + '" >> ' + loc + new_filter_par)
            if noisefits:
                os.system('echo "input.noise                =  ' + noisefits + '" >> ' + loc + new_filter_par)
            os.system('sofia ' + loc + new_filter_par + ' >> test.log')
        else:
            if b == 40:
                print("\tFull mosaic for field {} Cube {} is not present in this directory.".format(taskid, c))
            else:
                print("\tBeam {:02} Cube {} is not present in this directory.".format(b, c))
            continue

        # Check to see if the spline fitted file exists.  If not, make it from filtered file.
        if (not overwrite) & os.path.isfile(splinefits):
            print("[SOURCEFINDING] Spline fitted file exists and will not be overwritten.")
        elif os.path.isfile(filteredfits) & (not args.nospline):
            print("[SOURCEFINDING] Spline fitting this file: {}.".format(filteredfits))
            print(" - Loading the input cube")
            os.system('cp {} {}'.format(filteredfits, splinefits))              #DELETE THIS WHEN NECESSARY *******************************
            # os.system('cp {} {}'.format(sourcefits, splinefits))              #
            splinecube = fits.open(splinefits, mode='update')
            orig = fits.open(filteredfits)                                 #DELETE THIS WHEN NECESSARY *******************************
            # orig = fits.open(sourcefits)                                 #
            orig_data = orig[0].data
            splinecube_data = splinecube[0].data

            # Try masking strong sources to not bias fit
            print(" - Masking strong sources to not bias fit")
            mask = 2.5 * da.nanstd(orig_data)
            splinecube_data[np.abs(splinecube_data) >= mask] = np.nan

            # Defining the cases to analyse
            print(" - Defining the cases to analyse")
            xx = range(orig_data.shape[1])
            yy = range(orig_data.shape[2])
            x, y = np.meshgrid(xx, yy)
            x = x.ravel()
            y = y.ravel()
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
            splinecube.data = splinecube_data
            splinecube.flush()

            # Closing files
            print(" - Closing files")
            orig.close()
            splinecube.close()

        elif os.path.isfile(sourcefits) & args.nospline:
            print("\tWill not perform spline fitting.  Do source finding on just continuum filtered file.")
            print("\t[WARNING]: this is not the default but the file names are the SAME!"
                  " Keep track of what you're doing for future steps !!!")
        else:
            if b == 40:
                print("\tFull mosaic for field {} Cube {} is not present in this directory.".format(taskid, c))
            else:
                print("\tBeam {:02} Cube {} is not present in this directory.".format(b, c))
            continue

        # Reads in range of Galactic HI emission channels to exclude, calculates if unspecified:
        if args.chanrange != None:
            str_chan_range = args.chanrange.split('-')
            chan_range = int(str_chan_range[0]), int(str_chan_range[1])
        else:
            chan_range = find_rms_range(loc_dir=loc, field_name=taskid, splinefits=splinefits)
        
        new_paramfile, outlog = make_param_file(loc_dir=loc, cube_name=cube_name, cube=c, mosaic=args.mosaic, gal_em_range=chan_range)
        try:
            print("[SOURCEFINDING] Cleaning up old mask and catalog files before source finding.")
            # os.system('rm -rf ' + loc + cube_name + '*_sofia_*.fits ' + loc + cube_name + '*_sofia_cat.txt')
        except:
            pass
        print("[SOURCEFINDING] Doing source finding with SoFiA parameter file {}.".format(new_paramfile))
        tic1 = testtime.perf_counter()
        os.system('sofia ' + new_paramfile + ' >> ' + outlog)
        toc1 = testtime.perf_counter()
        print(f"Do sofia: {toc1 - tic1:0.4f} seconds")

    # After all cubes are done, run checkmasks to get summary plots for cleaning:
    checkmasks.main(loc, taskid, [b], cubes=cubes, nospline=args.nospline, mosaic=args.mosaic)
