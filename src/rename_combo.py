import os
from argparse import ArgumentParser, RawTextHelpFormatter
from astropy.table import Table

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
    filename = taskid + '_HIcube' + str(c) + '_clean_image'
    catalog = Table.read(mos_loc + filename + '_cat.xml')

    if args.sources == 'all':
        sources = [str(s + 1) for s in range(len(catalog))]
    elif '-' in args.sources:
        src_range = args.sources.split('-')
        sources = [str(s + int(src_range[0])) for s in range(int(src_range[1]) - int(src_range[0]) + 1)]
    else:
        sources = [str(s) for s in args.sources.split(',')]

    for s in sources:
        cat = catalog[catalog['id'] == int(s)]
        infile = f"{mos_loc}{filename}_figures/{filename}_{s}_combo.png"
        outfile = f"{mos_loc}{filename}_figures/AHC{cat['name'][0].split(' ')[1]}_{taskid}_{s}_combo.png"
        os.system(f"mv {infile} {outfile}")
