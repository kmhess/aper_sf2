'''
TNH, created 10/21/2025
Script to make and run and slurm script to run an input field through the Apertif data reduction pipeline (up until first source finding).
Assumes the data is staged and stored in the /project/apdw/Data/ directory.

Usage: python make_pipeline_slurm.py <FIELD_NAME> <CHAN_RANGE>
Example: python make_pipeline_slurm.py S1021+5815 None
    or:  python make_pipeline_slurm.py S1021+5815 480-490
'''

import os, sys
from datetime import datetime
import subprocess

date_string = datetime.today().strftime('%Y-%m-%d')

def write_slurm(field_name, chan_range):
    with open("/project/apdw/Software/aper_sf2/slurm_out/slurm_"+field_name+"_"+date_string+".sh", "w") as f:
        f.write("#!/bin/bash\n")
        f.write("#SBATCH --job-name=aper_"+field_name+"          # Assign an short name to your job\n")
        f.write("#SBATCH --nodes=1                     # Number of nodes you require\n")
        f.write("#SBATCH --ntasks=1             # Total # of tasks across all nodes\n")
        f.write("#SBATCH --cpus-per-task=50             # Cores per task (>1 if multithread tasks)\n")
        f.write("#SBATCH --time=8:00:00               # Total run time limit (HH:MM:SS)\n")
        f.write("#SBATCH --output=/project/apdw/Software/aper_sf2/slurm_out/aper_"+field_name+"_"+date_string+".out    # STDOUT output file\n")
        f.write("#SBATCH --error=/project/apdw/Software/aper_sf2/slurm_out/aper_"+field_name+"_"+date_string+".err     # STDERR output file (optional)\n")
        f.write("#SBATCH --export=ALL                  # Export you current env to the job env\n")
        f.write("#SBATCH --mail-type=END\n")
        f.write("#SBATCH --mail-user=tnh57@rutgers.edu\n\n")

        f.write("conda deactivate\n")
        f.write("conda activate snakemake\n")
        f.write("cd /project/apdw/Software/aper_sf2\n\n")
        
        f.write("snakemake -s Snakefile_slurm --cores 50 --config FIELD="+field_name+" CHAN_RANGE="+chan_range+"  --rerun-incomplete")

def run_slurm(field_name):
    subprocess.run(["sbatch", "/project/apdw/Software/aper_sf2/slurm_out/slurm_"+field_name+"_"+date_string+".sh"])

if __name__ == "__main__":
    field_name = str(sys.argv[1])
    chan_range = str(sys.argv[2])

    #making slurm script
    write_slurm(field_name, chan_range)
    print("Slurm script written!")

    #running slurm script
    run_slurm(field_name)
    print("Slurm job run!")
