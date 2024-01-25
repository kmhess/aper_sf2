# https://stackoverflow.com/questions/66466082/how-do-i-write-a-snakemake-input-when-not-all-jobs-successfully-output-files-fro

import os

configfile: "config.yaml"
FIELD = config['FIELD']
CUBE = str(config['CUBE'])

wildcard_constraints:
    beam="\d+"

rule all:
    input:
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_filtered.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.txt",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.xml",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask-2d.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_rel.eps",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_skellam.eps",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_filtered_spline.fits"

checkpoint stack_obs:
    output:
        directory(FIELD)
    shell:
        "python3 /mnt/scratch/stuff/aper_cube_stack/cube_stack.py -f {output} -b 0-39 -c "+CUBE

#run a separate job from each output of rule a
rule generate_pb:
    input:
        FIELD+"/HI_B0{bm}_cube"+CUBE+"_image.fits",
        "/mnt/scratch/apertif/cbeams/{bm}_gp_avg_orig.fits"
    output:
        FIELD+"/HI_B0{bm}_cube"+CUBE+"_pb.fits"
    shell:
        "python3 /mnt/scratch/stuff/aper_cube_stack/modules/regrid_aperpb.py -t {input[0]} -b {wildcards.bm} -p /mnt/scratch/apertif/cbeams/"

# input function for the rule aggregate
def aggregate_input(wildcards):
    # say this rule depends on a checkpoint.  DAG evaulation pauses here
    checkpoint_output = checkpoints.stack_obs.get(**wildcards).output[0]
    found_files = expand(FIELD+"/HI_B0{beam}_cube"+CUBE+"_pb.fits", 
                         beam=glob_wildcards(os.path.join(checkpoint_output, "HI_B0{beam}_cube"+CUBE+"_image.fits")).beam)
    found_files2 = expand(FIELD+"/HI_B0{beam}_cube"+CUBE+"_image.fits", 
                         beam=glob_wildcards(os.path.join(checkpoint_output, "HI_B0{beam}_cube"+CUBE+"_image.fits")).beam)
    return found_files+found_files2

rule make_mosaic:
    input:
        pb = aggregate_input,
    output:
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_noise.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_weights.fits"
    run:
        input = list(input)
        mos_params = " ".join([i.split("/")[1] for i in input if "image" in i])
        os.system('mosaic-queen -mc 0.1 -n '+FIELD+'_HIcube'+CUBE+' -i '+FIELD+' -o mos_'+FIELD+' -t '+mos_params+' -r')
        os.system('rm -rf mos_'+FIELD+'/*imageR*')
        os.system('rm -rf mos_'+FIELD+'/*pbR*')

rule run_sofia:
    input:
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_noise.fits"
    output:
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_filtered.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.txt",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.xml",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask-2d.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask.fits",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_rel.eps",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_skellam.eps",
        "mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_filtered_spline.fits"
    shell:
        "python3 /mnt/scratch/stuff/aper_sf2/sourcefinding.py -t "+FIELD+" -c "+CUBE+" -m"

