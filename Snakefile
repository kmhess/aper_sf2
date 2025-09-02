# https://stackoverflow.com/questions/66466082/how-do-i-write-a-snakemake-input-when-not-all-jobs-successfully-output-files-fro

import os

configfile: "config.yaml"
FIELD = config['FIELD']
CUBE = str(config['CUBE'])
DATA = str(config['PATH_TO_DATA'])
SOFTWARE = str(config['PATH_TO_SOFTWARE'])
CBEAMS = str(config['PATH_TO_CBEAMS'])

wildcard_constraints:
    beam="\d+"

rule all:
    input:
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.txt",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.xml",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask-2d.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_rel.eps",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_skellam.eps"

checkpoint stack_obs:
    output:
        directory(DATA+"/"+FIELD)
    shell:
        "python3 "+SOFTWARE+"/aper_cube_stack/cube_stack.py -f "+FIELD+" -b 0-39 -c "+CUBE+" -d "+DATA

#run a separate job from each output of rule a
rule generate_pb:
    input:
        DATA+"/"+FIELD+"/HI_B0{bm}_cube"+CUBE+"_image.fits",
        CBEAMS+"/{bm}_gp_avg_orig.fits"
    output:
        DATA+"/"+FIELD+"/HI_B0{bm}_cube"+CUBE+"_pb.fits"
    shell:
        "python3 "+SOFTWARE+"/aper_cube_stack/modules/regrid_aperpb.py -t {input[0]} -b {wildcards.bm} -p "+CBEAMS+"/"

# input function for the rule aggregate
def aggregate_input(wildcards):
    # say this rule depends on a checkpoint.  DAG evaulation pauses here
    checkpoint_output = checkpoints.stack_obs.get(**wildcards).output[0]
    found_files = expand(DATA+"/"+FIELD+"/HI_B0{beam}_cube"+CUBE+"_pb.fits", 
                         beam=glob_wildcards(os.path.join(checkpoint_output, "HI_B0{beam}_cube"+CUBE+"_image.fits")).beam)
    found_files2 = expand(DATA+"/"+FIELD+"/HI_B0{beam}_cube"+CUBE+"_image.fits", 
                         beam=glob_wildcards(os.path.join(checkpoint_output, "HI_B0{beam}_cube"+CUBE+"_image.fits")).beam)
    return found_files+found_files2

rule make_mosaic:
    input:
        pb = aggregate_input,
    output:
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_noise.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_weights.fits"
    run:
        input = list(input)
        mos_params = " ".join([i.split("/")[-1] for i in input if "image" in i])
        os.system('mosaic-queen -mc 0.1 -n '+FIELD+'_HIcube'+CUBE+' -i '+FIELD+' -o mos_'+FIELD+' -t '+mos_params+' -r')
        os.system('rm -rf '+DATA+'/mos_'+FIELD+'/*imageR*')
        os.system('rm -rf '+DATA+'/mos_'+FIELD+'/*pbR*')
        os.system('rm -rf '+DATA+'/mos_'+FIELD+'/*image_fields*')
        os.system('rm -rf '+DATA+'/'+FIELD+'/HI_B0{bm}_cube'+CUBE+'_pb.fits')

rule run_sofia:
    input:
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_noise.fits"
    output:
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.txt",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_cat.xml",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask-2d.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_mask.fits",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_rel.eps",
        DATA+"/mos_"+FIELD+"/"+FIELD+"_HIcube"+CUBE+"_image_sofiaFS_skellam.eps"
    shell:
        "python3 "+SOFTWARE+"/aper_sf2/sourcefinding.py -t "+FIELD+" -c "+CUBE+" -d "+DATA+" -m"

