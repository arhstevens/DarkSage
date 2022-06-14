# Run Dark Sage many times over a prescribed range of input parameters

import subprocess
import os
import numpy as np

ds_dir = '/Users/adam/DarkSage/'
base_param_file = '/Users/adam/DarkSage/input/genesis_calibration.par'
base_outdir = '/Users/adam/DarkSage_runs/Genesis/L75n324/ParameterSweep/'

fmove_grid = np.arange(0.1,0.91,0.1)
eta_grid = np.arange(0.05,0.41,0.05) # expand later

for fmove in fmove_grid:
    for eta in eta_grid:
    
        # create directory for next run of Dark Sage
        fmove_str = str(int(fmove*100))
        eta_str = str(int(eta*100))
        if len(eta_str)==1: eta_str = '0' + eta_str

        outdir = base_outdir + 'f0p'+ fmove_str + '_e0p' + eta_str +'/'
        if not os.path.exists(outdir): os.makedirs(outdir)
        if os.path.isfile(outdir+'model_z0.000_0'): continue

        # make new parameter file
        fmove_str = str(round(fmove,2))
        eta_str = str(round(eta,2))
        
        print("\n\n\n====================================================")
        print("====================================================")
        print("====================================================")
        print("Running Dark Sage for fmove="+fmove_str+" and eta="+eta_str)
        print("====================================================")
        print("====================================================")
        print("====================================================\n\n\n")
        
        print("outdir", outdir, "\nbase_param_file", base_param_file)
        
        outdir_w_bslash = outdir
        i = 0
        for j in range(len(outdir_w_bslash)):
            print(i)
            i = outdir_w_bslash.find('/', i)
            if i==-1: break
            outdir_w_bslash = outdir_w_bslash[:i] + "\\" + outdir_w_bslash[i:]
            i += 2
            
#        outdir_w_bslash += "\/"
        print(outdir_w_bslash)
        
        os.system("sed \'s/^OutputDir.*/OutputDir              "+outdir_w_bslash+"/\' "+base_param_file+" >"+outdir+"tempfile1.par")
        os.system("sed \'s/^RadiativeEfficiency.*/RadiativeEfficiency         "+eta_str+"/\' "+outdir+"tempfile1.par"+" >"+outdir+"tempfile2.par")
        os.system("sed \'s/^GasSinkRate.*/GasSinkRate                 "+fmove_str+"/\' "+outdir+"tempfile2.par"+" >"+outdir+"genesis_calibration.par")
        
        subprocess.call(["rm", outdir+"tempfile1.par"])
        subprocess.call(["rm", outdir+"tempfile2.par"])
        
        print(outdir+"genesis_calibration.par")
        subprocess.call(["/opt/local/bin/mpirun-openmpi-mp", "-np", "8", "./darksage", outdir+"genesis_calibration.par"])

