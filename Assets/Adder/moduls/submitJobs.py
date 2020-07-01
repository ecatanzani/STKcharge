import subprocess
import shutil
import os


def resubmit_condor_jobs(skipped_dirs, opts):
    for dir in skipped_dirs:
        os.chdir(dir)
        outCondor = [filename for filename in os.listdir('.') if filename.startswith("output")]
        
        # Clean the job dir
        if outCondor:
            for elm in outCondor:
                if os.path.isdir(elm):
                    shutil.rmtree(elm)
                if os.path.isfile(elm):
                    os.remove(elm)
        
        # Submit HTCondor job
        if opts.verbose:
            print('Resubmitting job from folder: {}'.format(dir))

        subprocess.run("condor_submit -name sn-01.cr.cnaf.infn.it -spool crawler.sub", shell=True, check=True)
