import os
import shutil
from ROOT import TFile


def getListOfFiles(condor_wd):
    list_dir = []
    skipped_dirs = []

    skipped_file_notFinalDir = 0
    skipped_file_notROOTfile = 0
    skipped_file_notReadable = 0
    skipped_file_noKeys = 0

    # Starting loop on output condor dirs
    for tmp_dir in os.listdir(condor_wd):
        if tmp_dir.startswith('job_'):
            full_dir_path = condor_wd + "/" + tmp_dir
            expected_condor_outDir = full_dir_path + "/outFiles"
            
            # Check if 'outFiles' dir exists
            if os.path.isdir(expected_condor_outDir):
                _list_dir = os.listdir(expected_condor_outDir)
                tmp_acc_full_path = ""
                for file in _list_dir:
                    if file.endswith(".root"):
                        tmp_acc_full_path = expected_condor_outDir + "/" + file
                        break
                
                # Check if output ROOT file exists
                if os.path.isfile(tmp_acc_full_path):
                    tmp_acc_file = TFile.Open(tmp_acc_full_path, "READ")
                    
                    # Check if output ROOT file is redable
                    if tmp_acc_file.IsOpen():
                        
                        # Check if output ROOT file has keys
                        outKeys = tmp_acc_file.GetNkeys()
                        
                        if outKeys:
                            list_dir.append(full_dir_path)

                        else:
                            # output ROOT file has been open but has not keys
                            skipped_dirs.append(full_dir_path)
                            skipped_file_noKeys += 1

                    else:
                        # output ROOT file has not been opened correctly
                        skipped_dirs.append(full_dir_path)
                        skipped_file_notReadable += 1

                else:
                    # output ROOT file does not exist
                    skipped_dirs.append(full_dir_path)
                    skipped_file_notROOTfile += 1
            
            else:
                # 'outFiles' dir does not exists
                skipped_dirs.append(full_dir_path)
                skipped_file_notFinalDir += 1

    return list_dir, skipped_dirs, skipped_file_notFinalDir, skipped_file_notROOTfile, skipped_file_notReadable, skipped_file_noKeys