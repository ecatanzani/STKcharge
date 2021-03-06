from ROOT import TFile, TH1D, TH2D
import sys
import os


def compute_final_histos(condor_dir_list, opts):
    
    # Charge histos
    h_chargeX = TH1D()
    h_chargeY = TH1D()
    h_charge = TH1D()
    h_charge2D = TH2D()

    for dIdx, tmp_dir in enumerate(condor_dir_list):
        tmp_dir += "/outFiles"
        tmp_dir_list = os.listdir(tmp_dir)
        for elm in tmp_dir_list:
            if elm.startswith("STKchargeOutFile_"):
                rFile_path = tmp_dir + "/" + elm

        # Open ROOT output file
        rFile = TFile.Open(rFile_path, "READ")
        if rFile.IsOpen():
            if opts.verbose:
                if dIdx == 0:
                    print('\nReading file {}: {}'.format((dIdx+1), rFile_path))
                else:
                    print('Reading file {}: {}'.format((dIdx+1), rFile_path))
        else:
            print('Error reading file {}: {}'.format((dIdx+1), rFile_path))
            sys.exit()
        
        # Reading histos
        h_chargeX_tmp = rFile.Get("h_chargeX")
        h_chargeY_tmp = rFile.Get("h_chargeY")
        h_charge_tmp = rFile.Get("h_charge")
        h_charge2D_tmp = rFile.Get("h_charge2D")
        
        # Unlink histos
        h_chargeX_tmp.SetDirectory(0)
        h_chargeY_tmp.SetDirectory(0)
        h_charge_tmp.SetDirectory(0)
        h_charge2D_tmp.SetDirectory(0)

        # Clone output file
        rFile.Close()

        # Add histos
        if dIdx == 0:
            
            h_chargeX = h_chargeX_tmp.Clone("h_chargeX")
            h_chargeY = h_chargeY_tmp.Clone("h_chargeY")
            h_charge = h_charge_tmp.Clone("h_charge")
            h_charge2D = h_charge2D_tmp.Clone("h_charge2D")

        else:
            
            h_chargeX.Add(h_chargeX_tmp)
            h_chargeY.Add(h_chargeY_tmp)
            h_charge.Add(h_charge_tmp)
            h_charge2D.Add(h_charge2D_tmp)
    
    # Create output file for full histos
    fOut = TFile.Open(opts.output, "RECREATE")
    if fOut.IsOpen():
        if opts.verbose:
            print('Output TFile has been created: {}'.format(opts.output))
    else:
        print('Error creating output TFile: {}'.format(opts.output))
        sys.exit()

    # Writing final histos to file
    h_chargeX.Write()
    h_chargeY.Write()  
    h_charge.Write()
    h_charge2D.Write()

    # Closing output file
    fOut.Close()