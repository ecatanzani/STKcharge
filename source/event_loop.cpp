#include "data_loop.h"
#include "aggregate_events.h"
#include "data_cuts.h"
#include "charge.h"
#include "DAMPE_geo_structure.h"

#include "DmpFilterOrbit.h"
#include "DmpEvtHeader.h"
#include "DmpIOSvc.h"
#include "DmpCore.h"

inline void updateProcessStatus(const int evIdx, int &kStep, const int nevents)
{
    auto percentage = ((evIdx + 1) / (double)nevents) * 100;
    if (floor(percentage) != 0 && ((int)floor(percentage) % kStep) == 0)
    {
        std::cout << "\n"
                  << floor(percentage) << " %\t | \tProcessed " << evIdx + 1 << " events / " << nevents;
        kStep += 10;
    }
}

void evLoop(
    const std::string inputPath,
    TFile &outFile,
    const bool verbose,
    const std::string wd)
{
    //auto dmpch = aggregateEventsDmpChain(inputPath,verbose);
    auto dmpch = aggregateDataEventsTChain(inputPath, verbose);

    // Register Header container
    std::shared_ptr<DmpEvtHeader> evt_header = std::make_shared<DmpEvtHeader>();
    dmpch->SetBranchAddress("EventHeader", &evt_header);
    
    // Register BGO constainer
    std::shared_ptr<DmpEvtBgoHits> bgohits = std::make_shared<DmpEvtBgoHits>();
    dmpch->SetBranchAddress("DmpEvtBgoHits", &bgohits);

    // Register BGO REC constainer
    std::shared_ptr<DmpEvtBgoRec> bgorec = std::make_shared<DmpEvtBgoRec>();
    dmpch->SetBranchAddress("DmpEvtBgoRec", &bgorec);

    // Register STK container
    std::shared_ptr<TClonesArray> stkclusters = std::make_shared<TClonesArray>("DmpStkSiCluster");
    dmpch->SetBranchAddress("StkClusterCollection", &stkclusters);

    // Check if STK tracks collection exists
    bool fStkKalmanTracksFound = false;
    for (int brIdx = 0; brIdx < dmpch->GetListOfBranches()->GetEntries(); ++brIdx)
        if (strcmp(dmpch->GetListOfBranches()->At(brIdx)->GetName(), "StkKalmanTracks"))
        {
            fStkKalmanTracksFound = true;
            break;
        }

    // Register STK tracks collection
    std::shared_ptr<TClonesArray> stktracks = std::make_shared<TClonesArray>("DmpStkTrack", 200);
    if (fStkKalmanTracksFound)
        dmpch->SetBranchAddress("StkKalmanTracks", &stktracks);

    // Register PSD container
    std::shared_ptr<DmpEvtPsdHits> psdhits = std::make_shared<DmpEvtPsdHits>();
    dmpch->SetBranchAddress("DmpPsdHits", &psdhits);
    
    // Orbit filter
    // Set gIOSvc
    gIOSvc->Set("OutData/NoOutput", "True");
    gIOSvc->Initialize();
    // Create orbit filter
    std::unique_ptr<DmpFilterOrbit> pFilter = std::make_unique<DmpFilterOrbit>("EventHeader");
    // Activate orbit filter
    pFilter->ActiveMe(); // Call this function to calculate SAA through House Keeping Data

    // Event loop
    auto nevents = dmpch->GetEntries();
    if (verbose)
        std::cout << "\n\nTotal number of events: " << nevents << "\n\n";
    
    // STK charge histos
    TH1D h_chargeX("h_chargeX", "Charge distribution X", 100, 0, 1000);
    TH1D h_chargeY("h_chargeY", "Charge distribution Y", 100, 0, 1000);

    // Create and load acceptance events cuts from config file
    cuts_conf charge_cuts;
    data_active_cuts active_cuts;
    load_charge_struct(charge_cuts, active_cuts, wd);

    double _GeV = 0.001;
    int kStep = 10;
    unsigned int ev_counter = 0;

    for (unsigned int evIdx = 0; evIdx < nevents; ++evIdx)
    {
        // Get chain event
        dmpch->GetEvent(evIdx);

        if (pFilter->IsInSAA(evt_header->GetSecond()))
            continue;

        // Event printout
        if (verbose)
            updateProcessStatus(evIdx, kStep, nevents);

        // Get event total energy
        double bgoTotalE_raw = bgorec->GetTotalEnergy(); // Energy in MeV - not corrected
        double bgoTotalE = bgorec->GetElectronEcor();    // Returns corrected energy assuming this was an electron (MeV)
        
        // Don't accept events outside the selected energy window
        if (bgoTotalE * _GeV < charge_cuts.min_event_energy || bgoTotalE * _GeV > charge_cuts.max_event_energy)
            continue;

        // Read trigger status
        // For MC events triggers 1 and 2 are always disabled
        bool unbiased_tr = evt_header->GeneratedTrigger(0) && evt_header->EnabledTrigger(0);
        bool mip1_tr = evt_header->GeneratedTrigger(1);
        bool mip2_tr = evt_header->GeneratedTrigger(2);
        bool HET_tr = evt_header->GeneratedTrigger(3) && evt_header->EnabledTrigger(3);
        bool LET_tr = evt_header->GeneratedTrigger(4) && evt_header->EnabledTrigger(4);
        bool MIP_tr = mip1_tr || mip2_tr;

        bool general_trigger = MIP_tr || HET_tr || LET_tr;

        // Check if the event has been triggered or not
        if (general_trigger)
        {
            if (!checkBGOreco_data(bgorec))
                continue;
        }
        else
            continue;
        
        ++ev_counter;   
        best_track event_best_track;

        // Load BGO event class
        DmpBgoContainer bgoVault(DAMPE_bgo_nLayers);
        bgoVault.scanBGOHits(
            bgohits,
            bgoTotalE,
            DAMPE_bgo_nLayers);
        
        auto filter_BGO_fiducial_cut = false;
        auto filter_BGO_fiducial_maxElayer_cut = false;
        auto filter_BGO_fiducial_maxBarLayer_cut = false;
        auto filter_BGO_fiducial_BGOTrackContainment_cut = false;
        auto filter_nBarLayer13_cut = false;
        auto filter_maxRms_cut = false;
        auto filter_track_selection_cut = false;
        auto filter_all_cut = true;

        // Cut check...

        // **** BGO Fiducial Volume ****
        if (active_cuts.BGO_fiducial)
        {
            // maxElayer_cut
            filter_BGO_fiducial_maxElayer_cut = maxElayer_cut(
                bgorec,
                charge_cuts,
                bgoTotalE);

            // maxBarLayer_cut
            filter_BGO_fiducial_maxBarLayer_cut = maxBarLayer_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetiMaxLayer(),
                bgoVault.GetIdxBarMaxLayer());

            // BGOTrackContainment_cut
            filter_BGO_fiducial_BGOTrackContainment_cut = BGOTrackContainment_cut(
                bgorec,
                charge_cuts);

            filter_BGO_fiducial_cut = filter_BGO_fiducial_maxElayer_cut && filter_BGO_fiducial_maxBarLayer_cut && filter_BGO_fiducial_BGOTrackContainment_cut;
            filter_all_cut *= filter_BGO_fiducial_cut;
        }

        // **** nBarLayer13 cut ****
        if (active_cuts.nBarLayer13)
        {
            filter_nBarLayer13_cut = nBarLayer13_cut(
                bgohits,
                bgoVault.GetSingleLayerBarNumber(13),
                bgoTotalE);
            filter_all_cut *= filter_nBarLayer13_cut;
        }

        // **** maxRms cut ****
        if (active_cuts.maxRms)
        {
            filter_maxRms_cut = maxRms_cut(
                bgoVault.GetLayerBarNumber(),
                bgoVault.GetRmsLayer(),
                bgoTotalE,
                charge_cuts);
            filter_all_cut *= filter_maxRms_cut;
        }

        // **** track_selection cut ****
        filter_track_selection_cut = 
            track_selection_cut(
                bgorec,
                bgohits,
                stkclusters,
                stktracks,
                charge_cuts,
                event_best_track);
        filter_all_cut *= filter_track_selection_cut;
        
        if (active_cuts.nActiveCuts)
            if (filter_all_cut)
                fillChargeHistos(
                    h_chargeX, 
                    h_chargeY, 
                    event_best_track,
                    stkclusters);

    }

    if (verbose)
    {
        std::cout << "\n\n ****** \n\n";
        std::cout << "Triggered events: " << ev_counter << std::endl;

        auto refEntries = ev_counter;
        
        if (h_chargeX.GetEntries())
            std::cout << "X charge measurements: " << h_chargeX.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_chargeX.GetEntries()) / refEntries << std::endl;

        if (h_chargeY.GetEntries())
            std::cout << "Y charge measurements: " << h_chargeY.GetEntries() << "/" << refEntries << " | statistic efficiency: " << static_cast<double>(h_chargeY.GetEntries()) / refEntries << std::endl;

        std::cout << "\n\n ****** \n\n";
    }

    outFile.cd();

    h_chargeX.Write();
    h_chargeY.Write();
       
}