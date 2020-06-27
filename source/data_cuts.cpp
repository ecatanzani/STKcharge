#include <data_cuts.h>

#include "DmpStkTrack.h"
#include "DmpStkSiCluster.h"
#include "DmpStkTrackHelper.h"

bool checkBGOreco_data(const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
    bool passed_BGO_check = true;
    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    if ((bgoRec_slope[0] == 0 && bgoRec_intercept[0] == 0) || (bgoRec_slope[1] == 0 && bgoRec_intercept[1] == 0))
        passed_BGO_check = false;
    
    return passed_BGO_check;
}

bool geometric_cut_data(const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
    bool passed_geometric_cut = false;

    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    double actual_topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double actual_topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    double actual_bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
    double actual_bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

    if (fabs(actual_topX) < BGO_SideXY && fabs(actual_topY) < BGO_SideXY &&
        fabs(actual_bottomX) < BGO_SideXY && fabs(actual_bottomY) < BGO_SideXY)
        passed_geometric_cut = true;
    
    return passed_geometric_cut;
}

bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts,
    const double bgoTotalE)
{
    bool passed_maxELayerTotalE_cut = true;
    auto layer_energies = bgorec->GetLayerEnergy();

    for (int idxLy = 0; idxLy < DAMPE_bgo_nLayers; ++idxLy)
    {
        auto tmp_ratio = layer_energies[idxLy] / bgoTotalE;
        if (tmp_ratio > data_cuts.energy_lRatio)
        {
            passed_maxELayerTotalE_cut = false;
            break;
        }
    }

    return passed_maxELayerTotalE_cut;
}

bool maxBarLayer_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<int> iMaxLayer,
    const std::vector<int> idxBarMaxLayer)
{
    bool passed_maxBarLayer_cut = true;

    for (auto lIdx = 1; lIdx <= 3; ++lIdx)
    {
        if (layerBarNumber[lIdx].size() == 0)
        {
            passed_maxBarLayer_cut = false;
            break;
        }
        if (iMaxLayer[lIdx] > -1)
            if (idxBarMaxLayer[lIdx] == 0 || idxBarMaxLayer[lIdx] == 21)
            {
                passed_maxBarLayer_cut = false;
                break;
            }
    }

    return passed_maxBarLayer_cut;
}

bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts)
{
    bool passed_bgo_containment_cut = false;

    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    double bottomX = bgoRec_slope[0] * BGO_BottomZ + bgoRec_intercept[0];
    double bottomY = bgoRec_slope[1] * BGO_BottomZ + bgoRec_intercept[1];

    if (
        fabs(topX) < data_cuts.shower_axis_delta &&
        fabs(topY) < data_cuts.shower_axis_delta &&
        fabs(bottomX) < data_cuts.shower_axis_delta &&
        fabs(bottomY) < data_cuts.shower_axis_delta)
        passed_bgo_containment_cut = true;

    return passed_bgo_containment_cut;
}

bool BGOTrackContainment_top_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts)
{
    bool passed_bgo_containment_cut = false;

    std::vector<double> bgoRec_slope(2);
    std::vector<double> bgoRec_intercept(2);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    if (fabs(topX) < data_cuts.shower_axis_delta && fabs(topY) < data_cuts.shower_axis_delta)
        passed_bgo_containment_cut = true;

    return passed_bgo_containment_cut;
}

bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> layerBarNumber,
    const double bgoTotalE)
{
    bool passed_nBarLayer13_cut = false;
    unsigned int nTriggeredBGO_13_bars = 0;
    double _GeV = 0.001;

    for (auto it = layerBarNumber.begin(); it != layerBarNumber.end(); ++it)
    {
        auto ihit = *it;
        auto hitE = (bgohits->fEnergy)[ihit];
        if (hitE > 10)
            ++nTriggeredBGO_13_bars;
    }
    double nBar13_threshold = 8 * log10(bgoTotalE * _GeV) - 5;
    if (nTriggeredBGO_13_bars < nBar13_threshold)
        passed_nBarLayer13_cut = true;

    return passed_nBarLayer13_cut;
}

bool maxRms_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<double> rmsLayer,
    const double bgoTotalE,
    const cuts_conf data_cuts)
{
    bool passed_maxRms_cut = false;
    auto max_rms = rmsLayer[0];
    auto eCut = bgoTotalE / 100.;

    for (auto lIdx = 0; lIdx < DAMPE_bgo_nLayers; ++lIdx)
    {
        double layerTotEnergy = 0;
        std::accumulate(layerBarNumber[lIdx].begin(), layerBarNumber[lIdx].end(), layerTotEnergy);
        if (layerTotEnergy > eCut)
            if (rmsLayer[lIdx] > max_rms)
                max_rms = rmsLayer[lIdx];
    }
    if (max_rms < data_cuts.max_rms_shower_width)
        passed_maxRms_cut = true;

    return passed_maxRms_cut;
}

inline void link_ladders(std::vector<int> &LadderToLayer)
{
    for (int ilad = 0; ilad < nSTKladders; ++ilad)
    {
        int iTRB = ilad / 24;
        int iladTRB = ilad % 24;
        int iPlane = 5 - iladTRB / 4;   // 0-5 for each double layer
        // int isX = (iTRB / 2) % 2;    // measuring X, TRB2, 3, 6, 7
        int isY = (iTRB / 2 + 1) % 2;   // measuring Y, TRB0, 1, 4, 5
        int iLay = iPlane * 2 + isY;    // 0-11 for each layer, starting from X (measuring X)
        LadderToLayer[ilad] = iLay;
    }
}

inline void fill_BGO_vectors(
    TVector3 &bgoRecEntrance,
    TVector3 &bgoRecDirection,
    const std::shared_ptr<DmpEvtBgoRec> bgorec)
{
    std::vector<double> bgoRec_slope(2, 0);
    std::vector<double> bgoRec_intercept(2, 0);

    bgoRec_slope[0] = bgorec->GetSlopeXZ();
    bgoRec_slope[1] = bgorec->GetSlopeYZ();
    bgoRec_intercept[0] = bgorec->GetInterceptXZ();
    bgoRec_intercept[1] = bgorec->GetInterceptYZ();

    // Build bgoRecDirection TVector3
    TVector3 vec_s0_a(bgoRec_intercept[0], bgoRec_intercept[1], 0.);
    TVector3 vec_s1_a(bgoRec_intercept[0] + bgoRec_slope[0], bgoRec_intercept[1] + bgoRec_slope[1], 1.);
    bgoRecDirection = (vec_s1_a - vec_s0_a).Unit(); //uni vector pointing from front to back

    // Build bgoRecEntrance TVector3
    double topZ = BGO_TopZ;
    double topX = bgoRec_slope[0] * BGO_TopZ + bgoRec_intercept[0];
    double topY = bgoRec_slope[1] * BGO_TopZ + bgoRec_intercept[1];

    if (fabs(topX) > BGO_SideXY || fabs(topY) > BGO_SideXY)
    {
        // possibly enter from the x-sides
        if (fabs(topX) > BGO_SideXY)
        {
            if (topX > 0)
                topX = BGO_SideXY;
            else
                topX = -BGO_SideXY;
            topZ = (topX - bgoRec_intercept[0]) / bgoRec_slope[0];
            topY = bgoRec_slope[1] * topZ + bgoRec_intercept[1];
            // possibly enter from the y-sides
            if (fabs(topY) > BGO_SideXY)
            {
                if (topY > 0)
                    topY = BGO_SideXY;
                else
                    topY = -BGO_SideXY;
                topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
                topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
            }
        }
        //enter from the y-sides
        else if (fabs(topY) > BGO_SideXY)
        {
            if (topY > 0)
                topY = BGO_SideXY;
            else
                topY = -BGO_SideXY;
            topZ = (topY - bgoRec_intercept[1]) / bgoRec_slope[1];
            topX = bgoRec_slope[0] * topZ + bgoRec_intercept[0];
        }
    }

    bgoRecEntrance[0] = topX;
    bgoRecEntrance[1] = topY;
    bgoRecEntrance[2] = topZ;
}

inline void get_track_points(
    DmpStkTrack *track,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::vector<int> LadderToLayer,
    std::vector<unsigned int> &track_nHoles,
    best_track &event_best_track,
    const bool best_track = false)
{
    std::vector<int> prevHole(2, - 2);
    std::vector<int> firstLayer(2, -1);
    std::vector<int> lastLayer(2, -1);
    std::vector<int> lastPoint(2, -1);
    std::vector<unsigned int> track_nHoles_cont(2, 0);
    
    // Loop on track points to find last layer values
    for (int ip = track->GetNPoints() - 1; ip >= 0; --ip)
    {
        if (lastLayer[0] == -1)
        {
            if (track->getHitMeasX(ip) > -99999)
            {
                lastPoint[0] = ip;
                DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                lastLayer[0] = LadderToLayer[hardID];
            }
        }
        if (lastLayer[1] == -1)
        {
            if (track->getHitMeasY(ip) > -99999)
            {
                lastPoint[1] = ip;
                DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
                auto hardID = cluster->getLadderHardware();
                lastLayer[1] = LadderToLayer[hardID];
            }
        }
    }
    
    // Found the number of holes on both X and Y
    for (int ip = 0; ip <= lastPoint[0]; ++ip)
    {
        if (track->getHitMeasX(ip) > -99999)
        {
            DmpStkSiCluster *cluster = track->GetClusterX(ip, stkclusters.get());
            auto hardID = cluster->getLadderHardware();
            if (firstLayer[0] == -1)
                firstLayer[0] = LadderToLayer[hardID];
        }
        else
        {
            if (firstLayer[0] != -1)
                ++track_nHoles[0];
            if (ip == prevHole[0] + 1)
                ++track_nHoles_cont[0];
            prevHole[0] = ip;
        }
    }
    
    for (int ip = 0; ip <= lastPoint[1]; ++ip)
    {
        if (track->getHitMeasY(ip) > -99999)
        {
            DmpStkSiCluster *cluster = track->GetClusterY(ip, stkclusters.get());
            auto hardID = cluster->getLadderHardware();
            if (firstLayer[1] == -1)
                firstLayer[1] = LadderToLayer[hardID];
        }
        else
        {
            if (firstLayer[1] != -1)
                ++track_nHoles[1];
            if (ip == prevHole[1] + 1)
                ++track_nHoles_cont[1];
            prevHole[1] = ip;
        }
    }

    if (best_track)
    {
        event_best_track.n_points = track->GetNPoints();
        for (int idx = 0; idx < 2; ++idx)
            event_best_track.n_holes[idx] = track_nHoles[idx];
        event_best_track.track_slope[0] = track->getTrackParams().getSlopeX();
        event_best_track.track_slope[1] = track->getTrackParams().getSlopeY();
        event_best_track.track_intercept[0] = track->getTrackParams().getInterceptX();
        event_best_track.track_intercept[1] = track->getTrackParams().getInterceptY();
        event_best_track.track_direction = (track->getDirection()).Unit();
    }
}

bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const cuts_conf data_cuts,
    best_track &event_best_track)
{
    bool passed_track_selection_cut = false;

    TVector3 bgoRecEntrance;
    TVector3 bgoRecDirection;
    std::vector<int> LadderToLayer(nSTKladders, -1);
    std::vector<DmpStkTrack *> selectedTracks;

    link_ladders(LadderToLayer);
    fill_BGO_vectors(
        bgoRecEntrance,
        bgoRecDirection,
        bgorec);
    
    // Loop on the tracks
    for (int trIdx = 0; trIdx < stktracks->GetLast() + 1; ++trIdx)
    {
        std::vector<unsigned int> track_nHoles(2, 0);
        std::vector<double> track_slope(2, 0);
        std::vector<double> track_intercept(2, 0);
        std::vector<double> extr_BGO_top(2, 0);

        // *********************

        // Get the track
        auto track = static_cast<DmpStkTrack *>(stktracks->ConstructedAt(trIdx));
        
        // Reject tracks with not enough X and Y clusters
        if (track->getNhitX() < data_cuts.track_X_clusters || track->getNhitY() < data_cuts.track_Y_clusters)
            continue;
        
        get_track_points(
            track,
            stkclusters,
            LadderToLayer,
            track_nHoles,
            event_best_track);
        
        if (track_nHoles[0] > 1 || track_nHoles[1] > 1)
            continue;

        // Find slope and intercept
        track_slope[0] = track->getTrackParams().getSlopeX();
        track_slope[1] = track->getTrackParams().getSlopeY();
        track_intercept[0] = track->getTrackParams().getInterceptX();
        track_intercept[1] = track->getTrackParams().getInterceptY();
        TVector3 trackDirection = (track->getDirection()).Unit();

        // Extrapolate to the top of BGO
        for (int coord = 0; coord < 2; ++coord)
            extr_BGO_top[coord] = track_slope[coord] * BGO_TopZ + track_intercept[coord];

        // Evaluate distance between Top STK and BGO points
        double dxTop = extr_BGO_top[0] - bgoRecEntrance[0];
        double dyTop = extr_BGO_top[1] - bgoRecEntrance[1];
        double drTop = sqrt(pow(dxTop, 2) + pow(dyTop, 2));
        // Evaluate angular distance between STK track and BGO Rec track
        double dAngleTrackBgoRec = trackDirection.Angle(bgoRecDirection) * TMath::RadToDeg();

        if (drTop > data_cuts.STK_BGO_delta_position)
            continue;
        if (dAngleTrackBgoRec > data_cuts.STK_BGO_delta_track)
            continue;

        selectedTracks.push_back(track);
    }
    
    // Sort selected tracks vector
    DmpStkTrackHelper tHelper(stktracks.get(), true, bgorec.get(), bgohits.get());
    tHelper.MergeSort(selectedTracks, &DmpStkTrackHelper::TracksCompare);

    if (selectedTracks.size() > 0)
    {
        DmpStkTrack *selected_track = static_cast<DmpStkTrack *>(selectedTracks[0]);
        std::vector<unsigned int> track_nHoles(2, 0);

        // Fill best track structure
        get_track_points(
            selected_track,
            stkclusters,
            LadderToLayer,
            track_nHoles,
            event_best_track,
            true);

        event_best_track.extr_BGO_topX = event_best_track.track_slope[0] * BGO_TopZ + event_best_track.track_intercept[0];
        event_best_track.extr_BGO_topY = event_best_track.track_slope[1] * BGO_TopZ + event_best_track.track_intercept[1];

        event_best_track.STK_BGO_topX_distance = event_best_track.extr_BGO_topX - bgoRecEntrance[0];
        event_best_track.STK_BGO_topY_distance = event_best_track.extr_BGO_topY - bgoRecEntrance[1];
        event_best_track.angular_distance_STK_BGO = event_best_track.track_direction.Angle(bgoRecDirection) * TMath::RadToDeg();
        event_best_track.STK_BGO_topY_distance = sqrt(pow(event_best_track.STK_BGO_topX_distance, 2) + pow(event_best_track.STK_BGO_topY_distance, 2));
        
        event_best_track.myBestTrack = static_cast<DmpStkTrack>(*selected_track);

        passed_track_selection_cut = true;
    }
    
    return passed_track_selection_cut;
}