#ifndef DATA_CUTS_H
#define DATA_CUTS_H

#include <memory>
#include <vector>
#include <numeric>

#include "DAMPE_geo_structure.h"

#include "TClonesArray.h"
#include "TH1D.h"
#include "TH2D.h"

// DAMPESW includes
#include "DmpEvtBgoRec.h"
#include "DmpEvtBgoHits.h"
#include "DmpEvtPsdHits.h"
#include "DmpStkTrack.h"

struct cuts_conf
{
    double min_event_energy;
    double max_event_energy;
    double energy_lRatio;
    int shower_axis_delta;
    double vertex_radius;
    int max_rms_shower_width;
    int track_X_clusters;
    int track_Y_clusters;
    int track_missingHit_X;
    int track_missingHit_Y;
    int STK_BGO_delta_track;
    int STK_BGO_delta_position;
};

struct data_active_cuts
{
    bool BGO_fiducial = false;
    bool nBarLayer13 = false;
    bool maxRms = false;
    bool track_selection = false;
    unsigned int nActiveCuts = 0;
};

struct best_track
{
    unsigned int n_points = 0;
    std::vector<unsigned int> n_holes = {0, 0};
    std::vector<double> track_slope = {-999, -999};
    std::vector<double> track_intercept = {-999, -999};
    TVector3 track_direction;
    double extr_BGO_topX = -999;
    double extr_BGO_topY = -999;
    double STK_BGO_topX_distance = -999;
    double STK_BGO_topY_distance = -999;
    double angular_distance_STK_BGO = -999;
    DmpStkTrack myBestTrack;
};

extern bool checkBGOreco_data(const std::shared_ptr<DmpEvtBgoRec> bgorec);
extern bool geometric_cut_data(const std::shared_ptr<DmpEvtBgoRec> bgorec);

extern void evaluateTopBottomPosition_data(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    TH1D &h_BGOrec_slopeX,
    TH1D &h_BGOrec_slopeY,
    TH1D &h_BGOrec_interceptX,
    TH1D &h_BGOrec_interceptY,
    TH2D &h_BGOreco_topMap,
    TH2D &h_BGOreco_bottomMap);

extern bool maxElayer_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf acceptance_cuts,
    const double bgoTotalE);

extern bool maxBarLayer_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<int> iMaxLayer,
    const std::vector<int> idxBarMaxLayer);

extern bool BGOTrackContainment_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts);

extern bool BGOTrackContainment_top_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf data_cuts);

extern bool nBarLayer13_cut(
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::vector<short> layerBarNumber,
    const double bgoTotalE);

extern bool maxRms_cut(
    const std::vector<std::vector<short>> layerBarNumber,
    const std::vector<double> rmsLayer,
    const double bgoTotalE,
    const cuts_conf data_cuts);

extern bool track_selection_cut(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const std::shared_ptr<DmpEvtBgoHits> bgohits,
    const std::shared_ptr<TClonesArray> stkclusters,
    const std::shared_ptr<TClonesArray> stktracks,
    const cuts_conf data_cuts,
    best_track &event_best_track);

#endif