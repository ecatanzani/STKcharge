#include "charge.h"

extern void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY, 
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters)
{
    double cluster_chargeX = -999;
    double cluster_chargeY = -999;
    auto btrack = track.myBestTrack;
    auto track_correction = btrack.getDirection().CosTheta();

    for(auto clIdx=0; clIdx<track.n_points; ++clIdx)
    {
        auto cluster_x = btrack.GetClusterX(clIdx, stkclusters.get());
        auto cluster_y = btrack.GetClusterY(clIdx, stkclusters.get());
        if (cluster_x && !cluster_x->getPlane())
            h_chargeX.Fill(cluster_x->getEnergy()*track_correction);
        if (cluster_y && !cluster_y->getPlane())
            h_chargeY.Fill(cluster_y->getEnergy()*track_correction);
    }
}