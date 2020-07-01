#include "charge.h"

void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY,
    TH1D &h_charge,
    TH2D &h_charge2D,
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters)
{
    double cluster_chargeX = -999;
    double cluster_chargeY = -999;
    auto btrack = track.myBestTrack;
    
    // Charge correction 
    auto track_correction = btrack.getDirection().CosTheta();
    
    // Define charge variables
    double chargeX = 0;
    double chargeY = 0;

    // Compute charges
    for(auto clIdx=0; clIdx<track.n_points; ++clIdx)
    {
        auto cluster_x = btrack.GetClusterX(clIdx, stkclusters.get());
        auto cluster_y = btrack.GetClusterY(clIdx, stkclusters.get());
        if (cluster_x && !cluster_x->getPlane())
            chargeX = cluster_x->getEnergy()*track_correction;
        if (cluster_y && !cluster_y->getPlane())
            chargeY = cluster_y->getEnergy()*track_correction;
    }

    // Compute mean charge
    auto mean_charge = 0.5*(chargeX + chargeY);
    
    // Fill histos
    h_chargeX.Fill(chargeX);
    h_chargeY.Fill(chargeY);
    h_charge.Fill(mean_charge);
    h_charge2D.Fill(chargeX, chargeY);

}