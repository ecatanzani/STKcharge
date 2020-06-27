#include "BGO_energy_cuts.h"

void evaluateEnergyRatio(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf acceptance_cuts,
    const double bgoTotalE,
    TH1D &h_layer_max_energy_ratio,
    std::vector<TH1D> &h_layer_energy_ratio)
{
    int iMaxELayer = -1;  // Index of the layer corresponding to the max energy
    double MaxELayer = 0; // Value of the max energy
    double _GeV = 0.001;

    // Found the max energy value and layer
    for (int idxLy = 0; idxLy < DAMPE_bgo_nLayers; ++idxLy)
    {
        auto layer_energy = static_cast<double>((bgorec->GetLayerEnergy())[idxLy]);
        h_layer_energy_ratio[idxLy].Fill(layer_energy / bgoTotalE);
        if (layer_energy > MaxELayer)
        {
            MaxELayer = layer_energy;
            iMaxELayer = idxLy;
        }
    }

    auto rMaxELayerTotalE = MaxELayer / bgoTotalE;
    h_layer_max_energy_ratio.Fill(rMaxELayerTotalE);
}