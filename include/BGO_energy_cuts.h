#ifndef BGO_ENERGY_CUTS_H
#define BGO_ENERGY_CUTS_H

#include <memory>
#include <vector>

#include "DmpEvtBgoRec.h"

#include "TH1D.h"

#include "data_cuts.h"
#include "DAMPE_geo_structure.h"

extern void evaluateEnergyRatio(
    const std::shared_ptr<DmpEvtBgoRec> bgorec,
    const cuts_conf acceptance_cuts,
    const double bgoTotalE,
    TH1D &h_layer_max_energy_ratio,
    std::vector<TH1D> &h_layer_energy_ratio);

#endif