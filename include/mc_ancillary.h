#ifndef MC_ANCILLARY_H
#define MC_ANCILLARY_H

#include <vector>

#include "data_cuts.h"

extern void compute_proton_background(
    const double sumRms,
    const std::vector<double> fracLayer,
    const cuts_conf data_cuts,
    const double energy,
    TH1D &h_background_under_xtrl_cut,
    TH1D &h_background_over_xtrl_cut);

#endif