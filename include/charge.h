#ifndef CHARGE_H
#define CHARGE_H

#include "data_cuts.h"

#include "TClonesArray.h"
#include "TH1D.h"

extern void load_charge_struct(
    cuts_conf &charge_cuts,
    data_active_cuts &active_cuts,
    const std::string wd);

extern void fillChargeHistos(
    TH1D &h_chargeX, 
    TH1D &h_chargeY, 
    const best_track track,
    const std::shared_ptr<TClonesArray> stkclusters);

#endif