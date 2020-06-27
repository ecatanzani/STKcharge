#ifndef AGGREGATE_EVENTS_H
#define AGGREGATE_EVENTS_H

#include <memory>
#include <string>
#include <fstream>
#include <sstream>

#include "DmpChain.h"
#include "DmpIOSvc.h"
#include "DmpCore.h"
#include "DmpFilterOrbit.h"

#include "TChain.h"

extern std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string listInputPath,
    const bool verbose);

extern std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string listInputPath,
    const bool verbose);

extern std::shared_ptr<TChain> aggregateDataEventsTChain(
    const std::string listInputPath,
    const bool verbose);

#endif