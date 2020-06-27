#include "myHeader.h"

#include <ctime>
#include <sstream>

const std::string uniqueOutFile(
    const std::string outputPath, 
    AnyOption &opt)
{
    std::time_t ctime = std::time(0);
    std::stringstream fPath;
    if (opt.getValue("outputDir") || opt.getValue('d'))
        fPath << outputPath << "/STKchargeOutFile_" << ctime << ".root";
    else if (opt.getValue("output") || opt.getValue('o'))
        fPath << outputPath;
    else
        fPath << "analysisOutFile_" << ctime << ".root";

    return fPath.str();
}