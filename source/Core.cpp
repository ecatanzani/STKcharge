#include "myHeader.h"
#include "data_loop.h"

#include "TFile.h"

void STKCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    AnyOption &opt,
    const std::string wd)
{
    // Create output TFile
    const char* outFilePath = static_cast<const char*>(uniqueOutFile(outputPath, opt).c_str());
    
    TFile outFile(outFilePath, "NEW", "Analysis Output File");
    if (!outFile.IsOpen())
    {
        std::cerr << "\n\nError writing output TFile: " << outFilePath << std::endl;
        exit(123);
    }

    evLoop(
        inputPath,
        outFile,
        verbose,
        wd);

    // Close output file ...
    outFile.Close();
}