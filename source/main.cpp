#include "myHeader.h"
#include "working_dir.h"

#include <sstream>

int main(int argc, char **argv)
{

    AnyOption opt;

    opt.addUsage("Usage: ");
    opt.addUsage("");
    opt.addUsage(" -h  --help                                                   Prints this help");
    opt.addUsage(" -i  --input          <path_to_input_DATA_list>           (*) Input DATA list");
    opt.addUsage(" -o  --output         <path_to_output_TFile>                  Output ROOT TFile");
    opt.addUsage(" -d  --outputDir      <path_to_output_TFile_dir>              Output ROOT TFile directory");
    opt.addUsage(" -v  --verbose                                                Verbose output");
    opt.addUsage("");

    opt.setFlag("help", 'h');
    opt.setOption("input", 'i');
    opt.setOption("output", 'o');
    opt.setOption("outputDir", 'd');
    opt.setFlag("verbose", 'v');

    opt.processCommandArgs(argc, argv);

    /*
        Input variables

    */

    std::string inputPath;
    std::string outputPath;
    std::string wd = getWorkingDir(argv[0]);

    bool verbose = false;
    
    if (!opt.hasOptions())
        opt.printUsage();

    if (opt.getFlag("help") || opt.getFlag('h'))
    {
        opt.printUsage();
        return 0;
    }
    if (opt.getValue("input") || opt.getValue('i'))
        inputPath = opt.getValue('i');
    if (opt.getValue("output") || opt.getValue('o'))
        outputPath = opt.getValue('o');
    if (opt.getValue("outputDir") || opt.getValue('d'))
        outputPath = opt.getValue('d');
    if (opt.getFlag("verbose") || opt.getFlag('v'))
        verbose = opt.getFlag('v');
    
    STKCore(
        inputPath,
        outputPath,
        verbose,
        opt,
        wd);

    return 0;
}