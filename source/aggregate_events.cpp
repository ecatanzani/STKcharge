#include "aggregate_events.h"

std::shared_ptr<DmpChain> aggregateEventsDmpChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using DAMPE Chain ******

    // Create DmpChain object
    std::shared_ptr<DmpChain> dmpch = std::make_shared<DmpChain>("CollectionTree");

    // Add MC file list to DmpChain
    //dmpch->AddFromList(getListPath(accInputPath, true).c_str());
    dmpch->AddFromList(listInputPath.c_str());
    if (verbose)
        dmpch->GetListOfFiles()->Print();

    return dmpch;
}

std::shared_ptr<TChain> aggregateEventsTChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    //TChain* dmpch = new TChain("CollectionTree");
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");
    //std::shared_ptr<TChain> dmpch( new TChain("CollectionTree") );

    // Reading list of MC files
    //std::ifstream input_file(getListPath(accInputPath, true).c_str());
    std::ifstream input_file(listInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << listInputPath << std::endl;
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    while (input_stream >> tmp_str)
    {
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }

    return dmpch;
}

std::shared_ptr<TChain> aggregateDataEventsTChain(
    const std::string listInputPath,
    const bool verbose)
{
    // ****** Access data using ROOT TChain ******

    // Create TChain object
    std::shared_ptr<TChain> dmpch = std::make_shared<TChain>("CollectionTree");
    
    // Reading list of MC files
    //std::ifstream input_file(getListPath(accInputPath, true).c_str());
    std::ifstream input_file(listInputPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << listInputPath << std::endl;
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    while (input_stream >> tmp_str)
    {
        // Setting up gIOSvc
        gIOSvc->Set("InData/Read", tmp_str.c_str());

        // Add file to the TChain
        dmpch->Add(tmp_str.c_str());
        if (verbose)
            std::cout << "\nAdding " << tmp_str << " to the chain ...";
    }

    return dmpch;
}