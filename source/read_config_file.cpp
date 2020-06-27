#include "data_cuts.h"

#include <fstream>
#include <sstream>

void load_charge_struct(
    cuts_conf &charge_cuts,
    data_active_cuts &active_cuts,
    const std::string wd)
{
    std::string cwd = wd;
    std::size_t index = cwd.find("STKcharge");
    std::string configPath = cwd.substr(0, index + 9);
    configPath += "/config/chargeConfig.txt";
    std::ifstream input_file(configPath.c_str());
    if (!input_file.is_open())
    {
        std::cerr << "\nERROR 100! File not open " << configPath << "\n\n";
        exit(100);
    }
    std::string input_string((std::istreambuf_iterator<char>(input_file)), (std::istreambuf_iterator<char>()));
    input_file.close();
    std::string tmp_str;
    std::istringstream input_stream(input_string);
    std::string::size_type sz;
    while (input_stream >> tmp_str)
    {
        // Load cuts variables
        if (!strcmp(tmp_str.c_str(), "min_event_energy"))
        {
            input_stream >> tmp_str;
            charge_cuts.min_event_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "max_event_energy"))
        {
            input_stream >> tmp_str;
            charge_cuts.max_event_energy = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "energy_lRatio"))
        {
            input_stream >> tmp_str;
            charge_cuts.energy_lRatio = stod(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "shower_axis_delta"))
        {
            input_stream >> tmp_str;
            charge_cuts.shower_axis_delta = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "max_rms_shower_width"))
        {
            input_stream >> tmp_str;
            charge_cuts.max_rms_shower_width = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_X_clusters"))
        {
            input_stream >> tmp_str;
            charge_cuts.track_X_clusters = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_Y_clusters"))
        {
            input_stream >> tmp_str;
            charge_cuts.track_Y_clusters = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_missingHit_X"))
        {
            input_stream >> tmp_str;
            charge_cuts.track_missingHit_X = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "track_missingHit_Y"))
        {
            input_stream >> tmp_str;
            charge_cuts.track_missingHit_Y = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_track"))
        {
            input_stream >> tmp_str;
            charge_cuts.STK_BGO_delta_track = stoi(tmp_str, &sz);
        }
        if (!strcmp(tmp_str.c_str(), "STK_BGO_delta_position"))
        {
            input_stream >> tmp_str;
            charge_cuts.STK_BGO_delta_position = stoi(tmp_str, &sz);
        }

        // Load cuts
        if (!strcmp(tmp_str.c_str(), "BGO_fiducial"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.BGO_fiducial = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "nBarLayer13"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.nBarLayer13 = true;
                ++active_cuts.nActiveCuts;
            }
        }
        if (!strcmp(tmp_str.c_str(), "maxRms"))
        {
            input_stream >> tmp_str;
            if (!strcmp(tmp_str.c_str(), "YES") || !strcmp(tmp_str.c_str(), "yes"))
            {
                active_cuts.maxRms = true;
                ++active_cuts.nActiveCuts;
            }
        }
    }
}
