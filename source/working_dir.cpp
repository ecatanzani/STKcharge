#include "working_dir.h"

std::string getWorkingDir(const char* exePath)
{
    std::string tmpPath(exePath);
    std::size_t index = tmpPath.find("STKcharge");
    auto wd = tmpPath.substr(0, index + 9);
    return wd;
}

std::string GetCurrentWorkingDir(void)
{
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}