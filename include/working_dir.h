#ifndef WORKING_DIR_H
#define WORKING_DIR_H

#include <string>
#include <unistd.h>

#define GetCurrentDir getcwd

extern std::string getWorkingDir(const char *exePath);
extern std::string GetCurrentWorkingDir(void);

#endif