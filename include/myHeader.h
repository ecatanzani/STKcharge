#ifndef MYHEADER_H
#define MYHEADER_H

#include <iostream>
#include <string>
#include <vector>

#include "anyoption.h"

#pragma once

extern void STKCore(
    const std::string inputPath,
    const std::string outputPath,
    const bool verbose,
    AnyOption &opt,
    const std::string wd);

extern const std::string uniqueOutFile(
    const std::string outputPath,
    AnyOption &opt);

#endif