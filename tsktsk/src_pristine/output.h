#ifndef OUTPUT_H
#define OUTPUT_H

#include "hdf5.h"
#include <vector>

#define H5FILE_NAME        "thermalOutput.h5"
#define NUMERICALNAME "NumericalSolution"
#define ANALYTICALNAME "AnalyticalSolution"

#define RANK   2

int populate_file(int xmesh, int ymesh, std::vector<double> n, std::vector<double> a = std::vector<double>());

#endif