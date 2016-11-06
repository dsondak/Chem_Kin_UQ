#ifndef __WRITE_DATA_H__
#define __WRITE_DATA_H__

#include <vector>

int create_file(const char* fname, int rows, int cols, int num_scen);
int write_file(std::vector<double>& time, std::vector<double>& sols, const char* fname, int rows, int cols, int scen, double phi, double scen_par);


#endif
