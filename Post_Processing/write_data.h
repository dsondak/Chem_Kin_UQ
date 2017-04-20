#ifndef __WRITE_DATA_H__
#define __WRITE_DATA_H__

#include <vector>

int create_file(const char* fname, int rows, int cols);
int write_file(std::vector<double>& fx, std::vector<double>& fT, const char* fname, int rows, int cols);


#endif
