#ifndef TRUTH_DATA_H
#define TRUTH_DATA_H

#include <hdf5.h>

class truth_data
{
  public:
    truth_data (const char* fname, int n_phis, int n_scen, int n_times, int num_fields);
    double * observation_data;
    double * ignition_data;
    double * sample_points;
    double * scenario_params;
};

#endif
