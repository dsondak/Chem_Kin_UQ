#ifndef LINEARINTERPOLATION_H
#define LINEARINTERPOLATION_H

#include <map>

class LinearInterpolation
{
  public:
    LinearInterpolation (std::map<double, double> Tdata_from_user);
    double Interpolate (double xstar_from_user);
    std::map<double, double> Tdata;
    double xstar;
};

#endif
