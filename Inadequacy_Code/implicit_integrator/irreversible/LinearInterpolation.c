#include "LinearInterpolation.h"
#include <map>

LinearInterpolation::LinearInterpolation (std::map<double, double> Tdata_from_user) 
:
    Tdata(Tdata_from_user) // Initialize data container
{}

double LinearInterpolation::Interpolate (double xstar)
{
  /* Linearly interpolate function between two points
 *** Data is stored in map s.t. Tdata[xi] = Ti
 *** Iterate over the map to see if the point we're 
 *** interested in is exactly represented.  Otherwise 
 *** determine bounds around the point.  Then just 
 *** do the interpolation. */

  std::map<double, double>::iterator it = Tdata.begin(); // Put iterator at first entry
  bool found = false; 
  double ystar = it->second; // Initialize output

  // Loop over each entry in the data
  while (it != Tdata.end() && !found) {
      // Check if we're at the desired point
      if (it->first >= xstar) {
         found = true;
         break;
      }
      it++;
  }
  // If the desired point is the first entry, then just return
  if (it == Tdata.begin()) {
     ystar = Tdata.begin()->second;
  }
  // If the desired point was not found, then just return last point
  else if (it == Tdata.end()) {
     it--;
     ystar = it->second;
  }
  // If desired point was exactly found then no interpolation necessary
  else if (it->first == xstar) {
     ystar = it->second;
  }
  // Otherwise do the interpolation
  else {
     double x2 = it->first;   // Second point
     double y2 = it->second;  // Function at second point

     it--;
     double x1 = it->first;  // First point
     double y1 = it->second; // Function at first point

     double slope = (y2 - y1) / (x2 - x1); // Slope of line
     ystar = slope * (xstar - x1) + y1;    // Interpolated value
  }

  return ystar;
}
