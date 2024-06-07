#ifndef HH_Parameters_HH
#define HH_Parameters_HH
#include <iosfwd>
#include <string>


//struct that defines the structure of the parameters

struct parameters
{ //analytical expression of the function R^2 --> R
  std::string f="8*_pi*_pi*sin(2*_pi*x)*sin(2*_pi*y)"; 
  //analytical expression of the boundary condition R^2 --> R
  std::string boundary="0";
  //number of subintervals in the grid for each dimension
  std::size_t n = 100;
  // tolerance for convergence 
  double eps = 1.e-6;
  // max number of iterations
  std::size_t itermax = 1000000;
};
//! Prints parameters
std::ostream &operator<<(std::ostream &, const parameters &);
#endif