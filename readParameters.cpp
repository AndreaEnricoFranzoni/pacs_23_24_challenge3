#include <fstream>
#include "GetPot"
#include "readParameters.hpp"


parameters
readParameters(std::string const &filename, bool verbose)
{
  //#include "GetPot"
  // Parameter default constructor fills it with the defaults values
  parameters defaults;
  // checks if file exixts and is readable
  std::ifstream check(filename);
  if(!check)
    {
      std::cerr << "ERROR: Parameter file " << filename << " does not exist"
                << std::endl;
      std::cerr << "Reverting to default values." << std::endl;
      if(verbose)
        std::cout << defaults;
      check.close();
      return defaults;
    }
  else
    check.close();

  GetPot     ifile(filename.c_str());
  parameters values;
  // Read parameters from getpot ddata base
  values.f = ifile("f", defaults.f);
  values.boundary = ifile("boundary", defaults.boundary);
  values.n = ifile("n", defaults.n);
  values.eps = ifile("eps", defaults.eps);
  values.itermax = ifile("itermax", defaults.itermax);
  
  if(verbose)
    {
      std::cout << "PARAMETER VALUES IN GETPOT FILE"
                << "\n";
      ifile.print();
      std::cout << std::endl;
      std::cout << "ACTUAL VALUES"
                << "\n"
                << values;
    }
  return values;
  

}