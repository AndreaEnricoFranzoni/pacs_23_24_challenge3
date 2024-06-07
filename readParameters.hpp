#ifndef HH_READPARAMETERS_HH
#define HH_READPARAMETERS_HH
#include <string>
#include "parameters.hpp"

/*!
  Reads problem parameters from GetPot file
  @param filename The getopot file with the new values
  @param verbose Prints some information on the parameters
 */
parameters readParameters(std::string const &filename, bool verbose = false);

#endif