#include "parameters.hpp"
#include <iostream>
std::ostream &
operator<<(std::ostream &out, const parameters &p)
{
  out << "PARAMETER VALUES:"
      << "\n";
  out << "f= " << p.f << "\n";
  out << "boundary= " << p.boundary << "\n";
  out << "n= " << p.n << "\n";
  out << "eps= " << p.eps << "\n";
  out << "itermax= " << p.itermax << "\n";
  return out;
}