#include "definitions.h"

void CustomExactSolution::derivatives(double x, double y, double& dx, double& dy) const
{
  dx = 0.;
  dy = 0.;
}

double CustomExactSolution::value(double x, double y) const
{
  return 1.;
}

Ord CustomExactSolution::ord(double x, double y) const
{
  return Ord(7);
}