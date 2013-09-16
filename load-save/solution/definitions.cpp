#include "definitions.h"

void CustomInitialCondition::derivatives (double x, double y, double& dx, double& dy) const 
{
  dx = (y + 10)/100.;
  dy = (x + 10)/100.;
};

double CustomInitialCondition::value (double x, double y) const 
{
  return (x + 10) * (y + 10) / 100.;
}

Ord CustomInitialCondition::ord(double x, double y) const 
{
  return Ord(x*y);
}
  
MeshFunction<double>* CustomInitialCondition::clone() const
{
  return new CustomInitialCondition(mesh);
}

CustomInitialCondition::~CustomInitialCondition()
{
}