#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomWeakFormPoissonNewton : public WeakForm<double>
{
public:
  CustomWeakFormPoissonNewton(std::string mat_al, Hermes::Hermes1DFunction<double>* lambda_al,
                              std::string mat_cu, Hermes::Hermes1DFunction<double>* lambda_cu,
                              Hermes::Hermes2DFunction<double>* vol_src_term, std::string bdy_heat_flux,
                              double alpha, double t_exterior);
};

/* Custom non-constant Dirichlet condition */

class CustomDirichletCondition : public EssentialBoundaryCondition<double>
{
public:
  CustomDirichletCondition(Hermes::vector<std::string> markers, double A, double B, double C);

  virtual EssentialBCValueType get_value_type() const;

  virtual double value(double x, double y) const;

  protected:
    double A, B, C;
};

