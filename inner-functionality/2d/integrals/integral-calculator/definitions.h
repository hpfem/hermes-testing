#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class CustomExactSolution : public ExactSolutionScalar<double>
{
public:
  CustomExactSolution(MeshSharedPtr mesh) : ExactSolutionScalar<double>(mesh) {};

  virtual void derivatives(double x, double y, double& dx, double& dy) const;

  virtual double value(double x, double y) const;

  virtual Ord ord(double x, double y) const;

  MeshFunction<double>* clone() const {
    return new CustomExactSolution(this->mesh);
  }
};

class CustomVolumeIntegralCalculator : public Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>
{
public:
  CustomVolumeIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_function, number_of_integrals) {} ;
  CustomVolumeIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_functions, number_of_integrals) {} ;
  
  /// The integral description.
  /// \param[in] n - number of integration points.
  /// \param[in] result - preallocated (see number_of_integrals in the constructor) and zeroed array for the results.
  virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
  {
    for (int i = 0; i < n; i++)
    {
      result[0] = wt[i] * fns[0]->val[i];
      result[1] = wt[i] * fns[0]->val[i];
      result[2] = wt[i] * fns[0]->val[i];
    }
  };

  /// The integration order calculation.
  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result)
  {
    result[0] = fns[0]->val[0];
    result[1] = fns[0]->val[0];
    result[2] = fns[0]->val[0];
  }
};

class CustomSurfaceIntegralCalculator : public Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>
{
public:
  CustomSurfaceIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>(source_function, number_of_integrals) {} ;
  CustomSurfaceIntegralCalculator(Hermes::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>(source_functions, number_of_integrals) {} ;

  /// The integral description.
  /// \param[in] n - number of integration points.
  /// \param[in] result - preallocated (see number_of_integrals in the constructor) and zeroed array for the results.
  virtual void integral(int n, double* wt, Func<double> **fns, Geom<double> *e, double* result)
  {
    for (int i = 0; i < n; i++)
      result[0] = wt[i] * fns[0]->val[i] * e->nx[i];
  };

  /// The integration order calculation.
  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result)
  {
    result[0] = fns[0]->val[0];
  }
};