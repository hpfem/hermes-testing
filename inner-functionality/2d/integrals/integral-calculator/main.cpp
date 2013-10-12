#include "definitions.h"
#include "../../../../testing-core/testing-core.h"
#include <algorithm>

//  This is a test of integral calculator.
int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh));

#ifdef _SHOW_OUTPUT
  ScalarView s;
  MeshView m;
  s.show(exact_sln);
  m.show(mesh);
  View::wait();
#endif

  int random_number = std::max(2, std::rand());

  CustomVolumeIntegralCalculator volume_integral(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln, exact_sln), random_number);
  CustomSurfaceIntegralCalculator surface_integral(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln, exact_sln), random_number);

  double volume = volume_integral.calculate("Copper")[random_number / 2];
  double surface = surface_integral.calculate("Outer")[random_number / 2];

  bool success = Testing::test_value(volume, 1.3927, "volumetric", 1e-3);
  success = Testing::test_value(surface, 2., "surface", 1e-3) && success;

  if (success)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }
}