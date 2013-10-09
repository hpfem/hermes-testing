#include "definitions.h"
#include "../../../../testing-core/testing-core.h"

//  This is a test of integral calculator.
int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh));

  CustomVolumeIntegralCalculator volume_integral(exact_sln, 1);
  CustomSurfaceIntegralCalculator surface_integral(exact_sln, 1);

  double volume = volume_integral.calculate("Copper")[0];
  double surface = surface_integral.calculate("Outer")[0];

  bool success = Testing::test_value(volume, 1.01039, "volumetric", 1e-3) && Testing::test_value(surface, -2, "surface", 1e-3);

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