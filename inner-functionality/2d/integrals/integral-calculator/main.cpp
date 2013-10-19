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

  int random_number = std::min<int>(10000, std::max<int>(2, std::rand()));

  CustomVolumeIntegralCalculator volume_integral(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln, exact_sln), random_number);
  CustomSurfaceIntegralCalculator surface_integral(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln, exact_sln), random_number);

  double volume = volume_integral.calculate("Copper")[random_number / 2];
  double surface = surface_integral.calculate("Outer")[random_number / 2];

  bool success = Testing::test_value(volume, 1.3927, "volumetric", 1e-3);
  success = Testing::test_value(surface, 2., "surface", 1e-3) && success;

  mloader.load("mesh_marker3.xml", mesh);

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln2(new CustomExactSolution(mesh));
  CustomVolumeIntegralCalculator volume_integral2(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln2, exact_sln2), random_number);
  CustomSurfaceIntegralCalculator surface_integral2(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln2, exact_sln2), random_number);

  Hermes::vector<double> volumetric_int;
  Hermes::vector<double> surface_int;
  try
  {
    volumetric_int.push_back(volume_integral2.calculate(Hermes::vector<std::string>("1", "2"))[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

  try
  {
    volumetric_int.push_back(volume_integral2.calculate(Hermes::vector<std::string>("1", "3"))[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

  try
  {
    volumetric_int.push_back(volume_integral2.calculate(Hermes::vector<std::string>("2", "4"))[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

  try
  {
    volumetric_int.push_back(volume_integral2.calculate(Hermes::vector<std::string>("3", "4", "5"))[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

  try
  {
    Hermes::vector<std::string> markers;
    char t[2];
    for(int i = 0; i < 29; i++)
      markers.push_back(itoa(i, t, 10));
    surface_int.push_back(surface_integral2.calculate(markers)[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }
  try
  {
    surface_int.push_back(surface_integral2.calculate("2")[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }
  try
  {
    surface_int.push_back(surface_integral2.calculate("3")[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }
  try
  {
    surface_int.push_back(surface_integral2.calculate("4")[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }
  try
  {
    surface_int.push_back(surface_integral2.calculate("5")[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

  // Load the mesh.
  MeshSharedPtr mesh_whole_domain(new Mesh), mesh_with_hole(new Mesh);
  Hermes::vector<MeshSharedPtr> meshes (mesh_whole_domain, mesh_with_hole);
  mloader.load("subdomains.xml", meshes);

  mesh_with_hole->refine_towards_boundary("Inner Wall");
  mesh_whole_domain->refine_towards_boundary("Inner Wall");

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln_whole(new CustomExactSolution(mesh_whole_domain));
  MeshFunctionSharedPtr<double> exact_sln_hole(new CustomExactSolution(mesh_with_hole));
  CustomVolumeIntegralCalculator volume_integral3(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln_whole, exact_sln_hole), random_number);
  CustomSurfaceIntegralCalculator surface_integral3(Hermes::vector<MeshFunctionSharedPtr<double> >(exact_sln_hole, exact_sln_whole), random_number);

  try
  {
    surface_int.push_back(surface_integral3.calculate("Inner Wall")[random_number / 2]);
    surface_int.push_back(surface_integral3.calculate("Inlet")[random_number / 2]);
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
    success = false;
  }

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