#include "hermes2d.h"
#include "../../../../testing-core/testing-core.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

class MyVolumetricIntegralCalculator : public Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>
{
public:
  MyVolumetricIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MyVolumetricIntegralCalculator(std::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Func<double> **fns, GeomVol<double> *e, double* result)
  {
    for(int j = 0; j < this->number_of_integrals; j++)
      for (int i = 0; i < n; i++)
        result[j] += wt[i];
  };

  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
    for(int j = 0; j < this->number_of_integrals; j++)
      result[j] = Hermes::Ord(j);
  }
};


class MyVolumetricIntegralCalculatorL2 : public Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>
{
public:
  MyVolumetricIntegralCalculatorL2(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MyVolumetricIntegralCalculatorL2(std::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : Hermes2D::PostProcessing::VolumetricIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Func<double> **fns, GeomVol<double> *e, double* result)
  {
    for(int j = 0; j < this->number_of_integrals; j++)
      for (int i = 0; i < n; i++)
        result[j] += wt[i] * fns[j]->val[i] * fns[j]->val[i];
  };

  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
    for(int j = 0; j < this->number_of_integrals; j++)
      result[j] = fns[j]->val[0];
  }
};

class MySurfaceIntegralCalculator : public Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>
{
public:
  MySurfaceIntegralCalculator(MeshFunctionSharedPtr<double> source_function, int number_of_integrals) : Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>(source_function, number_of_integrals)
  {
  }

  MySurfaceIntegralCalculator(std::vector<MeshFunctionSharedPtr<double> > source_functions, int number_of_integrals) : Hermes2D::PostProcessing::SurfaceIntegralCalculator<double>(source_functions, number_of_integrals)
  {
  }

  virtual void integral(int n, double* wt, Func<double> **fns, GeomSurf<double> *e, double* result)
  {
    for(int j = 0; j < this->number_of_integrals; j++)
      for (int i = 0; i < n; i++)
        result[j] += wt[i];
  };

  virtual void order(Func<Hermes::Ord> **fns, Hermes::Ord* result) {
    for(int j = 0; j < this->number_of_integrals; j++)
      result[j] = Hermes::Ord(j);
  }
};

void handle_success(bool success)
{
  if(success == true)
    printf("Success!\n");
  else
    printf("Failure!\n");
}

#ifdef _SHOW_OUTPUT
Views::MeshView m;
#endif
std::vector<int> measured_elems;
std::vector<int> elems;
std::vector<double> measured_integrals;
std::vector<double> integrals;

bool part2OneCase(int i, MeshSharedPtr mesh, std::vector<std::string> markers, int layers, int poly_degree)
{
  bool success = true;
  try
  {
    Hermes::Hermes2D::MeshSharedPtr egg_shell_mesh = Hermes::Hermes2D::EggShell::get_egg_shell(mesh, markers, layers);
    elems.push_back(egg_shell_mesh->get_num_active_elements());
    success = Hermes::Testing::test_value(egg_shell_mesh->get_num_active_elements(), measured_elems[i], "Num elems", 1) && success;

#ifdef _SHOW_OUTPUT
    m.show(egg_shell_mesh);
    m.wait_for_keypress();
#endif
    Hermes::Hermes2D::MeshFunctionSharedPtr<double> eggShell(new Hermes::Hermes2D::ExactSolutionEggShell(egg_shell_mesh, poly_degree));
#ifdef _SHOW_OUTPUT
    m.show(egg_shell_mesh);
    m.wait_for_keypress();
#endif

    Hermes::Hermes2D::MeshFunctionSharedPtr<double> zeroSlnWholeMesh(new Hermes::Hermes2D::ConstantSolution<double>(mesh, 1341.13213));
    Hermes::Hermes2D::MeshFunctionSharedPtr<double> zeroSlnEggShellMesh(new Hermes::Hermes2D::ConstantSolution<double>(egg_shell_mesh, 1.575674));

    MyVolumetricIntegralCalculatorL2 integral(std::vector<MeshFunctionSharedPtr<double> >({ eggShell, zeroSlnEggShellMesh, zeroSlnWholeMesh }), 3);
    double* result = integral.calculate(HERMES_ANY);
    success = Hermes::Testing::test_value(result[0], measured_integrals[i*3 + 0], "Integral-0") && success;
    success = Hermes::Testing::test_value(result[1], measured_integrals[i*3 + 1], "Integral-1") && success;
    success = Hermes::Testing::test_value(result[2], measured_integrals[i*3 + 2], "Integral-2") && success;
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    std::cout << e.what();
    return false;
  }
  return success;
}

bool part2OneCase(int i, MeshSharedPtr mesh, std::string marker, int layers, int poly_degree)
{
  std::vector<std::string> markers;
  markers.push_back({marker});
  return part2OneCase(i, mesh, markers, layers, poly_degree);
}

void part2FillData()
{
  measured_elems.push_back(82);
  measured_elems.push_back(130);
  measured_elems.push_back(88);
  measured_elems.push_back(57);
  measured_elems.push_back(82);
  measured_elems.push_back(425);
  measured_elems.push_back(379);
  measured_elems.push_back(271);
  measured_elems.push_back(293);
  measured_elems.push_back(244);

  measured_integrals.push_back(1.167631128e-005);
  measured_integrals.push_back(0.001440312967);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.0002445754128);
  measured_integrals.push_back(0.002081465591);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.0001259342804);
  measured_integrals.push_back(0.001343521083);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(6.952705348e-005);
  measured_integrals.push_back(0.0007275969098);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(7.275436683e-006);
  measured_integrals.push_back(0.001175314496);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.001016728677);
  measured_integrals.push_back(0.009732970719);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.001070668838);
  measured_integrals.push_back(0.009362072818);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.0005411781906);
  measured_integrals.push_back(0.004710499049);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.0006080980846);
  measured_integrals.push_back(0.005094778685);
  measured_integrals.push_back(134897.6543);
  measured_integrals.push_back(0.0003635958382);
  measured_integrals.push_back(0.00377022343);
  measured_integrals.push_back(134897.6543);
}

int main(int argc, char* argv[])
{
  Hermes::Hermes2D::EggShell::egg_shell_verbose = false;
  // Part1: eggShellCreation + Integral (both Vol. + Surf.) comparison to calculation with for_all_active_elements.
  try
  {
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("mesh_marker3.xml", mesh);

    Hermes::Hermes2D::MeshSharedPtr egg_shell_mesh = Hermes::Hermes2D::EggShell::get_egg_shell(mesh,"3", 2);

#ifdef SHOW_OUTPUT
    Views::MeshView m;
    m.show(egg_shell_mesh);
    m.wait_for_keypress();
#endif

    Hermes::Hermes2D::MeshFunctionSharedPtr<double> eggShell(new Hermes::Hermes2D::ExactSolutionEggShell(egg_shell_mesh, 3));

    MyVolumetricIntegralCalculator volume_integrator(eggShell, 15);
    double* area_int = volume_integrator.calculate(HERMES_ANY);

    MySurfaceIntegralCalculator surface_integrator(eggShell, 34);
    double* surface_1_int = surface_integrator.calculate(EggShell::eggShell1Marker);

    Element* e;
    double area = 0;
    double surface_1 = 0;
    Hermes::Hermes2D::Mesh::MarkersConversion::IntValid internalMarker = egg_shell_mesh->get_boundary_markers_conversion().get_internal_marker(EggShell::eggShell1Marker);
    if(!internalMarker.valid)
    {
      handle_success(false);
      return -1;
    }
    int innerEggShell1Marker = internalMarker.marker;
    for_all_active_elements(e, egg_shell_mesh)
    {
      e->calc_area(true);
      area += e->area;
      for(int i = 0; i < e->get_nvert(); i++)
      {
        if(e->en[i]->marker == innerEggShell1Marker)
          surface_1 += std::sqrt(((e->vn[(i + 1) % e->get_nvert()]->x - e->vn[i]->x) * (e->vn[(i + 1) % e->get_nvert()]->x - e->vn[i]->x)) + ((e->vn[(i + 1) % e->get_nvert()]->y - e->vn[i]->y) * (e->vn[(i + 1) % e->get_nvert()]->y - e->vn[i]->y)));
      }
    }

    bool success = true;
    for(int i = 0; i < 15; i++)
      success = Testing::test_value(area_int[i], area, "Area") && success;
    for(int i = 0; i < 34; i++)
      success = Testing::test_value(surface_1_int[i], surface_1, "Surface") && success;

    handle_success(success);
    if(!success)
      return -1;

#ifdef SHOW_OUTPUT
    Views::ScalarView s;
    s.show(eggShell);
    s.wait_for_keypress();
#endif
  }
  catch(std::exception& e)
  {
    handle_success(false);
    return -1;
  }

  // Part2: eggShellCreation with various settings.
  part2FillData();

  {
    bool success = true;

    EggShell::egg_shell_verbose = false;

    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("pokus.xml", mesh);

    success = part2OneCase(1, mesh, "2", 2, 3) && success;

    success = part2OneCase(2, mesh, "3", 3, 4) && success;

    success = part2OneCase(3, mesh, "4", 2, 5) && success;

    success = part2OneCase(5, mesh, std::vector<std::string>({ "2", "4" }), 5, 6) && success;

    success = part2OneCase(6, mesh, std::vector<std::string>({ "3", "5", "1" }), 4, 4) && success;

    success = part2OneCase(7, mesh, std::vector<std::string>({ "3", "2", "2" }), 3, 2) && success;

    success = part2OneCase(8, mesh, std::vector<std::string>({ "1", "2" }), 2, 3) && success;

    success = part2OneCase(9, mesh, std::vector<std::string>({ "3", "4" }), 4, 1) && success;

    if(!success)
    {
      handle_success(success);
      return (!success);
    }
  }
  // Part 3.
  {
#ifdef _SHOW_OUTPUT
    MeshView m;
    ScalarView s;
#endif
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.load("pokus.xml", mesh);

    for(int k = 0; k < 5; k++)
    {
      std::cout << "k == " << k << std::endl;

      mesh->refine_towards_boundary("3");

      for(int i = 2; i < 5; i++)
      {
        std::cout << "i == " << i << std::endl;
        Hermes::Hermes2D::MeshSharedPtr egg_shell_mesh = Hermes::Hermes2D::EggShell::get_egg_shell(mesh, "4", i);

#ifdef _SHOW_OUTPUT
        m.show(egg_shell_mesh);
          View::wait(HERMES_WAIT_KEYPRESS);
#endif

        for(int j = 2; j < 4; j++)
        {
          std::cout << "j == " << j << std::endl;
          Hermes::Hermes2D::MeshFunctionSharedPtr<double> eggShell(new Hermes::Hermes2D::ExactSolutionEggShell(egg_shell_mesh, j));

#ifdef _SHOW_OUTPUT
          s.show(eggShell);
          View::wait(HERMES_WAIT_KEYPRESS);
#endif
        }
      }
    }
  }

  handle_success(true);
  return 0;
}
