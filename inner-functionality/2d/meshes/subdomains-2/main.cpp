#define HERMES_REPORT_INFO
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Hermes2D::Views;

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 4;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr meshVertical(new Mesh), meshHorizontal(new Mesh);
  Hermes::vector<MeshSharedPtr> meshes (meshVertical, meshHorizontal);
  MeshReaderH2DXML mloader;
	mloader.set_validation(false);
  mloader.load("subdomains.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  DefaultEssentialBCConst<double> bc_essentialTemperature(HERMES_ANY, 0.0);
  EssentialBCs<double> bcsTemperature(&bc_essentialTemperature);
  DefaultEssentialBCConst<double> bc_essentialAdvectionDiffusion(HERMES_ANY, 0.0);
  EssentialBCs<double> bcsAdvectionDiffusion(&bc_essentialAdvectionDiffusion);

  SpaceSharedPtr<double> spaceTemperature(new H1Space<double>(meshVertical, &bcsTemperature, P_INIT));
  SpaceSharedPtr<double> spaceAdvectionDiffusion(new H1Space<double>(meshHorizontal, &bcsAdvectionDiffusion, P_INIT));
  
  CustomWeakForm wf("Top Left");

  MeshFunctionSharedPtr<double> slnTemp(new ZeroSolution<double>(meshVertical));
  MeshFunctionSharedPtr<double> slnAdv(new Solution<double>(meshHorizontal));

  wf.set_ext(slnTemp);

  ScalarView sTemp("Temperature");
  ScalarView sAdv("Advection");

  Hermes::vector<SpaceSharedPtr<double> > spaces(spaceTemperature, spaceAdvectionDiffusion);
  LinearSolver<double> solver(&wf, spaces);
  for(int step = 0; step <= 1; step++)
  {
    try
    {
      solver.solve();
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), Hermes::vector<SpaceSharedPtr<double> >(spaceTemperature, spaceAdvectionDiffusion), Hermes::vector<MeshFunctionSharedPtr<double> >(slnTemp, slnAdv));
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
      return -1;
    }
  }


  // Values as of 04022013
  /*
    9.7117986964428873
    -21.846810183868349
    -1.1517121406884645
    0.070238847349849406
    0.069616820552871647
    0.067679010586017496
    0.064375179813698405
    -0.25168293544547149
  */

  double testValues[8] = { 9.7117986964428873,
    -21.846810183868349,
    -1.1517121406884645,
    0.070238847349849406,
    0.069616820552871647,
    0.067679010586017496,
    0.064375179813698405,
    -0.25168293544547149
  };

  double test_point_1 = slnTemp->get_pt_value(0.5, 1.5)->val[0];
  double test_point_2 = slnTemp->get_pt_value(0.75, 1.23)->dx[0];
  double test_point_3 = slnTemp->get_pt_value(0.87, 1.24)->dy[0];

  double test_point_4 = slnAdv->get_pt_value(1.17, 1.8)->val[0];
  double test_point_5 = slnAdv->get_pt_value(1.27, 1.8)->val[0];
  double test_point_6 = slnAdv->get_pt_value(1.37, 1.8)->val[0];
  double test_point_7 = slnAdv->get_pt_value(1.47, 1.8)->val[0];
  double test_point_8 = slnAdv->get_pt_value(1.81, 1.67)->dx[0] + slnAdv->get_pt_value(1.81, 1.67)->dy[0];

  if(std::abs(test_point_1 - testValues[0]) > 1e-6)
    return -1;
  if(std::abs(test_point_2 - testValues[1]) > 1e-6)
    return -1;
  if(std::abs(test_point_3 - testValues[2]) > 1e-6)
    return -1;
  if(std::abs(test_point_4 - testValues[3]) > 1e-6)
    return -1;
  if(std::abs(test_point_5 - testValues[4]) > 1e-6)
    return -1;
  if(std::abs(test_point_6 - testValues[5]) > 1e-6)
    return -1;
  if(std::abs(test_point_7 - testValues[6]) > 1e-6)
    return -1;
  if(std::abs(test_point_8 - testValues[7]) > 1e-6)
    return -1;

  std::cout << "Success!" << std::endl;
  return 0;
}