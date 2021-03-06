#include "definitions.h"
#include "../../../../testing-core/testing-core.h"

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
  std::vector<MeshSharedPtr > meshes ({meshVertical, meshHorizontal});
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
  
  WeakFormSharedPtr<double> wf(new CustomWeakForm("Top Left"));

  MeshFunctionSharedPtr<double> slnTemp(new ZeroSolution<double>(meshVertical));
  MeshFunctionSharedPtr<double> slnAdv(new Solution<double>(meshHorizontal));

  wf->set_ext(slnTemp);

  std::vector<SpaceSharedPtr<double> > spaces({spaceTemperature, spaceAdvectionDiffusion});
  LinearSolver<double> solver(wf, spaces);
  for(int step = 0; step <= 1; step++)
  {
    try
    {
      solver.solve();
      Solution<double>::vector_to_solutions(solver.get_sln_vector(), { spaceTemperature, spaceAdvectionDiffusion }, {slnTemp, slnAdv});
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

  // Values as of 11112013
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

  bool success = true;

  success = Testing::test_value(test_point_1, testValues[0], "test_point_1") && success;
  success = Testing::test_value(test_point_2, testValues[1], "test_point_2") && success;
  success = Testing::test_value(test_point_3, testValues[2], "test_point_3") && success;

  success = Testing::test_value(test_point_4, testValues[3], "test_point_4") && success;
  success = Testing::test_value(test_point_5, testValues[4], "test_point_5") && success;
  success = Testing::test_value(test_point_6, testValues[5], "test_point_6") && success;
  success = Testing::test_value(test_point_7, testValues[6], "test_point_7") && success;
  success = Testing::test_value(test_point_8, testValues[7], "test_point_8") && success;

#ifdef WITH_BSON
  int first_mesh_elem_count = meshes[0]->get_num_active_elements();
  int second_mesh_elem_count = meshes[1]->get_num_active_elements();
  meshes[0]->refine_all_elements();
  meshes[1]->refine_all_elements(1);
  meshes[1]->refine_all_elements(2);

  Hermes2D::MeshReaderH2DBSON bson_loader;
  bson_loader.save("meshes.bson", meshes);
  MeshSharedPtr meshVertical2(new Mesh), meshHorizontal2(new Mesh);
  std::vector<MeshSharedPtr > meshes2({meshVertical2, meshHorizontal2});
  bson_loader.load("meshes.bson", meshes2);
  success = (4 * first_mesh_elem_count == meshes2[0]->get_num_active_elements()) && success;
  success = (4 * second_mesh_elem_count == meshes2[1]->get_num_active_elements()) && success;
#endif

  if(success)
  {
    std::cout << "Success!" << std::endl;
    return 0;
  }
  else
  {
    std::cout << "Failure!" << std::endl;
    return -1;
  }
}