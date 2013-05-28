#define HERMES_REPORT_ALL
#include "definitions.h"

const int P_INIT = 2;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.set_validation(false);
  mloader.load("domain.xml", mesh);

  // Perform initial mesh refinements (optional).
  mesh->refine_in_areas(Hermes::vector<std::string>("Aluminum", "Copper"), INIT_REF_NUM);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize linear solver.
  LinearSolver<double> linear_solver(&wf, space);

  // Solve the linear problem.
  linear_solver.solve();

  // Actual test. The values of 'sum' depend on the
  // current shapeset. If you change the shapeset,
  // you need to correct these numbers.
  double sum = 0;
  for (int i = 0; i < space->get_num_dofs(); i++)
    sum += linear_solver.get_sln_vector()[i];
  printf("coefficient sum = %f\n", sum);

  bool success = true;
  if(std::abs(sum - 3200.074172) > 1e-5)
    success = false;

  if(success == true)
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
