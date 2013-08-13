#include "definitions.h"
#include "../../../../testing-core/testing-core.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;
using namespace Hermes::Algebra::DenseMatrixOperations;
using namespace Hermes::Solvers;

//#define SHOW_OUTPUT
const int P_INIT = 3;                     // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 3;               // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

class CustomExternalSolver : public SimpleExternalSolver<double>
{
public:
  CustomExternalSolver(CSCMatrix<double> *m, SimpleVector<double> *rhs) : SimpleExternalSolver<double>(m, rhs) {};

  virtual std::string command()
  {
    SimpleVector<double> vec(this->m->get_size());
    vec.alloc(this->m->get_size());
    double* raw = new double[this->m->get_size()];
    memset(raw, 0, sizeof(double) * this->m->get_size());
    raw[1] = 500.;
    raw[3] = 300.;
    raw[5] = 100.;
    vec.set_vector(raw);
    vec.export_to_file("temp-vector-export", "a", Algebra::EXPORT_FORMAT_PLAIN_ASCII);
    return "temp-vector-export";
  }
};

Hermes::Solvers::ExternalSolver<double>* my_fn(CSCMatrix<double> *m, SimpleVector<double> *rhs)
{
  return new CustomExternalSolver(m, rhs);
}

int main(int argc, char* argv[])
{
  HermesCommonApi.set_integral_param_value(matrixSolverType, Hermes::SOLVER_EXTERNAL);
  ExternalSolver<double>::create_external_solver = my_fn;

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  // Refine all elements, do it INIT_REF_NUM-times.
  for(unsigned int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"), 0.);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize space->
  SpaceSharedPtr<double> space( new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT));

  std::cout << "Ndofs: " << space->get_num_dofs() << std::endl;

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, space);

  // Solve the linear problem.
  try
  {
    linear_solver.solve();
    linear_solver.solve(linear_solver.get_sln_vector());

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, space, sln);

#ifdef SHOW_OUTPUT
    {
      // Visualize the solution.
      Hermes::Hermes2D::Views::ScalarView viewS("Solution", new Hermes::Hermes2D::Views::WinGeom(50, 50, 1000, 800));

      viewS.show(sln, Hermes::Hermes2D::Views::HERMES_EPS_LOW);

      viewS.wait_for_close();
    }
#endif
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  std::cout << sln->get_pt_value(.1,.1)->val[0] << std::endl;
  std::cout << sln->get_pt_value(.2,.1)->val[0] << std::endl;
  std::cout << sln->get_pt_value(.3,.1)->val[0] << std::endl;

  bool success = true;
  success = Testing::test_value(sln->get_pt_value(.1,.1)->val[0], 434.315, "get_pt_value(.1,.1)", 1e-3) && success;
  success = Testing::test_value(sln->get_pt_value(.2,.1)->val[0], 286.978, "get_pt_value(.2,.1)", 1e-3) && success;
  success = Testing::test_value(sln->get_pt_value(.3,.1)->val[0], 101.786, "get_pt_value(.3,.1)", 1e-3) && success;

  if(success)
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