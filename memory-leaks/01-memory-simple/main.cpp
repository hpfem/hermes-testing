#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
int P_INIT = 2;                             // Uniform polynomial degree of mesh elements.
int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  if(argc > 1)
    {
      if(argc > 2)
      {
        P_INIT = atoi(argv[1]);
        INIT_REF_NUM = atoi(argv[2]);
      }
      else
      {
        std::cout << (std::string)"Wrong number of parameters.";
        return -1;
      }
    }

  Hermes::Hermes2D::H1Space<double>* space = NULL;
  Hermes::Hermes2D::Mesh* mesh = new Mesh();

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

wf.set_verbose_output(false);
  
	Hermes2DApi.set_text_param_value(xmlSchemasDirPath, "asfd");

  // This is in a block to test that the instances mesh and space can be deleted after being copied with no harm.
  {
    // Load the mesh.
    Hermes::Hermes2D::MeshReaderH2DXML mloader;
    mloader.set_validation(false);
    mloader.load("domain.xml", mesh);

    // Perform initial mesh refinements (optional).
    mesh->refine_in_areas(Hermes::vector<std::string>("Aluminum", "Copper"), INIT_REF_NUM);
    mesh->refine_in_area("Aluminum");

    // Create an H1 space with default shapeset.
    space = new Hermes::Hermes2D::H1Space<double>(mesh, &bcs, P_INIT);
  }

  Mesh* new_mesh = new Mesh();
  H1Space<double>* new_space = new H1Space<double>();
  new_space->copy(space, new_mesh);

  delete space;
  delete mesh;

  Hermes::Hermes2D::Element* e;
  int i = 1;
  for_all_active_elements(e, new_mesh)
  {
    new_space->set_element_order(e->id, i++ % 4 + 1);
  }

  new_space->assign_dofs();

  // Initialize the solution.
  Hermes::Hermes2D::Solution<double>* sln = new Hermes::Hermes2D::Solution<double>();

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver(&wf, new_space);
  linear_solver.set_verbose_output(false);

  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Hermes::Hermes2D::Solution<double>::vector_to_solution(sln_vector, new_space, sln);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  Views::Linearizer* lin = new Views::Linearizer();
  lin->process_solution(sln);
  lin->free();
  lin->process_solution(sln);
  lin->free();

  delete sln;
  delete new_space;
  delete new_mesh;

  std::cout << "OK";
  return 0;
}
