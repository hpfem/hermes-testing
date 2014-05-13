#include "definitions.h"
#include "../../../../../testing-core/testing-core.h"

// This example shows how to run adaptive hp-FEM, h-FEM and p-FEM with
// basic control parameters. The underlying problem is a planar model
// of an electrostatic micromotor (MEMS). You may want to experiment with
// various types of adaptivity via the options H2D_P_ISO, H2D_P_ANISO,
// H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H, H2D_HP_ANISO_P,
// and H2D_HP_ANISO. See the User Documentation for more details.
//
// Uniform initial polynomial degree of mesh elements can be set using
// the variable P_INIT. Before using adaptivity, you have to define a refinement
// selector as shown below. The function adapt() takes the selector
// as a parameter, along with THRESHOLD, STRATEGY, and MESH_REGULARITY.
//   
// Additional control parameters are available, these will be demonstrated
// in the following tutorial examples. In this example, two types of convergence
// graphs are created -- error estimate wrt. the number of degrees of freedom
// (DOF), and error estimate wrt. CPU time. Later we will show how to output
// the error wrt. exact solution when exact solution is available.
//   
// This example also demonstrates how to define different material parameters
// in various parts of the computational domain, and how to measure time.
//
// PDE: -div[eps_r(x,y) grad phi] = 0
//      eps_r = EPS_1 in Omega_1 (surrounding air)
//      eps_r = EPS_2 in Omega_2 (moving part of the motor)
//
// BC: phi = 0 V on Gamma_1 (left edge and also the rest of the outer boundary
//     phi = VOLTAGE on Gamma_2 (boundary of stator)
//
// The following parameters can be changed:

// Set to "false" to suppress Hermes OpenGL visualization. 
const bool HERMES_VISUALIZATION = true;           
// Set to "true" to enable VTK output.
const bool VTK_VISUALIZATION = false;             
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.8;                     
// Predefined list of element refinement candidates. Possible values are
// H2D_P_ISO, H2D_P_ANISO, H2D_H_ISO, H2D_H_ANISO, H2D_HP_ISO, H2D_HP_ANISO_H
// H2D_HP_ANISO_P, H2D_HP_ANISO.
const CandList CAND_LIST = H2D_HP_ANISO_H;        
// Stopping criterion for adaptivity.
const double ERR_STOP = 10.0;
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);

// Problem parameters.
const double EPS0 = 8.863e-12;
const double VOLTAGE = 50.0;
const double EPS_MOTOR = 10.0 * EPS0;
const double EPS_AIR = 1.0 * EPS0;

int main(int argc, char* argv[])
{
#ifdef WITH_PARALUTION
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_AMG);
  HermesCommonApi.set_integral_param_value(numThreads,1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakFormPoisson("Motor", EPS_MOTOR, "Air", EPS_AIR));

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential_out("Outer", 0.0);
  DefaultEssentialBCConst<double> bc_essential_stator("Stator", VOLTAGE);
  EssentialBCs<double> bcs(std::vector<EssentialBoundaryCondition<double> *>(&bc_essential_out, &bc_essential_stator));

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  adaptivity.set_space(space);

  // Initialize coarse and fine mesh solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>()), ref_sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize views.
#ifdef SHOW_OUTPUT
  Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 410, 600));
  sview.fix_scale_width(50);
  sview.show_mesh(false);
  Views::OrderView  oview("Polynomial orders", new Views::WinGeom(420, 0, 400, 600));
#endif

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  // Time measurement.
  Hermes::Mixins::TimeMeasurable cpu_time;

  // solver
  LinearSolver<double> solver(wf, space);
  solver.get_linear_matrix_solver()->as_AMGSolver()->set_tolerance(1e-1, RelativeTolerance);

  // Test values.
  int linear_iterations = 0;

  // Adaptivity loop:
  int as = 1;
  bool done = false;
  do
  {
    Hermes::Mixins::Loggable::Static::info("---- Adaptivity step %d:", as);

    // Time measurement.
    cpu_time.tick();

    // Ref. space creator.
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    // Construct globally refined mesh and setup fine mesh space.
    Space<double>::ReferenceSpaceCreator ref_space_creator;
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space(space, ref_mesh, true);

    // Initial vector.
    double* coeff_vec = new double[ref_space->get_num_dofs()];
    if(as > 1)
      OGProjection<double>::project_global(ref_space, ref_sln, coeff_vec);
    else
      memset(coeff_vec, 0, sizeof(double) * ref_space->get_num_dofs());

    // Initialize fine mesh problem.
    Hermes::Mixins::Loggable::Static::info("Solving on fine mesh.");


    // Perform solver's iteration.
    try
    {
      solver.set_space(ref_space);
      solver.solve(coeff_vec);
      linear_iterations += solver.get_linear_matrix_solver()->as_AMGSolver()->get_num_iters();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(solver.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    Hermes::Mixins::Loggable::Static::info("Projecting fine mesh solution on coarse mesh.");
    OGProjection<double>::project_global(space, ref_sln, sln);

    // Cleanup.
    delete [] coeff_vec;

    // Time measurement.
    cpu_time.tick();

#ifdef SHOW_OUTPUT
    // VTK output.
    if (VTK_VISUALIZATION) 
    {
      // Output solution in VTK format.
      Views::Linearizer lin(FileExport);
      char* title = new char[100];
      sprintf(title, "sln-%d.vtk", as);
      lin.save_solution_vtk(sln, title, "Potential", false);
      Hermes::Mixins::Loggable::Static::info("Solution in VTK format saved to file %s.", title);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      sprintf(title, "ord-%d.vtk", as);
      ord.save_orders_vtk(space, title);
      Hermes::Mixins::Loggable::Static::info("Element orders in VTK format saved to file %s.", title);
    }

    // View the coarse mesh solution and polynomial orders.
    if (HERMES_VISUALIZATION) 
    {
      sview.show(sln);
      oview.show(space);
    }
#endif

    // Skip visualization time.
    cpu_time.tick();

    // Calculate element errors and total error estimate.
    Hermes::Mixins::Loggable::Static::info("Calculating error estimate.");

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);
    std::cout << "Relative error: " << errorCalculator.get_total_error_squared() * 100. << '%' << std::endl;

    // Add entry to DOF and CPU convergence graphs.
    cpu_time.tick();

    // If err_est too large, adapt the mesh.
    // If err_est too large, adapt the mesh->
    if(errorCalculator.get_total_error_squared()  * 100. < ERR_STOP)
      done = true;
    else
    {
      std::cout << "Adapting..." << std::endl;
      adaptivity.adapt(&selector);
    }
    as++;
  }
  while (done == false);

  Hermes::Mixins::Loggable::Static::info("Total running time: %g s", cpu_time.accumulated());
  Hermes::Mixins::Loggable::Static::info("Number of adaptivity steps: %d.", as);
  Hermes::Mixins::Loggable::Static::info("Number of linear iterations: %d.", linear_iterations);

  bool success = true;
  success = Hermes::Testing::test_value(as, 6, "Adaptivity steps", 1) && success;
  success = Hermes::Testing::test_value(linear_iterations, 22, "Linear iterations", 1) && success;

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
#endif
  return 0;
}