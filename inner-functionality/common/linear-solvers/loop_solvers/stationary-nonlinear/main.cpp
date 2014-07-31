#include "definitions.h"
#include "../../../../../testing-core/testing-core.h"

// Initial polynomial degree.
const int P_INIT = 1;                             
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-4;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 2;                  
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 3;                   

// Problem parameters.
double heat_src = 1.0;
double alpha = 7.0;

int main(int argc, char* argv[])
{
#ifdef WITH_PARALUTION
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_GLOB_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Bdy", INIT_BDY_REF_NUM);

  // Initialize boundary conditions.
  CustomEssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> src(-heat_src);
  WeakFormSharedPtr<double> wf(new DefaultWeakFormPoisson<double>(HERMES_ANY, &lambda, &src));
#ifdef SHOW_OUTPUT
  ScalarView s_view("Solution");
#endif

  double* coeff_vec = new double[ndof];
  MeshFunctionSharedPtr<double> init_sln(new CustomInitialCondition(mesh));
  OGProjection<double> ogProjection; ogProjection.project_global(space, init_sln, coeff_vec); 
  MeshFunctionSharedPtr<double> sln(new Solution<double>);

  // Testing is happening here.
  std::vector<int> numbers_of_nonlinear_iterations;
  std::vector<int> numbers_of_linear_iterations_in_last_nonlinear_step;
  // Iterative - original initial guess.
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_ITERATIVE);
  {
    // Initialize Newton solver.
    NewtonSolver<double> newton(wf, space);
    newton.set_tolerance(NEWTON_TOL, ResidualNormAbsolute);
    newton.set_max_steps_with_reused_jacobian(0);

    newton.get_linear_matrix_solver()->as_IterSolver()->set_solver_type(Solvers::GMRES);
    newton.get_linear_matrix_solver()->as_IterSolver()->set_tolerance(1e-5);

    // Perform Newton's iteration.
    newton.solve(coeff_vec);
    numbers_of_nonlinear_iterations.push_back(newton.get_num_iters());
    numbers_of_linear_iterations_in_last_nonlinear_step.push_back(newton.get_linear_matrix_solver()->as_IterSolver()->get_num_iters());

    // Translate the resulting coefficient vector into a Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);
#ifdef SHOW_OUTPUT
    s_view.show(sln);
#endif
  }

  // Iterative - "exact" initial guess.
  {
    // Initialize Newton solver.
    NewtonSolver<double> newton(wf, space);
    newton.set_tolerance(NEWTON_TOL, ResidualNormAbsolute);

    newton.get_linear_matrix_solver()->as_IterSolver()->set_solver_type(Solvers::GMRES);

    // Perform Newton's iteration.
    newton.solve(sln);
    numbers_of_nonlinear_iterations.push_back(newton.get_num_iters());
    numbers_of_linear_iterations_in_last_nonlinear_step.push_back(newton.get_linear_matrix_solver()->as_IterSolver()->get_num_iters());

    // Translate the resulting coefficient vector into a Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);
#ifdef SHOW_OUTPUT
    s_view.show(sln);
#endif
  }

  bool success = true;
  
  success = Hermes::Testing::test_value(numbers_of_nonlinear_iterations[0], 10, "Nonlinear iterations[0]", 1) && success;
  success = Hermes::Testing::test_value(numbers_of_nonlinear_iterations[1], 1, "Nonlinear iterations[1]", 1) && success;
  success = Hermes::Testing::test_value(numbers_of_linear_iterations_in_last_nonlinear_step[0], 5, "Linear iterations[0]", 1) && success;
  success = Hermes::Testing::test_value(numbers_of_linear_iterations_in_last_nonlinear_step[1], 7, "Linear iterations[1]", 1) && success;
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