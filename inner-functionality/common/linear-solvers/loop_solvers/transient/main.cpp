#include "definitions.h"
#include "../../../../../testing-core/testing-core.h"

//  This example shows the simplest way to solve a linear time-dependent
//  PDE in Hermes using the implicit Euler method in time. The model describes 
//  in a naive approximation how the St. Vitus Cathedral in Prague 
//  (http://en.wikipedia.org/wiki/St._Vitus_Cathedral) responds to changes 
//  in the surrounding air temperature during one 24-hour cycle. 
//
//  PDE: non-stationary heat transfer equation
//  dT/dt - LAMBDA / (HEATCAP * RHO) * Laplace T = 0.
//
//  Domain: St. Vitus cathedral (file cathedral.mesh).
//
//  IC:  T = TEMP_INIT.
//  BC:  T = TEMP_INIT on the bottom edge ... Dirichlet,
//       LAMBDA * dT/dn = ALPHA*(t_exterior(time) - T) ... Newton, time-dependent.
//
//  Time-stepping: implicit Euler method.
//
//  The following parameters can be changed:

// Polynomial degree of mesh elements.
const int P_INIT = 2;                             
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;                       
// Number of initial uniform mesh refinements towards the boundary.
const int INIT_REF_NUM_BDY = 3;                   
// Time step in seconds.
const double time_step = 300.0;                   
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;  

// Problem parameters.
// Temperature of the ground (also initial temperature).
const double TEMP_INIT = 10;       
// Heat flux coefficient for Newton's boundary condition.
const double ALPHA = 10;           
// Thermal conductivity of the material.
const double LAMBDA = 1e2;         
// Heat capacity.
const double HEATCAP = 1e2;        
// Material density.
const double RHO = 3000;           
// Length of time interval (24 hours) in seconds.
const double T_FINAL = 15 * time_step;      

int main(int argc, char* argv[])
{
#ifdef WITH_PARALUTION
  HermesCommonApi.set_integral_param_value(matrixSolverType, SOLVER_PARALUTION_AMG);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();
  mesh->refine_towards_boundary("Boundary air", INIT_REF_NUM_BDY);
  mesh->refine_towards_boundary("Boundary ground", INIT_REF_NUM_BDY);

  // Previous time level solution (initialized by the external temperature).
  MeshFunctionSharedPtr<double> tsln(new ConstantSolution<double> (mesh, TEMP_INIT));

  // Initialize the weak formulation.
  double current_time = 0;
  WeakFormSharedPtr<double> wf(new CustomWeakFormHeatRK1("Boundary air", ALPHA, LAMBDA, HEATCAP, RHO, time_step, &current_time, TEMP_INIT, T_FINAL, tsln));
  
  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Boundary ground", TEMP_INIT);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);
 
  // Initialize Newton solver.
  NewtonSolver<double> newton(wf, space);
#ifdef SHOW_OUTPUT
  newton.set_verbose_output(true);
#else
  newton.set_verbose_output(false);
#endif
  newton.set_jacobian_constant();
  newton.get_linear_matrix_solver()->as_AMGSolver()->set_smoother(Solvers::GMRES, Preconditioners::ILU);
  newton.get_linear_matrix_solver()->as_LoopSolver()->set_tolerance(1e-1, RelativeTolerance);

#ifdef SHOW_OUTPUT
  // Initialize views.
  ScalarView Tview("Temperature", new WinGeom(0, 0, 450, 600));
  Tview.set_min_max_range(0,20);
  Tview.fix_scale_width(30);
#endif

  // Time stepping:
  int ts = 1;
  do 
  {
    Hermes::Mixins::Loggable::Static::info("---- Time step %d, time %3.5f s", ts, current_time);

    newton.solve();

    // Translate the resulting coefficient vector into the Solution sln.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), space, tsln);

#ifdef SHOW_OUTPUT
    // Visualize the solution.
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    Tview.set_title(title);
    Tview.show(tsln);
#endif

    // Increase current time and time step counter.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  // Wait for the view to be closed.
#ifdef SHOW_OUTPUT
  View::wait();
#endif
  return 0;
#endif
  return 0;
}
