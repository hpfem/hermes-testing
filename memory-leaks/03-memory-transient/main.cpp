#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

// This example solves a time-domain resonator problem for the Maxwell's equation. 
// It is very similar to resonator-time-domain-I but B is eliminated from the 
// equations, thus converting the first-order system into one second -order
// equation in time. The second-order equation in time is decomposed back into 
// a first-order system in time in the standard way (see example wave-1). Time 
// discretization is performed using the implicit Euler method.
//
// PDE: \frac{1}{SPEED_OF_LIGHT**2}\frac{\partial^2 E}{\partial t^2} + curl curl E = 0,
// converted into
//
//      \frac{\partial E}{\partial t} - F = 0,
//      \frac{\partial F}{\partial t} + SPEED_OF_LIGHT**2 * curl curl E = 0.
//
// Approximated by
// 
//      \frac{E^{n+1} - E^{n}}{tau} - F^{n+1} = 0,
//      \frac{F^{n+1} - F^{n}}{tau} + SPEED_OF_LIGHT**2 * curl curl E^{n+1} = 0.
//
// Domain: Square (-pi/2, pi/2) x (-pi/2, pi/2)... See mesh file domain.mesh.
//
// BC:  E \times \nu = 0 on the boundary (perfect conductor),
//      F \times \nu = 0 on the boundary (E \times \nu = 0 => \partial E / \partial t \times \nu = 0).
//
// IC:  Prescribed wave for E, zero for F.
//
// The following parameters can be changed:

// Initial polynomial degree of mesh elements.
const int P_INIT = 3;                              
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 3;
// Time step.
const double time_step = 0.05;                     
// Final time.
const double T_FINAL = 0.3;                       
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-8;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Matrix solver: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
MatrixSolverType matrix_solver = SOLVER_UMFPACK;   

// Problem parameters.
// Square of wave speed.  
const double C_SQUARED = 1;                                         

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize solutions.
  CustomInitialConditionWave E_sln(&mesh);
  ZeroSolutionVector<double> F_sln(&mesh);
  Hermes::vector<Solution<double>*> slns(&E_sln, &F_sln);

  // Initialize the weak formulation.
  CustomWeakFormWaveIE wf(time_step, C_SQUARED, &E_sln, &F_sln);
  
  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_essential("Perfect conductor", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create x- and y- displacement space using the default H1 shapeset.
  HcurlSpace<double> E_space(&mesh, &bcs, P_INIT);
  HcurlSpace<double> F_space(&mesh, &bcs, P_INIT);
  Hermes::vector<const Space<double> *> spaces = Hermes::vector<const Space<double> *>(&E_space, &F_space);
  int ndof = HcurlSpace<double>::get_num_dofs(spaces);
  Hermes::Mixins::Loggable::Static::info("ndof = %d.", ndof);

  // Initialize the FE problem.
  DiscreteProblem<double>* dp = new DiscreteProblem<double>(&wf, spaces);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  // NOTE: If you want to start from the zero vector, just define 
  // coeff_vec to be a vector of ndof zeros (no projection is needed).
  double* coeff_vec = new double[ndof];
  OGProjection<double> ogProjection; ogProjection.project_global(spaces, slns, coeff_vec); 

  // Initialize Newton solver.
  NewtonSolver<double> newton(dp);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Perform Newton's iteration.
    try
    {
      newton.set_newton_max_iter(NEWTON_MAX_ITER);
      newton.set_newton_tol(NEWTON_TOL);
      newton.solve_keep_jacobian(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    }

    // Translate the resulting coefficient vector into Solutions.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), spaces, slns);

    // Update time.
    current_time += time_step;

  } while (current_time < T_FINAL);

  delete dp;

  return 0;
}
