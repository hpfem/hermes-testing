#define HERMES_REPORT_WARN
#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

const bool STOKES = false;                        // For application of Stokes flow (creeping flow).

// If this is defined, the pressure is approximated using
// discontinuous L2 elements (making the velocity discreetely
// divergence-free, more accurate than using a continuous
// pressure approximation). Otherwise the standard continuous
// elements are used. The results are striking - check the
// tutorial for comparisons.
#define PRESSURE_IN_L2

const int P_INIT_VEL = 2;                         // Initial polynomial degree for velocity components.

// Initial polynomial degree for pressure.
// Note: P_INIT_VEL should always be greater than
// P_INIT_PRESSURE because of the inf-sup condition.
const int P_INIT_PRESSURE = 1;

const double RE = 200.0;                          // Reynolds number.
const double VEL_INLET = 1.0;                     // Inlet velocity (reached after STARTUP_TIME).

// During this time, inlet velocity increases gradually
// from 0 to VEL_INLET, then it stays constant.
const double STARTUP_TIME = 1.0;

const double TAU = 0.1;                           // Time step.
const double T_FINAL = 0.21;                      // Time interval length.
const double NEWTON_TOL = 1e-3;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 10;                   // Maximum allowed number of Newton iterations.

// Domain height (necessary to define the parabolic
// velocity profile at inlet).
const double H = 5;

// Boundary markers.
const std::string BDY_BOTTOM = "1";
const std::string BDY_RIGHT = "2";
const std::string BDY_TOP = "3";
const std::string BDY_LEFT = "4";
const std::string BDY_OBSTACLE = "5";

// Current time (used in weak forms).
double current_time = 0;

// Weak forms.
#include "definitions.cpp"

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Initial mesh refinements.
  //mesh->refine_all_elements();
  mesh->refine_towards_boundary(BDY_OBSTACLE, 4, false);
  mesh->refine_towards_boundary(BDY_TOP, 4, true);     // '4' is the number of levels,
  mesh->refine_towards_boundary(BDY_BOTTOM, 4, true);  // 'true' stands for anisotropic refinements.

  // Initialize boundary conditions.
  EssentialBCNonConst bc_left_vel_x(BDY_LEFT, VEL_INLET, H, STARTUP_TIME);
  DefaultEssentialBCConst<double> bc_other_vel_x(Hermes::vector<std::string>(BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_x(Hermes::vector<EssentialBoundaryCondition<double> *>(&bc_left_vel_x, &bc_other_vel_x));
  DefaultEssentialBCConst<double> bc_vel_y(Hermes::vector<std::string>(BDY_LEFT, BDY_BOTTOM, BDY_TOP, BDY_OBSTACLE), 0.0);
  EssentialBCs<double> bcs_vel_y(&bc_vel_y);
  EssentialBCs<double> bcs_pressure;

  // Spaces for velocity components and pressure.
  SpaceSharedPtr<double> xvel_space(new H1Space<double>(mesh, &bcs_vel_x, P_INIT_VEL));
  SpaceSharedPtr<double> yvel_space(new H1Space<double>(mesh, &bcs_vel_y, P_INIT_VEL));
#ifdef PRESSURE_IN_L2
  SpaceSharedPtr<double> p_space(new L2Space<double> (mesh, P_INIT_PRESSURE));
#else
  SpaceSharedPtr<double> p_space(new H1Space<double> (mesh, &bcs_pressure, P_INIT_PRESSURE));
#endif

  // Calculate and report the number of degrees of freedom.
  int ndof = Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space));

  // Define projection norms.
  ProjNormType vel_proj_norm = HERMES_H1_NORM;
#ifdef PRESSURE_IN_L2
  ProjNormType p_proj_norm = HERMES_L2_NORM;
#else
  ProjNormType p_proj_norm = HERMES_H1_NORM;
#endif

  // Solutions for the Newton's iteration and time stepping.
  MeshFunctionSharedPtr<double> xvel_prev_time(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> yvel_prev_time(new ConstantSolution<double> (mesh, 0.0));
  MeshFunctionSharedPtr<double> p_prev_time(new ConstantSolution<double> (mesh, 0.0));

  // Initialize weak formulation.
  WeakForm<double>* wf = new WeakFormNSNewton(STOKES, RE, TAU, xvel_prev_time, yvel_prev_time);

  wf->set_ext(Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_prev_time, yvel_prev_time));

	// Initialize the Newton solver.
	NewtonSolver<double> newton;
	newton.set_weak_formulation(wf);
	newton.set_spaces(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space));

  // Project the initial condition on the FE space to obtain initial
  // coefficient vector for the Newton's method.
  double* coeff_vec = new double[Space<double>::get_num_dofs(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space))];
  OGProjection<double> ogProjection;

  ogProjection.project_global(Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space),
    Hermes::vector<MeshFunctionSharedPtr<double> >(xvel_prev_time, yvel_prev_time, p_prev_time),
    coeff_vec, Hermes::vector<ProjNormType>(vel_proj_norm, vel_proj_norm, p_proj_norm));

  newton.set_newton_max_iter(NEWTON_MAX_ITER);
  newton.set_newton_tol(NEWTON_TOL);

  // Time-stepping loop:
  int num_time_steps = T_FINAL / TAU;
  for (int ts = 1; ts <= num_time_steps; ts++)
  {
    current_time += TAU;

    // Update time-dependent essential BCs.
    if(current_time <= STARTUP_TIME)
      newton.set_time(current_time);

    // Perform Newton's iteration and translate the resulting coefficient vector into previous time level solutions.
    try
    {
			newton.set_weak_formulation(wf);
      newton.solve_keep_jacobian(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }
    Hermes::vector<MeshFunctionSharedPtr<double> > tmp(xvel_prev_time, yvel_prev_time, p_prev_time);
    Hermes::Hermes2D::Solution<double>::vector_to_solutions(newton.get_sln_vector(), 
      Hermes::vector<SpaceSharedPtr<double> >(xvel_space, yvel_space, p_space), tmp);

  }
  delete [] coeff_vec;

  int success = 1;
  double eps = 1e-5;
  if(fabs(xvel_prev_time->get_pt_value(0.0, 2.5)->val[0] - 0.200000) > eps) {
    printf("Coordinate (   0, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(0.0, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(xvel_prev_time->get_pt_value(5, 2.5)->val[0] - 0.134291) > eps) {
    printf("Coordinate (   5, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(xvel_prev_time->get_pt_value(7.5, 2.5)->val[0] - 0.135088) > eps) {
    printf("Coordinate ( 7.5, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(7.5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(xvel_prev_time->get_pt_value(10, 2.5)->val[0] - 0.134944) > eps) {
    printf("Coordinate (  10, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(10, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(xvel_prev_time->get_pt_value(12.5, 2.5)->val[0] - 0.134888) > eps) {
    printf("Coordinate (12.5, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(12.5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(xvel_prev_time->get_pt_value(15, 2.5)->val[0] - 0.134864) > eps) {
    printf("Coordinate (  15, 2.5)->val[0] xvel value is %g\n", xvel_prev_time->get_pt_value(15, 2.5)->val[0]);
    success = 0;
  }

  if(fabs(yvel_prev_time->get_pt_value(0.0, 2.5)->val[0] - 0.000000) > eps) {
    printf("Coordinate (   0, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(0.0, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(yvel_prev_time->get_pt_value(5, 2.5)->val[0] - 0.000493) > eps) {
    printf("Coordinate (   5, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(yvel_prev_time->get_pt_value(7.5, 2.5)->val[0] - 0.000070) > eps) {
    printf("Coordinate ( 7.5, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(7.5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(yvel_prev_time->get_pt_value(10, 2.5)->val[0] - 0.000008) > eps) {
    printf("Coordinate (  10, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(10, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(yvel_prev_time->get_pt_value(12.5, 2.5)->val[0] + 0.000003) > eps) {
    printf("Coordinate (12.5, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(12.5, 2.5)->val[0]);
    success = 0;
  }
  if(fabs(yvel_prev_time->get_pt_value(15, 2.5)->val[0] + 0.000006) > eps) {
    printf("Coordinate (  15, 2.5)->val[0] yvel value is %g\n", yvel_prev_time->get_pt_value(15, 2.5)->val[0]);
    success = 0;
  }

  if(success == 1) {
    printf("Success!\n");
    return 0;
  }
  else {
    printf("Failure!\n");
    return -1;
  }
}
