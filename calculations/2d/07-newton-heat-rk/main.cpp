#include "definitions.h"
#include "../../../testing-core/testing-core.h"

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 1;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 3e3;                       // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;         // Thermal conductivity of the material.
const double HEATCAP = 1e2;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 86400;      // Length of time interval (24 hours) in seconds.

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Initialize the time.
  double current_time = 0;

  // mesh->
  MeshSharedPtr mesh(new Mesh);

  // Init mesh->
  MeshReaderH2D mloader;
  mloader.load("cathedral.mesh", mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  mesh->refine_towards_boundary("Boundary_air", INIT_REF_NUM_BDY);
  mesh->refine_towards_boundary("Boundary_ground", INIT_REF_NUM_BDY);

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential("Boundary_ground", TEMP_INIT);
  Hermes::Hermes2D::EssentialBCs<double> bcs(&bc_essential);

  // space->
  SpaceSharedPtr<double> space1(new H1Space<double>(mesh, &bcs, P_INIT));
  SpaceSharedPtr<double> space2(new H1Space<double>(mesh, &bcs, P_INIT + 1));
  std::vector<SpaceSharedPtr<double> > spaces({ space1, space2 });

  // Solution pointer.
  MeshFunctionSharedPtr<double> sln_time_prev1(new ConstantSolution<double>(mesh, TEMP_INIT));
  MeshFunctionSharedPtr<double> sln_time_prev2(new ConstantSolution<double>(mesh, TEMP_INIT));
  std::vector<MeshFunctionSharedPtr<double> > sln_time_prev({sln_time_prev1, sln_time_prev2});

  MeshFunctionSharedPtr<double> sln_time_new1(new Solution<double>(mesh));
  MeshFunctionSharedPtr<double> sln_time_new2(new Solution<double>(mesh));
  std::vector<MeshFunctionSharedPtr<double> > sln_time_new({sln_time_new1, sln_time_new2});

  CustomWeakFormHeatRK wf("Boundary_air", ALPHA, LAMBDA, HEATCAP, RHO,
    &current_time, TEMP_INIT, T_FINAL);

#ifdef _SHOW_OUTPUT
  // Initialize views.
  Hermes::Hermes2D::Views::ScalarView Tview("Temperature", new Hermes::Hermes2D::Views::WinGeom(0, 0, 450, 600));
  Tview.set_min_max_range(0, 20);
  Tview.fix_scale_width(30);
#endif

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, spaces, &bt);
  runge_kutta.set_tolerance(NEWTON_TOL);
  runge_kutta.set_verbose_output(true);
  runge_kutta.set_time_step(time_step);

  // Iteration number.
  int iteration = 0;

  // Time stepping loop:
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    try
    {
      runge_kutta.set_time(current_time);
      runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
    }
    catch(Exceptions::Exception& e)
    {
      e.print_msg();
    }

    // Show the new time level solution.
#ifdef _SHOW_OUTPUT
    char title[100];
    sprintf(title, "Time %3.2f s", current_time);
    Tview.set_title(title);
    Tview.show(sln_time_new2);
#endif

    // Copy solution for the new time step.
    sln_time_prev1->copy(sln_time_new1);
    sln_time_prev2->copy(sln_time_new2);

    // Increase current time and time step counter.
    current_time += time_step;
    iteration++;
  }
  while (iteration < 10);

  /* Begin test */

  bool success = true;

  success = Testing::test_value(sln_time_new1->get_pt_value(-3.5, 17.0)->val[0], 17.305684213090167, "sln_time_new1->get_pt_value(-3.5, 17.0)->val[0]", 1e-6) && success;
  success = Testing::test_value(sln_time_new2->get_pt_value(-4.9, 2.0)->val[0], 11.762112192614151, "sln_time_new1->get_pt_value(-3.5, 17.0)->val[0]", 1e-6) && success;
  success = Testing::test_value(sln_time_new2->get_pt_value(0.0, 9.95)->val[0], 15.216437017723038, "sln_time_new1->get_pt_value(-3.5, 17.0)->val[0]", 1e-6) && success;
  success = Testing::test_value(sln_time_new1->get_pt_value(4.9, 2.0)->val[0], 11.756504696911865, "sln_time_new1->get_pt_value(-3.5, 17.0)->val[0]", 1e-6) && success;
  success = Testing::test_value(sln_time_new1->get_pt_value(3.5, 17.0)->val[0], 17.305684213090171, "sln_time_new1->get_pt_value(-3.5, 17.0)->val[0]", 1e-6) && success;

#ifdef _SHOW_OUTPUT
  View::wait();
#endif

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
