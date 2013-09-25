#include "definitions.h"
#include "../../../testing-core/testing-core.h"

const int P_INIT = 2;                             // Polynomial degree of all mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.
const int INIT_REF_NUM_BDY = 3;                   // Number of initial uniform mesh refinements towards the boundary.
const double time_step = 1;                       // Time step in seconds.
const double NEWTON_TOL = 1e-5;                   // Stopping criterion for the Newton's method.
const int NEWTON_MAX_ITER = 100;                  // Maximum allowed number of Newton iterations.

ButcherTableType butcher_table_type = Implicit_SDIRK_2_2;

// Problem parameters.
const double TEMP_INIT = 10;       // Temperature of the ground (also initial temperature).
const double ALPHA = 10;           // Heat flux coefficient for Newton's boundary condition.
const double LAMBDA = 1e2;         // Thermal conductivity of the material.
const double HEATCAP = 1e2;        // Heat capacity.
const double RHO = 3000;           // Material density.
const double T_FINAL = 5*time_step;

int main(int argc, char* argv[])
{
  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("cathedral.mesh", mesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();
  mesh->refine_towards_boundary("Boundary_air", INIT_REF_NUM_BDY);
  mesh->refine_towards_boundary("Boundary_ground", INIT_REF_NUM_BDY);

  // Previous and next time level solutions.
  MeshFunctionSharedPtr<double> sln_time_prev(new ConstantSolution<double>(mesh, TEMP_INIT));
  MeshFunctionSharedPtr<double> sln_time_new(new Solution<double>(mesh));

  // Initialize the weak formulation.
	double current_time = 0;

  CustomWeakFormHeatRK wf("Boundary_air", ALPHA, LAMBDA, HEATCAP, RHO,
                          &current_time, TEMP_INIT, T_FINAL);
  wf.set_global_integration_order(10);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Boundary_ground", TEMP_INIT);
  EssentialBCs<double>bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  int ndof = space->get_num_dofs();

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(&wf, space, &bt);
  runge_kutta.set_verbose_output(true);

  // Iteration number.
  int iteration = 0;
    
  // Time stepping loop:
  do
  {
    // Perform one Runge-Kutta time step according to the selected Butcher's table.
    try
    {
      runge_kutta.set_space(space);
      runge_kutta.set_time(current_time);
      runge_kutta.set_time_step(time_step);
      runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
    }
    catch(Exceptions::Exception& e)
    {
      e.print_msg();
    }

    // Copy solution for the new time step.
    sln_time_prev->copy(sln_time_new);

    // Increase current time and time step counter.
    current_time += time_step;
  }
  while (current_time < T_FINAL);

  /* Begin test */

  bool success = true;

  if(fabs(sln_time_new->get_pt_value(-3.5, 17.0)->val[0] - 10.00271206) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(-1.0, 2.0)->val[0] - 10.0) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(0.0, 9.5)->val[0] - 10.00005812) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value( 1.0, 2.0)->val[0] - 10.0) > 1E-6) success = false;
  if(fabs(sln_time_new->get_pt_value(3.5, 17.0)->val[0] - 10.00271206) > 1E-6) success = false;

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
