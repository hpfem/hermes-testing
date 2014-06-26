#include "definitions.h"
#include "../../../testing-core/testing-core.h"

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 3;                       
// Initial polynomial degree of all mesh elements.
const int P_INIT = 1;                             
// Time step. 
double time_step = 0.05;                           
// Time interval length.
const double T_FINAL = 0.3;                       

// Error calculation & adaptivity.
// Every UNREF_FREQth time step the mesh is derefined.
const int UNREF_FREQ = 2;
// 1... mesh reset to basemesh and poly degrees to P_INIT.   
// 2... one ref. layer shaved off, poly degrees reset to P_INIT.
// 3... one ref. layer shaved off, poly degrees decreased by one. 
const int UNREF_METHOD = 3;                       
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 5e-2;

// Newton's method
// Stopping criterion for Newton on fine mesh->
const double NEWTON_TOL = 1e-5;                   
// Maximum allowed number of Newton iterations.
const int max_allowed_iterations = 20;

ButcherTableType butcher_table_type = Implicit_RK_1;

// Problem parameters.
// Parameter for nonlinear thermal conductivity.
const double alpha = 4.0;                         
const double heat_src = 1.0;

int main(int argc, char* argv[])
{
  // Set just one thread.
  Hermes::HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Choose a Butcher's table or define your own.
  ButcherTable bt(butcher_table_type);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh), basemesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square.mesh", basemesh);

  // Perform initial mesh refinements.
  for(int i = 0; i < INIT_REF_NUM; i++) 
    basemesh->refine_all_elements(0, true);
  mesh->copy(basemesh);

  // Initialize boundary conditions.
  EssentialBCNonConst bc_essential("Bdy");
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));
  adaptivity.set_space(space);
  int ndof_coarse = space->get_num_dofs();

  // Previous time level solution (initialized by initial condition).
  MeshFunctionSharedPtr<double> sln_time_prev(new CustomInitialCondition (mesh));

  // Initialize the weak formulation
  CustomNonlinearity lambda(alpha);
  Hermes2DFunction<double> f(heat_src);
  WeakFormSharedPtr<double> wf(new CustomWeakFormPoisson(&lambda, &f));

  // Next time level solution.
  MeshFunctionSharedPtr<double> sln_time_new(new Solution<double>(mesh));

  // Create a refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Initialize Runge-Kutta time stepping.
  RungeKutta<double> runge_kutta(wf, space, &bt);

  // Time stepping loop.
  double current_time = 0; int ts = 1;
  do
  {
    // Periodic global derefinement.
    if (ts > 1 && ts % UNREF_FREQ == 0)
    {
      switch (UNREF_METHOD) 
      {
      case 1: mesh->copy(basemesh);
        space->set_uniform_order(P_INIT);
        break;
      case 2: mesh->unrefine_all_elements();
        space->set_uniform_order(P_INIT);
        break;
      case 3: mesh->unrefine_all_elements();
        space->adjust_element_order(-1, -1, P_INIT, P_INIT);
        break;
      }

      space->assign_dofs();
      ndof_coarse = Space<double>::get_num_dofs(space);
    }

    // Spatial adaptivity loop. Note: sln_time_prev must not be changed
    // during spatial adaptivity.
    bool done = false; int as = 1;
    do
    {
      // Construct globally refined reference mesh and setup reference space->
      Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
      MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
      int ndof_ref = Space<double>::get_num_dofs(ref_space);

      // Perform one Runge-Kutta time step according to the selected Butcher's table.
      try
      {
        runge_kutta.set_space(ref_space);
        runge_kutta.set_verbose_output(true);
        runge_kutta.set_time(current_time);
        runge_kutta.set_time_step(time_step);
        runge_kutta.set_tolerance(NEWTON_TOL);
        runge_kutta.rk_time_step_newton(sln_time_prev, sln_time_new);
      }
      catch(Exceptions::Exception& e)
      {
        std::cout << e.what();
      }

      // Project the fine mesh solution onto the coarse mesh.
      MeshFunctionSharedPtr<double> sln_coarse(new Solution<double>());
      OGProjection<double> ogProjection; ogProjection.project_global(space, sln_time_new, sln_coarse);

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(sln_coarse, sln_time_new);
      double err_est_rel = errorCalculator.get_total_error_squared() * 100;

      // If err_est too large, adapt the mesh->
      if (err_est_rel < ERR_STOP) done = true;
      else
      {
        done = adaptivity.adapt(&selector);
        as++;
      }
    }
    while (done == false);

    sln_time_prev->copy(sln_time_new);

    // Increase current time and counter of time steps.
    current_time += time_step;
    ts++;
  }
  while (current_time < T_FINAL);

  double* result = new double[runge_kutta.get_space(0)->get_num_dofs()];
  OGProjection<double>::project_global(runge_kutta.get_space(0), sln_time_prev, result);
  double sum = 0;
  for (int i = 0; i < space->get_num_dofs(); i++)
    sum += result[i];

  bool success = Testing::test_value(sum, 349.579072949291, "coefficient sum", 2.);

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

  return 0;
}
