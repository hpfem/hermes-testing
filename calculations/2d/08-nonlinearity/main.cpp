#include "definitions.h"

using namespace RefinementSelectors;

// Initial polynomial degree.
const int P_INIT = 2;                             
// Stopping criterion for the Newton's method.
const double NEWTON_TOL = 1e-8;                   
// Maximum allowed number of Newton iterations.
const int NEWTON_MAX_ITER = 100;                  
// Number of initial uniform mesh refinements.
const int INIT_GLOB_REF_NUM = 3;                  
// Number of initial refinements towards boundary.
const int INIT_BDY_REF_NUM = 4;                   

// Problem parameters.
double heat_src = 1.0;
double alpha = 7.0;

int main(int argc, char* argv[])
{
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
  DefaultWeakFormPoisson<double> wf(HERMES_ANY, &lambda, &src);

  // Initialize the FE problem.
  DiscreteProblem<double> dp(&wf, space);

  // Project the initial condition on the FE space to obtain initial 
  // coefficient vector for the Newton's method.
  // NOTE: If you want to start from the zero vector, just define 
  // coeff_vec to be a vector of ndof zeros (no projection is needed).
  double* coeff_vec = new double[ndof];
  MeshFunctionSharedPtr<double> init_sln(new CustomInitialCondition(mesh));
  OGProjection<double> ogProjection; ogProjection.project_global(space, init_sln, coeff_vec); 

  // Initialize Newton solver.
  NewtonSolver<double> newton(&dp);
  newton.set_tolerance(NEWTON_TOL, ResidualNormAbsolute);
  newton.set_max_allowed_residual_norm(1e99);
  newton.set_max_allowed_iterations(NEWTON_MAX_ITER);
  
  // 1st - OK
  newton.set_sufficient_improvement_factor_jacobian(0.5);
  newton.set_max_steps_with_reused_jacobian(5);
  newton.set_initial_auto_damping_coeff(0.95);
  newton.set_sufficient_improvement_factor(1.1);
  newton.set_min_allowed_damping_coeff(1e-10);

  // 2nd - OK
  // newton.set_sufficient_improvement_factor_jacobian(0.9);
  // newton.set_max_steps_with_reused_jacobian(10);
  // newton.set_manual_damping_coeff(0.1);
    
  // Perform Newton's iteration.
  try
  {
    newton.solve(coeff_vec);
  }
  catch(NewtonSolver<double>::NewtonException& e)
  {
    switch(e.get_exception_state())
    {
       case NewtonSolver<double>::AboveMaxIterations:
         std::cout << std::endl << "\t\t\tAboveMaxIterations" << std::endl;
         break;
      case NewtonSolver<double>::BelowMinDampingCoeff:
        std::cout << std::endl << "\t\t\tBelowMinDampingCoeff" << std::endl;
         break;
      case NewtonSolver<double>::AboveMaxAllowedResidualNorm:
        std::cout << std::endl << "\t\t\tAboveMaxAllowedResidualNorm" << std::endl;
         break;
    }
  }
  
  catch(std::exception& e)
  {
    std::cout << e.what();
    printf("Failure!\n");
    return -1;
  }

  // Translate the resulting coefficient vector into a Solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>);
  Solution<double>::vector_to_solution(newton.get_sln_vector(), space, sln);

  // Clean up.
  delete [] coeff_vec;

  printf("Success!\n");
  return 0;
}

