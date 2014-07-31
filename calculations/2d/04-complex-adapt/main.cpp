#include "definitions.h"
#include "../../../testing-core/testing-core.h"

//  The following parameters can be changed:
const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;                     // This is a quantitative parameter of Adaptivity.

// Error calculation & adaptivity.
DefaultErrorCalculator<::complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionCumulative<::complex> stoppingCriterion(0.3);
// Adaptivity processor class.
Adapt<::complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e1;

// Problem parameters.
const double MU_0 = 4.0*M_PI*1e-7;
const double MU_IRON = 1e3 * MU_0;
const double GAMMA_IRON = 6e6;
const double J_EXT = 1e6;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<::complex> bc_essential("Dirichlet", ::complex (0.0, 0.0));
  EssentialBCs<::complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<::complex> space(new H1Space<::complex>(mesh, &bcs, P_INIT));

  // Initialize the weak formulation.
  WeakFormSharedPtr<::complex> wf(new CustomWeakForm("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON,
    "Wire", MU_0, ::complex (J_EXT, 0.0), OMEGA));

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<::complex> sln(new Hermes::Hermes2D::Solution<::complex>());
  MeshFunctionSharedPtr<::complex> ref_sln(new Hermes::Hermes2D::Solution<::complex>());

  // Initialize refinement selector.
  H1ProjBasedSelector<::complex> selector(CAND_LIST);

  // DOF and CPU convergence graphs initialization.
  SimpleGraph graph_dof, graph_cpu;

  DiscreteProblem<::complex> dp(wf, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<::complex> newton(&dp);
  
  // Adaptivity loop:
  int as = 1; bool done = false;
  adaptivity.set_space(space);
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<::complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<::complex> ref_space = ref_space_creator.create_ref_space();
    
    newton.set_space(ref_space);

    int ndof_ref = ref_space->get_num_dofs();

    // Initialize reference problem.

    // Initial coefficient vector for the Newton's method.
    ::complex * coeff_vec = new ::complex [ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(::complex));

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }

    Hermes::Hermes2D::Solution<::complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<::complex> ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);

    // If err_est too large, adapt the mesh->
    if(errorCalculator.get_total_error_squared()  * 100. < ERR_STOP)
      done = true;
    else
    {
      adaptivity.adapt(&selector);
    }

    // Clean up.
    delete [] coeff_vec;

    // Increase counter.
    as++;
  }
  while (done == false);

  ::complex  sum = 0;
  for (int i = 0; i < space->get_num_dofs(); i++)
    sum += newton.get_sln_vector()[i];
  printf("coefficient sum = %f\n", sum);

  ::complex  expected_sum;
  expected_sum.real(1.4685364e-005);
  expected_sum.imag(-5.45632171e-007);

  bool success = true;
  if(std::abs(sum - expected_sum) > 1e-6)
    success = false;

  int ndof = space->get_num_dofs();
  if(ndof != 82) // Tested value as of May 2013.
    success = false;

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