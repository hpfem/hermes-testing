#include "definitions.h"
#include "../../../../testing-core/testing-core.h"

//  This is a test of adaptivity.
//
//  PDE: -Laplace u - f = 0.
//
//  Known exact solution, see class CustomExactSolution in definitions.cpp.
//
//  Domain: square domain (0, pi) x (0, pi), mesh file square_quad.mesh->
//
//  BC:  Dirichlet, given by exact solution.
//
//  The following parameters can be changed:

int P_INIT = 1;                                   // Initial polynomial degree of all mesh elements.
const double THRESHOLD = 0.3;
const CandList CAND_LIST = H2D_HP_ANISO;          // Predefined list of element refinement candidates.
const double ERR_STOP = 1e-4;                     // Stopping criterion for adaptivity (rel. error tolerance between the
                                                  // reference mesh and coarse mesh solution in percent).
// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionCumulative<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("square_quad.mesh", mesh);

  // Avoid zero ndof situation.
  if(P_INIT == 1) {
    if(is_hp(CAND_LIST)) P_INIT++;
    else mesh->refine_element_id(0, 0);
  }

  // Define exact solution.
  MeshFunctionSharedPtr<double> exact_sln(new CustomExactSolution(mesh));

  // Initialize the weak formulation.
  Hermes1DFunction<double> lambda(1.0);
  CustomFunction f;
  DefaultWeakFormPoisson<double> wf(Hermes::HERMES_ANY, &lambda, &f);

  // Initialize boundary conditions.
  DefaultEssentialBCConst<double> bc_essential("Bdy", 0.0);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  // Initialize approximate solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  // Assemble the discrete problem.
  DiscreteProblem<double> dp(&wf, space);

  // Adaptivity loop:
  int as = 1; bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

    int ndof_ref = ref_space->get_num_dofs();

    dp.set_space(ref_space);

    // Initial coefficient vector for the Newton's method.
    double* coeff_vec = new double[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(double));

    NewtonSolver<double> newton(&dp);
    newton.set_verbose_output(false);

    MeshFunctionSharedPtr<double> ref_sln(new Solution<double>());

    try{
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }
    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    adaptivity.set_space(space);

    // If err_est too large, adapt the mesh-> The NDOF test must be here, so that the solution may be visualized
    // after ending due to this criterion.
    if(err_est_rel_total < ERR_STOP)
      done = true;
    else
      done = adaptivity.adapt(&selector);

    // Increase the counter of adaptivity steps.
    if(done == false)
      as++;

    // Clean up.
    delete [] coeff_vec;
  }
  while (done == false);

  bool success = Hermes::Testing::test_value(space->get_num_dofs(), 225, "DOFs", 1);

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
}