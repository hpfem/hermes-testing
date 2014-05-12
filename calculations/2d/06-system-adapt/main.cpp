#include "definitions.h"
#include "../../../testing-core/testing-core.h"

// Initial polynomial degree for u.
const int P_INIT_U = 2;
// Initial polynomial degree for v.
const int P_INIT_V = 1;
// Number of initial boundary refinements
const int INIT_REF_BDY = 5;
// MULTI = true  ... use multi-mesh,
// MULTI = false ... use single-mesh->
// Note: In the single mesh option, the meshes are
// forced to be geometrically the same but the
// polynomial degrees can still vary.
const bool MULTI = true;
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.5;

// Error calculation & adaptivity.
DefaultErrorCalculator<double, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 2);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionLevels<double> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 5.;

// Problem parameters.
const double D_u = 1;
const double D_v = 1;
const double SIGMA = 1;
const double LAMBDA = 1;
const double KAPPA = 1;
const double K = 100.;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr u_mesh(new Mesh), v_mesh(new Mesh);
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", u_mesh);
  if (MULTI == false)
    u_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Create initial mesh (master mesh).
  v_mesh->copy(u_mesh);

  // Initial mesh refinements in the v_mesh towards the boundary.
  if (MULTI == true)
    v_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Set exact solutions.
  MeshFunctionSharedPtr<double> exact_u(new ExactSolutionFitzHughNagumo1 (u_mesh));
  MeshFunctionSharedPtr<double> exact_v(new ExactSolutionFitzHughNagumo2 (MULTI ? v_mesh : u_mesh, K));

  // Define right-hand sides.
  CustomRightHandSide1 g1(K, D_u, SIGMA);
  CustomRightHandSide2 g2(K, D_v);

  // Initialize the weak formulation.
  CustomWeakForm wf(&g1, &g2);

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u("Bdy", 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v("Bdy", 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(u_mesh, &bcs_u, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(MULTI ? v_mesh : u_mesh, &bcs_v, P_INIT_V));

  // Initialize coarse and reference mesh solutions.
  MeshFunctionSharedPtr<double> u_sln(new Solution<double>()), v_sln(new Solution<double>()), u_ref_sln(new Solution<double>()), v_ref_sln(new Solution<double>());
  std::vector<MeshFunctionSharedPtr<double> > slns({u_sln, v_sln});
  std::vector<MeshFunctionSharedPtr<double> > ref_slns({u_ref_sln, v_ref_sln});
  std::vector<MeshFunctionSharedPtr<double> > exact_slns({exact_u, exact_v});

  // Initialize refinement selector.
  H1ProjBasedSelector<double> selector(CAND_LIST);

  NewtonSolver<double> newton;
  newton.set_max_steps_with_reused_jacobian(0);
  newton.set_weak_formulation(&wf);

  // Adaptivity loop:
  int as = 1;
  bool done = false;

  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator u_ref_mesh_creator(u_mesh);
    MeshSharedPtr u_ref_mesh = u_ref_mesh_creator.create_ref_mesh();
    Mesh::ReferenceMeshCreator v_ref_mesh_creator(v_mesh);
    MeshSharedPtr v_ref_mesh = v_ref_mesh_creator.create_ref_mesh();
    Space<double>::ReferenceSpaceCreator u_ref_space_creator(u_space, u_ref_mesh);
    SpaceSharedPtr<double> u_ref_space = u_ref_space_creator.create_ref_space();
    Space<double>::ReferenceSpaceCreator v_ref_space_creator(v_space, MULTI ? v_ref_mesh : u_ref_mesh);
    SpaceSharedPtr<double> v_ref_space = v_ref_space_creator.create_ref_space();

    std::vector<SpaceSharedPtr<double> > ref_spaces_const({u_ref_space, v_ref_space});

    newton.set_spaces(ref_spaces_const);

    int ndof_ref = Space<double>::get_num_dofs(ref_spaces_const);

    // Perform Newton's iteration.
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      std::cout << e.what();
    }
    catch(std::exception& e)
    {
      std::cout << e.what();
    }

    // Translate the resulting coefficient vector into the instance of Solution.
    Solution<double>::vector_to_solutions(newton.get_sln_vector(), ref_spaces_const, std::vector<MeshFunctionSharedPtr<double> >(u_ref_sln, v_ref_sln));

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<double> ogProjection; ogProjection.project_global(std::vector<SpaceSharedPtr<double> >(u_space, v_space), ref_slns, slns);

    // Calculate element errors.
    errorCalculator.calculate_errors(slns, exact_slns, false);
    double err_exact_rel_total = errorCalculator.get_total_error_squared() * 100;

    errorCalculator.calculate_errors(slns, ref_slns, true);
    double err_est_rel_total = errorCalculator.get_total_error_squared() * 100;

    adaptivity.set_spaces(std::vector<SpaceSharedPtr<double> >(u_space, v_space));

    // If err_est too large, adapt the mesh->
    if (err_est_rel_total < ERR_STOP)
      done = true;
    else
    {
      std::vector<RefinementSelectors::Selector<double> *> selectors({&selector, &selector});
      done = adaptivity.adapt(selectors);
    }

    // Increase counter.
    as++;
  }

  while (done == false);

  bool success = true;

  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, -0.98)->val[0], 0.000986633, "u_ref_sln->get_pt_value(-0.98, -0.98)->val[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, -0.98)->val[0], .747675, "v_ref_sln->get_pt_value(-0.98, -0.98)->val[0]") && success;
  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, 0.98)->val[0], .000986633, "u_ref_sln->get_pt_value(-0.98, 0.98)->val[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, 0.98)->val[0], .747675, "v_ref_sln->get_pt_value(-0.98, 0.98)->val[0]") && success;

  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, -0.98)->dx[0], 0.0493155, "u_ref_sln->get_pt_value(-0.98, -0.98)->dx[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, -0.98)->dx[0], 11.701, "v_ref_sln->get_pt_value(-0.98, -0.98)->dx[0]") && success;
  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, 0.98)->dx[0], 0.0493155, "u_ref_sln->get_pt_value(-0.98, 0.98)->dx[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, 0.98)->dx[0], 11.701, "v_ref_sln->get_pt_value(-0.98, 0.98)->dx[0]") && success;

  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, -0.98)->dy[0], 0.0493155, "u_ref_sln->get_pt_value(-0.98, -0.98)->dy[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, -0.98)->dy[0], 11.701, "v_ref_sln->get_pt_value(-0.98, -0.98)->dy[0]") && success;
  success = Testing::test_value(u_ref_sln->get_pt_value(-0.98, 0.98)->dy[0], -0.0493155, "u_ref_sln->get_pt_value(-0.98, 0.98)->dy[0]") && success;
  success = Testing::test_value(v_ref_sln->get_pt_value(-0.98, 0.98)->dy[0], -11.701, "v_ref_sln->get_pt_value(-0.98, 0.98)->dy[0]") && success;

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
