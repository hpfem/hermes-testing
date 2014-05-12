#include "../../../testing-core/testing-core.h"

// Initial polynomial degree. NOTE: The meaning is different from
// standard continuous elements in the space H1. Here, P_INIT refers
// to the maximum poly order of the tangential component, and polynomials
// of degree P_INIT + 1 are present in element interiors. P_INIT = 0
// is for Whitney elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 1;

// Error calculation & adaptivity.
DefaultErrorCalculator<complex, HERMES_HCURL_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(0.3);
// Adaptivity processor class.
Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
double ERR_STOP = 10.;

// Problem parameters.
const double MU_R   = 1.0;
const double KAPPA  = 1.0;
const double LAMBDA = 1.0;

// Bessel functions, exact solution, and weak forms.
#include "definitions.cpp"

int main_element_type_spec(int argc, char* argv[], bool use_triangles)
{
  // Set just one thread.
  Hermes::HermesCommonApi.set_integral_param_value(numThreads, 1);

  if(!use_triangles)
    ERR_STOP = 1.0;
  else
    ERR_STOP = 10.0;

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2D mloader;
  if(use_triangles)
    mloader.load("lshape3t.mesh", mesh);  // triangles
  else
    mloader.load("lshape3q.mesh", mesh);    // quadrilaterals

  // Perform initial mesh refinemets.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential(std::vector<std::string>({"Corner_horizontal",
    "Corner_vertical"}), 0);
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an Hcurl space with default shapeset.
  SpaceSharedPtr<complex> space(new HcurlSpace<complex>(mesh, &bcs, P_INIT));

  // Initialize the weak formulation.
  WeakFormSharedPtr<complex> wf(new CustomWeakForm(MU_R, KAPPA));

  // Initialize coarse and reference mesh solutions.
  MeshFunctionSharedPtr<complex> sln(new Hermes::Hermes2D::Solution<complex>());
  MeshFunctionSharedPtr<complex> ref_sln(new Hermes::Hermes2D::Solution<complex>());

  // Initialize exact solution.
  MeshFunctionSharedPtr<complex> sln_exact(new CustomExactSolution(mesh));

  // Initialize refinement selector.
  HcurlProjBasedSelector<complex> selector(CAND_LIST);

  DiscreteProblem<complex> dp(wf, space);

  // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
  Hermes::Hermes2D::NewtonSolver<complex> newton(&dp);

  // Adaptivity loop:
  int as = 1; bool done = false;
  adaptivity.set_space(space);
  do
  {
    // Construct globally refined reference mesh and setup reference space->
    Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
    MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
    Space<complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = ref_space_creator.create_ref_space();

    newton.set_space(ref_space);
    int ndof_ref = ref_space->get_num_dofs();

    // Initial coefficient vector for the Newton's method.
    complex* coeff_vec = new complex[ndof_ref];
    memset(coeff_vec, 0, ndof_ref * sizeof(complex));
    
    try
    {
      newton.solve(coeff_vec);
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
    }
    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<complex> ogProjection;
    ogProjection.project_global(space, ref_sln, sln);

    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, sln_exact, false);
    double err_exact_rel = errorCalculator.get_total_error_squared() * 100;
    // Calculate exact error.
    errorCalculator.calculate_errors(sln, ref_sln);
    double err_est_rel = errorCalculator.get_total_error_squared() * 100;

    // If err_est_rel too large, adapt the mesh->
    if(err_est_rel < ERR_STOP)
      done = true;
    else
    {
      done = adaptivity.adapt(&selector);
      // Increase the counter of performed adaptivity steps.
      if(done == false)  as++;
    }

    // Clean up.
    delete [] coeff_vec;
  }
  while (done == false);

  complex sum = 0.;
  for (int i = 0; i < space->get_num_dofs(); i++)
    sum += newton.get_sln_vector()[i];

  complex expected_sum;
  bool success = true;
  if(use_triangles)
  {
    success = Testing::test_value(sum.real(), -3.4370807963401635, "sum - triangles - real", 0.5) && success; // Tested value as of May 2013.
    success = Testing::test_value(sum.imag(), -0.0065074801813934085, "sum - triangles - complex", 0.1) && success; // Tested value as of May 2013.
  }
  else
  {
    success = Testing::test_value(sum.real(), -0.414547, "sum - quads - real", 0.1) && success; // Tested value as of May 2013.
    success = Testing::test_value(sum.imag(), -0.0005088685276074, "sum - quads - complex", 0.1) && success; // Tested value as of May 2013.
  }

  int ndof = space->get_num_dofs();

  if(use_triangles)
    return Testing::test_value(703, ndof, "ndof - triangles") && success; // Tested value as of May 2013.
  else
    return Testing::test_value(988, ndof, "ndof - quads") && success; // Tested value as of May 2013.
}

int main(int argc, char* argv[])
{
  bool success = main_element_type_spec(argc, argv, false) && main_element_type_spec(argc, argv, true);

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