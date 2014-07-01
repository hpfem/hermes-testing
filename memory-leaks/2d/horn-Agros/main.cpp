#include "definitions.h"

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;                             
// This is a quantitative parameter of the adapt(...) function and
// it has different meanings for various adaptive strategies.
const double THRESHOLD = 0.3;                     
// Error calculation & adaptivity.
DefaultErrorCalculator<complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
// Stopping criterion for an adaptivity step.
AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(THRESHOLD);
// Adaptivity processor class.
Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);
// Predefined list of element refinement candidates.
const CandList CAND_LIST = H2D_HP_ANISO;
// Stopping criterion for adaptivity.
const double ERR_STOP = 1e-1;

// Problem parameters.
const double RHO = 1.25;
const double FREQ = 5e3;
const double OMEGA = 2 * M_PI * FREQ;
const double SOUND_SPEED = 353.0;
const std::complex<double>  P_SOURCE(1.0, 0.0);

int main(int argc, char* argv[])
{
#ifdef THREAD_TESTING
  HermesCommonApi.set_integral_param_value(numThreads, 8);
#endif
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  std::vector<MeshSharedPtr > meshes;
  meshes.push_back({mesh});
  MeshReaderH2DXML mloader;
  mloader.load("agrosMesh.msh", meshes);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++)
    mesh->refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<complex> bc_essential("4", P_SOURCE);
  EssentialBCs<complex> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<complex> space(new H1Space<complex> (mesh, &bcs, P_INIT));
  adaptivity.set_space(space);

  // Initialize the weak formulation.
  WeakFormSharedPtr<complex> wf(new CustomWeakFormAcoustics("0", RHO, SOUND_SPEED, OMEGA));

  // Initialize coarse and reference mesh solution.
  MeshFunctionSharedPtr<complex>  sln(new Solution<complex>), ref_sln(new Solution<complex>);

  // Initialize refinement selector.
  H1ProjBasedSelector<complex> selector(CAND_LIST);

  Hermes::Hermes2D::NewtonSolver<complex> newton;
  newton.set_weak_formulation(wf);

  // 2 Adaptivity steps:
  int as = 1;
  bool done = false;
  do
  {
    // Construct globally refined reference mesh and setup reference space.
    Mesh::ReferenceMeshCreator refMeshCreator(mesh);
    MeshSharedPtr ref_mesh = refMeshCreator.create_ref_mesh();

    Space<complex>::ReferenceSpaceCreator refSpaceCreator(space, ref_mesh);
    SpaceSharedPtr<complex> ref_space = refSpaceCreator.create_ref_space();

    // Perform Newton's iteration.
    try
    {
      newton.set_space(ref_space);
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception& e)
    {
      e.print_msg();
      throw Hermes::Exceptions::Exception("Newton's iteration failed.");
    };

    // Translate the resulting coefficient vector into the Solution<complex> sln->
    Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

    // Project the fine mesh solution onto the coarse mesh.
    OGProjection<complex> ogProjection; ogProjection.project_global(space, ref_sln, sln);
    
    // Calculate element errors and total error estimate.
    errorCalculator.calculate_errors(sln, ref_sln);
    adaptivity.adapt(&selector);
  }
  while (as++ < 2);
  return 0;
}
