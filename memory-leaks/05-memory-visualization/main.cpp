#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "definitions.h"

//  This problem describes the distribution of the vector potential in
//  a 2D domain comprising a wire carrying electrical current, air, and
//  an iron which is not under voltage.
//
//  PDE: -div(1/rho grad p) - omega**2 / (rho c**2) * p = 0.
//
//  Domain: Floor plan of an existing apartment.
//
//  BC: Prescribed pressure at the source.
//      Newton matched boundary on the walls.
//
//  The following parameters can be changed:

// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 0;                       
// Initial polynomial degree of mesh elements.
const int P_INIT = 2;

// Problem parameters.
const double RHO = 1.25;
const double FREQ = 5e2;
const double OMEGA = 2 * M_PI * FREQ;
const double SOUND_SPEED = 353.0;
const std::complex<double> P_SOURCE(1.0, 0.0);

int main(int argc, char* argv[])
{
  // Load the mesh.
  Mesh mesh;
  MeshReaderH2D mloader;
  mloader.load("domain.mesh", &mesh);

  // Perform initial mesh refinements.
  for (int i = 0; i < INIT_REF_NUM; i++) mesh.refine_all_elements();

  // Initialize boundary conditions.
  DefaultEssentialBCConst<std::complex<double> > bc_essential("Source", P_SOURCE);
  EssentialBCs<std::complex<double> > bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  H1Space<std::complex<double> > space(&mesh, &bcs, P_INIT);
  int ndof = Space<std::complex<double> >::get_num_dofs(&space);
  Hermes::Mixins::Loggable::Static::info("ndof = %d", ndof);

  // Initialize the weak formulation.
  CustomWeakFormAcoustics wf("Wall", RHO, SOUND_SPEED, OMEGA);

  // Assemble the reference problem.
  DiscreteProblem<std::complex<double> > dp(&wf, &space);

  // Perform Newton's iteration.
  Hermes::Hermes2D::NewtonSolver<std::complex<double> > newton(&dp);

  //Tests.
  Linearizer* lin = new Linearizer();
  Orderizer* ord = new Orderizer();
  Vectorizer* vec = new Vectorizer();

  // Adaptivity loop:
  bool done = false;
  for(int i = 0; i < 2; i++)
  {
    try
    {
      newton.solve();
    }
    catch(Hermes::Exceptions::Exception e)
    {
      e.print_msg();
    };
    
    // Initialize coarse and reference mesh solution.
    Solution<std::complex<double> >* sln = new Solution<std::complex<double> >();

    // Translate the resulting coefficient vector into the Solution<std::complex<double> > sln.
    Hermes::Hermes2D::Solution<std::complex<double> >::vector_to_solution(newton.get_sln_vector(), &space, sln);

    RealFilter* ref_mag = new RealFilter(sln);
    AbsFilter* abs_ref_mag = new AbsFilter(ref_mag);
    TopValFilter* top_abs_ref_mag = new TopValFilter(Hermes::vector<MeshFunction<double>*>(ref_mag, abs_ref_mag), Hermes::vector<double>(0.4, 0.5));

    std::cout << "Linearizer #1" << std::endl;
    lin->process_solution(ref_mag);
    std::cout << "Linearizer #2" << std::endl;
    lin->process_solution(abs_ref_mag);
    std::cout << "Linearizer #3" << std::endl;
    lin->process_solution(top_abs_ref_mag);

    std::cout << "Orderizer #1" << std::endl;
    ord->process_space(&space);
    Mesh::ReferenceMeshCreator ref_mesh_creator(&mesh);
		Mesh* ref_mesh = ref_mesh_creator.create_ref_mesh();
		Space<std::complex<double> >::ReferenceSpaceCreator ref_space_creator(&space, ref_mesh);
		Space<std::complex<double> >* ref_space = ref_space_creator.create_ref_space();
    std::cout << "Orderizer #2" << std::endl;
    ord->process_space(ref_space);
    delete ref_space;
    delete ref_mesh;

    std::cout << "Vectorizer #1" << std::endl;
    vec->process_solution(ref_mag, abs_ref_mag, H2D_FN_VAL_0, H2D_FN_VAL_0);
    std::cout << "Vectorizer #2" << std::endl;
    vec->process_solution(top_abs_ref_mag, abs_ref_mag, H2D_FN_VAL_0, H2D_FN_VAL_0);
  
    delete top_abs_ref_mag;
    delete abs_ref_mag;
    delete ref_mag;
    delete sln;
  }

  delete lin;
  delete ord;
  delete vec;

  // Wait for all views to be closed.
  View::wait();
  return 0;
}
