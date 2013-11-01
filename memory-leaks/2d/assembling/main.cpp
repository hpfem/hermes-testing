#include "definitions.h"

const int P_INIT_U = 2;
const int P_INIT_V = 3;
const int INIT_REF_BDY = 4;

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
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", u_mesh);
  u_mesh->refine_all_elements();
  
  v_mesh->copy(u_mesh);
  v_mesh->refine_towards_boundary("Bdy", INIT_REF_BDY);

  // Define right-hand sides.
  CustomRightHandSide1* g1 = new CustomRightHandSide1(K, D_u, SIGMA);
  CustomRightHandSide2* g2 = new CustomRightHandSide2(K, D_v);

  // Initialize the weak formulation.
  CustomWeakForm wf(g1, g2);

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u("Bdy", 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v("Bdy", 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(u_mesh, &bcs_u, P_INIT_U));
  SpaceSharedPtr<double> v_space(new H1Space<double>(v_mesh, &bcs_v, P_INIT_V));
  Hermes::vector<SpaceSharedPtr<double> > spaces(u_space, v_space);
  NewtonSolver<double> newton(&wf, spaces);

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

  return 0;
}
