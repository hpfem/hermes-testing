#include "../../../../testing-core/testing-core.h"

// This test makes sure that subdomains work correctly.

// Uniform polynomial degree of mesh elements.
const int P_INIT = 2;
// Number of initial uniform mesh refinements.
const int INIT_REF_NUM = 2;
// Possibilities: SOLVER_AMESOS, SOLVER_AZTECOO, SOLVER_MUMPS,
// SOLVER_PETSC, SOLVER_SUPERLU, SOLVER_UMFPACK.
Hermes::MatrixSolverType matrix_solver_type = Hermes::SOLVER_UMFPACK;

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh_whole_domain(new Mesh), mesh_bottom_left_corner(new Mesh), mesh_complement(new Mesh);
  std::vector<MeshSharedPtr > meshes({ mesh_bottom_left_corner, mesh_whole_domain, mesh_complement });
  MeshReaderH2DXML mloader;
  mloader.set_validation(false);
  mloader.load("subdomains.xml", meshes);

  // Perform initial mesh refinements (optional).
  for(int i = 0; i < INIT_REF_NUM; i++)
    for(unsigned int meshes_i = 0; meshes_i < meshes.size(); meshes_i++)
      meshes[meshes_i]->refine_all_elements();

  mloader.save("subdomains2.xml", meshes);
  mloader.load("subdomains2.xml", meshes);

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential_whole_domain(std::vector<std::string>({"Bottom Left", "Bottom Right", "Top Left", "Top Right"}), 0.0);
  EssentialBCs<double> bcs_whole_domain(&bc_essential_whole_domain);

  DefaultEssentialBCConst<double> bc_essential_bottom_left_corner(std::vector<std::string>({"Bottom Left", "Horizontal Left"}), 0.0);
  EssentialBCs<double> bcs_bottom_left_corner(&bc_essential_bottom_left_corner);

  DefaultEssentialBCConst<double> bc_essential_complement(std::vector<std::string>({"Bottom Right", "Top Right", "Top Left", "Horizontal Left", "Vertical Bottom"}), 0.0);
  EssentialBCs<double> bcs_complement(&bc_essential_complement);

  // Create H1 spaces with default shapeset.
  SpaceSharedPtr<double> space_whole_domain(new H1Space<double>(mesh_whole_domain, &bcs_whole_domain, P_INIT));
  int ndof_whole_domain = space_whole_domain->get_num_dofs();

  SpaceSharedPtr<double> space_bottom_left_corner(new H1Space<double>(mesh_bottom_left_corner, &bcs_bottom_left_corner, P_INIT));
  int ndof_bottom_left_corner = space_bottom_left_corner->get_num_dofs();

  SpaceSharedPtr<double> space_complement(new H1Space<double>(mesh_complement, &bcs_complement, P_INIT));
  int ndof_complement = space_complement->get_num_dofs();

  if(ndof_whole_domain == 225 && ndof_bottom_left_corner == 56 && ndof_complement == 161)
  {
    return 0;
  }
  else
  {
    return -1;
  }

  return 0;
}
