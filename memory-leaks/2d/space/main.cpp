#include "../../../testing-core/testing-core.h"

namespace L2_real
{
  // Number of initial uniform mesh refinements.
  const int INIT_REF = 2;

  int main(int argc, char* args[])
  {
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("square.mesh", mesh);

    // Perform initial mesh refinement.
    for (int i = 0; i < INIT_REF; i++)
      mesh->refine_all_elements();

    int orders[2] = { 3, 2 };
    for (int k = 0; k < 2; k++)
    {
      // Create an L2 space.
      SpaceSharedPtr<double> space(new L2Space<double>(mesh, orders[k]));

      int as = 1; bool done = false;
      do
      {
        // Construct globally refined reference mesh
        // and setup reference space->
        Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
        MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
        Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
        SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

        if (as == 3)
        {
          space->copy(ref_space, ref_mesh);

          space->save("space-real.xml");
          space->free();
          space->load("space-real.xml");
#ifdef WITH_BSON
          space->save_bson("space-real.bson");
          space->free();
          space->load_bson("space-real.bson");
#endif
        }
        else
        {
          ref_space->save("space-real.xml");
          ref_space->free();
          ref_space->load("space-real.xml");
#ifdef WITH_BSON
          ref_space->save_bson("space-real.bson");
          ref_space->free();
          ref_space->load_bson("space-real.bson");
#endif
        }
      } while (as++ < 3);
    }
    return 0;
  }
}

namespace H1_complex
{
  typedef std::complex<double> complex;
  const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.

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
    Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential("Dirichlet", complex(0.0, 0.0));
    EssentialBCs<complex> bcs(&bc_essential);

    int orders[2] = { 1, 2 };
    for (int k = 0; k < 2; k++)
    {
      // Create an H1 space with default shapeset.
      SpaceSharedPtr<complex> space(new H1Space<complex>(mesh, &bcs, orders[k]));
      space->update_essential_bc_values();
      int ndof = space->get_num_dofs();

      // Adaptivity loop:
      int as = 1; bool done = false;
      do
      {
        // Construct globally refined reference mesh and setup reference space->
        Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
        MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
        Space<complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
        SpaceSharedPtr<complex> ref_space = ref_space_creator.create_ref_space();
        ref_space->update_essential_bc_values();

        ref_space->save("space-complex.xml");
        ref_space->free();
        ref_space->load("space-complex.xml");
#ifdef WITH_BSON
        ref_space->save_bson("space-complex.bson");
        ref_space->free();
        ref_space->load_bson("space-complex.bson");
#endif

        int ndof_ref = ref_space->get_num_dofs();

        // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
        SpaceSharedPtr<complex> space_test = Space<complex>::load("space-complex.xml", ref_mesh, false, &bcs);

#ifdef WITH_BSON
        space->save_bson("space-complex-coarse.bson");
        SpaceSharedPtr<complex> space_test2 = Space<complex>::load_bson("space-complex-coarse.bson", mesh, &bcs);
#else
        space->save("space-complex-coarse.xml2");
        SpaceSharedPtr<complex> space_test2 = Space<complex>::load("space-complex-coarse.xml2", mesh, false, &bcs);
#endif
      } while (as++ < 4);
    }

    return 0;
  }
}

int main(int argc, char* args[])
{
  L2_real::main(argc, args);
  H1_complex::main(argc, args);

  return 0;
}