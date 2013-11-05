#include "definitions.h"
#include "../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
  Hermes::Mixins::TimeMeasurable time;
  time.tick();
  // Load the mesh.
  MeshSharedPtr mesh1(new Mesh);
  MeshSharedPtr mesh2(new Mesh);
  Hermes::vector<MeshSharedPtr> meshes;
  Hermes::vector<MeshSharedPtr> meshes1;
  meshes.push_back(mesh1);
  meshes1.push_back(mesh1);
  meshes1.push_back(mesh2);
  MeshReaderH2DXML mloader;
  bool exceptionCaughtCorrectly = false;
  mloader.set_validation(true);
  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-1.xml", meshes);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "instance document parsing failed\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  mloader.set_validation(false);

  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-1.xml", meshes);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Number of subdomains( = 2) does not equal the number of provided meshes in the vector( = 1).\n"))
      exceptionCaughtCorrectly = true;
  }
  
  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-2.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Some of the vertices of element #2 are identical which is not right.\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }
  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-3.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Some of the vertices of element #10 are identical which is not right.\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-4.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Boundary data #8: edge 156-5 does not exist.\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-5.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Curve #2: edge 111-6 does not exist.\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }
  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-6.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Wrong element number:123 in subdomain 0.\n"))
      exceptionCaughtCorrectly = true;
  }

  if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  try
  {
    exceptionCaughtCorrectly = false;
    mloader.load("bad-7.xml", meshes1);
  }
  catch(Hermes::Exceptions::MeshLoadFailureException& e)
  {
    std::cout << e.what();
    if(!strcmp(e.what(), "Boundary data error (edge 1 does not exist).\n"))
      exceptionCaughtCorrectly = true;
  }

	if(!exceptionCaughtCorrectly)
  {
    std::cout << "Exceptions not caught correctly";
    return -1;
  }

  try
  {
    mloader.load("initial.xml", meshes);
    mloader.save("saved.xml", meshes);
    mloader.load("saved.xml", meshes);

    mloader.load("bad-1.xml", meshes1);
    int nelem = mesh1->get_num_active_elements();
    mesh1->refine_all_elements();
    
    if(mesh1->get_num_active_elements() != 4 * nelem)
      throw Hermes::Exceptions::Exception("fail");
    mesh1->unrefine_all_elements();
    mloader.save("saved2.xml", meshes);
    mloader.load("saved2.xml", meshes);
    if(mesh1->get_num_active_elements() != nelem)
      throw Hermes::Exceptions::Exception("fail");
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    return -1;
  }

  try
  {
    mloader.set_validation(false);
    MeshSharedPtr mesh_agros1(new Mesh), mesh_agros2(new Mesh);
    Hermes::vector<MeshSharedPtr> meshes_agros(mesh_agros1, mesh_agros2);
    mloader.load("agros-test.msh", meshes_agros);
    DefaultEssentialBCConst<double> bc(HERMES_ANY, 1.0);
    EssentialBCs<double> bcs(&bc);
    SpaceSharedPtr<double> space_1(new H1Space<double>(mesh_agros1, &bcs, 1));
    SpaceSharedPtr<double> space_2(new H1Space<double>(mesh_agros2, &bcs, 1));
  }
  catch (std::exception& e)
  {
    std::cout << e.what();
    return -1;
  }

  time.tick();

  std::cout << "Time taken: " << time.last() << std::endl;
  return 0;
}
