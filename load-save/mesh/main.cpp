#define HERMES_REPORT_ALL
#include "definitions.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

int main(int argc, char* argv[])
{
  Hermes::Mixins::TimeMeasurable time;
  time.tick();
  // Load the mesh.
  Hermes::Hermes2D::Mesh mesh;
  Hermes::vector<Mesh*> meshes;
  meshes.push_back(&mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.set_validation(false);
  try
  {
    mloader.load("initial.xml", meshes);
    int nelem = mesh.get_num_active_elements();
    mesh.refine_all_elements();
    mloader.save("saved.xml", meshes);
    mloader.load("saved.xml", meshes);
    if(mesh.get_num_active_elements() != 4 * nelem)
      throw Hermes::Exceptions::Exception("fail");
    mesh.unrefine_all_elements();
    mloader.save("saved2.xml", meshes);
    mloader.load("saved2.xml", meshes);
    if(mesh.get_num_active_elements() != nelem)
      throw Hermes::Exceptions::Exception("fail");
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
    return -1;
  }

  time.tick();

  std::cout << "Time taken: " << time.last() << std::endl;
  return 0;
}
