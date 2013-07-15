#include "definitions.h"
#include "../../testing-core/testing-core.h"

using namespace Hermes;
using namespace Hermes::Hermes2D;

int main(int argc, char* argv[])
{
#ifdef WITH_BSON
  Hermes::Mixins::TimeMeasurable time;
  time.tick();
  // Load the mesh.
  MeshSharedPtr mesh1(new Mesh);
  MeshReaderH2DXML mloader;
  MeshReaderH2DBSON mloaderb;
  Hermes::vector<MeshSharedPtr> meshes;
  meshes.push_back(mesh1);

  std::stringstream ss;
  ss << "saved" << std::rand() << ".bson";
  mloader.load("initial.xml", meshes);
  mloaderb.save(ss.str().c_str(), mesh1);
  mloaderb.load(ss.str().c_str(), mesh1);

  SpaceSharedPtr<double> space(new L2Space<double>(mesh1, 0));
  space->assign_dofs();

  MeshFunctionSharedPtr<double> fn(new ConstantSolution<double>(mesh1, 1.123124));
  Views::Linearizer linearizer;
  linearizer.save_solution_vtk(fn, "fn.vtk", "sln", true, 1, 1.0);

  bool success = Testing::compare_files(ss.str().c_str(), "reference.bson");

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
#endif      
  return 0;
}
