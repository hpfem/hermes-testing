#include "definitions.h"
#include "../../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);
  mloader.save("domain2.xml", mesh);

  return 0;
}
