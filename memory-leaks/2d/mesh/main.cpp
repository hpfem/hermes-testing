#include "../../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshSharedPtr mesh2(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);
  mloader.save("domain2.xml", mesh);
  mloader.load("domain.xml", mesh2);
  mesh->free();
  mloader.load("domain.xml", mesh);

#ifdef WITH_BSON
  MeshReaderH2DBSON mloader_bson;

  MeshSharedPtr mesh_bson(new Mesh);
  MeshSharedPtr mesh_bson2(new Mesh);
  mloader_bson.save("domain.bson", mesh);
  mloader_bson.load("domain.bson", mesh_bson);
  mloader_bson.save("domain2.bson", mesh_bson2);
  mesh_bson->free();
  mloader_bson.load("domain2.bson", mesh_bson);
#endif

  return 0;
}
