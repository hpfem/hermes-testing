#include "../../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  RefMap rm;
  MeshFunctionSharedPtr<double> sln(new Solution<double>(mesh));
  
  Element* e;
  for(int ref = 0; ref < 5; ref++)
  {
    mesh->refine_all_elements();

    MeshFunctionSharedPtr<double> zsln(new ZeroSolution<double>(mesh));

    SpaceSharedPtr<double> space(new H1Space<double>(mesh, 3));
    double* coeff = (double*)calloc(space->get_num_dofs(), sizeof(double));
    Solution<double>::vector_to_solution(coeff, space, sln);

    for_all_active_elements(e, mesh)
    {
      rm.set_active_element(e);
      zsln->set_active_element(e);
      sln->set_active_element(e);
    }

    ::free(coeff);
  }
  return 0;
}