#include "../../../testing-core/testing-core.h"

using namespace Hermes::Solvers;
using namespace Hermes::Algebra;

// This is a test of excptions.

int main(int argc, char* argv[])
{
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  //NullException test
  try
  {
    ((Solution<double>*)sln.get())->get_ref_value(nullptr,0,0,0,0);
    std::cout << "Failure - get_ref_value!";
    return -1;
  }
  catch(Exceptions::NullException& e)
  {
    e.print_msg();
  }

  //LengthException test
  double solution_vector[3];
  std::vector<SpaceSharedPtr<double> > spaces({ nullptr, nullptr, nullptr, nullptr });
  std::vector<MeshFunctionSharedPtr<double> > solutions({ nullptr, nullptr, nullptr });
  try
  {
    Solution<double>::vector_to_solutions(solution_vector,spaces,solutions);
    std::cout << "Failure - vector_to_solutions!";
    return -1;
  }
  catch(Exceptions::LengthException& e)
  {
    e.print_msg();
  }

  //1/2Exception test

  CSCMatrix<double> mat;
  int ap[]={0,1,1};
  int ai[]={0};
  double ax[]={0.0};
  mat.create(2,1,ap,ai,ax);
  SimpleVector<double> vec(2);

  UMFPackLinearMatrixSolver<double> linsolv(&mat,&vec);
  try
  {
    linsolv.solve();
    std::cout << "Failure - algebra!";
    return -1;
  }
  catch(Exceptions::LinearMatrixSolverException& e)
  {
    e.print_msg();
  }

  //ValueException test
  std::vector<SpaceSharedPtr<double> > spaces2;
  std::vector<Hermes2D::NormType> proj_norms;
  for (int i=0;i>H2D_MAX_COMPONENTS+1;i++)
  {
    spaces2.push_back(nullptr);
    proj_norms.push_back(Hermes2D::HERMES_UNSET_NORM);
  }

  try
  {
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2DXML reader;
    reader.load("domain.xml", mesh);
    std::cout << "Failure - mesh!";
    return -1;
  }
  catch(Exceptions::MeshLoadFailureException& e)
  {
    e.print_msg();
  }

  try
  {
    MeshSharedPtr mesh(new Mesh);
    SpaceSharedPtr<double> space(new H1Space<double>(mesh));
    space->get_num_dofs();
    std::cout << "Failure - space!";
    return -1;
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }

  try
  {
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("domain.mesh", mesh);

    // Create an H1 space with default shapeset.
    SpaceSharedPtr<double> space(new L2Space<double>(mesh, 3));

    LinearSolver<double> ls;
    ls.set_space(space);
    ls.solve();
    std::cout << "Failure - solver!";
    return -1;
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    e.print_msg();
  }

  std::cout << "Success!";
  return 0;
}
