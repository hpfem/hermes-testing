#include "definitions.h"
#include "../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
#ifdef WITH_BSON
#pragma region Initialization
  MeshSharedPtr mesh1(new Mesh);
  MeshSharedPtr mesh2(new Mesh);
  MeshSharedPtr mesh3(new Mesh);
  MeshSharedPtr mesh4(new Mesh);
  MeshReaderH2DXML mloader;
  MeshReaderH2DBSON mloaderb;
  std::vector<MeshSharedPtr > meshes;
  std::vector<MeshSharedPtr > meshes1;
  std::vector<MeshSharedPtr > meshes2;
  meshes.push_back({ mesh1 });
  meshes1.push_back(mesh2);
  meshes1.push_back(mesh3);
  meshes1.push_back(mesh4);
  meshes2.push_back(mesh3);
  meshes2.push_back(mesh4);
  std::stringstream ss;
  std::stringstream ss1;
  std::stringstream ss2;
  ss << "saved" << std::rand() << ".bson";
  ss1 << "saved" << (2 + std::rand()) << ".bson";
  ss2 << "saved" << (12 + std::rand()) << ".bson";
  Views::Linearizer linearizer(FileExport);
#pragma endregion
#pragma region Single mesh
  mloader.load("initial.xml", meshes);
  mloaderb.save(ss.str().c_str(), mesh1);
  mloaderb.load(ss.str().c_str(), mesh1);

  SpaceSharedPtr<double> space(new L2Space<double>(mesh1, 0));
  space->assign_dofs();

  MeshFunctionSharedPtr<double> fn(new ConstantSolution<double>(mesh1, 1.123124));
  linearizer.save_solution_vtk(fn, "fn.vtk", "sln", true, 1);
#pragma endregion

#pragma region Multi-mesh 1
  mloader.load("subdomains.xml", meshes1);
  mloaderb.save(ss1.str().c_str(), meshes1);
  mloaderb.load(ss1.str().c_str(), meshes1);

  SpaceSharedPtr<double> space2(new L2Space<double>(mesh2, 0));
  SpaceSharedPtr<double> space21(new L2Space<double>(mesh3, 0));
  SpaceSharedPtr<double> space22(new L2Space<double>(mesh4, 0));

  MeshFunctionSharedPtr<double> fn2(new ConstantSolution<double>(mesh2, 1.123124));
  linearizer.save_solution_vtk(fn2, "fn2.vtk", "sln", true, 1);

  MeshFunctionSharedPtr<double> fn3(new ConstantSolution<double>(mesh3, 1.123124));
  linearizer.save_solution_vtk(fn3, "fn3.vtk", "sln", true, 1);

  MeshFunctionSharedPtr<double> fn4(new ConstantSolution<double>(mesh4, 1.123124));
  linearizer.save_solution_vtk(fn4, "fn4.vtk", "sln", true, 1);
#pragma endregion

#pragma region Multi-mesh 2
  mloader.load("subdomains2.xml", meshes2);
  mloaderb.save(ss2.str().c_str(), meshes2);
  mloaderb.load(ss2.str().c_str(), meshes2);

  MeshFunctionSharedPtr<double> fn33(new ConstantSolution<double>(mesh3, 1.123124));
  linearizer.save_solution_vtk(fn33, "fn33.vtk", "sln", true, 1);

  MeshFunctionSharedPtr<double> fn44(new ConstantSolution<double>(mesh4, 1.123124));
  linearizer.save_solution_vtk(fn44, "fn44.vtk", "sln", true, 1);

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf(new CustomWeakForm);

  // Initialize boundary conditions
  DefaultEssentialBCConst<double> bc_u(HERMES_ANY, 0.0);
  EssentialBCs<double> bcs_u(&bc_u);
  DefaultEssentialBCConst<double> bc_v(HERMES_ANY, 0.0);
  EssentialBCs<double> bcs_v(&bc_v);

  // Create H1 spaces with default shapeset for both displacement components.
  SpaceSharedPtr<double> u_space(new H1Space<double>(mesh3, &bcs_u, 2));
  SpaceSharedPtr<double> v_space(new H1Space<double>(mesh4, &bcs_v, 2));
  LinearSolver<double> solver(wf, { u_space, v_space });
  solver.solve();
  double norm = get_l2_norm(solver.get_sln_vector(), Space<double>::get_num_dofs({u_space, v_space}));

#ifdef SHOW_OUTPUT
  MeshFunctionSharedPtr<double> sln1(new Solution<double>());
  MeshFunctionSharedPtr<double> sln2(new Solution<double>());
  Solution<double>::vector_to_solutions(solver.get_sln_vector(), {u_space, v_space}, {sln1, sln2});

  Views::ScalarView s;
  s.show(sln1);
  s.wait_for_keypress();
  s.show(sln2);
  s.wait_for_keypress();
#endif
#pragma endregion

  bool success = Testing::test_value(norm, 25.002510204555620, "solution-norm", 1e-6);
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
