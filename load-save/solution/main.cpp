#include "definitions.h"
#include "../../testing-core/testing-core.h"

int main(int argc, char* argv[])
{
  bool success = true;
  HermesCommonApi.set_integral_param_value(numThreads, 1);

  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  Hermes::Hermes2D::MeshReaderH2DXML mloader;
  mloader.load("domain.xml", mesh);

  mesh->refine_element_id(0);
  mesh->refine_element_id(2, 1);
  mesh->refine_all_elements();
  mesh->refine_all_elements();
  mesh->refine_all_elements();

  // Initialize essential boundary conditions.
  Hermes::Hermes2D::DefaultEssentialBCConst<double> bc_essential_r("Bdy", 0.);
  Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential_c("Bdy", complex(.132, -.12));
  Hermes::Hermes2D::EssentialBCs<double> bcs_r(&bc_essential_r);
  Hermes::Hermes2D::EssentialBCs<complex> bcs_c(&bc_essential_c);

  // Initialize space.
  SpaceSharedPtr<double> space_r( new Hermes::Hermes2D::H1Space<double>(mesh, &bcs_r, 2));
  SpaceSharedPtr<double> space_l2( new Hermes::Hermes2D::L2Space<double>(mesh, 0));
  SpaceSharedPtr<complex> space_c( new Hermes::Hermes2D::H1Space<complex>(mesh, &bcs_c, 2));

  // Initialize the weak formulation.
  WeakFormSharedPtr<double> wf_r(new WeakFormsH1::DefaultWeakFormPoissonLinear<double>(HERMES_ANY, nullptr));
  WeakFormSharedPtr<complex> wf_c(new WeakFormsH1::DefaultWeakFormPoissonLinear<complex>(HERMES_ANY, nullptr));

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln_r(new Solution<double>);
  MeshFunctionSharedPtr<double> sln1_r(new Solution<double>);
  MeshFunctionSharedPtr<double> sln2_r(new Solution<double>);

  MeshFunctionSharedPtr<complex> sln_c(new Solution<complex>);
  MeshFunctionSharedPtr<complex> sln1_c(new Solution<complex>);
  MeshFunctionSharedPtr<complex> sln2_c(new Solution<complex>);

  // Initialize linear solver.
  Hermes::Hermes2D::LinearSolver<double> linear_solver_r(wf_r, space_r);
  Hermes::Hermes2D::LinearSolver<complex> linear_solver_c(wf_c, space_c);

#ifdef SHOW_OUTPUT
  Views::ScalarView s;
#endif

  // 1st - save & load real solution.
  linear_solver_r.solve();
  Solution<double>::vector_to_solution(linear_solver_r.get_sln_vector(), space_r, sln_r);
  sln_r.get_solution()->save("saved_sln_r.xml");
  sln1_r.get_solution()->load("saved_sln_r.xml", space_r);
  sln1_r.get_solution()->save("saved_sln_r-final.xml");
#ifdef SHOW_OUTPUT
  s.set_title("Real - Original");
  s.show(sln_r);
  s.wait_for_keypress();
  s.set_title("Real - XML");
  s.show(sln1_r);
  s.wait_for_keypress();
#endif

#ifdef WITH_BSON
  sln1_r.get_solution()->save_bson("saved_sln_r.bson");
  sln2_r.get_solution()->load_bson("saved_sln_r.bson", space_r);
  sln2_r.get_solution()->save_bson("saved_sln_r-final.bson");
#ifdef SHOW_OUTPUT
  s.set_title("Real - BSON");
  s.show(sln2_r);
  s.wait_for_keypress();
#endif
#endif

  // 2nd - save & load complex one.
  linear_solver_c.solve();
  Solution<complex>::vector_to_solution(linear_solver_c.get_sln_vector(), space_c, sln_c);
  sln_c.get_solution()->save("saved_sln_c.xml");
  sln1_c.get_solution()->load("saved_sln_c.xml", space_c);
  sln1_c.get_solution()->save("saved_sln_c-final.xml");
#ifdef SHOW_OUTPUT
  MeshFunctionSharedPtr<double> filter(new RealFilter(sln_c));
  MeshFunctionSharedPtr<double> filter_1(new RealFilter(sln1_c));

  s.set_title("Complex - Original");
  s.show(filter);
  s.wait_for_keypress();
  s.set_title("Complex - XML");
  s.show(filter_1);
  s.wait_for_keypress();
#endif

#ifdef WITH_BSON
  sln1_c.get_solution()->save_bson("saved_sln_c.bson");
  sln2_c.get_solution()->load_bson("saved_sln_c.bson", space_c);
  sln2_c.get_solution()->save_bson("saved_sln_c-final.bson");
  MeshFunctionSharedPtr<double> filter_2(new RealFilter(sln2_c));
#ifdef SHOW_OUTPUT
  s.set_title("Complex - BSON");
  s.show(filter_2);
  s.wait_for_keypress();
#endif
#endif

  // 3th - save & load constant one.
  // 3.1 - bad one.
  MeshFunctionSharedPtr<double> initial_condition_bad(new CustomInitialCondition(mesh));
  bool local_success = false;
  try
  {
    initial_condition_bad.get_solution()->save("this-will-fail");
  }
  catch(Hermes::Exceptions::Exception& e)
  {
    local_success = true;
    std::string s = "Arbitrary exact solution can not be saved to disk. Only constant one can. Project to a space to get a saveable solution.";
    std::string s1 = e.info();
    if(s.compare(s1))
      success = false;
  }
  if(!local_success)
  {
    throw Exceptions::Exception("Exception not caught correctly.");
    success = false;
  }

  // 3.2 - good one.
  MeshFunctionSharedPtr<double> initial_condition_good(new ConstantSolution<double>(mesh, 3.1415926));
  MeshFunctionSharedPtr<double> test_solution_1(new Solution<double>());
  MeshFunctionSharedPtr<double> test_solution_2(new Solution<double>());
  initial_condition_good.get_solution()->save("constant_sln.xml");
  test_solution_1.get_solution()->load("constant_sln.xml", space_l2);
  test_solution_1.get_solution()->save("constant_sln-final.xml");
#ifdef SHOW_OUTPUT
  s.set_title("Exact - Original");
  s.show(initial_condition_good);
  s.wait_for_keypress();
  s.set_title("Exact - XML");
  s.show(test_solution_1);
  s.wait_for_keypress();
#endif
#ifdef WITH_BSON
  initial_condition_good.get_solution()->save_bson("constant_sln.bson");
  test_solution_2.get_solution()->load_bson("constant_sln.bson", space_l2);
  test_solution_2.get_solution()->save_bson("constant_sln-final.bson");
#ifdef SHOW_OUTPUT
  s.set_title("Exact - BSON");
  s.show(test_solution_2);
  s.wait_for_keypress();
#endif
#endif

/// \todo re-do the test files, but with some higher precision of stored numbers, this way it usually crashes on different truncation
/*
#if defined (_WINDOWS) || defined (WIN32) || defined (_MSC_VER)
  success = Testing::compare_files("saved_sln_r-final.xml", "win\\saved_sln_r-template.xml") && success;
  success = Testing::compare_files("saved_sln_r-final.bson", "win\\saved_sln_r-template.bson") &&  success;
  success = Testing::compare_files("saved_sln_c-final.xml", "win\\saved_sln_c-template.xml") &&  success;
  success = Testing::compare_files("saved_sln_c-final.bson", "win\\saved_sln_c-template.bson") &&  success;
  success = Testing::compare_files("constant_sln-final.xml", "win\\constant_sln-template.xml") &&  success;
  success = Testing::compare_files("constant_sln-final.bson", "win\\constant_sln-template.bson") && success;
#else
  success = Testing::compare_files("saved_sln_r-final.xml", "linux/saved_sln_r-template.xml") && success;
  success = Testing::compare_files("saved_sln_r-final.bson", "linux/saved_sln_r-template.bson") &&  success;
  success = Testing::compare_files("saved_sln_c-final.xml", "linux/saved_sln_c-template.xml") &&  success;
  success = Testing::compare_files("saved_sln_c-final.bson", "linux/saved_sln_c-template.bson") &&  success;
  success = Testing::compare_files("constant_sln-final.xml", "linux/constant_sln-template.xml") &&  success;
  success = Testing::compare_files("constant_sln-final.bson", "linux/constant_sln-template.bson") && success;
#endif
*/
  if(success)
  {
    printf("Success!");
    return 0;
  }
  else
  {
    printf("Failure!");
    return -1;
  }
}
