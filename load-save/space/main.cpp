#include "definitions.h"
#include "../../testing-core/testing-core.h"

namespace L2_real
{
  // Number of initial uniform mesh refinements.
  const int INIT_REF = 2;
  // Initial polynomial degrees of mesh elements in vertical and horizontal directions.
  const int P_INIT = 1;
  // This is a quantitative parameter of the adapt(...) function and
  // it has different meanings for various adaptive strategies.
  const double THRESHOLD = 0.9;

  // Error calculation & adaptivity.
  DefaultErrorCalculator<double, HERMES_L2_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<double> stoppingCriterion(THRESHOLD);
  // Adaptivity processor class.
  Adapt<double> adaptivity(&errorCalculator, &stoppingCriterion);
  // Predefined list of element refinement candidates.
  const CandList CAND_LIST = H2D_HP_ANISO;
  // Stopping criterion for adaptivity.
  const double ERR_STOP = 2e-1;

  int main(int argc, char* args[])
  {
    HermesCommonApi.set_integral_param_value(numThreads, 1);

    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("square.mesh", mesh);

    // Perform initial mesh refinement.
    for (int i=0; i<INIT_REF; i++)
      mesh->refine_all_elements();

    // Create an L2 space->
    SpaceSharedPtr<double> space(new L2Space<double>(mesh, P_INIT));

    // Initialize refinement selector.
    L2ProjBasedSelector<double> selector(CAND_LIST);

    // Display the mesh.
#ifdef SHOW_OUTPUT
    OrderView oview("Coarse mesh", new WinGeom(0, 0, 440, 350));
    oview.show(space);
#endif

    MeshFunctionSharedPtr<double> sln(new Solution<double>);
    MeshFunctionSharedPtr<double> ref_sln(new Solution<double>);

    // Initialize the weak formulation.
    WeakFormSharedPtr<double> wf(new CustomWeakForm("Bdy_bottom_left", mesh));

#ifdef SHOW_OUTPUT
    ScalarView view1("Solution", new WinGeom(900, 0, 450, 350));
    view1.fix_scale_width(60);
#endif

    // Initialize linear solver.
    Hermes::Hermes2D::LinearSolver<double> linear_solver(wf, space);

    int as = 1; bool done = false;
    do
    {
      // Construct globally refined reference mesh
      // and setup reference space->
      Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
      MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();

      ref_space->save("space-real.xml");
      ref_space->free();
      ref_space->load("space-real.xml");
#ifdef WITH_BSON
      ref_space->save_bson("space-real.bson");
      ref_space->free();
      ref_space->load_bson("space-real.bson");
#endif

      linear_solver.set_space(ref_space);

      // Solve the linear system. If successful, obtain the solution.
      linear_solver.solve();
      Solution<double>::vector_to_solution(linear_solver.get_sln_vector(), ref_space, ref_sln);
      
      // Project the fine mesh solution onto the coarse mesh.
      OGProjection<double> ogProjection;
      ogProjection.project_global(space, ref_sln, sln, HERMES_L2_NORM);

#ifdef SHOW_OUTPUT
      MeshFunctionSharedPtr<double> val_filter(new ValFilter(ref_sln, 0.0, 1.0));

      // View the coarse mesh solution.
      view1.show(val_filter);
      oview.show(space);
#endif

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(sln, ref_sln);
      double err_est_rel = errorCalculator.get_total_error_squared() * 100;

      adaptivity.set_space(space);
#ifdef SHOW_OUTPUT
      std::cout << "Error: " << err_est_rel << "%." << std::endl;
#endif
      // If err_est_rel too large, adapt the mesh->
      if(err_est_rel < ERR_STOP)
        done = true;
      else
        done = adaptivity.adapt(&selector);
      as++;
    }
    while (done == false);

    // Wait for keyboard or mouse input.
#ifdef SHOW_OUTPUT
    View::wait();
#endif
    return as;
  }
}

namespace H1_complex
{
typedef std::complex<double> complex;
  const int INIT_REF_NUM = 0;                       // Number of initial uniform mesh refinements.
  const int P_INIT = 1;                             // Initial polynomial degree of all mesh elements.
  const double THRESHOLD = 0.95;                    // This is a quantitative parameter of Adaptivity.

  // Error calculation & adaptivity.
  DefaultErrorCalculator<complex, HERMES_H1_NORM> errorCalculator(RelativeErrorToGlobalNorm, 1);
  // Stopping criterion for an adaptivity step.
  AdaptStoppingCriterionSingleElement<complex> stoppingCriterion(0.9);
  // Adaptivity processor class.
  Adapt<complex> adaptivity(&errorCalculator, &stoppingCriterion);
  // Predefined list of element refinement candidates.
  const CandList CAND_LIST = H2D_HP_ANISO;
  // Stopping criterion for adaptivity.
  const double ERR_STOP = 0.8;

  // Problem parameters.
  const double MU_0 = 4.0*M_PI*1e-7;
  const double MU_IRON = 1e3 * MU_0;
  const double GAMMA_IRON = 6e6;
  const double J_EXT = 1e6;
  const double FREQ = 5e3;
  const double OMEGA = 2 * M_PI * FREQ;


  int main(int argc, char* argv[])
  {
    // Load the mesh.
    MeshSharedPtr mesh(new Mesh);
    MeshReaderH2D mloader;
    mloader.load("domain.mesh", mesh);

    // Perform initial mesh refinements.
    for (int i = 0; i < INIT_REF_NUM; i++) mesh->refine_all_elements();

    // Initialize boundary conditions.
    Hermes::Hermes2D::DefaultEssentialBCConst<complex> bc_essential("Dirichlet", complex(0.0, 0.0));
    EssentialBCs<complex> bcs(&bc_essential);

    // Create an H1 space with default shapeset.
    SpaceSharedPtr<complex> space(new H1Space<complex>(mesh, &bcs, P_INIT));
    int ndof = space->get_num_dofs();

    // Initialize the weak formulation.
    WeakFormSharedPtr<complex> wf(new CustomWeakForm("Air", MU_0, "Iron", MU_IRON, GAMMA_IRON,
      "Wire", MU_0, complex(J_EXT, 0.0), OMEGA));

    // Initialize coarse and reference mesh solution.
    MeshFunctionSharedPtr<complex> sln(new Hermes::Hermes2D::Solution<complex>());
    MeshFunctionSharedPtr<complex> ref_sln(new Hermes::Hermes2D::Solution<complex>());

    // Initialize refinement selector.
    H1ProjBasedSelector<complex> selector(CAND_LIST);

    // Initialize views.
#ifdef SHOW_OUTPUT
    Views::ScalarView sview("Solution", new Views::WinGeom(0, 0, 600, 350));
    Views::ScalarView sview2("Ref. Solution", new Views::WinGeom(0, 0, 600, 350));
    Views::OrderView oview("Polynomial orders", new Views::WinGeom(610, 0, 520, 350));
#endif

    DiscreteProblem<complex> dp(wf, space);

    // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
    Hermes::Hermes2D::NewtonSolver<complex> newton(&dp);

    // Adaptivity loop:
    int as = 1; bool done = false;
    adaptivity.set_space(space);
    do
    {
      // Construct globally refined reference mesh and setup reference space->
      Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
      MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<complex>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
      SpaceSharedPtr<complex> ref_space = ref_space_creator.create_ref_space();

      newton.set_space(ref_space);
      ref_space->save("space-complex.xml");
      ref_space->free();
      ref_space->load("space-complex.xml");
#ifdef WITH_BSON
      ref_space->save_bson("space-complex.bson");
      ref_space->free();
      ref_space->load_bson("space-complex.bson");
#endif

      int ndof_ref = ref_space->get_num_dofs();

      // Initialize reference problem.

      // Initial coefficient vector for the Newton's method.
      complex* coeff_vec = new complex[ndof_ref];
      memset(coeff_vec, 0, ndof_ref * sizeof(complex));

      // Perform Newton's iteration and translate the resulting coefficient vector into a Solution.
      SpaceSharedPtr<complex> space_test = Space<complex>::load("space-complex.xml", ref_mesh, false, &bcs);
      newton.set_space(space_test);
      newton.solve(coeff_vec);

      Hermes::Hermes2D::Solution<complex>::vector_to_solution(newton.get_sln_vector(), ref_space, ref_sln);

      // Project the fine mesh solution onto the coarse mesh.
      OGProjection<complex> ogProjection;

#ifdef WITH_BSON
      space->save_bson("space-complex-coarse.bson");
      SpaceSharedPtr<complex> space_test2 = Space<complex>::load_bson("space-complex-coarse.bson", mesh, &bcs);
      ogProjection.project_global(space_test2, ref_sln, sln);
#else
      space->save("space-complex-coarse.xml2");
      SpaceSharedPtr<complex> space_test2 = Space<complex>::load("space-complex-coarse.xml2", mesh, false, &bcs);
      ogProjection.project_global(space_test2, ref_sln, sln);
#endif
      // View the coarse mesh solution and polynomial orders.
#ifdef SHOW_OUTPUT
      MeshFunctionSharedPtr<double> real_filter(new RealFilter(sln));
      MeshFunctionSharedPtr<double> rreal_filter(new RealFilter(ref_sln));
      sview2.show(rreal_filter);
      oview.show(space);
#endif

      // Calculate element errors and total error estimate.
      errorCalculator.calculate_errors(sln, ref_sln);

#ifdef SHOW_OUTPUT
      std::cout << "Relative error: " << errorCalculator.get_total_error_squared() * 100. << '%' << std::endl;
#endif

      // Add entry to DOF and CPU convergence graphs.
#ifdef SHOW_OUTPUT
      sview.show(errorCalculator.get_errorMeshFunction());
#endif

      // If err_est too large, adapt the mesh->
      if(errorCalculator.get_total_error_squared()  * 100. < ERR_STOP)
        done = true;
      else
      {
        std::cout << "Adapting..." << std::endl << std::endl;
        adaptivity.adapt(&selector);
      }

      // Clean up.
      delete [] coeff_vec;

      // Increase counter.
      as++;
    }
    while (done == false);

#ifdef SHOW_OUTPUT
    // Show the reference solution - the final result.
    sview.set_title("Fine mesh solution");
    MeshFunctionSharedPtr<double> real_filter(new RealFilter(ref_sln));
    sview.show(real_filter);

    // Wait for all views to be closed.
    Views::View::wait();
#endif
    return as;
  }
}

int main(int argc, char* args[])
{
  // Measure total adaptivity steps.
  // If something goes wrong, this will probably differ.
  // But mostly, if anything goes really wrong, there will be an exception, or a segfault.

  int total_check = 0;
  total_check += L2_real::main(argc, args);
  total_check += H1_complex::main(argc, args);

  bool success = Hermes::Testing::test_value(total_check, 18, "Adaptivity steps", 1);
  if (total_check == 18)
  {
    printf("Success!\n");
    return 0;
  }
  else
  {
    printf("Failure!\n");
    return -1;
  }

  return 0;
}