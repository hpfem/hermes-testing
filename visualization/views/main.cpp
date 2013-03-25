#define HERMES_REPORT_ALL
#include "definitions.h"

// This example shows how to solve a simple PDE that describes stationary
// heat transfer in an object consisting of two materials (aluminum and
// copper). The object is heated by constant volumetric heat sources
// (generated, for example, by a DC electric current). The temperature
// on the boundary is fixed. We will learn how to:
//
//   - load the mesh,
//   - perform initial refinements,
//   - create a H1 space over the mesh,
//   - define weak formulation,
//   - initialize matrix solver,
//   - assemble and solve the matrix system,
//   - output the solution and element orders in VTK format
//     (to be visualized, e.g., using Paraview),
//   - visualize the solution using Hermes' native OpenGL-based functionality.
//
// PDE: Poisson equation -div(LAMBDA grad u) - VOLUME_HEAT_SRC = 0.
//
// Boundary conditions: Dirichlet u(x, y) = FIXED_BDY_TEMP on the boundary.
//
// Geometry: L-Shape domain (see file domain.mesh).
//
// The following parameters can be changed:

const bool HERMES_VISUALIZATION = true;           // Set to "false" to suppress Hermes OpenGL visualization.
const bool VTK_VISUALIZATION = false;              // Set to "true" to enable VTK output.
const int P_INIT = 1;                             // Uniform polynomial degree of mesh elements.
const int INIT_REF_NUM = 1;                       // Number of initial uniform mesh refinements.

// Problem parameters.
const double LAMBDA_AL = 236.0;            // Thermal cond. of Al for temperatures around 20 deg Celsius.
const double LAMBDA_CU = 386.0;            // Thermal cond. of Cu for temperatures around 20 deg Celsius.
const double VOLUME_HEAT_SRC = 5e2;        // Volume heat sources generated (for example) by electric current.
const double FIXED_BDY_TEMP = 20.0;        // Fixed temperature on the boundary.

int main(int argc, char* argv[])
{
  // Load the mesh.
  MeshSharedPtr mesh(new Mesh);
  MeshReaderH2DXML mloader;
  mloader.set_validation(false);
  try
  {
    mloader.load("domain.xml", mesh);
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  // Perform initial mesh refinements (optional).
  mesh->refine_in_areas(Hermes::vector<std::string>("Aluminum", "Copper"), INIT_REF_NUM);
  mesh->refine_in_area("Aluminum");

  // Initialize the weak formulation.
  CustomWeakFormPoisson wf("Aluminum", new Hermes::Hermes1DFunction<double>(LAMBDA_AL), "Copper",
    new Hermes::Hermes1DFunction<double>(LAMBDA_CU), new Hermes::Hermes2DFunction<double>(-VOLUME_HEAT_SRC));

  // Initialize essential boundary conditions.
  DefaultEssentialBCConst<double> bc_essential(Hermes::vector<std::string>("Bottom", "Inner", "Outer", "Left"),
    FIXED_BDY_TEMP);
  EssentialBCs<double> bcs(&bc_essential);

  // Create an H1 space with default shapeset.
  SpaceSharedPtr<double> space(new H1Space<double>(mesh, &bcs, P_INIT));

  Element* e;
  int i = 1;
  for_all_active_elements(e, mesh)
  {
    space->set_element_order(e->id, i++ % 9 + 1);
  }

  space->assign_dofs();

  // Initialize the solution.
  MeshFunctionSharedPtr<double> sln(new Solution<double>());

  // Initialize linear solver.
  LinearSolver<double> linear_solver(&wf, space);

  // Solve the linear problem.
  try
  {
    linear_solver.solve();

    // Get the solution vector.
    double* sln_vector = linear_solver.get_sln_vector();

    // Translate the solution vector into the previously initialized Solution.
    Solution<double>::vector_to_solution(sln_vector, space, sln);

    // VTK output.
    if(VTK_VISUALIZATION)
    {
      // Output solution in VTK format.
      Views::Linearizer lin;
      bool mode_3D = false;
      lin.save_solution_vtk(sln, "sln.vtk", "Temperature", mode_3D, 1, Views::HERMES_EPS_LOW);

      // Output mesh and element orders in VTK format.
      Views::Orderizer ord;
      ord.save_mesh_vtk(space, "mesh->vtk");
      ord.save_orders_vtk(space, "ord.vtk");
    }

    // Visualize the solution.
    Views::ScalarView viewS("Solution", new Views::WinGeom(50, 50, 1000, 800));
    Views::OrderView viewO("Orders", new Views::WinGeom(50, 50, 1000, 800));
    Views::VectorView viewV("Vectors", new Views::WinGeom(50, 50, 1000, 800));
    Views::MeshView viewM("Mesh", new Views::WinGeom(50, 50, 1000, 800));

    if(HERMES_VISUALIZATION)
    {
      viewS.show(sln);
      viewS.save_screenshot("000-base.bmp");
      viewS.show_contours(1.0);
      viewS.save_screenshot("001-contours1.bmp");
      viewS.show_contours(3.0);
      viewS.save_screenshot("002-contours2.bmp");
      viewS.show_contours(5.0);
      viewS.save_screenshot("003-contours3.bmp");
      viewS.set_3d_mode(true);
      viewS.save_screenshot("004-3D.bmp");
      viewS.set_3d_mode(false);
      viewS.save_screenshot("005-backTo2D.bmp");
      viewS.hide_contours();
      viewS.save_screenshot("006-noContours.bmp");
      viewS.set_palette_filter(false);
      viewS.save_screenshot("007-paletteSmooth.bmp");
      viewS.set_palette_filter(true);
      viewS.save_screenshot("008-paletteBack.bmp");
      viewS.show_mesh(false);
      viewS.save_screenshot("009-nomesh.bmp");
      viewS.set_scale_size(40, 1230, 14);
      viewS.save_screenshot("010-scaleSize.bmp");

      MeshFunctionSharedPtr<double> slnZero(new ZeroSolution<double>(mesh));
      //viewS.set_3d_mode(false);
      viewS.show_mesh(true);
      viewS.set_scale_size(30, 120, 14);
      viewS.show(slnZero);
      viewS.save_screenshot("011-zeroSolution.bmp");

      Mesh::ReferenceMeshCreator ref_mesh_creator(mesh);
      MeshSharedPtr ref_mesh = ref_mesh_creator.create_ref_mesh();
      Space<double>::ReferenceSpaceCreator ref_space_creator(space, ref_mesh);
      SpaceSharedPtr<double> ref_space = ref_space_creator.create_ref_space();
      MeshFunctionSharedPtr<double> slnZeroReference(new ZeroSolution<double>(ref_space->get_mesh()));
      viewS.show(slnZeroReference);
      viewS.save_screenshot("012-zeroSolutionRef1.bmp");
      for(int i = 0; i < 4; i++)
      {
        int number = ref_space->get_mesh()->get_max_element_id();
        for(int j = 0; j < number; j++)
          if((j % 5) < 3)
            if(ref_space->get_mesh()->get_element(j)->active)
              ref_space->get_mesh()->refine_element_id(j);
      }
      MeshFunctionSharedPtr<double> slnConstRefined(new ConstantSolution<double>(ref_space->get_mesh(), 1.234567));
      MeshFunctionSharedPtr<double> c1s(new TestExactSolution1(ref_space->get_mesh()));
      MeshFunctionSharedPtr<double> c2s(new TestExactSolution2(ref_space->get_mesh()));
      
      viewS.show(slnConstRefined, Views::HERMES_EPS_NORMAL, 1, c1s, c2s, 0.4);
      //viewS.set_3d_mode(true);
      viewS.save_screenshot("013-constSolutionRefWithDisplacement.bmp");
      viewO.show(space); 
      viewO.save_screenshot("100-space->bmp");
      viewO.set_b_orders(true);
      viewO.save_screenshot("101-spaceBOrders.bmp");
      viewM.show(mesh);
      viewM.save_screenshot("200-mesh->bmp");
      viewM.set_b_elem_mrk(true);
      viewM.save_screenshot("201-meshBElemMrk.bmp");
      viewV.show(sln, sln);
      viewV.save_screenshot("300-vectorizer.bmp");
      viewV.show(sln, sln, Views::HERMES_EPS_VERYHIGH);
      viewV.save_screenshot("301-vectorizerFiner.bmp");
      viewV.get_vectorizer()->set_curvature_epsilon(1e-5);
      viewV.show(sln, sln);
      viewV.save_screenshot("302-vectorizerFinerWithFinerCurves.bmp");
      viewV.get_vectorizer()->set_curvature_epsilon(1e-3);
      viewV.set_mode(1);
      viewV.show(sln, sln);
      viewV.save_screenshot("303-vectorizerArrowsMode.bmp");
      
      MeshFunctionSharedPtr<double> c1v(new TestExactSolution1(mesh));
      MeshFunctionSharedPtr<double> c2v(new TestExactSolution2(mesh));
      
      viewV.show(sln, sln, Views::HERMES_EPS_NORMAL, 1, 1, c1v, c2v, 0.5);
      viewV.save_screenshot("303-vectorizerArrowsModeWithDisplacement.bmp");
    }
  }
  catch(std::exception& e)
  {
    std::cout << e.what();
  }

  return 0;
}
