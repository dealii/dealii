/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2021 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * <br>
 *
 * <i>
 * This program was contributed by Peter Munch. This work and the required
 * generalizations of the internal data structures of deal.II form part of the
 * project "Virtual Materials Design" funded by the Helmholtz Association of
 * German Research Centres.
 * </i>
 *
 *
 * <a name="Intro"></a>
 * <h1>Introduction</h1>
 *
 * <h3>Motivation</h3>
 *
 * Many freely available mesh-generation tools produce meshes that consist of
 * simplices (triangles in 2d; tetrahedra in 3d). The reason for this is that
 * generating such kind of meshes for complex geometries is simpler than the
 * generation of hex-only meshes. This tutorial shows how to work on such kind
 * of meshes with the experimental simplex features in deal.II. For this
 * purpose, we solve the Poisson problem from step-3 in 2d with a mesh only
 * consisting of triangles.
 *
 *
 * <h3>Working on simplex meshes</h3>
 *
 * To be able to work on simplex meshes, one has to select appropriate finite
 * elements, quadrature rules, and mapping objects. In step-3, we used FE_Q,
 * QGauss, and (implicitly by not specifying a mapping) MappingQ1. The
 * equivalent classes for the first two classes in the context of simplices are
 * FE_SimplexP and QGaussSimplex, which we will utilize here. For mapping
 * purposes, we use the class MappingFE, which implements an isoparametric
 * mapping. We initialize it with an FE_SimplexP object so that it can be
 * applied on simplex meshes.
 *
 *
 * <h3>Mesh generation</h3>
 *
 * In contrast to step-3, we do not use a function from the GridGenerator
 * namespace, but rather read an externally generated mesh. For this tutorial,
 * we have created the mesh (square with width and height of one) with Gmsh with
 * the following journal file "box_2D_tri.geo":
 *
 * @code
 * Rectangle(1) = {0, 0, 0, 1, 1, 0};
 * Mesh 2;
 * Save "box_2D_tri.msh";
 * @endcode
 *
 * The journal file can be processed by Gmsh generating the actual mesh with the
 * ending ".geo":
 *
 * @code
 * gmsh box_2D_tri.geo
 * @endcode
 *
 * We have included in the tutorial folder both the journal file and the mesh
 * file in the event that one does not have access to Gmsh.
 *
 * The mesh can be simply read by deal.II with methods provided by the GridIn
 * class, as shown below.
 *
 */


// @sect3{Include files}

// Include files, as used in step-3:
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

// Include files that contain appropriate quadrature rules, finite elements,
// and mapping objects for simplex meshes.
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

// The following file contains the class GridIn, which allows us to read
// external meshes.
#include <deal.II/grid/grid_in.h>

using namespace dealii;

// @sect3{The <code>Step3</code> class}
//
// This is the main class of the tutorial. Since it is very similar to the
// version from step-3, we will only point out and explain the relevant
// differences that allow to perform simulations on simplex meshes.

class Step3
{
public:
  Step3();

  void run();

private:
  void make_grid();
  void setup_system();
  void assemble_system();
  void solve();
  void output_results() const;

  Triangulation<2> triangulation;

  // Here, we select a mapping object, a finite element, and a quadrature rule
  // that are compatible with simplex meshes.
  const MappingFE<2>     mapping;
  const FE_SimplexP<2>   fe;
  const QGaussSimplex<2> quadrature_formula;

  DoFHandler<2> dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};


// @sect4{Step3::Step3}
//
// In the constructor, we set the polynomial degree of the finite element and
// the number of quadrature points. Furthermore, we initialize the MappingFE
// object with a (linear) FE_SimplexP object so that it can work on simplex
// meshes.
Step3::Step3()
  : mapping(FE_SimplexP<2>(1))
  , fe(2)
  , quadrature_formula(3)
  , dof_handler(triangulation)
{}


// @sect4{Step3::make_grid}
//
// Read the external mesh file "box_2D_tri.msh" as in step-3.
void Step3::make_grid()
{
  GridIn<2>(triangulation).read("box_2D_tri.msh");

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}


// @sect4{Step3::setup_system}
//
// From here on, nothing has changed. Not even, the
// cell integrals have been changed depending on whether one operates on
// hypercube or simplex meshes. This is astonishing and is possible due to the
// design of the following two classes:
//  - DoFHandler: this class stores degrees of freedom in a flexible way and
//    allows simple access to them depending on the element type independent of
//    the cell type.
//  - FEValues: this class hides the details of finite element, quadrature rule,
//    and mapping (even if the implementations are inherently different)
//    behind a unified interface.
void Step3::setup_system()
{
  dof_handler.distribute_dofs(fe);
  std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
            << std::endl;
  DynamicSparsityPattern dsp(dof_handler.n_dofs());
  DoFTools::make_sparsity_pattern(dof_handler, dsp);
  sparsity_pattern.copy_from(dsp);

  system_matrix.reinit(sparsity_pattern);

  solution.reinit(dof_handler.n_dofs());
  system_rhs.reinit(dof_handler.n_dofs());
}


// @sect4{Step3::assemble_system}
//
// Nothing has changed here.
void Step3::assemble_system()
{
  FEValues<2> fe_values(mapping,
                        fe,
                        quadrature_formula,
                        update_values | update_gradients | update_JxW_values);

  const unsigned int dofs_per_cell = fe.n_dofs_per_cell();
  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);
  std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);

      cell_matrix = 0;
      cell_rhs    = 0;

      for (const unsigned int q_index : fe_values.quadrature_point_indices())
        {
          for (const unsigned int i : fe_values.dof_indices())
            for (const unsigned int j : fe_values.dof_indices())
              cell_matrix(i, j) +=
                (fe_values.shape_grad(i, q_index) * // grad phi_i(x_q)
                 fe_values.shape_grad(j, q_index) * // grad phi_j(x_q)
                 fe_values.JxW(q_index));           // dx

          for (const unsigned int i : fe_values.dof_indices())
            cell_rhs(i) += (fe_values.shape_value(i, q_index) * // phi_i(x_q)
                            1. *                                // f(x_q)
                            fe_values.JxW(q_index));            // dx
        }
      cell->get_dof_indices(local_dof_indices);

      for (const unsigned int i : fe_values.dof_indices())
        for (const unsigned int j : fe_values.dof_indices())
          system_matrix.add(local_dof_indices[i],
                            local_dof_indices[j],
                            cell_matrix(i, j));

      for (const unsigned int i : fe_values.dof_indices())
        system_rhs(local_dof_indices[i]) += cell_rhs(i);
    }


  std::map<types::global_dof_index, double> boundary_values;
  VectorTools::interpolate_boundary_values(
    mapping, dof_handler, 0, Functions::ZeroFunction<2>(), boundary_values);
  MatrixTools::apply_boundary_values(boundary_values,
                                     system_matrix,
                                     solution,
                                     system_rhs);
}


// @sect4{Step3::solve}
//
// Nothing has changed here.
void Step3::solve()
{
  SolverControl            solver_control(1000, 1e-12);
  SolverCG<Vector<double>> solver(solver_control);
  solver.solve(system_matrix, solution, system_rhs, PreconditionIdentity());

  std::cout << solver_control.last_step()
            << " CG iterations needed to obtain convergence." << std::endl;
}


// @sect4{Step3::output_results}
//
// Nothing has changed here.
void Step3::output_results() const
{
  DataOut<2> data_out;

  DataOutBase::VtkFlags flags;
  flags.write_higher_order_cells = true;
  data_out.set_flags(flags);

  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(solution, "solution");
  data_out.build_patches(mapping, 2);
  std::ofstream output("solution.vtk");
  data_out.write_vtk(output);
}


// @sect4{Step3::run}
//
// Nothing has changed here.
void Step3::run()
{
  make_grid();
  setup_system();
  assemble_system();
  solve();
  output_results();
}


// @sect3{The <code>main</code> function}
//
// Nothing has changed here.
int main()
{
  Step3 laplace_problem;
  laplace_problem.run();

  return 0;
}

/*
 * <h1>Results</h1>
 *
 * The following figures show the mesh and the result obtained by executing this
 * program:
 *
 * <table align="center" class="doxtable" style="width:65%">
 *   <tr>
 *     <td>
 *         @image html step_3_simplex_0.png
 *     </td>
 *     <td>
 *         @image html step_3_simplex_1.png
 *     </td>
 *   </tr>
 * </table>
 *
 * Not surprisingly, the result looks as expected.
 *
 *
 * <h3>Possibilities for extensions</h3>
 *
 * In this tutorial, we presented how to use the deal.II simplex infrastructure
 * to solve a simple Poisson problem on a simplex mesh in 2d. In this scope, we
 * could only present a small section of the capabilities. In the following, we
 * point out further capabilities briefly.
 *
 *
 * <h4>3d meshes and codim-1 meshes in 3d</h4>
 *
 * An extension to 3d is quite straightforward. Both FE_SimplexP and
 * QGaussSimplex are implemented in a dimensional-independent way so that simply
 * replacing everywhere dim=2 with dim=3 should work out of the box.
 *
 * Furthermore, embedding of a 2d mesh consisting of triangles in 3d space is
 * possible.
 *
 *
 * <h4>Mixed meshes</h4>
 *
 * In step-3, we considered meshes only consisting of quadrilaterals. In this
 * tutorial, we took a look at the case that the mesh only consists of
 * triangles. In the general case (also known as mixed mesh), the mesh consists
 * of both cell types. In 3d, meshes might even consist of more cell types, like
 * wedges/prisms and pyramids. We consider such meshes in the tutorial
 * step-3mixed.
 *
 *
 * <h4>Alternative finite elements, quadrature rules, and mapping objects</h4>
 *
 * In this tutorial, we used the most basic finite-element, quadrature-rule, and
 * mapping classes. However, more classes are compatible with simplices. The
 * following list gives an overview of these classes:
 * - finite elements: FE_SimplexP, FE_SimplexDGP, FE_SimplexP_Bubbles
 * - quadrature rules: QGaussSimplex, QWitherdenVincentSimplex, QDuffy
 * - mapping objects: MappingFE, MappingFEField
 *
 * It should be also pointed out that FESystems can also handle simplex finite
 * elements which is crucial to solve vector-valued problems, as needed, e.g.,
 * to solve elasticity and fluid problems (see also step-17).
 *
 *
 * <h4>Alternative mesh generation approaches</h4>
 *
 * In this tutorial, we have created the mesh externally and read it with the
 * help of GridIn. Since we believe that the main motivation to work on simplex
 * meshes is that one has a complex geometry that can only be meshed with
 * an external tool with simplices, deal.II does not have too many functions in
 * the GridGenerator namespace, targeting simplex meshes. However, we would like
 * to point out the following functions:
 *  - GridGenerator::subdivided_hyper_cube_with_simplices() and
 *    GridGenerator::subdivided_hyper_rectangle_with_simplices(), which fill a
 *    hypercube and a hyperrectangle domain with simplices
 *  - GridGenerator::convert_hypercube_to_simplex_mesh(), which converts meshes
 *    consisting of hypercube cells to simplex meshes by replacing one
 *    quadrilateral with 4 triangles and one hexahedron with 24 tetrahedrons
 *
 *
 * <h4>hp-adaptivity</h4>
 *
 * Here, we considered a mesh without refinements and with all cells assigned
 * the same type of element with the same polynomial degree. However, one is not
 * restricted to this. For further details on hp-methods, see step-27.
 *
 *
 * <h4>Parallelization</h4>
 *
 * To parallelize the code, one needs to replace the Triangulation object either
 * with parallel::shared::Triangulation or
 * parallel::fullydistributed::Triangulation and make some minor adjustments, as
 * discussed in step-6.
 *
 *
 * <h4>Face integrals and discontinuous Galerkin methods</h4>
 *
 * The classes FEFaceValues and FEInterfaceValues are also compatible with
 * simplex meshes.
 *
 *
 * <h4>Matrix-free operator evaluation</h4>
 *
 * In this tutorial, we showed a matrix-based approach. However, one could also
 * rewrite the code using MatrixFree, FEEvaluation, and FEFaceEvaluation, which
 * are also compatible with simplex meshes.
 *
 */
