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
 * The motivation for using simplex meshes (as done in step-3simplex) is
 * straightforward: many freely available mesh-generation tools are very good in
 * creating good-quality meshes in such a format, while they struggle with
 * hex-only meshes. Hex-only meshes, on the other hand, are characterized with
 * better numerical properties (e.g., less degrees of freedom for the same
 * degree of accuracy and possibly better performance since the tensor-product
 * structure can be exploited) and are, as a consequence, the natural choice for
 * rather simple geometries and for meshes described by a coarse mesh with a few
 * cells (like a hyperball) and obtained in their final form through iterative
 * local refinement.
 *
 * Mixed meshes try to combine the best of both worlds by partitioning the
 * geometry in parts that can be easily meshed by hypercube cells
 * (quadrilaterals in 2d, hexahedrons in 3d) and in parts that can not be meshed
 * easily, requiring simplices (triangles in 2d, tetrahedrons in 3d). Since one
 * assumes that the region requiring simplices is rather small compared to the
 * rest of the domain where more efficient and accurate methods can be applied,
 * one can expect that the overall efficiency is hardly impacted by such an
 * approach.
 *
 * One should note that in 3d, one also needs a transition region between
 * hypercube and simplex regions. Here, one can use wedges/prisms and/or
 * pyramids.
 *
 *
 * <h3>Working with mixed meshes</h3>
 *
 * <i>
 * In the following, we concentrate, for the sake of simplicity, on 2d meshes:
 * they can only contain triangles and quadrilaterals. However, as detailed in
 * the outlook, an extension of the presented approach to the 3d case is
 * straightforward.
 * </i>
 *
 * The complexity of working with mixed meshes in 2d results from the fact
 * that it contains of two
 * types of geometrical objects: quadrilaterals and triangles. How to deal with
 * quadrilaterals, we have discussed in step-3: we selected an appropriate
 * finite element, quadrature rule and mapping object, e.g., FE_Q, QGauss, and
 * MappingFE (initialized with FE_Q). For simplex meshes, we selected in
 * step-3simplex FE_SimplexP, QGaussSimplex, and MappingFE (initialized with
 * FE_SimplexP).
 *
 * For mixed meshes, we need multiple finite elements, quadrature rules, and
 * mapping objects (one set for triangles and one set for quadrilaterals) in the
 * same program. To ease the work with the multitude of objects (in particular
 * in 3d, we need at least four of each), you can collect the objects and group
 * them together in hp::FECollection, hp::QCollection, and
 * hp::MappingCollection.
 *
 * Just like in the context of finite elements, quadrature rules, and mapping
 * objects, we need multiple FEValues objects: the collection of FEValues is
 * called hp::FEValues. It returns the FEValues object needed for the current
 * cell via the method hp::FEValues::get_present_fe_values().
 *
 * For hp::FEValues, to be able to select the right finite element/quadrature
 * rule/mapping object set, it queries the active FE index of the given cell
 * during hp::FEValues::reinit(). The indices have to be set - as shown below -
 * before calling DoFHandler::distribute_dofs() by the user.
 *
 * <i>
 * The namespace name of hp::FECollection, hp::QCollection,
 * hp::MappingCollection, and hp::FEValues indicates that these classes have not
 * been written for mixed meshes in the first place, but for problems where each
 * (hypercube) cell could have a different type of finite element assigned - in
 * the simplest case, all cells have the same element type but different
 * polynomial degrees p (the reason for the letter "p" in "hp"). An extension of
 * this infrastructure to work not only on different element types but also on
 * different geometrical objects was a natural choice. For further details on
 * hp-methods, see step-27.
 * </i>
 *
 * <h3>Mesh generation</h3>
 *
 * Just like in step-3simplex, we read an externally generated mesh. For this
 * tutorial, we have created the mesh (square with width and height of one;
 * quadrilaterals on the left half and triangles on the right half) with Gmsh
 * with the following journal file "box_2D_mixed.geo":
 *
 * @code
 * Rectangle(1) = {0.0, 0, 0, 0.5, 1, 0};
 * Rectangle(2) = {0.5, 0, 0, 0.5, 1, 0};
 * Recombine Surface{1};
 * Physical Surface("All") = {1, 2};
 * Mesh 2;
 * Coherence Mesh;
 * Save "box_2D_mixed.msh";
 * @endcode
 *
 * The journal file can be processed by Gmsh generating the actual mesh with the
 * ending ".msh":
 *
 * @code
 * gmsh box_2D_mixed.geo
 * @endcode
 *
 * We have included in the tutorial folder both the journal file and the mesh
 * file in the event that one does not have access to Gmsh.
 *
 */


// @sect3{Include files}

// Include files, as used in step-3:
#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
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

// Include files, as added in step-3simplex:
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_in.h>

// Include files that we need in this tutorial to be able to deal with
// collections of finite elements, quadrature rules, mapping objects, and
// FEValues.
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

using namespace dealii;

// @sect3{The <code>Step3</code> class}
//
// This is the main class of the tutorial. Since it is very similar to the
// version from step-3 and step-3simplex, we will only point out and explain
// the relevant differences that allow to perform simulations on mixed meshes.

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

  // As already explained, we are not working with mapping objects, finite
  // elements, and quadrature rules directly but with collections of them.
  const hp::MappingCollection<2> mapping;
  const hp::FECollection<2>      fe;
  const hp::QCollection<2>       quadrature_formula;

  DoFHandler<2> dof_handler;

  SparsityPattern      sparsity_pattern;
  SparseMatrix<double> system_matrix;

  Vector<double> solution;
  Vector<double> system_rhs;
};


// @sect4{Step3::Step3}
//
// In the constructor of the Step3 class, we fill the collections. Here, we
// position the objects related to triangles in the first place (index 0) and
// the ones related to quadrilaterals in the second place (index 1).
Step3::Step3()
  : mapping(MappingFE<2>(FE_SimplexP<2>(1)), MappingFE<2>(FE_Q<2>(1)))
  , fe(FE_SimplexP<2>(2), FE_Q<2>(2))
  , quadrature_formula(QGaussSimplex<2>(3), QGauss<2>(3))
  , dof_handler(triangulation)
{}


// @sect4{Step3::make_grid}
//
// Read the external mesh file "box_2D_mixed.msh" as in step-3simplex.
void Step3::make_grid()
{
  GridIn<2>(triangulation).read("box_2D_mixed.msh");

  std::cout << "Number of active cells: " << triangulation.n_active_cells()
            << std::endl;
}


// @sect4{Step3::setup_system}
//
// In contrast to step-3 and step-3simplex, we need here a preprocessing step
// that assigns an active FE index to each cell consistently according to the
// indices in the collections and the cell type.
void Step3::setup_system()
{
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      if (cell->reference_cell() == ReferenceCells::Triangle)
        cell->set_active_fe_index(0);
      else if (cell->reference_cell() == ReferenceCells::Quadrilateral)
        cell->set_active_fe_index(1);
      else
        DEAL_II_NOT_IMPLEMENTED();
    }

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
// The following function looks similar to the version in step-3 and
// step-3simplex with the following two differences:
//  - We do not work with FEValues directly but with the collection class
//    hp::FEValues. It gives us - after it has been initialized with the current
//    cell - a reference to the right FEValues object (constructed
//    with the correct mapping object, finite element, and quadrature rule),
//    which can be used as usual to compute the cell integrals.
//  - The cell-local @ref GlossStiffnessMatrix "stiffness matrix" and the right-hand-side vector have
//    different sizes depending on the cell type (6 DoFs vs. 9 DoFs) so that
//    they might need to be resized for each cell.
//
// Apart from these two changes, the code has not changed. Not even, the
// cell integrals have been changed depending on whether one operates on
// hypercube, simplex, or mixed meshes.
void Step3::assemble_system()
{
  hp::FEValues<2> hp_fe_values(mapping,
                               fe,
                               quadrature_formula,
                               update_values | update_gradients |
                                 update_JxW_values);

  FullMatrix<double>                   cell_matrix;
  Vector<double>                       cell_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      hp_fe_values.reinit(cell);

      const auto &fe_values = hp_fe_values.get_present_fe_values();

      const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
      cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
      cell_rhs.reinit(dofs_per_cell);
      local_dof_indices.resize(dofs_per_cell);

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

/**
 * <h1>Results</h1>
 *
 * The following figures show the mesh and the result obtained by executing this
 * program:
 *
 * <table align="center" class="doxtable" style="width:65%">
 *   <tr>
 *     <td>
 *         @image html step_3_mixed_0.png
 *     </td>
 *     <td>
 *         @image html step_3_mixed_1.png
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
 * to solve a simple Poisson problem on a mixed mesh in 2d. In this scope, we
 * could only present a small section of the capabilities. In the following, we
 * point out further capabilities briefly.
 *
 *
 * <h4>Pure hypercube and simplex meshes</h4>
 *
 * In this tutorial, we worked on a mesh consisting both of quadrilaterals and
 * triangles. However, the program and the underlying concepts also work if the
 * mesh only contains either quadrilaterals or triangles. Interested users can
 * try this out: we have provided appropriate journal files and meshes for such
 * cases.
 *
 *
 * <h4>3d meshes</h4>
 *
 * In 3d, meshes might also consist of wedges/prisms and pyramids. Therefore,
 * the above introduced collections might consist of four components.
 *
 * For wedge/prism and pyramid cell types, following finite-element and
 * quadrature-rule classes are available:
 *  - wedge: FE_WedgeP, FE_WedgeDGP, QGaussWedge, MappingFE
 *  - pyramid: FE_PyramidP, FE_PyramidDGP, QGaussPyramid, MappingFE
 *
 * <h4>Parallelization, face integrals, discontinuous Galerkin methods, and
 * matrix-free operator evaluation</h4>
 *
 * Regarding these aspects, the same comments are valid that are described in
 * step-3simplex for pure simplex meshes.
 *
 */
