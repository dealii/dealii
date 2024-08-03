/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2023 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 */

// Make sure that we get the correct convergence order for a combination of P3
// and Q3 elements and and two reversed quadrilateral line orientations.

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_simplex_p.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_fe.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools_topology.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools_boundary.h>
#include <deal.II/numerics/vector_tools_integrate_difference.h>
#include <deal.II/numerics/vector_tools_interpolate.h>

#include "../tests.h"

// u
class Exact : public Function<2>
{
public:
  double
  value(const Point<2> &p, unsigned int /*component*/ = 0) const
  {
    return std::sin(p[0]) * std::sin(p[1]);
  }
};

// -Laplacian(u) + u
class Forcing : public Function<2>
{
public:
  double
  value(const Point<2> &p, unsigned int /*component*/ = 0) const
  {
    return 3.0 * std::sin(p[0]) * std::sin(p[1]);
  }
};

int
main()
{
  initlog();

  const std::vector<Point<2>> vertices{{0.0, 0.0},
                                       {1.0, 0.0},
                                       {0.0, 2.0},
                                       {3.0, 1.0},
                                       {3.0, 3.0},
                                       {2.0, 3.0},
                                       {0.0, 3.0},
                                       {3.0, 0.0}};

  // faces 1 and 3 of the quadrilateral will have reversed orientations
  std::vector<CellData<2>> cells(5);
  cells[0].vertices = {0u, 1u, 2u};
  cells[1].vertices = {3u, 4u, 5u};
  cells[2].vertices = {1u, 3u, 7u};
  cells[3].vertices = {5u, 2u, 6u};
  cells[4].vertices = {1u, 3u, 2u, 5u};

  const hp::MappingCollection<2> mapping(MappingFE<2>(FE_SimplexP<2>(1)),
                                         MappingFE<2>(FE_Q<2>(1)));

  const hp::FECollection<2> fe(FE_SimplexP<2>(3),
                               FE_Q<2>(QIterated<1>(QTrapezoid<1>(), 3)));
  const hp::QCollection<2>  quadrature_formula(QGaussSimplex<2>(4),
                                              QGauss<2>(4));

  Triangulation<2> tria;
  GridTools::invert_cells_with_negative_measure(vertices, cells);
  tria.create_triangulation(vertices, cells, SubCellData());

  for (unsigned int r = 0; r < 5; ++r)
    {
      tria.refine_global(1);
      DoFHandler<2> dof_handler(tria);
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          if (cell->reference_cell() == ReferenceCells::Triangle)
            cell->set_active_fe_index(0);
          else if (cell->reference_cell() == ReferenceCells::Quadrilateral)
            cell->set_active_fe_index(1);
        }
      dof_handler.distribute_dofs(fe);

      // set up matrices and vectors:
      AffineConstraints<double> constraints(dof_handler.locally_owned_dofs(),
                                            dof_handler.locally_owned_dofs());
      VectorTools::interpolate_boundary_values(
        mapping, dof_handler, 0, Exact(), constraints);
      constraints.close();
      SparsityPattern        sparsity_pattern;
      DynamicSparsityPattern dsp(dof_handler.n_dofs());
      DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints, false);
      sparsity_pattern.copy_from(dsp);

      SparseMatrix<double> system_matrix(sparsity_pattern);
      Vector<double>       system_rhs(dof_handler.n_dofs());
      Vector<double>       solution(dof_handler.n_dofs());

      hp::FEValues<2> hp_fe_values(mapping,
                                   fe,
                                   quadrature_formula,
                                   update_values | update_gradients |
                                     update_quadrature_points |
                                     update_JxW_values);

      FullMatrix<double>                   cell_matrix;
      Vector<double>                       cell_rhs;
      std::vector<types::global_dof_index> local_dof_indices;

      Forcing forcing;
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          hp_fe_values.reinit(cell);
          const auto &fe_values = hp_fe_values.get_present_fe_values();

          const unsigned int dofs_per_cell = cell->get_fe().n_dofs_per_cell();
          cell_matrix.reinit(dofs_per_cell, dofs_per_cell);
          cell_rhs.reinit(dofs_per_cell);
          local_dof_indices.resize(dofs_per_cell);
          cell->get_dof_indices(local_dof_indices);

          cell_matrix = 0;
          cell_rhs    = 0;

          for (const unsigned int q_index :
               fe_values.quadrature_point_indices())
            {
              for (const unsigned int i : fe_values.dof_indices())
                for (const unsigned int j : fe_values.dof_indices())
                  cell_matrix(i, j) += (fe_values.shape_grad(i, q_index) *
                                          fe_values.shape_grad(j, q_index) +
                                        fe_values.shape_value(i, q_index) *
                                          fe_values.shape_value(j, q_index)) *
                                       fe_values.JxW(q_index);

              for (const unsigned int i : fe_values.dof_indices())
                cell_rhs(i) +=
                  fe_values.shape_value(i, q_index) *
                  forcing.value(fe_values.quadrature_point(q_index)) *
                  fe_values.JxW(q_index);
            }

          constraints.distribute_local_to_global(cell_matrix,
                                                 cell_rhs,
                                                 local_dof_indices,
                                                 system_matrix,
                                                 system_rhs);
        }

      SolverControl            solver_control(1000, 1e-12);
      SolverCG<Vector<double>> solver(solver_control);

      // Cut down on a lot of solver iterations by interpolating the solution as
      // the initial guess:
      VectorTools::interpolate(mapping, dof_handler, Exact(), solution);

      PreconditionSSOR<> precondition;
      precondition.initialize(system_matrix,
                              PreconditionSSOR<>::AdditionalData());
      deallog << "rhs l2 = " << system_rhs.l2_norm() << std::endl;
      solver.solve(system_matrix, solution, system_rhs, precondition);
      constraints.distribute(solution);

      Vector<double> cell_errors(tria.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        solution,
                                        Exact(),
                                        cell_errors,
                                        quadrature_formula,
                                        VectorTools::L2_norm);

      deallog << "global error = "
              << VectorTools::compute_global_error(tria,
                                                   cell_errors,
                                                   VectorTools::L2_norm)
              << std::endl;

      std::ofstream out("out-" + std::to_string(r) + ".vtu");
      DataOut<2>    data_out;
      data_out.attach_dof_handler(dof_handler);
      data_out.add_data_vector(solution, "solution");
      data_out.add_data_vector(cell_errors, "cell_errors");
      data_out.build_patches(mapping);
      data_out.write_vtu(out);
    }
}
