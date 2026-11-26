// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------



// Check that FEEvaluation::read_dof_values() reads the unknowns of the right
// cell in every lane of a cell batch, and that a matrix-free mass operator
// equals the action of the assembled cell matrices, when MatrixFree is set up
// with face data and the unknowns of each cell are contiguous but are not
// disjoint blocks numbered in the order of the cells of the first
// DoFHandler. This happens in 1d for FE_Q, whose neighboring cells share an
// unknown, and for two DG DoFHandlers that number the cells in different
// orders. The non-buffering shared-memory indices of MatrixFree assumed such
// blocks, and lanes 1, ..., n_lanes-1 then read the unknowns of another cell.
// The results must not depend on the SIMD width.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <algorithm>
#include <map>

#include "../tests.h"


using VectorType = LinearAlgebra::distributed::Vector<double>;



// Read a vector that holds its own indices on every cell batch and compare,
// for each lane, the values read with the unknowns of that lane's cell.
template <int fe_degree>
void
check_read(const MatrixFree<1, double> &matrix_free, const unsigned int dof_no)
{
  VectorType x;
  matrix_free.initialize_dof_vector(x, dof_no);
  for (unsigned int i = 0; i < x.size(); ++i)
    x(i) = i;

  FEEvaluation<1, fe_degree>           phi(matrix_free, dof_no);
  std::vector<types::global_dof_index> dof_indices(phi.dofs_per_cell);

  // collect by active cell index, so that the output does not depend on the
  // composition of the cell batches
  std::map<unsigned int, std::pair<std::vector<double>, std::vector<double>>>
    cells;
  for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
    {
      phi.reinit(cell);
      phi.read_dof_values(x);
      for (unsigned int v = 0;
           v < matrix_free.n_active_entries_per_cell_batch(cell);
           ++v)
        {
          const auto cell_iterator =
            matrix_free.get_cell_iterator(cell, v, dof_no);
          cell_iterator->get_dof_indices(dof_indices);
          auto &entry = cells[cell_iterator->active_cell_index()];
          for (unsigned int i = 0; i < phi.dofs_per_cell; ++i)
            {
              entry.first.push_back(dof_indices[i]);
              entry.second.push_back(phi.get_dof_value(i)[v]);
            }
          // FEEvaluation uses the lexicographic order of the unknowns,
          // get_dof_indices() the hierarchical one
          std::sort(entry.first.begin(), entry.first.end());
          std::sort(entry.second.begin(), entry.second.end());
        }
    }

  unsigned int n_correct = 0;
  for (const auto &[index, entry] : cells)
    if (entry.first == entry.second)
      ++n_correct;
    else
      {
        deallog << "cell " << index << " has unknowns";
        for (const double d : entry.first)
          deallog << ' ' << d;
        deallog << " but reads";
        for (const double d : entry.second)
          deallog << ' ' << d;
        deallog << std::endl;
      }
  deallog << "DoFHandler " << dof_no << ": " << n_correct << " of "
          << cells.size() << " cells read their own unknowns" << std::endl;
}



// Compare a matrix-free mass operator with the action of the assembled cell
// matrices.
template <int fe_degree>
void
check_mass(const MatrixFree<1, double> &matrix_free, const unsigned int dof_no)
{
  const DoFHandler<1> &dof_handler = matrix_free.get_dof_handler(dof_no);

  VectorType x, y, y_ref;
  matrix_free.initialize_dof_vector(x, dof_no);
  matrix_free.initialize_dof_vector(y, dof_no);
  matrix_free.initialize_dof_vector(y_ref, dof_no);
  for (unsigned int i = 0; i < x.size(); ++i)
    x(i) = 1. + 0.25 * i;

  matrix_free.cell_loop<VectorType, VectorType>(
    [&](const MatrixFree<1, double>                 &data,
        VectorType                                  &dst,
        const VectorType                            &src,
        const std::pair<unsigned int, unsigned int> &range) {
      FEEvaluation<1, fe_degree> phi(data, dof_no);
      for (unsigned int cell = range.first; cell < range.second; ++cell)
        {
          phi.reinit(cell);
          phi.read_dof_values(src);
          phi.evaluate(EvaluationFlags::values);
          for (const unsigned int q : phi.quadrature_point_indices())
            phi.submit_value(phi.get_value(q), q);
          phi.integrate(EvaluationFlags::values);
          phi.distribute_local_to_global(dst);
        }
    },
    y,
    x,
    true);

  const FiniteElement<1> &fe = dof_handler.get_fe();
  const MappingQ1<1>      mapping;
  FEValues<1>             fe_values(mapping,
                                    fe,
                                    QGauss<1>(fe_degree + 1),
                                    update_values | update_JxW_values);
  std::vector<types::global_dof_index> dof_indices(fe.n_dofs_per_cell());
  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      fe_values.reinit(cell);
      cell->get_dof_indices(dof_indices);
      for (unsigned int i = 0; i < fe.n_dofs_per_cell(); ++i)
        for (unsigned int j = 0; j < fe.n_dofs_per_cell(); ++j)
          for (const unsigned int q : fe_values.quadrature_point_indices())
            y_ref(dof_indices[i]) += fe_values.shape_value(i, q) *
                                     fe_values.shape_value(j, q) *
                                     fe_values.JxW(q) * x(dof_indices[j]);
    }

  const double norm = y_ref.linfty_norm();
  y -= y_ref;
  const double error = y.linfty_norm() / norm;
  deallog << "DoFHandler " << dof_no << ": mass operator relative error ";
  if (error < 1e-12)
    deallog << "below 1e-12" << std::endl;
  else
    deallog << error << std::endl;
}



// Report whether the non-buffering shared-memory indices are set up as
// expected: only for disjoint blocks of unknowns, and only if the cell batches
// are stored as 'contiguous', which requires more than one lane.
void
check_shared_memory_indices(const MatrixFree<1, double> &matrix_free,
                            const bool                   disjoint_blocks)
{
  bool expected = false;
#ifdef DEAL_II_WITH_MPI
  expected = disjoint_blocks && VectorizedArray<double>::size() > 1;
#else
  (void)disjoint_blocks;
#endif
  const bool set_up =
    matrix_free.get_dof_info(0).dof_indices_contiguous_sm[2].size() > 0;
  deallog << "shared-memory cell indices "
          << (set_up == expected ? "as expected" : "NOT as expected")
          << std::endl;
}



MatrixFree<1, double>::AdditionalData
additional_data(const bool inner_faces)
{
  MatrixFree<1, double>::AdditionalData data;
  data.mapping_update_flags                = update_values | update_JxW_values;
  data.mapping_update_flags_boundary_faces = update_JxW_values;
  if (inner_faces)
    data.mapping_update_flags_inner_faces = update_JxW_values;
  return data;
}



// FE_Q(1) and FE_Q(2): neighboring cells share an unknown. FE_Q(2) is
// renumbered along x, so that the unknowns of each cell are contiguous.
template <int fe_degree>
void
test_fe_q(const bool inner_faces)
{
  deallog.push("FE_Q" + std::to_string(fe_degree) +
               (inner_faces ? "_all_faces" : "_boundary_faces"));

  Triangulation<1> tria;
  GridGenerator::subdivided_hyper_cube(tria, 7);
  FE_Q<1>       fe(fe_degree);
  DoFHandler<1> dof_handler(tria);
  dof_handler.distribute_dofs(fe);
  if (fe_degree > 1)
    {
      Tensor<1, 1> direction;
      direction[0] = 1.;
      DoFRenumbering::downstream(dof_handler, direction, true);
    }

  AffineConstraints<double> constraints;
  constraints.close();

  const MappingQ1<1>    mapping;
  MatrixFree<1, double> matrix_free;
  matrix_free.reinit(mapping,
                     dof_handler,
                     constraints,
                     QGauss<1>(fe_degree + 1),
                     additional_data(inner_faces));

  check_shared_memory_indices(matrix_free, false);
  check_read<fe_degree>(matrix_free, 0);
  check_mass<fe_degree>(matrix_free, 0);

  deallog.pop();
}



// FE_DGQ(1) with one DoFHandler: disjoint blocks in the order of the cells,
// for which the shared-memory indices remain in use. With a second DoFHandler
// whose cells are numbered in the reverse order, they must not be used.
void
test_fe_dgq(const bool second_handler_reversed)
{
  deallog.push(second_handler_reversed ? "FE_DGQ1_two_handlers" :
                                         "FE_DGQ1_one_handler");

  Triangulation<1> tria;
  GridGenerator::subdivided_hyper_cube(tria, 7);
  FE_DGQ<1>     fe(1);
  DoFHandler<1> dof_handler_0(tria);
  dof_handler_0.distribute_dofs(fe);
  DoFHandler<1> dof_handler_1(tria);
  dof_handler_1.distribute_dofs(fe);
  std::vector<DoFHandler<1>::active_cell_iterator> reversed;
  for (const auto &cell : dof_handler_1.active_cell_iterators())
    reversed.push_back(cell);
  std::reverse(reversed.begin(), reversed.end());
  DoFRenumbering::cell_wise(dof_handler_1, reversed);

  AffineConstraints<double> constraints;
  constraints.close();

  std::vector<const DoFHandler<1> *> dof_handlers{&dof_handler_0};
  if (second_handler_reversed)
    dof_handlers.push_back(&dof_handler_1);
  const std::vector<const AffineConstraints<double> *> all_constraints(
    dof_handlers.size(), &constraints);

  const MappingQ1<1>    mapping;
  MatrixFree<1, double> matrix_free;
  matrix_free.reinit(mapping,
                     dof_handlers,
                     all_constraints,
                     QGauss<1>(2),
                     additional_data(true));

  check_shared_memory_indices(matrix_free, !second_handler_reversed);
  for (unsigned int dof_no = 0; dof_no < dof_handlers.size(); ++dof_no)
    {
      check_read<1>(matrix_free, dof_no);
      check_mass<1>(matrix_free, dof_no);
    }

  deallog.pop();
}



int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_init(argc, argv, 1);
  initlog();

  test_fe_q<1>(false);
  test_fe_q<1>(true);
  test_fe_q<2>(false);
  test_fe_dgq(false);
  test_fe_dgq(true);
}
