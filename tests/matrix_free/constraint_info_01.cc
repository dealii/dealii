// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test internal::MatrixFreeFunctions::ConstraintInfo and compare the
// result with FEEvaluation.

#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/matrix_free/constraint_info.h>
#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"


int
main()
{
  initlog();

  const int dim    = 2;
  const int degree = 1;

  using Number              = double;
  using VectorizedArrayType = VectorizedArray<Number>;
  using VectorType          = LinearAlgebra::distributed::Vector<Number>;

  Triangulation<dim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 2);
  tria.begin()->set_refine_flag();
  tria.execute_coarsening_and_refinement();

  QGauss<dim>    quad(degree + 1);
  FE_Q<dim>      fe(degree);
  MappingQ1<dim> mapping;

  DoFHandler<dim> dof_handler;
  dof_handler.reinit(tria);
  dof_handler.distribute_dofs(fe);

  AffineConstraints<Number> constraints;

  DoFTools::make_zero_boundary_constraints(dof_handler, constraints);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  typename MatrixFree<dim, Number, VectorizedArrayType>::AdditionalData
    additional_data;
  additional_data.mapping_update_flags = update_values | update_gradients |
                                         update_JxW_values |
                                         dealii::update_quadrature_points;

  MatrixFree<dim, Number, VectorizedArrayType> matrix_free;
  matrix_free.reinit(mapping, dof_handler, constraints, quad, additional_data);

  internal::MatrixFreeFunctions::ConstraintInfo<dim, VectorizedArrayType>
    constraint_info;
  constraint_info.reinit(dof_handler, matrix_free.n_physical_cells(), true);

  for (unsigned int cell = 0, cell_no = 0; cell < matrix_free.n_cell_batches();
       ++cell)
    for (unsigned int v = 0;
         v < matrix_free.n_active_entries_per_cell_batch(cell);
         ++v, ++cell_no)
      constraint_info.read_dof_indices(cell_no,
                                       -1,
                                       matrix_free.get_cell_iterator(cell, v),
                                       constraints,
                                       matrix_free.get_vector_partitioner());

  constraint_info.finalize();

  const auto print_stat = [](const auto &dof_info) {
    for (const auto i : dof_info.dof_indices)
      deallog << i << " ";
    deallog << std::endl;

    for (const auto i : dof_info.constraint_indicator)
      deallog << "(" << i.first << "," << i.second << ") ";
    deallog << std::endl;

    // make sure to filter out contributions that only get filled up to make
    // all SIMD lanes populated
    for (unsigned int i = 0; i < dof_info.row_starts.size(); ++i)
      if (i == 0 ||
          (dof_info.row_starts[i].first > dof_info.row_starts[i - 1].first))
        deallog << "(" << dof_info.row_starts[i].first << ","
                << dof_info.row_starts[i].second << ") ";
    deallog << std::endl;

    deallog << std::endl;
  };

  print_stat(constraint_info);
  print_stat(matrix_free.get_dof_info());

  VectorType src, dst, dst_mf;

  const auto initialize_dof_vectors = [&]() {
    matrix_free.initialize_dof_vector(src);
    matrix_free.initialize_dof_vector(dst);
    matrix_free.initialize_dof_vector(dst_mf);

    for (const auto i : src.locally_owned_elements())
      src[i] = i;
  };

  initialize_dof_vectors();

  AlignedVector<VectorizedArrayType> local_vector(fe.n_dofs_per_cell());

  const auto read_dof_values = [&](const unsigned int first_cell,
                                   const unsigned int n_cells,
                                   const unsigned int n_dofs_per_cell) {
    internal::VectorReader<Number, VectorizedArrayType> reader;
    constraint_info.read_write_operation(reader,
                                         src,
                                         local_vector.data(),
                                         first_cell,
                                         n_cells,
                                         n_dofs_per_cell,
                                         true);
    constraint_info.apply_hanging_node_constraints(first_cell,
                                                   n_cells,
                                                   false,
                                                   local_vector);
  };

  const auto distribute_local_to_global =
    [&](const unsigned int first_cell,
        const unsigned int n_cells,
        const unsigned int n_dofs_per_cell) {
      internal::VectorDistributorLocalToGlobal<Number, VectorizedArrayType>
        writer;
      constraint_info.apply_hanging_node_constraints(first_cell,
                                                     n_cells,
                                                     true,
                                                     local_vector);
      constraint_info.read_write_operation(writer,
                                           dst,
                                           local_vector.data(),
                                           first_cell,
                                           n_cells,
                                           n_dofs_per_cell,
                                           true);
    };

  const auto read_dof_values_plain = [&](const unsigned int first_cell,
                                         const unsigned int n_cells,
                                         const unsigned int n_dofs_per_cell) {
    internal::VectorReader<Number, VectorizedArrayType> reader;
    constraint_info.read_write_operation(reader,
                                         src,
                                         local_vector.data(),
                                         first_cell,
                                         n_cells,
                                         n_dofs_per_cell,
                                         false);
  };

  const auto set_dof_values_plain = [&](const unsigned int first_cell,
                                        const unsigned int n_cells,
                                        const unsigned int n_dofs_per_cell) {
    internal::VectorSetter<Number, VectorizedArrayType> writer;
    constraint_info.read_write_operation(writer,
                                         dst,
                                         local_vector.data(),
                                         first_cell,
                                         n_cells,
                                         n_dofs_per_cell,
                                         false);
  };

  unsigned int       first_cell      = 0;
  const unsigned int n_dofs_per_cell = fe.n_dofs_per_cell();

  FEEvaluation<dim, -1, 0, 1, Number, VectorizedArrayType> fe_eval(matrix_free);

  for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
    {
      {
        const unsigned int n_cells =
          matrix_free.n_active_entries_per_cell_batch(cell);

        read_dof_values(first_cell, n_cells, n_dofs_per_cell);
        distribute_local_to_global(first_cell, n_cells, n_dofs_per_cell);

        first_cell += n_cells;
      }

      {
        fe_eval.reinit(cell);

        fe_eval.read_dof_values(src);
        fe_eval.distribute_local_to_global(dst_mf);
      }
    }

  dst.print(deallog.get_file_stream());
  dst_mf.print(deallog.get_file_stream());
  deallog << std::endl;

  initialize_dof_vectors();
  first_cell = 0;

  for (unsigned int cell = 0; cell < matrix_free.n_cell_batches(); ++cell)
    {
      {
        const unsigned int n_cells =
          matrix_free.n_active_entries_per_cell_batch(cell);

        read_dof_values_plain(first_cell, n_cells, n_dofs_per_cell);
        set_dof_values_plain(first_cell, n_cells, n_dofs_per_cell);

        first_cell += n_cells;
      }

      {
        fe_eval.reinit(cell);

        fe_eval.read_dof_values_plain(src);
        fe_eval.set_dof_values_plain(dst_mf);
      }
    }

  dst.print(deallog.get_file_stream());
  dst_mf.print(deallog.get_file_stream());
}
