// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

// check FEEvaluation with manual set/get via RowsBlockAccessor
// same setup as in bcsr_13.cc
// for BCSR sparse vectors. To that end,
// 1) reorder DoFs based on support point locations
// 2) use 1 cell and quadratic FE
// 3) block by number of DoFs in one direction == 3
// 4) add couple of columns

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_csr_matrix.h>

#include <deal.II/matrix_free/fe_evaluation.h>
#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "bcsr_helper.h"

using namespace dealii;



template <int dim,
          int fe_degree,
          int n_q_points_1d = fe_degree + 1,
          typename Number   = double,
          int n_components  = 1>
class MatrixFreeTest
{
public:
  MatrixFreeTest(const std::shared_ptr<const MatrixFree<dim, Number>> &data)
    : data(data)
  {}

  void
  vmult(BlockCSRMatrix<Number> &dst, const BlockCSRMatrix<Number> &src) const
  {
    dst = 0;
    data->cell_loop(&MatrixFreeTest::local_apply_cell,
                    this,
                    dst,
                    src,
                    /*zero_dst_vector*/ false);
  }

private:
  void
  local_apply_cell(
    const MatrixFree<dim, Number> & /*data*/,
    BlockCSRMatrix<Number> &                     dst,
    const BlockCSRMatrix<Number> &               src,
    const std::pair<unsigned int, unsigned int> &cell_range) const
  {
    FEEvaluation<dim, fe_degree, n_q_points_1d, n_components, Number> phi(
      *data);

    const auto &              dof_info = data->get_dof_info();
    std::vector<unsigned int> my_rows;
    my_rows.reserve(phi.dofs_per_component *
                    VectorizedArray<Number>::n_array_elements);

    const VectorizedArray<Number> zero = make_vectorized_array<Number>(0);

    std::vector<bool> touched;

    unsigned int nonzero_columns = 0;
    for (unsigned int cell = cell_range.first; cell < cell_range.second; ++cell)
      {
        // collect DoFs on all cell blocks
        dof_info.get_dof_indices_on_cell_batch(my_rows,
                                               cell,
                                               false /*apply_constraints*/);

        deallog << "Rows on cell: " << cell << std::endl;
        for (auto r : my_rows)
          deallog << " " << r;
        deallog << std::endl;

        DoFInfo dof_info;
        dof_info.initialize(my_rows, src.get_row_blocks());

        BlockCSRMatrixIterators::RowsBlockAccessor<Number, true>
          src_row_accessor(&src, dof_info);
        BlockCSRMatrixIterators::RowsBlockAccessor<Number, false>
          dst_row_accessor(&dst, dof_info);


        AssertDimension(my_rows.size(),
                        data->n_components_filled(cell) *
                          phi.dofs_per_component);

        types::global_dof_index src_col = src_row_accessor.reinit(0);
        types::global_dof_index dst_col = dst_row_accessor.reinit(0);

        phi.reinit(cell);

        // sanity/understanding check. Output JxW to see what happens
        // with un-filled cells
        deallog << "JxW at first quadrature point:";
        for (unsigned int i = 0; i < VectorizedArray<Number>::n_array_elements;
             ++i)
          deallog << " " << phi.JxW(0)[i];
        deallog << std::endl;

        while (src_col != numbers::invalid_dof_index)
          {
            const auto N = src_row_accessor.get_col_block_size();
            deallog << "Column block: " << src_col << " of size " << N
                    << std::endl;
            while (dst_col != numbers::invalid_dof_index && dst_col != src_col)
              dst_col = dst_row_accessor.advance();

            // at this point we should really have matching columns:
            AssertDimension(src_col, dst_col);
            AssertDimension(N, dst_row_accessor.get_col_block_size());

            // set touched based on sparsity of source for this column block
            touched.assign(my_rows.size(), false);

            // now do the standard matrix-free operations using accessors
            for (unsigned int c = 0; c < N; ++c)
              {
                // for scalar-valued operators submit_dof_value() takes
                // VectorizedArray<Number> val_in for each
                // phi.dofs_per_component
                VectorizedArray<Number> val_in = zero;

                src_row_accessor.process_active_rows_vectorized(
                  [&](const ArrayView<
                        const std::pair<unsigned int, unsigned int>> &dof_view,
                      typename BlockCSRMatrixIterators::RowsBlockAccessor<
                        double,
                        true>::vectorized_pointer const val,
                      const unsigned int                stride) {
                    for (unsigned int i = 0; i < dof_view.size(); ++i)
                      {
                        // we know that we have one cell only, so we will set
                        // first component only.
                        // FIXME: adjust loop above and use val directly
                        val_in[0] = *(&val[dof_view[i].first * stride][0] + c);
                        phi.submit_dof_value(val_in, dof_view[i].second);
                        touched[dof_view[i].second] = true;
                      }
                  });

                // reset all non-touched dofs to zero:
                for (unsigned int dof = 0; dof < phi.dofs_per_component; ++dof)
                  if (!touched[dof])
                    phi.submit_dof_value(zero, dof);

                if (src_col == 2 && c == 1)
                  {
                    deallog << "Block " << src_col << " subcol " << c
                            << " dof values:";
                    for (unsigned int dof = 0; dof < phi.dofs_per_component;
                         ++dof)
                      deallog << " " << phi.get_dof_value(dof)[0];
                    deallog << std::endl;
                  }

                // standard matrix-free kernel
                phi.evaluate(true, true, false);
                for (unsigned int q = 0; q < phi.n_q_points; ++q)
                  {
                    phi.submit_gradient(0.5 * phi.get_gradient(q), q);
                    phi.submit_value(phi.get_value(q), q);
                  }
                phi.integrate(true, true);

                if (src_col == 2 && c == 1)
                  {
                    deallog << "Block " << src_col << " subcol " << c
                            << " dof values:";
                    for (unsigned int dof = 0; dof < phi.dofs_per_component;
                         ++dof)
                      deallog << " " << phi.get_dof_value(dof)[0];
                    deallog << std::endl;
                  }

                // distribute manually
                dst_row_accessor.process_active_rows_vectorized(
                  [&](const ArrayView<
                        const std::pair<unsigned int, unsigned int>> &dof_view,
                      typename BlockCSRMatrixIterators::RowsBlockAccessor<
                        double,
                        false>::vectorized_pointer const val,
                      const unsigned int                 stride) {
                    for (unsigned int i = 0; i < dof_view.size(); ++i)
                      {
                        // we know that we have one cell only, so we will set
                        // first component only
                        // FIXME: adjust loop over c and use val directly
                        *(&val[dof_view[i].first * stride][0] + c) +=
                          phi.get_dof_value(dof_view[i].second)[0];
                      }
                  });
              }

            // finally advance the src accessor:
            src_col = src_row_accessor.advance();
            ++nonzero_columns;
          }

      } // loop over all cells

    deallog << "Nonzero column blocks over "
            << cell_range.second - cell_range.first
            << " cells: " << nonzero_columns << std::endl;
  }

  const std::shared_ptr<const MatrixFree<dim, Number>> data;
};



template <int dim,
          typename Number   = double,
          int fe_degree     = 2,
          int n_q_points_1d = fe_degree + 1>
void
test(const unsigned int n_cells = 1)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_core =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
  parallel::distributed::Triangulation<dim> triangulation(
    mpi_communicator,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  // Setup system
  {
    std::vector<unsigned int> repetitions(dim, 1);
    repetitions[0] = n_cells;
    const Point<dim> p1;
    Point<dim>       p2;
    for (unsigned int d = 1; d < dim; ++d)
      p2[d] = 1.;
    p2[0] = n_cells;
    GridGenerator::subdivided_hyper_rectangle(
      triangulation, repetitions, p1, p2, true);
  }

  DoFHandler<dim> dh(triangulation);

  FE_Q<dim> fe(fe_degree);
  dh.distribute_dofs(fe);

  const std::vector<unsigned int> row_blocks = renumber_based_on_nodes(dh);

  // now test with evaluating
  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dh, locally_relevant_dofs);

  AffineConstraints<double> constraints;
  constraints.reinit(locally_relevant_dofs);
  DoFTools::make_hanging_node_constraints(dh, constraints);

  // FIXME: figure out how to work with Dirichlet constraints
  // or how to distribute them on src vector efficiently
  // VectorTools::interpolate_boundary_values(
  //   dh, 0 /*left side*/, ZeroFunction<dim>(), constraints);

  constraints.close();

  deallog << "Constraints:" << std::endl;
  constraints.print(deallog.get_file_stream());

  const auto n_row_blocks = row_blocks.size();

  const std::vector<unsigned int> col_blocks = {{2, 1, 3}};

  DynamicSparsityPattern sp_src(n_row_blocks, col_blocks.size());
  DynamicSparsityPattern sp_dst(n_row_blocks, col_blocks.size());

  // after renumbering we know that we have fe_degree+1 nodes with the same x
  // coordinate in one row block.
  // FIXME: use setup_1d_sparsity() to get more interesting blocking.

  // make 3 blocks so that on one cell we have
  // 1: both dst/src non-empty
  {
    const unsigned int col = 0;
    for (unsigned int i = 0; i < n_row_blocks; ++i)
      {
        sp_src.add(i, col);
        sp_dst.add(i, col);
      }
  }

  // 2: both dst/src empty (let this be column 1)

  // 3: some blocks in src are empty, whereas dst is of course full
  {
    const unsigned int col = 2;
    sp_src.add(0, col);
    sp_src.add(1, col);
    for (unsigned int i = 0; i < n_row_blocks; ++i)
      sp_dst.add(i, col);
  }

  std::shared_ptr<BlockIndices> rows =
    std::make_shared<BlockIndices>(row_blocks);
  std::shared_ptr<BlockIndices> cols =
    std::make_shared<BlockIndices>(col_blocks);

  const types::global_dof_index full_rows = dh.n_dofs();
  const types::global_dof_index full_cols =
    std::accumulate(col_blocks.begin(),
                    col_blocks.end(),
                    types::global_dof_index(0));

  deallog << "Sparsity src:" << std::endl;
  sp_src.print(deallog.get_file_stream());
  deallog << "Sparsity dst:" << std::endl;
  sp_dst.print(deallog.get_file_stream());

  // prepare MF data
  std::shared_ptr<MatrixFree<dim, Number>> fine_level_data =
    std::make_shared<MatrixFree<dim, Number>>();
  {
    typename MatrixFree<dim, Number>::AdditionalData fine_level_additional_data;
    fine_level_additional_data.mapping_update_flags =
      update_values | update_gradients | update_JxW_values |
      update_quadrature_points;
    fine_level_additional_data.overlap_communication_computation = false;
    fine_level_additional_data.tasks_parallel_scheme =
      MatrixFree<dim, Number>::AdditionalData::none;
    fine_level_data->reinit(dh,
                            constraints,
                            QGauss<1>(n_q_points_1d),
                            fine_level_additional_data);
  }

  const std::shared_ptr<const dealii::Utilities::MPI::Partitioner>
    bcsr_row_part = fine_level_data->get_vector_partitioner();

  // setup matrices
  BlockCSRMatrix<Number> src, dst;
  src.reinit(sp_src, rows, cols, bcsr_row_part);
  dst.reinit(sp_dst, rows, cols, bcsr_row_part);

  deallog << "Internal sparsity src:" << std::endl;
  src.get_sparsity_pattern().print(deallog.get_file_stream());

  deallog << "Internal sparsity dst:" << std::endl;
  dst.get_sparsity_pattern().print(deallog.get_file_stream());

  // randomize content of matrices
  const auto randomize_mat = [](BlockCSRMatrix<Number> &mat) {
    for (unsigned int i = 0; i < mat.n_local_row_blocks(); ++i)
      {
        const auto M = mat.get_row_blocks()->block_size(i);
        for (auto it = mat.begin_local(i); it != mat.end_local(i); ++it)
          {
            const auto j = it->column();
            const auto N = mat.get_col_blocks()->block_size(j);
            for (unsigned int jj = 0; jj < N; ++jj)
              for (unsigned int ii = 0; ii < M; ++ii)
                *(it->data() +
                  BlockCSRMatrix<double>::local_index(ii, jj, M, N)) =
                  Utilities::generate_normal_random_number(0, 0.2);
          }
      }
  };

  randomize_mat(src);

  MatrixFreeTest<dim, fe_degree, n_q_points_1d, Number> mf_test(
    fine_level_data);
  mf_test.vmult(dst, src);


  // now do the same using full serial vectors
  LAPACKFullMatrix<Number> full_src(full_rows, full_cols),
    full_dst(full_rows, full_cols), full_diff(full_rows, full_cols);

  FEEvaluation<dim, fe_degree, n_q_points_1d, 1, double> phi(*fine_level_data);

  src.copy_to(full_src);

  LinearAlgebra::distributed::Vector<Number> src_col(bcsr_row_part);
  LinearAlgebra::distributed::Vector<Number> dst_col(bcsr_row_part);

  const auto loc_w_ghost =
    bcsr_row_part->local_size() + bcsr_row_part->n_ghost_indices();

  full_dst = 0.;

  // understanding check:
  // set src_col
  {
    deallog << "FEEvaluation dof relationship check:" << std::endl;
    for (unsigned int i = 0; i < loc_w_ghost; ++i)
      src_col.local_element(i) = i;

    const unsigned int        cell = 0;
    std::vector<unsigned int> my_rows;
    fine_level_data->get_dof_info().get_dof_indices_on_cell_batch(
      my_rows, cell, false /*apply_constraints*/);
    deallog << "get_dof_indices_on_cell_batch:" << std::endl;
    for (const auto &row : my_rows)
      deallog << " " << row;
    deallog << std::endl;
    phi.reinit(cell);
    phi.read_dof_values(src_col);
    deallog << "get_dof_value(dof):" << std::endl;
    for (unsigned int dof = 0; dof < phi.dofs_per_component; ++dof)
      deallog << dof << " " << phi.get_dof_value(dof)[0] << std::endl;
    deallog << std::endl;
  }

  for (unsigned int cell = 0; cell < fine_level_data->n_macro_cells(); ++cell)
    {
      phi.reinit(cell);
      for (unsigned int j = 0; j < full_cols; ++j)
        {
          for (unsigned int i = 0; i < loc_w_ghost; ++i)
            src_col.local_element(i) =
              full_src(bcsr_row_part->local_to_global(i), j);

          dst_col = 0.;

          phi.read_dof_values(src_col);

          if (j == 4)
            {
              deallog << "Column " << j << " dof values:";
              for (unsigned int dof = 0; dof < phi.dofs_per_component; ++dof)
                deallog << " " << phi.get_dof_value(dof)[0];
              deallog << std::endl;
            }

          phi.evaluate(true, true, false);
          for (unsigned int q = 0; q < phi.n_q_points; ++q)
            {
              phi.submit_gradient(0.5 * phi.get_gradient(q), q);
              phi.submit_value(phi.get_value(q), q);
            }
          phi.integrate(true, true);

          if (j == 4)
            {
              deallog << "Column " << j << " dof values:";
              for (unsigned int dof = 0; dof < phi.dofs_per_component; ++dof)
                deallog << " " << phi.get_dof_value(dof)[0];
              deallog << std::endl;
            }

          phi.distribute_local_to_global(dst_col);

          for (unsigned int i = 0; i < loc_w_ghost; ++i)
            full_dst(bcsr_row_part->local_to_global(i), j) +=
              dst_col.local_element(i);
        } // loop over all columns
    }     // loop over all cells

  Utilities::MPI::sum(full_dst, mpi_communicator, full_dst);

  dst.copy_to(full_diff);
  full_diff.add(-1, full_dst);

  deallog << "Diff:" << std::endl;
  full_diff.print_formatted(deallog.get_file_stream(), 3, true, 0, "0");

  deallog << "frobenius_norm src: " << full_src.frobenius_norm() << std::endl
          << "frobenius_norm dst: " << full_dst.frobenius_norm() << std::endl
          << "linfty_norm diff:   " << full_diff.linfty_norm() << std::endl;

  // print grid and DoFs for visual inspection
  if (dim == 2)
    {
      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping, dh, support_points);

      const std::string base_filename =
        "grid" + dealii::Utilities::int_to_string(dim) + "_p" +
        dealii::Utilities::int_to_string(this_mpi_core);

      const std::string filename = base_filename + ".gp";
      std::ofstream     f(filename.c_str());

      f << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
        << std::endl
        << "set output \"" << base_filename << ".png\"" << std::endl
        << "set size square" << std::endl
        << "set view equal xy" << std::endl
        << "unset xtics" << std::endl
        << "unset ytics" << std::endl
        << "unset grid" << std::endl
        << "unset border" << std::endl
        << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 "
           "offset 1,1 notitle"
        << std::endl;
      GridOut().write_gnuplot(triangulation, f);
      f << "e" << std::endl;

      /*
      // to make life easier, filter out all support points which are not in
      // the sets:
      for (auto it = support_points.cbegin(); it != support_points.cend(); )
        {
          if (local_support[i].is_element(it->first))
            ++it;
          else
            support_points.erase(it++);
        }
      */

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

      f << "e" << std::endl;
    }

  dh.clear();
  deallog << "Ok" << std::endl;
}

int
main(int argc, char **argv)
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  const unsigned int n_procs =
    dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  std::string   deallogname = "output" + dealii::Utilities::int_to_string(myid);
  std::ofstream logfile(deallogname);
  deallog.attach(logfile, /*do not print job id*/ false);
  deallog.depth_console(0);

  test<2>();

  logfile.close();

  MPI_Barrier(MPI_COMM_WORLD);

  if (myid == 0)
    for (unsigned int p = 0; p < n_procs; ++p)
      {
        std::string deallogname =
          "output" + dealii::Utilities::int_to_string(p);
        std::ifstream f(deallogname);
        std::string   line;
        while (std::getline(f, line))
          std::cout << p << ":" << line << std::endl;
      }

  return 0;
}
