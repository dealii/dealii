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



// Tests DoFInfo::get_dof_indices()

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/matrix_free/matrix_free.h>

#include <deal.II/numerics/vector_tools.h>

#include <iostream>

#include "../tests.h"



template <int dim, int fe_degree, typename number = double>
void
test(const bool adaptive_ref = true)
{
  MPI_Comm           mpi_communicator(MPI_COMM_WORLD);
  const unsigned int this_mpi_core =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);

  parallel::distributed::Triangulation<dim> tria(mpi_communicator);
  GridGenerator::hyper_cube(tria, 0, 1, true);
  tria.refine_global(1);
  if (adaptive_ref)
    {
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center().norm() < 0.5)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
      for (auto &cell : tria.active_cell_iterators())
        if (cell->is_locally_owned())
          if (cell->center()[0] < 0.2)
            cell->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }
  else
    {
      tria.refine_global(1);
    }

  FE_Q<dim>       fe(fe_degree);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe);

  IndexSet owned_set = dof.locally_owned_dofs();
  IndexSet relevant_set;
  DoFTools::extract_locally_relevant_dofs(dof, relevant_set);

  AffineConstraints<double> constraints(relevant_set);
  DoFTools::make_hanging_node_constraints(dof, constraints);
  // constrain bottom part of the boundary (lower in y direction)
  VectorTools::interpolate_boundary_values(dof,
                                           2,
                                           Functions::ZeroFunction<dim>(),
                                           constraints);
  constraints.close();

  deallog << "Testing " << dof.get_fe().get_name() << std::endl;

  std::shared_ptr<MatrixFree<dim, number>> mf_data(
    new MatrixFree<dim, number>());
  {
    const QGauss<1>                                  quad(fe_degree + 1);
    typename MatrixFree<dim, number>::AdditionalData data;
    data.tasks_parallel_scheme = MatrixFree<dim, number>::AdditionalData::none;
    data.tasks_block_size      = 7;
    mf_data->reinit(dof, constraints, quad, data);
  }

  const unsigned int     n_cells         = mf_data->n_macro_cells();
  const auto &           dof_info        = mf_data->get_dof_info();
  constexpr unsigned int n_vectorization = VectorizedArray<number>::size();

  std::vector<unsigned int> my_rows;
  my_rows.reserve(fe.dofs_per_cell * n_vectorization);

  constexpr auto dof_access_index =
    internal::MatrixFreeFunctions::DoFInfo::dof_access_cell;

  bool checked_interleaved = false;
  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      deallog << "Cell: " << cell << std::endl;
      checked_interleaved |=
        (dof_info.index_storage_variants[dof_access_index][cell] ==
         internal::MatrixFreeFunctions::DoFInfo::IndexStorageVariants::
           interleaved);
      auto get_and_log = [&](const bool apply_constraints) {
        dof_info.get_dof_indices_on_cell_batch(my_rows,
                                               cell,
                                               apply_constraints);
        // sort and make unique to make visual inspection easier:
        std::sort(my_rows.begin(), my_rows.end());
        my_rows.erase(std::unique(my_rows.begin(), my_rows.end()),
                      my_rows.end());
        for (auto el : my_rows)
          deallog << " " << el;
        deallog << std::endl;
      };
      get_and_log(true);
      get_and_log(false);
    }

  // if we don't have adaptive refinement, we should have
  // some cells interleaved (at least that is the intention of this
  // part of the test)
  Assert(checked_interleaved || adaptive_ref, ExcInternalError());

  // output in Gnuplot
  if (dim == 2)
    {
      std::map<types::global_dof_index, Point<dim>> support_points;
      MappingQ1<dim>                                mapping;
      DoFTools::map_dofs_to_support_points(mapping, dof, support_points);

      const std::string prefix =
        std::is_same<float, number>::value ? "float_" : "double_";
      const std::string href = (adaptive_ref ? "" : "global_");
      const std::string base_filename =
        prefix + href + "grid" + dealii::Utilities::int_to_string(dim) + "_" +
        dealii::Utilities::int_to_string(fe_degree) + "_p" +
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
        << "plot '-' using 1:2 with lines notitle, '-' with labels tc rgb 'red' nopoint notitle, '-' with labels point pt 4 offset 1,1 notitle"
        << std::endl;
      GridOut().write_gnuplot(tria, f);
      f << "e" << std::endl;

      // output cell blocks:
      for (unsigned int cell = 0; cell < n_cells; ++cell)
        for (unsigned int c = 0; c < mf_data->n_components_filled(cell); ++c)
          {
            const auto dof_cell = mf_data->get_cell_iterator(cell, c);
            f << dof_cell->center() << " \"" << cell << "\"\n";
          }

      f << std::flush;
      f << "e" << std::endl << std::endl;

      DoFTools::write_gnuplot_dof_support_point_info(f, support_points);

      f << "e" << std::endl;
    }
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  MPILogInitAll mpi_init_log;

  deallog.push("2d");
  test<2, 1>();
  test<2, 1, float>();
  test<2, 1>(false);
  deallog.pop();
}
