// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

// Check that we can create a sparsity pattern from two DoFHandlers with the
// same distributed mesh. This previously failed since the implementation
// called a GridTools function that assumed all cells without children were
// active (which is not the case for distributed meshes).

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_cartesian.h>

#include <deal.II/grid/grid_generator.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include "../tests.h"


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;

  parallel::distributed::Triangulation<2> triangulation(MPI_COMM_WORLD);
  GridGenerator::subdivided_hyper_rectangle(triangulation,
                                            {4u, 1u},
                                            {0.0, 0.0},
                                            {4.0, 1.0});
  triangulation.refine_global(1);
  deallog << "global number of active cells: "
          << triangulation.n_global_active_cells() << std::endl;
  deallog << "local number of active cells: " << triangulation.n_active_cells()
          << std::endl;
  deallog << "number of locally owned active cells: "
          << triangulation.n_locally_owned_active_cells() << std::endl;

  FE_DGQ<2>     fe_0(0);
  DoFHandler<2> dof_handler_0;
  FE_Q<2>       fe_1(1);
  DoFHandler<2> dof_handler_1;
  dof_handler_0.initialize(triangulation, fe_0);
  dof_handler_1.initialize(triangulation, fe_1);

  IndexSet locally_relevant_dofs_0;
  IndexSet locally_relevant_dofs_1;
  DoFTools::extract_locally_relevant_dofs(dof_handler_0,
                                          locally_relevant_dofs_0);
  DoFTools::extract_locally_relevant_dofs(dof_handler_1,
                                          locally_relevant_dofs_1);
  deallog << "locally owned dofs 0: ";
  dof_handler_0.locally_owned_dofs().print(deallog);
  deallog << std::endl;
  deallog << "locally owned dofs 1: ";
  dof_handler_1.locally_owned_dofs().print(deallog);
  deallog << std::endl;
  deallog << "locally relevant dofs 0: ";
  locally_relevant_dofs_0.print(deallog);
  deallog << std::endl;
  deallog << "locally relevant dofs 1: ";
  locally_relevant_dofs_1.print(deallog);
  deallog << std::endl;

#if 0 // reenable this to see where the DoFs are
  MappingCartesian<2> mapping;
  {
    deallog << "support points of dofs 0:" << std::endl;
    std::vector<types::global_dof_index> cell_dofs(fe_0.dofs_per_cell);
    std::size_t cell_n = 0;
    for (const auto &cell : dof_handler_0.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            deallog << "cell_n = " << cell_n << std::endl;
            const std::vector<Point<2>> &support_points = fe_0.get_unit_support_points();
            cell->get_dof_indices(cell_dofs);
            for (std::size_t dof_n = 0; dof_n < fe_0.dofs_per_cell; ++dof_n)
              {
                deallog << "dof " << std::setw(2) << cell_dofs[dof_n] << ": "
                        << mapping.transform_unit_to_real_cell(cell, support_points[dof_n])
                        << std::endl;
              }
            ++cell_n;
          }
      }
  }

  {
    deallog << "support points of dofs 1:" << std::endl;
    std::vector<types::global_dof_index> cell_dofs(fe_1.dofs_per_cell);
    std::size_t cell_n = 0;
    for (const auto &cell : dof_handler_1.active_cell_iterators())
      {
        if (cell->is_locally_owned())
          {
            deallog << "cell_n = " << cell_n << std::endl;
            const std::vector<Point<2>> &support_points = fe_1.get_unit_support_points();
            cell->get_dof_indices(cell_dofs);
            for (std::size_t dof_n = 0; dof_n < fe_1.dofs_per_cell; ++dof_n)
              {
                deallog << "dof " << std::setw(2) << cell_dofs[dof_n] << ": "
                        << mapping.transform_unit_to_real_cell(cell, support_points[dof_n])
                        << std::endl;
              }
            ++cell_n;
          }
      }
  }
#endif

  DynamicSparsityPattern dynamic_sparsity_pattern(
    dof_handler_0.locally_owned_dofs().size(),
    dof_handler_1.locally_owned_dofs().size(),
    dof_handler_0.locally_owned_dofs());
  DoFTools::make_sparsity_pattern(dof_handler_0,
                                  dof_handler_1,
                                  dynamic_sparsity_pattern);
  dynamic_sparsity_pattern.print(deallog.get_file_stream());
  SparsityTools::distribute_sparsity_pattern(
    dynamic_sparsity_pattern,
    dof_handler_0.locally_owned_dofs(),
    MPI_COMM_WORLD,
    dof_handler_0.locally_owned_dofs());
  dynamic_sparsity_pattern.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}
