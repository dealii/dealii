// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test that DofRenumbering::Cuthill_McKee works in parallel also when
// a set of starting indices is given.
//
// this set of starting indices may also contain locally active DoFs,
// even though they will not be part of the renumbering on each
// cell. verify that this works by using as starting indices the ones
// that are located on processor boundaries


#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid   = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int nprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr, -1.0, 1.0);
  tr.refine_global(2);

  const FE_Q<dim> fe(1);
  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);

  deallog << "Before:" << std::endl;
  for (const auto &cell :
       dofh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      deallog << "locally owned cell: " << cell << std::endl;
      deallog << "       dof indices: ";

      std::vector<types::global_dof_index> cell_dofs(
        cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(cell_dofs);

      for (auto i : cell_dofs)
        deallog << i << ' ';
      deallog << std::endl;
    }

  std::set<types::global_dof_index> starting_indices;
  for (const auto &cell :
       dofh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    for (const unsigned int f : GeometryInfo<dim>::face_indices())
      if (!cell->at_boundary(f) && cell->neighbor(f)->is_ghost())
        {
          // we've identified a subdomain interface. use these DoFs
          // as starting indices
          std::vector<types::global_dof_index> face_dofs(
            cell->get_fe().dofs_per_face);
          cell->face(f)->get_dof_indices(face_dofs);
          for (auto i : face_dofs)
            starting_indices.insert(i);
        }

  DoFRenumbering::Cuthill_McKee(
    dofh,
    false,
    false,
    std::vector<types::global_dof_index>(starting_indices.begin(),
                                         starting_indices.end()));

  // output the renumbered DoF indices
  deallog << "After:" << std::endl;
  for (const auto &cell :
       dofh.active_cell_iterators() | IteratorFilters::LocallyOwnedCell())
    {
      deallog << "locally owned cell: " << cell << std::endl;
      deallog << "       dof indices: ";

      std::vector<types::global_dof_index> cell_dofs(
        cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(cell_dofs);

      for (auto i : cell_dofs)
        deallog << i << ' ';
      deallog << std::endl;
    }

  std::map<types::global_dof_index, Point<dim>> support_points;
  DoFTools::map_dofs_to_support_points(MappingQ1<dim>(), dofh, support_points);
  DoFTools::write_gnuplot_dof_support_point_info(deallog.get_file_stream(),
                                                 support_points);
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
