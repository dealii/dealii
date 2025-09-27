// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Read in a large grid from a file and distribute hp-DoFs on it using
// FE_Q elements of different orders on different cells. The
// active_fe_index on each cell is determined in a mostly random way,
// but so that it is the same regardless of the number of processors.
//
// We used to treat hp-DoF unification on vertices and faces
// differently depending on whether we are in the interior of a
// subdomain or at a processor boundary. But later versions of the
// code did away with this distinction, and now the total number of
// DoFs, must be the same regardless of the number of subdomains.
//
// This test checks this on a large 2d mesh (~30k cells) and a large
// 3d mesh (~13k cells).


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include <numeric>

#include "../tests.h"


template <int dim>
void
test()
{
  parallel::distributed::Triangulation<dim> triangulation(
    MPI_COMM_WORLD, Triangulation<dim>::limit_level_difference_at_vertices);

  // First, read a complicated mesh
  GridIn<dim> gi;
  gi.attach_triangulation(triangulation);
  if (dim == 2)
    {
      std::ifstream in(SOURCE_DIR
                       "/../grid/grid_in_02/2d.xda"); // ~29k 2d cells
      gi.read_xda(in);
    }
  else
    {
      std::ifstream in(SOURCE_DIR "/../grid/grid_in_3d/4.in"); // ~14k 3d cells
      gi.read_xda(in);
    }

  // Then build a collection of FE_Q objects; to make things a bit
  // more interesting, duplicate each element twice in the collection
  // so that different active_fe_index values may correspond to
  // different FE objects but the same underlying FE
  hp::FECollection<dim> fe;
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int p = 1; p <= 5; ++p)
      fe.push_back(FE_Q<dim>(p));


  // Then more or less randomly assign elements to cells. We use a
  // coarse mesh, so the cell->active_cell_index() is globally unique
  // regardless of the number of processors involved, and we can use
  // that to build a hash value from it that is then used to assign an
  // active_fe_index
  DoFHandler<dim> dof_handler(triangulation);
  for (const auto &cell : dof_handler.active_cell_iterators() |
                            IteratorFilters::LocallyOwnedCell())
    cell->set_active_fe_index(
      (cell->active_cell_index() +
       13 * cell->active_cell_index() * cell->active_cell_index()) %
      fe.size());
  dof_handler.distribute_dofs(fe);

  deallog << "n_globally_active_cells: "
          << triangulation.n_global_active_cells() << std::endl;
  deallog << "n_locally_owned_dofs: " << dof_handler.n_locally_owned_dofs()
          << std::endl;
  deallog << "n_global_dofs: " << dof_handler.n_dofs() << std::endl;
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
