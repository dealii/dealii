// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test Tria::(level_)ghost_owners()

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/grid_refinement.h>
#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



// make sure if i is in s on proc j, j is in s on proc i
void
mpi_check(const std::set<types::subdomain_id> &s)
{
  MPI_Barrier(MPI_COMM_WORLD);
  unsigned int tag = 1234;
  for (std::set<types::subdomain_id>::iterator it = s.begin(); it != s.end();
       ++it)
    MPI_Send(nullptr, 0, MPI_INT, *it, tag, MPI_COMM_WORLD);

  for (unsigned int i = 0; i < s.size(); ++i)
    {
      MPI_Status status;
      MPI_Recv(
        nullptr, 0, MPI_INT, MPI_ANY_SOURCE, tag, MPI_COMM_WORLD, &status);
      if ((s.end() == s.find(status.MPI_SOURCE)))
        deallog << status.MPI_SOURCE << " NOTOKAY" << std::endl;
      Assert(s.end() != s.find(status.MPI_SOURCE), ExcInternalError());
    }
  MPI_Barrier(MPI_COMM_WORLD);
}


template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::limit_level_difference_at_vertices,
    parallel::distributed::Triangulation<dim>::construct_multigrid_hierarchy);

  std::vector<unsigned int> sub(dim, 1);
  sub[0] = 5;
  Point<dim> p;
  p[0] = 5;
  for (unsigned int i = 1; i < dim; ++i)
    p[i] = 1;
  GridGenerator::subdivided_hyper_rectangle(tr, sub, Point<dim>(), p);
  tr.refine_global(2);


  for (unsigned int ref = 0; ref <= 3; ++ref)
    {
      deallog << "* cycle " << ref << std::endl;

      deallog << "ghost owners: ";
      std::set<types::subdomain_id> ghost_owners = tr.ghost_owners();
      for (std::set<types::subdomain_id>::iterator it = ghost_owners.begin();
           it != ghost_owners.end();
           ++it)
        deallog << *it << ' ';
      deallog << std::endl;

      mpi_check(ghost_owners);

      deallog << "level ghost owners: ";
      std::set<types::subdomain_id> level_ghost_owners =
        tr.level_ghost_owners();
      for (std::set<types::subdomain_id>::iterator it =
             level_ghost_owners.begin();
           it != level_ghost_owners.end();
           ++it)
        deallog << *it << ' ';
      deallog << std::endl;

      mpi_check(level_ghost_owners);

      // owners need to be a subset of level owners:
      bool is_subset = std::includes(level_ghost_owners.begin(),
                                     level_ghost_owners.end(),
                                     ghost_owners.begin(),
                                     ghost_owners.end());
      Assert(is_subset, ExcInternalError());

      Vector<float>                 indicators(tr.n_active_cells());
      std::set<types::subdomain_id> neighbors;
      {
        for (typename Triangulation<dim>::active_cell_iterator cell =
               tr.begin_active();
             cell != tr.end();
             ++cell)
          if (!cell->is_artificial())
            {
              if (cell->is_ghost())
                neighbors.insert(cell->subdomain_id());
              indicators[cell->index()] = cell->index();
            }
      }

      Assert(neighbors == ghost_owners, ExcInternalError());

      parallel::distributed::GridRefinement::refine_and_coarsen_fixed_number(
        tr, indicators, 0.3, 0.0);
      tr.execute_coarsening_and_refinement();
      if (myid == 0)
        deallog << "total active cells = " << tr.n_global_active_cells()
                << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  MPILogInitAll all;
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
}
