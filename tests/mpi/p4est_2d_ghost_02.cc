// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2020 by the deal.II authors
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



// check existence of ghost layer in 2d with global refinement

#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (true)
    {
      if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
        deallog << "hyper_cube" << std::endl;

      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      for (int i = 0; i < 5; ++i)
        {
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            deallog << "refine loop:" << i << std::endl;

          tr.refine_global(1);

          if (myid == 0)
            {
              std::vector<types::subdomain_id> cell_subd(tr.n_active_cells());

              GridTools::get_subdomain_association(tr, cell_subd);
              for (unsigned int i = 0; i < tr.n_active_cells(); ++i)
                deallog << cell_subd[i] << " ";
              deallog << std::endl;
            }

          // check that all local
          // neighbors have the
          // correct level
          typename Triangulation<dim, dim>::active_cell_iterator cell;

          for (cell = tr.begin_active(); cell != tr.end(); ++cell)
            {
              if (cell->subdomain_id() != (unsigned int)myid)
                {
                  AssertThrow(cell->is_ghost() || cell->is_artificial(),
                              ExcInternalError());
                  continue;
                }

              for (const unsigned int n : GeometryInfo<dim>::face_indices())
                {
                  if (cell->at_boundary(n))
                    continue;
                  AssertThrow(cell->neighbor(n).state() == IteratorState::valid,
                              ExcInternalError());

                  AssertThrow(cell->neighbor(n)->level() == cell->level(),
                              ExcInternalError());

                  AssertThrow(!cell->neighbor(n)->has_children(),
                              ExcInternalError());
                }
            }

          const unsigned int checksum = tr.get_checksum();
          if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
            deallog << "Checksum: " << checksum << std::endl;

          AssertThrow(
            tr.n_global_active_cells() ==
              static_cast<unsigned int>(
                std::pow(1. * GeometryInfo<dim>::max_children_per_cell, i + 1)),
            ExcInternalError());
        }
    }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();

      deallog.push("2d");
      test<2>();
      deallog.pop();
    }
  else
    test<2>();
}
