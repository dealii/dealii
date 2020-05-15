// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2018 by the deal.II authors
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



// Check register_data_attach and notify_ready_to_unpack
// for variable size transfer. In particular check that the
// transfer succeeds if one process does not get any data.

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
std::vector<char>
pack_function(
  const typename parallel::distributed::Triangulation<dim, dim>::cell_iterator
    &cell,
  const typename parallel::distributed::Triangulation<dim, dim>::CellStatus
    status)
{
  // Only pack data for one cell. The test is run for 3 ranks. The cell is
  // chosen to be one that is transferred from rank 0 to rank 1 so that we
  // test both that the transfer from 0 to 1 is successful, and that rank
  // 3 does not complain about not receiving any data.
  unsigned int      some_number = (cell->id().to_string() == "0_1:2") ? 5 : 0;
  char              some_char   = 'a';
  std::vector<char> buffer(some_number, some_char);

  deallog << "packing cell " << cell->id()
          << " with data size =" << buffer.size() << " status=";
  if (status == parallel::distributed::Triangulation<dim, dim>::CELL_PERSIST)
    deallog << "PERSIST";
  else if (status ==
           parallel::distributed::Triangulation<dim, dim>::CELL_REFINE)
    deallog << "REFINE";
  else if (status ==
           parallel::distributed::Triangulation<dim, dim>::CELL_COARSEN)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status == parallel::distributed::Triangulation<dim, dim>::CELL_COARSEN)
    {
      Assert(cell->has_children(), ExcInternalError());
    }
  else
    {
      Assert(!cell->has_children(), ExcInternalError());
    }

  return buffer;
}

template <int dim>
void
unpack_function(
  const typename parallel::distributed::Triangulation<dim, dim>::cell_iterator
    &cell,
  const typename parallel::distributed::Triangulation<dim, dim>::CellStatus
                                                                  status,
  const boost::iterator_range<std::vector<char>::const_iterator> &data_range)
{
  const unsigned int data_in_bytes =
    std::distance(data_range.begin(), data_range.end());

  deallog << "unpacking cell " << cell->id() << " with data size="
          << std::distance(data_range.begin(), data_range.end())
          << " and data: '" << std::string(data_range.begin(), data_range.end())
          << "' status=";
  if (status == parallel::distributed::Triangulation<dim, dim>::CELL_PERSIST)
    deallog << "PERSIST";
  else if (status ==
           parallel::distributed::Triangulation<dim, dim>::CELL_REFINE)
    deallog << "REFINE";
  else if (status ==
           parallel::distributed::Triangulation<dim, dim>::CELL_COARSEN)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status == parallel::distributed::Triangulation<dim, dim>::CELL_REFINE)
    {
      Assert(cell->has_children(), ExcInternalError());
    }
  else
    {
      Assert(!cell->has_children(), ExcInternalError());
    }
}


template <int dim>
void
test()
{
  unsigned int myid     = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  unsigned int numprocs = Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::subdivided_hyper_cube(tr, 2);
      tr.refine_global(1);
      deallog << "cells before: " << tr.n_global_active_cells() << std::endl;

      typename Triangulation<dim, dim>::active_cell_iterator cell;

      for (cell = tr.begin_active(); cell != tr.end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              if (cell->id().to_string() == "0_1:0")
                cell->set_refine_flag();
              else if (cell->parent()->id().to_string() == "3_0:")
                cell->set_coarsen_flag();

              deallog << "myid=" << myid << " cellid=" << cell->id();
              if (cell->coarsen_flag_set())
                deallog << " coarsening" << std::endl;
              else if (cell->refine_flag_set())
                deallog << " refining" << std::endl;
              else
                deallog << std::endl;
            }
        }

      unsigned int handle =
        tr.register_data_attach(pack_function<dim>,
                                /*returns_variable_size_data=*/true);

      deallog << "handle=" << handle << std::endl;
      tr.execute_coarsening_and_refinement();

      deallog << "cells after: " << tr.n_global_active_cells() << std::endl;

      /*
      for (cell = tr.begin_active();
           cell != tr.end();
           ++cell)
        {
      if (cell->is_locally_owned())
      deallog << "myid=" << myid << " cellid=" << cell->id() << std::endl;
      }*/

      tr.notify_ready_to_unpack(handle, unpack_function<dim>);

      const unsigned int checksum = tr.get_checksum();
      deallog << "Checksum: " << checksum << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
