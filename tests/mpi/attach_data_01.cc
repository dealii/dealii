// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check register_data_attach and notify_ready_to_unpack

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
  const CellStatus status)
{
  static int some_number = 0;
  deallog << "packing cell " << cell->id() << " with data=" << some_number
          << " status=";
  if (status == CellStatus::cell_will_persist)
    deallog << "PERSIST";
  else if (status == CellStatus::cell_will_be_refined)
    deallog << "REFINE";
  else if (status == CellStatus::children_will_be_coarsened)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status == CellStatus::children_will_be_coarsened)
    {
      Assert(cell->has_children(), ExcInternalError());
    }
  else
    {
      Assert(!cell->has_children(), ExcInternalError());
    }

  return Utilities::pack(some_number++, /*allow_compression=*/false);
}

template <int dim>
void
unpack_function(
  const typename parallel::distributed::Triangulation<dim, dim>::cell_iterator
                                                                 &cell,
  const CellStatus                                                status,
  const boost::iterator_range<std::vector<char>::const_iterator> &data_range)
{
  const int intdata = Utilities::unpack<int>(data_range.begin(),
                                             data_range.end(),
                                             /*allow_compression=*/false);

  deallog << "unpacking cell " << cell->id() << " with data=" << intdata
          << " status=";
  if (status == CellStatus::cell_will_persist)
    deallog << "PERSIST";
  else if (status == CellStatus::cell_will_be_refined)
    deallog << "REFINE";
  else if (status == CellStatus::children_will_be_coarsened)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status == CellStatus::cell_will_be_refined)
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
                                /*returns_variable_size_data=*/false);

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
