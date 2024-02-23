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



// save and load a triangulation with a different number of cpus
// with variable size data attach
// this is a combination of tests p4est_save_04 and attach_data_02

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
  static unsigned int       some_number = 1;
  std::vector<unsigned int> some_vector(some_number);
  for (unsigned int i = 0; i < some_number; ++i)
    some_vector[i] = i;

  std::vector<char> buffer;
  buffer.reserve(some_number * sizeof(unsigned int));
  for (auto vector_it = some_vector.cbegin(); vector_it != some_vector.cend();
       ++vector_it)
    {
      Utilities::pack(*vector_it, buffer, /*allow_compression=*/false);
    }

  deallog << "packing cell " << cell->id()
          << " with data size=" << buffer.size() << " accumulated data="
          << std::accumulate(some_vector.begin(), some_vector.end(), 0)
          << std::endl;

  Assert((status == CellStatus::cell_will_persist), ExcInternalError());

  ++some_number;
  return buffer;
}



template <int dim>
void
unpack_function(
  const typename parallel::distributed::Triangulation<dim, dim>::cell_iterator
                                                                 &cell,
  const CellStatus                                                status,
  const boost::iterator_range<std::vector<char>::const_iterator> &data_range)
{
  const unsigned int data_in_bytes =
    std::distance(data_range.begin(), data_range.end());

  std::vector<unsigned int> intdatavector(data_in_bytes / sizeof(unsigned int));

  auto vector_it = intdatavector.begin();
  auto data_it   = data_range.begin();
  for (; data_it != data_range.end();
       ++vector_it, data_it += sizeof(unsigned int))
    {
      *vector_it =
        Utilities::unpack<unsigned int>(data_it,
                                        data_it + sizeof(unsigned int),
                                        /*allow_compression=*/false);
    }

  deallog << "unpacking cell " << cell->id() << " with data size="
          << std::distance(data_range.begin(), data_range.end())
          << " accumulated data="
          << std::accumulate(intdatavector.begin(), intdatavector.end(), 0)
          << std::endl;

  Assert((status == CellStatus::cell_will_persist), ExcInternalError());
}



template <int dim>
void
test()
{
  unsigned int myid    = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  MPI_Comm     com_all = MPI_COMM_WORLD;
  MPI_Comm     com_small;

  // split the communicator in proc 0,1,2 and 3,4
  MPI_Comm_split(com_all, (myid < 3) ? 0 : 1, myid, &com_small);

  // write with small com
  if (myid < 3)
    {
      deallog << "writing with " << Utilities::MPI::n_mpi_processes(com_small)
              << std::endl;

      parallel::distributed::Triangulation<dim> tr(com_small);
      GridGenerator::subdivided_hyper_cube(tr, 2);
      tr.refine_global(1);

      typename Triangulation<dim, dim>::active_cell_iterator cell;
      for (cell = tr.begin_active(); cell != tr.end(); ++cell)
        {
          if (cell->is_locally_owned())
            {
              if (cell->id().to_string() == "0_1:0")
                cell->set_refine_flag();
              else if (cell->parent()->id().to_string() == "3_0:")
                cell->set_coarsen_flag();
            }
        }
      tr.execute_coarsening_and_refinement();

      unsigned int handle =
        tr.register_data_attach(pack_function<dim>,
                                /*returns_variable_size_data=*/true);

      tr.save("file");
      deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
      deallog << "Checksum: " << tr.get_checksum() << std::endl;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  deallog << "reading with " << Utilities::MPI::n_mpi_processes(com_all)
          << std::endl;

  {
    parallel::distributed::Triangulation<dim> tr(com_all);

    GridGenerator::subdivided_hyper_cube(tr, 2);
    tr.load("file");

    unsigned int handle =
      tr.register_data_attach(pack_function<dim>,
                              /*returns_variable_size_data=*/true);

    tr.notify_ready_to_unpack(handle, unpack_function<dim>);

    deallog << "#cells = " << tr.n_global_active_cells() << std::endl;
    deallog << "Checksum: " << tr.get_checksum() << std::endl;
  }

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  test<2>();
}
