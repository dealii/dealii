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
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



// test manual repartitioning and transferring data

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/base/utilities.h>



template <int dim>
void
pack_function (const typename parallel::distributed::Triangulation<dim,dim>::cell_iterator &cell,
               const typename parallel::distributed::Triangulation<dim,dim>::CellStatus status,
               void *data)
{
  static int some_number = cell->index();
  deallog << "packing cell " << cell->id() << " with data=" << some_number << " status=";
  if (status==parallel::distributed::Triangulation<dim,dim>::CELL_PERSIST)
    deallog << "PERSIST";
  else if (status==parallel::distributed::Triangulation<dim,dim>::CELL_REFINE)
    deallog << "REFINE";
  else if (status==parallel::distributed::Triangulation<dim,dim>::CELL_COARSEN)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status==parallel::distributed::Triangulation<dim,dim>::CELL_COARSEN)
    {
      Assert(cell->has_children(), ExcInternalError());
    }
  else
    {
      Assert(!cell->has_children(), ExcInternalError());
    }

  int *intdata = reinterpret_cast<int *>(data);
  *intdata = some_number;

  ++some_number;
}

template <int dim>
void
unpack_function (const typename parallel::distributed::Triangulation<dim,dim>::cell_iterator &cell,
                 const typename parallel::distributed::Triangulation<dim,dim>::CellStatus status,
                 const void *data)
{
  const int *intdata = reinterpret_cast<const int *>(data);

  deallog << "unpacking cell " << cell->id() << " with data=" << (*intdata) << " status=";
  if (status==parallel::distributed::Triangulation<dim,dim>::CELL_PERSIST)
    deallog << "PERSIST";
  else if (status==parallel::distributed::Triangulation<dim,dim>::CELL_REFINE)
    deallog << "REFINE";
  else if (status==parallel::distributed::Triangulation<dim,dim>::CELL_COARSEN)
    deallog << "COARSEN";
  deallog << std::endl;

  if (status==parallel::distributed::Triangulation<dim,dim>::CELL_REFINE)
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
  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);

  if (true)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD,
                                                   dealii::Triangulation<dim,dim>::none,
                                                   parallel::distributed::Triangulation<dim>::no_automatic_repartitioning);

      GridGenerator::hyper_cube(tr);
      tr.refine_global(1);

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;

      deallog << "* global refine:" << std::endl;

      unsigned int offset = tr.register_data_attach(sizeof(int), pack_function<dim>);

      tr.refine_global(1);

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;

      tr.notify_ready_to_unpack(offset, unpack_function<dim>);

      //tr.write_mesh_vtk("a");




      deallog << "* repartition:" << std::endl;

      offset = tr.register_data_attach(sizeof(int), pack_function<dim>);

      tr.repartition();

      //tr.write_mesh_vtk("b");

      deallog << "locally owned cells: " << tr.n_locally_owned_active_cells()
              << " / "
              << tr.n_global_active_cells()
              << std::endl;

      tr.notify_ready_to_unpack(offset, unpack_function<dim>);

      const unsigned int checksum = tr.get_checksum ();
      if (myid == 0)
        deallog << "Checksum: "
                << checksum
                << std::endl;
    }

  if (Utilities::MPI::this_mpi_process (MPI_COMM_WORLD) == 0)
    deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);
  MPILogInitAll log;
  test<2>();
}
