// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2017 by the deal.II authors
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


// create a shared tria mesh and distribute with a custom function

#include "../tests.h"
#include <deal.II/base/tensor.h>
#include <deal.II/grid/tria.h>
#include <deal.II/distributed/tria.h>
#include <deal.II/distributed/shared_tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>


template <int dim>
void
mypartition(parallel::shared::Triangulation<dim> &tria)
{
  const unsigned int nproc = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  typename Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active();
  for (unsigned int i=0; i < tria.n_active_cells(); ++i)
    {
      unsigned int j=0;
      if (i<4)
        j=1;
      else if (i<8)
        j=2;

      cell->set_subdomain_id(j % nproc);
      ++cell;
    }
}

template <int dim>
void
test()
{
  parallel::shared::Triangulation<dim> shared_tria(MPI_COMM_WORLD,typename Triangulation<dim>::MeshSmoothing
                                                   (Triangulation<dim>::limit_level_difference_at_vertices), true,
                                                   parallel::shared::Triangulation<dim>::partition_custom_signal);
  shared_tria.signals.post_refinement.connect (std::bind(&mypartition<dim>, std::ref(shared_tria)));


  GridGenerator::hyper_L(shared_tria);
  shared_tria.refine_global();
  typename Triangulation<dim>::active_cell_iterator it = shared_tria.begin_active();
  it->set_refine_flag();
  shared_tria.execute_coarsening_and_refinement ();

  deallog << "(CellId,subdomain_id) for each active cell:" << std::endl;
  typename Triangulation<dim>::active_cell_iterator
  cell = shared_tria.begin_active(),
  endc = shared_tria.end();
  for (; cell!=endc; ++cell)
    if (cell->subdomain_id() != numbers::artificial_subdomain_id)
      deallog << "(" << cell->id().to_string() << "," << cell->subdomain_id() << ")" << std::endl;
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

  deallog << "OK" << std::endl;
}

