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



// create a parallel DoFHandler and output data on a single
// cell. DataOut was not prepared to handle situations where a
// processor has no active cells at all.
// We are checking the output by proc 0, who has 0 own cells to write.
// Note that we can not write an empty file, because the collection
// needs to be readable by paraview and visit. We check that
// the format stays as it is right now.

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test()
{
  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    deallog << "hyper_cube" << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dofh(tr);

  static const FE_Q<dim> fe(2);
  dofh.distribute_dofs(fe);

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dofh);

  TrilinosWrappers::MPI::Vector x;
  x.reinit(dofh.locally_owned_dofs(), MPI_COMM_WORLD);
  x = 2.0;

  data_out.add_data_vector(x, "x");
  data_out.build_patches();

  std::vector<types::global_dof_index> n_locally_owned_dofs_per_processor =
    Utilities::MPI::all_gather(MPI_COMM_WORLD, dofh.n_locally_owned_dofs());
  if (myid == 0)
    {
      for (unsigned int i = 0; i < n_locally_owned_dofs_per_processor.size();
           ++i)
        deallog << n_locally_owned_dofs_per_processor[i] << std::endl;
      data_out.write_vtu(deallog.get_file_stream());
    }
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
