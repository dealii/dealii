// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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

// Test DataOut::write_deal_II_intermediate_in_parallel()

#include <deal.II/base/mpi.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/filtered_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
check()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);

  Vector<double> cell_data(tria.n_active_cells());
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
    cell_data(i) = i * 1.0;

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  Vector<double> x(dof_handler.n_dofs());

  DataOut<dim> data_out;
  data_out.attach_dof_handler(dof_handler);
  data_out.add_data_vector(cell_data,
                           "cell_data",
                           DataOut<dim>::type_cell_data);
  data_out.add_data_vector(x, "solution");

  data_out.build_patches();

  data_out.write_deal_II_intermediate_in_parallel(
    "test.d2p", MPI_COMM_WORLD, DataOutBase::VtkFlags::no_compression);

  const unsigned int my_rank =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (my_rank == 0)
    {
      // can't checksum as the file contains "[written by deal.II
      // 9.5.0-pre]" which would change often, so we just look at the
      // size for now. This will of course fail once we go to
      // 10.0. :-)
      std::ifstream in("test.d2p", std::ifstream::ate | std::ifstream::binary);
      Assert(in, dealii::ExcIO());
      deallog << "size: " << in.tellg() << std::endl;
    }

  deallog << "OK" << std::endl;
}

int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    all;

  check<2>();

  return 0;
}
