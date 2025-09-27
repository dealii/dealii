// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test DataOut::write_deal_II_intermediate_in_parallel() and
// DataOutReader::read_whole_parallel_file()

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
    "test.pd2", MPI_COMM_WORLD, DataOutBase::CompressionLevel::no_compression);

  const unsigned int my_rank =
    dealii::Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  if (my_rank == 0)
    {
      // Read the data back in and dump it into the deallog:
      std::ifstream in("test.pd2");
      Assert(in, dealii::ExcIO());
      DataOutReader<dim, dim> reader;
      reader.read_whole_parallel_file(in);
      reader.write_deal_II_intermediate(deallog.get_file_stream());
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
