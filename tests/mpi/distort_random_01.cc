// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check GridTools::distort_random


#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"



template <int dim>
void
test1(const bool keep_boundary)
{
  const unsigned int my_id = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  parallel::distributed::Triangulation<dim> tria(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  GridTools::distort_random(0.1, tria, keep_boundary);

  deallog << "dim=" << dim << ", keep_boundary=" << keep_boundary << std::endl;
  std::string filename;
  if (keep_boundary)
    filename = "keep_true-";
  else
    filename = "keep_false-";
  filename += Utilities::int_to_string(dim);

  std::ofstream logfile(
    (filename + "-" + Utilities::int_to_string(my_id, 2)).c_str());

  GridOut().write_gnuplot(tria, logfile);
  MPI_Barrier(MPI_COMM_WORLD);

  if (my_id == 0)
    for (unsigned int i = 0;
         i < Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD);
         ++i)
      {
        deallog << "Process " << i << ':' << std::endl;
        cat_file((filename + "-" + Utilities::int_to_string(i, 2)).c_str());
      }

  deallog << "OK" << std::endl;
}



int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test1<2>(true);
      test1<2>(false);
      test1<3>(true);
      test1<3>(false);
    }
  else
    {
      test1<2>(true);
      test1<2>(false);
      test1<3>(true);
      test1<3>(false);
    }

  return 0;
}
