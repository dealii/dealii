// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// like the second part of 3d_refinement_06 but use only four coarse
// grid cells. we want to make sure that the 2:1 relationship holds
// across an edge


#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim>
void
test(std::ostream & /*out*/)
{
  {
    parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

    std::vector<unsigned int> subdivisions(3, 2);
    subdivisions[2] = 1;
    GridGenerator::subdivided_hyper_rectangle(tr,
                                              subdivisions,
                                              Point<3>(0, 0, 0),
                                              Point<3>(2, 2, 1));
    tr.begin_active()->set_refine_flag();
    tr.execute_coarsening_and_refinement();
    for (unsigned int c = 0; c < 8; ++c)
      tr.begin(0)->child(c)->set_refine_flag();
    tr.execute_coarsening_and_refinement();

    //    write_vtk(tr, "2");
    deallog << "cells test2: " << tr.n_active_cells() << std::endl;

    Assert(tr.n_active_cells() == 88, ExcInternalError());
  }
}


int
main(int argc, char *argv[])
{
  initlog();
#ifdef DEAL_II_WITH_MPI
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
#else
  (void)argc;
  (void)argv;
#endif

  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
