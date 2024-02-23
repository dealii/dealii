// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Test Tria<2,3>

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim, int spacedim>
void
test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim, spacedim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim, spacedim>::none,
    parallel::distributed::Triangulation<dim, spacedim>::
      communicate_vertices_to_p4est);

  GridGenerator::torus(tr, 1, 0.2);
  tr.reset_all_manifolds();
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  write_vtk(tr, "file");
}


int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  initlog();

  deallog.push("2-3");
  test<2, 3>(deallog.get_file_stream());
  deallog.pop();
}
