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



// Test interaction with p4est with a complicated 2d grid read from file. the
// grid describes a cross-section of an airfoil with flaps at the front and
// back. it has some 30,000 cells

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

#include "coarse_grid_common.h"



template <int dim>
void
test(std::ostream & /*out*/)
{
  parallel::distributed::Triangulation<dim> tr(
    MPI_COMM_WORLD,
    Triangulation<dim>::none,
    parallel::distributed::Triangulation<dim>::communicate_vertices_to_p4est);

  GridIn<dim> gi;
  gi.attach_triangulation(tr);
  std::ifstream in(SOURCE_DIR "/../grid/grid_in_02/2d.xda");
  try
    {
      gi.read_xda(in);
    }
  catch (const typename Triangulation<dim>::DistortedCellList &distorted_cells)
    {
      // ignore distorted cells
      deallog << distorted_cells.distorted_cells.size()
              << " distorted cells after creating mesh." << std::endl;
    }

  write_vtk(tr, "1");
}


int
main(int argc, char *argv[])
{
  initlog();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
}
