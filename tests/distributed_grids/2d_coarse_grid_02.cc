// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2020 by the deal.II authors
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
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

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
