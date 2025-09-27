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



// Test interaction with p4est with a simple grid in 2d. here, we
// check that if we refine a square once, and then one of the children
// once again, that we get 7 cells
//
// at the time of writing this test, the results for this testcase
// were erratic and apparently non-deterministic. the actual cause was
// an uninitialized variable, fixed in revision 16414.

#include <deal.II/base/tensor.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
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
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

      GridGenerator::hyper_cube(tr);

      deallog << i << ' ' << tr.n_active_cells() << std::endl;

      tr.refine_global(1);

      deallog << i << ' ' << tr.n_active_cells() << std::endl;

      Assert(tr.n_active_cells() == 4, ExcInternalError());

      typename Triangulation<dim>::active_cell_iterator cell =
        tr.begin_active();
      std::advance(cell, i);
      cell->set_refine_flag();
      tr.execute_coarsening_and_refinement();

      deallog << i << ' ' << tr.n_active_cells() << std::endl;

      Assert(tr.n_active_cells() == 7, ExcInternalError());
    }
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
