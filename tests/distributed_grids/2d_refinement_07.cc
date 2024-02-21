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



// Test that we abort if you refine further than p4est::max_level

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
  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  for (unsigned int i = 0; i < 30; ++i)
    {
      deallog << "cells: " << tr.n_active_cells() << " level:" << tr.n_levels()
              << std::endl;

      typename parallel::distributed::Triangulation<dim>::cell_iterator it;
      it = tr.begin_active();
      while (it->level() < static_cast<int>(tr.n_levels() - 1))
        ++it;


      it->set_refine_flag();
      try
        {
          tr.execute_coarsening_and_refinement();
        }
      catch (ExceptionBase &e)
        {
          deallog << e.get_exc_name() << std::endl;
        }
    }
}


int
main(int argc, char *argv[])
{
  initlog();
  deal_II_exceptions::disable_abort_on_exception();
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  deallog.push("2d");
  test<2>(deallog.get_file_stream());
  deallog.pop();
}
