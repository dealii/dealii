// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2022 by the deal.II authors
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

#include <deal.II/distributed/p4est_wrappers.h>
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
  const unsigned int max_level = internal::p4est::functions<dim>::max_level;

  deallog << "The maximal level of p4est refinements is " << max_level
          << std::endl;

  parallel::distributed::Triangulation<dim> tr(MPI_COMM_WORLD);

  GridGenerator::hyper_cube(tr);
  // The first refinement level is "1". So in order to refine past the
  // limit we have to attempt to refine max_level times:
  for (unsigned int i = 0; i < max_level; ++i)
    {
      deallog << "cells: " << tr.n_active_cells() << " level:" << tr.n_levels()
              << std::endl;

      typename parallel::distributed::Triangulation<dim>::cell_iterator it;
      it = tr.begin_active();
      while (it->level() < static_cast<int>(i))
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

  deallog.push("3d");
  test<3>(deallog.get_file_stream());
  deallog.pop();
}
