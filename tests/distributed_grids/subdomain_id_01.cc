// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check (level)subdomain_id(). Should be =0.

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include "../tests.h"



template <class TRIA>
void
check(TRIA &tr)
{
  typename TRIA::cell_iterator cell = tr.begin(), endc = tr.end();

  for (; cell != endc; ++cell)
    {
      deallog << cell->level_subdomain_id() << ' ';
      try
        {
          deallog << cell->subdomain_id();
        }
      catch (...)
        {
          deallog << '.';
        }
      deallog << std::endl;
    }

  deallog << "OK" << std::endl;
}


int
main(int argc, char *argv[])
{
  deal_II_exceptions::disable_abort_on_exception();

  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  initlog();

  Triangulation<2>                        tria;
  parallel::distributed::Triangulation<2> tria2(MPI_COMM_WORLD);
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);
  GridGenerator::hyper_cube(tria2);
  tria2.refine_global(2);
  check(tria);
  check(tria2);
}
