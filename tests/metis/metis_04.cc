// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check GridTools::partition_triangulation with too many subdomains

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include "../tests.h"


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(1);
  {
    triangulation.begin_active()->set_refine_flag();
    triangulation.execute_coarsening_and_refinement();
  }

  for (unsigned int procs = 1; procs < 20; ++procs)
    {
      deallog << " procs = " << procs << std::endl;
      GridTools::partition_triangulation(procs, triangulation);
      {
        typename Triangulation<dim>::active_cell_iterator cell =
          triangulation.begin_active();
        for (; cell != triangulation.end(); ++cell)
          {
            deallog << cell->subdomain_id() << ' ';
          }
        deallog << std::endl;
      }
    }
  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  test<2>();
}
