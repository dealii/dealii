// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check DoFHandler::has_level_dofs and DoFHandler::has_active_dofs

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <algorithm>

#include "../tests.h"



template <int dim>
void
check()
{
  FE_Q<dim> fe(1);

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  // check the two functions mentioned above
  // in their natural order
  {
    DoFHandler<dim> dof(tr);
    deallog << "check " << dim << " before distribute " << dof.has_active_dofs()
            << ' ' << dof.has_level_dofs() << std::endl;

    dof.distribute_dofs(fe);
    deallog << "check " << dim << " after  distribute " << dof.has_active_dofs()
            << ' ' << dof.has_level_dofs() << std::endl;


    dof.distribute_dofs(fe);
    dof.distribute_mg_dofs();
    deallog << "check " << dim << " level  distribute " << dof.has_active_dofs()
            << ' ' << dof.has_level_dofs() << std::endl;

    dof.clear();
    deallog << "check " << dim << " after  clear      " << dof.has_active_dofs()
            << ' ' << dof.has_level_dofs() << std::endl;
  }
}

int
main()
{
  initlog(__FILE__);
  check<1>();
  check<2>();
  check<3>();
}
