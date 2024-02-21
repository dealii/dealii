// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test internal::extract_dofs_by_component for some corner cases that
// I was unsure about when refactoring some code in there
//
// this particular test checks the call path to
// internal::extract_dofs_by_component from
// DoFTools::count_dofs_per_fe_component with argument only_once=false


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);

  FESystem<dim>   element(FE_Q<dim>(1),
                        1,
                        FE_RaviartThomas<dim>(0),
                        1,
                        FE_Q<dim>(1),
                        1,
                        FE_Nedelec<dim>(0),
                        1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  const std::vector<types::global_dof_index> count =
    DoFTools::count_dofs_per_fe_component(dof, false);

  for (unsigned int d = 0; d < count.size(); ++d)
    deallog << count[d] << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
