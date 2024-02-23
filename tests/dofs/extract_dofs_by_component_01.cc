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
// internal::extract_dofs_by_component from DoFTools::extract_dofs via
// the component_select flag


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

  // try all possible component
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask = 0; int_mask < (1U << element.n_components());
       ++int_mask)
    {
      std::vector<bool> component_mask(element.n_components());
      for (unsigned int c = 0; c < element.n_components(); ++c)
        component_mask[c] = (int_mask & (1 << c));

      const IndexSet dofs =
        DoFTools::extract_dofs(dof, ComponentMask(component_mask));

      for (unsigned int d = 0; d < dof.n_dofs(); ++d)
        deallog << dofs.is_element(d);
      deallog << std::endl;
    }
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
