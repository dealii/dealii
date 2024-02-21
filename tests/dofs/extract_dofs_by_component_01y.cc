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
//
// this is a variant of the _01x test that shows that we were doing something
// wrong all along even on mainline but that nobody apparently ever realized
// this. the bug is that internal::extract_dofs_by_component got things wrong
// when a non-primitive element was part of an FESystem but it was
// accidentally correct whenever there was only one such element; this test
// verifies the same with two non-primitive elements in one FESystem


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
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

  FESystem<dim>   element(FE_Nedelec<dim>(0), 2);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);

  // use a mask that only has the one
  // component set and cycle over
  // which one that is
  for (unsigned int comp = 0; comp < element.n_components(); ++comp)
    {
      std::vector<bool> component_mask(element.n_components(), false);
      component_mask[comp] = true;

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
}
