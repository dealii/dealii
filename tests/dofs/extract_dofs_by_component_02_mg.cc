// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2022 by the deal.II authors
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
// the BlockMask argument


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
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);
  tr.begin_active()->set_refine_flag();
  tr.execute_coarsening_and_refinement();

  FESystem<dim>   element(FE_Q<dim>(2), 1, FE_Nedelec<dim>(0), 1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(element);
  dof.distribute_mg_dofs();

  // try all possible block
  // masks, which we encode as bit
  // strings
  for (unsigned int int_mask = 0; int_mask < (1U << element.n_blocks());
       ++int_mask)
    {
      std::vector<bool> component_mask(element.n_blocks());
      for (unsigned int c = 0; c < element.n_blocks(); ++c)
        component_mask[c] = (int_mask & (1 << c));

      for (unsigned int level = 0; level < tr.n_levels(); ++level)
        {
          deallog << "level=" << level << std::endl;

          std::vector<bool> dofs(dof.n_dofs(level));
          DoFTools::extract_level_dofs(level,
                                       dof,
                                       BlockMask(component_mask),
                                       dofs);

          for (unsigned int d = 0; d < dofs.size(); ++d)
            deallog << dofs[d];
          deallog << std::endl;
        }
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
