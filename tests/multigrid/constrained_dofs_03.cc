// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// check mg constrained dofs for primitive and non-primitive FiniteElements

#include <deal.II/base/function.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_tools.h>

#include <algorithm>

#include "../tests.h"



template <int dim>
void
check_fe(FiniteElement<dim> &fe, ComponentMask &component_mask)
{
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global(1);

  DoFHandler<dim> dofh(tr);
  dofh.distribute_dofs(fe);
  dofh.distribute_mg_dofs();

  MGConstrainedDoFs                  mg_constrained_dofs;
  const std::set<types::boundary_id> boundary_indicators = {0};
  mg_constrained_dofs.initialize(dofh);
  mg_constrained_dofs.make_zero_boundary_constraints(dofh,
                                                     boundary_indicators,
                                                     component_mask);

  const unsigned int n_levels = tr.n_global_levels();
  for (unsigned int level = 0; level < n_levels; ++level)
    {
      deallog << "Level " << level << ':' << std::endl;
      IndexSet boundary_indices =
        mg_constrained_dofs.get_boundary_indices(level);
      boundary_indices.print(deallog);
    }
}

template <int dim>
void
check()
{
  // All primitive
  {
    FE_Q<dim>     q1(2);
    FE_Q<dim>     q2(1);
    FESystem<dim> s1(q1, dim, q2, 1);

    // All selected
    ComponentMask component_mask1(dim + 1, true);
    // Partially selected
    ComponentMask component_mask2(dim + 1, true);
    component_mask2.set(dim, false);

    check_fe(s1, component_mask1);
    check_fe(s1, component_mask2);
  }

  // Non-primitive
  {
    FE_RaviartThomas<dim> q1(1);
    FE_DGQ<dim>           q2(1);
    FESystem<dim>         s1(q1, 1, q2, 1);

    // All selected
    ComponentMask component_mask1(dim + 1, true);
    // Partially selected
    ComponentMask component_mask2(dim + 1, true);
    component_mask2.set(dim, false);

    check_fe(s1, component_mask1);
    check_fe(s1, component_mask2);
  }
}

int
main(int argc, char *argv[])
{
  initlog();

  check<2>();
  check<3>();
}
