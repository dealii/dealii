// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// the DoFTools::count_dofs_per_{component,block} functions resized the output
// array to fe.n_components or fe.n_blocks even if a grouping argument was
// given. this would appear wrong and can lead to all sorts of interesting
// behavior if not caught by an assertion early enough


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <string>

#include "../tests.h"



void
print(const std::vector<types::global_dof_index> &v)
{
  deallog << v.size();
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << ' ' << v[i];
  deallog << std::endl;
}



template <int dim>
void
check()
{
  // create tria and dofhandler
  // objects. set different boundary
  // and sub-domain ids
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0., 1.);
  tria.refine_global(1);
  for (int i = 0; i < 2; ++i)
    {
      tria.begin_active()->set_refine_flag();
      tria.execute_coarsening_and_refinement();
    }

  // Taylor-Hood element
  FESystem<dim>   fe(FE_Q<dim>(1), dim, FE_DGQ<dim>(0), 1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // no grouping
  {
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_component(dof_handler);
    print(dpc);
  }

  {
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_block(dof_handler);
    print(dpc);
  }


  // grouping into less groups than
  // components
  {
    std::vector<unsigned int> group(dim + 1, 0);
    group[dim] = 1;
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_component(dof_handler, false, group);
    Assert(dpc.size() == 2, ExcInternalError());
    print(dpc);
  }

  {
    std::vector<unsigned int> group(dim + 1, 0);
    group[dim] = 1;
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_block(dof_handler, group);
    Assert(dpc.size() == 2, ExcInternalError());
    print(dpc);
  }

  // grouping into more groups than
  // components
  {
    std::vector<unsigned int> group(dim + 1, 2 * dim);
    group[dim] = 0;
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_component(dof_handler, false, group);
    Assert(dpc.size() == 2 * dim + 1, ExcInternalError());
    print(dpc);
  }

  {
    std::vector<unsigned int> group(dim + 1, 2 * dim);
    group[dim] = 0;
    const std::vector<types::global_dof_index> dpc =
      DoFTools::count_dofs_per_fe_block(dof_handler, group);
    Assert(dpc.size() == 2 * dim + 1, ExcInternalError());
    print(dpc);
  }
}



int
main()
{
  initlog();
  deallog << std::setprecision(2);

  check<1>();
  check<2>();
  check<3>();
}
