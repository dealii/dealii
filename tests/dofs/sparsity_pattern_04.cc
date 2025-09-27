// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2007 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// use the version of DoFTools::make_sparsity_pattern that takes two
// DoFHandler arguments for two DoFHandlers that are actually from different
// meshes (though with the same base)
//
// like sparsity_pattern_03, but use different finite elements and use a mesh
// with a different number of cells


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"



template <int dim>
void
check()
{
  // create two different triangulations
  Triangulation<dim> triangulation_1;
  if (dim == 2)
    GridGenerator::hyper_ball(triangulation_1, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(triangulation_1, -1, 1);
  triangulation_1.refine_global(1);
  triangulation_1.begin_active()->set_refine_flag();
  triangulation_1.execute_coarsening_and_refinement();
  triangulation_1.begin_active(2)->set_refine_flag();
  triangulation_1.execute_coarsening_and_refinement();
  if (dim == 1)
    triangulation_1.refine_global(2);


  Triangulation<dim> triangulation_2;
  if (dim == 2)
    GridGenerator::hyper_ball(triangulation_2, Point<dim>(), 1);
  else
    GridGenerator::hyper_cube(triangulation_2, -1, 1);
  triangulation_2.refine_global(1);
  (std::next(triangulation_2.begin_active()))->set_refine_flag();
  triangulation_2.execute_coarsening_and_refinement();
  (std::next(triangulation_2.begin_active(2)))->set_refine_flag();
  triangulation_2.execute_coarsening_and_refinement();
  if (dim == 1)
    triangulation_2.refine_global(2);
  triangulation_2.refine_global(1);



  FESystem<dim>   element_1(FE_Q<dim>(1), 2, FE_Q<dim>(2), 1);
  FESystem<dim>   element_2(FE_Q<dim>(3), 1, FE_DGQ<dim>(0), 2);
  DoFHandler<dim> dof_1(triangulation_1);
  DoFHandler<dim> dof_2(triangulation_2);
  dof_1.distribute_dofs(element_1);
  dof_2.distribute_dofs(element_2);

  SparsityPattern sparsity(dof_1.n_dofs(),
                           dof_2.n_dofs(),
                           std::max(dof_1.n_dofs(), dof_2.n_dofs()));
  DoFTools::make_sparsity_pattern(dof_1, dof_2, sparsity);
  sparsity.compress();

  sparsity.print(deallog.get_file_stream());
}



int
main()
{
  initlog();
  deallog << std::setprecision(2) << std::fixed;

  deallog.push("1d");
  check<1>();
  deallog.pop();
  deallog.push("2d");
  check<2>();
  deallog.pop();
  deallog.push("3d");
  check<3>();
  deallog.pop();
}
