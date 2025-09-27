// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// document a bug in hp::DoFHandler and copy_triangulation

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test1()
{
  Triangulation<dim> tr1;
  DoFHandler<dim>    dofh(tr1);

  Triangulation<dim> tr2;
  GridGenerator::hyper_ball(tr2);
  tr2.refine_global(2);
  tr1.copy_triangulation(tr2);

  dofh.begin_active()->set_active_fe_index(0);
}

int
main()
{
  initlog();

  test1<2>();
  return 0;
}
