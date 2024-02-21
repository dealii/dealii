// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// verify that iterators can be copied. when you just copy the
// iterator itself, everything is alright: the copy operator of the
// iterator calls Accessor::copy_from.
//
// the test was originally written because at that time it was
// possible to do something like *it1 = *it2 for DoF iterators --
// likely not what the author intended since it does not copy the
// cell2 into the cell1, only the accessor. furthermore, the operator=
// was wrongly implemented, so it was removed in the process


#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  typename DoFHandler<dim>::active_cell_iterator cell =
    dof_handler.begin_active();
  ++cell;

  typename DoFHandler<dim>::face_iterator face =
    dof_handler.begin_active()->face(0);
  face = cell->face(0);
  deallog << cell->face(0) << ' ' << face << std::endl;
  Assert(cell->face(0) == face, ExcInternalError());
  Assert(!(cell->face(0) != face), ExcInternalError());
}



int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
