// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test for equality of accessor objects

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << dim << 'd' << std::endl;

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  tria.refine_global(1);

  deallog << tria.begin_active() << std::endl;

  // test a few comparison operators
  {
    const typename Triangulation<dim>::active_cell_iterator cell =
      tria.begin_active();
    AssertThrow(cell == cell, ExcInternalError());
    AssertThrow(!(cell != cell), ExcInternalError());
    AssertThrow(!(cell < cell), ExcInternalError());
  }

  // same with non-active iterators
  {
    const typename Triangulation<dim>::cell_iterator cell = tria.begin();
    AssertThrow(cell == cell, ExcInternalError());
    AssertThrow(!(cell != cell), ExcInternalError());
    AssertThrow(!(cell < cell), ExcInternalError());
  }

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  // now do the same tests with the
  // DoFHandler iterators
  {
    const typename DoFHandler<dim>::active_cell_iterator cell =
      dof_handler.begin_active();
    AssertThrow(cell == cell, ExcInternalError());
    AssertThrow(!(cell != cell), ExcInternalError());
    AssertThrow(!(cell < cell), ExcInternalError());
  }
  {
    const typename DoFHandler<dim>::cell_iterator cell = dof_handler.begin();
    AssertThrow(cell == cell, ExcInternalError());
    AssertThrow(!(cell != cell), ExcInternalError());
    AssertThrow(!(cell < cell), ExcInternalError());
  }

  // finally check that two iterators
  // pointing into different DoFHandlers
  // can't be the same
  //
  // this test accidentally failed before a
  // fix was checked in on 2005-08-08 that
  // completely forbids comparing iterators
  // into different DoFHandlers at all
  {
    DoFHandler<dim> dof_handler2(tria);
    dof_handler2.distribute_dofs(fe);
    try
      {
        bool is_same =
          (dof_handler.begin_active() != dof_handler2.begin_active());
        (void)is_same;
      }
    catch (ExceptionBase &e)
      {
        deallog << e.get_exc_name() << std::endl;
      }
  }
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  initlog();

  test<1>();
  test<2>();
  test<3>();
}
