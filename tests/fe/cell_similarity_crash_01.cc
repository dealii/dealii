// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2011 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// The FEValues class stores a copy of the current cell so that when we move
// on to the next cell we can compare the two for shape similarity and decide
// which parts of the data we need to re-compute. But this runs into trouble
// if the lifetime of FEValues object extends beyond a refinement or
// coarsening step in the triangulation since in that case the stored copy of
// the cell previously worked on may no longer be valid. The solution is to
// invalidate the stored cell once mesh refinement happens.

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"



template <int dim>
void
test()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);

  FE_Q<dim>       fe(1);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  const QGauss<dim> quadrature(2);
  FEValues<dim>     fe_values(fe, quadrature, update_values);

  // initialize FEValues with the first cell
  fe_values.reinit(dof.begin_active());

  // then make this first cell invalid by
  // coarsening it away
  for (unsigned int c = 0; c < GeometryInfo<dim>::max_children_per_cell; ++c)
    tr.begin(1)->child(c)->set_coarsen_flag();
  tr.execute_coarsening_and_refinement();
  dof.distribute_dofs(fe);

  // initialize FEValues with what is now the
  // first cell. this used to fail (see the
  // reason given at the top)
  fe_values.reinit(dof.begin_active());

  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(4);

  test<1>();
  test<2>();
  test<3>();
}
