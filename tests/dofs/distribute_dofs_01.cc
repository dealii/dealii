// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_nedelec.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Check DoFHandler::distribute_dofs with multiple calls to reinit

void
test()
{
  dealii::Triangulation<3> tria;
  dealii::GridGenerator::hyper_cube(tria, -1, 1);

  dealii::FE_Nedelec<3> fe_nedelec(0);

  dealii::DoFHandler<3> dh_nedelec;

  dh_nedelec.reinit(tria);
  dh_nedelec.distribute_dofs(fe_nedelec);

  dh_nedelec.reinit(tria);
  dh_nedelec.distribute_dofs(fe_nedelec);

  for (const auto &cell : dh_nedelec.active_cell_iterators())
    for (unsigned int i = 0; i < dealii::GeometryInfo<3>::lines_per_cell; ++i)
      cell->line(i)->dof_index(0);

  deallog << "ok" << std::endl;
}

int
main()
{
  initlog();

  test();

  return 0;
}
