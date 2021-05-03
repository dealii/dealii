// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
