// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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


// Test optional arguments of MGLevelObject.

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>

#include "../tests.h"

using namespace dealii;

int
main()
{
  initlog();

  Triangulation<2> tria;

  MGLevelObject<DoFHandler<2>> dof_handlers(0, 3, tria);

  dof_handlers.resize(0, 5, tria);

  deallog << "OK!" << std::endl;
}
