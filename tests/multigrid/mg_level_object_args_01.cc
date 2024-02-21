// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Test optional arguments of MGLevelObject.

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/grid/tria.h>

#include "../tests.h"


int
main()
{
  initlog();

  Triangulation<2> tria;

  MGLevelObject<DoFHandler<2>> dof_handlers(0, 3, tria);

  dof_handlers.resize(0, 5, tria);

  deallog << "OK!" << std::endl;
}
