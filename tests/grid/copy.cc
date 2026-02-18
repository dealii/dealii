// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

// Test Triangulation<dim, spacedim>::operator=(Triangulation<dim, spacedim>
// &&tria)

dealii::Triangulation<2>
createDummyTria()
{
  dealii::Triangulation<2, 2> patch;
  dealii::GridGenerator::hyper_cube<2>(patch);
  dealii::Triangulation<2> tria;
  tria.copy_triangulation(patch); // here the policies are copied (just a check)
  return patch;                   // patch has the policies
}

int
main()
{
  initlog();

  dealii::Triangulation<2> tria;
  tria.copy_triangulation(createDummyTria());

  tria.refine_global();

  deallog << "OK" << std::endl;
}
