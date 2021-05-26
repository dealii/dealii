// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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
