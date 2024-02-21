// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test conversion from deal.II Point to VTK double array

#include <deal.II/base/point.h>

#include <deal.II/vtk/utilities.h>

#include "../tests.h"

int
main()
{
  initlog();

  // Test conversion from deal.II Point to VTK array.
  {
    const Point<3> point_dealii(1.0, 2.0, 3.0);
    const auto point_vtk = VTKWrappers::dealii_point_to_vtk_array(point_dealii);
    deallog << "VTK Point: " << point_vtk->GetComponent(0, 0) << ", "
            << point_vtk->GetComponent(0, 1) << ", "
            << point_vtk->GetComponent(0, 2) << std::endl;
  }

  return 0;
}
