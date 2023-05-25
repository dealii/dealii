// ---------------------------------------------------------------------
//
// Copyright (C) 2023 by the deal.II authors
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
