// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

// Check that VectorTools::interpolate correctly recovers a
// constant vector field for H1, Hdiv and Hcurl conforming elements.

#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>

#include "../tests.h"

#include "generalized_interpolation.h"

int
main()
{
  initlog();

  test<2>(FE_Q<2>(1), ConstantFunction<2>(1.0, 1), 1, false, true);
  test<2>(FE_Q<2>(2), ConstantFunction<2>(1.0, 1), 2, false, true);
  test<3>(FE_Q<3>(3), ConstantFunction<3>(1.0, 1), 3, false, true);

  test<2>(FE_Q<2>(1), ConstantFunction<2>(1.0, 1), 1, true, true);
  test<2>(FE_Q<2>(2), ConstantFunction<2>(1.0, 1), 2, true, true);
  test<3>(FE_Q<3>(3), ConstantFunction<3>(1.0, 1), 3, true, true);

  test<2>(FE_RaviartThomas<2>(0), ConstantFunction<2>(1.0, 2), 1, false, true);
  test<2>(FE_RaviartThomas<2>(1), ConstantFunction<2>(1.0, 2), 2, false, true);
  test<2>(FE_RaviartThomas<2>(2), ConstantFunction<2>(1.0, 2), 3, false, true);
  test<3>(FE_RaviartThomas<3>(0), ConstantFunction<3>(1.0, 3), 1, false, true);
  test<3>(FE_RaviartThomas<3>(1), ConstantFunction<3>(1.0, 3), 2, false, true);

  test<2>(FE_RaviartThomas<2>(0), ConstantFunction<2>(1.0, 2), 1, true, true);
  test<2>(FE_RaviartThomas<2>(1), ConstantFunction<2>(1.0, 2), 2, true, true);
  test<2>(FE_RaviartThomas<2>(2), ConstantFunction<2>(1.0, 2), 3, true, true);
  // lowest order RT in 3D does not contain constant 1 function on a
  // distorted mesh.
  test<3>(FE_RaviartThomas<3>(1), ConstantFunction<3>(1.0, 3), 2, true, true);

  test<2>(FE_Nedelec<2>(0), ConstantFunction<2>(1.0, 2), 1, false, true);
  test<2>(FE_Nedelec<2>(1), ConstantFunction<2>(1.0, 2), 2, false, true);
  test<2>(FE_Nedelec<2>(2), ConstantFunction<2>(1.0, 2), 3, false, true);
  test<2>(FE_Nedelec<2>(0), ConstantFunction<2>(1.0, 2), 1, true, true);
  test<2>(FE_Nedelec<2>(1), ConstantFunction<2>(1.0, 2), 2, true, true);
  test<2>(FE_Nedelec<2>(2), ConstantFunction<2>(1.0, 2), 3, true, true);
}
