// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GeometryInfo<dim>::unit_normal_vector

#include <deal.II/base/geometry_info.h>

#include "../tests.h"

const double eps = 1e-10;
DeclException3(ExcWrongValue,
               double,
               double,
               double,
               << arg1 << " != " << arg2 << " with delta = " << arg3);

template <int dim>
void
test()
{
  for (const unsigned int i : GeometryInfo<dim>::face_indices())
    {
      const Tensor<1, dim> unit_normal_vector =
        GeometryInfo<dim>::unit_normal_vector[i];
      deallog << "Direction " << i << " = " << unit_normal_vector << std::endl;
      // check consistency with other arrays within GeometryInfo
      for (unsigned int j = 0; j < dim; ++j)
        {
          const double obtained = unit_normal_vector[j];
          const double expected =
            (j == GeometryInfo<dim>::unit_normal_direction[i]) ?
              static_cast<double>(
                GeometryInfo<dim>::unit_normal_orientation[i]) :
              0;
          AssertThrow(obtained == expected,
                      ExcWrongValue(obtained, expected, obtained - expected));
        }
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  deallog.push("1d");
  test<1>();
  deallog.pop();
  deallog.push("2d");
  test<2>();
  deallog.pop();
  deallog.push("3d");
  test<3>();
  deallog.pop();
  return 0;
}
