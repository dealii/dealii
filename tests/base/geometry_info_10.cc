// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test GeometryInfo<dim>::unit_tangential_vectors

#include <deal.II/base/geometry_info.h>

#include "../tests.h"


template <int dim>
void
test()
{
  for (const unsigned int i : GeometryInfo<dim>::face_indices())
    {
      const Tensor<1, dim> unit_normal_vector =
        GeometryInfo<dim>::unit_normal_vector[i];
      deallog << "Direction " << i << " = " << unit_normal_vector << std::endl;

      deallog << "Tangential vectors: ";
      for (const auto t : GeometryInfo<dim>::unit_tangential_vectors[i])
        deallog << t << " ; ";
      deallog << std::endl;
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
  deallog.push("4d");
  test<4>();
  deallog.pop();
  return 0;
}
