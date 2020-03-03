// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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
