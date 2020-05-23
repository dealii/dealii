// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
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



// check GeometryInfo::d_linear_shape_function_gradient

#include <deal.II/base/geometry_info.h>

#include "../tests.h"



template <int dim>
void
test()
{
  deallog << "Checking in " << dim << "d" << std::endl;

  // check phi_i(v_j) = delta_{ij}
  for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
    for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
      {
        const Tensor<1, dim> phi_i_grad =
          GeometryInfo<dim>::d_linear_shape_function_gradient(
            GeometryInfo<dim>::unit_cell_vertex(v), i);

        deallog << phi_i_grad << std::endl;
      }

  // check that
  //    sum_i phi_i(x) == const
  // at all points by verifying that the
  // gradient of the sum of shape functions
  // is zero. do so at every vertex, and then
  // at the center
  for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
    {
      Tensor<1, dim> s;
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        s += GeometryInfo<dim>::d_linear_shape_function_gradient(
          GeometryInfo<dim>::unit_cell_vertex(v), i);
      AssertThrow(s.norm() == 0, ExcInternalError());

      deallog << "Sum of shape functions: " << s << std::endl;
    }
  {
    Point<dim> center;
    for (unsigned int i = 0; i < dim; ++i)
      center[i] = 0.5;

    Tensor<1, dim> s;
    for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
      s += GeometryInfo<dim>::d_linear_shape_function_gradient(center, i);
    AssertThrow(s.norm() == 0, ExcInternalError());

    deallog << "Sum of shape functions: " << s << std::endl;
  }
}


int
main()
{
  initlog();

  test<1>();
  test<2>();
  test<3>();

  return 0;
}
