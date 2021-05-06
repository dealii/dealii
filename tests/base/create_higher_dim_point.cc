// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2019 by the deal.II authors
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


#include <deal.II/base/function_restriction.h>

#include "../tests.h"

namespace
{
  using namespace dealii;

  // Test the function create_higher_dim_point by making sure that we
  // can construct the point (x, y, z),
  const double x = 1, y = 2, z = 3;
  // for each possible component. That is, we test that we get the following
  // results:
  // {point, component, coordinate} ->  result
  //        {(y, z), 0, x}          -> (x, y, z)
  //        {(z, x), 1, y}          -> (x, y, z)
  //        {(y, z), 2, x}          -> (x, y, z)
  //
  // The same thing is done for the 2D case.
  void
  test3D()
  {
    const int dim = 3;
    deallog << "dim=" << dim << std::endl;
    {
      const unsigned int   component = 0;
      const Point<dim - 1> low_dim_point(y, z);
      const Point<dim>     point =
        internal::create_higher_dim_point(low_dim_point, component, x);
      deallog << "Point: " << point << std::endl;
    }
    {
      const unsigned int component = 1;
      // The order (z,x) here follows convention in
      // coordinate_to_one_dim_higher.
      const Point<dim - 1> low_dim_point(z, x);
      const Point<dim>     point =
        internal::create_higher_dim_point(low_dim_point, component, y);
      deallog << "Point: " << point << std::endl;
    }
    {
      const unsigned int component = 2;
      Point<dim - 1>     low_dim_point(x, y);
      Point<dim>         point =
        internal::create_higher_dim_point(low_dim_point, component, z);
      deallog << "Point: " << point << std::endl;
    }
  }



  void
  test2D()
  {
    const int dim = 2;
    deallog << "dim=" << dim << std::endl;
    {
      const unsigned int   component = 0;
      const Point<dim - 1> low_dim_point(y);
      const Point<dim>     point =
        internal::create_higher_dim_point(low_dim_point, component, x);
      deallog << "Point: " << point << std::endl;
    }
    {
      const unsigned int   component = 1;
      const Point<dim - 1> low_dim_point(x);
      const Point<dim>     point =
        internal::create_higher_dim_point(low_dim_point, component, y);
      deallog << "Point: " << point << std::endl;
    }
  }



  void
  test1D()
  {
    const int dim = 1;
    deallog << "dim=" << dim << std::endl;
    const unsigned int   component = 0;
    const Point<dim - 1> low_dim_point;
    const Point<dim>     point =
      internal::create_higher_dim_point(low_dim_point, component, x);
    deallog << "Point: " << point << std::endl;
  }
} // namespace



int
main()
{
  initlog();
  test1D();
  test2D();
  test3D();
}
