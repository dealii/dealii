// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// Check ChartManifold::get_intermediate_point

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/grid/manifold_lib.h>

template <int dim>
void
print_intermediate_point(const Manifold<dim>& manifold,
                         const std::string&   manifold_name,
                         const Point<dim>&    p1,
                         const Point<dim>&    p2,
                         const double         weight)
{
  const std::vector<Point<dim>> points({p1, p2});
  const std::vector<double>     weights({1 - weight, weight});
  deallog.precision(3);
  deallog << manifold_name << " between points [" << p1 << "] and [" << p2
          << "] with weight " << weight << std::endl;
  deallog.precision(12);
  deallog << "Intermediate point: "
          << manifold.get_intermediate_point(p1, p2, weight) << std::endl
          << "get_new_point:      " << manifold.get_new_point(points, weights)
          << std::endl;
}

int
main()
{
  initlog();

  const PolarManifold<2> polar;
  {
    Point<2> p1(1., 0.1);
    Point<2> p2(0.1, 1);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.1);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.5);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.9);
  }
  {
    Point<2> p1(1., 0.1);
    Point<2> p2(1., -0.1);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.1);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.5);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.9);
  }
  {
    Point<2> p1(1., 0.1);
    Point<2> p2(-1., -0.099);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.5);
  }
  {
    Point<2> p1(1., 0.1);
    Point<2> p2(-1., -0.101);
    print_intermediate_point(polar, "PolarManifold", p1, p2, 0.5);
  }

  const CylindricalManifold<3> cylindrical(2);
  {
    Point<3> p1(1., 0.1, 0.2);
    Point<3> p2(0.1, 1, 0.2);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.1);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.5);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.9);
  }
  {
    Point<3> p1(1., 0.1, 0.2);
    Point<3> p2(-1, -0.099, 0.2);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.5);
  }
  {
    Point<3> p1(1., 0.1, 0.2);
    Point<3> p2(-1, -0.101, 0.2);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.5);
  }
  {
    Point<3> p1(1., 0.1, 0.2);
    Point<3> p2(-1, -0.101, 0.5);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.1);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.5);
    print_intermediate_point(cylindrical, "CylindricalManifold", p1, p2, 0.9);
  }

  return 0;
}
