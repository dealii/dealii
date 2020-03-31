// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2018 by the deal.II Authors
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

// Test all member functions of the HyperCube<dim> class. Each function is
// tested in the function named test_{member_function_name}.

#include <deal.II/base/geometry_info.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"


using namespace dealii;
using NonMatching::internal::QuadratureGeneratorImplementation::Hypercube;

// Print the center and the side length of the incoming hypercube, which are
// its degrees of freedom.
template <int dim>
void
print_cube(const Hypercube<dim> &hyper_cube)
{
  deallog << "center = " << hyper_cube.center() << ", ";
  deallog << "side = " << hyper_cube.side_length() << std::endl;
}



// Construct the cross section orthogonal to all axes and print them.
template <int dim>
void
test_cross_section()
{
  deallog << "test_cross_section" << std::endl;

  const double width = 1;
  Point<dim>   center;
  for (unsigned int d = 0; d < dim; ++d)
    center[d] = d;

  const Hypercube<dim> hypercube(center, width);

  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << "orthogonal to " << d << std::endl;
      const Hypercube<dim - 1> down = hypercube.cross_section(d);
      print_cube(down);
    }
}



// 1D-specialization of the above test. We can't take a cross section in 1D, so
// we do nothing.
template <>
void
test_cross_section<1>()
{}



// Construct all possible children of a hypercube and print them.
template <int dim>
void
test_child()
{
  deallog << "test_child" << std::endl;

  const double width = 4;
  Point<dim>   center;

  const Hypercube<dim> hypercube(center, width);

  // Take the cross section in all directions.
  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      deallog << "child " << i << std::endl;
      const Hypercube<dim> child = hypercube.child(i);
      print_cube(child);
    }
}



// Comptue and print the volume of a hypercube.
template <int dim>
void
test_volume()
{
  deallog << "test_volume" << std::endl;

  const double width = 2;
  Point<dim>   center;

  const Hypercube<dim> hypercube(center, width);
  deallog << "volume = " << hypercube.volume() << std::endl;
}



// Print the result of the lower_bound and upper_bound functions in all
// coordinate directions.
template <int dim>
void
test_lower_upper_bound()
{
  deallog << "test_lower_upper_bound" << std::endl;

  const double width = 2;
  Point<dim>   center;
  for (unsigned int d = 0; d < dim; ++d)
    center[d] = d + 1;

  const Hypercube<dim> hypercube(center, width);
  for (int d = 0; d < dim; ++d)
    {
      deallog << "[" << hypercube.lower_bound(d) << ", "
              << hypercube.upper_bound(d) << "]";
      if (dim - 1 != d)
        {
          deallog << "x";
        }
    }
}



// Extract and print all the vertices of a hypercube.
template <int dim>
void
test_vertex()
{
  deallog << "test_vertex" << std::endl;

  const Hypercube<dim> hypercube;
  for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      deallog << "vertex " << i << " = " << hypercube.vertex(i) << std::endl;
    }
}



// Print the results of the bounds function for all coordinate directions.
template <int dim>
void
test_bounds()
{
  deallog << "test_bounds" << std::endl;

  const double width = 2;
  Point<dim>   center;
  for (unsigned int d = 0; d < dim; ++d)
    center[d] = d + 1;

  const Hypercube<dim> hypercube(center, width);

  for (unsigned int i = 0; i < dim; ++i)
    {
      const Hypercube<1> bounds = hypercube.bounds(i);
      deallog << "Bounds in direction " << i << std::endl;
      print_cube(bounds);
    }
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;

  test_cross_section<dim>();
  deallog << std::endl;

  test_child<dim>();
  deallog << std::endl;

  test_volume<dim>();
  deallog << std::endl;

  test_lower_upper_bound<dim>();
  deallog << std::endl;

  test_vertex<dim>();
  deallog << std::endl;

  test_bounds<dim>();
  deallog << std::endl;

  deallog << std::endl;
}



int
main()
{
  initlog();

  run_test<1>();
  run_test<2>();
  run_test<3>();
}
