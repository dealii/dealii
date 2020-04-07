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

// Test the following functions of the BoundingBox class
// cross_section
// center
// side_length
// child
// vertex
// bounds
// and the non-member but related function create_unit_bounding_box
// Each function is tested in the function named test_{member_function_name}.

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/geometry_info.h>

#include "../tests.h"

using namespace dealii;

// Print the bounds of the incoming box to deallog.
template <int dim>
void
print_box(const BoundingBox<dim> &box)
{
  for (unsigned int d = 0; d < dim; ++d)
    {
      if (d > 0)
        deallog << "x";
      deallog << "[" << box.lower_bound(d) << ", " << box.upper_bound(d) << "]";
    }
  deallog << std::endl;
}



// Construct the cross section orthogonal to all axes and print them.
template <int dim>
void
test_cross_section()
{
  deallog << "test_cross_section" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    {
      lower_upper_corners.first[d]  = -(d + 1);
      lower_upper_corners.second[d] = d + 1;
    }

  const BoundingBox<dim> box(lower_upper_corners);

  for (unsigned int d = 0; d < dim; ++d)
    {
      deallog << "orthogonal to " << d << std::endl;
      const BoundingBox<dim - 1> cross_section = box.cross_section(d);
      print_box(cross_section);
    }
}



// Compute and print the center of the box.
template <int dim>
void
test_center()
{
  deallog << "test_center" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    lower_upper_corners.second[d] = 1;

  const BoundingBox<dim> box(lower_upper_corners);

  deallog << "center " << box.center() << std::endl;
}



// Print all side lengths of a box.
template <int dim>
void
test_side_length()
{
  deallog << "test_side_length" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    lower_upper_corners.second[d] = d + 1;

  const BoundingBox<dim> box(lower_upper_corners);

  for (unsigned int d = 0; d < dim; ++d)
    deallog << box.side_length(d) << " ";
  deallog << std::endl;
}



// Construct all possible children of a box and print them.
template <int dim>
void
test_child()
{
  deallog << "test_child" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    {
      lower_upper_corners.first[d]  = -(d + 1);
      lower_upper_corners.second[d] = (d + 1);
    }

  const BoundingBox<dim> box(lower_upper_corners);

  for (unsigned int i = 0; i < GeometryInfo<dim>::max_children_per_cell; ++i)
    {
      deallog << "child " << i << std::endl;
      const BoundingBox<dim> child = box.child(i);
      print_box(child);
    }
}



// Print all the vertices of a box.
template <int dim>
void
test_vertex()
{
  deallog << "test_vertex" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    {
      lower_upper_corners.first[d]  = -(d + 1);
      lower_upper_corners.second[d] = d + 1;
    }

  const BoundingBox<dim> box(lower_upper_corners);

  for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      deallog << "vertex " << i << " = " << box.vertex(i) << std::endl;
    }
}



// Print the results of the bounds function for all coordinate directions.
template <int dim>
void
test_bounds()
{
  deallog << "test_bounds" << std::endl;

  std::pair<Point<dim>, Point<dim>> lower_upper_corners;
  for (int d = 0; d < dim; ++d)
    {
      lower_upper_corners.first[d]  = -(d + 1);
      lower_upper_corners.second[d] = d + 1;
    }

  const BoundingBox<dim> box(lower_upper_corners);

  for (unsigned int i = 0; i < dim; ++i)
    {
      const BoundingBox<1> bounds = box.bounds(i);
      deallog << "Bounds in direction " << i << std::endl;
      print_box(bounds);
    }
}



// Test that we can call the create_unit_box function.
template <int dim>
void
test_create_unit_bounding_box()
{
  deallog << "test_create_unit_bounding_box" << std::endl;
  BoundingBox<dim> box = create_unit_bounding_box<dim>();
  print_box(box);
}



template <int dim>
void
run_test()
{
  deallog << "dim = " << dim << std::endl;

  test_cross_section<dim>();
  deallog << std::endl;

  test_center<dim>();
  deallog << std::endl;

  test_side_length<dim>();
  deallog << std::endl;

  test_child<dim>();
  deallog << std::endl;

  test_vertex<dim>();
  deallog << std::endl;

  test_bounds<dim>();
  deallog << std::endl;

  test_create_unit_bounding_box<dim>();
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
