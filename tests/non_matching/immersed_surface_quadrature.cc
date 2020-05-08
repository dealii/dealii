// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include "../tests.h"


// Test that an ImmersedSurfaceQuadrature can be constructed for each dimension
// and that quadrature points can be added to it.



template <int dim>
void
print_quadrature(const NonMatching::ImmersedSurfaceQuadrature<dim> &quadrature)
{
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    {
      deallog << quadrature.point(i) << ", " << quadrature.weight(i) << ", "
              << quadrature.normal_vector(i) << std::endl;
    }
}



// Check that get_normals() are callable and are of the same size as
// points and weights.
template <int dim>
void
check_get_normals(const NonMatching::ImmersedSurfaceQuadrature<dim> &quadrature)
{
  const std::vector<Point<dim>> &    points  = quadrature.get_points();
  const std::vector<Tensor<1, dim>> &normals = quadrature.get_normal_vectors();
  AssertThrow(points.size() == normals.size(), ExcInternalError())
}



template <int dim>
void
test_non_default_constructor()
{
  deallog << "Using constructor" << std::endl;
  std::vector<Point<dim>>     points(1);
  std::vector<double>         weights(1, 1);
  std::vector<Tensor<1, dim>> normals;
  normals.push_back(Point<dim>::unit_vector(dim - 1));
  NonMatching::ImmersedSurfaceQuadrature<dim> quadrature(points,
                                                         weights,
                                                         normals);

  print_quadrature(quadrature);
}



template <int dim>
void
test_push_back()
{
  deallog << "Using push_back" << std::endl;
  const Point<dim>     point;
  const double         weight = 1;
  const Tensor<1, dim> normal = Point<dim>::unit_vector(dim - 1);

  NonMatching::ImmersedSurfaceQuadrature<dim> quadrature;
  quadrature.push_back(point, weight, normal);

  print_quadrature(quadrature);
}



template <int dim>
void
construct_quadrature_and_print_points()
{
  test_push_back<dim>();
  test_non_default_constructor<dim>();
}



int
main()
{
  initlog();
  construct_quadrature_and_print_points<1>();
  construct_quadrature_and_print_points<2>();
  construct_quadrature_and_print_points<3>();
}
