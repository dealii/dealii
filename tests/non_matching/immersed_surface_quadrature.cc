// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/non_matching/immersed_surface_quadrature.h>

#include "../tests.h"


// Test that an ImmersedSurfaceQuadrature can be constructed for each dimension
// and that quadrature points can be added to it.



template <int dim, int spacedim>
void
print_quadrature(
  const NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> &quadrature)
{
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    {
      if (dim > 0) // Can't print 0-dim points.
        deallog << quadrature.point(i);
      deallog << ", " << quadrature.weight(i) << ", "
              << quadrature.normal_vector(i) << std::endl;
    }
}



// Check that get_normals() are callable and are of the same size as
// points and weights.
template <int dim, int spacedim>
void
check_get_normals(
  const NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> &quadrature)
{
  const std::vector<Point<dim>>     &points  = quadrature.get_points();
  const std::vector<Tensor<1, dim>> &normals = quadrature.get_normal_vectors();
  AssertThrow(points.size() == normals.size(), ExcInternalError());
}



template <int dim, int spacedim>
void
test_non_default_constructor()
{
  deallog << "Using constructor" << std::endl;
  std::vector<Point<dim>>          points(1);
  std::vector<double>              weights(1, 1);
  std::vector<Tensor<1, spacedim>> normals;
  normals.push_back(Point<spacedim>::unit_vector(spacedim - 1));
  NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> quadrature(points,
                                                                   weights,
                                                                   normals);

  print_quadrature(quadrature);
}



template <int dim, int spacedim>
void
test_push_back()
{
  deallog << "Using push_back" << std::endl;
  const Point<dim>          point;
  const double              weight = 1;
  const Tensor<1, spacedim> normal = Point<spacedim>::unit_vector(spacedim - 1);

  NonMatching::ImmersedSurfaceQuadrature<dim, spacedim> quadrature;
  quadrature.push_back(point, weight, normal);

  print_quadrature(quadrature);
}



template <int dim, int spacedim>
void
construct_quadrature_and_print_points()
{
  test_push_back<dim, spacedim>();
  test_non_default_constructor<dim, spacedim>();
}



int
main()
{
  initlog();
  construct_quadrature_and_print_points<1, 1>();
  construct_quadrature_and_print_points<2, 2>();
  construct_quadrature_and_print_points<3, 3>();

  // Face quadrature
  construct_quadrature_and_print_points<0, 1>();
  construct_quadrature_and_print_points<1, 2>();
  construct_quadrature_and_print_points<2, 3>();
}
