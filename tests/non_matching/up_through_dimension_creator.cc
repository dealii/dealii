// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/**
 * Test the class UpThroughDimensionCreator
 * in NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"

#include "quadrature_printing.h"

using namespace NonMatching::internal::QuadratureGeneratorImplementation;


// Print each quadrature in the incoming QPartitioning.
template <int dim>
void
print(const QPartitioning<dim> &q_partitioning)
{
  deallog << "Negative" << std::endl;
  print_quadrature(q_partitioning.negative);
  deallog << "Positive" << std::endl;
  print_quadrature(q_partitioning.positive);
  deallog << "Indefinite" << std::endl;
  print_quadrature(q_partitioning.indefinite);
  deallog << "Surface" << std::endl;
  print_surface_quadrature(q_partitioning.surface);
}



// Let the height function direction be dim - 1 and let the lower dimensional
// quadrature contain a single point. Call UpThroughDimensionCreator with the
// incoming level set function over the unit box to generate a QPartitioning.
// Print the quadratures in the partitioning to make sure it is correct.
template <int dim>
void
create_and_print_partitioning(const Function<dim> &level_set)
{
  const hp::QCollection<1>                    q_collection1D(QGauss<1>(2));
  const NonMatching::AdditionalQGeneratorData additional_data;

  UpThroughDimensionCreator<dim, dim> up_through_dimension_creator(
    q_collection1D, additional_data);

  std::vector<std::reference_wrapper<const Function<dim>>> level_sets;
  level_sets.push_back(level_set);
  const BoundingBox<dim>   box = create_unit_bounding_box<dim>();
  const QMidpoint<dim - 1> low_dim_quadrature;
  const unsigned int       height_function_direction = dim - 1;

  QPartitioning<dim> q_partitioning;
  up_through_dimension_creator.generate(level_sets,
                                        box,
                                        low_dim_quadrature,
                                        height_function_direction,
                                        q_partitioning);
  print(q_partitioning);
}



// Set up a level set function with the zero contour along, x_{dim-1} = 0.5
// Call create_and_print_partitioning to test that the points are added as
// expected:
// "negative" points should have x_{dim-1} \in  (0, 0.5)
// "positive" points should have x_{dim-1} \in  (0.5, 1)
// surface points should have x_{dim-1} = 0.5
template <int dim>
void
test_cut_through_center()
{
  deallog << "test_cut_through_center" << std::endl;
  deallog << std::endl;

  Point<dim>     point_through_plane = .5 * Point<dim>::unit_vector(dim - 1);
  Tensor<1, dim> plane_normal        = Point<dim>::unit_vector(dim - 1);
  const Functions::SignedDistance::Plane<dim> level_set(point_through_plane,
                                                        plane_normal);

  create_and_print_partitioning(level_set);
}



// Fabricate the case when we have missed roots when creating the quadrature in
// the lower dimensions. See the comment in the implementation of
// UpThroughDimensionCreator::create_surface_point(..).
//
// In this test, the zero contour goes outside the cell but close to the "bottom
// face" at x_{dim-1} = 0. Check that the surface quadrature points gets placed
// on x_{dim-1} = 0. This is the "least bad" option we have.
template <int dim>
void
test_missed_roots_on_bottom_face()
{
  deallog << "test_missed_roots_on_bottom_face" << std::endl;
  deallog << std::endl;

  const Tensor<1, dim> plane_normal = Point<dim>::unit_vector(dim - 1);
  Point<dim>           point_in_plane;
  point_in_plane[dim - 1] = -.1;
  const Functions::SignedDistance::Plane<dim> level_set(point_in_plane,
                                                        plane_normal);

  create_and_print_partitioning(level_set);
}



// Same test as above, but the zero contour just outside the "top face" at
// x_{dim-1} = 1. Check that the surface quadrature points are placed on
// x_{dim-1} = 1.
template <int dim>
void
test_missed_roots_on_top_face()
{
  deallog << "test_missed_roots_on_top_face" << std::endl;
  deallog << std::endl;

  const Tensor<1, dim> plane_normal = Point<dim>::unit_vector(dim - 1);
  Point<dim>           point_in_plane;
  point_in_plane[dim - 1] = 1.1;
  const Functions::SignedDistance::Plane<dim> level_set(point_in_plane,
                                                        plane_normal);

  create_and_print_partitioning(level_set);
}



int
main()
{
  initlog();

  const int dim = 2;

  test_cut_through_center<dim>();
  deallog << std::endl;

  test_missed_roots_on_bottom_face<dim>();
  deallog << std::endl;

  test_missed_roots_on_top_face<dim>();
  deallog << std::endl;
}
