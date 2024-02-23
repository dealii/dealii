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

/*
 * Test the QuadratureGenerator class, by setting up a few simple cuts over the
 * unit box and writing the generated quadrature rules to the output file.
 *
 * Each function beginning with "test_" sets up a level set function and then
 * calls the function create_and_print_quadratures() to generate the
 * quadratures.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/point.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/hp/q_collection.h>

#include <deal.II/lac/vector.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <vector>

#include "../tests.h"

#include "quadrature_printing.h"

using NonMatching::QuadratureGenerator;

/*
 * Create immersed quadrature rules over a unit box intersected by the
 * incoming level set function. Use a 1D-Gauss quadrature with n_1D_points
 * points as a base. Print the constructed quadrature rules to deallog.
 */
template <int dim>
void
create_and_print_quadratures(
  const Function<dim>                                     &level_set,
  const unsigned int                                       n_1D_points = 2,
  const typename QuadratureGenerator<dim>::AdditionalData &additional_data =
    typename QuadratureGenerator<dim>::AdditionalData())
{
  deallog << "dim=" << dim << std::endl;

  hp::QCollection<1> q_collection;
  q_collection.push_back(QGauss<1>(n_1D_points));

  QuadratureGenerator<dim> quadrature_generator(q_collection, additional_data);

  const BoundingBox<dim> box = create_unit_bounding_box<dim>();
  quadrature_generator.generate(level_set, box);

  deallog << "Inside quadrature" << std::endl;
  print_quadrature(quadrature_generator.get_inside_quadrature());
  deallog << "Outside quadrature" << std::endl;
  print_quadrature(quadrature_generator.get_outside_quadrature());
  deallog << "Surface quadrature" << std::endl;
  print_surface_quadrature(quadrature_generator.get_surface_quadrature());
}



/*
 * Construct level set with a zero contour as a plane cutting straight
 * through the unit box. Create and print the constructed quadratures. Do this
 * for all unit normals aligned with the coordinate directions. We expect that
 * the constructed quadrature has equal number of points in the inside/outside
 * region and that they are tensor products.
 */
template <int dim>
void
test_vertical_cuts_through_center()
{
  deallog << "test_vertical_cuts_through_center" << std::endl;

  Point<dim> center;
  for (int d = 0; d < dim; ++d)
    {
      center[d] = .5;
    }
  for (int direction = 0; direction < dim; ++direction)
    {
      deallog << "direction=" << direction << std::endl;
      const Tensor<1, dim> normal = Point<dim>::unit_vector(direction);
      const Functions::SignedDistance::Plane<dim> level_set(center, normal);
      create_and_print_quadratures(level_set);
    }
}



/*
 * Set up a constant negative/positive level set function. Check that
 * the constructed quadratures only the inside/outside have points and that
 * this quadrature is a tensor product.
 */
template <int dim>
void
test_constant_level_sets_both_signs()
{
  const Functions::ConstantFunction<dim> constant_positive(1);
  const Functions::ConstantFunction<dim> constant_negative(-1);
  deallog << std::endl;

  deallog << "constant_positive" << std::endl;
  create_and_print_quadratures(constant_positive);

  deallog << std::endl;

  deallog << "constant_negative" << std::endl;
  create_and_print_quadratures(constant_negative);
}



// Set up a level set function corresponding to a plane with normal (1,1) in 2D
// and (1,1,1) in 3D. This makes the inside region a simplex, with vertices
// (0, 0), (0, l), (l, 0), in 2D,
// (0, 0 ,0), (0, 0, l), (0, l, 0), (l, 0, 0), in 3D.
// where l is the edge length.
// This is a good test because we know that the inside weights should sum to the
// area/volume: $V = l^{dim}/dim!$, and that the surface weights should sum
// to $S = \sqrt(2) l$ in 2D and $S = \sqrt(3)/2 l^2$ in 3D.
template <int dim>
void
test_simplex_cut()
{
  deallog << "test_simplex_cut" << std::endl;

  const double edge_length = 1. / std::sqrt(2);

  Tensor<1, dim> normal;
  for (int i = 0; i < dim; ++i)
    normal[i] = 1;

  Point<dim> point_in_plane;
  point_in_plane[0] = edge_length;

  const Functions::SignedDistance::Plane<dim> level_set(point_in_plane, normal);

  create_and_print_quadratures(level_set);
}



// Set up a level set function with a zero contour being a plane in the
// direction (1,1) in 2D and (1,1,1) in 3D, such that it cuts the bottom corner
// of the reference cell with a cut of size epsilon. Test that the epsilon cut
// is ignored and we get a tensor product quadrature over the outside region.
template <int dim>
void
test_epsilon_cut_at_bottom_corner()
{
  deallog << "test_epsilon_cut_at_bottom_corner" << std::endl;
  const double   epsilon = 1e-15;
  Tensor<1, dim> normal;
  Point<dim>     center;
  for (int i = 0; i < dim; ++i)
    {
      normal[i] = 1;
      center[i] += epsilon;
    }
  const Functions::SignedDistance::Plane<dim> level_set(center, normal);

  create_and_print_quadratures(level_set);
}



/*
 * Set up a spherical level set with radius R centered in (0, R) in 2D and
 * (0, 0, R) in 3D. The result of this is that the zero contour of the level set
 * function cuts exactly through vertex 0 of the unit box
 *
 * When we choose to split the cell in 4. This test case is difficult for the
 * algorithm in 3D. We first get dim - 1 as the first height-direction. But
 * after restricting once we get L_a = 0 (where L_a is defined by |\partial_i
 * psi| > L_a) for both i = 1,2. Thus we can not choose a second height function
 * direction. The results is that the cell is split several times until the
 * maximum recursion is reached. When this happens the algorithm uses the
 * midpoint method as a fallback.
 */
template <int dim>
void
test_sphere_cutting_corner_exactly()
{
  deallog << "test_sphere_cutting_corner_exactly" << std::endl;
  const double radius = 4;
  Point<dim>   center;
  center[dim - 1] = radius;
  const Functions::SignedDistance::Sphere<dim> level_set(center, radius);

  typename QuadratureGenerator<dim>::AdditionalData data;
  data.split_in_half  = false;
  data.max_box_splits = 2;

  const unsigned int n_1D_points = 2;

  create_and_print_quadratures(level_set, n_1D_points, data);
}



// A "fake" function used in test_splitting(). This function is constant 1,
// except close to the unit box center, x_i = 0.5, where it has a very large
// Hessian.
template <int dim>
class ConstantOneButLargeHessianInCenter
  : public Functions::ConstantFunction<dim>
{
public:
  ConstantOneButLargeHessianInCenter()
    : Functions::ConstantFunction<dim>(1)
  {
    for (int d = 0; d < dim; ++d)
      unit_box_center[d] = .5;
  }

  SymmetricTensor<2, dim>
  hessian(const Point<dim> &point, const unsigned int) const override
  {
    SymmetricTensor<2, dim> hessian;

    const double max_distance = 1e-3;
    const double diagonal_value =
      point.distance(unit_box_center) < max_distance ? 1E3 : 0;

    for (int d = 0; d < dim; ++d)
      hessian[d][d] = diagonal_value;

    return hessian;
  }

private:
  Point<dim> unit_box_center;
};



// Test the box splitting. Call QuadratureGenerator with a function that is
// constant 1, but has a large Hessian close to the center of the box. This
// should make the algorithm split the box, since the function bounds will be
// large.
template <int dim>
void
test_splitting()
{
  deallog << "test_splitting" << std::endl;

  const ConstantOneButLargeHessianInCenter<dim> level_set;
  create_and_print_quadratures(level_set);
}



// Some of the test cases only make sense for a given dimension,
// so we list the cases for each dimension.
int
main()
{
  initlog();
  // 1D
  test_vertical_cuts_through_center<1>();
  deallog << std::endl;
  // 2D
  test_vertical_cuts_through_center<2>();
  deallog << std::endl;
  test_constant_level_sets_both_signs<2>();
  deallog << std::endl;
  test_simplex_cut<2>();
  deallog << std::endl;
  test_epsilon_cut_at_bottom_corner<2>();
  deallog << std::endl;
  test_sphere_cutting_corner_exactly<2>();
  deallog << std::endl;
  test_splitting<2>();
  deallog << std::endl;
  // 3D
  test_vertical_cuts_through_center<3>();
  deallog << std::endl;
  test_simplex_cut<3>();
  deallog << std::endl;
  test_epsilon_cut_at_bottom_corner<3>();
  deallog << std::endl;
  test_sphere_cutting_corner_exactly<3>();
}
