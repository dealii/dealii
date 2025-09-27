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
 * Test the RootFinder class in internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/function_signed_distance.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/point.h>
#include <deal.II/base/tensor.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"


using namespace NonMatching::internal::QuadratureGeneratorImplementation;

// Use RootFinder to find the roots of the incoming functions over the interval
// [0, 1]. Print the roots to deallog.
void
find_and_print_roots(
  const std::vector<std::reference_wrapper<const Function<1>>> &functions)
{
  const BoundingBox<1> interval = create_unit_bounding_box<1>();

  std::vector<double> roots;
  RootFinder          root_finder;
  root_finder.find_roots(functions, interval, roots);

  for (unsigned int i = 0; i < roots.size(); ++i)
    {
      deallog << roots[i];
      if (i < roots.size() - 1)
        deallog << ", ";
    }
  deallog << std::endl;
}



// Test that the roots we get back from RootFinder are sorted and
// not duplicated.
//
// Call find_roots with 3 linear functions f_i(x) = x - x_i,
// where x_0 = 0.75, x_1 = 0.25, x_2 = 0.25
// and check that RootFinder gives back the vector {0.25, 0.75}.
void
test_roots_sorted_not_duplicated()
{
  deallog << "test_roots_sorted_not_duplicated" << std::endl;

  std::vector<Functions::SignedDistance::Plane<1>> linear_functions;

  const std::vector<double> roots = {.75, .25, .25};
  for (unsigned int i = 0; i < roots.size(); ++i)
    {
      Tensor<1, 1> normal;
      normal[0] = 1;
      const Point<1> point(roots.at(i));
      linear_functions.push_back(
        Functions::SignedDistance::Plane<1>(point, normal));
    }

  const std::vector<std::reference_wrapper<const Function<1>>> functions(
    linear_functions.begin(), linear_functions.end());

  find_and_print_roots(functions);
}



/*
 * The function:
 * f(x) = C(x - x_0)^2 + y_0
 */
class QuadraticFunction : public Function<1>
{
public:
  QuadraticFunction(const double C, const double x_0, const double y_0)
    : C(C)
    , x_0(x_0)
    , y_0(y_0)
  {}

  double
  value(const Point<1> &point, const unsigned int component = 0) const override
  {
    return C * std::pow(point[0] - x_0, 2) + y_0;
  };

  Tensor<1, 1>
  gradient(const Point<1>    &point,
           const unsigned int component = 0) const override
  {
    Tensor<1, 1> grad;
    grad[0] = 2 * C * (point[0] - x_0);

    return grad;
  };

  SymmetricTensor<2, 1>
  hessian(const Point<1>    &point,
          const unsigned int component = 0) const override
  {
    SymmetricTensor<2, 1> grad;
    grad[0][0] = 2 * C;

    return grad;
  };

private:
  const double C;
  const double x_0;
  const double y_0;
};



// Test that RootFinder can find both roots of the function
// f(x) = 4(x-0.5)^2 - 0.25
// which are x_0 = 0.25 and x_1 = 0.75.
void
test_find_both_roots()
{
  deallog << "test_find_both_roots" << std::endl;

  const QuadraticFunction function(4, 0.5, -0.25);

  std::vector<std::reference_wrapper<const Function<1>>> functions;
  functions.push_back(function);

  find_and_print_roots(functions);
}



int
main()
{
  initlog();
  test_roots_sorted_not_duplicated();
  deallog << std::endl;
  test_find_both_roots();
}
