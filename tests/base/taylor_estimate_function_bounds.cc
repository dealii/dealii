// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2020 by the deal.II authors
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

#include <deal.II/base/function_tools.h>

#include "../tests.h"

/*
 * Test taylor_estimate_function_bounds, by calling it with a function and a box
 * such that the function bounds are simple to compute.
 */

namespace
{
  using namespace dealii;

  /**
   * Returns a bounding box with different side lengths, 2*h_i, in all
   * directions. Here, h_i = 1/(1 + i), i.e. the box is
   * [-1, 1] x [-1/2, 1/2] x [-1/3, 1/3] in 3D.
   */
  template <int dim>
  BoundingBox<dim>
  create_box()
  {
    std::pair<Point<dim>, Point<dim>> lower_upper_corner;
    for (unsigned int i = 0; i < dim; i++)
      {
        lower_upper_corner.first[i]  = -1.0 / (i + 1);
        lower_upper_corner.second[i] = 1.0 / (i + 1);
      }

    return BoundingBox<dim>(lower_upper_corner);
  }



  /*
   * The Taylor expansion multiplies the gradient entries with h_i,
   * and the Hessian entries with h_i * h_j, where h_i is the half the
   * side length of the of the box in the ith coordinate direction. Choose the
   * gradient and Hessian as
   *
   * grad[i] = 1 / h_i,
   * hessian[i][j] = 1 / (h_i * h_j),
   *
   * to get function bounds that are simple to compute.
   */
  template <int dim>
  class TestFunction : public Function<dim>
  {
  public:
    double
    value(const Point<dim> &, const unsigned int) const override
    {
      return 1;
    }

    Tensor<1, dim>
    gradient(const Point<dim> &, const unsigned int) const override
    {
      Tensor<1, dim> grad;
      for (unsigned int i = 0; i < dim; i++)
        grad[i] = i + 1;

      return grad;
    }

    SymmetricTensor<2, dim>
    hessian(const Point<dim> &, const unsigned int) const
    {
      SymmetricTensor<2, dim> hess;
      for (unsigned int i = 0; i < dim; i++)
        for (unsigned int j = 0; j < dim; j++)
          hess[i][j] = (i + 1) * (j + 1);

      return hess;
    }
  };



  void
  print_bounds(const std::pair<double, double> &bounds)
  {
    deallog << "[" << bounds.first << ", " << bounds.second << "]" << std::endl;
  }



  /**
   * Call taylor_estimate_function_bounds with the above TestFunction and
   * the box returned by create_box. Print the computed bounds to deallog.
   */
  template <int dim>
  void
  run_test()
  {
    deallog << dim << "D" << std::endl;

    const TestFunction<dim> function;
    const BoundingBox<dim>  box = create_box<dim>();

    std::pair<double, double>                  value_bounds;
    std::array<std::pair<double, double>, dim> gradient_bounds;
    FunctionTools::taylor_estimate_function_bounds<dim>(function,
                                                        box,
                                                        value_bounds,
                                                        gradient_bounds);

    deallog << "value: ";
    print_bounds(value_bounds);

    for (unsigned int i = 0; i < dim; i++)
      {
        deallog << "gradient[" << i << "]: ";
        print_bounds(gradient_bounds[i]);
      }
    deallog << std::endl;
  }
} // namespace



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
}
