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


/*
 * Test that PointRestriction returns the derivatives of its underlying function
 * truncated correctly.
 */


#include <deal.II/base/function_restriction.h>

#include "../tests.h"

namespace
{
  using namespace dealii;

  /*
   * A simple test-function which is constant in space, but varies
   * in each gradient component, and on the diagonal of the Hessian.
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
      Tensor<1, dim> tensor;
      for (unsigned int i = 0; i < dim; ++i)
        tensor[i] = i + 1;
      return tensor;
    }

    SymmetricTensor<2, dim>
    hessian(const Point<dim> &, const unsigned int) const override
    {
      SymmetricTensor<2, dim> tensor;
      for (unsigned int i = 0; i < dim; ++i)
        tensor[i][i] = i + 1;

      return tensor;
    }
  };


  /*
   * Create a PointRestriction that has the incoming direction open and has
   * the above TestFunction as underlying function. Check that the values that
   *
   * PointRestriction::values
   * PointRestriction::gradient
   * PointRestriction::hessian
   *
   * return are
   *
   * value
   * gradient[open_direction]
   * hessian[open_direction][open_direction]
   *
   * where value, gradient and hessian are the values from TestFunction.
   */
  template <int dim>
  void
  test_truncates_derivatives_correctly(const unsigned int open_direction)
  {
    deallog << "open direction = " << open_direction << std::endl;
    const Point<dim - 1>    locked_point;
    const TestFunction<dim> function;

    const Functions::PointRestriction<dim - 1> restriction(function,
                                                           open_direction,
                                                           locked_point);

    const Point<1>     point;
    const unsigned int component = 0;

    deallog << restriction.value(point, component) << std::endl;
    deallog << restriction.gradient(point, component) << std::endl;
    deallog << restriction.hessian(point, component) << std::endl;

    deallog << std::endl;
  }



  template <int dim>
  void
  test_all_coordinate_directions()
  {
    deallog << "dim = " << dim << std::endl;

    for (unsigned int direction = 0; direction < dim; direction++)
      test_truncates_derivatives_correctly<dim>(direction);
  }

} // namespace

int
main()
{
  initlog();
  test_all_coordinate_directions<2>();
  test_all_coordinate_directions<3>();
}
