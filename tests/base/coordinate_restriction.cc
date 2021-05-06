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
 * Test that CoordinateRestriction returns the derivatives of its underlying
 * function truncated correctly.
 */


#include <deal.II/base/function_restriction.h>

#include "../tests.h"

namespace
{
  using namespace dealii;

  /*
   * A simple test-function which is constant in space, but varies
   * in each component.
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
        for (unsigned int j = 0; j < dim; ++j)
          tensor[i][j] = 1 + i + j;

      return tensor;
    }
  };


  /*
   * Create a CoordinateRestriction that has the incoming direction of the
   * above TestFunction restricted. Write value, gradient & hessian
   * to deallog to check that they are truncated as expected.
   */
  template <int dim>
  void
  test_truncates_derivatives_correctly(const unsigned int restricted_direction)
  {
    deallog << "Restricted direction = " << restricted_direction << std::endl;

    const TestFunction<dim> function;
    const double            coordinate_value = 0;

    const Functions::CoordinateRestriction<dim - 1> restriction(
      function, restricted_direction, coordinate_value);

    const Point<dim - 1> point;
    const unsigned int   component = 0;

    deallog << "value = " << restriction.value(point, component) << std::endl;
    deallog << "gradient = " << restriction.gradient(point, component)
            << std::endl;
    deallog << "Hessian = " << restriction.hessian(point, component)
            << std::endl;

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
