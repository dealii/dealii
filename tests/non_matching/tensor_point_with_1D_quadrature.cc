// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2021 by the deal.II authors
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

/**
 * Test the function tensor_point_with_1D_quadrature()
 * in NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"

#include "quadrature_printing.h"

using namespace dealii;
using namespace NonMatching::internal::QuadratureGeneratorImplementation;


/*
 * Set up a (dim-1)-dimensional point and a 1D-quadrature. Call
 * tensor_point_with_1D_quadrature for each possible coordinate direction
 * and print the resulting dim-dimensional quadrature.
 */
template <int dim>
void
create_and_output_quadrature_for_each_direction()
{
  deallog << "dim=" << dim << std::endl;

  const unsigned int     n_quadrature_points = 2;
  const QGaussLobatto<1> quadrature1D(n_quadrature_points);
  // Choose the points coordinates to something
  // easily distinguished.
  Point<dim - 1> point;
  for (int i = 0; i < dim - 1; ++i)
    {
      point(i) = 10 * (i + 1);
    }
  // Both points in the 1D-quadrature have weight 1/2 so
  // this should also be the weight of the points in the final
  // quadrature.
  const double weight = 5;
  const double start = -1, end = 1;

  for (int direction = 0; direction < dim; ++direction)
    {
      deallog << "direction=" << direction << std::endl;
      ExtendableQuadrature<dim> result;
      tensor_point_with_1D_quadrature(
        point, weight, quadrature1D, start, end, direction, result);
      print_quadrature(result);
      deallog << std::endl;
    }
}



int
main()
{
  initlog();
  create_and_output_quadrature_for_each_direction<1>();
  create_and_output_quadrature_for_each_direction<2>();
  create_and_output_quadrature_for_each_direction<3>();
}
