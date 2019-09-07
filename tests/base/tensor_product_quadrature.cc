// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2018 by the deal.II authors
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

#include <deal.II/base/quadrature.h>

#include "../tests.h"

// Write points and weights of the incoming quadrature to deallog.
template <int dim>
void
print_quadrature(const Quadrature<dim> &quadrature)
{
  deallog << "points" << std::endl;
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    deallog << quadrature.point(i) << ", ";
  deallog << std::endl;

  deallog << "weights" << std::endl;
  for (unsigned int i = 0; i < quadrature.size(); ++i)
    deallog << quadrature.weight(i) << ", ";
  deallog << std::endl;
}



// Write points weights and 1D base quadratures to deallog
template <int dim>
void
print_tensor_product_quadrature(const TensorProductQuadrature<dim> &quadrature)
{
  print_quadrature(quadrature);

  const std::array<Quadrature<1>, dim> &tensor_basis =
    quadrature.get_tensor_basis();
  for (unsigned int i = 0; i < dim; ++i)
    {
      deallog << "base quadrature " << i << std::endl;
      print_quadrature(tensor_basis[i]);
    }
}



template <int dim>
void
test_constructor_taking_1D_quadratures()
{
  deallog << "Constructor taking 1D Quadratures" << std::endl;

  const std::vector<Point<1>> points  = {Point<1>(.2), Point<1>(.8)};
  const std::vector<double>   weights = {.1, .9};
  const Quadrature<1>         quadrature_1D(points, weights);

  const TensorProductQuadrature<dim> tensor_product_quadrature(quadrature_1D);
  print_tensor_product_quadrature(tensor_product_quadrature);
}



template <int dim>
void
test_constructor_taking_1D_and_subquadrature()
{
  deallog << "Constructor taking 1D and SubQuadrature" << std::endl;

  const std::vector<Point<1>> points1D  = {Point<1>(.2), Point<1>(.8)};
  const std::vector<double>   weights1D = {.1, .9};
  const Quadrature<1>         quadrature_1D(points1D, weights1D);

  // Need to create a SubTensorProductQuadrautre to give the constructor. Create
  // a different 1-dimensional quadrature to create the
  // SubTensorProductQuadrautre from.
  const std::vector<Point<1>> points1D_for_sub  = {Point<1>(.4), Point<1>(.6)};
  const std::vector<double>   weights1D_for_sub = {.3, .7};
  const Quadrature<1>         quadrature_1D_for_sub(points1D_for_sub,
                                            weights1D_for_sub);

  const TensorProductQuadrature<dim - 1> subquadrature(quadrature_1D_for_sub);

  const TensorProductQuadrature<dim> tensor_product_quadrature(subquadrature,
                                                               quadrature_1D);
  print_tensor_product_quadrature(tensor_product_quadrature);
}



template <int dim>
void
test_constructor_taking_n_points()
{
  deallog << "Constructor taking n_points" << std::endl;

  const unsigned int                 n_points = 1;
  const TensorProductQuadrature<dim> tensor_product_quadrature(n_points);
  AssertThrow(n_points == tensor_product_quadrature.size(), ExcInternalError());
  deallog << "OK" << std::endl;
}



template <int dim>
void
test()
{
  deallog << "dim = " << dim << std::endl;

  test_constructor_taking_1D_quadratures<dim>();
  deallog << std::endl;

  test_constructor_taking_1D_and_subquadrature<dim>();
  deallog << std::endl;

  test_constructor_taking_n_points<dim>();
}



int
main()
{
  initlog();
  test<2>();
}
