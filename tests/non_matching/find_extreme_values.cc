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

/*
 * Test the function find_extreme_values in
 * NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"


using namespace dealii;
using namespace NonMatching::internal::QuadratureGeneratorImplementation;


/**
 * Send in a vector with two different function bounds to
 * find_extreme_values(). Check that what we get back are the actual extreme
 * values.
 */
template <int dim>
void
test_extreme_values_are_found()
{
  deallog << "test_extreme_values_are_found" << std::endl;

  std::vector<FunctionBounds<dim>> bounds(2);
  bounds[0].value.first  = 1;
  bounds[0].value.second = 2;
  bounds[1].value.first  = -1;
  bounds[1].value.second = 3;

  const std::pair<double, double> extremes = find_extreme_values(bounds);

  deallog << "min = " << extremes.first << std::endl;
  deallog << "max = " << extremes.second << std::endl;
}



/**
 * Since the implementation of find_extreme_values() treats the 0th entry
 * differently, we check that we get the same entry back if we send in a
 * vector with only one entry.
 */
template <int dim>
void
test_extreme_values_initialized_to_first()
{
  deallog << "test_extreme_values_initialized_to_first" << std::endl;

  std::vector<FunctionBounds<dim>> bounds(1);
  bounds[0].value.first  = 1;
  bounds[0].value.second = 2;

  const std::pair<double, double> extremes = find_extreme_values(bounds);

  deallog << "min = " << extremes.first << std::endl;
  deallog << "max = " << extremes.second << std::endl;
}



int
main()
{
  initlog();
  test_extreme_values_initialized_to_first<1>();
  deallog << std::endl;
  test_extreme_values_are_found<1>();
}
