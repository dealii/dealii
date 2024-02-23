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
 * Test the function find_extreme_values in
 * NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"


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
