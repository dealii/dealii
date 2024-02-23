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
 * Test the function find_best_height_direction in
 * NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/non_matching/quadrature_generator.h>

#include "../tests.h"


using namespace NonMatching::internal::QuadratureGeneratorImplementation;


// Return a pair with both entries equal to value.
std::pair<double, double>
pair_with_equal_entries(const double value)
{
  return std::pair<double, double>(value, value);
}



/*
 * Test that find_best_height_direction returns an unset optional if the
 * incoming bounds correspond to negative/positive definite functions.
 */
void
test_ignores_definite_functions()
{
  const int dim = 2;
  deallog << "test_ignores_definite_functions" << std::endl;

  // Bounds corresponding to one negative and one positive definite function.
  std::vector<FunctionBounds<dim>> bounds(2);
  bounds[0].value = pair_with_equal_entries(-1);
  bounds[1].value = pair_with_equal_entries(1);

  const std::optional<HeightDirectionData> data =
    find_best_height_direction(bounds);

  if (!data)
    deallog << "OK" << std::endl;
}



/**
 * Create a vector containing two FunctionBounds, set them up so that
 * there is one height function direction that is the best. Test that this is
 * the direction returned from find_best_height_direction().
 */
void
test_find_best_height_direction()
{
  deallog << "test_find_best_height_direction" << std::endl;

  const int dim = 2;

  std::vector<FunctionBounds<dim>> bounds(2);
  // Set up so that the bounds correspond to indefinite functions.
  for (unsigned int i = 0; i < bounds.size(); ++i)
    {
      bounds[i].value.first  = -1;
      bounds[i].value.second = 1;
    }

  // Set up the bounds so that the componenetwise min (over function bounds)
  // of the gradient is [3, 5]. This makes 1 the best direction.
  bounds[0].gradient[0] = pair_with_equal_entries(3);
  bounds[1].gradient[0] = pair_with_equal_entries(4);
  bounds[0].gradient[1] = pair_with_equal_entries(6);
  bounds[1].gradient[1] = pair_with_equal_entries(5);

  const std::optional<HeightDirectionData> data =
    find_best_height_direction(bounds);

  deallog << "height direction = " << data->direction << std::endl;
  deallog << "min_abs_dfdx = " << data->min_abs_dfdx << std::endl;
}



int
main()
{
  initlog();
  test_ignores_definite_functions();
  deallog << std::endl;
  test_find_best_height_direction();
}
