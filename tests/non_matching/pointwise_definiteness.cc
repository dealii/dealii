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
 * Test the function pointwise_definiteness in
 * NonMatching::internal::QuadratureGeneratorImplementation.
 */

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/non_matching/quadrature_generator.h>

#include <vector>

#include "../tests.h"

using namespace NonMatching::internal::QuadratureGeneratorImplementation;


/**
 * Call pointwise_definiteness with two positive Functions,
 * check that it returns Definiteness::positive.
 */
template <int dim>
void
test_with_positive_functions()
{
  std::vector<Functions::ConstantFunction<dim>> functions;
  functions.push_back(Functions::ConstantFunction<dim>(1));
  functions.push_back(Functions::ConstantFunction<dim>(1));

  const std::vector<std::reference_wrapper<const Function<dim>>> function_refs(
    functions.begin(), functions.end());

  const Definiteness definiteness =
    pointwise_definiteness(function_refs, Point<dim>());

  AssertThrow(definiteness == Definiteness::positive, ExcInternalError());
}



/**
 * Call pointwise_definiteness with two negative Functions,
 * check that it returns Definiteness::negative.
 */
template <int dim>
void
test_with_negative_functions()
{
  std::vector<Functions::ConstantFunction<dim>> functions;
  functions.push_back(Functions::ConstantFunction<dim>(-1));
  functions.push_back(Functions::ConstantFunction<dim>(-1));

  const std::vector<std::reference_wrapper<const Function<dim>>> function_refs(
    functions.begin(), functions.end());

  const Definiteness definiteness =
    pointwise_definiteness(function_refs, Point<dim>());

  AssertThrow(definiteness == Definiteness::negative, ExcInternalError());
}



/**
 * Call pointwise_definiteness with with one positive and one negative Function,
 * check that it returns Definiteness::indefinite.
 */
template <int dim>
void
test_with_functions_of_different_sign()
{
  std::vector<Functions::ConstantFunction<dim>> functions;
  functions.push_back(Functions::ConstantFunction<dim>(-1));
  functions.push_back(Functions::ConstantFunction<dim>(1));

  const std::vector<std::reference_wrapper<const Function<dim>>> function_refs(
    functions.begin(), functions.end());

  const Definiteness definiteness =
    pointwise_definiteness(function_refs, Point<dim>());

  AssertThrow(definiteness == Definiteness::indefinite, ExcInternalError());
}



/**
 * Call pointwise_definiteness with a single Function which is zero,
 * check that it returns Definiteness::indefinite.
 *
 * This is a special case in the implementation.
 */
template <int dim>
void
test_first_function_zero()
{
  Functions::ZeroFunction<dim> zero_function;

  std::vector<std::reference_wrapper<const Function<dim>>> function_refs;
  function_refs.push_back(zero_function);

  const Definiteness definiteness =
    pointwise_definiteness(function_refs, Point<dim>());

  AssertThrow(definiteness == Definiteness::indefinite, ExcInternalError());
}



template <int dim>
void
run_test()
{
  test_with_positive_functions<dim>();
  test_with_negative_functions<dim>();
  test_with_functions_of_different_sign<dim>();
  test_first_function_zero<dim>();
}



int
main()
{
  initlog();
  run_test<1>();
  run_test<2>();
  run_test<3>();
  deallog << "OK" << std::endl;
}
