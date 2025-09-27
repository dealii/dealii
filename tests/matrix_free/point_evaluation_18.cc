// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// Check that all template specializations of FEPointEvaluation are
// compiling.

#include <deal.II/matrix_free/fe_point_evaluation.h>

#include "../tests.h"

template <int n_components, int dim, int spacedim, typename Number>
void
test()
{
  return; // nothing to do, since we are only interested if the code
          // compiles

  std::unique_ptr<Mapping<dim, spacedim>>       mapping;
  std::unique_ptr<FiniteElement<dim, spacedim>> fe;

  FEPointEvaluation<n_components, dim, spacedim, Number> fpe(
    *mapping, *fe, UpdateFlags::update_default);

  Triangulation<dim, spacedim> tria;

  fpe.reinit(tria.begin(), ArrayView<const Point<dim>>());

  fpe.evaluate(ArrayView<const Number>(), EvaluationFlags::values);

  fpe.integrate(ArrayView<Number>(), EvaluationFlags::values);
}

int
main()
{
  initlog();

  test<1, 1, 1, double>();
  test<2, 1, 1, double>();

  test<1, 2, 2, double>();
  test<2, 2, 2, double>();
  test<3, 2, 2, double>();

  test<1, 3, 3, double>();
  test<2, 3, 3, double>();
  test<3, 3, 3, double>();
  test<4, 3, 3, double>();

  test<1, 1, 1, float>();
  test<2, 1, 1, float>();

  test<1, 2, 2, float>();
  test<2, 2, 2, float>();
  test<3, 2, 2, float>();

  test<1, 3, 3, float>();
  test<2, 3, 3, float>();
  test<3, 3, 3, float>();
  test<4, 3, 3, float>();

  deallog << "OK!" << std::endl;
}
