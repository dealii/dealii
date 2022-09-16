// ---------------------------------------------------------------------
//
// Copyright (C) 2022 by the deal.II authors
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
