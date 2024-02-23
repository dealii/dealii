// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include "../tests.h"

#include "fe_tools_common.h"

// check
//   FETools::get_projection_matrix



template <int dim>
void
check_this(const FiniteElement<dim> &fe1, const FiniteElement<dim> &fe2)
{
  // use a higher output accuracy for this test. the reason is that many of
  // the constraints are negative powers of 2, which have exact
  // representations with 3 or 4 digits of accuracy, but not with the usual
  // 2 digits (for example, 0.375, which sometimes rounds to 0.38 and
  // sometimes to 0.37, depending on how intermediate errors have
  // accumulated)
  deallog << std::setprecision(8);

  if (fe1.n_components() != 1)
    return;
  if (fe1.n_components() != fe2.n_components())
    return;

  FullMatrix<double> X(fe2.dofs_per_cell, fe1.dofs_per_cell);
  FETools::get_projection_matrix(fe1, fe2, X);

  output_matrix(X);
}
