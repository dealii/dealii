// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


#include "../tests.h"
#include "fe_tools_common.h"

// check
//   FETools::get_projection_matrix


std::string output_file_name = "output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  // use a higher output accuracy for this test. the reason is that many of
  // the constraints are negative powers of 2, which have exact
  // representations with 3 or 4 digits of accuracy, but not with the usual
  // 2 digits (for example, 0.375, which sometimes rounds to 0.38 and
  // sometimes to 0.37, depending on how intermediate errors have
  // accumulated)
  deallog << std::setprecision (8);

  if (fe1.n_components() != 1)
    return;
  if (fe1.n_components() != fe2.n_components())
    return;

  FullMatrix<double> X (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FETools::get_projection_matrix (fe1, fe2, X);

  output_matrix (X);
}

