//----------------------------  fe_tools_12.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_12.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.cc"

// check
//   FETools::get_projection_matrix


std::string output_file_name = "fe_tools_12.output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
  if (fe1.n_components() != 1)
    return;
  if (fe1.n_components() != fe2.n_components())
    return;
  
  FullMatrix<double> X (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FETools::get_projection_matrix (fe1, fe2, X);

  output_matrix (X);
}

