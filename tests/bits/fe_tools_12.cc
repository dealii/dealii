//----------------------------  fe_tools_12.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2005, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_12.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.h"

// check
//   FETools::get_projection_matrix


std::string output_file_name = "fe_tools_12/output";



template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
				   // use a higher output accuracy for this
				   // test. the reason is that many of the
				   // constraints are negative powers of 2,
				   // which have exact representations with 3
				   // or 4 digits of accuracy, but not with
				   // the usual 2 digits (for example, 0.375,
				   // which sometimes rounds to 0.38 and
				   // sometimes to 0.37, depending on how
				   // intermediate errors have accumulated)
  deallog << std::setprecision (3);

  if (fe1.n_components() != 1)
    return;
  if (fe1.n_components() != fe2.n_components())
    return;
  
  FullMatrix<double> X (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FETools::get_projection_matrix (fe1, fe2, X);

  output_matrix (X);
}

