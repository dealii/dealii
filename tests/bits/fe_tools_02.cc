//----------------------------  fe_tools_02.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools_02.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.cc"
#include <lac/sparsity_pattern.h>

// check
//   FETools::get_interpolation_matrix


std::string output_file_name = "fe_tools_02.output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe1,
            const FiniteElement<dim> &fe2)
{
                                   // only check if both elements have
                                   // support points. otherwise,
                                   // interpolation doesn't really
                                   // work.
  if ((fe1.get_unit_support_points().size() == 0) ||
      (fe2.get_unit_support_points().size() == 0))
    return;
                                   //  likewise for non-primitive elements
  if (!fe1.is_primitive() || !fe2.is_primitive())
    return;
  
  FullMatrix<double> m (fe2.dofs_per_cell,
                        fe1.dofs_per_cell);
  FETools::get_interpolation_matrix (fe1, fe2, m);

  output_matrix (m);
}

