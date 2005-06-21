//----------------------------  fe_tools_cpfqpm_03.cc  ---------------------------
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
//----------------------------  fe_tools_cpfqpm_03.cc  ---------------------------

#include "../tests.h"
#include "fe_tools_common.cc"
#include <base/quadrature_lib.h>

// check
//   FETools::compute_projection_from_quadrature_points_matrix
// we put this into the fe_tools_common framework for simplicity, but
// in fact we ignore the second FE it passes to the check_this() function
// and we can only test as well for scalar elements, since this is all
// the function presently supports.
//
// this test simply computes the matrix and outputs some of its
// characteristics. in contrast to the _01 test, we choose a FE that has less
// DoFs than there are quadrature points -- not an uncommon case, and one that
// needs to work


std::string output_file_name = "fe_tools_cpfqpm_03.output";


template <int dim>
void
check_this (const FiniteElement<dim> &fe,
            const FiniteElement<dim> &/*fe2*/)
{
                                   // only check if both elements have
                                   // support points. otherwise,
                                   // interpolation doesn't really
                                   // work
  if (fe.n_components() != 1)
    return;

                                   // ignore this check if this fe has already
                                   // been treated
  static std::set<std::string> already_checked;
  if (already_checked.find(fe.get_name()) != already_checked.end())
    return;
  already_checked.insert (fe.get_name());
  

                                   // test with different quadrature formulas
  QGauss<dim> q_lhs(fe.degree+1);
  QGauss<dim> q_rhs(fe.degree+3);
      
  FullMatrix<double> X (fe.dofs_per_cell,
                        q_rhs.n_quadrature_points);

  FETools::compute_projection_from_quadrature_points_matrix (fe,
                                                             q_lhs, q_rhs,
                                                             X);
      
  output_matrix (X);
}

