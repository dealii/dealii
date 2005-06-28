//----------------------------  fe_tools_cpfqpm_04.cc  ---------------------------
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
//----------------------------  fe_tools_cpfqpm_04.cc  ---------------------------

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
// the test makes sure that if the FE has support points and we use
// these support points as quadrature points, that the resulting
// matrix is the unit matrix


std::string output_file_name = "fe_tools_cpfqpm_04.output";


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

				   // only test elements with support
				   // points
  if (fe.has_support_points() == false)
    return;

                                   // test with different quadrature formulas
  Quadrature<dim> q_rhs(fe.get_unit_support_points(),
			std::vector<double> (fe.dofs_per_cell,
					     1./fe.dofs_per_cell));
      
  FullMatrix<double> X (fe.dofs_per_cell,
                        q_rhs.n_quadrature_points);

  Assert (X.m() == X.n(), ExcInternalError());
  
  FETools::compute_projection_from_quadrature_points_matrix (fe,
                                                             q_rhs, q_rhs,
                                                             X);

  for (unsigned int i=0; i<X.m(); ++i)
    X(i,i) -= 1;

  Assert (X.frobenius_norm() < 1e-12, ExcInternalError());
}

