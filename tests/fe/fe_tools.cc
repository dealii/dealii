//--------------------------------------------------------------------------
//    fe_tools.cc,v 1.1 2003/11/28 15:03:26 guido Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------------

// Test the cell matrices generated in FETools.

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <base/logstream.h>

#include <fe/fe_q.h>
#include <fe/fe_dgq.h>
#include <fe/fe_dgp.h>
#include <fe/fe_raviart_thomas.h>
#include <fe/fe_system.h>

#include <fe/fe_tools.h>


template<int dim>
void test_projection (const FiniteElement<dim>& fe1,
		      const FiniteElement<dim>& fe2,
		      std::ostream& out)
{
  out << fe1.get_name() << " -> "
      << fe2.get_name() << std::endl;
  
  const unsigned int n1 = fe1.dofs_per_cell;
  const unsigned int n2 = fe2.dofs_per_cell;

  FullMatrix<double> P(n2,n1);

  FETools::get_projection_matrix(fe1, fe2, P);
  P.print_formatted(out, 2, false, 5);
}

  
template<int dim>
void test_projection (std::ostream& out)
{
  FE_DGQ<dim> q0(0);
  FE_DGQ<dim> q1(1);
  FE_DGQ<dim> q2(2);
  FE_DGQ<dim> q3(3);
  FE_DGQ<dim> q4(4);

  FE_DGP<dim> p0(0);
  FE_DGP<dim> p1(1);
  FE_DGP<dim> p2(2);
  FE_DGP<dim> p3(3);
  FE_DGP<dim> p4(4);

  test_projection(p1,p0, out);
  test_projection(p0,p1, out);
  test_projection(p2,p1, out);
  test_projection(p1,p2, out);
  test_projection(p2,p0, out);
  test_projection(p0,p2, out);
  test_projection(p3,p2, out);
  test_projection(p2,p3, out);
  test_projection(p4,p3, out);
  test_projection(p3,p4, out);

  test_projection(q1,q0, out);
  test_projection(q2,q0, out);
  test_projection(q3,q0, out);
  test_projection(q4,q0, out);
  test_projection(q2,q1, out);
  test_projection(q1,q2, out);
  test_projection(q3,q2, out);
  test_projection(q4,q3, out);
}


int main()
{
  std::ofstream logfile("fe_tools.output");
  
  test_projection<1>(logfile);
  test_projection<2>(logfile);
  test_projection<3>(logfile);
}
