//----------------------------------------------------------------------
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
//----------------------------------------------------------------------

#include "interpolate_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/fe/fe_raviart_thomas.h>

#include <fstream>

// FE_RaviartThomas<dim>::interpolate(...)

template <int dim>
void check1(const Function<dim>& f,
	    const unsigned int degree)
{
  FE_RaviartThomas<dim> fe(degree);
  deallog << fe.get_name() << ' ';
  deallog << fe.get_generalized_support_points().size() << ' ';
  
  std::vector<double> dofs(fe.dofs_per_cell);
  
  std::vector<std::vector<double> >
    values(dim, std::vector<double>(fe.get_generalized_support_points().size()));
  std::vector<Vector<double> >
    vectors(fe.get_generalized_support_points().size(),
	    Vector<double>(dim));
  f.vector_value_list(fe.get_generalized_support_points(), vectors);

  for (unsigned int c=0;c<values.size();++c)
    for (unsigned int k=0;k<values[c].size();++k)
      values[c][k] = vectors[k](c);
  
  fe.interpolate(dofs, values);
  deallog << " vector " << vector_difference(fe,dofs,f,0);
  
  fe.interpolate(dofs, vectors, 0);
  deallog << " Vector " << vector_difference(fe,dofs,f,0) << std::endl;
}

int main()
{
  std::ofstream logfile ("interpolate_rt/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-15);

//   Q1WedgeFunction<1,1,2> w1;
//   check1(w1,1,2);
//   check1(w1,2,2);
//   check1(w1,3,2);
  Q1WedgeFunction<2,1,2> w21;
  check1(w21,1);
  check1(w21,2);
  Q1WedgeFunction<2,2,2> w22;
  check1(w22,2);
  Q1WedgeFunction<2,3,2> w23;
  check1(w23,3);
//  Q1WedgeFunction<3,1,3> w3;
//  check1(w3,1);
//  check1(w3,2);
}
