//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

#include <iostream>
#include <fstream>
#include <cmath>

#include <base/logstream.h>
#include <base/polynomial.h>
#include <base/quadrature_lib.h>


double scalar_product (const Polynomial<double>& p1,
		       const Polynomial<double>& p2)
{
  unsigned int degree = (p1.degree() + p2.degree())/2 + 1;
  QGauss<1> gauss(degree);

  double sum = 0.;
  for (unsigned int i=0;i<gauss.n_quadrature_points;++i)
    {
      double x = 2.*gauss.point(i)(0)-1.;
      double P1 = p1.value(x);
      double P2 = p2.value(x);
      sum += 2.*gauss.weight(i) * P1 * P2;
    }
  return sum;
}

int main ()
{
  std::ofstream logfile("polynomial1d.output");
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

  std::vector<Polynomial<double> > p (15);

  deallog << "Legendre" << endl;
  
  for (unsigned int i=0;i<p.size();++i)
    p[i] = Legendre<double>(i);
  
  for (unsigned int i=0;i<p.size();++i)
    for (unsigned int j=0;j<=i;++j)
      deallog << 'P' << i << " * P" << j
	      << " =" << scalar_product(p[i], p[j]) << std::endl;


  deallog << "LagrangeEquidistant" << endl;
  
  p.resize(6);
  for (unsigned int i=0;i<p.size();++i)
    p[i] = LagrangeEquidistant(p.size(), i);

				   // We add 1.0001 bacuse of bugs in
				   // the ostream classes
  for (unsigned int i=0;i<p.size();++i)
    for (unsigned int j=0;j<p.size();++j)
      deallog << 'P' << i << "(x" << j
	      << ") =" << p[i].value((double) j/p.size())+1.0001 << std::endl;
}

