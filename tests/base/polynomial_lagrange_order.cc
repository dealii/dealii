//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// Similar to polyomial_lagrange, but test Lagrange interpolation for high
// order with tighter tolerances, in particular the effect of stability of the
// polynomial evaluation at random points

#include "../tests.h"
#include <iomanip>
#include <fstream>
#include <cmath>

#include <deal.II/base/logstream.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>


using namespace Polynomials;


void check_interpolation (const std::vector<Polynomial<double> >& p,
			  const std::vector<Point<1> >& x)
{
  for (unsigned int i=0;i<p.size();++i)
    {
      deallog << i;
      for (unsigned int k=0;k<x.size();++k)
	{
	  deallog << '.';
	  const double y = p[i].value(x[k](0));
	  if (i == k)
	    {
	      if (std::fabs(y-1.) > 1e-13)
		deallog << "Error1  lg y=" << std::log10(std::fabs(y-1.))
			<< std::endl;
	    }
	  else
	    {
	      if (std::fabs(y) > 1e-13)
		deallog << "Error0  lg y=" << std::log10(std::fabs(y))
			<< std::endl;
	    }
	}
      deallog << std::endl;
    }
}



void check_constant (const std::vector<Polynomial<double> >& p)
{
				// check whether the sum of all polynomials in
				// the basis gives one for a given point
  deallog << "Representation of one at random points";
  for (unsigned int j=0; j<12; ++j)
    {
      double x = (double)rand()/RAND_MAX;
      double value = 0;
      for (unsigned int i=0;i<p.size();++i)
	{
	  value += p[i].value(x);
	}
      deallog << ".";
      if (std::fabs (1.-value) > 1e-13)
	deallog << "Error1  lg y=" << std::log10(std::fabs(1.-value))
		<< std::endl;
    }
  deallog << std::endl;

}



void
check_poly (const Quadrature<1>& q)
{
  deallog << "Points: " << q.size() << std::endl;
  std::vector<Polynomial<double> > p = generate_complete_Lagrange_basis(q.get_points());
  check_interpolation(p, q.get_points());
  check_constant (p);
  deallog << std::endl;
}



void
check_lge (unsigned int n)
{
  deallog << "Points: " << n+1 << std::endl;
  std::vector<Polynomial<double> > p = LagrangeEquidistant::generate_complete_basis(n);
  std::vector<Point<1> > x(n+1);
  const double h = 1./n;
  for (unsigned int i=0;i<=n;++i)
    x[i](0) = h*i;
  check_interpolation(p, x);
  check_constant(p);
  deallog << std::endl;
}



int main()
{
  std::ofstream logfile("polynomial_lagrange_order/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  deallog.push("LagrangeEquidistant");
  for (unsigned i=8; i<18; i+=2)
    check_lge (i);
  deallog.pop();
  deallog << std::endl;

				// Lagrange elements on GL points have good
				// conditioning, so test to some very high
				// orders
  deallog.push("GaussLobatto");
  for (unsigned i=8; i<40; i+=3)
    check_poly (QGaussLobatto<1>(i+1));
  deallog.pop();
}
