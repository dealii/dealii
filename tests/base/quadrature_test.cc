//----------------------------  quadrature_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature_test.cc  ---------------------------


#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <cmath>

int main(int,char)
{
  ofstream logfile("quadrature_test.output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  vector<Quadrature<2> *> quadratures;
  quadratures.push_back (new QGauss2<2>());
  quadratures.push_back (new QGauss3<2>());
  quadratures.push_back (new QGauss4<2>());
  quadratures.push_back (new QGauss5<2>());
  quadratures.push_back (new QGauss6<2>());
  quadratures.push_back (new QGauss7<2>());
  quadratures.push_back (new QMidpoint<2>());
  quadratures.push_back (new QTrapez<2>());
  quadratures.push_back (new QSimpson<2>());
  quadratures.push_back (new QMilne<2>());
  quadratures.push_back (new QWeddle<2>());
  
  for (unsigned int n=0; n<quadratures.size(); ++n)
    {
      Quadrature<2> *quadrature=quadratures[n];
      const vector<Point<2> > &points=quadrature->get_points();
      const vector<double> &weights=quadrature->get_weights();

      deallog << "Quadrature no." << n
	      << " (" << typeid(*quadrature).name() << ")";
      
      unsigned int i=0;
      double quadrature_int=0;
      double exact_int=0;
      double err = 0;
      do
	{
	  ++i;

	  quadrature_int=0;

					   // Check the polynomial x^i*y^i

	  for (unsigned int x=0; x<quadrature->n_quadrature_points; ++x)
	    quadrature_int+=pow(points[x](0), i)*pow(points[x](1), i)*weights[x];

					   // the exact integral is 1/(i+1)
	  exact_int=1./(i+1)/(i+1);
	  err = fabs(quadrature_int-exact_int);
	}
      while (err<1e-15);

				       // Uncomment here for testing
//      deallog << " (Error " << err << ")";
      deallog << " is exact for polynomials of degree " << i-1 << endl;
    }
}


