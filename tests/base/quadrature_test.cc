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

  vector<Quadrature<1> *> quadratures(14);
  quadratures[2]=new QGauss2<1>();
  quadratures[3]=new QGauss3<1>();
  quadratures[4]=new QGauss4<1>();
  quadratures[5]=new QGauss5<1>();
  quadratures[6]=new QGauss6<1>();
  quadratures[7]=new QGauss7<1>();
  quadratures[8]=new QGauss8<1>();
  quadratures[9] = new QMidpoint<1>();
  quadratures[10] = new QTrapez<1>();
  quadratures[11] = new QSimpson<1>();
  quadratures[12] = new QMilne<1>();
  quadratures[13] = new QWeddle<1>();
  
  for (unsigned int n=2; n<14; ++n)
    {
      Quadrature<1> *quadrature=quadratures[n];
      const vector<Point<1> > &points=quadrature->get_points();
      const vector<double> &weights=quadrature->get_weights();

      deallog << "Quadrature no." << n;
      
      unsigned int i=0;
      double quadrature_int=0;
      double exact_int=0;
      double err = 0;
      do
	{
	  ++i;

	  quadrature_int=0;
	  for (unsigned int x=0; x<quadrature->n_quadrature_points; ++x)
	    quadrature_int+=pow(points[x](0), i)*weights[x];

					   // the exact integral is 1/(i+1)
	  exact_int=1./(i+1);
	  err = fabs(quadrature_int-exact_int);
	}
      while (err<1e-15);

				       // Uncomment here for testing
//      deallog << " (Error " << err << ")";
      deallog << " is exact for polynomials of degree " << i-1 << endl;
    }
}


