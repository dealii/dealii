/* $Id$   */
/*            Ralf Hartmann, University of Heidelberg, Feb 99         */

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

  vector<Quadrature<1> *> gauss_quadratures(9);
  gauss_quadratures[2]=new QGauss2<1>();
  gauss_quadratures[3]=new QGauss3<1>();
  gauss_quadratures[4]=new QGauss4<1>();
  gauss_quadratures[5]=new QGauss5<1>();
  gauss_quadratures[6]=new QGauss6<1>();
  gauss_quadratures[7]=new QGauss7<1>();
  gauss_quadratures[8]=new QGauss8<1>();
  
  for (unsigned int n=2; n<9; ++n)
    {
      Quadrature<1> *quadrature=gauss_quadratures[n];
      const vector<Point<1> > &points=quadrature->get_quad_points();
      const vector<double> &weights=quadrature->get_weights();

      deallog << "Gauss" << n;
      
      unsigned int i=0;
      double quadrature_int=0, exact_int=0;
      do
	{
	  ++i;

	  quadrature_int=0;
	  for (unsigned int x=0; x<quadrature->n_quadrature_points; ++x)
	    quadrature_int+=pow(points[x](0), i)*weights[x];

					   // the exact integral is 1/(i+1)
	  exact_int=1./(i+1);
	}
      while (fabs(quadrature_int-exact_int)<1e-15);

      deallog << " is exact for polynomials of degree " << i-1 << endl;
    }
}



