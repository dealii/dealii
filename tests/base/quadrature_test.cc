//----------------------------  quadrature_test.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
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

template <int dim>
void
fill_vector (vector<Quadrature<dim> *>& quadratures)
{
  quadratures.push_back (new QGauss2<dim>());
  quadratures.push_back (new QGauss3<dim>());
  quadratures.push_back (new QGauss4<dim>());
  quadratures.push_back (new QGauss5<dim>());
  quadratures.push_back (new QGauss6<dim>());
  quadratures.push_back (new QGauss7<dim>());
  quadratures.push_back (new QMidpoint<dim>());
  quadratures.push_back (new QTrapez<dim>());
  quadratures.push_back (new QSimpson<dim>());
  quadratures.push_back (new QMilne<dim>());
  quadratures.push_back (new QWeddle<dim>());
  for (unsigned int i=1;i<9;++i)
    {
      quadratures.push_back (new QGauss<dim>(i));
    }
}

template <int dim>
void
check_cells (vector<Quadrature<dim>*>& quadratures)
{
  for (unsigned int n=0; n<quadratures.size(); ++n)
    {
      Quadrature<dim>& quadrature = *quadratures[n];
      const vector<Point<dim> > &points=quadrature.get_points();
      const vector<double> &weights=quadrature.get_weights();

      deallog << "Quadrature no." << n
	      << " (" << typeid(*quadratures[n]).name() << ")";

      unsigned int i=0;
      double quadrature_int=0;
      double exact_int=0;
      double err = 0;
      
      do
	{
	  ++i;

	  quadrature_int=0;
					   // Check the polynomial x^i*y^i

	  for (unsigned int x=0; x<quadrature.n_quadrature_points; ++x)
	    {
	      double f=1.;
	      switch (dim)
		{
		case 3:
		  f *= pow(points[x](2), i);
		case 2:
		  f *= pow(points[x](1), i);
		case 1:
		  f *= pow(points[x](0), i);
		}
	      quadrature_int+=f*weights[x];
	    }
	  
					   // the exact integral is 1/(i+1)
	  exact_int=1./pow(i+1,dim);	  
	  err = fabs(quadrature_int-exact_int);
	}
      while (err<1e-15);
				       // Uncomment here for testing
//      deallog << " (Int " << quadrature_int << ',' << exact_int << ")";
      deallog << " is exact for polynomials of degree " << i-1 << endl;

      if (dim==1)
	{
					   // check the ordering of
					   // the quadrature points
	  bool in_order=true;
	  for (unsigned int x=1; x<quadrature.n_quadrature_points; ++x)
	    {
	      if (points[x](0)<=points[x-1](0))
		in_order=false;
	    }
	  if (!in_order)
	    for (unsigned int x=0; x<quadrature.n_quadrature_points; ++x)
	      deallog << points[x] << endl;
	}
    }  
}


template <int dim>
void
check_faces (vector<Quadrature<dim-1>*>& quadratures, bool sub)
{
  if (sub)
    deallog.push("subfaces");
  else
    deallog.push("faces");

  for (unsigned int n=0; n<quadratures.size(); ++n)
    {
      QProjector<dim> quadrature(*quadratures[n], sub);
      const vector<Point<dim> > &points=quadrature.get_points();
      const vector<double> &weights=quadrature.get_weights();

      deallog << "Quadrature no." << n
	      << " (" << typeid(*quadratures[n]).name() << ")";

      unsigned int i=0;
      double quadrature_int=0;
      double exact_int=0;
      double err = 0;

      do
	{
	  ++i;

	  quadrature_int=0;
					   // Check the polynomial
	                                   // x^i*y^i*z^i

	  for (unsigned int x=0; x<quadrature.n_quadrature_points; ++x)
	    {
	      double f=1.;
	      switch (dim)
		{
		case 3:
		  f *= pow(points[x](2), i);
		case 2:
		  f *= pow(points[x](1), i);
		case 1:
		  f *= pow(points[x](0), i);
		}
	      quadrature_int+=f*weights[x];
	    }
	  
					   // the exact integral is
	                                   // 1/(i+1)^(dim-1)
	  switch (dim)
	    {
	    case 2:
	      exact_int = 2 * (sub ? 2:1) / (double) (i+1);
	      break;
	    case 3:
	      exact_int = 3 * (sub ? 4:1) / (double) (i+1)/(i+1);
	      break;
	    }
      
	  err = fabs(quadrature_int-exact_int);
	}
      while (err<5e-15);
				       // Uncomment here for testing
      //      deallog << " (Int " << quadrature_int << '-' << exact_int << '=' << err << ")";
      deallog << " is exact for polynomials of degree " << i-1 << endl;
    }
  deallog.pop();
}

int main(int,char)
{
  ofstream logfile("quadrature_test.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  vector<Quadrature<1> *> q1;
  vector<Quadrature<2> *> q2;
  vector<Quadrature<3> *> q3;
  fill_vector (q1);
  fill_vector (q2);
  fill_vector (q3);

  deallog.push("1d");
  check_cells(q1);
  deallog.pop();

  deallog.push("2d");
  check_cells(q2);
  check_faces(q1,false);
  check_faces(q1,true);
  deallog.pop();

  deallog.push("3d");
  check_cells(q3);
  check_faces(q2,false);
  check_faces(q2,true);
  deallog.pop();

				   // delete objects again to avoid
				   // messages about memory leaks
  for (unsigned int i=0; i<q1.size(); ++i)
    delete q1[i];
  for (unsigned int i=0; i<q2.size(); ++i)
    delete q2[i];
  for (unsigned int i=0; i<q3.size(); ++i)
    delete q3[i];  
}


