//----------------------------  gauss_lobatto.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  gauss_lobatto.cc  ---------------------------


// check accuracy of various quadrature formulas by using them to
// integrate polynomials of increasing degree, and finding the degree
// until which they integrate exactly


#include "../tests.h"
#include <iomanip>
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/qprojector.h>
#include <cmath>

template <int dim>
void
fill_vector (std::vector<Quadrature<dim> *>& quadratures)
{
  quadratures.push_back (new QMidpoint<dim>());
  quadratures.push_back (new QTrapez<dim>());
  quadratures.push_back (new QSimpson<dim>());
  quadratures.push_back (new QMilne<dim>());
  quadratures.push_back (new QWeddle<dim>());
  for (unsigned int i=0;i<9;++i)
    {
      quadratures.push_back (new QGauss<dim>(i));
    }
  QMilne<1> q1d;
  quadratures.push_back (new Quadrature<dim>(q1d));
  for (unsigned int i=2;i<8;++i)
    {
      quadratures.push_back (new QGaussLobatto<dim>(i));
    }
}

template <int dim>
void
check_cells (std::vector<Quadrature<dim>*>& quadratures)
{
  Quadrature<dim> quadrature;
  for (unsigned int n=0; n<quadratures.size(); ++n)
    {
      quadrature = *quadratures[n];
      const std::vector<Point<dim> > &points=quadrature.get_points();
      const std::vector<double> &weights=quadrature.get_weights();

      deallog << "Quadrature no." << n;

      unsigned int i=0;
      double quadrature_int=0;
      double exact_int=0;
      double err = 0;

      do
	{
	  ++i;

	  quadrature_int=0;
					   // Check the polynomial x^i*y^i

	  for (unsigned int x=0; x<quadrature.size(); ++x)
	    {
	      double f=1.;
	      switch (dim)
		{
		case 3:
		  f *= std::pow(static_cast<double>(points[x](2)), i*1.0);
		case 2:
		  f *= std::pow(static_cast<double>(points[x](1)), i*1.0);
		case 1:
		  f *= std::pow(static_cast<double>(points[x](0)), i*1.0);
		}
	      quadrature_int+=f*weights[x];
	    }

					   // the exact integral is 1/(i+1)
	  exact_int=1./std::pow(static_cast<double>(i+1),dim);
	  err = std::fabs(quadrature_int-exact_int);
	}
      while (err<1e-14);
				       // Uncomment here for testing
//      deallog << " (Int " << quadrature_int << ',' << exact_int << ")";
      deallog << " is exact for polynomials of degree " << i-1 << std::endl;

      if (dim==1)
	{
					   // check the ordering of
					   // the quadrature points
	  bool in_order=true;
	  for (unsigned int x=1; x<quadrature.size(); ++x)
	    {
	      if (points[x](0)<=points[x-1](0))
		in_order=false;
	    }
	  if (!in_order)
	    for (unsigned int x=0; x<quadrature.size(); ++x)
	      deallog << points[x] << std::endl;
	}
    }
}


template <int dim>
void
check_faces (const std::vector<Quadrature<dim-1>*>& quadratures, const bool sub)
{
  if (sub)
    deallog.push("subfaces");
  else
    deallog.push("faces");

  for (unsigned int n=0; n<quadratures.size(); ++n)
    {
      Quadrature<dim> quadrature (sub == false?
				  QProjector<dim>::project_to_all_faces(*quadratures[n]) :
				  QProjector<dim>::project_to_all_subfaces(*quadratures[n]));
      const std::vector<Point<dim> > &points=quadrature.get_points();
      const std::vector<double> &weights=quadrature.get_weights();

      deallog << "Quadrature no." << n;

      unsigned int i=0;
      long double quadrature_int=0;
      double exact_int=0;
      double err = 0;

      do
	{
	  ++i;

	  quadrature_int=0;
					   // Check the polynomial
	                                   // x^i*y^i*z^i

	  for (unsigned int x=0; x<quadrature.size(); ++x)
	    {
	      long double f=1.;
	      switch (dim)
		{
		case 3:
		  f *= std::pow((long double) points[x](2), i*1.0L);
		case 2:
		  f *= std::pow((long double) points[x](1), i*1.0L);
		case 1:
		  f *= std::pow((long double) points[x](0), i*1.0L);
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
	      exact_int = 3 * (sub ? (4+2+2):1)*8 / (double) (i+1)/(i+1);
	      break;
	    }

	  err = std::fabs(quadrature_int-exact_int);
	}
				       // for comparison: use factor 8 in case
				       // of dim==3, as we integrate 8 times
				       // over the whole surface (all
				       // combinations of face_orientation,
				       // face_flip and face_rotation)
      while (err < (dim==3 ? 8 : 1) * 2e-14);
				       // Uncomment here for testing
      //      deallog << " (Int " << quadrature_int << '-' << exact_int << '=' << err << ")";
      deallog << " is exact for polynomials of degree " << i-1 << std::endl;
    }
  deallog.pop();
}

int main()
{
  std::ofstream logfile("gauss_lobatto/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int n=2; n<20; ++n)
    {
      deallog << "QGaussLobatto(" << n << ")" << std::endl;

      QGaussLobatto<1> q(n);
      for (unsigned int i=0; i<q.size(); ++i)
	deallog << q.point(i) << ' ' << q.weight(i) << std::endl;

				       // the points must be
				       // symmetrically located around
				       // 0.5
      double p = 0;
      for (unsigned int i=0; i<q.size(); ++i)
	p += (q.point(i)[0] - 0.5);
      Assert (std::fabs(p) < 1e-12, ExcInternalError());

				       // the sum of weights must be
				       // one
      double w = 0;
      for (unsigned int i=0; i<q.size(); ++i)
	w += q.weight(i);
      Assert (std::fabs(w-1) < 1e-12, ExcInternalError());
    }
}


