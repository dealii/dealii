//-----------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2010, 2011, 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//-----------------------------------------------------------------------------

// synthesize the laplacian_list function out of repeated calls to
// laplacian

#include "../tests.h"
#include <deal.II/base/function.h>
#include <deal.II/lac/vector.h>


template <int dim>
class F : public Function<dim>
{
  public:
    double laplacian (const Point<dim> &p,
		      const unsigned int c) const
      {
	Assert (c == 0, ExcInternalError());
	return p.norm();
      }
};


template <int dim>
void check ()
{
  std::vector<Point<dim> > points;
  for (unsigned int i=0; i<10; ++i)
    {
      Point<dim> p;
      for (unsigned int d=0; d<dim; ++d)
	p[d] = d;
      points.push_back (p);
    }

  F<dim> f;
  std::vector<double> laplacians(10);
  f.laplacian_list (points, laplacians);

  for (unsigned int i=0; i<10; ++i)
    Assert (points[i].norm() == laplacians[i],
	    ExcInternalError());

  deallog << "OK" << std::endl;
}




int main()
{
  std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}


