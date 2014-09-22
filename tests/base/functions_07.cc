// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------


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
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check<1> ();
  check<2> ();
  check<3> ();
}


