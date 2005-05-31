//----------------------------  quadrature_test.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature_test.cc  ---------------------------

// Test projection onto lines

#include "../tests.h"
#include <iostream>
#include <fstream>

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/qprojector.h>
#include <cmath>

template<int dim>
void check_line(Quadrature<1>& quadrature)
{
  Point<dim> p1;
  Point<dim> p2;
  p1(0) = 1.;
  p2(0) = 7.;
  if (dim>1)
    {
      p1(1) = 3;
      p2(1) = -5.;
    }
  if (dim>2)
    {
      p1(2) = 0;
      p2(2) = 10.;
    }
  Quadrature<dim> q = QProjector<dim>::project_to_line(quadrature, p1, p2);
  double s = 0.;
  
  for (unsigned int k=0;k<q.n_quadrature_points;++k)
    {
      deallog << k << '\t' << q.point(k) << std::endl;
      s += q.weight(k);
    }
  deallog << "length: " << s << std::endl;
}


int main()
{
  std::ofstream logfile("qprojector.output");
  deallog.attach(logfile);
  deallog.depth_console(10);
  deallog.threshold_double(1.e-10);
  
  QTrapez<1> trapez;
  check_line<1> (trapez);
  check_line<2> (trapez);
  check_line<3> (trapez);
  
  QSimpson<1> simpson;
  check_line<1> (simpson);
  check_line<2> (simpson);
  check_line<3> (simpson);

  QMilne<1> milne;
  check_line<1> (milne);
  check_line<2> (milne);
  check_line<3> (milne);

  QWeddle<1> weddle;
  check_line<1> (weddle);
  check_line<2> (weddle);
  check_line<3> (weddle);
}
