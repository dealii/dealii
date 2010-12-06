//--------------------------  quadrature_qiterated.cc  --------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------  quadrature_qiterated.cc  --------------------

// check serialization for QIterated

#include "serialization.h"
#include <base/quadrature.h>
#include <boost/serialization/vector.hpp>

void test ()
{ 
  const unsigned int dim = 2;
  unsigned int n_copies = 3;
  
  std::vector<Point <1> > points1;
  points1.push_back(Point<1>(0.));
  points1.push_back(Point<1>(1.));
  double w1[2] = {0.5, 0.5};
  std::vector<double> weights1(w1, &w1[2]);
  
  std::vector<Point <1> > points2;
  points2.push_back(Point<1>(0.25));
  points2.push_back(Point<1>(0.75));
  double w2[2] = {0.4, 0.6}; 
  std::vector<double> weights2(w2, &w2[2]);
  
  Quadrature<1> qx(points1, weights1);
  Quadrature<1> qy(points2, weights2);
  
  QIterated<dim> q1(qx, n_copies);
   
  QIterated<dim> q2(qy, n_copies);

  verify (q1, q2);
}


int main ()
{
  std::ofstream logfile("quadrature_qiterated/output");
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test ();

  deallog << "OK" << std::endl;
}
