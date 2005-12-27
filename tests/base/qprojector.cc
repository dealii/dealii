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
#include <base/geometry_info.h>
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

template<int dim>
void check_face(Quadrature<1>& q1)
{
  deallog << "Checking dim " << dim
	  << " 1d-points " << q1.n_quadrature_points
	  << std::endl;
  
  Quadrature<dim-1> subquadrature(q1);
  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    {
      deallog << "Face " << f
	      << std::endl;
      
      Quadrature<dim> quadrature
	= QProjector<dim>::project_to_face(subquadrature, f);
      for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
	deallog << quadrature.point(k) << std::endl;
    }

  for (unsigned int f=0;f<GeometryInfo<dim>::faces_per_cell;++f)
    for (unsigned int s=0;s<GeometryInfo<dim>::subfaces_per_face;++s)
      {
	deallog << "Face " << f << " subface " << s
		<< std::endl;
	
	Quadrature<dim> quadrature
	  = QProjector<dim>::project_to_face(subquadrature, f);
	for (unsigned int k=0;k<quadrature.n_quadrature_points;++k)
	  deallog << quadrature.point(k) << std::endl;
      }
}

void check(Quadrature<1>& q)
{
  check_line<1> (q);
  check_line<2> (q);
  check_line<3> (q);
  
  check_face<2>(q);
  check_face<3>(q);  
}

int main()
{
  std::ofstream logfile("qprojector/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  logfile.precision(2);
  
  deallog.threshold_double(1.e-10);

  Quadrature<1> none(0);
  check(none);
  
  QTrapez<1> trapez;
  check(trapez);
  
  QSimpson<1> simpson;
  check(simpson);

  QMilne<1> milne;
  check(milne);

  QWeddle<1> weddle;
  check(weddle);
}
