//----------------------------  quadrature_selector.cc  ---------------------------
//    quadrature_selector.cc,v 1.18 2003/01/08 17:58:18 wolf Exp
//    Version: 
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  quadrature_selector.cc  ---------------------------


// make sure that the QuadratureSelector works for a selection of
// arguments


#include "../tests.h"
#include <fstream>

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <base/quadrature_selector.h>
#include <string>
#include <cmath>


template <int dim>
void check (const std::string     &name,
	    const unsigned int     order,
	    const Quadrature<dim> &q)
{
  Assert (QuadratureSelector<dim>(name, order).get_points() ==
	  q.get_points(),
	  ExcInternalError());
  deallog << name << ' ' << order << " ok" << std::endl;
}


int main()
{
  std::ofstream logfile("quadrature_selector.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  check ("gauss", 2, QGauss2<1>());
  check ("gauss", 2, QGauss2<2>());
  check ("gauss", 2, QGauss2<3>());
  
  check ("gauss", 2, QGauss2<3>());
  check ("gauss", 6, QGauss6<3>());
  check ("gauss", 10, QGauss<3>(10));

  check ("weddle", 0, QWeddle<2>());
}


