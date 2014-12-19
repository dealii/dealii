// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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



// make sure that the QuadratureSelector works for a selection of
// arguments


#include "../tests.h"
#include <fstream>

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/quadrature_selector.h>
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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  check ("gauss", 2, QGauss<1>(2));
  check ("gauss", 2, QGauss<2>(2));
  check ("gauss", 2, QGauss<3>(2));

  check ("gauss", 2, QGauss<3>(2));
  check ("gauss", 6, QGauss<3>(6));
  check ("gauss", 10, QGauss<3>(10));

  check ("weddle", 0, QWeddle<2>());
}


