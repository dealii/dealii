// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// Check the documented property of the DGP element that the first shape
// function is constant in space

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>
#include <string>

#define PRECISION 3



template<int dim>
void
test(const unsigned int degree)
{
  FE_DGP<dim> fe(degree);
  deallog << fe.get_name() << std::endl;

  Triangulation<dim> tr;
  if (dim > 1)
    GridGenerator::hyper_ball (tr);
  else
    GridGenerator::hyper_cube (tr);

  QGauss<dim> q(degree+1);
  FEValues<dim> fe_values (fe, q, update_values);
  for (typename Triangulation<dim>::active_cell_iterator
       cell = tr.begin_active(); cell != tr.end(); ++cell)
    {
      fe_values.reinit (cell);
      for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
        Assert (fe_values.shape_value (0,q) == 1,
                ExcInternalError());
    }

  deallog << "OK" << std::endl;
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  for (unsigned int degree=0; degree<=4; ++degree)
    test<1>(degree);

  for (unsigned int degree=0; degree<=4; ++degree)
    test<2>(degree);

  for (unsigned int degree=0; degree<=4; ++degree)
    test<3>(degree);

  return 0;
}



