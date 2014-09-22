// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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


#include "../tests.h"
#include "shapes.h"
#include <deal.II/fe/fe_dgp_monomial.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_DGPMonomial_shape_functions()
{
  MappingQ1<dim> m;

  FE_DGPMonomial<dim> p1(1);
  plot_shape_functions(m, p1, "DGPMonomial1");
  plot_face_shape_functions(m, p1, "DGPMonomial1");
  test_compute_functions(m, p1, "DGPMonomial1");

  FE_DGPMonomial<dim> p2(2);
  plot_shape_functions(m, p2, "DGPMonomial2");
  plot_face_shape_functions(m, p2, "DGPMonomial2");
  test_compute_functions(m, p2, "DGPMonomial2");

  if (dim<3)
    {
      FE_DGPMonomial<dim> p3(3);
      plot_shape_functions(m, p3, "DGPMonomial3");
      plot_face_shape_functions(m, p3, "DGPMonomial3");
      test_compute_functions(m, p3, "DGPMonomial3");
    }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_DGPMonomial_shape_functions<1>();
  plot_FE_DGPMonomial_shape_functions<2>();
  plot_FE_DGPMonomial_shape_functions<3>();

  return 0;
}



