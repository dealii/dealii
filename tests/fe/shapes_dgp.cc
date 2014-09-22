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
#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_DGP_shape_functions()
{
  MappingQ1<dim> m;

  FE_DGP<dim> p1(1);
  plot_shape_functions(m, p1, "DGP1");
  plot_face_shape_functions(m, p1, "DGP1");
  test_compute_functions(m, p1, "DGP1");

  FE_DGP<dim> p2(2);
  plot_shape_functions(m, p2, "DGP2");
  plot_face_shape_functions(m, p2, "DGP2");
  test_compute_functions(m, p2, "DGP2");

  FE_DGP<dim> p3(3);
  plot_shape_functions(m, p3, "DGP3");
  plot_face_shape_functions(m, p3, "DGP3");
  test_compute_functions(m, p3, "DGP3");

//    FE_DGP<dim> p4(4);
//    plot_shape_functions(m, p4, "DGP4");
//    plot_face_shape_functions(m, p4, "DGP4");
//    test_compute_functions(m, p4, "DGP4");

//    FE_DGP<dim> p5(5);
//    plot_shape_functions(m, p5, "DGP5");
//    FE_DGP<dim> p6(6);
//    plot_shape_functions(m, p6, "DGP6");
//    FE_DGP<dim> p7(7);
//    plot_shape_functions(m, p7, "DGP7");
//    FE_DGP<dim> p8(8);
//    plot_shape_functions(m, p8, "DGP8");
//    FE_DGP<dim> p9(9);
//    plot_shape_functions(m, p9, "DGP9");
//    FE_DGP<dim> p10(10);
//    plot_shape_functions(m, p10, "DGP10");
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_DGP_shape_functions<1>();
  plot_FE_DGP_shape_functions<2>();
  plot_FE_DGP_shape_functions<3>();

  return 0;
}



