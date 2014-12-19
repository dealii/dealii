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
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_DGQ_shape_functions()
{
  MappingQ1<dim> m;

  FE_DGQ<dim> q1(1);
  plot_shape_functions(m, q1, "DGQ1");
  plot_face_shape_functions(m, q1, "DGQ1");
  test_compute_functions(m, q1, "DGQ1");

  FE_DGQ<dim> q2(2);
  plot_shape_functions(m, q2, "DGQ2");
  plot_face_shape_functions(m, q2, "DGQ2");
  test_compute_functions(m, q2, "DGQ2");

  FE_DGQ<dim> q3(3);
  plot_shape_functions(m, q3, "DGQ3");
  plot_face_shape_functions(m, q3, "DGQ3");
  test_compute_functions(m, q3, "DGQ3");

  QGaussLobatto<1> quadrature_gl(5);
  FE_DGQArbitraryNodes<dim> qgl(quadrature_gl);
  plot_shape_functions(m, qgl, "DGQGL");
  plot_face_shape_functions(m, qgl, "DGQGL");
  test_compute_functions(m, qgl, "DGQGL");

  QGauss<1> quadrature_g(5);
  FE_DGQArbitraryNodes<dim> qg(quadrature_g);
  plot_shape_functions(m, qg, "DGQG");
  plot_face_shape_functions(m, qg, "DGQG");
  test_compute_functions(m, qg, "DGQG");

//    FE_DGQ<dim> q4(4);
//    plot_shape_functions(m, q4, "DGQ4");
//    plot_face_shape_functions(m, q4, "DGQ4");
//    test_compute_functions(m, q4, "DGQ4");

//    FE_DGQ<dim> q5(5);
//    plot_shape_functions(m, q5, "DGQ5");
//    FE_DGQ<dim> q6(6);
//    plot_shape_functions(m, q6, "DGQ6");
//    FE_DGQ<dim> q7(7);
//    plot_shape_functions(m, q7, "DGQ7");
//    FE_DGQ<dim> q8(8);
//    plot_shape_functions(m, q8, "DGQ8");
//    FE_DGQ<dim> q9(9);
//    plot_shape_functions(m, q9, "DGQ9");
//    FE_DGQ<dim> q10(10);
//    plot_shape_functions(m, q10, "DGQ10");
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_DGQ_shape_functions<1>();
  plot_FE_DGQ_shape_functions<2>();
  plot_FE_DGQ_shape_functions<3>();

  return 0;
}



