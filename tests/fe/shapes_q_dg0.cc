// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2014 by the deal.II authors
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
#include <deal.II/fe/fe_q_dg0.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 2


template<int dim>
void plot_FE_Q_DG0_shape_functions()
{
  MappingQ1<dim> m;

  FE_Q_DG0<dim> q1(1);
  plot_shape_functions(m, q1, "Q1_DG0");
  plot_face_shape_functions(m, q1, "Q1_DG0");
  test_compute_functions(m, q1, "Q1_DG0");

  FE_Q_DG0<dim> q2(2);
  plot_shape_functions(m, q2, "Q2_DG0");
  plot_face_shape_functions(m, q2, "Q2_DG0");
  test_compute_functions(m, q2, "Q2_DG0");

  // skip the following tests to
  // reduce run-time
  if (dim < 3)
    {
      FE_Q_DG0<dim> q3(3);
      plot_shape_functions(m, q3, "Q3_DG0");
      plot_face_shape_functions(m, q3, "Q3_DG0");
      test_compute_functions(m, q3, "Q3_DG0");

      FE_Q_DG0<dim> q4(4);
      plot_shape_functions(m, q4, "Q4_DG0");
      plot_face_shape_functions(m, q4, "Q4_DG0");
      test_compute_functions(m, q4, "Q4_DG0");
    };

//    FE_Q<dim> q5(5);
//    plot_shape_functions(m, q5, "Q5");
//    FE_Q<dim> q6(6);
//    plot_shape_functions(m, q6, "Q6");
//    FE_Q<dim> q7(7);
//    plot_shape_functions(m, q7, "Q7");
//    FE_Q<dim> q8(8);
//    plot_shape_functions(m, q8, "Q8");
//    FE_Q<dim> q9(9);
//    plot_shape_functions(m, q9, "Q9");
//    FE_Q<dim> q10(10);
//    plot_shape_functions(m, q10, "Q10");
}


int
main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  plot_FE_Q_DG0_shape_functions<1>();
  plot_FE_Q_DG0_shape_functions<2>();
  plot_FE_Q_DG0_shape_functions<3>();

  return 0;
}



