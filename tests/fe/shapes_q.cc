// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2015 by the deal.II authors
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
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q1.h>
#include <fstream>
#include <string>

#define PRECISION 8


template<int dim>
void plot_FE_Q_shape_functions()
{
  MappingQGeneric<dim> m(1);

  FE_Q<dim> q1(1);
  plot_shape_functions(m, q1, "Q1");
  plot_face_shape_functions(m, q1, "Q1");
  test_compute_functions(m, q1, "Q1");

  FE_Q<dim> q2(2);
  plot_shape_functions(m, q2, "Q2");
  plot_face_shape_functions(m, q2, "Q2");
  test_compute_functions(m, q2, "Q2");

  // skip the following tests to
  // reduce run-time
  if (dim < 3)
    {
      FE_Q<dim> q3(QIterated<1>(QTrapez<1>(),3));
      plot_shape_functions(m, q3, "Q3");
      plot_face_shape_functions(m, q3, "Q3");
      test_compute_functions(m, q3, "Q3");

      FE_Q<dim> q4(QIterated<1>(QTrapez<1>(),4));
      plot_shape_functions(m, q4, "Q4");
      plot_face_shape_functions(m, q4, "Q4");
      test_compute_functions(m, q4, "Q4");
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
  deallog.threshold_double(1.e-10);

  plot_FE_Q_shape_functions<1>();
  plot_FE_Q_shape_functions<2>();
  plot_FE_Q_shape_functions<3>();

  return 0;
}
