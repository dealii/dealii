// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_trace.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_TraceQ_shape_functions()
{
  MappingQGeneric<dim> m(1);

  FE_TraceQ<dim> tq1(1);
  FE_TraceQ<dim> tq2(2);
  FE_TraceQ<dim> tq3(3);
  FE_TraceQ<dim> tq4(4);
  FE_Q<dim>      q1(1);
  FE_Q<dim>      q2(2);
  FESystem<dim>  sys1(tq1, 1);
  FESystem<dim>  sys2(tq2, 1);
  FESystem<dim>  sys11(tq1, 1, q1, 1);
  FESystem<dim>  sys22(tq2, 1, q2, 1);
  plot_face_shape_functions(m, tq1, "TraceQ1", update_values);
  plot_face_shape_functions(m, tq2, "TraceQ2", update_values);
  plot_face_shape_functions(m, tq3, "TraceQ3", update_values);
  plot_face_shape_functions(m, tq4, "TraceQ4", update_values);
  plot_face_shape_functions(m, sys1, "System1", update_values);
  plot_face_shape_functions(m, sys2, "System2", update_values);
  plot_face_shape_functions(m, sys11, "System1-1", update_values);
  plot_face_shape_functions(m, sys22, "System2-2", update_values);
}


int
main()
{
  const std::string logname = "output";
  std::ofstream     logfile(logname.c_str());
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_TraceQ_shape_functions<2>();
  plot_FE_TraceQ_shape_functions<3>();

  return 0;
}
