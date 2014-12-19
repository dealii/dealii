// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


#include "interpolate_common.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>

#include <fstream>

// FE_Q_DG0<dim>::interpolate(...)

template <int dim>
void check_q_dg0(const Function<dim> &f,
                 const unsigned int degree)
{
  FE_Q_DG0<dim> fe(degree);
  deallog << fe.get_name() << ' ';

  std::vector<double> dofs(fe.dofs_per_cell);

  std::vector<std::vector<double> >
  values(1, std::vector<double>(fe.get_unit_support_points().size()));
  f.value_list(fe.get_unit_support_points(), values[0]);
  fe.interpolate(dofs, values[0]);
  deallog << " value " << difference(fe,dofs,f);
  fe.interpolate(dofs, values);
  deallog << " vector " << difference(fe,dofs,f);

  std::vector<Vector<double> >
  vectors(fe.get_unit_support_points().size(), Vector<double>(1));
  f.vector_value_list(fe.get_unit_support_points(), vectors);
  fe.interpolate(dofs, vectors, 0);
  deallog << " Vector " << difference(fe,dofs,f) << std::endl;
}

int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(2.e-15);

  Q1WedgeFunction<1,1> w1;
  check_q_dg0(w1,1);
  check_q_dg0(w1,2);
  check_q_dg0(w1,3);
  Q1WedgeFunction<2,1> w2;
  check_q_dg0(w2,1);
  check_q_dg0(w2,2);
  check_q_dg0(w2,3);
  Q1WedgeFunction<2,2> w22;
  check_q_dg0(w22,2);
  check_q_dg0(w22,3);
  Q1WedgeFunction<2,3> w23;
  check_q_dg0(w23,3);
  Q1WedgeFunction<3,1> w3;
  check_q_dg0(w3,1);
  check_q_dg0(w3,2);
  check_q_dg0(w3,3);
}
