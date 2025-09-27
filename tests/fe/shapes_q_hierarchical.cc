// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/fe/fe_q_hierarchical.h>
#include <deal.II/fe/mapping_q1.h>

#include <string>

#include "../tests.h"

#include "shapes.h"

#define PRECISION 8


template <int dim>
void
plot_FE_Q_Hierarchical_shape_functions()
{
  MappingQ<dim> m(1);

  FE_Q_Hierarchical<dim> q1(1);
  plot_shape_functions(m, q1, "QHierarchical1");
  plot_face_shape_functions(m, q1, "QHierarchical1");
  test_compute_functions(m, q1, "QHierarchical1");

  FE_Q_Hierarchical<dim> q2(2);
  plot_shape_functions(m, q2, "QHierarchical2");
  plot_face_shape_functions(m, q2, "QHierarchical2");
  test_compute_functions(m, q2, "QHierarchical2");

  // skip the following tests to
  // reduce run-time
  if (dim < 3)
    {
      FE_Q_Hierarchical<dim> q3(3);
      plot_shape_functions(m, q3, "QHierarchical3");
      plot_face_shape_functions(m, q3, "QHierarchical3");
      test_compute_functions(m, q3, "QHierarchical3");

      FE_Q_Hierarchical<dim> q4(4);
      plot_shape_functions(m, q4, "QHierarchical4");
      plot_face_shape_functions(m, q4, "QHierarchical4");
      test_compute_functions(m, q4, "QHierarchical4");
    }
}


int
main()
{
  std::ofstream logfile("output");
  deallog << std::setprecision(PRECISION) << std::fixed;
  deallog.attach(logfile);

  plot_FE_Q_Hierarchical_shape_functions<1>();
  plot_FE_Q_Hierarchical_shape_functions<2>();
  plot_FE_Q_Hierarchical_shape_functions<3>();

  return 0;
}
