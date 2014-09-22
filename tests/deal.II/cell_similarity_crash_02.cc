// ---------------------------------------------------------------------
//
// Copyright (C) 2009 - 2013 by the deal.II authors
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



// The first attempt to fix the _01 test had a bug in it in that we
// managed to invalidate the cell stored in FEValues but forgot to
// disconnect from the triangulation's signal. We were therefore still
// getting reminders from the triangulation even though we weren't
// doing anything with it any more

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>

#include <fstream>



template <int dim>
void test()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global (2);

  FE_Q<dim> fe(1);
  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (fe, quadrature, update_values);

  // initialize FEValues with the first cell
  fe_values.reinit (tr.begin_active());

  // then invalidate the cell iterator
  tr.refine_global (1);
  // and invalidate it again. this
  // shouldn't do any further harm
  tr.refine_global (1);

  deallog << "OK" << std::endl;
}


int main()
{
  initlog();
  deallog << std::setprecision (4);
  deallog.threshold_double(1.e-7);

  test<1>();
  test<2>();
  test<3>();
}
