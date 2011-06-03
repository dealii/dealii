//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


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
  std::ofstream logfile ("cell_similarity_crash_02/output");
  deallog << std::setprecision (4);

  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-7);

  test<1>();
  test<2>();
  test<3>();
}
