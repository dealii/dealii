//----------------------------  negative_jac_det.cc  -------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  negative_jac_det.cc  -------------------------

// Test that the computation of the Jacobian determinant is correct also when
// the mesh is in a non-standard orientation

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/tria.h>


void test()
{
  const int dim = 2;
  Triangulation<dim> tria;
  std::vector<Point<dim> > points (4);
  points[0] = (Point<dim> (0, 0));
  points[1] = (Point<dim> (0, 1));
  points[2] = (Point<dim> (1, 0));
  points[3] = (Point<dim> (1, 1));

  std::vector<CellData<dim> > cells (1);
  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 2;
  cells[0].vertices[3] = 3;
  cells[0].material_id = 0;

  tria.create_triangulation(points, cells, SubCellData());

  FE_Nothing<dim> dummy;
  QGauss<dim> quad(2);
  FEValues<dim> fe_val (dummy, quad, update_JxW_values);
  double integral = 0.;
  /*typename*/ Triangulation<dim>::active_cell_iterator
    cell = tria.begin_active(), endc = tria.end();
  for ( ; cell != endc; ++cell)
    {
      fe_val.reinit (cell);
      for (unsigned int q=0; q<quad.size(); ++q)
        integral += fe_val.JxW(q);
    }
  deallog << "Integral = " << integral << ", should be 1." << std::endl;
}


int
main()
{
  std::ofstream logfile ("negative_jac_det/output");
  deallog << std::setprecision(4) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();

  return 0;
}



