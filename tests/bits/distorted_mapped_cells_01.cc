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


// Test that an assertion is thrown when MappingQEulerian produces a cell with
// negative volume

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/grid/tria.h>


void test()
{
  const int dim = 2;
  // create a dummy triangulation with no extension and set the geometry
  // through MappingQEulerian
  Triangulation<dim> tria;
  std::vector<Point<dim> > points (4);

  std::vector<CellData<dim> > cells (1);
  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 2;
  cells[0].vertices[3] = 3;
  cells[0].material_id = 0;

  tria.create_triangulation(points, cells, SubCellData());

  FE_Q<dim> fe(1);
  FESystem<dim> fe_sys(fe, dim);
  DoFHandler<dim> dof_h(tria);
  dof_h.distribute_dofs(fe_sys);
  Vector<double> displacements(dof_h.n_dofs());
  displacements(2) = -1.;
  displacements(5) = 1.;
  displacements(6) = -1.;
  displacements(7) = 1.;

  // this gives a Cartesian cell but in non-standard orientation (x-coordinate
  // is gone through backwards)
  MappingQEulerian<dim> mapping(1,displacements, dof_h);
  QGauss<dim> quad(1);
  FEValues<dim> fe_val (mapping, fe, quad, update_JxW_values);
  double integral = 0.;
  /*typename*/ Triangulation<dim>::active_cell_iterator
  cell = tria.begin_active(), endc = tria.end();
  for ( ; cell != endc; ++cell)
    {
      try
        {
          fe_val.reinit (cell);
          for (unsigned int q=0; q<quad.size(); ++q)
            integral += fe_val.JxW(q);
        }
      catch (ExceptionBase &e)
        {
          deallog << e.get_exc_name() << std::endl;
        }
    }
  deallog << "Integral = " << integral << std::endl;
}


int
main()
{
  deal_II_exceptions::disable_abort_on_exception();

  std::ofstream logfile ("output");
  deallog << std::setprecision(4) << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test();

  return 0;
}



