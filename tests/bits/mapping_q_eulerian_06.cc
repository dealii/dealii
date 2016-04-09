// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2015 by the deal.II authors
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



// computes points in real space for 1D Eulerian mapping where middle points
// are moved

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files you need here

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/mapping_q_eulerian.h>
#include <deal.II/base/quadrature_lib.h>

#include <fstream>
#include <string>

std::ofstream logfile("output");

void test(unsigned int degree)
{
  const unsigned int dim = 1;
  const unsigned int spacedim = 1;
  Triangulation<dim, spacedim> tria;
  GridGenerator::hyper_cube(tria, 0, 1);
  FE_Q<dim, spacedim> fe(QIterated<1>(QTrapez<1>(),degree));

  DoFHandler<dim, spacedim> shift_dh(tria);

  shift_dh.distribute_dofs(fe);

  // Shift just the interior points but not the boundary points
  Vector<double> shift(shift_dh.n_dofs());
  for (unsigned int i=2; i<=degree; ++i)
    shift(i) = 0.1;

  QGauss<dim> quad(degree+1);
  MappingQEulerian<dim,Vector<double>,spacedim> mapping(degree,shift, shift_dh);

  Triangulation<dim,spacedim>::active_cell_iterator cell=tria.begin_active(),
                                                    endc=tria.end() ;
  Point<spacedim> real;
  Point<dim> unit;
  double eps = 1e-10;
  for (; cell!=endc; ++cell)
    {
      deallog<<cell<< std::endl;
      for (unsigned int q=0; q<quad.size(); ++q)
        {
          real = mapping.transform_unit_to_real_cell(cell, quad.point(q));
          unit = mapping.transform_real_to_unit_cell(cell, real);
          deallog<<quad.point(q)<< " -> " << real << std::endl;
          if ( (unit-quad.point(q)).norm()>eps)
            deallog<<quad.point(q)<< " != " << unit << std::endl;
        }
    }

}

int main ()
{
  deallog.attach(logfile);

  test(1);
  test(2);
  test(3);
  return 0;
}
