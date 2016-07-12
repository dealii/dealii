// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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


// Test that the flat manifold does what it should on a sphere surface.

#include "../tests.h"

#include <fstream>
#include <deal.II/base/logstream.h>


// all include files you need here
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

// Helper function
template <int dim>
void test(unsigned int degree)
{
  deallog << "Testing dim=" << dim <<", degree="
          << degree << std::endl;

  PolarManifold<dim> manifold;
  Triangulation<dim> tria;
  GridGenerator::hyper_shell(tria, Point<dim>(), .4, .6, 6);
  typename Triangulation<dim>::active_cell_iterator cell;

  for (cell = tria.begin_active(); cell != tria.end(); ++cell)
    cell->set_all_manifold_ids(1);

  tria.set_manifold(1, manifold);

  MappingQ<dim> mapping(degree, true);
  FE_Q<dim> fe(degree);
  DoFHandler<dim> dh(tria);
  dh.distribute_dofs(fe);
  QGaussLobatto<dim> quad(degree+1);

  FEValues<dim> fe_v(mapping, fe, quad, update_quadrature_points);

  // char fname[50];
  // sprintf(fname, "out_%d_%d.gpl", dim, degree);
  // std::ofstream ofile(fname);
  std::ostream &ofile = deallog.get_file_stream();

  for (typename DoFHandler<dim>::active_cell_iterator cell = dh.begin_active();
       cell != dh.end(); ++cell)
    {
      fe_v.reinit(cell);
      for (unsigned int q=0; q<quad.size(); ++q)
        {
          ofile << fe_v.get_quadrature_points()[q] << std::endl;
        }
      ofile << std::endl;
    }
  ofile << std::endl;
}

int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.threshold_double(1.e-10);

  for (unsigned int d=1; d<5; ++d)
    {
      test<2>(d);
      test<3>(d);
    }
  return 0;
}

