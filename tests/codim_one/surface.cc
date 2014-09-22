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



// calculates the measure of the surface of a hypersphere

#include "../tests.h"
#include <fstream>
#include <deal.II/base/logstream.h>

// all include files needed for the program

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>

#include <fstream>
#include <string>

std::ofstream logfile("output");

// Computes the area and the outer normals of circles and spheres
// which are more and more refined, and prints the error on the
// output.
template <int dim, int spacedim>
void test(std::string filename)
{

  Triangulation<dim, spacedim> triangulation;
  GridIn<dim, spacedim> gi;

  gi.attach_triangulation (triangulation);
  std::ifstream in (filename.c_str());
  gi.read_ucd (in);

  const QGauss<dim> quadrature(2);
  const FE_Q<dim,spacedim> dummy_fe (1);
  DoFHandler<dim,spacedim> dof_handler (triangulation);

  FEValues<dim,spacedim> fe_values (dummy_fe, quadrature,
                                    update_JxW_values |
                                    update_cell_normal_vectors |
                                    update_quadrature_points);

  dof_handler.distribute_dofs (dummy_fe);

  double area = 0;
  double normals = 0;

  typename DoFHandler<dim,spacedim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();

  std::vector<Point<spacedim> > expectedcellnormals(fe_values.n_quadrature_points);

  for (; cell!=endc; ++cell)
    {
      fe_values.reinit (cell);
      const std::vector<Point<spacedim> > &cellnormals = fe_values.get_cell_normal_vectors();
      const std::vector<Point<spacedim> > &quad_points = fe_values.get_quadrature_points();

      for (unsigned int i=0; i<fe_values.n_quadrature_points; ++i)
        {
          expectedcellnormals[i] = quad_points[i]/quad_points[i].norm();
          area += fe_values.JxW (i);
          normals += (expectedcellnormals[i]-cellnormals[i]).norm();
        }
    };

  deallog<<"Approximate measure of hyper sphere = "<<area<<std::endl;
  deallog<<"Error = "<<std::fabs(dim*2*numbers::PI-area)<<std::endl;
  deallog << "Average error in norms: "
          << ( normals/dof_handler.get_tria().n_active_cells()
               /fe_values.n_quadrature_points)
          << std::endl;

}



int main ()
{
  deallog.attach(logfile);
  deallog.depth_console(0);

  deallog<<"Test <1,2>"<<std::endl;
  test<1,2>(SOURCE_DIR "/grids/circle_1.inp");
  test<1,2>(SOURCE_DIR "/grids/circle_2.inp");
  test<1,2>(SOURCE_DIR "/grids/circle_3.inp");
  test<1,2>(SOURCE_DIR "/grids/circle_4.inp");

  deallog<<"Test <2,3>"<<std::endl;
  test<2,3>(SOURCE_DIR "/grids/sphere_1.inp");
  test<2,3>(SOURCE_DIR "/grids/sphere_2.inp");
  test<2,3>(SOURCE_DIR "/grids/sphere_3.inp");
  test<2,3>(SOURCE_DIR "/grids/sphere_4.inp");

  return 0;
}
