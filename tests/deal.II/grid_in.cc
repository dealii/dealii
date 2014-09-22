// ---------------------------------------------------------------------
//
// Copyright (C) 2002 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_boundary.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <iomanip>
#include <string>

std::ofstream logfile("output");


template <int dim>
void test1 ()
{
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (SOURCE_DIR "/grid_in/2d.inp");
  gi.read_ucd (in);

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);
}


template <int dim>
void test2 ()
{
  // read a much larger grid (30k
  // cells). with the old grid
  // reordering scheme, this took >90
  // minutes (exact timing not
  // available, program was killed
  // before), with the new one it
  // takes less than 8 seconds
  //
  // note that the input file is no good: it
  // contains two pairs of two cells, where
  // the two members of each pair share 3
  // vertices (in 2d) -- this can of course
  // not work properly. it makes the
  // grid_in_02 testcase fail when using this
  // input file, but grid_in_02/2d.xda is a
  // corrected input file.
  Triangulation<dim> tria (Triangulation<dim>::none, true);
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  std::ifstream in (SOURCE_DIR "/grid_in/2d.xda");
  try
    {
      gi.read_xda (in);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      // ignore the exception that we
      // get because the mesh has
      // distorted cells
      deallog << dcv.distorted_cells.size() << " cells are distorted."
              << std::endl;
    }


  int hash = 0;
  int index = 0;
  for (typename Triangulation<dim>::active_cell_iterator c=tria.begin_active();
       c!=tria.end(); ++c, ++index)
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      hash += (index * i * c->vertex_index(i)) % (tria.n_active_cells()+1);
  deallog << hash << std::endl;
}


template <int dim>
void test3 ()
{
  Triangulation<dim> tria;
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  gi.read (SOURCE_DIR "/grid_in/2d.nc");

  GridOut grid_out;
  std::ofstream gnufile("grid_in_2d.gnuplot");
  grid_out.write_gnuplot (tria, gnufile);
}


template<int dim>
void check_file (const std::string name,
                 typename GridIn<dim>::Format format)
{
  Triangulation<dim> tria (Triangulation<dim>::none, true);
  GridIn<dim> gi;
  gi.attach_triangulation (tria);
  try
    {
      gi.read(name, format);
    }
  catch (typename Triangulation<dim>::DistortedCellList &dcv)
    {
      // ignore the exception
      deallog << dcv.distorted_cells.size() << " cells are distorted."
              << std::endl;
    }

  deallog << '\t' << tria.n_vertices()
          << '\t' << tria.n_cells()
          << std::endl;
}

void filename_resolution()
{
  check_file<2> (std::string(SOURCE_DIR "/grid_in/2d"), GridIn<2>::ucd);
  check_file<2> (std::string(SOURCE_DIR "/grid_in/2d"), GridIn<2>::xda);
}


int main ()
{
  deallog << std::setprecision (2);
  logfile << std::setprecision (2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test1<2> ();
  test2<2> ();
  // test3 needs NetCDF
//    test3<2> ();

  filename_resolution();
}

