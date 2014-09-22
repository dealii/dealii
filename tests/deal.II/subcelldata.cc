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



#include "../tests.h"
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/tria_accessor.h>

#include <fstream>


static unsigned subcells[6][4] = {{0, 1, 2, 3},
  {4, 5, 6, 7},
  {0, 1, 5, 4},
  {1, 5, 6, 2},
  {3, 2, 6, 7},
  {0, 4, 7, 3}
};



template <int dim>
void test()
{
  Assert(dim==2 || dim==3, ExcNotImplemented());

  std::vector<Point<dim> > vertices (GeometryInfo<dim>::vertices_per_cell);
  vertices[0](0)=0;
  vertices[0](1)=0;
  vertices[1](0)=2;
  vertices[1](1)=1;
  vertices[2](0)=3;
  vertices[2](1)=3;
  vertices[3](0)=0;
  vertices[3](1)=1;
  if (dim==3)
    {
      // for the new numbering
//       for (unsigned int i=0; i<4; ++i)
//  {
//    vertices[i+4]=vertices[i];
//    vertices[i+4](2)=1;
//  }
      // for the old numbering
      for (unsigned int i=0; i<4; ++i)
        {
          std::swap(vertices[i](1),vertices[i](2));
          vertices[i+4]=vertices[i];
          vertices[i+4](1)=1;
        }

    }

  std::vector<CellData<dim> > cells (1);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    cells[0].vertices[i] = i;
  cells[0].material_id = 0;

  SubCellData subcelldata;
  if (dim==2)
    {
      subcelldata.boundary_lines.resize(GeometryInfo<dim>::faces_per_cell);
      for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
        {
          subcelldata.boundary_lines[i].vertices[0]=i;
          subcelldata.boundary_lines[i].vertices[1]=(i+1)%4;
          subcelldata.boundary_lines[i].material_id=10*i+1;
        }
    }
  else if (dim==3)
    {
      subcelldata.boundary_quads.resize(GeometryInfo<dim>::faces_per_cell);
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        {
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_face; ++i)
            subcelldata.boundary_quads[f].vertices[i]=subcells[f][i];
          subcelldata.boundary_quads[f].material_id=10*f+1;
        }
    }

  Triangulation<dim> tria;
  tria.create_triangulation_compatibility (vertices, cells, subcelldata);

  GridOutFlags::Ucd ucd_flags(true,true);
  GridOut grid_out;
  grid_out.set_flags(ucd_flags);
  grid_out.write_ucd(tria, deallog.get_file_stream());

//   std::ofstream gnuplot_file("subcelldata.gnuplot");
//   grid_out.write_gnuplot(tria, gnuplot_file);
//   std::ofstream ucd_file("subcelldata.inp");
//   grid_out.write_ucd(tria, ucd_file);
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  test<2>();
  test<3>();

  return 0;
}
