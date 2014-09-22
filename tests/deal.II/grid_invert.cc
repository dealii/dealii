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
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/base/logstream.h>

#include <fstream>
#include <string>

std::ofstream logfile("output");

template <int dim>
void test (bool second_case=false)
{
  std::vector<Point<dim> > vertices(GeometryInfo<dim>::vertices_per_cell);
  vertices[1](1)=1;
  vertices[2](0)=1;
  vertices[2](1)=1;
  vertices[3](0)=1;
  if (dim==3)
    {
      for (unsigned int i=4; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i](2)=-1;
      vertices[5](1)=1;
      vertices[6](0)=1;
      vertices[6](1)=1;
      vertices[7](0)=1;
    }
  std::vector<CellData<dim> > cells(1);
  for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    cells[0].vertices[i]=i;

  if (dim==3 && second_case)
    {
      std::swap(cells[0].vertices[1],cells[0].vertices[3]);
      std::swap(cells[0].vertices[5],cells[0].vertices[7]);
      for (unsigned int i=4; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        vertices[i](2)=1;
    }

  SubCellData subcelldata;
  GridReordering<dim>::invert_all_cells_of_negative_grid(vertices,cells);

  Triangulation<dim> tria;
  tria.create_triangulation_compatibility(vertices,cells,subcelldata);

  logfile<<"---------------------------------------------"
         <<std::endl<<"dim="<<dim
         <<(second_case ? ", second case" : ", first case")
         <<std::endl<<std::endl;

  GridOut grid_out;
  grid_out.set_flags (GridOutFlags::Ucd(true));
  grid_out.write_ucd (tria, logfile);
}

int main ()
{
  test<2> ();
  test<3> (false);
  test<3> (true);
}

