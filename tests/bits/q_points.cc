// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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



// make sure that if we loop over all quadrature points on a face and
// over the same quadrature points on the same face of the neighboring
// cell, that we then visit them in the same order. this test is
// basically the same as the mesh_3d_7 test, except that we leave out
// the complication of mis-oriented faces

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>


void create_two_cubes (Triangulation<3> &coarse_grid)
{
  const Point<3> points[6] = { Point<3>(0,0,0),
                               Point<3>(1,0,0),
                               Point<3>(1,1,0),
                               Point<3>(0,1,0),
                               Point<3>(2,0,0),
                               Point<3>(2,1,0)
                             };
  std::vector<Point<3> > vertices;
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i]);
  for (unsigned int i=0; i<6; ++i)
    vertices.push_back (points[i] + Point<3>(0,0,-1));

  std::vector<CellData<3> > cells;

  const unsigned int connectivity[2][4]
  = { { 0,1,2,3 },
    { 1,4,5,2 }
  };
  for (unsigned int i=0; i<2; ++i)
    {
      CellData<3> cell;
      for (unsigned int j=0; j<4; ++j)
        {
          cell.vertices[j]   = connectivity[i][j];
          cell.vertices[j+4] = connectivity[i][j]+6;
        }
      cells.push_back (cell);
    }

  // finally generate a triangulation
  // out of this
  GridReordering<3>::reorder_cells (cells);
  coarse_grid.create_triangulation_compatibility (vertices, cells, SubCellData());
}


void check (Triangulation<3> &tria)
{
  QGauss<2> quadrature(3);
  FE_Q<3> fe(1);
  FEFaceValues<3> fe_face_values1 (fe, quadrature,
                                   update_q_points | update_JxW_values);
  FEFaceValues<3> fe_face_values2 (fe, quadrature,
                                   update_q_points | update_JxW_values);

  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  for (DoFHandler<3>::cell_iterator cell=dof_handler.begin();
       cell != dof_handler.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      if (!cell->at_boundary(f))
        {
          const unsigned int nn
            = cell->neighbor_of_neighbor (f);

          // this test is about
          // properly oriented
          // faces. mesh_3d_7 does it
          // for mis-oriented faces
          Assert (cell->face_orientation(f) == true,
                  ExcInternalError());
          Assert (cell->neighbor(f)->face_orientation(nn) == true,
                  ExcInternalError());

          fe_face_values1.reinit (cell, f);
          fe_face_values2.reinit (cell->neighbor(f), nn);

          deallog << "Cell " << cell << ", face " << f
                  << std::endl;

          for (unsigned int q=0; q<quadrature.size(); ++q)
            {
              deallog << "  " << fe_face_values1.quadrature_point(q)
                      << std::endl;
              Assert ((fe_face_values1.quadrature_point(q)-
                       fe_face_values2.quadrature_point(q)).square()
                      < 1e-20,
                      ExcInternalError());
            }
        }
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  Triangulation<3> coarse_grid;
  create_two_cubes (coarse_grid);
  check (coarse_grid);

  coarse_grid.refine_global (1);
  check (coarse_grid);
}



