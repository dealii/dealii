//----------------------------  mesh_3d_7.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_7.cc  ---------------------------


// check that face orientation flags work by looping over all cells
// and check on all faces that quadrature points match up

#include "mesh_3d.h"

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <fstream>


void check_this (Triangulation<3> &tria)
{
  QGauss3<2> quadrature;
  FE_Q<3> fe(1);
  FEFaceValues<3> fe_face_values1 (fe, quadrature,
                                   update_q_points | update_JxW_values);
  FEFaceValues<3> fe_face_values2 (fe, quadrature,
                                   update_q_points | update_JxW_values);
  
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

                                   // look at all faces, not only
                                   // active ones
  for (DoFHandler<3>::cell_iterator cell=dof_handler.begin();
       cell != dof_handler.end(); ++cell)
    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
      if (!cell->at_boundary(f))
        {
          const unsigned int nn
            = cell->neighbor_of_neighbor (f);
          fe_face_values1.reinit (cell, f);
          fe_face_values2.reinit (cell->neighbor(f), nn);

          deallog << "Cell " << cell << ", face " << f
                  << std::endl;

          for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
            {
              deallog << "  " << fe_face_values1.quadrature_point(q)
                      << ", " << fe_face_values1.JxW(q)
                      << std::endl;
              
              Assert (fe_face_values1.quadrature_point(q) ==
                      fe_face_values2.quadrature_point(q),
                      ExcInternalError());

              Assert (fe_face_values1.JxW(q) ==
                      fe_face_values2.JxW(q),
                      ExcInternalError());
            }
        }
}


void check (Triangulation<3> &tria)
{
  deallog << "Initial check" << std::endl;
  check_this (tria);
  
  for (unsigned int r=0; r<3; ++r)
    {
      tria.refine_global (1);
      deallog << "Check " << r << std::endl;
      check_this (tria);
    }

  coarsen_global (tria);
  deallog << "Check " << 1 << std::endl;
  check_this (tria);
  
  tria.refine_global (1);
  deallog << "Check " << 2 << std::endl;
  check_this (tria);
}


int main () 
{
  std::ofstream logfile("mesh_3d_7.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {  
    Triangulation<3> coarse_grid;
    create_two_cubes (coarse_grid);
    check (coarse_grid);
  }
  
  {  
    Triangulation<3> coarse_grid;
    create_L_shape (coarse_grid);
    check (coarse_grid);
  }
  
  {  
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    check (coarse_grid);
  }
  
}

  
  
