//----------------------------  mesh_3d_6.cc  ---------------------------
//    mesh_3d_6.cc,v 1.2 2003/09/29 15:16:48 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_6.cc  ---------------------------


// check that face orientation flags work by looping over all cells
// and check on all faces that if we look from both sides that normal
// vectors point in opposite directions

#include "../tests.h"
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
  QMidpoint<2> q;
  FE_Q<3> fe(1);
  FEFaceValues<3> fe_face_values1 (fe, q, update_normal_vectors);
  FEFaceValues<3> fe_face_values2 (fe, q, update_normal_vectors);

  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  unsigned int global_face = 0;
  
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

                                           // in order to reduce
                                           // output file size, only
                                           // write every seventeenth
                                           // normal vector. if the
                                           // normals differ anyway,
                                           // then the assertion below
                                           // will catch this, and if
                                           // we compute _all_ normals
                                           // wrongly, then outputting
                                           // some will be ok, I guess
          if (global_face++ % 17 == 0)
            deallog << "Cell " << cell << ", face " << f
                    << " n=" << fe_face_values1.normal_vector(0)
                    << std::endl;

                                           // normal vectors should be
                                           // in opposite directions,
                                           // so their sum should be
                                           // close to zero
          Assert ((fe_face_values1.normal_vector(0) +
                  fe_face_values2.normal_vector(0)).square()
                  < 1e-20,
                  ExcInternalError());
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
  std::ofstream logfile("mesh_3d_6.output");
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

  
  
