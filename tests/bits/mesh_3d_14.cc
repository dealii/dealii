//----------------------------  mesh_3d_14.cc  ---------------------------
//    $Id$
//    Version: 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_14.cc  ---------------------------


// similar to mesh_3d_13, but check that quadrature points are
// correct. it is thus also similar to mesh_3d_7, but checks for faces
// and subfaces

#include "../tests.h"
#include "mesh_3d.h"

#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_reordering.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>

#include <fstream>


void check_this (Triangulation<3> &tria)
{
  QTrapez<2> quadrature;
  FE_Q<3> fe(1);
  FEFaceValues<3> fe_face_values1 (fe, quadrature,
                                   update_q_points | update_JxW_values |
                                   update_normal_vectors);
  FESubfaceValues<3> fe_face_values2 (fe, quadrature,
                                      update_q_points | update_JxW_values |
                                      update_normal_vectors);

  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);
  
  DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
  for (; cell!=dof_handler.end(); ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<3>::faces_per_cell;
         ++face_no)
      if (!cell->at_boundary(face_no)
          &&
          cell->neighbor(face_no)->has_children())
        for (unsigned int subface_no=0;
             subface_no<GeometryInfo<3>::subfaces_per_face;
             ++subface_no)
          {
                                             // get an iterator
                                             // pointing to the cell
                                             // behind the present
                                             // subface
            const DoFHandler<3>::cell_iterator
              neighbor = cell->neighbor(face_no);
            const unsigned int neighbor_neighbor
              = cell->neighbor_of_neighbor (face_no);

                                             // see whether face and
                                             // the neighbor's
                                             // counterface share the
                                             // same indexing of
                                             // children. if not so,
                                             // translate child
                                             // indices
            const bool face_orientations_match
              = (neighbor->face_orientation(neighbor_neighbor) ==
                 cell->face_orientation(face_no));
            static const unsigned int subface_translation[4]
              = { 0, 3, 2, 1 };
            const unsigned int neighbor_child_index
              = (GeometryInfo<3>::
                 child_cell_on_face(neighbor_neighbor,
                                    (face_orientations_match ?
                                     subface_no :
                                     subface_translation[subface_no])));
            const DoFHandler<3>::active_cell_iterator neighbor_child
              = neighbor->child(neighbor_child_index);

            fe_face_values1.reinit (neighbor_child, neighbor_neighbor);
            fe_face_values2.reinit (cell, face_no, subface_no);

            for (unsigned int q=0; q<quadrature.n_quadrature_points; ++q)
              {
                Assert ((fe_face_values1.quadrature_point(q)-
                         fe_face_values2.quadrature_point(q)).square()
                        < 1e-20,
                        ExcInternalError());

                Assert (std::fabs(fe_face_values1.JxW(q)-
                                  fe_face_values2.JxW(q)) < 1e-15,
                        ExcInternalError());
                Assert ((fe_face_values1.normal_vector(q) +
                         fe_face_values2.normal_vector(q)).square()
                        < 1e-20,
                        ExcInternalError());
              }            
          }
}



void check (Triangulation<3> &tria)
{
  (++tria.begin_active())->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  
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
  std::ofstream logfile("mesh_3d_14.output");
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

  
  
