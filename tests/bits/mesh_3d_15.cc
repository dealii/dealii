//----------------------------  mesh_3d_15.cc  ---------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mesh_3d_15.cc  ---------------------------


// make sure TriaCellAccessor::neighbor_child_on_subface does what it
// is supposed to do. check it for dof accessors, since they simply
// call the tria accessors, and this way we catch both cases at the
// same time

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

#include <fstream>


void check_this (Triangulation<3> &tria)
{
  FE_Q<3> fe(1);
  
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

	    const DoFHandler<3>::face_iterator subface=
	      cell->face(face_no)->child(subface_no);
	    unsigned int i=0;
	    for (; i<GeometryInfo<3>::children_per_cell; ++i)
	      if (neighbor->child(i)->face(neighbor_neighbor)==subface)
		break;

            Assert (neighbor->child(i)==
                    cell->neighbor_child_on_subface (face_no,
                                                     subface_no),
                    ExcInternalError());
          }
}



void check (Triangulation<3> &tria)
{
  (++tria.begin_active())->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();
  
  deallog << "Initial check" << std::endl;
  check_this (tria);
//TODO:[WB] Is there a reason to do this three times?
// Changed to two. Guido
  for (unsigned int r=0; r<2; ++r)
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
  std::ofstream logfile("mesh_3d_15.output");
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

  
  
