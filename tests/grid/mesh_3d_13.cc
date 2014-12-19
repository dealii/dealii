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



// this tests the assertion that is actually triggered in mesh_3d_12,
// isolated from the rest of the code in which it sits in that test
//
// actually, the computation of the orientation flag was wrong as we
// were also considering the orientation flag of the present cell,
// which we shouldn't have

#include "../tests.h"
#include "mesh_3d.h"

#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>

#include <fstream>


void check_this (Triangulation<3> &tria)
{
  Triangulation<3>::active_cell_iterator cell = tria.begin_active();
  for (; cell!=tria.end(); ++cell)
    for (unsigned int face_no=0; face_no<GeometryInfo<3>::faces_per_cell;
         ++face_no)
      if (!cell->at_boundary(face_no)
          &&
          cell->neighbor(face_no)->has_children())
        for (unsigned int subface_no=0;
             subface_no<GeometryInfo<3>::max_children_per_face;
             ++subface_no)
          {
            // get an iterator
            // pointing to the cell
            // behind the present
            // subface
            const Triangulation<3>::cell_iterator
            neighbor = cell->neighbor(face_no);
            const unsigned int neighbor_neighbor
              = cell->neighbor_of_neighbor (face_no);
            const bool orientation_flag
              = (neighbor->face_orientation(neighbor_neighbor) == true);
            static const unsigned int subface_translation[4]
              = { 0, 2, 1, 3 };
            const unsigned int neighbor_child_index
              = (GeometryInfo<3>::
                 child_cell_on_face(RefinementCase<3>::isotropic_refinement,neighbor_neighbor,
                                    (orientation_flag ?
                                     subface_no :
                                     subface_translation[subface_no])));
            const Triangulation<3>::active_cell_iterator neighbor_child
              = neighbor->child(neighbor_child_index);

            Assert (neighbor_child->face(neighbor_neighbor) ==
                    cell->face(face_no)->child(subface_no),
                    ExcInternalError());
            Assert (!neighbor->child(neighbor_child_index)->has_children(),
                    ExcInternalError());
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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

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



