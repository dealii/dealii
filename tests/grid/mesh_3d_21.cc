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



// check, that we find our way back from fine cells over coarser neighbors to
// the fine cell agaion and vice versa

#include "../tests.h"
#include "mesh_3d.h"

#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

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
      if (!cell->at_boundary(face_no))
        {
          if (cell->neighbor(face_no)->has_children())
            // we are coarser than the neighbor
            for (unsigned int subface_no=0;
                 subface_no<cell->face(face_no)->n_children();
                 ++subface_no)
              {
                // get an iterator
                // pointing to the cell
                // behind the present
                // subface
                const unsigned int neighbor_neighbor
                  = cell->neighbor_of_neighbor (face_no);

                const DoFHandler<3>::active_cell_iterator neighbor_child
                  = cell->neighbor_child_on_subface(face_no,subface_no);

                // make sure, that we find the
                // way back
                const unsigned int our_face_no=neighbor_child->neighbor_of_coarser_neighbor(neighbor_neighbor).first;
                const unsigned int our_subface_no=neighbor_child->neighbor_of_coarser_neighbor(neighbor_neighbor).second;

                Assert (our_face_no==face_no, ExcInternalError());
                Assert (our_subface_no==subface_no, ExcInternalError());
                deallog << "from coarse to fine and back: OK" <<std::endl;
              }
          else if (cell->neighbor(face_no)->level()<cell->level())
            // the neighbor is coarser
            {
              const unsigned int neighbor_face_no=cell->neighbor_of_coarser_neighbor(face_no).first;
              const unsigned int neighbor_subface_no=cell->neighbor_of_coarser_neighbor(face_no).second;

              // try to find the way back to our cell
              const DoFHandler<3>::active_cell_iterator our_cell=cell->neighbor(face_no)->neighbor_child_on_subface(neighbor_face_no,
                                                                 neighbor_subface_no);
              Assert (our_cell==cell, ExcInternalError());
              deallog << "from fine to coarse and back: OK" <<std::endl;
            }
        }


}



void check (Triangulation<3> &tria)
{
  (++tria.begin_active())->set_refine_flag ();
  tria.execute_coarsening_and_refinement ();

  deallog << "Initial check" << std::endl;
  check_this (tria);

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
    create_two_cubes_rotation (coarse_grid,1);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes_rotation (coarse_grid,2);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    create_two_cubes_rotation (coarse_grid,3);
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



