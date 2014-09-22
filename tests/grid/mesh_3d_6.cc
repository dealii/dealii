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



// check that face orientation flags work by looping over all cells
// and check on all faces that if we look from both sides that normal
// vectors point in opposite directions

#include "../tests.h"
#include "mesh_3d.h"

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



