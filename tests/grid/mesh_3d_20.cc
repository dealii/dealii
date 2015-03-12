// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2015 by the deal.II authors
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



// check that face_rotation and face_flip flags work by looping over all cells
// and check on all faces that quadrature points match up. based on mesh_3d_7

#include "../tests.h"
#include "../grid/mesh_3d.h"

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
#include <deal.II/numerics/data_out.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

void check_this (Triangulation<3> &tria)
{
  QTrapez<2> quadrature;
  FE_Q<3> fe(1);
  FEFaceValues<3> fe_face_values1 (fe, quadrature,
                                   update_q_points | update_JxW_values);
  FEFaceValues<3> fe_face_values2 (fe, quadrature,
                                   update_q_points | update_JxW_values);

  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  unsigned int global_datum = 0;

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

          for (unsigned int q=0; q<quadrature.size(); ++q)
            {
              // in order to reduce
              // output file size,
              // only write every
              // 289th datum. if the
              // values differ
              // anyway, then the
              // assertion below will
              // catch this, and if
              // we compute _all_
              // values wrongly, then
              // outputting some will
              // be ok, I guess
              if (global_datum++ % 17*17 == 0)
                deallog << "Cell " << cell << ", face " << f
                        << std::endl
                        << "  " << fe_face_values1.quadrature_point(q)
                        << ", " << fe_face_values1.JxW(q)
                        << std::endl;

              Assert ((fe_face_values1.quadrature_point(q)-
                       fe_face_values2.quadrature_point(q)).norm_square()
                      < 1e-20,
                      ExcInternalError());

              Assert (std::fabs(fe_face_values1.JxW(q)-
                                fe_face_values2.JxW(q)) < 1e-15,
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
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

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

}



