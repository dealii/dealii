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



// like normals_1, but with the type of cells we have in the mesh_3d_*
// tests: integrating the normals over the surface of any cell
// (distorted or not) should yield a zero vector. (To prove this, use
// the divergence theorem.)

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




void check (Triangulation<3> &tria)
{
  FE_Q<3> fe(1);
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss<2> q_face(3);

  FEFaceValues<3>    fe_face_values (fe, q_face,
                                     update_normal_vectors | update_JxW_values);
  FESubfaceValues<3> fe_subface_values (fe, q_face,
                                        update_normal_vectors | update_JxW_values);

  for (DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      Tensor<1,3> n1, n2;

      // first integrate over faces
      // and make sure that the
      // result of the integration is
      // close to zero
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
        {
          fe_face_values.reinit (cell, f);
          for (unsigned int q=0; q<q_face.size(); ++q)
            n1 += fe_face_values.normal_vector(q) *
                  fe_face_values.JxW(q);
        }
      Assert (n1*n1 < 1e-24, ExcInternalError());
      deallog << cell << " face integration is ok: "
              << std::sqrt (n1*n1)
              << std::endl;

      // now same for subface
      // integration
      for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f)
        for (unsigned int sf=0; sf<GeometryInfo<3>::max_children_per_face; ++sf)
          {
            fe_subface_values.reinit (cell, f, sf);
            for (unsigned int q=0; q<q_face.size(); ++q)
              n2 += fe_subface_values.normal_vector(q) *
                    fe_subface_values.JxW(q);
          }
      Assert (n2*n2 < 1e-24, ExcInternalError());
      deallog << cell << " subface integration is ok: "
              << std::sqrt (n2*n2)
              << std::endl;
    }
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



