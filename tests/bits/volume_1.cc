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



// like a union of the normals_* and q_point_sum_* tests: integrating
// \vec x times the normal over the surface of any body should yield
// the volume of this body times the space dimension

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>



template <int dim>
void check (const Triangulation<dim> &tria)
{
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss<dim-1> q_face(3);

  FEFaceValues<dim>    fe_face_values (fe, q_face,
                                       update_normal_vectors |
                                       update_q_points |
                                       update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (fe, q_face,
                                          update_normal_vectors |
                                          update_q_points |
                                          update_JxW_values);

  double v1=0, v2=0;
  for (typename DoFHandler<dim>::active_cell_iterator
       cell = dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {

      // first integrate over faces
      // and make sure that the
      // result of the integration is
      // close to zero
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          {
            fe_face_values.reinit (cell, f);
            for (unsigned int q=0; q<q_face.size(); ++q)
              v1 += (fe_face_values.normal_vector(q) *
                     fe_face_values.quadrature_point(q)) *
                    fe_face_values.JxW(q);
          }

      // now same for subface
      // integration
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          for (unsigned int sf=0; sf<GeometryInfo<dim>::max_children_per_face; ++sf)
            {
              fe_subface_values.reinit (cell, f, sf);
              for (unsigned int q=0; q<q_face.size(); ++q)
                v2 += (fe_subface_values.normal_vector(q) *
                       fe_subface_values.quadrature_point(q)) *
                      fe_subface_values.JxW(q);
            }
    }

  Assert (std::fabs(v1-v2) < 1e-12, ExcInternalError());
  deallog << " face integration: "
          << v1 / dim
          << std::endl;
  deallog << " subface integration: "
          << v2 / dim
          << std::endl;
}


int main ()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid);
    check (coarse_grid);
  }


  {
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    check (coarse_grid);
  }

  {
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_ball (coarse_grid);
    check (coarse_grid);
  }

}



