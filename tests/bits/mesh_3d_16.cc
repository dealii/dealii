//----------------------------  mesh_3d_16.cc  ---------------------------
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
//----------------------------  mesh_3d_16.cc  ---------------------------


// like normals_1, but with the type of cells we have in the mesh_3d_*
// tests: integrating the normals over the surface of any cell
// (distorted or not) should yield a zero vector. (To prove this, use
// the divergence theorem.)

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




void check (Triangulation<3> &tria)
{
  FE_Q<3> fe(1);
  DoFHandler<3> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss3<2> q_face;
  
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
          for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
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
        for (unsigned int sf=0; sf<GeometryInfo<3>::subfaces_per_face; ++sf)
        {
          fe_subface_values.reinit (cell, f, sf);
          for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
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
  std::ofstream logfile("mesh_3d_16.output");
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

  
  
