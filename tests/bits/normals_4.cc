//----------------------------  normals_4.cc  ---------------------------
//    normals_4.cc,v 1.1 2003/10/19 22:29:38 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  normals_4.cc  ---------------------------


// integrating the normals over the surface of any cell (distorted or
// not) should yield a zero vector. (To prove this, use the divergence
// theorem.) check that this is indeed so for the hypercube in 2d and 3d
//
// in contrast to normals_1, do this with a cartesian mapping
// associated with the geometry

#include "../tests.h"
#include <base/logstream.h>
#include <base/quadrature_lib.h>
#include <grid/tria.h>
#include <grid/tria_boundary_lib.h>
#include <grid/tria_accessor.h>
#include <grid/tria_iterator.h>
#include <grid/grid_generator.h>
#include <dofs/dof_handler.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_cartesian.h>

#include <fstream>



template <int dim>
void check (const Triangulation<dim> &tria)
{
  MappingCartesian<dim> mapping;
  FE_Q<dim> fe(1);
  DoFHandler<dim> dof_handler (tria);
  dof_handler.distribute_dofs (fe);

  QGauss3<dim-1> q_face;
  
  FEFaceValues<dim>    fe_face_values (mapping, fe, q_face,
                                       update_normal_vectors | update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (mapping, fe, q_face,
                                          update_normal_vectors | update_JxW_values);

  for (typename DoFHandler<dim>::active_cell_iterator
         cell = dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
      Tensor<1,dim> n1, n2;

                                       // first integrate over faces
                                       // and make sure that the
                                       // result of the integration is
                                       // close to zero
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
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
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        for (unsigned int sf=0; sf<GeometryInfo<dim>::subfaces_per_face; ++sf)
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
  std::ofstream logfile("normals_4.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

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
}

  
  
