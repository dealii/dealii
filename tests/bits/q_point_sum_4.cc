//----------------------------  q_point_sum_4.cc  ---------------------------
//    q_point_sum_4.cc,v 1.1 2003/10/19 22:29:39 wolf Exp
//    Version: 
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  q_point_sum_4.cc  ---------------------------


// integrating \vec x over the surface of the [-1,1] hypercube and
// hyperball in 2d and 3d should yield zero
//
// same as q_point_sum_1, but with a Cartesian mapping


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
#include <fe/mapping_cartesian.h>
#include <fe/fe_values.h>

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
                                       update_q_points | update_JxW_values);
  FESubfaceValues<dim> fe_subface_values (mapping, fe, q_face,
                                          update_q_points | update_JxW_values);

  Point<dim> n1, n2;
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
            for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
              n1 += fe_face_values.quadrature_point(q) *
                    fe_face_values.JxW(q);
          }
      
                                       // now same for subface
                                       // integration
      for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
        if (cell->at_boundary(f))
          for (unsigned int sf=0; sf<GeometryInfo<dim>::subfaces_per_face; ++sf)
            {
              fe_subface_values.reinit (cell, f, sf);
              for (unsigned int q=0; q<q_face.n_quadrature_points; ++q)
                n2 += fe_subface_values.quadrature_point(q) *
                      fe_subface_values.JxW(q);
            }
    }
  
  Assert (n1*n1 < 1e-24, ExcInternalError());
  deallog << " face integration is ok: "
          << std::sqrt (n1*n1)
          << std::endl;
  Assert (n2*n2 < 1e-24, ExcInternalError());
  deallog << " subface integration is ok: "
          << std::sqrt (n2*n2)
          << std::endl;
}


int main () 
{
  std::ofstream logfile("q_point_sum_4.output");
  deallog.attach(logfile);
  deallog.depth_console(0);

  {  
    Triangulation<2> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid, -1, 1);
    check (coarse_grid);
  }
  
  {  
    Triangulation<3> coarse_grid;
    GridGenerator::hyper_cube (coarse_grid, -1, 1);
    check (coarse_grid);
  }
}

  
  
