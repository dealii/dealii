// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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


// like fe_face_values_1d but testing that one can use MappingCartesian

#include "../tests.h"
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe.h>
#include <vector>
#include <fstream>
#include <iomanip>
#include <string>
#include <sstream>

#define PRECISION 4



template<int dim>
inline void
plot_faces(Mapping<dim> &mapping,
           FiniteElement<dim> &fe,
           typename DoFHandler<dim>::cell_iterator &cell)
{
  // create a QGauss<0>(4), which should
  // still only have 1 quadrature point
  QGauss<dim-1> q(4);
  Assert (q.size() == 1, ExcInternalError());

  FEFaceValues<dim> fe_values(mapping, fe, q,
                              UpdateFlags(update_q_points
                                          | update_JxW_values
                                          | update_values
                                          | update_gradients
                                          | update_hessians
                                          | update_normal_vectors));

  for (unsigned int face_nr=0;
       face_nr < GeometryInfo<dim>::faces_per_cell;
       ++ face_nr)
    {
      deallog << "Face=" << face_nr << std::endl;
      fe_values.reinit(cell, face_nr);

      // print some data on the location
      deallog << "x=" << fe_values.quadrature_point(0) << std::endl;
      deallog << "n=" << fe_values.normal_vector(0) << std::endl;
      deallog << "JxW=" << fe_values.JxW(0) << std::endl;

      // now print some data on the shape
      // functions
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        deallog << "shape_function " << i << ":" << std::endl
                << "  phi=" << fe_values.shape_value(i,0) << std::endl
                << "  grad phi=" << fe_values.shape_grad(i,0) << std::endl
                << "  grad^2 phi=" << fe_values.shape_hessian(i,0) << std::endl;

      deallog << std::endl;
    }
}



template<int dim>
inline void
plot_subfaces(Mapping<dim> &mapping,
              FiniteElement<dim> &fe,
              typename DoFHandler<dim>::cell_iterator &cell)
{

  // create a QGauss<0>(4), which should
  // still only have 1 quadrature point
  QGauss<dim-1> q(4);
  Assert (q.size() == 1, ExcInternalError());

  FESubfaceValues<dim> fe_values(mapping, fe, q,
                                 UpdateFlags(update_q_points
                                             | update_JxW_values
                                             | update_values
                                             | update_gradients
                                             | update_hessians
                                             | update_normal_vectors));
  for (unsigned int face_nr=0;
       face_nr < GeometryInfo<dim>::faces_per_cell;
       ++ face_nr)
    for (unsigned int sub_nr=0;
         sub_nr < GeometryInfo<dim>::max_children_per_face;
         ++ sub_nr)
      {
        deallog << "Face=" << face_nr
                << ", subface=" << sub_nr
                << std::endl;

        fe_values.reinit(cell, face_nr, sub_nr);

        // print some data on the location
        deallog << "x=" << fe_values.quadrature_point(0) << std::endl;
        deallog << "n=" << fe_values.normal_vector(0) << std::endl;
        deallog << "JxW=" << fe_values.JxW(0) << std::endl;

        // now print some data on the shape
        // functions
        for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
          deallog << "shape_function " << i << ":" << std::endl
                  << "  phi=" << fe_values.shape_value(i,0) << std::endl
                  << "  grad phi=" << fe_values.shape_grad(i,0) << std::endl
                  << "  grad^2 phi=" << fe_values.shape_hessian(i,0) << std::endl;

        deallog << std::endl;
      }
}


template<int dim>
void mapping_test()
{
  deallog << "dim=" << dim << std::endl;

  std::vector<Mapping<dim> *> mapping_ptr;
  std::vector<std::string> mapping_strings;

  MappingCartesian<dim> mapping;
  std::string mapping_name = "MappingCartesian";

  Triangulation<dim> tria;
  GridGenerator::hyper_cube (tria);

  FE_Q<dim> fe_q4(2);
  DoFHandler<dim> dof(tria);
  dof.distribute_dofs(fe_q4);

  typename DoFHandler<dim>::cell_iterator cell = dof.begin_active();
  {
    std::ostringstream ost;
    ost << "MappingFace" << dim << "d-"
        << mapping_name;
    deallog << ost.str() << std::endl;
    plot_faces(mapping, fe_q4, cell);
  }

  {
    std::ostringstream ost;
    ost << "MappingSubface" << dim << "d-"
        << mapping_name;
    deallog << ost.str() << std::endl;
    plot_subfaces(mapping, fe_q4, cell);
  }
}



int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision(PRECISION);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // -----------------------
  // Tests for dim=1
  // -----------------------
  mapping_test<1>();
}

