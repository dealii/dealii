// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
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


// FE_RaviartThomasNodal
// Show the shape functions in support points

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/fe_values.h>

#include <vector>
#include <fstream>
#include <string>

#include <cmath>


template <int dim>
void
check_support_points (const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl;

  const std::vector< Point<dim > > &points = fe.get_generalized_support_points();
  std::vector<double> weights(points.size());

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  MappingCartesian<dim> mapping;

  Quadrature<dim> support_quadrature(points, weights);
  UpdateFlags flags = update_values | update_gradients;
  FEValues<dim> vals(mapping, fe, support_quadrature, flags);
  vals.reinit(dof.begin_active());

  for (unsigned int k=0; k<points.size(); ++k)
    {
      const Point<dim> &p = points[k];
      deallog << p;
      for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
        {
          deallog << '\t';
          for (unsigned int d=0; d<dim; ++d)
            {
              const double sf = fe.shape_value_component(i,p,d);
              deallog << ' ' << (int) rint(60*sf);

              const double diff = std::abs(sf - vals.shape_value_component(i,k,d));
              if (diff > 1.e-12)
                deallog << "Error values" << std::endl;
              Tensor<1,dim> grad = fe.shape_grad_component(i,p,d);
              grad -= vals.shape_grad_component(i,k,d);
              if (grad.norm() > 1.e-12)
                deallog << "Error grads" << std::endl;
            }
        }
      deallog << std::endl;
    }
  deallog << std::endl;
}


template <int dim>
void
check_face_support_points (const FiniteElement<dim> &fe)
{
  deallog << "Face " << fe.get_name() << std::endl;

  const std::vector< Point<dim-1> > & sub_points = fe.get_generalized_face_support_points();
  std::vector<double> weights(sub_points.size());
  Quadrature<dim-1> sub_quadrature(sub_points, weights);
  std::vector< Point<dim> > points(sub_points.size());

  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  MappingCartesian<dim> mapping;

  UpdateFlags flags = update_values | update_gradients;
  FEFaceValues<dim> vals(mapping, fe, sub_quadrature, flags);

  for (unsigned int face=0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    {
      QProjector<dim>::project_to_face(sub_quadrature, face, points);
      vals.reinit(dof.begin_active(), face);

      for (unsigned int k=0; k<points.size(); ++k)
        {
          const Point<dim> &p = points[k];
          deallog << p;
          for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
            {
              deallog << '\t';
              for (unsigned int d=0; d<dim; ++d)
                {
                  const double sf = fe.shape_value_component(i,p,d);
                  deallog << ' ' << (int) rint(60*sf);

                  const double diff = std::abs(sf - vals.shape_value_component(i,k,d));
                  if (diff > 1.e-12)
                    deallog << "Error values" << std::endl;
                  Tensor<1,dim> grad = fe.shape_grad_component(i,p,d);
                  grad -= vals.shape_grad_component(i,k,d);
                  if (grad.norm() > 1.e-12)
                    deallog << "Error grads" << std::endl;
                }
            }
          deallog << std::endl;
        }
      deallog << std::endl;
    }
}


int
main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  FE_RaviartThomasNodal<2> e20(0);
  check_support_points(e20);
  FE_RaviartThomasNodal<3> e30(0);
  check_support_points(e30);
  check_face_support_points(e20);
  FE_RaviartThomasNodal<2> e21(1);
  check_support_points(e21);
  check_face_support_points(e21);
}
