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


// compare the FE_RaviartThomasNodal and FE_RaviartThomas elements

// compare the shape funcions and shape values after converting to the
// same basis.

// Summary: the different Raviart-Thomas implementations use the same
// polynomial spaces, but different basis functions. Here, we convert
// between the bases and test if the resulting functions are the same
// point-wise.

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/fe_values.h>

#include <fstream>

// This function copied from FERaviartThomasNodal. nodes is the
// element having the support points and the value of other in these
// points is computed.
template <int dim>
void
initialize_node_matrix (const FiniteElement<dim> &other,
                        const FiniteElement<dim> &nodes,
                        FullMatrix<double> &N)
{
  const unsigned int n_dofs = other.dofs_per_cell;
  Assert (n_dofs == nodes.dofs_per_cell,
          ExcDimensionMismatch(n_dofs, nodes.dofs_per_cell));

  N.reinit(n_dofs, n_dofs);

  const std::vector<Point<dim> > &unit_support_points = nodes.get_generalized_support_points();

  // The curent node functional index
  unsigned int current = 0;

  // For each face and all quadrature
  // points on it, the node value is
  // the normal component of the
  // shape function, possibly
  // pointing in negative direction.
  for (unsigned int face = 0; face<GeometryInfo<dim>::faces_per_cell; ++face)
    for (unsigned int k=0; k<other.dofs_per_face; ++k)
      {
        for (unsigned int i=0; i<n_dofs; ++i)
          N(current,i) = other.shape_value_component(
                           i, unit_support_points[current],
                           GeometryInfo< dim >::unit_normal_direction[face]);
        ++current;
      }
  // Interior degrees of freedom in each direction
  const unsigned int n_cell = (n_dofs - current) / dim;

  for (unsigned int d=0; d<dim; ++d)
    for (unsigned int k=0; k<n_cell; ++k)
      {
        for (unsigned int i=0; i<n_dofs; ++i)
          N(current,i) = other.shape_value_component(i, unit_support_points[current], d);
        ++current;
      }
  Assert (current == n_dofs, ExcInternalError());
}


template <int dim>
void
compare_shapes (const FiniteElement<dim> &other,
                const FiniteElement<dim> &nodes,
                FullMatrix<double> &M)
{
  QGauss<dim> quadrature(other.degree+1);
  Table<3,double> other_values(quadrature.size(), other.dofs_per_cell, dim);
  Table<3,double> nodes_values(quadrature.size(), other.dofs_per_cell, dim);
  Table<3,Tensor<1,dim> > other_grads(quadrature.size(), other.dofs_per_cell, dim);
  Table<3,Tensor<1,dim> > nodes_grads(quadrature.size(), other.dofs_per_cell, dim);
  for (unsigned int k=0; k<quadrature.size(); ++k)
    for (unsigned int i=0; i<other.dofs_per_cell; ++i)
      for (unsigned int d=0; d<dim; ++d)
        {
          other_values[k][i][d] = other.shape_value_component(i,quadrature.point(k),d);
          nodes_values[k][i][d] = nodes.shape_value_component(i,quadrature.point(k),d);
          other_grads[k][i][d] = other.shape_grad_component(i,quadrature.point(k),d);
          nodes_grads[k][i][d] = nodes.shape_grad_component(i,quadrature.point(k),d);
        }

  for (unsigned int k=0; k<quadrature.size(); ++k)
    {
      for (unsigned int i=0; i<other.dofs_per_cell; ++i)
        for (unsigned int d=0; d<dim; ++d)
          {
            double value = other_values[k][i][d];
            Tensor<1,dim> grad = other_grads[k][i][d];
            for (unsigned int j=0; j<other.dofs_per_cell; ++j)
              {
                value -= M(j,i) * nodes_values[k][j][d];
                grad -= M(j,i) * nodes_grads[k][j][d];
              }
            deallog << '.';
            if (std::fabs(value) > 1.e-12)
              deallog << "Error value\t" << k << '\t' << i << '\t' << d << '\t' << value
                      << std::endl;
            if (grad.norm() > 1.e-12)
              deallog << "Error grad\t" << k << '\t' << i << '\t' << d << '\t' << grad
                      << '\t' << other_grads[k][i][d]
                      << std::endl;
          }
      deallog << std::endl;
    }
}


template<int dim>
void
test (unsigned int degree)
{
  FE_RaviartThomas<dim> rt1(degree);
  FE_RaviartThomasNodal<dim> rtn1(degree);
  FullMatrix<double> N;
  initialize_node_matrix(rt1, rtn1, N);
  compare_shapes(rt1, rtn1, N);
}

int main()
{
  std::ofstream logfile ("output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test<2>(0);
  test<2>(1);
  test<2>(2);
}
