// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2014 by the deal.II authors
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


// Check Singularity functions and SingularityGrad functions for

#include "../tests.h"
#include <deal.II/base/data_out_base.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/job_identifier.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function_lib.h>
#include <deal.II/base/auto_derivative_function.h>
#include <deal.II/base/flow_function.h>
#include <deal.II/lac/vector.h>

#include "functions.h"

#include <vector>
#include <iomanip>
#include <fstream>
#include <string>

template<int dim>
void
check_function_consistency(
  const Function<dim> &f,
  const Function<dim> &gradf,
  unsigned int sub)
{
  QMidpoint<1> mid;
  QIterated<dim> quadrature(mid, sub);

  std::vector<Tensor<1,dim> > fg1(quadrature.size());
  std::vector<std::vector<Tensor<1,dim> > > fg2(quadrature.size(), std::vector<Tensor<1,dim> >(1));

  std::vector<double> g1(quadrature.size());
  std::vector<Vector<double> > g2(quadrature.size(), Vector<double>(dim));

  // Derivative values are in fg1,
  // fg2, g1, and g2
  f.gradient_list(quadrature.get_points(), fg1, 0);
  f.vector_gradient_list(quadrature.get_points(), fg2);
  gradf.vector_value_list(quadrature.get_points(), g2);
  deallog << "Gradient consistency ";

  for (unsigned int d=0; d<dim; ++d)
    {
      gradf.value_list(quadrature.get_points(), g1, d);
      for (unsigned int i=0; i<g1.size(); ++i)
        {
          if (std::fabs(g1[i]-g2[i](d)) > 1.e-14)
            deallog << ' ' << i << ':' << g1[i]-g2[i](d);
          if (std::fabs(g1[i]-fg2[i][0][d]) > 1.e-14)
            deallog << ' ' << i << ';' << g1[i]-fg2[i][0][d];
          if (std::fabs(g1[i]-fg1[i][d]) > 1.e-14)
            deallog << ' ' << i << '#' << g1[i]-fg1[i][d];
        }
    }
  deallog << std::endl;
}


template<int dim>
void
check_function_derivative(const Functions::FlowFunction<dim> &f,
                          unsigned int sub,
                          std::ostream &out)
{
  DerivativeTestFunction<dim> dtest1(f, 1.e-2);
  DerivativeTestFunction<dim> dtest2(f, 2.e-2);
  // Prepare a vector with a single
  // patch stretching over the cube
  // [-1,1]^dim
  std::vector<DataOutBase::Patch<dim,dim> > patches(1);
  unsigned int vertex_number = 0;
  for (unsigned int iz=0; iz < ((dim>2) ? 2 : 1) ; ++iz)
    for (unsigned int iy=0; iy < ((dim>1) ? 2 : 1) ; ++iy)
      for (unsigned int ix=0; ix < 2 ; ++ix)
        {
          if (dim>0) patches[0].vertices[vertex_number](0) = -1. + 2.*ix;
          if (dim>1) patches[0].vertices[vertex_number](1) = -1. + 2.*iy;
          if (dim>2) patches[0].vertices[vertex_number](2) = -1. + 2.*iz;
          ++vertex_number;
        }
  for (unsigned int i=0; i<GeometryInfo<dim>::faces_per_cell; ++i)
    patches[0].neighbors[i] = numbers::invalid_unsigned_int;
  patches[0].patch_index = 0;
  patches[0].n_subdivisions = sub;
  patches[0].points_are_available = false;

  vertex_number = 1;
  for (unsigned int d=0; d<dim; ++d)
    vertex_number *= (sub+1);
  patches[0].data.reinit(f.n_components, vertex_number);

  // Build the vector of quadrature points;
  std::vector<Point<dim> > points(vertex_number);
  const double h = 2./sub;
  vertex_number = 0;
  for (unsigned int iz=0; iz <= ((dim>2) ? sub : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? sub : 0) ; ++iy)
      for (unsigned int ix=0; ix <= sub ; ++ix)
        {
          if (dim>0) points[vertex_number](0) = -1.+ix*h;
          if (dim>1) points[vertex_number](1) = -1.+iy*h;
          if (dim>2) points[vertex_number](2) = -1.+iz*h;
          ++vertex_number;
        }

  std::vector<Vector<double> > values(points.size(), Vector<double>(f.n_components));
  std::vector<std::vector<double> > values2(f.n_components,
                                            std::vector<double>(points.size()));
  f.vector_value_list(points, values);
  f.vector_values(points, values2);
  for (unsigned int i=0; i<values.size(); ++i)
    for (unsigned int j=0; j<values[i].size(); ++j)
      {
        // generate data, but
        // truncate too small values
        // to avoid output that
        // depends on round-off
        if (std::fabs (values[i](j)) > 1e-10)
          patches[0].data(j, i) = values[i](j);
        else
          patches[0].data(j, i) = 0;
        if (values[i](j) != values2[j][i])
          deallog << "Error values (" << i << ',' << j << ") : "
                  << values[i](j) << " != " << values2[j][i] << std::endl;
      }

  deallog << "Gradients ";
  // Compute gradients and difference approximations
  std::vector<std::vector<Tensor<1,dim> > >
  gradients(f.n_components, std::vector<Tensor<1,dim> >(points.size()));
  std::vector<std::vector<Tensor<1,dim> > >
  gradients1(points.size(), std::vector<Tensor<1,dim> >(f.n_components));
  std::vector<std::vector<Tensor<1,dim> > >
  gradients2(points.size(), std::vector<Tensor<1,dim> >(f.n_components));

  f.vector_gradients(points, gradients);
  dtest1.vector_gradient_list(points, gradients1);
  dtest2.vector_gradient_list(points, gradients2);

  // Compare gradients and difference quotients
  for (unsigned int k=0; k<gradients.size(); ++k)
    for (unsigned int i=0; i<gradients[k].size(); ++i)
      {
        // Compute difference
        Tensor<1,dim> d1 = gradients1[i][k] - gradients[k][i];
        Tensor<1,dim> d2 = gradients2[i][k] - gradients[k][i];

        // If the difference is
        // already small, we are fine
        if (d1.norm() > 1.e-13)
          {
            // Check for
            // convergence. For full
            // 4th order, gradients2
            // should be 16 times as
            // large, so let's be a
            // bit generous
            if (d2.norm() < 12.* d1.norm())
              {
                deallog << "Gradient error: point " << i
                        << " (" << points[i] << " )"
                        << " comp " << k
//      << " norms " << d1.norm() << " " << d2.norm()
                        << std::endl;
                for (unsigned int d=0; d<dim; ++d)
                  deallog
                      << " " << gradients[k][i][d]
                      << " " << gradients1[i][k][d]
                      << std::endl;
              }
          }
      }
  deallog << "tested" << std::endl;

  // Check if divergence is zero
  deallog << "Divergence ";

  for (unsigned int k=0; k<points.size(); ++k)
    {
      double div = 0.;
      for (unsigned int d=0; d<dim; ++d)
        div += gradients[d][k][d];
      if (std::fabs(div)> 1.e-13)
        deallog << "Divergence " << k << " "
                << div << std::endl;
    }
  deallog << "tested" << std::endl;

  f.vector_laplacian_list(points, values);
  f.vector_laplacians(points, values2);
  double sum = 0.;
  for (unsigned int i=0; i<values.size(); ++i)
    for (unsigned int j=0; j<values[i].size(); ++j)
      {
        sum += values[i](j)*values[i](j);
        if (values[i](j) != values2[j][i])
          deallog << "Error values (" << i << ',' << j << ") : "
                  << values[i](j) << " != " << values2[j][i] << std::endl;
      }
  deallog << "Laplacians " << std::sqrt(sum)/points.size() << std::endl;

  std::vector<std::string> names(f.n_components);
  for (unsigned int i=0; i<names.size(); ++i)
    {
      names[i] = std::string("comp");
    }

  DataOutBase::DXFlags dxflags;
  DataOutBase::GnuplotFlags gflags;
  std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> > vectors;
  if (dim==2)
    DataOutBase::write_gnuplot(patches, names, vectors, gflags, out);
  else
    DataOutBase::write_dx(patches, names, vectors, dxflags, out);
}


int main()
{
  std::string logname = "output";
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console(0);

  Functions::LSingularityFunction fl;
  Functions::LSingularityGradFunction flg;
  // Use odd number of points to
  // avoid lines with
  // discontinuous derivatives.
  deallog << "Functions::LSingularityFunction" << std::endl;
  check_function_value_consistency(fl, 5);
  check_function_gradient_consistency(fl, 5);
  deallog << "Functions::LSingularityGradFunction" << std::endl;
  check_function_value_consistency(flg, 5);
  check_function_consistency(fl, flg, 5);
}

