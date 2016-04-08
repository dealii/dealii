// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2015 by the deal.II authors
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



// test the FEValues views and extractor classes. this test is for
// get_function_hessians for vector components and a non-primitive element

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria_boundary_lib.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <fstream>

template<int dim>
class VectorFunction : public Function<dim>
{
public:
  VectorFunction() : Function<dim>(dim) {}
  virtual double value (const Point<dim> &p, const unsigned int component) const;
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template<>
double VectorFunction<2>::value(const Point<2> &p, const unsigned int component) const
{
  Assert (component < 2,  ExcIndexRange (component, 0, 1));

  const double PI = numbers::PI;
  double val = 0.0;
  switch (component)
    {
    case 0:
      val = pow(p(0),3);
      break;
    case 1:
      val = pow(p(1),2)*p(0);
      break;
    }
  return val;
}

template<>
double VectorFunction<3>::value(const Point<3> &p, const unsigned int component) const
{
  Assert (component < 3, ExcIndexRange (component, 0, 2));

  const double PI = numbers::PI;
  double val = 0.0;
  switch (component)
    {
    case 0:
      val = pow(p(0),3);
      break;
    case 1:
      val = pow(p(1),2)*p(0);
      break;
    case 2:
      val = p(2)*p(1)*p(0);
      break;
    }
  return val;
}

template<int dim>
void VectorFunction<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (int i = 0; i < dim; ++i)
    values(i) = value(p, i);
}

template<int dim>
void test (const Triangulation<dim> &tr,
           const FiniteElement<dim> &fe)
{
  deallog << "FE=" << fe.get_name()
          << std::endl;

  DoFHandler<dim> dof_handler(tr);
  dof_handler.distribute_dofs(fe);

  VectorFunction<dim> fe_function;

  MappingQGeneric<dim> mapping(1);

  const QGauss<dim> quadrature(2);
  FEValues<dim> fe_values (mapping, fe, quadrature,
                           update_values | update_gradients | update_3rd_derivatives);

  ConstraintMatrix constraints;
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);
  constraints.close();

  Vector<double> function_vals(dof_handler.n_dofs());
  VectorTools::project(mapping, dof_handler, constraints, quadrature, fe_function, function_vals);

  fe_values.reinit (dof_handler.begin_active());

  std::vector<Tensor<4,dim> > selected_vector_values (quadrature.size());
  std::vector<std::vector<Tensor<3,dim> > >
  vector_values (quadrature.size(),
                 std::vector<Tensor<3,dim> >(fe.n_components()));

  fe_values.get_function_third_derivatives (function_vals, vector_values);

  for (unsigned int c=0; c<fe.n_components(); ++c)
    // use a vector extractor if there
    // are sufficiently many components
    // left after the current component
    // 'c'
    if (c+dim <= fe.n_components())
      {
        FEValuesExtractors::Vector vector_components (c);
        fe_values[vector_components].get_function_third_derivatives (function_vals,
            selected_vector_values);
        deallog << "component=" << c << std::endl;

        for (unsigned int q=0; q<fe_values.n_quadrature_points; ++q)
          for (unsigned int d=0; d<dim; ++d)
            {
              for (unsigned int e=0; e<dim; ++e)
                for (unsigned int f=0; f<dim; ++f)
                  for (unsigned int g=0; g<dim; ++g)
                    deallog << selected_vector_values[q][d][e][f][g] << (e<dim-1 && f<dim-1 && g<dim-1? ", " : "; ");
              deallog << std::endl;
              Assert ((selected_vector_values[q][d] - vector_values[q][c+d]).norm()
                      <= 1e-12 * selected_vector_values[q][d].norm(),
                      ExcInternalError());
            }
      }
}



template<int dim>
void test_hyper_sphere()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_ball(tr);

  static const HyperBallBoundary<dim> boundary;
  tr.set_boundary (0, boundary);


  const unsigned int order = 3;

  test(tr, FESystem<dim> (FE_Q<dim>(QIterated<1>(QTrapez<1>(),order)), dim) );
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (3);

  deallog.attach(logfile);
  deallog.threshold_double(1.e-5);

  test_hyper_sphere<2>();
  test_hyper_sphere<3>();
}
