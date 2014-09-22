// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/lac/vector.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

#include <fstream>


/**
 * A vector-valued polynomial for testing RT elements.
 */
template<int dim>
class TestFunction : public Function<dim>
{
public:

  /// Construct a polynomial of degree p
  TestFunction(unsigned int degree);

  virtual void vector_value_list (const std::vector<Point<dim> > &points,
                                  std::vector<Vector<double> >   &values) const;
private:
  unsigned int degree;
};


template<int dim>
TestFunction<dim>::TestFunction (unsigned int p)
  :
  Function<dim>(dim), degree(p)
{}


template<int dim>
void
TestFunction<dim>::vector_value_list (const std::vector<Point<dim> > &points,
                                      std::vector<Vector<double> >   &values) const
{
  for (unsigned int k=0; k<points.size(); ++k)
    {
      if (degree < 2)
        {
          for (unsigned int d=0; d<dim; ++d)
            values[k](d) = points[k](d) - d;
        }
      else
        {
          // Base of the function is
          // the distance to a
          // different point in each
          // component
          for (unsigned int d=0; d<dim; ++d)
            {
              Point<dim> p = points[k];
              for (unsigned int dd=0; dd<dim; ++dd)
                p(dd) -= d;
              const double r2 = p.square();
              values[k](d) = std::pow(r2, (int) degree/2);
            }
        }
    }
}



template<int dim>
double integrate_error(const DoFHandler<dim> &dof,
                       FEFaceValues<dim> &fe,
                       const Vector<double> &u,
                       const Function<dim> &f)
{
  double result = 0.;
  std::vector<Vector<double> > f_values(fe.n_quadrature_points, Vector<double>(dim));
  std::vector<Vector<double> > fe_values(fe.n_quadrature_points, Vector<double>(dim));

  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end(); ++cell)
    {
      for (unsigned int face=0 ; face != GeometryInfo<dim>::faces_per_cell; ++face)
        {
          if (!cell->at_boundary(face)) continue;

          fe.reinit(cell, face);
          f.vector_value_list(fe.get_quadrature_points(), f_values);
          fe.get_function_values(u, fe_values);
          for (unsigned int k=0; k<fe.n_quadrature_points; ++k)
            {
              double diff = 0.;
              for (unsigned int d=0; d<dim; ++d)
                diff += fe.normal_vector(k)(d) * (f_values[k](d) - fe_values[k](d));
              result += fe.JxW(k) * diff * diff;
            }
        }
    }
  return result;
}


template<int dim>
void test_projection (const Triangulation<dim> &tr,
                      const FiniteElement<dim> &fe)
{
  deallog << fe.get_name() << std::endl << "Cells: " << tr.n_active_cells() << std::endl;

  const unsigned int degree = fe.tensor_degree();

  DoFHandler<dim> dof(tr);
  dof.distribute_dofs(fe);

  QGauss<dim-1> quadrature(degree+2);
  MappingQ1<dim> mapping;

  TestFunction<dim> f(degree-1);
  std::map<types::global_dof_index, double> boundary_constraints;
  typename FunctionMap<dim>::type boundary_map;
  for (types::boundary_id i=0; i<255; ++i)
    boundary_map[i] = &f;
  VectorTools::project_boundary_values(mapping, dof, boundary_map, quadrature,
                                       boundary_constraints);

  deallog << "Constraints: " << boundary_constraints.size() << std::endl;

  // Fill a vector with the projected
  // boundary values
  Vector<double> u(dof.n_dofs());
  u = -1.;
  for (typename std::map<types::global_dof_index, double>::const_iterator
       i = boundary_constraints.begin(); i != boundary_constraints.end(); ++i)
    u(i->first) = i->second;

  FEFaceValues<dim> feval(mapping, fe, quadrature,
                          update_quadrature_points
                          | update_normal_vectors
                          | update_JxW_values
                          | update_values);
  double err = integrate_error(dof, feval, u, f);
  deallog << err << std::endl;
}


template<int dim>
void test_hyper_cube(const FiniteElement<dim> &fe)
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  test_projection(tr, fe);
}


int main()
{
  std::ofstream logfile ("output");
  deallog << std::setprecision (2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);
  deallog.threshold_double(1.e-12);

  FE_RaviartThomasNodal<2> rt21(1);
  test_hyper_cube(rt21);
  FE_RaviartThomasNodal<2> rt22(2);
  test_hyper_cube(rt22);
  FE_RaviartThomasNodal<2> rt23(3);
  test_hyper_cube(rt23);
  FE_RaviartThomasNodal<2> rt24(4);
  test_hyper_cube(rt24);
}
