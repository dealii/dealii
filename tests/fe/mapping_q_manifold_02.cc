// ---------------------------------------------------------------------
//
// Copyright (C) 2016 by the deal.II authors
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

// Test that high order MappingQ on a manifold puts the intermediate points
// along the manifold

#include "../tests.h"

#include <fstream>
#include <deal.II/base/function.h>
#include <deal.II/base/logstream.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/manifold.h>
#include <deal.II/grid/manifold_lib.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/fe_nothing.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;

double f_x(double x_m)
{
  double x = x_m*1000.0;
  double y_m = 0.0;

  if (x <= 9.0)
    y_m = 0.001*(-28.0 + std::min(28.0, 2.8e1 + 6.775070969851e-3*x*x - 2.124527775800e-3*x*x*x));

  else if (x > 9.0 && x <= 14.0)
    y_m = 0.001*(-28.0 + 2.507355893131e1 + 9.754803562315e-1*x - 1.016116352781e-1*x*x + 1.889794677828e-3*x*x*x);

  else if (x > 14.0 && x <= 20.0)
    y_m = 0.001*(-28.0 + 2.579601052357e1 + 8.206693007457e-1*x - 9.055370274339e-2*x*x + 1.626510569859e-3*x*x*x);

  else if (x > 20.0 && x <= 30.0)
    y_m = 0.001*(-28.0 + 4.046435022819e1 - 1.379581654948*x + 1.945884504128e-2*x*x - 2.070318932190e-4*x*x*x);

  else if (x > 30.0 && x <= 40.0)
    y_m = 0.001*(-28.0 + 1.792461334664e1 + 8.743920332081e-1*x - 5.567361123058e-2*x*x + 6.277731764683e-4*x*x*x);

  else if (x > 40.0 && x <= 54.0)
    y_m = 0.001*(-28.0 + std::max(0.0, 5.639011190988e1 - 2.010520359035*x + 1.644919857549e-2*x*x + 2.674976141766e-5*x*x*x));

  else if (x > 54.0)
    y_m = 0.001*(-28.0);
  return y_m;
}

template<int dim>
class PushForward : public Function<dim>
{
public:
  PushForward ()
    :
    Function<dim>(dim, 0.),
    h (0.028),
    x_max (4.5*h),
    y_max (2.036*h),
    z_max (4.5*h),
    y_FoR (h)
  {}

  virtual ~PushForward() {};

  virtual double value (const Point<dim> &p,const unsigned int component = 0) const;

private:
  const double h;

  // data from initial block
  const double x_max;
  const double y_max;
  const double z_max;

  const double y_FoR;


};


template<int dim>
double PushForward<dim>::value(const Point<dim> &p,const unsigned int component) const
{
  double result = 0;

  // x component
  if (component == 0)
    result = p[0];

  // y component
  else if (component == 1)
    {
      if (p[0] <= x_max/2.0)
        result = p[1] + (1 - (p[1] - y_FoR)/y_max)*f_x(p[0]);

      else if (p[0] > x_max/2.0)
        result = p[1] + (1 - (p[1] - y_FoR)/y_max)*f_x(x_max - p[0]);
    }

  // z component
  else if (component == 2)
    result = p[2];

  return result;
}

template<int dim>
class PullBack : public Function<dim>
{
public:
  PullBack ()
    :
    Function<dim>(dim, 0.),
    h (0.028),
    x_max (4.5*h),
    y_max (2.036*h),
    z_max (4.5*h),
    y_FoR (h)
  {}

  virtual ~PullBack() {};

  virtual double value (const Point<dim> &p,const unsigned int component = 0) const;

private:
  const double h;

  const double x_max;
  const double y_max;
  const double z_max;

  const double y_FoR;
};

template<int dim>
double PullBack<dim>::value(const Point<dim> &p,const unsigned int component) const
{
  double result = 0;

  // x component
  if (component == 0)
    result = p[0];

  // y component
  else if (component == 1)
    {
      if (p[0] <= x_max/2.0)
        result = (p[1] - f_x(p[0])*(1 + y_FoR/y_max))/(1.0 - f_x(p[0])/y_max);
      else if (p[0] > x_max/2.0)
        result = (p[1] - f_x(x_max - p[0])*(1 + y_FoR/y_max))/(1.0 - f_x(x_max - p[0])/y_max);
    }

  // z component
  else if (component == 2)
    result = p[2];

  return result;
}

template <int dim>
void create_tria(Triangulation<dim> &triangulation, const Manifold<dim> &manifold)
{
  const double h = 0.028;
  std::vector<unsigned int> refinements(dim,1);
  refinements[1] = 2;

  Point<dim> p1, p2;
  p2[0] = 4.5*h;//9.0*h;
  p1[1] = h;
  p2[1] = 2.018*h;//2.036*h;
  if (dim == 3)
    {
      p1[2] = -2.25*h;
      p2[2] = 2.25*h;
    }
  GridGenerator::hyper_rectangle(triangulation, p1, p2);

  for (typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
       cell != triangulation.end(); ++cell)
    cell->set_all_manifold_ids(111);
  triangulation.set_manifold(111, manifold);

  triangulation.refine_global(1);
}

template <int dim>
void test()
{
  deallog << "dim: " << dim << std::endl;

  PushForward<dim> push_forward;
  PullBack<dim> pull_back;
  FunctionManifold<dim,dim,dim> manifold(push_forward, pull_back);

  Triangulation<dim> triangulation;
  create_tria(triangulation, manifold);

  for (unsigned mapping_p = 4; mapping_p <= 5; ++mapping_p)
    {
      deallog << "Mapping degree: " << mapping_p << std::endl;
      MappingQ<dim> mapping(mapping_p, true);
      std::vector<Point<dim> > points(Utilities::fixed_power<dim>(2));
      for (unsigned int i=0, c=0; i<(dim==2?1:2); ++i)
        for (unsigned int j=0; j<2; ++j)
          for (unsigned int k=0; k<2; ++k, ++c)
            {
              points[c][0] = 0.1 + 0.8*k;
              points[c][1] = 0.1 + 0.8*j;
              if (dim == 3)
                points[c][2] = 0.1 + 0.8*i;
            }
      FE_Nothing<dim> dummy;
      Quadrature<dim> quad(points);
      FEValues<dim> fe_values(mapping, dummy, quad, update_quadrature_points);

      for (typename Triangulation<dim>::active_cell_iterator
           cell=triangulation.begin_active(); cell != triangulation.end(); ++cell)
        {
          fe_values.reinit(cell);
          deallog << "Points for cell with first point: " << cell->vertex(0)
                  << std::endl;
          for (unsigned int q=0; q<quad.size(); ++q)
            deallog << fe_values.quadrature_point(q) << "   ";
          deallog << std::endl;
        }

      deallog << std::endl;
    }
}

int main ()
{
  deallog << std::setprecision (5);
  deallog.attach (std::cout);
  deallog.depth_console (0);
  deallog.threshold_double (1e-12);

  test<2>();
  test<3>();
}
