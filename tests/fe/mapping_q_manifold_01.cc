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

// Test high order MappingQ on a ChartManifold.

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
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>

using namespace dealii;

double zvalue (const double x, const double y)
{
  double xh = x * 5., yh = y * 5.;
  return (xh * exp(-xh*xh - yh*yh)) / 10.;
}

template<int dim>
class Geometry: public ChartManifold<dim>
{
public:
  virtual Point<dim> pull_back(const Point<dim> &space_point) const;
  virtual Point<dim> push_forward(const Point<dim> &chart_point) const;
};

template<int dim>
Point<dim> Geometry<dim>::pull_back(const Point<dim> &space_point) const
{
  const double d = space_point[dim - 1];
  const double z = zvalue(space_point[0], dim == 3 ? space_point[1] : 0);

  double d_hat = 0.;
  if ((d - z) <= 0)
    d_hat = (d - z) / (1. + z);
  else
    d_hat = (d - z) / (1. - z);

  Point<dim> p;
  for (unsigned i = 0; i < dim - 1; ++i)
    p[i] = space_point[i];
  p[dim - 1] = d_hat;

  return p;
}

template<int dim>
Point<dim> Geometry<dim>::push_forward(const Point<dim> &chart_point) const
{
  const double d_hat = chart_point[dim - 1];
  const double z = zvalue(chart_point[0], dim == 3 ? chart_point[1] : 0);

  double d = 0.;
  if (d_hat <= 0)
    d = d_hat + (d_hat + 1.) * z;
  else
    d = d_hat - (d_hat - 1.) * z;

  Point<dim> p;
  for (unsigned i = 0; i < dim - 1; ++i)
    p[i] = chart_point[i];
  p[dim - 1] = d;

  return p;
}

template<int dim>
class VectorFunction : public Function<dim>
{
public:
  VectorFunction(unsigned n_components):
    Function<dim>(n_components)
  {}
  virtual double value (const Point<dim> &p, const unsigned int component) const;
  virtual void value (const Point<dim> &p, Vector<double> &values) const;
  virtual void vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template<>
double VectorFunction<3>::value(const Point<3> &p, const unsigned int component) const
{
  return (1 - p(0)*p(0)) * (1 - p(1)*p(1)) * (1 - p(2)*p(2));
}

template<>
double VectorFunction<2>::value(const Point<2> &p, const unsigned int component) const
{
  return (1 - p(0)*p(0)) * (1 - p(1)*p(1));
}

template<int dim>
void VectorFunction<dim>::value(const Point<dim> &p, Vector<double> &values) const
{
  for (unsigned i = 0; i < values.size(); ++i)
    values[i] = value(p, i);
}

template<int dim>
void VectorFunction<dim>::vector_value(const Point<dim> &p, Vector<double> &values) const
{
  for (int i = 0; i < dim; ++i)
    values(i) = value(p, i);
}

template <int dim>
void create_tria(Triangulation<dim> &triangulation, const Geometry<dim> &geometry)
{
  GridGenerator::hyper_cube (triangulation, -1., 1.);
  triangulation.refine_global(3);

  GridTools::transform (std_cxx11::bind(&Geometry<dim>::push_forward,
                                        std_cxx11::cref(geometry),
                                        std_cxx11::_1),
                        triangulation);

  triangulation.set_manifold(0, geometry);
  for (Triangulation<3>::active_cell_iterator cell=triangulation.begin_active();
       cell!=triangulation.end(); ++cell)
    cell->set_all_manifold_ids(0);
}

template <int dim>
void test(const FiniteElement<dim> &fe)
{
  deallog << "dim: " << dim << "\t" << fe.get_name() << std::endl;
  deallog << "Mapping degree\t||u-u_h||" << std::endl;

  Geometry<dim> geometry;
  Triangulation<dim> triangulation;
  create_tria(triangulation, geometry);

  DoFHandler<dim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  VectorFunction<dim> fe_function(fe.n_components());

  for (unsigned mapping_p = 1; mapping_p <= 5; ++mapping_p)
    {
      MappingQ<dim> mapping(mapping_p, true);

      ConstraintMatrix constraints;
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      constraints.close();

      Vector<double> v(dof_handler.n_dofs());
      VectorTools::project(mapping, dof_handler, constraints,
                           QGauss<dim>(fe.degree + 1), fe_function, v);

      Vector<float> diff(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping, dof_handler, v, fe_function, diff,
                                        QGauss<dim>(fe.degree + 2), VectorTools::L2_norm);

      deallog << mapping.get_degree() << "\t"
              << diff.l2_norm() << std::endl;

      DataOut<dim> data_out;
      std::ofstream output("output_" + Utilities::int_to_string(mapping_p) + ".vtk");
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector(v, "v");
      data_out.add_data_vector(diff, "e");
      data_out.build_patches(mapping, fe.degree);
      data_out.write_vtk(output);
      output.close();
    }
}

int main ()
{
  deallog << std::setprecision (5);
  deallog.attach (std::cout);
  deallog.depth_console (0);
  deallog.threshold_double (1e-12);

  const static unsigned dim = 3;

  for (unsigned p = 1; p < 5; ++p)
    {
      test<dim>(FE_Q<dim>(QGaussLobatto<1>(p+1)));

      //test<dim>(FE_RaviartThomas<dim> (p));
      //test<dim>(FE_Nedelec<dim> (p));
      //test<dim>(FESystem<dim> (FE_Q<dim>(p), dim));
    }
}
