// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------


#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_abf.h>
#include <deal.II/fe/fe_nedelec.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/vector_tools.h>

#include "../tests.h"

/*
 * This program projects a function into FE spaces defined on meshes
 * consisting of: rectangular cells, affine cells, non-affine cells
 *
 * The error, curl and divergence are then numerically calculated on a
 * series of globally refined meshes and get output.
 *
 * Among FE spaces tested are: FE_ABF, FE_Nedelec, FE_RaviartThomas,
 * FE_Q^dim (via FESystem)
 *
 * Alexander Grayver
 */


static const Point<2> vertices_nonaffine[] = {
  Point<2>(-1., -1.),
  Point<2>(0., -1.),
  Point<2>(1., -1.),

  Point<2>(-1., 0.),
  Point<2>(0.3, 0.3),
  Point<2>(1., 0.),

  Point<2>(-1., 1.),
  Point<2>(0., 1.),
  Point<2>(1., 1.),
};

static const Point<2> vertices_affine[] = {
  Point<2>(-1.4, -1.),
  Point<2>(-0.4, -1.),
  Point<2>(0.6, -1.),

  Point<2>(-1.2, 0.),
  Point<2>(-0.2, 0.),
  Point<2>(0.8, 0.),

  Point<2>(-1., 1.),
  Point<2>(0., 1.),
  Point<2>(1., 1.),
};

static const Point<2> vertices_rectangular[] = {
  Point<2>(-1., -1.),
  Point<2>(0., -1.),
  Point<2>(1., -1.),

  Point<2>(-1., 0.),
  Point<2>(0., 0.),
  Point<2>(1., 0.),

  Point<2>(-1., 1.),
  Point<2>(0., 1.),
  Point<2>(1., 1.),
};

static const unsigned n_vertices =
  sizeof(vertices_rectangular) / sizeof(vertices_rectangular[0]);
static const unsigned n_cells = 4;

template <int dim>
class VectorFunction : public Function<dim>
{
public:
  VectorFunction()
    : Function<dim>(dim)
  {}
  virtual double
  value(const Point<dim> &p, const unsigned int component) const;
  virtual void
  vector_value(const Point<dim> &p, Vector<double> &values) const;
};

template <>
double
VectorFunction<2>::value(const Point<2> &p, const unsigned int component) const
{
  Assert(component < 2, ExcIndexRange(component, 0, 1));

  const double PI  = numbers::PI;
  double       val = 0.0;
  switch (component)
    {
      case 0:
        val = cos(PI * p(0)) * sin(PI * p(1));
        break;
      case 1:
        val = -sin(PI * p(0)) * cos(PI * p(1));
        break;
    }
  return val;
}

template <int dim>
void
VectorFunction<dim>::vector_value(const Point<dim> &p,
                                  Vector<double> &  values) const
{
  for (int i = 0; i < dim; ++i)
    values(i) = value(p, i);
}

void create_tria(Triangulation<2> &triangulation,
                 const Point<2> *  vertices_parallelograms)
{
  const std::vector<Point<2>> vertices(&vertices_parallelograms[0],
                                       &vertices_parallelograms[n_vertices]);

  static const int cell_vertices[][GeometryInfo<2>::vertices_per_cell] = {
    {0, 1, 3, 4}, {1, 2, 4, 5}, {3, 4, 6, 7}, {4, 5, 7, 8}};

  std::vector<CellData<2>> cells(n_cells, CellData<2>());
  for (unsigned i = 0; i < cells.size(); ++i)
    {
      for (const unsigned int j : GeometryInfo<2>::vertex_indices())
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  triangulation.create_triangulation(vertices, cells, SubCellData());
  triangulation.refine_global(1);
}

template <int dim>
void
test(const FiniteElement<dim> &fe,
     unsigned                  n_cycles,
     bool                      global,
     const Point<dim> *        vertices_parallelograms)
{
  deallog << "dim: " << dim << "\t" << fe.get_name() << std::endl;
  deallog << "DoFs\t\t||u-u_h||\tcurl(u_h)\tdiv(u_h)" << std::endl;

  Triangulation<dim> triangulation;
  create_tria(triangulation, vertices_parallelograms);

  DoFHandler<dim> dof_handler(triangulation);

  VectorFunction<dim>              fe_function;
  const FEValuesExtractors::Vector vec(0);
  const QGauss<dim>                quadrature(fe.degree + 1);
  const unsigned int               n_q_points = quadrature.size();
  MappingQ<dim>                    mapping(1);
  // MappingQGeneric<dim> mapping(1);
  std::vector<double>                                         div_v(n_q_points);
  std::vector<typename FEValuesViews::Vector<dim>::curl_type> curl_v(
    n_q_points);

  for (unsigned cycle = 0; cycle < n_cycles; ++cycle)
    {
      dof_handler.distribute_dofs(fe);

      AffineConstraints<double> constraints;
      DoFTools::make_hanging_node_constraints(dof_handler, constraints);
      constraints.close();

      Vector<double> v(dof_handler.n_dofs());
      VectorTools::project(
        mapping, dof_handler, constraints, quadrature, fe_function, v);

      Vector<float> diff(triangulation.n_active_cells());
      VectorTools::integrate_difference(mapping,
                                        dof_handler,
                                        v,
                                        fe_function,
                                        diff,
                                        QGauss<dim>(fe.degree + 2),
                                        VectorTools::L2_norm);

      typename FEValuesViews::Vector<dim>::curl_type total_curl;
      double                                         total_div = 0;
      total_curl *= 0.;

      FEValues<dim> fe_values(mapping,
                              fe,
                              quadrature,
                              update_JxW_values | update_quadrature_points |
                                update_values | update_gradients);
      unsigned int  cell_index = 0;

      for (typename DoFHandler<dim>::active_cell_iterator cell =
             dof_handler.begin_active();
           cell != dof_handler.end();
           ++cell, ++cell_index)
        {
          fe_values.reinit(cell);
          const std::vector<double> &JxW_values = fe_values.get_JxW_values();
          fe_values[vec].get_function_divergences(v, div_v);
          fe_values[vec].get_function_curls(v, curl_v);
          for (const auto q_point : fe_values.quadrature_point_indices())
            {
              total_div += JxW_values[q_point] * div_v[q_point];
              total_curl += JxW_values[q_point] * curl_v[q_point];
            }
        }

      deallog << dof_handler.n_dofs() << "\t\t" << diff.l2_norm() << "\t"
              << total_curl << "\t" << total_div << std::endl;

      if (global)
        triangulation.refine_global();
      else
        {
          GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                          diff,
                                                          0.3,
                                                          0.0);
          triangulation.prepare_coarsening_and_refinement();
          triangulation.execute_coarsening_and_refinement();
        }
    }
}

int
main()
{
  initlog();
  deallog << std::setprecision(7) << std::fixed;

  const static unsigned dim      = 2;
  unsigned              order    = 1;
  unsigned              n_cycles = 4;

  deallog << "2d\nRectangular grid:\n";

  const Point<dim> *vertices = &vertices_rectangular[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);
  test<dim>(FESystem<dim>(FE_Q<dim>(order), dim), n_cycles, true, vertices);
  test<dim>(FE_ABF<dim>(order), n_cycles, true, vertices);

  deallog << "\nAffine grid:\n";

  vertices = &vertices_affine[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);
  test<dim>(FESystem<dim>(FE_Q<dim>(order), dim), n_cycles, true, vertices);
  test<dim>(FE_ABF<dim>(order), n_cycles, true, vertices);

  deallog << "\nNon-affine grid:\n";

  vertices = &vertices_nonaffine[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);
  test<dim>(FESystem<dim>(FE_Q<dim>(order), dim), n_cycles, true, vertices);
  test<dim>(FE_ABF<dim>(order), n_cycles, true, vertices);

  deallog << std::endl;
}
