// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2018 by the deal.II authors
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
 * The error, curl and divergence are then numerically calculated
 * on a series of globally refined meshes and output.
 *
 * Among FE spaces tested are FE_Nedelec and FE_RaviartThomas
 *
 * Alexander Grayver, Maien Hamed
 */


static const Point<3> vertices_affine[] = {
  Point<3>(-1., -1., -1.), Point<3>(0., -1., -1.),  Point<3>(1., -1., -1.),

  Point<3>(-1.2, -1., 0.), Point<3>(-0.2, -1., 0.), Point<3>(0.8, -1., 0.),

  Point<3>(-1.4, -1., 1),  Point<3>(-0.4, -1., 1),  Point<3>(0.6, -1., 1),

  Point<3>(-1., 0., -1.),  Point<3>(0., 0., -1.),   Point<3>(1., 0., -1.),

  Point<3>(-1.2, 0., 0.),  Point<3>(-0.2, 0., 0.),  Point<3>(0.8, 0., 0.),

  Point<3>(-1.4, 0., 1),   Point<3>(-0.4, 0., 1),   Point<3>(0.6, 0., 1),

  Point<3>(-1., 1., -1.),  Point<3>(0., 1., -1.),   Point<3>(1., 1., -1.),

  Point<3>(-1.2, 1., 0.),  Point<3>(-0.2, 1., 0.),  Point<3>(0.8, 1., 0.),

  Point<3>(-1.4, 1., 1),   Point<3>(-0.4, 1., 1),   Point<3>(0.6, 1., 1)};

static const Point<3> vertices_nonaffine[] = {
  Point<3>(-1., -1., -1.), Point<3>(0., -1., -1.),  Point<3>(1., -1., -1.),

  Point<3>(-1., -1., 0.),  Point<3>(0., -1., 0.),   Point<3>(1., -1., 0.),

  Point<3>(-1., -1., 1),   Point<3>(0., -1., 1),    Point<3>(1., -1., 1),

  Point<3>(-1., 0., -1.),  Point<3>(0., 0., -1.),   Point<3>(1., 0., -1.),

  Point<3>(-1., 0., 0.),   Point<3>(0.2, 0.3, 0.1), Point<3>(1., 0., 0.),

  Point<3>(-1., 0., 1),    Point<3>(0., 0., 1),     Point<3>(1., 0., 1),

  Point<3>(-1., 1., -1.),  Point<3>(0., 1., -1.),   Point<3>(1., 1., -1.),

  Point<3>(-1., 1., 0.),   Point<3>(0., 1., 0.),    Point<3>(1., 1., 0.),

  Point<3>(-1., 1., 1.),   Point<3>(0., 1., 1.),    Point<3>(1., 1., 1.)};

static const Point<3> vertices_rectangular[] = {
  Point<3>(-1., -1., -1.), Point<3>(0., -1., -1.), Point<3>(1., -1., -1.),

  Point<3>(-1., -1., 0.),  Point<3>(0., -1., 0.),  Point<3>(1., -1., 0.),

  Point<3>(-1., -1., 1),   Point<3>(0., -1., 1),   Point<3>(1., -1., 1),

  Point<3>(-1., 0., -1.),  Point<3>(0., 0., -1.),  Point<3>(1., 0., -1.),

  Point<3>(-1., 0., 0.),   Point<3>(0., 0., 0.),   Point<3>(1., 0., 0.),

  Point<3>(-1., 0., 1),    Point<3>(0., 0., 1),    Point<3>(1., 0., 1),

  Point<3>(-1., 1., -1.),  Point<3>(0., 1., -1.),  Point<3>(1., 1., -1.),

  Point<3>(-1., 1., 0.),   Point<3>(0., 1., 0.),   Point<3>(1., 1., 0.),

  Point<3>(-1., 1., 1.),   Point<3>(0., 1., 1.),   Point<3>(1., 1., 1.)};

static const unsigned n_vertices =
  sizeof(vertices_rectangular) / sizeof(vertices_rectangular[0]);
static const unsigned n_cells = 8;

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
  virtual Tensor<1, dim>
  gradient(const Point<dim> &p, const unsigned int component = 0) const;
};

template <>
double
VectorFunction<3>::value(const Point<3> &p, const unsigned int component) const
{
  Assert(component < 3, ExcIndexRange(component, 0, 2));

  const double PI  = numbers::PI;
  double       val = 0.0;
  switch (component)
    {
      case 0:
        val = -sin(PI * p(0)) * cos(PI * p(1)) * cos(PI * p(2));
        break;
      case 1:
        val = -cos(PI * p(0)) * sin(PI * p(1)) * cos(PI * p(2));
        break;
      case 2:
        val = 2 * cos(PI * p(0)) * cos(PI * p(1)) * sin(PI * p(2));
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

template <>
Tensor<1, 3>
VectorFunction<3>::gradient(const Point<3> &   p,
                            const unsigned int component) const
{
  const double PI = numbers::PI;
  Tensor<1, 3> val;
  double       x = p(0), y = p(1), z = p(2);

  switch (component)
    {
      case 0:
        val[0] = -PI * cos(PI * x) * cos(PI * y) * cos(PI * z);
        val[1] = PI * cos(PI * z) * sin(PI * x) * sin(PI * y);
        val[2] = -2 * PI * cos(PI * y) * sin(PI * x) * sin(PI * z);
        break;
      case 1:
        val[0] = PI * cos(PI * z) * sin(PI * x) * sin(PI * y);
        val[1] = -PI * cos(PI * x) * cos(PI * y) * cos(PI * z);
        val[2] = -2 * PI * cos(PI * x) * sin(PI * y) * sin(PI * z);
        break;
      case 2:
        val[0] = PI * cos(PI * y) * sin(PI * x) * sin(PI * z);
        val[1] = PI * cos(PI * x) * sin(PI * y) * sin(PI * z);
        val[2] = 2 * PI * cos(PI * x) * cos(PI * y) * cos(PI * z);
        break;
    }
  return val;
}

void create_tria(Triangulation<3> &triangulation,
                 const Point<3> *  vertices_parallelograms)
{
  const std::vector<Point<3>> vertices(&vertices_parallelograms[0],
                                       &vertices_parallelograms[n_vertices]);

  // create grid with all possible combintations of face_flip, face_orientation
  // and face_rotation flags
  static const int cell_vertices[][GeometryInfo<3>::vertices_per_cell] = {
    {0, 1, 9, 10, 3, 4, 12, 13},  // cell 1 standard
    {1, 2, 10, 11, 4, 5, 13, 14}, // cell 2 standard
    //{10, 11, 13, 14, 1, 2, 4, 5},       // cell 2 rotated by 270 deg
    {9, 10, 18, 19, 12, 13, 21, 22},  // cell 3 standard
    {10, 11, 19, 20, 13, 14, 22, 23}, // cell 4 standard
    //{13, 14, 10, 11, 22, 23, 19, 20},   // cell 4 rotated by 90 deg
    {3, 4, 12, 13, 6, 7, 15, 16},     // cell 5 standard
    {4, 5, 13, 14, 7, 8, 16, 17},     // cell 6 standard
    {12, 13, 21, 22, 15, 16, 24, 25}, // cell 7 standard
    //{24, 25, 15, 16, 21, 22, 12, 13},   // cell 7 rotated by 180 deg
    {13, 14, 22, 23, 16, 17, 25, 26} // cell 8 standard
  };

  std::vector<CellData<3>> cells(n_cells, CellData<3>());
  for (unsigned i = 0; i < cells.size(); ++i)
    {
      for (const unsigned int j : GeometryInfo<3>::vertex_indices())
        cells[i].vertices[j] = cell_vertices[i][j];
      cells[i].material_id = 0;
    }

  triangulation.create_triangulation(vertices, cells, SubCellData());
}

template <int dim>
void
test(const FiniteElement<dim> &fe,
     unsigned                  n_cycles,
     bool                      global,
     const Point<dim> *        vertices_parallelograms)
{
  deallog << "dim: " << dim << "\t" << fe.get_name() << std::endl;
  deallog
    << "DoFs\t\t||u-u_h||_1\tcurl(u_h)\ttangentials\tcurl(curl(u_h))\tcurl_curl_traces\tdiv(u_h)\tboundary_flux"
    << std::endl;

  Triangulation<dim> triangulation;
  create_tria(triangulation, vertices_parallelograms);

  DoFHandler<dim> dof_handler(triangulation);

  VectorFunction<dim>              fe_function;
  const FEValuesExtractors::Vector vec(0);
  const QGauss<dim>                quadrature(fe.degree + 1);
  const QGauss<dim - 1>            face_quadrature(fe.degree + 1);
  const unsigned int               n_q_points      = quadrature.size();
  const unsigned int               n_face_q_points = face_quadrature.size();
  // MappingQ<dim> mapping(2);
  MappingQGeneric<dim>                                        mapping(1);
  std::vector<double>                                         div_v(n_q_points);
  std::vector<typename FEValuesViews::Vector<dim>::curl_type> curl_v(
    n_q_points);
  std::vector<Tensor<3, dim>> hessians(n_q_points);

  std::vector<Tensor<1, dim>> face_values(n_face_q_points);
  std::vector<typename FEValuesViews::Vector<dim>::curl_type> face_curls(
    n_face_q_points);

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
                                        VectorTools::L1_norm);

      typename FEValuesViews::Vector<dim>::curl_type total_curl,
        boundary_tangentials;
      Tensor<1, dim> total_curl_curl, boundary_curl_curl_traces;
      double         total_div     = 0;
      double         boundary_flux = 0;
      total_curl *= 0.;
      boundary_tangentials *= 0.;

      FEValues<dim>     fe_values(mapping,
                              fe,
                              quadrature,
                              update_JxW_values | update_quadrature_points |
                                update_values | update_gradients |
                                update_hessians);
      FEFaceValues<dim> fe_face_values(mapping,
                                       fe,
                                       face_quadrature,
                                       update_JxW_values |
                                         update_quadrature_points |
                                         update_values | update_gradients |
                                         update_normal_vectors);
      unsigned int      cell_index = 0;

      for (typename DoFHandler<dim>::active_cell_iterator cell =
             dof_handler.begin_active();
           cell != dof_handler.end();
           ++cell, ++cell_index)
        {
          fe_values.reinit(cell);
          const std::vector<double> &JxW_values = fe_values.get_JxW_values();
          fe_values[vec].get_function_divergences(v, div_v);
          fe_values[vec].get_function_curls(v, curl_v);
          fe_values[vec].get_function_hessians(v, hessians);
          for (unsigned int q_point = 0; q_point < n_q_points; ++q_point)
            {
              total_div += JxW_values[q_point] * div_v[q_point];
              total_curl += JxW_values[q_point] * curl_v[q_point];
              if (dim == 3)
                {
                  total_curl_curl[0] +=
                    JxW_values[q_point] *
                    (hessians[q_point][1][0][1] + hessians[q_point][2][0][2] -
                     hessians[q_point][0][1][1] - hessians[q_point][0][2][2]);
                  total_curl_curl[1] +=
                    JxW_values[q_point] *
                    (hessians[q_point][2][1][2] + hessians[q_point][0][0][1] -
                     hessians[q_point][1][2][2] - hessians[q_point][1][0][0]);
                  total_curl_curl[2] +=
                    JxW_values[q_point] *
                    (hessians[q_point][0][0][2] + hessians[q_point][1][1][2] -
                     hessians[q_point][2][0][0] - hessians[q_point][2][1][1]);
                }
            }

          for (const unsigned int face : GeometryInfo<dim>::face_indices())
            {
              fe_face_values.reinit(cell, face);
              const std::vector<double> &face_JxW_values =
                fe_face_values.get_JxW_values();
              fe_face_values[vec].get_function_values(v, face_values);
              if (dim == 3)
                fe_face_values[vec].get_function_curls(v, face_curls);
              for (unsigned int q_point = 0; q_point < n_face_q_points;
                   ++q_point)
                {
                  const Tensor<1, dim> &normal =
                    fe_face_values.normal_vector(q_point);

                  // boundary flux
                  if (cell->at_boundary(face))
                    boundary_flux += face_JxW_values[q_point] *
                                     (face_values[q_point] * normal);
                  else
                    total_div -= face_JxW_values[q_point] *
                                 (face_values[q_point] * normal);

                  // boundary tangentials (curl traces)
                  typename FEValuesViews::Vector<dim>::curl_type n_x_v;
                  if (dim == 2)
                    n_x_v[0] = (-normal[1] * face_values[q_point][0] +
                                normal[0] * face_values[q_point][1]);
                  else if (dim == 3)
                    n_x_v = cross_product_3d(normal, face_values[q_point]);

                  if (cell->at_boundary(face))
                    boundary_tangentials += face_JxW_values[q_point] * n_x_v;
                  else
                    total_curl -= face_JxW_values[q_point] * n_x_v;

                  // boundary curl curl traces
                  if (dim == 3)
                    {
                      Tensor<1, dim> n_x_curl_u =
                        cross_product_3d(normal,
                                         *reinterpret_cast<Tensor<1, dim> *>(
                                           &face_curls[q_point]));
                      if (cell->at_boundary(face))
                        boundary_curl_curl_traces +=
                          face_JxW_values[q_point] * n_x_curl_u;
                      else
                        total_curl_curl -=
                          face_JxW_values[q_point] * n_x_curl_u;
                    }
                }
            }
        }

      deallog << dof_handler.n_dofs() << "\t\t" << diff.l1_norm() << "\t"
              << total_curl.norm() << "\t" << boundary_tangentials.norm()
              << "\t" << total_curl_curl.norm() << "\t"
              << boundary_curl_curl_traces.norm() << "\t" << total_div << "\t"
              << boundary_flux << std::endl;

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

  const static unsigned dim      = 3;
  unsigned              order    = 1;
  unsigned              n_cycles = 2;

  deallog << "3d\nRectangular grid:\n";

  const Point<dim> *vertices = &vertices_rectangular[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);

  deallog << "\nAffine grid:\n";

  vertices = &vertices_affine[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);

  deallog << "\nNon-affine grid:\n";

  vertices = &vertices_nonaffine[0];
  test<dim>(FE_Nedelec<dim>(order), n_cycles, true, vertices);
  test<dim>(FE_RaviartThomas<dim>(order), n_cycles, true, vertices);

  deallog << std::endl;
}
