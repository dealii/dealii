// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

// Test VectorTools::project_boundary_values_div_conforming convergence rates
// for the case that the DoFHandler constains more than one
// FE_RaviartThomas element.

#include <deal.II/base/function.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_raviart_thomas.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_cartesian.h>
#include <deal.II/fe/mapping_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>

#include <fstream>
#include <iostream>

#include "../tests.h"


template <int dim>
class BoundaryFunctionDisp : public Function<dim>
{
public:
  BoundaryFunctionDisp()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> & point,
        const unsigned int component = 0) const override
  {
    Assert(component < dim, ExcNotImplemented());
    return std::sin(10 * point[component] +
                    point[component == 0 ? 1 : component - 1]);
  }
};

template <int dim>
class BoundaryFunctionVelo : public Function<dim>
{
public:
  BoundaryFunctionVelo()
    : Function<dim>(dim)
  {}

  virtual double
  value(const Point<dim> & point,
        const unsigned int component = 0) const override
  {
    Assert(component < dim, ExcNotImplemented());
    return std::cos(10 * point[component] +
                    point[component == 0 ? 1 : component - 1]);
  }
};

template <int dim>
class BoundaryFunctionPres : public Function<dim>
{
public:
  BoundaryFunctionPres()
    : Function<dim>(2 * dim + 1)
  {}

  // To make interpolate_boundary_values happy, the pressure is assigned to
  // the 2*dim + 1th component
  virtual double
  value(const Point<dim> & point,
        const unsigned int component = 0) const override
  {
    if (component != 2 * dim)
      return std::numeric_limits<double>::quiet_NaN();
    else
      return std::sin(7.0 * point[0]) + std::cos(5.0 * point[dim - 1]);
  }
};


template <int dim>
void
test_boundary_values(const FiniteElement<dim> &fe)
{
  Triangulation<dim> triangulation;

  GridGenerator::hyper_cube(triangulation);
  MappingQ<dim>      mapping(2);
  QIterated<dim - 1> face_quadrature(QTrapez<1>(), 5);

  double old_max_disp_error = 0.0;
  double old_max_velo_error = 0.0;
  double old_max_pres_error = 0.0;
  for (unsigned int refinement_n = 0; refinement_n < 7 - dim; ++refinement_n)
    {
      triangulation.refine_global(1);

      DoFHandler<dim> dof_handler(triangulation);
      dof_handler.distribute_dofs(fe);
      const FEValuesExtractors::Vector displacements(0);
      const FEValuesExtractors::Vector velocities(dim);
      const FEValuesExtractors::Scalar pressure(2 * dim);

      BoundaryFunctionDisp<dim> boundary_function_disp;
      BoundaryFunctionVelo<dim> boundary_function_velo;
      BoundaryFunctionPres<dim> boundary_function_pres;

      AffineConstraints<double> constraints;

      std::map<types::boundary_id, const Function<dim> *> boundary_functions;
      boundary_functions[0] = &boundary_function_pres;
      VectorTools::interpolate_boundary_values(mapping,
                                               dof_handler,
                                               boundary_functions,
                                               constraints,
                                               fe.component_mask(pressure));
      VectorTools::project_boundary_values_div_conforming(
        dof_handler,
        0, /*first_vector_component*/
        boundary_function_disp,
        0, /*bdry_id*/
        constraints,
        mapping);
      VectorTools::project_boundary_values_div_conforming(
        dof_handler,
        dim, /*first_vector_component*/
        boundary_function_velo,
        0, /*bdry_id*/
        constraints,
        mapping);
      constraints.close();
      Vector<double> solution(dof_handler.n_dofs());
      constraints.distribute(solution);

      FEFaceValues<dim> face_values(mapping,
                                    fe,
                                    face_quadrature,
                                    update_values | update_normal_vectors |
                                      update_quadrature_points);
      double            max_disp_error = 0.0;
      double            max_velo_error = 0.0;
      double            max_pres_error = 0.0;

      std::vector<Tensor<1, dim>> cell_disp_values(face_quadrature.size());
      std::vector<Tensor<1, dim>> cell_velo_values(face_quadrature.size());
      std::vector<double>         cell_pres_values(face_quadrature.size());
      for (const auto &cell : dof_handler.active_cell_iterators())
        {
          for (const unsigned int face_n : GeometryInfo<dim>::face_indices())
            {
              auto face = cell->face(face_n);
              if (face->at_boundary())
                {
                  face_values.reinit(cell, face_n);
                  face_values[displacements].get_function_values(
                    solution, cell_disp_values);
                  face_values[velocities].get_function_values(solution,
                                                              cell_velo_values);
                  face_values[pressure].get_function_values(solution,
                                                            cell_pres_values);

                  for (unsigned int q_point_n = 0;
                       q_point_n < face_quadrature.size();
                       ++q_point_n)
                    {
                      const Point<dim> q_point =
                        face_values.quadrature_point(q_point_n);
                      const Tensor<1, dim> n_vector =
                        face_values.normal_vector(q_point_n);
                      Tensor<1, dim> exact_disp;
                      Tensor<1, dim> exact_velo;
                      for (unsigned int dim_n = 0; dim_n < dim; ++dim_n)
                        {
                          exact_disp[dim_n] =
                            boundary_function_disp.value(q_point, dim_n);
                          exact_velo[dim_n] =
                            boundary_function_velo.value(q_point, dim_n);
                        }
                      const double exact_pres =
                        boundary_function_pres.value(q_point, 2 * dim);

                      max_disp_error =
                        std::max(std::abs(cell_disp_values[q_point_n] *
                                            n_vector -
                                          exact_disp * n_vector),
                                 max_disp_error);
                      max_velo_error =
                        std::max(std::abs(cell_velo_values[q_point_n] *
                                            n_vector -
                                          exact_velo * n_vector),
                                 max_velo_error);
                      max_pres_error =
                        std::max(std::abs(cell_pres_values[q_point_n] -
                                          exact_pres),
                                 max_pres_error);
                    }
                }
            }
        }

      if (refinement_n != 0)
        {
          deallog << std::endl;
          deallog << "disp ratio: " << old_max_disp_error / max_disp_error
                  << std::endl;
          deallog << "velo ratio: " << old_max_velo_error / max_velo_error
                  << std::endl;
          deallog << "pres ratio: " << old_max_pres_error / max_pres_error
                  << std::endl;
        }

      old_max_disp_error = max_disp_error;
      old_max_velo_error = max_velo_error;
      old_max_pres_error = max_pres_error;
    }
}

int
main()
{
  initlog();
  {
    deallog << "2d:" << std::endl;
    constexpr unsigned int dim = 2;
    FE_RaviartThomas<dim>  u(1);
    FE_Q<dim>              p(1);
    FESystem<dim>          fesys(u, 2, p, 1);
    test_boundary_values(fesys);
  }

  {
    deallog << "3d:" << std::endl;
    constexpr unsigned int dim = 3;
    FE_RaviartThomas<dim>  u(1);
    FE_Q<dim>              p(1);
    FESystem<dim>          fesys(u, 2, p, 1);
    test_boundary_values(fesys);
  }
}
