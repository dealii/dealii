// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_tools_rhs_templates_h
#define dealii_vector_tools_rhs_templates_h

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools_rhs.h>

#include <set>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Quadrature<dim - 1>                                 &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids)
  {
    const FiniteElement<dim> &fe = dof_handler.get_fe();
    AssertDimension(fe.n_components(), rhs_function.n_components);
    AssertDimension(rhs_vector.size(), dof_handler.n_dofs());

    rhs_vector = typename VectorType::value_type(0.0);

    UpdateFlags update_flags =
      UpdateFlags(update_values | update_quadrature_points | update_JxW_values);
    FEFaceValues<dim> fe_values(mapping, fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points,
                       n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs(dofs_per_cell);
    Vector<double>                       cell_vector(dofs_per_cell);

    if (n_components == 1)
      {
        std::vector<double> rhs_values(n_q_points);

        for (const auto &cell : dof_handler.active_cell_iterators())
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->at_boundary() &&
                (boundary_ids.empty() ||
                 (boundary_ids.find(cell->face(face)->boundary_id()) !=
                  boundary_ids.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights = fe_values.get_JxW_values();
                rhs_function.value_list(fe_values.get_quadrature_points(),
                                        rhs_values);

                cell_vector = 0;
                for (unsigned int point = 0; point < n_q_points; ++point)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_vector(i) += rhs_values[point] *
                                      fe_values.shape_value(i, point) *
                                      weights[point];

                cell->get_dof_indices(dofs);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
    else
      {
        std::vector<Vector<double>> rhs_values(n_q_points,
                                               Vector<double>(n_components));

        for (const auto &cell : dof_handler.active_cell_iterators())
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->at_boundary() &&
                (boundary_ids.empty() ||
                 (boundary_ids.find(cell->face(face)->boundary_id()) !=
                  boundary_ids.end())))
              {
                fe_values.reinit(cell, face);

                const std::vector<double> &weights = fe_values.get_JxW_values();
                rhs_function.vector_value_list(
                  fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;

                // Use the faster code if the
                // FiniteElement is primitive
                if (fe.is_primitive())
                  {
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          const unsigned int component =
                            fe.system_to_component_index(i).first;

                          cell_vector(i) += rhs_values[point](component) *
                                            fe_values.shape_value(i, point) *
                                            weights[point];
                        }
                  }
                else
                  {
                    // And the full featured
                    // code, if vector valued
                    // FEs are used
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int comp_i = 0; comp_i < n_components;
                             ++comp_i)
                          if (fe.get_nonzero_components(i)[comp_i])
                            {
                              cell_vector(i) +=
                                rhs_values[point](comp_i) *
                                fe_values.shape_value_component(i,
                                                                point,
                                                                comp_i) *
                                weights[point];
                            }
                  }

                cell->get_dof_indices(dofs);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Quadrature<dim - 1>                                 &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids)
  {
    create_boundary_right_hand_side(StaticMappingQ1<dim>::mapping,
                                    dof_handler,
                                    quadrature,
                                    rhs_function,
                                    rhs_vector,
                                    boundary_ids);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const hp::QCollection<dim - 1>                            &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids)
  {
    const hp::FECollection<dim> &fe = dof_handler.get_fe_collection();
    AssertDimension(fe.n_components(), rhs_function.n_components);
    AssertDimension(rhs_vector.size(), dof_handler.n_dofs());

    rhs_vector = typename VectorType::value_type(0.0);

    UpdateFlags update_flags =
      UpdateFlags(update_values | update_quadrature_points | update_JxW_values);
    hp::FEFaceValues<dim> x_fe_values(mapping, fe, quadrature, update_flags);

    const unsigned int n_components = fe.n_components();

    std::vector<types::global_dof_index> dofs(fe.max_dofs_per_cell());
    Vector<double>                       cell_vector(fe.max_dofs_per_cell());

    if (n_components == 1)
      {
        std::vector<double> rhs_values;

        for (const auto &cell : dof_handler.active_cell_iterators())
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->at_boundary() &&
                (boundary_ids.empty() ||
                 (boundary_ids.find(cell->face(face)->boundary_id()) !=
                  boundary_ids.end())))
              {
                x_fe_values.reinit(cell, face);

                const FEFaceValues<dim> &fe_values =
                  x_fe_values.get_present_fe_values();

                const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                   n_q_points = fe_values.n_quadrature_points;
                rhs_values.resize(n_q_points);

                const std::vector<double> &weights = fe_values.get_JxW_values();
                rhs_function.value_list(fe_values.get_quadrature_points(),
                                        rhs_values);

                cell_vector = 0;
                for (unsigned int point = 0; point < n_q_points; ++point)
                  for (unsigned int i = 0; i < dofs_per_cell; ++i)
                    cell_vector(i) += rhs_values[point] *
                                      fe_values.shape_value(i, point) *
                                      weights[point];

                dofs.resize(dofs_per_cell);
                cell->get_dof_indices(dofs);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
    else
      {
        std::vector<Vector<double>> rhs_values;

        for (const auto &cell : dof_handler.active_cell_iterators())
          for (const unsigned int face : cell->face_indices())
            if (cell->face(face)->at_boundary() &&
                (boundary_ids.empty() ||
                 (boundary_ids.find(cell->face(face)->boundary_id()) !=
                  boundary_ids.end())))
              {
                x_fe_values.reinit(cell, face);

                const FEFaceValues<dim> &fe_values =
                  x_fe_values.get_present_fe_values();

                const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                   n_q_points = fe_values.n_quadrature_points;
                rhs_values.resize(n_q_points, Vector<double>(n_components));

                const std::vector<double> &weights = fe_values.get_JxW_values();
                rhs_function.vector_value_list(
                  fe_values.get_quadrature_points(), rhs_values);

                cell_vector = 0;

                // Use the faster code if the
                // FiniteElement is primitive
                if (cell->get_fe().is_primitive())
                  {
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        {
                          const unsigned int component =
                            cell->get_fe().system_to_component_index(i).first;

                          cell_vector(i) += rhs_values[point](component) *
                                            fe_values.shape_value(i, point) *
                                            weights[point];
                        }
                  }
                else
                  {
                    // And the full featured
                    // code, if vector valued
                    // FEs are used
                    for (unsigned int point = 0; point < n_q_points; ++point)
                      for (unsigned int i = 0; i < dofs_per_cell; ++i)
                        for (unsigned int comp_i = 0; comp_i < n_components;
                             ++comp_i)
                          if (cell->get_fe().get_nonzero_components(i)[comp_i])
                            {
                              cell_vector(i) +=
                                rhs_values[point](comp_i) *
                                fe_values.shape_value_component(i,
                                                                point,
                                                                comp_i) *
                                weights[point];
                            }
                  }
                dofs.resize(dofs_per_cell);
                cell->get_dof_indices(dofs);

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  rhs_vector(dofs[i]) += cell_vector(i);
              }
      }
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const hp::QCollection<dim - 1>                            &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids)
  {
    create_boundary_right_hand_side(
      hp::StaticMappingQ1<dim>::mapping_collection,
      dof_handler,
      quadrature,
      rhs_function,
      rhs_vector,
      boundary_ids);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints)
  {
    using Number = typename VectorType::value_type;

    const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();
    AssertDimension(fe.n_components(), rhs_function.n_components);
    AssertDimension(rhs_vector.size(), dof_handler.n_dofs());
    rhs_vector = typename VectorType::value_type(0.0);

    UpdateFlags update_flags =
      UpdateFlags(update_values | update_quadrature_points | update_JxW_values);
    FEValues<dim, spacedim> fe_values(mapping, fe, quadrature, update_flags);

    const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                       n_q_points    = fe_values.n_quadrature_points,
                       n_components  = fe.n_components();

    std::vector<types::global_dof_index> dofs(dofs_per_cell);
    Vector<Number>                       cell_vector(dofs_per_cell);

    if (n_components == 1)
      {
        std::vector<Number> rhs_values(n_q_points);

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);

              const std::vector<double> &weights = fe_values.get_JxW_values();
              rhs_function.value_list(fe_values.get_quadrature_points(),
                                      rhs_values);

              cell_vector = 0;
              for (unsigned int point = 0; point < n_q_points; ++point)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_vector(i) += rhs_values[point] *
                                    fe_values.shape_value(i, point) *
                                    weights[point];

              cell->get_dof_indices(dofs);

              constraints.distribute_local_to_global(cell_vector,
                                                     dofs,
                                                     rhs_vector);
            }
      }
    else
      {
        std::vector<Vector<Number>> rhs_values(n_q_points,
                                               Vector<Number>(n_components));

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              fe_values.reinit(cell);

              const std::vector<double> &weights = fe_values.get_JxW_values();
              rhs_function.vector_value_list(fe_values.get_quadrature_points(),
                                             rhs_values);

              cell_vector = 0;
              // Use the faster code if the
              // FiniteElement is primitive
              if (fe.is_primitive())
                {
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      {
                        const unsigned int component =
                          fe.system_to_component_index(i).first;

                        cell_vector(i) += rhs_values[point](component) *
                                          fe_values.shape_value(i, point) *
                                          weights[point];
                      }
                }
              else
                {
                  // Otherwise do it the way
                  // proposed for vector valued
                  // elements
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      for (unsigned int comp_i = 0; comp_i < n_components;
                           ++comp_i)
                        if (fe.get_nonzero_components(i)[comp_i])
                          {
                            cell_vector(i) +=
                              rhs_values[point](comp_i) *
                              fe_values.shape_value_component(i,
                                                              point,
                                                              comp_i) *
                              weights[point];
                          }
                }
              cell->get_dof_indices(dofs);

              constraints.distribute_local_to_global(cell_vector,
                                                     dofs,
                                                     rhs_vector);
            }
      }

    rhs_vector.compress(VectorOperation::values::add);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints)
  {
    create_right_hand_side(get_default_linear_mapping(
                             dof_handler.get_triangulation()),
                           dof_handler,
                           quadrature,
                           rhs_function,
                           rhs_vector,
                           constraints);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const hp::QCollection<dim>                                &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints)
  {
    using Number = typename VectorType::value_type;

    const hp::FECollection<dim, spacedim> &fe = dof_handler.get_fe_collection();
    AssertDimension(fe.n_components(), rhs_function.n_components);
    AssertDimension(rhs_vector.size(), dof_handler.n_dofs());
    rhs_vector = typename VectorType::value_type(0.0);

    UpdateFlags update_flags =
      UpdateFlags(update_values | update_quadrature_points | update_JxW_values);
    hp::FEValues<dim, spacedim> x_fe_values(mapping,
                                            fe,
                                            quadrature,
                                            update_flags);

    const unsigned int n_components = fe.n_components();

    std::vector<types::global_dof_index> dofs(fe.max_dofs_per_cell());
    Vector<Number>                       cell_vector(fe.max_dofs_per_cell());

    if (n_components == 1)
      {
        std::vector<Number> rhs_values;

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              x_fe_values.reinit(cell);

              const FEValues<dim, spacedim> &fe_values =
                x_fe_values.get_present_fe_values();

              const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                 n_q_points    = fe_values.n_quadrature_points;
              rhs_values.resize(n_q_points);
              dofs.resize(dofs_per_cell);
              cell_vector.reinit(dofs_per_cell);

              const auto &weights = fe_values.get_JxW_values();
              rhs_function.value_list(fe_values.get_quadrature_points(),
                                      rhs_values);

              cell_vector = 0;
              for (unsigned int point = 0; point < n_q_points; ++point)
                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                  cell_vector(i) += rhs_values[point] *
                                    fe_values.shape_value(i, point) *
                                    weights[point];

              cell->get_dof_indices(dofs);

              constraints.distribute_local_to_global(cell_vector,
                                                     dofs,
                                                     rhs_vector);
            }
      }
    else
      {
        std::vector<Vector<Number>> rhs_values;

        for (const auto &cell : dof_handler.active_cell_iterators())
          if (cell->is_locally_owned())
            {
              x_fe_values.reinit(cell);

              const FEValues<dim, spacedim> &fe_values =
                x_fe_values.get_present_fe_values();

              const unsigned int dofs_per_cell = fe_values.dofs_per_cell,
                                 n_q_points    = fe_values.n_quadrature_points;
              rhs_values.resize(n_q_points, Vector<Number>(n_components));
              dofs.resize(dofs_per_cell);
              cell_vector.reinit(dofs_per_cell);

              const auto &weights = fe_values.get_JxW_values();
              rhs_function.vector_value_list(fe_values.get_quadrature_points(),
                                             rhs_values);

              cell_vector = 0;

              // Use the faster code if the
              // FiniteElement is primitive
              if (cell->get_fe().is_primitive())
                {
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      {
                        const unsigned int component =
                          cell->get_fe().system_to_component_index(i).first;

                        cell_vector(i) += rhs_values[point](component) *
                                          fe_values.shape_value(i, point) *
                                          weights[point];
                      }
                }
              else
                {
                  // Otherwise do it the way proposed
                  // for vector valued elements
                  for (unsigned int point = 0; point < n_q_points; ++point)
                    for (unsigned int i = 0; i < dofs_per_cell; ++i)
                      for (unsigned int comp_i = 0; comp_i < n_components;
                           ++comp_i)
                        if (cell->get_fe().get_nonzero_components(i)[comp_i])
                          {
                            cell_vector(i) +=
                              rhs_values[point](comp_i) *
                              fe_values.shape_value_component(i,
                                                              point,
                                                              comp_i) *
                              weights[point];
                          }
                }

              cell->get_dof_indices(dofs);

              constraints.distribute_local_to_global(cell_vector,
                                                     dofs,
                                                     rhs_vector);
            }
      }

    rhs_vector.compress(VectorOperation::values::add);
  }



  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof_handler,
    const hp::QCollection<dim>                                &quadrature,
    const Function<spacedim, typename VectorType::value_type> &rhs_function,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints)
  {
    create_right_hand_side(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof_handler,
      quadrature,
      rhs_function,
      rhs_vector,
      constraints);
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_rhs_templates_h
