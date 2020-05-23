// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2020 by the deal.II authors
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

#ifndef dealii_vector_tools_point_value_templates_h
#define dealii_vector_tools_point_value_templates_h


#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_tools.h>

#include <deal.II/hp/fe_values.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <deal.II/numerics/vector_tools_common.h>
#include <deal.II/numerics/vector_tools_point_value.h>


DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim> &                  point,
              Vector<typename VectorType::value_type> &value)
  {
    point_value(
      StaticMappingQ1<dim, spacedim>::mapping, dof, fe_function, point, value);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const hp::DoFHandler<dim, spacedim> &    dof,
              const VectorType &                       fe_function,
              const Point<spacedim> &                  point,
              Vector<typename VectorType::value_type> &value)
  {
    point_value(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
                dof,
                fe_function,
                point,
                value);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim> &          point)
  {
    return point_value(StaticMappingQ1<dim, spacedim>::mapping,
                       dof,
                       fe_function,
                       point);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const hp::DoFHandler<dim, spacedim> &dof,
              const VectorType &                   fe_function,
              const Point<spacedim> &              point)
  {
    return point_value(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
                       dof,
                       fe_function,
                       point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const Mapping<dim, spacedim> &           mapping,
              const DoFHandler<dim, spacedim> &        dof,
              const VectorType &                       fe_function,
              const Point<spacedim> &                  point,
              Vector<typename VectorType::value_type> &value)
  {
    using Number                 = typename VectorType::value_type;
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(value.size() == fe.n_components(),
           ExcDimensionMismatch(value.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
                    Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim> quadrature(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell_point.first);

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<Number>> u_value(1, Vector<Number>(fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const hp::DoFHandler<dim, spacedim> &       dof,
              const VectorType &                          fe_function,
              const Point<spacedim> &                     point,
              Vector<typename VectorType::value_type> &   value)
  {
    using Number                              = typename VectorType::value_type;
    const hp::FECollection<dim, spacedim> &fe = dof.get_fe_collection();

    Assert(value.size() == fe.n_components(),
           ExcDimensionMismatch(value.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<
      typename hp::DoFHandler<dim, spacedim>::active_cell_iterator,
      Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim> quadrature(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    hp::FEValues<dim, spacedim> hp_fe_values(mapping,
                                             fe,
                                             hp::QCollection<dim>(quadrature),
                                             update_values);
    hp_fe_values.reinit(cell_point.first);
    const FEValues<dim, spacedim> &fe_values =
      hp_fe_values.get_present_fe_values();

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<Number>> u_value(1, Vector<Number>(fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    value = u_value[0];
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const Mapping<dim, spacedim> &   mapping,
              const DoFHandler<dim, spacedim> &dof,
              const VectorType &               fe_function,
              const Point<spacedim> &          point)
  {
    Assert(dof.get_fe(0).n_components() == 1,
           ExcMessage(
             "Finite element is not scalar as is necessary for this function"));

    Vector<typename VectorType::value_type> value(1);
    point_value(mapping, dof, fe_function, point, value);

    return value(0);
  }


  template <int dim, typename VectorType, int spacedim>
  typename VectorType::value_type
  point_value(const hp::MappingCollection<dim, spacedim> &mapping,
              const hp::DoFHandler<dim, spacedim> &       dof,
              const VectorType &                          fe_function,
              const Point<spacedim> &                     point)
  {
    Assert(dof.get_fe(0).n_components() == 1,
           ExcMessage(
             "Finite element is not scalar as is necessary for this function"));

    Vector<typename VectorType::value_type> value(1);
    point_value(mapping, dof, fe_function, point, value);

    return value(0);
  }

  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_function,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim> &                                    point)
  {
    point_difference(StaticMappingQ1<dim>::mapping,
                     dof,
                     fe_function,
                     exact_function,
                     difference,
                     point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_difference(
    const Mapping<dim, spacedim> &                             mapping,
    const DoFHandler<dim, spacedim> &                          dof,
    const VectorType &                                         fe_function,
    const Function<spacedim, typename VectorType::value_type> &exact_function,
    Vector<typename VectorType::value_type> &                  difference,
    const Point<spacedim> &                                    point)
  {
    using Number                 = typename VectorType::value_type;
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(difference.size() == fe.n_components(),
           ExcDimensionMismatch(difference.size(), fe.n_components()));

    // first find the cell in which this point
    // is, initialize a quadrature rule with
    // it, and then a FEValues object
    const std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
                    Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof, point);

    AssertThrow(cell_point.first->is_locally_owned(),
                ExcPointNotAvailableHere());
    Assert(GeometryInfo<dim>::distance_to_unit_cell(cell_point.second) < 1e-10,
           ExcInternalError());

    const Quadrature<dim> quadrature(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));
    FEValues<dim> fe_values(mapping, fe, quadrature, update_values);
    fe_values.reinit(cell_point.first);

    // then use this to get at the values of
    // the given fe_function at this point
    std::vector<Vector<Number>> u_value(1, Vector<Number>(fe.n_components()));
    fe_values.get_function_values(fe_function, u_value);

    if (fe.n_components() == 1)
      difference(0) = exact_function.value(point);
    else
      exact_function.vector_value(point, difference);

    for (unsigned int i = 0; i < difference.size(); ++i)
      difference(i) -= u_value[0](i);
  }

  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &          p,
                             Vector<double> &                 rhs_vector)
  {
    Assert(rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert(dof_handler.get_fe(0).n_components() == 1,
           ExcMessage("This function only works for scalar finite elements"));

    rhs_vector = 0;

    std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
              Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof_handler, p);

    Quadrature<dim> q(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim, spacedim> fe_values(mapping,
                                      dof_handler.get_fe(),
                                      q,
                                      UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell_point.first->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) = fe_values.shape_value(i, 0);
  }



  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &          p,
                             Vector<double> &                 rhs_vector)
  {
    create_point_source_vector(StaticMappingQ1<dim, spacedim>::mapping,
                               dof_handler,
                               p,
                               rhs_vector);
  }


  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof_handler,
    const Point<spacedim> &                     p,
    Vector<double> &                            rhs_vector)
  {
    Assert(rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert(dof_handler.get_fe(0).n_components() == 1,
           ExcMessage("This function only works for scalar finite elements"));

    rhs_vector = 0;

    std::pair<typename hp::DoFHandler<dim, spacedim>::active_cell_iterator,
              Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof_handler, p);

    Quadrature<dim> q(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
                            cell_point.first->get_fe(),
                            q,
                            UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell_point.first->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) = fe_values.shape_value(i, 0);
  }



  template <int dim, int spacedim>
  void
  create_point_source_vector(const hp::DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &              p,
                             Vector<double> &                     rhs_vector)
  {
    create_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
                               dof_handler,
                               p,
                               rhs_vector);
  }



  template <int dim, int spacedim>
  void
  create_point_source_vector(const Mapping<dim, spacedim> &   mapping,
                             const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &          p,
                             const Point<dim> &               orientation,
                             Vector<double> &                 rhs_vector)
  {
    Assert(rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert(dof_handler.get_fe(0).n_components() == dim,
           ExcMessage(
             "This function only works for vector-valued finite elements."));

    rhs_vector = 0;

    const std::pair<typename DoFHandler<dim, spacedim>::active_cell_iterator,
                    Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof_handler, p);

    const Quadrature<dim> q(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    const FEValuesExtractors::Vector vec(0);
    FEValues<dim, spacedim>          fe_values(mapping,
                                      dof_handler.get_fe(),
                                      q,
                                      UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell_point.first->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =
        orientation * fe_values[vec].value(i, 0);
  }



  template <int dim, int spacedim>
  void
  create_point_source_vector(const DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &          p,
                             const Point<dim> &               orientation,
                             Vector<double> &                 rhs_vector)
  {
    create_point_source_vector(StaticMappingQ1<dim, spacedim>::mapping,
                               dof_handler,
                               p,
                               orientation,
                               rhs_vector);
  }


  template <int dim, int spacedim>
  void
  create_point_source_vector(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof_handler,
    const Point<spacedim> &                     p,
    const Point<dim> &                          orientation,
    Vector<double> &                            rhs_vector)
  {
    Assert(rhs_vector.size() == dof_handler.n_dofs(),
           ExcDimensionMismatch(rhs_vector.size(), dof_handler.n_dofs()));
    Assert(dof_handler.get_fe(0).n_components() == dim,
           ExcMessage(
             "This function only works for vector-valued finite elements."));

    rhs_vector = 0;

    std::pair<typename hp::DoFHandler<dim, spacedim>::active_cell_iterator,
              Point<spacedim>>
      cell_point =
        GridTools::find_active_cell_around_point(mapping, dof_handler, p);

    Quadrature<dim> q(
      GeometryInfo<dim>::project_to_unit_cell(cell_point.second));

    const FEValuesExtractors::Vector vec(0);
    FEValues<dim> fe_values(mapping[cell_point.first->active_fe_index()],
                            cell_point.first->get_fe(),
                            q,
                            UpdateFlags(update_values));
    fe_values.reinit(cell_point.first);

    const unsigned int dofs_per_cell = cell_point.first->get_fe().dofs_per_cell;

    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
    cell_point.first->get_dof_indices(local_dof_indices);

    for (unsigned int i = 0; i < dofs_per_cell; i++)
      rhs_vector(local_dof_indices[i]) =
        orientation * fe_values[vec].value(i, 0);
  }



  template <int dim, int spacedim>
  void
  create_point_source_vector(const hp::DoFHandler<dim, spacedim> &dof_handler,
                             const Point<spacedim> &              p,
                             const Point<dim> &                   orientation,
                             Vector<double> &                     rhs_vector)
  {
    create_point_source_vector(hp::StaticMappingQ1<dim>::mapping_collection,
                               dof_handler,
                               p,
                               orientation,
                               rhs_vector);
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_value_templates_h
