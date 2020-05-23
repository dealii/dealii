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

#ifndef dealii_vector_tools_point_gradient_templates_h
#define dealii_vector_tools_point_gradient_templates_h


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
#include <deal.II/numerics/vector_tools_point_gradient.h>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType &               fe_function,
    const Point<spacedim> &          point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>>
      &gradients)
  {
    point_gradient(StaticMappingQ1<dim, spacedim>::mapping,
                   dof,
                   fe_function,
                   point,
                   gradients);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const hp::DoFHandler<dim, spacedim> &dof,
    const VectorType &                   fe_function,
    const Point<spacedim> &              point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>>
      &gradients)
  {
    point_gradient(hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
                   dof,
                   fe_function,
                   point,
                   gradients);
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const DoFHandler<dim, spacedim> &dof,
                 const VectorType &               fe_function,
                 const Point<spacedim> &          point)
  {
    return point_gradient(StaticMappingQ1<dim, spacedim>::mapping,
                          dof,
                          fe_function,
                          point);
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const hp::DoFHandler<dim, spacedim> &dof,
                 const VectorType &                   fe_function,
                 const Point<spacedim> &              point)
  {
    return point_gradient(
      hp::StaticMappingQ1<dim, spacedim>::mapping_collection,
      dof,
      fe_function,
      point);
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const VectorType &               fe_function,
    const Point<spacedim> &          point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &gradient)
  {
    const FiniteElement<dim> &fe = dof.get_fe();

    Assert(gradient.size() == fe.n_components(),
           ExcDimensionMismatch(gradient.size(), fe.n_components()));

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

    FEValues<dim> fe_values(mapping, fe, quadrature, update_gradients);
    fe_values.reinit(cell_point.first);

    // then use this to get the gradients of
    // the given fe_function at this point
    using Number = typename VectorType::value_type;
    std::vector<std::vector<Tensor<1, dim, Number>>> u_gradient(
      1, std::vector<Tensor<1, dim, Number>>(fe.n_components()));
    fe_values.get_function_gradients(fe_function, u_gradient);

    gradient = u_gradient[0];
  }


  template <int dim, typename VectorType, int spacedim>
  void
  point_gradient(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof,
    const VectorType &                          fe_function,
    const Point<spacedim> &                     point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &gradient)
  {
    using Number                              = typename VectorType::value_type;
    const hp::FECollection<dim, spacedim> &fe = dof.get_fe_collection();

    Assert(gradient.size() == fe.n_components(),
           ExcDimensionMismatch(gradient.size(), fe.n_components()));

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
                                             update_gradients);
    hp_fe_values.reinit(cell_point.first);
    const FEValues<dim, spacedim> &fe_values =
      hp_fe_values.get_present_fe_values();

    std::vector<std::vector<Tensor<1, dim, Number>>> u_gradient(
      1, std::vector<Tensor<1, dim, Number>>(fe.n_components()));
    fe_values.get_function_gradients(fe_function, u_gradient);

    gradient = u_gradient[0];
  }


  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const Mapping<dim, spacedim> &   mapping,
                 const DoFHandler<dim, spacedim> &dof,
                 const VectorType &               fe_function,
                 const Point<spacedim> &          point)
  {
    Assert(dof.get_fe(0).n_components() == 1,
           ExcMessage(
             "Finite element is not scalar as is necessary for this function"));

    std::vector<Tensor<1, dim, typename VectorType::value_type>> gradient(1);
    point_gradient(mapping, dof, fe_function, point, gradient);

    return gradient[0];
  }



  template <int dim, typename VectorType, int spacedim>
  Tensor<1, spacedim, typename VectorType::value_type>
  point_gradient(const hp::MappingCollection<dim, spacedim> &mapping,
                 const hp::DoFHandler<dim, spacedim> &       dof,
                 const VectorType &                          fe_function,
                 const Point<spacedim> &                     point)
  {
    Assert(dof.get_fe(0).n_components() == 1,
           ExcMessage(
             "Finite element is not scalar as is necessary for this function"));

    std::vector<Tensor<1, dim, typename VectorType::value_type>> gradient(1);
    point_gradient(mapping, dof, fe_function, point, gradient);

    return gradient[0];
  }
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_gradient_templates_h
