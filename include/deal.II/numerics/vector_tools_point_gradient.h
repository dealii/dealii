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

#ifndef dealii_vector_tools_point_gradient_h
#define dealii_vector_tools_point_gradient_h


#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>

#include <vector>


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;

template <int dim, typename Number>
DEAL_II_CXX20_REQUIRES(dim >= 0)
class Point;

template <int rank_, int dim, typename Number>
class Tensor;
template <typename Number>
class Vector;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name Evaluation of functions and errors
   */
  /** @{ */

  /**
   * Evaluate a possibly vector-valued finite element function defined by the
   * given DoFHandler and nodal vector at the given point, and return the
   * (vector) gradient of this function through the last argument.
   *
   * This is a wrapper function using a Q1-mapping for cell boundaries to call
   * the other point_gradient() function.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  void point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  void point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * Evaluate a scalar finite element function defined by the given DoFHandler
   * and nodal vector at the given point, and return the gradient of this
   * function.
   *
   * Compared with the other function of the same name, this is a wrapper
   * function using a Q1-mapping for cells.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  Tensor<1, spacedim, typename VectorType::value_type> point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  Tensor<1, spacedim, typename VectorType::value_type> point_gradient(
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point);

  /**
   * Evaluate a possibly vector-valued finite element function defined by the
   * given DoFHandler and nodal vector at the given point, and return the
   * gradients of this function through the last argument.
   *
   * Compared with the other function of the same name, this function uses an
   * arbitrary mapping for evaluation.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  void point_gradient(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  void point_gradient(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const VectorType                           &fe_function,
    const Point<spacedim, double>              &point,
    std::vector<Tensor<1, spacedim, typename VectorType::value_type>> &value);

  /**
   * Evaluate a scalar finite element function defined by the given DoFHandler
   * and nodal vector at the given point, and return the gradient of this
   * function.
   *
   * Compared with the other function of the same name, this function uses an
   * arbitrary mapping for evaluation.
   *
   * This function is not particularly cheap. This is because it first
   * needs to find which cell a given point is in, then find the point
   * on the reference cell that matches the given evaluation point,
   * and then evaluate the shape functions there. You probably do not
   * want to use this function to evaluate the solution at <i>many</i>
   * points. For this kind of application, the FEFieldFunction class
   * offers at least some optimizations. On the other hand, if you
   * want to evaluate <i>many solutions</i> at the same point, you may
   * want to look at the VectorTools::create_point_source_vector()
   * function.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  Tensor<1, spacedim, typename VectorType::value_type> point_gradient(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const VectorType                &fe_function,
    const Point<spacedim, double>   &point);

  /**
   * Same as above for hp.
   *
   * @note If the cell in which the point is found is not locally owned, an
   * exception of type VectorTools::ExcPointNotAvailableHere is thrown.
   *
   * @note This function needs to find the cell within which a point lies,
   *   and this can only be done up to a certain numerical tolerance of course.
   *   Consequently, for points that are on, or close to, the boundary of
   *   a cell, you may get the gradient of the finite element field either
   *   here or there, depending on which cell the point is found in. Since
   *   the gradient is, for most elements, discontinuous from one cell or
   *   the other, you will get unpredictable values for
   *   points on or close to the boundary of the cell, as one would expect
   *   when trying to evaluate point values of discontinuous functions.
   *
   * @dealiiConceptRequires{concepts::is_dealii_vector_type<VectorType>}
   */
  template <int dim, typename VectorType, int spacedim>
  DEAL_II_CXX20_REQUIRES(concepts::is_dealii_vector_type<VectorType>)
  Tensor<1, spacedim, typename VectorType::value_type> point_gradient(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const VectorType                           &fe_function,
    const Point<spacedim, double>              &point);

  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_point_gradient_h
