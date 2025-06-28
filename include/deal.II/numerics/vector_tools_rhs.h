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

#ifndef dealii_vector_tools_rhs_h
#define dealii_vector_tools_rhs_h

#include <deal.II/base/config.h>

#include <deal.II/base/template_constraints.h>
#include <deal.II/base/types.h>

#include <set>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;

template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;
template <int dim>
class Quadrature;
namespace hp
{
  template <int dim, int spacedim>
  class MappingCollection;
  template <int dim>
  class QCollection;
} // namespace hp


namespace VectorTools
{
  /**
   * @name Assembling of right hand sides
   */
  /** @{ */

  /**
   * Create a right hand side vector. Prior content of the given @p rhs_vector
   * vector is deleted.
   *
   * See the general documentation of this namespace for further information.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const Quadrature<dim>                                     &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * Call the create_right_hand_side() function, see above, with
   * <tt>mapping=MappingQ@<dim@>(1)</tt>.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof,
    const Quadrature<dim>                                     &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * Like the previous set of functions, but for hp-objects.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const hp::QCollection<dim>                                &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * Like the previous set of functions, but for hp-objects.
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof,
    const hp::QCollection<dim>                                &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const AffineConstraints<typename VectorType::value_type>  &constraints =
      AffineConstraints<typename VectorType::value_type>());

  /**
   * Create a right hand side vector from boundary forces. Prior content of
   * the given @p rhs_vector vector is deleted.
   *
   * See the general documentation of this namespace for further information.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const Quadrature<dim - 1>                                 &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids =
      std::set<types::boundary_id>());

  /**
   * Call the create_boundary_right_hand_side() function, see above, with
   * <tt>mapping=MappingQ@<dim@>(1)</tt>.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof,
    const Quadrature<dim - 1>                                 &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids =
      std::set<types::boundary_id>());

  /**
   * Same as the set of functions above, but for hp-objects.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const hp::MappingCollection<dim, spacedim>                &mapping,
    const DoFHandler<dim, spacedim>                           &dof,
    const hp::QCollection<dim - 1>                            &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids =
      std::set<types::boundary_id>());

  /**
   * Call the create_boundary_right_hand_side() function, see above, with a
   * single Q1 mapping as collection. This function therefore will only work
   * if the only active FE index in use is zero.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   *
   * @dealiiConceptRequires{concepts::is_writable_dealii_vector_type<VectorType>}
   */
  template <int dim, int spacedim, typename VectorType>
  DEAL_II_CXX20_REQUIRES(concepts::is_writable_dealii_vector_type<VectorType>)
  void create_boundary_right_hand_side(
    const DoFHandler<dim, spacedim>                           &dof,
    const hp::QCollection<dim - 1>                            &q,
    const Function<spacedim, typename VectorType::value_type> &rhs,
    VectorType                                                &rhs_vector,
    const std::set<types::boundary_id>                        &boundary_ids =
      std::set<types::boundary_id>());
  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_rhs_h
