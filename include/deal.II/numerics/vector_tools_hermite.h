// ---------------------------------------------------------------------
//
// Copyright (C) 2023 - 2023 by the deal.II authors
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

#ifndef dealii_vector_tools_hermite_h
#define dealii_vector_tools_hermite_h

#include <deal.II/base/config.h>

#include <deal.II/base/function.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/types.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/mapping.h>

#include <map>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace VectorTools
{
  /**
   * @name Domain and boundary projection for Hermite elements
   * @{
   */

  /**
   * Enforces boundary conditions by projecting onto the Hermite finite element
   * space at the boundary.
   *
   * This is required even in the case of Dirichlet boundary conditions where
   * similar functions already exist for other bases, because Hermite assigns
   * degrees of freedom directly onto element boundaries that have a zero shape
   * value across entire element faces. For most existing methods this can
   * create a singular mass matrix unless special steps are taken to remove
   * those degrees of freedom from the boundary system.
   *
   * Note that only Dirichlet, Neumann and second degree Neumann
   * conditions are currently supported.
   */
  template <int dim, int spacedim = dim, typename Number = double>
  void
  project_hermite_boundary_values(
    const Mapping<dim, spacedim>    &mapping_h,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &boundary_functions,
    const Quadrature<dim - 1>                 &quadrature,
    const unsigned int                         boundary_norm_deriv_order,
    std::map<types::global_dof_index, Number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});



  /**
   * The following function is a specialisation of the above for the case where
   * Dirichlet boundary conditions should be enforced.
   */
  template <int dim, int spacedim = dim, typename Number = double>
  void
  project_hermite_boundary_values(
    const Mapping<dim, spacedim>    &mapping_h,
    const DoFHandler<dim, spacedim> &dof_handler,
    const std::map<types::boundary_id, const Function<spacedim, Number> *>
                                              &boundary_functions,
    const Quadrature<dim - 1>                 &quadrature,
    std::map<types::global_dof_index, Number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});



  /**
   * Projection function for use with Hermite finite elements. Usually the
   * existing VectorTools::project() function will work with Hermite, however
   * setting any of the flags for enforcing boundary conditions to true would
   * cause errors, either by setting all derivatives to 0 at the boundary or
   * creating a singular boundary mass matrix. Since the Hermite mass matrix can
   * become poorly conditioned for higher polynomial degrees it is often useful
   * to perform separate boundary and domain projection, which means a separate
   * projection function for Hermite elements is necessary for their effective
   * use.
   */
  template <int dim, typename VectorType, int spacedim>
  void
  project_hermite(
    const Mapping<dim, spacedim>                              &mapping,
    const DoFHandler<dim, spacedim>                           &dofhandler,
    const AffineConstraints<typename VectorType::value_type>  &constraints,
    const Quadrature<dim>                                     &quadrature,
    const Function<spacedim, typename VectorType::value_type> &function,
    VectorType                                                &vec,
    const bool                 enforce_zero_boundary = false,
    const Quadrature<dim - 1> &q_boundary = (dim > 1 ? QGauss<dim - 1>(2) :
                                                       Quadrature<dim - 1>(0)),
    const bool                 project_to_boundary_first = false);

  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif
