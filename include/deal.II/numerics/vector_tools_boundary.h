// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_vector_tools_boundary_h
#define dealii_vector_tools_boundary_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/hp/mapping_collection.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
template <int dim, typename Number>
class Function;
namespace hp
{
  template <int dim>
  class QCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name Interpolation and projection
   */
  /** @{ */

  /**
   * Compute constraints on the solution that corresponds to the imposition
   * of Dirichlet boundary conditions on parts of the boundary.
   * This function creates a map of
   * degrees of freedom subject to Dirichlet boundary conditions and the
   * corresponding values to be assigned to them, by interpolation around the
   * boundary. For each degree of freedom at the boundary, its boundary value
   * will be overwritten if its index already exists in @p boundary_values.
   * Otherwise, a new entry with proper index and boundary value for this
   * degree of freedom will be inserted into @p boundary_values.
   *
   * The parameter @p function_map provides a list of boundary indicators to
   * be handled by this function and corresponding boundary value functions.
   * The key of this map corresponds to the number @p boundary_id of the face.
   * numbers::internal_face_boundary_id is an illegal value for this key since
   * it is reserved for interior faces. For an example of how to use this
   * argument with a non-empty map, see the step-16 tutorial program.
   *
   * The flags in the last parameter, @p component_mask, denote which
   * components of the finite element space shall be interpolated. If it is
   * left as specified by the default value (i.e. an empty array), all
   * components are interpolated. If it is different from the default value,
   * it is assumed that the number of entries equals the number of components
   * in the boundary functions and the finite element, and those components in
   * the given boundary function will be used for which the respective flag
   * was set in the component mask. See also
   * @ref GlossComponentMask.
   * As an example, assume that you are solving the Stokes equations in 2d,
   * with variables $(u,v,p)$ and that you only want to interpolate boundary
   * values for the velocity, then the component mask should correspond to
   * <code>(true,true,false)</code>.
   *
   * @note Whether a component mask has been specified or not, the number of
   * components of the functions in @p function_map must match that of the
   * finite element used by @p dof. In other words, for the example above, you
   * need to provide a Function object that has 3 components (the two
   * velocities and the pressure), even though you are only interested in the
   * first two of them. interpolate_boundary_values() will then call this
   * function to obtain a vector of 3 values at each interpolation point but
   * only take the first two and discard the third. In other words, you are
   * free to return whatever you like in the third component of the vector
   * returned by Function::vector_value, but the Function object must state
   * that it has 3 components.
   *
   * If the finite element used has shape functions that are non-zero in more
   * than one component (in deal.II speak: they are non-primitive), then these
   * components can presently not be used for interpolating boundary values.
   * Thus, the elements in the component mask corresponding to the components
   * of these non-primitive shape functions must be @p false.
   *
   * @note This function applies the same ComponentMask to all parts of
   * the boundary indicated by keys of the `function_map` argument.
   * If you want to apply different component masks to parts of the boundary
   * represented by different boundary indicators, this function needs to be
   * called multiple times. For performance reasons, it might be reasonable to
   * use the present function by grouping together all boundary indicators with
   * the same ComponentMask. An alternative is to use one of the other
   * functions with this name, which take only one boundary indicator with
   * corresponding boundary function, to be called separately for every
   * boundary indicator.
   *
   * @note Mathematically, boundary conditions can only be applied to a
   *   part of the boundary that has a nonzero $(d-1)$-dimensional measure;
   *   in other words, it must be the union of *faces* of a mesh, rather than
   *   a set of edges in 3d, or even just a few vertices. That is because
   *   applying boundary conditions on individual vertices (rather than
   *   on the entire face of which this vertex might be a part) would
   *   correspond to using Dirac delta functions as boundary values, and
   *   this generally leads to singular solutions that can not adequately
   *   be resolved with the finite element method. These considerations
   *   notwithstanding, people often do apply boundary conditions at individual
   *   vertices -- in particular in solid mechanics, where one would then
   *   impose constraints on one or all components of the displacement at a
   *   vertex. This function does not support this operation: It works solely
   *   by looping over faces, checking whether the boundary indicator of the
   *   face is one of the ones of interest, and then considers all of the
   *   degrees of freedom on the face; it does not consider vertices (or,
   *   in 3d, edges) separately from the faces they are part of. But you can
   *   impose constraints on individual vertices by looping over all cells,
   *   over all vertices of each cell, and identifying whether this is the
   *   vertex you care about; then you use DoFAccessor::vertex_dof_index()
   *   to obtain the indices of the DoFs located on it. You can then
   *   entries for these degrees of freedom by hand to the `std::map`
   *   or AffineConstraints object you typically use to represent
   *   boundary value constraints.
   *
   * @note When solving a partial differential equation with boundary
   *   conditions $u|_{\partial\Omega}=g$ (on the entire boundary
   *   $\partial\Omega$, or perhaps only on parts $\Gamma\subset\partial\Omega$
   *   of the boundary),
   *   then this boundary condition is in general not satisfiable exactly
   *   using finite elements in the form $u_h|_{\partial\Omega}=g$. That is
   *   because the function $g$ is generally not a polynomial, whereas
   *   $u_h|_{\partial\Omega}$ *is* a polynomial on each face of the
   *   mesh that is located at the boundary. In other words, it is in
   *   general not possible to *impose* such boundary condition; what one
   *   *can* do, however, is to impose
   *     @f[ u_h|_{\partial\Omega}=I_h^{\partial\Omega} g, @f]
   *   where $I_h^{\partial\Omega} g$ is a function that equals $g$ at each node
   *   of the finite element space located on the boundary, and is piecewise
   *   polynomial in between. In other words, $I_h^{\partial\Omega}$ is an
   *   *interpolation operator* and $I_h^{\partial\Omega} g$ are the
   *   interpolated boundary values -- thus the name. The use of
   *   $I_h^{\partial\Omega} g$ instead of $g$ as boundary values imposes
   *   an additional error (in the same spirit as using quadrature introduces
   *   an additional error compared to being able to compute the integrals of
   *   the weak form exactly). In most cases, this additional error is of the
   *   same order as the other error terms in the finite element method,
   *   though there are some subtle differences when measuring the error in
   *   the $L^2$ norm. For some details, see @cite Bartels2004 .
   *
   * @note An alternative to using the interpolant,
   *     @f[ u_h|_{\partial\Omega}=I_h^{\partial\Omega} g @f]
   *   is to use the *projection* of the boundary values $g$ onto the
   *   finite element space on the boundary:
   *     @f[ u_h|_{\partial\Omega}=\Pi_h^{\partial\Omega} g. @f]
   *   The projection is available using the project_boundary_values()
   *   function. Using the projection may have some theoretical advantages
   *   (see again @cite Bartels2004) but has the practical disadvantage
   *   that computing the projection is far more expensive than computing
   *   the interpolation because the latter can be done one face at a time
   *   whereas the projection requires the solution of a problem on the entire
   *   boundary. On the other hand, interpolation is only possible for
   *   "nodal" finite element spaces (such as FE_Q, but not
   *   FE_Q_Hierarchical), whereas the projection is always possible.
   *
   * See the general documentation of this namespace for more information.
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask                       &component_mask = {});

  /**
   * Like the previous function, but take a mapping collection to go with
   * DoFHandler objects with hp-capabilities.
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask                       &component_mask = {});

  /**
   * Like the previous functions but without Mapping argument, using
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt> internally.
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask                       &component_mask = {});

  /**
   * Take only one boundary indicator with corresponding boundary function.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim>              &mapping,
    const DoFHandler<dim, spacedim>           &dof,
    const types::boundary_id                   boundary_indicator,
    const Function<spacedim, number>          &boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask                       &component_mask = {});

  /**
   * Like the previous function, but take a mapping collection to go with
   * DoFHandler objects with hp-capabilities.
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const types::boundary_id                    boundary_indicator,
    const Function<spacedim, number>           &boundary_function,
    std::map<types::global_dof_index, number>  &boundary_values,
    const ComponentMask                        &component_mask = {});

  /**
   * Like the previous functions but without Mapping argument, using
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt> internally.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim>           &dof,
    const types::boundary_id                   boundary_indicator,
    const Function<spacedim, number>          &boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask                       &component_mask = {});


  /**
   * Insert the (algebraic) constraints due to Dirichlet boundary conditions
   * into an AffineConstraints object. This function identifies the
   * degrees of freedom subject to Dirichlet boundary conditions, adds them to
   * the list of constrained DoFs in @p constraints and sets the respective
   * inhomogeneity to the value interpolated around the boundary. If this
   * routine encounters a DoF that already is constrained (for instance by a
   * hanging node constraint, see below, or any other type of constraint, e.g.
   * from periodic boundary conditions), the old setting of the constraint is
   * kept and nothing happens.
   *
   * @note When combining adaptively refined meshes with hanging node
   * constraints and boundary conditions like from the current function within
   * one AffineConstraints object, the hanging node constraints should always
   * be set first and then the boundary conditions, since boundary conditions
   * are not set in the second operation on degrees of freedom that are
   * already constrained. This makes sure that the discretization remains
   * conforming as is needed. See the discussion on conflicting constraints in
   * the topic on @ref constraints.
   *
   * For further information and details on the other function arguments, see
   * the interpolate_boundary_values() function with `std::map` arguments and
   * the general class documentation.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask       &component_mask = {});

  /**
   * Like the previous function, but take a mapping collection to go with
   * DoFHandler objects with hp-capabilities.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask       &component_mask = {});

  /**
   * Like the previous functions but without Mapping argument, using
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt> internally.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask       &component_mask = {});

  /**
   * Take only one boundary indicator with corresponding boundary function.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim>     &mapping,
    const DoFHandler<dim, spacedim>  &dof,
    const types::boundary_id          boundary_indicator,
    const Function<spacedim, number> &boundary_function,
    AffineConstraints<number>        &constraints,
    const ComponentMask              &component_mask = {});

  /**
   * Like the previous function, but take a mapping collection to go with
   * DoFHandler objects with hp-capabilities.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const types::boundary_id                    boundary_indicator,
    const Function<spacedim, number>           &boundary_function,
    AffineConstraints<number>                  &constraints,
    const ComponentMask                        &component_mask = {});

  /**
   * Like the previous functions but without Mapping argument, using
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt> internally.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const DoFHandler<dim, spacedim>  &dof,
    const types::boundary_id          boundary_indicator,
    const Function<spacedim, number> &boundary_function,
    AffineConstraints<number>        &constraints,
    const ComponentMask              &component_mask = {});


  /**
   * Project a function or a set of functions to the boundary of the domain.
   * In other words, compute the solution of the following problem: Find $u_h
   * \in V_h$ (where $V_h$ is the finite element space represented by the
   * DoFHandler argument of this function) so that
   * @f{align*}{
   * \int_{\Gamma} \varphi_i u_h
   * = \sum_{k \in {\cal K}} \int_{\Gamma_k} \varphi_i f_k,
   * \qquad \forall \varphi_i \in V_h
   * @f}
   * where $\Gamma = \bigcup_{k \in {\cal K}} \Gamma_k$, $\Gamma_k \subset
   * \partial\Omega$, $\cal K$ is the set of indices and $f_k$ the
   * corresponding boundary functions represented in the function map argument
   * @p boundary_values to this function, and the integrals are evaluated by
   * quadrature. This problem has a non-unique solution in the interior, but
   * it is well defined for the degrees of freedom on the part of the
   * boundary, $\Gamma$, for which we do the integration. The values of
   * $u_h|_\Gamma$, i.e., the nodal values of the degrees of freedom of this
   * function along the boundary, are then what is computed by this function.
   *
   * In case this function is used with $H_{div}$ conforming finite element
   * space, the solution of a different problem is computed, namely: Find
   * $\vec{u}_h \in V_h \subset H(\text{div}; \Omega)$ so that
   * @f{align*}{
   * \int_{\Gamma} (\vec{\varphi}_i \cdot \vec{n}) (\vec{u}_h \cdot \vec{n})
   * = \sum_{k \in {\cal K}} \int_{\Gamma_k} (\vec{\varphi}_i \cdot \vec{n})
   * (\vec{f}_k \cdot \vec{n}),
   * \qquad \forall \vec{\varphi_i} \in V_h,
   * @f}
   * where $\vec{n}$ is an outward normal vector.
   *
   * This function throws an exception if used with $H_\text{curl}$ conforming
   * elements, so the project_boundary_values_curl_conforming_l2() should be
   * used instead.
   *
   * @param[in] mapping The mapping that will be used in the transformations
   * necessary to integrate along the boundary.
   * @param[in] dof The DoFHandler that describes the finite element space and
   * the numbering of degrees of freedom.
   * @param[in] boundary_functions A map from boundary indicators to pointers
   * to functions that describe the desired values on those parts of the
   * boundary marked with this boundary indicator (see
   * @ref GlossBoundaryIndicator "Boundary indicator").
   * The projection happens on only those parts of the boundary whose
   * indicators are represented in this map.
   * @param[in] q The face quadrature used in the integration necessary to
   * compute the @ref GlossMassMatrix "mass matrix" and right hand side of the projection.
   * @param[out] boundary_values The result of this function. It is a map
   * containing all indices of degrees of freedom at the boundary (as covered
   * by the boundary parts in @p boundary_functions) and the computed dof
   * value for this degree of freedom. For each degree of freedom at the
   * boundary, if its index already exists in @p boundary_values then its
   * boundary value will be overwritten, otherwise a new entry with proper
   * index and boundary value for this degree of freedom will be inserted into
   * @p boundary_values.
   * @param[in] component_mapping It is sometimes convenient to project a
   * vector-valued function onto only parts of a finite element space (for
   * example, to project a function with <code>dim</code> components onto the
   * velocity components of a <code>dim+1</code> component DoFHandler for a
   * Stokes problem). To allow for this, this argument allows components to be
   * remapped. If the vector is not empty, it has to have one entry for each
   * vector component of the finite element used in @p dof. This entry is the
   * component number in @p boundary_functions that should be used for this
   * component in @p dof. By default, no remapping is applied.
   *
   * @note Using the *projection* rather than the *interpolation* of
   *   boundary values makes relatively little difference in
   *   practice. That said, it is far more computationally expensive
   *   to compute projections because the require the solution of a
   *   problem that couples all unknowns on the boundary, whereas
   *   interpolation works on one face at a time. On the other hand,
   *   interpolation is only possible for "nodal" finite element
   *   spaces (such as FE_Q, but not FE_Q_Hierarchical), whereas the
   *   projection is always possible. (For some more theoretical
   *   considerations, see the documentation of the first
   *   interpolate_boundary_values() function above.)
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &boundary_functions,
    const Quadrature<dim - 1>                 &q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt>.
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &boundary_function,
    const Quadrature<dim - 1>                 &q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Same as above, but with hp-capabilities.
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const DoFHandler<dim, spacedim>            &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &boundary_functions,
    const hp::QCollection<dim - 1>            &q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt>.
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                                              &boundary_function,
    const hp::QCollection<dim - 1>            &q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Project a function to the boundary of the domain, using the given
   * quadrature formula for the faces. This function identifies the degrees of
   * freedom subject to Dirichlet boundary conditions, adds them to the list
   * of constrained DoFs in @p constraints and sets the respective
   * inhomogeneity to the value resulting from the projection operation. If
   * this routine encounters a DoF that already is constrained (for instance
   * by a hanging node constraint, see below, or any other type of constraint,
   * e.g. from periodic boundary conditions), the old setting of the
   * constraint (dofs the entry is constrained to, inhomogeneities) is kept
   * and nothing happens.
   *
   * @note When combining adaptively refined meshes with hanging node
   * constraints and boundary conditions like from the current function within
   * one AffineConstraints object, the hanging node constraints should always
   * be set first, and then the boundary conditions since boundary conditions
   * are not set in the second operation on degrees of freedom that are
   * already constrained. This makes sure that the discretization remains
   * conforming as is needed. See the discussion on conflicting constraints in
   * the topic on
   * @ref constraints.
   *
   * If @p component_mapping is empty, it is assumed that the number of
   * components of @p boundary_function matches that of the finite element
   * used by @p dof.
   *
   * In 1d, projection equals interpolation. Therefore,
   * interpolate_boundary_values is called.
   *
   * @arg @p component_mapping: if the components in @p boundary_functions and
   * @p dof do not coincide, this vector allows them to be remapped. If the
   * vector is not empty, it has to have one entry for each component in @p
   * dof. This entry is the component number in @p boundary_functions that
   * should be used for this component in @p dof. By default, no remapping is
   * applied.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim>    &mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &boundary_functions,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQ@<dim,spacedim@>(1)</tt>.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
                              &boundary_function,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});

  /**
   * This function is an updated version of the
   * project_boundary_values_curl_conforming function. The intention is to fix
   * a problem when using the previous function in conjunction with
   * non-rectangular geometries (i.e. elements with non-rectangular faces). The
   * L2-projection method used has been taken from the paper "Electromagnetic
   * scattering simulation using an H (curl) conforming hp-finite element
   * method in three dimensions" by PD Ledger, K Morgan and O Hassan ( Int. J.
   * Num. Meth. Fluids, Volume 53, Issue 8, pages 1267-1296).
   *
   * This function will compute constraints that correspond to Dirichlet
   * boundary conditions of the form
   * $\vec{n}\times\vec{E}=\vec{n}\times\vec{F}$ i.e. the tangential
   * components of $\vec{E}$ and $f$ shall coincide.
   *
   * <h4>Computing constraints</h4>
   *
   * To compute the constraints we use a projection method based upon the
   * paper mentioned above. In 2d this is done in a single stage for the
   * edge-based shape functions, regardless of the order of the finite element.
   * In 3d this is done in two stages, edges first and then faces.
   *
   * For each cell, each edge, $e$, is projected by solving the linear system
   * $Ax=b$ where $x$ is the vector of constraints on degrees of freedom on the
   * edge and
   *
   * $A_{ij} = \int_{e} (\vec{s}_{i}\cdot\vec{t})(\vec{s}_{j}\cdot\vec{t}) dS$
   *
   * $b_{i} = \int_{e} (\vec{s}_{i}\cdot\vec{t})(\vec{F}\cdot\vec{t}) dS$
   *
   * with $\vec{s}_{i}$ the $i^{th}$ shape function and $\vec{t}$ the tangent
   * vector.
   *
   * Once all edge constraints, $x$, have been computed, we may compute the
   * face constraints in a similar fashion, taking into account the residuals
   * from the edges.
   *
   * For each face on the cell, $f$, we solve the linear system $By=c$ where
   * $y$ is the vector of constraints on degrees of freedom on the face and
   *
   * $B_{ij} = \int_{f} (\vec{n} \times \vec{s}_{i}) \cdot (\vec{n} \times
   * \vec{s}_{j}) dS$
   *
   * $c_{i} = \int_{f} (\vec{n} \times \vec{r}) \cdot (\vec{n} \times
   * \vec{s}_i) dS$
   *
   * and $\vec{r} = \vec{F} - \sum_{e \in f} \sum{i \in e} x_{i}\vec{s}_i$,
   * the edge residual.
   *
   * The resulting constraints are then given in the solutions $x$ and $y$.
   *
   * If the AffineConstraints @p constraints contained values or other
   * constraints before, the new ones are added or the old ones overwritten,
   * if a node of the boundary part to be used was already in the list of
   * constraints. This is handled by using inhomogeneous constraints. Please
   * note that when combining adaptive meshes and this kind of constraints,
   * the Dirichlet conditions should be set first, and then completed by
   * hanging node constraints, in order to make sure that the discretization
   * remains consistent. See the discussion on conflicting constraints in the
   * topic on
   * @ref constraints.
   *
   * <h4>Arguments to this function</h4>
   *
   * This function is explicitly for use with FE_Nedelec elements, or with
   * FESystem elements which contain FE_Nedelec elements. It will throw an
   * exception if called with any other finite element. The user must ensure
   * that FESystem elements are correctly setup when using this function as
   * this check not possible in this case.
   *
   * The second argument of this function denotes the first vector component
   * of the finite element which corresponds to the vector function that you
   * wish to constrain. For example, if we are solving Maxwell's equations in
   * 3d and have components $(E_x,E_y,E_z,B_x,B_y,B_z)$ and we want the
   * boundary conditions $\vec{n}\times\vec{B}=\vec{n}\times\vec{f}$, then @p
   * first_vector_component would be 3. The @p boundary_function must return 6
   * components in this example, with the first 3 corresponding to $\vec{E}$
   * and the second 3 corresponding to $\vec{B}$. Vectors are implicitly
   * assumed to have exactly <code>dim</code> components that are ordered in
   * the same way as we usually order the coordinate directions, i.e. $x$-,
   * $y$-, and finally $z$-component.
   *
   * The parameter @p boundary_component corresponds to the number @p
   * boundary_id of the face. numbers::internal_face_boundary_id is an illegal
   * value, since it is reserved for interior faces.
   *
   * The last argument is denoted to compute the normal vector $\vec{n}$ at
   * the boundary points.
   *
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const DoFHandler<dim, dim>  &dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, number> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<number>   &constraints,
    const Mapping<dim>          &mapping);


  /**
   * hp-namespace version of project_boundary_values_curl_conforming_l2
   * (above).
   *
   * @ingroup constraints
   */
  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const DoFHandler<dim, dim>            &dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, number>           &boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<number>             &constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);


  /**
   * Compute constraints that correspond to boundary conditions of the form
   * $\vec{n}^T\vec{u}=\vec{n}^T\vec{f}$, i.e. the normal components of the
   * solution $u$ and a given $f$ shall coincide. The function $f$ is given by
   * @p boundary_function and the resulting constraints are added to @p
   * constraints for faces with boundary indicator @p boundary_component.
   *
   * This function is explicitly written to use with the FE_RaviartThomas
   * elements. Thus it throws an exception, if it is called with other finite
   * elements.
   *
   * If the AffineConstraints object @p constraints contained values or other
   * constraints before, the new ones are added or the old ones overwritten,
   * if a node of the boundary part to be used was already in the list of
   * constraints. This is handled by using inhomogeneous constraints. Please
   * note that when combining adaptive meshes and this kind of constraints,
   * the Dirichlet conditions should be set first, and then completed by
   * hanging node constraints, in order to make sure that the discretization
   * remains consistent. See the discussion on conflicting constraints in the
   * topic on
   * @ref constraints.
   *
   * The argument @p first_vector_component denotes the first vector component
   * in the finite element that corresponds to the vector function $\vec{u}$
   * that you want to constrain. Vectors are implicitly assumed to have
   * exactly <code>dim</code> components that are ordered in the same way as
   * we usually order the coordinate directions, i.e., $x$-, $y$-, and finally
   * $z$-component.
   *
   * The parameter @p boundary_component corresponds to the @p boundary_id of
   * the faces where the boundary conditions are applied.
   * numbers::internal_face_boundary_id is an illegal value, since it is
   * reserved for interior faces. The @p mapping is used to compute the normal
   * vector $\vec{n}$ at the boundary points.
   *
   * <h4>Computing constraints</h4>
   *
   * To compute the constraints we use interpolation operator proposed in
   * Brezzi, Fortin (Mixed and Hybrid Finite Element Methods, Springer, 1991)
   * on every face located at the boundary.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, typename number, typename number2 = number>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim, dim>   &dof_handler,
    const unsigned int            first_vector_component,
    const Function<dim, number2> &boundary_function,
    const types::boundary_id      boundary_component,
    AffineConstraints<number>    &constraints,
    const Mapping<dim>           &mapping);

  /**
   * Same as above for the hp-namespace.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, typename number, typename number2 = number>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim, dim>            &dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, number2>          &boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<number>             &constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);

  /** @} */
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h
