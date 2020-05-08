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


#ifndef dealii_vector_tools_boundary_h
#define dealii_vector_tools_boundary_h

#include <deal.II/base/config.h>

#include <deal.II/hp/mapping_collection.h>

#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
class Function;
namespace hp
{
  template <int dim, int spacedim>
  class DoFHandler;
  template <int dim>
  class QCollection;
} // namespace hp

namespace VectorTools
{
  /**
   * @name Interpolation and projection
   */
  //@{

  /**
   * Compute Dirichlet boundary conditions.  This function makes up a map of
   * degrees of freedom subject to Dirichlet boundary conditions and the
   * corresponding values to be assigned to them, by interpolation around the
   * boundary. For each degree of freedom at the boundary, if its index
   * already exists in @p boundary_values then its boundary value will be
   * overwritten, otherwise a new entry with proper index and boundary value
   * for this degree of freedom will be inserted into @p boundary_values.
   *
   * The parameter @p function_map provides a list of boundary indicators to
   * be handled by this function and corresponding boundary value functions.
   * The keys of this map correspond to the number @p boundary_id of the face.
   * numbers::internal_face_boundary_id is an illegal value for this key since
   * it is reserved for interior faces. For an example of how to use this
   * argument with a non-empty map, see the step-16 tutorial program.
   *
   * The flags in the last parameter, @p component_mask denote which
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
   * See the general documentation of this namespace for more information.
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * Like the previous function, but take a mapping collection to go with the
   * hp::DoFHandler object.
   */
  template <int dim, int spacedim, typename number>
  void
  interpolate_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * Same function as above, but taking only one pair of boundary indicator
   * and corresponding boundary function. The same comments apply as for the
   * previous function, in particular about the use of the component mask and
   * the requires size of the function object.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &             mapping,
    const DoFHandlerType<dim, spacedim> &      dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * Call the other interpolate_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>. The same comments
   * apply as for the previous function, in particular about the use of the
   * component mask and the requires size of the function object.
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &      dof,
    const types::boundary_id                   boundary_component,
    const Function<spacedim, number> &         boundary_function,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());


  /**
   * Call the other interpolate_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>. The same comments
   * apply as for the previous function, in particular about the use of the
   * component mask and the requires size of the function object.
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        function_map,
    std::map<types::global_dof_index, number> &boundary_values,
    const ComponentMask &component_mask = ComponentMask());


  /**
   * Insert the (algebraic) constraints due to Dirichlet boundary conditions
   * into a AffineConstraints @p constraints. This function identifies the
   * degrees of freedom subject to Dirichlet boundary conditions, adds them to
   * the list of constrained DoFs in @p constraints and sets the respective
   * inhomogeneity to the value interpolated around the boundary. If this
   * routine encounters a DoF that already is constrained (for instance by a
   * hanging node constraint, see below, or any other type of constraint, e.g.
   * from periodic boundary conditions), the old setting of the constraint
   * (dofs the entry is constrained to, inhomogeneities) is kept and nothing
   * happens.
   *
   * @note When combining adaptively refined meshes with hanging node
   * constraints and boundary conditions like from the current function within
   * one AffineConstraints object, the hanging node constraints should always
   * be set first, and then the boundary conditions since boundary conditions
   * are not set in the second operation on degrees of freedom that are
   * already constrained. This makes sure that the discretization remains
   * conforming as is needed. See the discussion on conflicting constraints in
   * the module on
   * @ref constraints.
   *
   * The parameter @p boundary_component corresponds to the number @p
   * boundary_id of the face.
   *
   * The flags in the last parameter, @p component_mask denote which
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
   * values for the pressure, then the component mask should correspond to
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
   * See the general documentation of this namespace for more information.
   *
   * @ingroup constraints
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask = ComponentMask());

  /**
   * Same function as above, but taking only one pair of boundary indicator
   * and corresponding boundary function. The same comments apply as for the
   * previous function, in particular about the use of the component mask and
   * the requires size of the function object.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const Mapping<dim, spacedim> &       mapping,
    const DoFHandlerType<dim, spacedim> &dof,
    const types::boundary_id             boundary_component,
    const Function<spacedim, number> &   boundary_function,
    AffineConstraints<number> &          constraints,
    const ComponentMask &                component_mask = ComponentMask());

  /**
   * Call the other interpolate_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>. The same comments
   * apply as for the previous function, in particular about the use of the
   * component mask and the requires size of the function object.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const types::boundary_id             boundary_component,
    const Function<spacedim, number> &   boundary_function,
    AffineConstraints<number> &          constraints,
    const ComponentMask &                component_mask = ComponentMask());


  /**
   * Call the other interpolate_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>. The same comments
   * apply as for the previous function, in particular about the use of the
   * component mask and the requires size of the function object.
   *
   * @ingroup constraints
   */
  template <int dim,
            int spacedim,
            template <int, int> class DoFHandlerType,
            typename number>
  void
  interpolate_boundary_values(
    const DoFHandlerType<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        function_map,
    AffineConstraints<number> &constraints,
    const ComponentMask &      component_mask = ComponentMask());


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
   * This function throws exception if used with $H_{curl}$ conforming elements,
   * so the project_boundary_values_curl_conforming() should be used instead.
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
   * compute the mass matrix and right hand side of the projection.
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
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>.
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_function,
    const Quadrature<dim - 1> &                q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Same as above, but for objects of type hp::DoFHandler
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::MappingCollection<dim, spacedim> &mapping,
    const hp::DoFHandler<dim, spacedim> &       dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_functions,
    const hp::QCollection<dim - 1> &           q,
    std::map<types::global_dof_index, number> &boundary_values,
    std::vector<unsigned int>                  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>.
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const hp::DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                                        boundary_function,
    const hp::QCollection<dim - 1> &           q,
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
   * the module on
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
    const Mapping<dim, spacedim> &   mapping,
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_functions,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});

  /**
   * Call the project_boundary_values() function, see above, with
   * <tt>mapping=MappingQGeneric@<dim,spacedim@>(1)</tt>.
   *
   * @ingroup constraints
   */
  template <int dim, int spacedim, typename number>
  void
  project_boundary_values(
    const DoFHandler<dim, spacedim> &dof,
    const std::map<types::boundary_id, const Function<spacedim, number> *>
      &                        boundary_function,
    const Quadrature<dim - 1> &q,
    AffineConstraints<number> &constraints,
    std::vector<unsigned int>  component_mapping = {});


  /**
   * Compute constraints that correspond to boundary conditions of the form
   * $\vec{n}\times\vec{u}=\vec{n}\times\vec{f}$, i.e. the tangential
   * components of $u$ and $f$ shall coincide.
   *
   * If the AffineConstraints @p constraints contained values or other
   * constraints before, the new ones are added or the old ones overwritten,
   * if a node of the boundary part to be used was already in the list of
   * constraints. This is handled by using inhomogeneous constraints. Please
   * note that when combining adaptive meshes and this kind of constraints,
   * the Dirichlet conditions should be set first, and then completed by
   * hanging node constraints, in order to make sure that the discretization
   * remains consistent. See the discussion on conflicting constraints in the
   * module on
   * @ref constraints.
   *
   * This function is explicitly written to use with the FE_Nedelec elements.
   * Thus it throws an exception, if it is called with other finite elements.
   *
   * The second argument of this function denotes the first vector component
   * in the finite element that corresponds to the vector function that you
   * want to constrain. For example, if we want to solve Maxwell's equations
   * in 3d and the finite element has components $(E_x,E_y,E_z,B_x,B_y,B_z)$
   * and we want the boundary conditions
   * $\vec{n}\times\vec{B}=\vec{n}\times\vec{f}$, then @p
   * first_vector_component would be 3. Vectors are implicitly assumed to have
   * exactly <code>dim</code> components that are ordered in the same way as
   * we usually order the coordinate directions, i.e. $x$-, $y$-, and finally
   * $z$-component.
   *
   * The parameter @p boundary_component corresponds to the number @p
   * boundary_id of the face. numbers::internal_face_boundary_id is an illegal
   * value, since it is reserved for interior faces.
   *
   * The last argument is denoted to compute the normal vector $\vec{n}$ at
   * the boundary points.
   *
   * <h4>Computing constraints</h4>
   *
   * To compute the constraints we use projection-based interpolation as
   * proposed in Solin, Segeth and Dolezel (Higher order finite elements,
   * Chapman&amp;Hall, 2004) on every face located at the boundary.
   *
   * First one projects $\vec{f}$ on the lowest-order edge shape functions.
   * Then the remaining part $(I-P_0)\vec{f}$ of the function is projected on
   * the remaining higher-order edge shape functions. In the last step we
   * project $(I-P_0-P_e)\vec{f}$ on the bubble shape functions defined on the
   * face.
   *
   * @deprecated Use the project_boundary_values_curl_conforming_l2() function
   * instead of this one.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim>
  DEAL_II_DEPRECATED void
  project_boundary_values_curl_conforming(
    const DoFHandler<dim, dim> & dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, double> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<double> &  constraints,
    const Mapping<dim> &         mapping = StaticMappingQ1<dim>::mapping);

  /**
   * Same as above for the hp-namespace.
   *
   * @deprecated Use the project_boundary_values_curl_conforming_l2() function
   * instead of this one.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim>
  DEAL_II_DEPRECATED void
  project_boundary_values_curl_conforming(
    const hp::DoFHandler<dim, dim> &       dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, double> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<double> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);

  /**
   * This function is an updated version of the
   * project_boundary_values_curl_conforming function. The intention is to fix
   * a problem when using the previous function in conjunction with non-
   * rectangular geometries (i.e. elements with non-rectangular faces). The
   * L2-projection method used has been taken from the paper "Electromagnetic
   * scattering simulation using an H (curl) conforming hp finite element
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
   * paper mentioned above. In 2D this is done in a single stage for the edge-
   * based shape functions, regardless of the order of the finite element. In
   * 3D this is done in two stages, edges first and then faces.
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
   * module on
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
   * 3D and have components $(E_x,E_y,E_z,B_x,B_y,B_z)$ and we want the
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
    const DoFHandler<dim, dim> & dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, number> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<number> &  constraints,
    const Mapping<dim> &         mapping = StaticMappingQ1<dim>::mapping);


  /**
   * hp-namespace version of project_boundary_values_curl_conforming_l2
   * (above).
   *
   * @ingroup constraints
   */
  template <int dim, typename number>
  void
  project_boundary_values_curl_conforming_l2(
    const hp::DoFHandler<dim, dim> &       dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, number> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<number> &            constraints,
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
   * module on
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
   * Brezzi, Fortin (Mixed and Hybrid (Finite Element Methods, Springer, 1991)
   * on every face located at the boundary.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim>
  void
  project_boundary_values_div_conforming(
    const DoFHandler<dim, dim> & dof_handler,
    const unsigned int           first_vector_component,
    const Function<dim, double> &boundary_function,
    const types::boundary_id     boundary_component,
    AffineConstraints<double> &  constraints,
    const Mapping<dim> &         mapping = StaticMappingQ1<dim>::mapping);

  /**
   * Same as above for the hp-namespace.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim>
  void
  project_boundary_values_div_conforming(
    const hp::DoFHandler<dim, dim> &       dof_handler,
    const unsigned int                     first_vector_component,
    const Function<dim, double> &          boundary_function,
    const types::boundary_id               boundary_component,
    AffineConstraints<double> &            constraints,
    const hp::MappingCollection<dim, dim> &mapping_collection =
      hp::StaticMappingQ1<dim>::mapping_collection);

  // @}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_boundary_h
