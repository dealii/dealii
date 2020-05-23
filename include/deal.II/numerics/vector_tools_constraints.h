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


#ifndef dealii_vector_tools_constraints_h
#define dealii_vector_tools_constraints_h

#include <deal.II/base/config.h>

#include <map>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <typename number>
class AffineConstraints;
template <int dim, int spacedim>
struct StaticMappingQ1;
template <int dim, typename Number>
class Function;
template <int dim, int spacedim>
class Mapping;


namespace VectorTools
{
  /**
   * @name Interpolation and projection
   */
  //@{

  /**
   * This function computes the constraints that correspond to boundary
   * conditions of the form $\vec u \cdot \vec n=\vec u_\Gamma \cdot \vec n$,
   * i.e., normal flux constraints where $\vec u$ is a vector-valued solution
   * variable and $\vec u_\Gamma$ is a prescribed vector field whose normal
   * component we want to be equal to the normal component of the solution.
   * These conditions have exactly the form handled by the
   * AffineConstraints class, in that they relate a <i>linear
   * combination</i> of boundary degrees of freedom to a corresponding
   * value (the inhomogeneity of the constraint). Consequently, the current
   * function creates a list of constraints that are written into an
   * AffineConstraints container. This object may already have some
   * content, for example from hanging node constraints, that remains
   * untouched. These constraints have to be applied to the linear system
   * like any other such constraints, i.e., you have to condense the linear
   * system with the constraints before solving, and you have to distribute
   * the solution vector afterwards.
   *
   * This function treats a more general case than
   * VectorTools::compute_no_normal_flux_constraints() (which can only handle
   * the case where $\vec u_\Gamma \cdot \vec n = 0$, and is used in
   * step-31 and step-32). However, because everything that would apply
   * to that function also applies as a special case to the current
   * function, the following discussion is relevant to both.
   *
   * @note This function doesn't make much sense in 1d, so it throws an
   *   exception if @p dim equals one.
   *
   *
   * <h4>Arguments to this function</h4>
   *
   * The second argument of this function denotes the first vector component
   * in the finite element that corresponds to the vector function that you
   * want to constrain. For example, if we were solving a Stokes equation in
   * 2d and the finite element had components $(u,v,p)$, then @p
   * first_vector_component needs to be zero if you intend to constraint
   * the vector $(u,v)^T \cdot \vec n = \vec u_\Gamma \cdot \vec n$.
   * On the other hand, if we solved the
   * Maxwell equations in 3d and the finite element has components
   * $(E_x,E_y,E_z,B_x,B_y,B_z)$ and we want the boundary condition $\vec
   * B\cdot \vec n=\vec B_\Gamma\cdot \vec n$, then @p first_vector_component
   * would be 3. Vectors are implicitly assumed to have exactly
   * <code>dim</code> components that are ordered in the same way as we
   * usually order the coordinate directions, i.e. $x$-, $y$-, and finally
   * $z$-component. The function assumes, but can't check, that the vector
   * components in the range
   * <code>[first_vector_component,first_vector_component+dim)</code> come
   * from the same base finite element. For example, in the Stokes example
   * above, it would not make sense to use a
   * <code>FESystem@<dim@>(FE_Q@<dim@>(2), 1, FE_Q@<dim@>(1), dim)</code>
   * (note that the first velocity vector component is a $Q_2$ element,
   * whereas all the other ones are $Q_1$ elements) as there would be points
   * on the boundary where the $x$-velocity is defined but no corresponding
   * $y$- or $z$-velocities.
   *
   * The third argument denotes the set of boundary indicators on which the
   * boundary condition is to be enforced. Note that, as explained below, this
   * is one of the few functions where it makes a difference where we call the
   * function multiple times with only one boundary indicator, or whether we
   * call the function once with the whole set of boundary indicators at once.
   *
   * Argument four (@p function_map) describes the boundary function $\vec
   * u_\Gamma$ for each boundary id. The function <code>function_map[id]</code>
   * is used on boundary with id @p id taken from the set @p boundary_ids.
   * Each function in @p function_map is expected to have @p dim
   * components, which are used independent of @p first_vector_component.
   *
   * The mapping argument is used to compute the boundary points at which the
   * function needs to request the normal vector $\vec n$ from the boundary
   * description.
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
   *
   * <h4>Computing constraints in 2d</h4>
   *
   * Computing these constraints requires some smarts. The main question
   * revolves around the question what the normal vector is. Consider the
   * following situation:
   *
   * <p ALIGN="center">
   * @image html no_normal_flux_1.png
   * </p>
   *
   * Here, we have two cells that use a bilinear mapping (i.e.,
   * MappingQGeneric(1)). Consequently, for each of the cells, the normal
   * vector is perpendicular to the straight edge. If the two edges at the top
   * and right are meant to approximate a curved boundary (as indicated by the
   * dashed line), then neither of the two computed normal vectors are equal
   * to the exact normal vector (though they approximate it as the mesh is
   * refined further). What is worse, if we constrain $\vec u \cdot \vec n=
   * \vec u_\Gamma \cdot \vec n$ at the common vertex with the normal vector
   * from both cells, then we constrain the vector $\vec u$ with respect to
   * two linearly independent vectors; consequently, the constraint would be
   * $\vec u=\vec u_\Gamma$ at this point (i.e. <i>all</i> components of the
   * vector), which is not what we wanted.
   *
   * To deal with this situation, the algorithm works in the following way: at
   * each point where we want to constrain $\vec u$, we first collect all
   * normal vectors that adjacent cells might compute at this point. We then
   * do not constrain $\vec u \cdot \vec n=\vec u_\Gamma \cdot \vec n$ for
   * <i>each</i> of these normal vectors but only for the <i>average</i> of
   * the normal vectors. In the example above, we therefore record only a
   * single constraint $\vec u \cdot \vec {\bar n}=\vec u_\Gamma \cdot \vec
   * {\bar n}$, where $\vec {\bar n}$ is the average of the two indicated
   * normal vectors.
   *
   * Unfortunately, this is not quite enough. Consider the situation here:
   *
   * <p ALIGN="center">
   * @image html no_normal_flux_2.png
   * </p>
   *
   * If again the top and right edges approximate a curved boundary, and the
   * left boundary a separate boundary (for example straight) so that the
   * exact boundary has indeed a corner at the top left vertex, then the above
   * construction would not work: here, we indeed want the constraint that
   * $\vec u$ at this point (because the normal velocities with respect to
   * both the left normal as well as the top normal vector should be zero),
   * not that the velocity in the direction of the average normal vector is
   * zero.
   *
   * Consequently, we use the following heuristic to determine whether all
   * normal vectors computed at one point are to be averaged: if two normal
   * vectors for the same point are computed on <i>different</i> cells, then
   * they are to be averaged. This covers the first example above. If they are
   * computed from the same cell, then the fact that they are different is
   * considered indication that they come from different parts of the boundary
   * that might be joined by a real corner, and must not be averaged.
   *
   * There is one problem with this scheme. If, for example, the same domain
   * we have considered above, is discretized with the following mesh, then we
   * get into trouble:
   *
   * <p ALIGN="center">
   * @image html no_normal_flux_3.png
   * </p>
   *
   * Here, the algorithm assumes that the boundary does not have a corner at
   * the point where faces $F1$ and $F2$ join because at that point there are
   * two different normal vectors computed from different cells. If you intend
   * for there to be a corner of the exact boundary at this point, the only
   * way to deal with this is to assign the two parts of the boundary
   * different boundary indicators and call this function twice, once for each
   * boundary indicators; doing so will yield only one normal vector at this
   * point per invocation (because we consider only one boundary part at a
   * time), with the result that the normal vectors will not be averaged. This
   * situation also needs to be taken into account when using this function
   * around reentrant corners on Cartesian meshes. If normal-flux boundary
   * conditions are to be enforced on non-Cartesian meshes around reentrant
   * corners, one may even get cycles in the constraints as one will in
   * general constrain different components from the two sides. In that case,
   * set a no-slip constraint on the reentrant vertex first.
   *
   *
   * <h4>Computing constraints in 3d</h4>
   *
   * The situation is more complicated in 3d. Consider the following case
   * where we want to compute the constraints at the marked vertex:
   *
   * <p ALIGN="center">
   * @image html no_normal_flux_4.png
   * </p>
   *
   * Here, we get four different normal vectors, one from each of the four
   * faces that meet at the vertex. Even though they may form a complete set
   * of vectors, it is not our intent to constrain all components of the
   * vector field at this point. Rather, we would like to still allow
   * tangential flow, where the term "tangential" has to be suitably defined.
   *
   * In a case like this, the algorithm proceeds as follows: for each cell
   * that has computed two tangential vectors at this point, we compute the
   * unconstrained direction as the outer product of the two tangential
   * vectors (if necessary multiplied by minus one). We then average these
   * tangential vectors. Finally, we compute constraints for the two
   * directions perpendicular to this averaged tangential direction.
   *
   * There are cases where one cell contributes two tangential directions and
   * another one only one; for example, this would happen if both top and
   * front faces of the left cell belong to the boundary selected whereas only
   * the top face of the right cell belongs to it, maybe indicating that the
   * entire front part of the domain is a smooth manifold whereas the top
   * really forms two separate manifolds that meet in a ridge, and that
   * normal-flux boundary conditions are only desired on the front manifold
   * and the right one on top. In cases like these, it's difficult to define
   * what should happen. The current implementation simply ignores the one
   * contribution from the cell that only contributes one normal vector. In
   * the example shown, this is acceptable because the normal vector for the
   * front face of the left cell is the same as the normal vector provided by
   * the front face of the right cell (the surface is planar) but it would be
   * a problem if the front manifold would be curved. Regardless, it is
   * unclear how one would proceed in this case and ignoring the single cell
   * is likely the best one can do.
   *
   *
   * <h4>Results</h4>
   *
   * Because it makes for good pictures, here are two images of vector fields
   * on a circle and on a sphere to which the constraints computed by this
   * function have been applied (for illustration purposes, we enforce zero
   * normal flux, which can more easily be computed using
   * VectorTools::compute_no_normal_flux_constraints(), as this must
   * lead to a <i>tangential</i> vector field):
   *
   * <p ALIGN="center">
   * @image html no_normal_flux_5.png
   * @image html no_normal_flux_6.png
   * </p>
   *
   * The vectors fields are not physically reasonable but the tangentiality
   * constraint is clearly enforced. The fact that the vector fields are zero
   * at some points on the boundary is an artifact of the way it is created,
   * it is not constrained to be zero at these points.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  compute_nonzero_normal_flux_constraints(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_vector_component,
    const std::set<types::boundary_id> & boundary_ids,
    const std::map<types::boundary_id, const Function<spacedim, double> *>
      &                           function_map,
    AffineConstraints<double> &   constraints,
    const Mapping<dim, spacedim> &mapping =
      StaticMappingQ1<dim, spacedim>::mapping);

  /**
   * This function does the same as the
   * compute_nonzero_normal_flux_constraints() function (see there for more
   * information), but for the simpler case of homogeneous normal-flux
   * constraints, i.e., for imposing the condition
   * $\vec u \cdot \vec n= 0$. This function is used in step-31 and step-32.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  compute_no_normal_flux_constraints(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_vector_component,
    const std::set<types::boundary_id> & boundary_ids,
    AffineConstraints<double> &          constraints,
    const Mapping<dim, spacedim> &       mapping =
      StaticMappingQ1<dim, spacedim>::mapping);

  /**
   * Compute the constraints that correspond to boundary conditions of the
   * form $\vec u \times \vec n=\vec u_\Gamma \times \vec n$, i.e., tangential
   * flow constraints where $\vec u$ is a vector-valued solution
   * variable and $\vec u_\Gamma$ is prescribed vector field whose tangential
   * component(s) we want to be equal to the tangential component(s) of the
   * solution. This function constrains exactly those dim-1 vector-valued
   * components that are left unconstrained by
   * VectorTools::compute_no_normal_flux_constraints(), and leaves the one
   * component unconstrained that is constrained by that function.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  compute_nonzero_tangential_flux_constraints(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_vector_component,
    const std::set<types::boundary_id> & boundary_ids,
    const std::map<types::boundary_id, const Function<spacedim, double> *>
      &                           function_map,
    AffineConstraints<double> &   constraints,
    const Mapping<dim, spacedim> &mapping =
      StaticMappingQ1<dim, spacedim>::mapping);

  /**
   * Same as above for homogeneous tangential-flux constraints.
   *
   * @ingroup constraints
   *
   * @see
   * @ref GlossBoundaryIndicator "Glossary entry on boundary indicators"
   */
  template <int dim, int spacedim, template <int, int> class DoFHandlerType>
  void
  compute_normal_flux_constraints(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_vector_component,
    const std::set<types::boundary_id> & boundary_ids,
    AffineConstraints<double> &          constraints,
    const Mapping<dim, spacedim> &       mapping =
      StaticMappingQ1<dim, spacedim>::mapping);

  //@}
} // namespace VectorTools

DEAL_II_NAMESPACE_CLOSE

#endif // dealii_vector_tools_constraints_h
