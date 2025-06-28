// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_constrained_dofs_h
#define dealii_mg_constrained_dofs_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/affine_constraints.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif


/**
 * Collection of boundary constraints and refinement edge constraints for
 * level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public EnableObserverPointer
{
public:
  using size_dof = std::vector<std::set<types::global_dof_index>>::size_type;
  /**
   * Fill the internal data structures with hanging node constraints extracted
   * from the dof handler object. Works with natural boundary conditions only.
   * There exists a sister function setting up boundary constraints as well.
   *
   * This function ensures that on every level, degrees of freedom at interior
   * edges of a refinement level are treated corrected but leaves degrees of
   * freedom at the boundary of the domain untouched assuming that no
   * Dirichlet boundary conditions for them exist.
   *
   * Furthermore, this call sets up an AffineConstraints object on each
   * level that contains possible periodicity constraints in case those
   * have been added to the underlying triangulation. The AffineConstraints
   * object can be queried by get_level_constraints(level). Note that the
   * current implementation of periodicity constraints in this class does
   * not support rotation matrices in the periodicity definition, i.e., the
   * respective argument in the GridTools::collect_periodic_faces() may not
   * be different from the identity matrix.
   * If no level_relevant_dofs are passed as the second argument, the function
   * uses the locally relevant level DoFs, extracted by
   * DoFTools::extract_locally_relevant_level_dofs(). Otherwise, the
   * user-provided IndexSets, which should define a superset of locally relevant
   * DoFs, are used on each level to allow the user to add additional indices to
   * the set of constrained DoFs.
   */
  template <int dim, int spacedim>
  void
  initialize(const DoFHandler<dim, spacedim> &dof,
             const MGLevelObject<IndexSet>   &level_relevant_dofs =
               MGLevelObject<IndexSet>());

  /**
   * Fill the internal data structures with information
   * about Dirichlet boundary dofs.
   *
   * The initialize() function has to be called before
   * to set hanging node constraints.
   *
   * This function can be called multiple times to allow considering
   * different sets of boundary_ids for different components.
   */
  template <int dim, int spacedim>
  void
  make_zero_boundary_constraints(
    const DoFHandler<dim, spacedim>    &dof,
    const std::set<types::boundary_id> &boundary_ids,
    const ComponentMask                &component_mask = {});

  /**
   * Add Dirichlet boundary dofs to the internal data structures
   * on level @p level.
   * The indices are restricted to the set of locally relevant
   * level dofs.
   */
  template <int dim, int spacedim>
  void
  add_boundary_indices(const DoFHandler<dim, spacedim> &dof,
                       const unsigned int               level,
                       const IndexSet                  &boundary_indices);

  /**
   * Add user defined constraints to be used on level @p level.
   *
   * The user can call this function multiple times and any new,
   * conflicting constraints will overwrite the previous constraints
   * for that DoF.
   *
   * Before the transfer, the user defined constraints will be distributed
   * to the source vector, and then any DoF index set using
   * make_zero_boundary_constraints() will be overwritten with
   * value zero.
   *
   * @note This is currently only implemented for MGTransferMatrixFree.
   */
  void
  add_user_constraints(const unsigned int               level,
                       const AffineConstraints<double> &constraints_on_level);

  /**
   * Fill the internal data structures with information
   * about no normal flux boundary dofs.
   *
   * This function is limited to meshes whose no normal flux boundaries
   * have faces which are normal to the x-, y-, or z-axis. Also, for a
   * specific boundary id, all faces must be facing in the same direction,
   * i.e., a boundary normal to the x-axis must have a different boundary
   * id than a boundary normal to the y- or z-axis and so on. If the mesh
   * was produced, for example, using the <tt>GridGenerator::hyper_cube()</tt>
   * function, setting <tt>colorize=true</tt> during mesh generation and calling
   * make_no_normal_flux_constraints() for each no normal flux boundary is
   * sufficient.
   */
  template <int dim, int spacedim>
  void
  make_no_normal_flux_constraints(const DoFHandler<dim, spacedim> &dof,
                                  const types::boundary_id         bid,
                                  const unsigned int first_vector_component);

  /**
   * Clear the user constraints on all levels.
   */
  void
  clear_user_constraints();

  /**
   * Reset the data structures.
   */
  void
  clear();

  /**
   * Determine whether a dof index is subject to a boundary constraint.
   */
  bool
  is_boundary_index(const unsigned int            level,
                    const types::global_dof_index index) const;

  /**
   * Determine whether a dof index is at the refinement edge.
   */
  bool
  at_refinement_edge(const unsigned int            level,
                     const types::global_dof_index index) const;


  /**
   * Determine whether the (i,j) entry of the interface matrix
   * on a given level should be set. This is taken in terms of
   * dof i, that is, return true if i is at a refinement edge,
   * j is not, and both are not on the external boundary.
   */
  bool
  is_interface_matrix_entry(const unsigned int            level,
                            const types::global_dof_index i,
                            const types::global_dof_index j) const;

  /**
   * Return the indices of level dofs on the given level that are subject to
   * Dirichlet boundary conditions (as set by the @p function_map parameter in
   * initialize()).  The indices are restricted to the set of locally relevant
   * level dofs.
   */
  const IndexSet &
  get_boundary_indices(const unsigned int level) const;


  /**
   * Return the indices of dofs on the given level that lie on an refinement
   * edge (dofs on faces to neighbors that are coarser).
   */
  const IndexSet &
  get_refinement_edge_indices(unsigned int level) const;


  /**
   * Return if Dirichlet boundary indices are set in initialize().
   */
  bool
  have_boundary_indices() const;

  /**
   * Return the AffineConstraints object for a given level, containing
   * periodicity constraints (if enabled on the triangulation).
   */
  const AffineConstraints<double> &
  get_level_constraints(const unsigned int level) const;

  /**
   * Return the user defined constraint matrix for a given level. These
   * constraints are set using the function add_user_constraints() and
   * should not contain constraints for DoF indices set in
   * make_zero_boundary_constraints() as they will be overwritten during
   * the transfer.
   */
  const AffineConstraints<double> &
  get_user_constraint_matrix(const unsigned int level) const;

  /**
   * Merge selected constraints of a specified @p level into a given single
   * AffineConstraints object.
   *
   * @param constraints AffineConstraints object to be filled.
   * @param level Refinement to be considered.
   * @param add_boundary_indices Add boundary indices.
   * @param add_refinement_edge_indices Add refinement-edge indices.
   * @param add_level_constraints Add level constraints including the one passed
   *   during initialize() and periodicity constraints.
   * @param add_user_constraints Add user constraints.
   */
  template <typename Number>
  void
  merge_constraints(AffineConstraints<Number> &constraints,
                    const unsigned int         level,
                    const bool                 add_boundary_indices,
                    const bool                 add_refinement_edge_indices,
                    const bool                 add_level_constraints,
                    const bool                 add_user_constraints) const;

private:
  /**
   * The indices of boundary dofs for each level.
   */
  std::vector<IndexSet> boundary_indices;

  /**
   * The degrees of freedom on a given level that live on the refinement edge
   * between the level and cells on a coarser level.
   */
  std::vector<IndexSet> refinement_edge_indices;

  /**
   * Constraint matrices containing information regarding potential
   * periodic boundary conditions for each level .
   */
  std::vector<AffineConstraints<double>> level_constraints;

  /**
   * Constraint matrices defined by user.
   */
  std::vector<AffineConstraints<double>> user_constraints;
};



inline bool
MGConstrainedDoFs::is_boundary_index(const unsigned int            level,
                                     const types::global_dof_index index) const
{
  if (boundary_indices.empty())
    return false;

  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level].is_element(index);
}



inline bool
MGConstrainedDoFs::at_refinement_edge(const unsigned int            level,
                                      const types::global_dof_index index) const
{
  AssertIndexRange(level, refinement_edge_indices.size());

  return refinement_edge_indices[level].is_element(index);
}



inline bool
MGConstrainedDoFs::is_interface_matrix_entry(
  const unsigned int            level,
  const types::global_dof_index i,
  const types::global_dof_index j) const
{
  const IndexSet &interface_dofs_on_level =
    this->get_refinement_edge_indices(level);

  return interface_dofs_on_level.is_element(i)     // at_refinement_edge(i)
         && !interface_dofs_on_level.is_element(j) // !at_refinement_edge(j)
         && !this->is_boundary_index(level, i)     // !on_boundary(i)
         && !this->is_boundary_index(level, j);    // !on_boundary(j)
}



inline const IndexSet &
MGConstrainedDoFs::get_boundary_indices(const unsigned int level) const
{
  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level];
}



inline const IndexSet &
MGConstrainedDoFs::get_refinement_edge_indices(unsigned int level) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  return refinement_edge_indices[level];
}



inline bool
MGConstrainedDoFs::have_boundary_indices() const
{
  return boundary_indices.size() != 0;
}



inline const AffineConstraints<double> &
MGConstrainedDoFs::get_level_constraints(const unsigned int level) const
{
  AssertIndexRange(level, level_constraints.size());
  return level_constraints[level];
}



inline const AffineConstraints<double> &
MGConstrainedDoFs::get_user_constraint_matrix(const unsigned int level) const
{
  AssertIndexRange(level, user_constraints.size());
  return user_constraints[level];
}



DEAL_II_NAMESPACE_CLOSE

#endif
