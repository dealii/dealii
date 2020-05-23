// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_mg_constrained_dofs_h
#define dealii_mg_constrained_dofs_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/multigrid/mg_tools.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
#endif


/**
 * Collection of boundary constraints and refinement edge constraints for
 * level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public Subscriptor
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
   */
  template <int dim, int spacedim>
  void
  initialize(const DoFHandler<dim, spacedim> &dof);

  /**
   * Fill the internal data structures with values extracted from the dof
   * handler object and apply the boundary values provided.
   *
   * This function internally calls the initialize() function above and the
   * constrains degrees on the external boundary of the domain by calling
   * MGTools::make_boundary_list() with the given second and third argument.
   *
   * @deprecated Use initialize() followed by make_zero_boundary_constraints() instead
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED void
  initialize(const DoFHandler<dim, spacedim> &dof,
             const std::map<types::boundary_id, const Function<spacedim> *>
               &                  function_map,
             const ComponentMask &component_mask = ComponentMask());

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
    const DoFHandler<dim, spacedim> &   dof,
    const std::set<types::boundary_id> &boundary_ids,
    const ComponentMask &               component_mask = ComponentMask());

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
   * Return the AffineConstraints object for a given level, containing
   * periodicity constraints (if enabled on the triangulation).
   *
   * @deprecated Use get_level_constraints instead, which has a more descriptive name.
   */
  DEAL_II_DEPRECATED
  const AffineConstraints<double> &
  get_level_constraint_matrix(const unsigned int level) const;

  /**
   * Return the user defined constraint matrix for a given level. These
   * constraints are set using the function add_user_constraints() and
   * should not contain constraints for DoF indices set in
   * make_zero_boundary_constraints() as they will be overwritten during
   * the transfer.
   */
  const AffineConstraints<double> &
  get_user_constraint_matrix(const unsigned int level) const;

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


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::initialize(const DoFHandler<dim, spacedim> &dof)
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  level_constraints.clear();
  user_constraints.clear();

  const unsigned int nlevels = dof.get_triangulation().n_global_levels();

  // At this point level_constraint and refinement_edge_indices are empty.
  refinement_edge_indices.resize(nlevels);
  level_constraints.resize(nlevels);
  user_constraints.resize(nlevels);
  for (unsigned int l = 0; l < nlevels; ++l)
    {
      IndexSet relevant_dofs;
      DoFTools::extract_locally_relevant_level_dofs(dof, l, relevant_dofs);
      level_constraints[l].reinit(relevant_dofs);

      // Loop through relevant cells and faces finding those which are periodic
      // neighbors.
      typename DoFHandler<dim, spacedim>::cell_iterator cell = dof.begin(l),
                                                        endc = dof.end(l);
      for (; cell != endc; ++cell)
        if (cell->level_subdomain_id() != numbers::artificial_subdomain_id)
          {
            for (auto f : GeometryInfo<dim>::face_indices())
              if (cell->has_periodic_neighbor(f) &&
                  cell->periodic_neighbor(f)->level() == cell->level())
                {
                  if (cell->is_locally_owned_on_level())
                    {
                      Assert(
                        cell->periodic_neighbor(f)->level_subdomain_id() !=
                          numbers::artificial_subdomain_id,
                        ExcMessage(
                          "Periodic neighbor of a locally owned cell must either be owned or ghost."));
                    }
                  // Cell is a level-ghost and its neighbor is a
                  // level-artificial cell nothing to do here
                  else if (cell->periodic_neighbor(f)->level_subdomain_id() ==
                           numbers::artificial_subdomain_id)
                    {
                      Assert(cell->is_locally_owned_on_level() == false,
                             ExcInternalError());
                      continue;
                    }

                  const unsigned int dofs_per_face =
                    cell->face(f)->get_fe(0).dofs_per_face;
                  std::vector<types::global_dof_index> dofs_1(dofs_per_face);
                  std::vector<types::global_dof_index> dofs_2(dofs_per_face);

                  cell->periodic_neighbor(f)
                    ->face(cell->periodic_neighbor_face_no(f))
                    ->get_mg_dof_indices(l, dofs_1, 0);
                  cell->face(f)->get_mg_dof_indices(l, dofs_2, 0);
                  // Store periodicity information in the level
                  // AffineConstraints object. Skip DoFs for which we've
                  // previously entered periodicity constraints already; this
                  // can happen, for example, for a vertex dof at a periodic
                  // boundary that we visit from more than one cell
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    if (level_constraints[l].can_store_line(dofs_2[i]) &&
                        level_constraints[l].can_store_line(dofs_1[i]) &&
                        !level_constraints[l].is_constrained(dofs_2[i]) &&
                        !level_constraints[l].is_constrained(dofs_1[i]))
                      {
                        level_constraints[l].add_line(dofs_2[i]);
                        level_constraints[l].add_entry(dofs_2[i],
                                                       dofs_1[i],
                                                       1.);
                      }
                }
          }
      level_constraints[l].close();

      // Initialize with empty IndexSet of correct size
      refinement_edge_indices[l] = IndexSet(dof.n_dofs(l));
    }

  MGTools::extract_inner_interface_dofs(dof, refinement_edge_indices);
}


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::initialize(
  const DoFHandler<dim, spacedim> &                               dof,
  const std::map<types::boundary_id, const Function<spacedim> *> &function_map,
  const ComponentMask &component_mask)
{
  initialize(dof);

  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  // At this point boundary_indices is empty.
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list(dof,
                              function_map,
                              boundary_indices,
                              component_mask);
}


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::make_zero_boundary_constraints(
  const DoFHandler<dim, spacedim> &   dof,
  const std::set<types::boundary_id> &boundary_ids,
  const ComponentMask &               component_mask)
{
  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  Assert(boundary_indices.size() == 0 || boundary_indices.size() == n_levels,
         ExcInternalError());
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list(dof,
                              boundary_ids,
                              boundary_indices,
                              component_mask);
}


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::make_no_normal_flux_constraints(
  const DoFHandler<dim, spacedim> &dof,
  const types::boundary_id         bid,
  const unsigned int               first_vector_component)
{
  // For a given boundary id, find which vector component is on the boundary
  // and set a zero boundary constraint for those degrees of freedom.
  const unsigned int n_components = DoFTools::n_components(dof);
  AssertIndexRange(first_vector_component + dim - 1, n_components);

  ComponentMask comp_mask(n_components, false);


  typename Triangulation<dim>::face_iterator
    face = dof.get_triangulation().begin_face(),
    endf = dof.get_triangulation().end_face();
  for (; face != endf; ++face)
    if (face->at_boundary() && face->boundary_id() == bid)
      for (unsigned int d = 0; d < dim; ++d)
        {
          Tensor<1, dim, double> unit_vec;
          unit_vec[d] = 1.0;

          const Tensor<1, dim> normal_vec =
            face->get_manifold().normal_vector(face, face->center());

          if (std::abs(std::abs(unit_vec * normal_vec) - 1.0) < 1e-10)
            comp_mask.set(d + first_vector_component, true);
          else
            Assert(
              std::abs(unit_vec * normal_vec) < 1e-10,
              ExcMessage(
                "We can currently only support no normal flux conditions "
                "for a specific boundary id if all faces are normal to the "
                "x, y, or z axis."));
        }

  Assert(comp_mask.n_selected_components() == 1,
         ExcMessage(
           "We can currently only support no normal flux conditions "
           "for a specific boundary id if all faces are facing in the "
           "same direction, i.e., a boundary normal to the x-axis must "
           "have a different boundary id than a boundary normal to the "
           "y- or z-axis and so on. If the mesh here was produced using "
           "GridGenerator::..., setting colorize=true during mesh generation "
           "and calling make_no_normal_flux_constraints() for each no normal "
           "flux boundary will fulfill the condition."));

  this->make_zero_boundary_constraints(dof, {bid}, comp_mask);
}


inline void
MGConstrainedDoFs::add_user_constraints(
  const unsigned int               level,
  const AffineConstraints<double> &constraints_on_level)
{
  AssertIndexRange(level, user_constraints.size());

  // Get the relevant DoFs from level_constraints if
  // the user constraint matrix has not been initialized
  if (user_constraints[level].get_local_lines().size() == 0)
    user_constraints[level].reinit(level_constraints[level].get_local_lines());

  user_constraints[level].merge(
    constraints_on_level,
    AffineConstraints<double>::MergeConflictBehavior::right_object_wins);
  user_constraints[level].close();
}


inline void
MGConstrainedDoFs::clear()
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  user_constraints.clear();
}


inline bool
MGConstrainedDoFs::is_boundary_index(const unsigned int            level,
                                     const types::global_dof_index index) const
{
  if (boundary_indices.size() == 0)
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
MGConstrainedDoFs::get_level_constraint_matrix(const unsigned int level) const
{
  return get_level_constraints(level);
}



inline const AffineConstraints<double> &
MGConstrainedDoFs::get_user_constraint_matrix(const unsigned int level) const
{
  AssertIndexRange(level, user_constraints.size());
  return user_constraints[level];
}



DEAL_II_NAMESPACE_CLOSE

#endif
