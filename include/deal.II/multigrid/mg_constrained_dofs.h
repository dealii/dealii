// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_mg_constrained_dofs_h
#define dealii_mg_constrained_dofs_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <deal.II/multigrid/mg_tools.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
class DoFHandler;
template <int dim, typename Number>
struct FunctionMap;


/**
 * Collection of boundary constraints and refinement edge constraints for
 * level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public Subscriptor
{
public:
  typedef std::vector<std::set<types::global_dof_index>>::size_type size_dof;
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
   * Furthermore, this call sets up a ConstraintMatrix on each level that
   * contains possible periodicity constraints in case those have been added to
   * the underlying triangulation. The constraint matrix can be queried by
   * get_level_constraint_matrix(level). Note that the current implementation of
   * periodicity constraints in this class does not support rotation matrices in
   * the periodicity definition, i.e., the respective argument in the
   * GridTools::collect_periodic_faces() may not be different from the identity
   * matrix.
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
   * @deprecated Use initialize() followed by make_zero_boundary_constraints()
   * instead
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED void
  initialize(const DoFHandler<dim, spacedim> &      dof,
             const typename FunctionMap<dim>::type &function_map,
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
   * Return the level constraint matrix for a given level, containing
   * periodicity constraints (if enabled on the triangulation).
   */
  const ConstraintMatrix &
  get_level_constraint_matrix(const unsigned int level) const;

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
  std::vector<ConstraintMatrix> level_constraints;
};


template <int dim, int spacedim>
inline void
MGConstrainedDoFs::initialize(const DoFHandler<dim, spacedim> &dof)
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
  level_constraints.clear();

  const unsigned int nlevels = dof.get_triangulation().n_global_levels();

  // At this point level_constraint and refinement_edge_indices are empty.
  level_constraints.resize(nlevels);
  refinement_edge_indices.resize(nlevels);
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
            for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
              if (cell->has_periodic_neighbor(f))
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
                  // Store periodicity information in the level constraint
                  // matrix Skip DoFs for which we've previously entered
                  // periodicity constraints already; this can happen, for
                  // example, for a vertex dof at a periodic boundary that we
                  // visit from more than one cell
                  for (unsigned int i = 0; i < dofs_per_face; ++i)
                    if (level_constraints[l].can_store_line(dofs_2[i]) &&
                        level_constraints[l].can_store_line(dofs_1[i]) &&
                        !level_constraints[l].is_constrained(dofs_2[i]) &&
                        !level_constraints[l].is_constrained(dofs_1[i]))
                      {
                        level_constraints[l].add_line(dofs_2[i]);
                        level_constraints[l].add_entry(
                          dofs_2[i], dofs_1[i], 1.);
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
  const DoFHandler<dim, spacedim> &      dof,
  const typename FunctionMap<dim>::type &function_map,
  const ComponentMask &                  component_mask)
{
  initialize(dof);

  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  // At this point boundary_indices is empty.
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list(
    dof, function_map, boundary_indices, component_mask);
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

  MGTools::make_boundary_list(
    dof, boundary_ids, boundary_indices, component_mask);
}


inline void
MGConstrainedDoFs::clear()
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
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



inline const ConstraintMatrix &
MGConstrainedDoFs::get_level_constraint_matrix(const unsigned int level) const
{
  AssertIndexRange(level, level_constraints.size());
  return level_constraints[level];
}



DEAL_II_NAMESPACE_CLOSE

#endif
