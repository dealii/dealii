// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

#ifndef dealii__mg_constrained_dofs_h
#define dealii__mg_constrained_dofs_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/multigrid/mg_tools.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class DoFHandler;
template <int dim, typename Number> struct FunctionMap;


/**
 * Collection of boundary constraints and refinement edge constraints for
 * level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public Subscriptor
{
public:

  typedef std::vector<std::set<types::global_dof_index> >::size_type size_dof;
  /**
   * Fill the internal data structures with hanging node constraints extracted
   * from the dof handler object. Works with natural boundary conditions only.
   * There exists a sister function setting up boundary constraints as well.
   *
   * This function ensures that on every level, degrees of freedom at interior
   * edges of a refinement level are treated corrected but leaves degrees of
   * freedom at the boundary of the domain untouched assuming that no
   * Dirichlet boundary conditions for them exist.
   */
  template <int dim, int spacedim>
  void initialize (const DoFHandler<dim,spacedim> &dof);

  /**
   * Fill the internal data structures with values extracted from the dof
   * handler object and apply the boundary values provided.
   *
   * This function internally calls the initialize() function above and the
   * constrains degrees on the external boundary of the domain by calling
   * MGTools::make_boundary_list() with the given second and third argument.
   */
  template <int dim, int spacedim>
  void initialize(const DoFHandler<dim,spacedim> &dof,
                  const typename FunctionMap<dim>::type &function_map,
                  const ComponentMask &component_mask = ComponentMask());

  /**
   * Reset the data structures.
   */
  void clear();

  /**
   * Determine whether a dof index is subject to a boundary constraint.
   */
  bool is_boundary_index (const unsigned int level,
                          const types::global_dof_index index) const;

  /**
   * Determine whether a dof index is at the refinement edge.
   */
  bool at_refinement_edge (const unsigned int level,
                           const types::global_dof_index index) const;

  /**
   * Return the indices of level dofs on the given level that are subject to
   * Dirichlet boundary conditions (as set by the @p function_map parameter in
   * initialize()).  The indices are restricted to the set of locally relevant
   * level dofs.
   */
  const IndexSet &
  get_boundary_indices (const unsigned int level) const;


  /**
   * Return the indices of dofs on the given level that lie on an refinement
   * edge (dofs on faces to neighbors that are coarser).
   */
  const IndexSet &
  get_refinement_edge_indices (unsigned int level) const;


  /**
   * Return if Dirichlet boundary indices are set in initialize().
   */
  bool have_boundary_indices () const;

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
};


template <int dim, int spacedim>
inline
void
MGConstrainedDoFs::initialize(const DoFHandler<dim,spacedim> &dof)
{
  boundary_indices.clear();

  const unsigned int nlevels = dof.get_triangulation().n_global_levels();

  refinement_edge_indices.resize(nlevels);
  for (unsigned int l=0; l<nlevels; ++l)
    refinement_edge_indices[l] = IndexSet(dof.n_dofs(l));

  MGTools::extract_inner_interface_dofs (dof, refinement_edge_indices);
}


template <int dim, int spacedim>
inline
void
MGConstrainedDoFs::initialize(const DoFHandler<dim,spacedim> &dof,
                              const typename FunctionMap<dim>::type &function_map,
                              const ComponentMask &component_mask)
{
  initialize (dof);

  // allocate an IndexSet for each global level. Contents will be
  // overwritten inside make_boundary_list.
  const unsigned int n_levels = dof.get_triangulation().n_global_levels();
  boundary_indices.resize(n_levels);

  MGTools::make_boundary_list (dof,
                               function_map,
                               boundary_indices,
                               component_mask);
}


inline
void
MGConstrainedDoFs::clear()
{
  boundary_indices.clear();
  refinement_edge_indices.clear();
}


inline
bool
MGConstrainedDoFs::is_boundary_index (const unsigned int level,
                                      const types::global_dof_index index) const
{
  if (boundary_indices.size() == 0)
    return false;

  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level].is_element(index);
}

inline
bool
MGConstrainedDoFs::at_refinement_edge (const unsigned int level,
                                       const types::global_dof_index index) const
{
  AssertIndexRange(level, refinement_edge_indices.size());

  return refinement_edge_indices[level].is_element(index);
}




inline
const IndexSet &
MGConstrainedDoFs::get_boundary_indices (const unsigned int level) const
{
  AssertIndexRange(level, boundary_indices.size());
  return boundary_indices[level];
}



inline
const IndexSet &
MGConstrainedDoFs::get_refinement_edge_indices (unsigned int level) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  return refinement_edge_indices[level];
}




inline
bool
MGConstrainedDoFs::have_boundary_indices () const
{
  return boundary_indices.size()!=0;
}


DEAL_II_NAMESPACE_CLOSE

#endif
