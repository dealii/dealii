// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2014 by the deal.II authors
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

#ifndef __deal2__mg_constraints_h
#define __deal2__mg_constraints_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>

#include <deal.II/multigrid/mg_tools.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class DoFHandler;
template <int dim, typename Number> struct FunctionMap;


/**
 * Collection of boundary constraints and refinement edge constraints
 * for level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public Subscriptor
{
public:

  typedef std::vector<std::set<types::global_dof_index> >::size_type size_dof;
  /**
   * Fill the internal data structures with hanging node constraints
   * extracted from the dof handler object. Works with natural
   * boundary conditions only. There exists a sister function setting
   * up boundary constraints as well.
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
   * Determine whether a dof index is subject to a boundary
   * constraint.
   */
  bool is_boundary_index (const unsigned int level,
                          const types::global_dof_index index) const;

  /**
   * @deprecated Determine whether a dof index is at an edge that is not a
   * refinement edge.
   */
  bool non_refinement_edge_index (const unsigned int level,
                                  const types::global_dof_index index) const DEAL_II_DEPRECATED;

  /**
   * Determine whether a dof index is at the refinement edge.
   */
  bool at_refinement_edge (const unsigned int level,
                           const types::global_dof_index index) const;

  /**
   * @deprecated Use is_boundary_index() instead. The logic behind
   * this function here is unclear and for practical purposes, the
   * other is needed.
   *
   * Determine whether a dof index is subject to a boundary
   * constraint.
   */
  bool at_refinement_edge_boundary (const unsigned int level,
                                    const types::global_dof_index index) const DEAL_II_DEPRECATED;

  /**
   * Return the indices of dofs for each level that are subject to
   * boundary constraints.
   */
  const std::vector<std::set<types::global_dof_index> > &
  get_boundary_indices () const;

  /**
   * @deprecated Use at_refinement_edge() if possible, else
   * get_refinement_edge_indices(unsigned int).
   *
   * Return the indices of dofs for each level that lie on the
   * refinement edge (i.e. are on faces between cells of this level
   * and cells on the level below).
   */
  const std::vector<std::vector<bool> > &
  get_refinement_edge_indices () const DEAL_II_DEPRECATED;

  /**
   * Return the indices of dofs on the given level that lie on an
   * refinement edge (dofs on faces to neighbors that are coarser)
   */
  const IndexSet &
  get_refinement_edge_indices (unsigned int level) const;

  /**
   * @deprecated Use at_refinement_edge_boundary() if possible, else
   * use get_refinement_edge_boundary_indices().
   *
   * Return the indices of dofs for each level that are in the
   * intersection of the sets returned by get_boundary_indices() and
   * get_refinement_edge_indices().
   */
  const std::vector<std::vector<bool> > &
  get_refinement_edge_boundary_indices () const DEAL_II_DEPRECATED;

  /**
   * @deprecated The function is_boundary_index() now returns false if
   * no boundary values are set.
   *
   * Return if boundary_indices need to be set or not.
   */
  bool set_boundary_values () const DEAL_II_DEPRECATED;

private:

  /**
   * The indices of boundary dofs for each level.
   */
  std::vector<std::set<types::global_dof_index> > boundary_indices;

  /**
   * The degrees of freedom on the refinement edge between a level and
   * coarser cells.
   */
  std::vector<IndexSet> refinement_edge_indices;

  /**
   * old data structure only filled on demand
   */
  mutable std::vector<std::vector<bool> > refinement_edge_boundary_indices_old;

  /**
   * old data structure only filled on demand
   */
  mutable std::vector<std::vector<bool> > refinement_edge_indices_old;
};


template <int dim, int spacedim>
inline
void
MGConstrainedDoFs::initialize(const DoFHandler<dim,spacedim> &dof)
{
  const unsigned int nlevels = dof.get_tria().n_global_levels();

  boundary_indices.resize(nlevels);

  refinement_edge_indices.resize(nlevels);
  refinement_edge_indices_old.clear();
  refinement_edge_boundary_indices_old.clear();
  for (unsigned int l=0; l<nlevels; ++l)
    {
      boundary_indices[l].clear();

      refinement_edge_indices[l] = IndexSet(dof.n_dofs(l));
    }

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

  MGTools::make_boundary_list (dof,
                               function_map,
                               boundary_indices,
                               component_mask);
}


inline
void
MGConstrainedDoFs::clear()
{
  for (unsigned int l=0; l<boundary_indices.size(); ++l)
    boundary_indices[l].clear();

  for (unsigned int l=0; l<refinement_edge_indices.size(); ++l)
    refinement_edge_indices[l].clear();

  refinement_edge_indices_old.clear();
  refinement_edge_boundary_indices_old.clear();
}


inline
bool
MGConstrainedDoFs::is_boundary_index (const unsigned int level,
                                      const types::global_dof_index index) const
{
  if (boundary_indices.size() == 0)
    return false;

  AssertIndexRange(level, boundary_indices.size());
  return (boundary_indices[level].find(index) != boundary_indices[level].end());
}

inline
bool
MGConstrainedDoFs::non_refinement_edge_index (const unsigned int level,
                                              const types::global_dof_index index) const
{
  return !at_refinement_edge (level, index);
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
bool
MGConstrainedDoFs::at_refinement_edge_boundary (const unsigned int level,
                                                const types::global_dof_index index) const
{
  return is_boundary_index(level, index);
}

inline
const std::vector<std::set<types::global_dof_index> > &
MGConstrainedDoFs::get_boundary_indices () const
{
  return boundary_indices;
}

inline
const std::vector<std::vector<bool> > &
MGConstrainedDoFs::get_refinement_edge_indices () const
{
  if (refinement_edge_indices_old.size()!=refinement_edge_indices.size())
    {
      unsigned int n_levels = refinement_edge_indices.size();
      refinement_edge_indices_old.resize(n_levels);
      for (unsigned int l=0; l<n_levels; ++l)
        {
          refinement_edge_indices_old[l].resize(refinement_edge_indices[l].size(), false);
          refinement_edge_indices[l].fill_binary_vector(refinement_edge_indices_old[l]);
        }
    }

  return refinement_edge_indices_old;
}

inline
const IndexSet &
MGConstrainedDoFs::get_refinement_edge_indices (unsigned int level) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  return refinement_edge_indices[level];
}


inline
const std::vector<std::vector<bool> > &
MGConstrainedDoFs::get_refinement_edge_boundary_indices () const
{
  if (refinement_edge_boundary_indices_old.size()==0)
    {
      unsigned int n_levels = refinement_edge_indices.size();
      refinement_edge_boundary_indices_old.resize(n_levels);
      for (unsigned int l=0; l<n_levels; ++l)
        {
          refinement_edge_boundary_indices_old[l].resize(refinement_edge_indices[l].size());
          for (types::global_dof_index idx=0; idx<refinement_edge_indices[l].size(); ++idx)
            refinement_edge_boundary_indices_old[l][idx] = this->is_boundary_index(l, idx);
        }
    }

  return refinement_edge_boundary_indices_old;
}


inline
bool
MGConstrainedDoFs::set_boundary_values () const
{
  const bool boundary_values_need_to_be_set
    = boundary_indices.size()!=0;
  return boundary_values_need_to_be_set;
}


DEAL_II_NAMESPACE_CLOSE

#endif
