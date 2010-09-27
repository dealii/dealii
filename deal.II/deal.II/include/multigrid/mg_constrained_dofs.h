//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_constraints_h
#define __deal2__mg_constraints_h

#include <base/config.h>
#include <base/subscriptor.h>

#include <multigrid/mg_tools.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class MGDoFHandler;
template <int dim> class FunctionMap;


/**
 * Collection of boundary constraints and refinement edge constraints
 * for level vectors.
 *
 * @ingroup mg
 */
class MGConstrainedDoFs : public Subscriptor
{
  public:
				     /**
				      * Fill the internal data
				      * structures with values
				      * extracted from the dof
				      * handler.
				      *
				      * This function leaves
				      * #boundary_indices empty, since
				      * no boundary values are
				      * provided.
				      */
    template <int dim, int spacedim>
    void initialize(const MGDoFHandler<dim,spacedim>& dof);

				     /**
				      * Fill the internal data
				      * structures with values
				      * extracted from the dof
				      * handler, applying the boundary
				      * values provided.
				      */
    template <int dim, int spacedim>
    void initialize(const MGDoFHandler<dim,spacedim>& dof,
		    const typename FunctionMap<dim>::type& function_map,
		    const std::vector<bool>& component_mask = std::vector<bool>());

				     /**
				      * Reset the data structures.
				      */
    void clear();

				     /**
				      * Determine whether a dof index
				      * is subject to a boundary
				      * constraint.
				      */
    bool is_boundary_index (const unsigned int level,
			    const unsigned int index) const;

				     /**
				      * Determine whether a dof index
				      * is at the refinement edge.
				      */
    bool at_refinement_edge (const unsigned int level,
			     const unsigned int index) const;

				     /**
				      * Determine whether a dof index
				      * is at the refinement edge and 
                                      * subject to a boundary
				      * constraint .
				      */
    bool at_refinement_edge_boundary (const unsigned int level,
				      const unsigned int index) const;

				     /**
				      * Return the indices of dofs for each
				      * level that lie on the boundary of the
				      * domain.
				      */
    const std::vector<std::set<unsigned int> > &
      get_boundary_indices () const;

				     /**
				      * Return the indices of dofs for each
				      * level that lie on the refinement edge
				      * (i.e. are on faces between cells of
				      * this level and cells on the level
				      * below).
				      */
    const std::vector<std::vector<bool> > &
      get_refinement_edge_indices () const;

				     /**
				      * Return the indices of dofs for each
				      * level that are in the intersection of
				      * the sets returned by
				      * get_boundary_indices() and
				      * get_refinement_edge_indices().
				      */
    const std::vector<std::vector<bool> > &
      get_refinement_edge_boundary_indices () const;

				     /**
				      * Return if boundary_indices need to 
                                      * be set or not.
				      */

    bool set_boundary_values () const;
    bool continuity_across_edges () const;
  private:

				     /**
				      * The indices of boundary dofs
				      * for each level.
				      */
    std::vector<std::set<unsigned int> > boundary_indices;

				     /**
				      * The degrees of freedom on the
				      * refinement edge between a
				      * level and coarser cells.
				      */
    std::vector<std::vector<bool> > refinement_edge_indices;

				     /**
				      * The degrees of freedom on the
				      * refinement edge between a
				      * level and coarser cells, which
				      * are also on the boundary.
				      *
				      * This is a subset of
				      * #refinement_edge_indices.
				      */
    std::vector<std::vector<bool> > refinement_edge_boundary_indices;
};


template <int dim, int spacedim>
inline
void
MGConstrainedDoFs::initialize(const MGDoFHandler<dim,spacedim>& dof)
{
  const unsigned int nlevels = dof.get_tria().n_levels();
  refinement_edge_indices.resize(nlevels);
  refinement_edge_boundary_indices.resize(nlevels);
  for(unsigned int l=0; l<nlevels; ++l)
    {
      refinement_edge_indices[l].resize(dof.n_dofs(l));
      refinement_edge_boundary_indices[l].resize(dof.n_dofs(l));
    }
  MGTools::extract_inner_interface_dofs (dof, refinement_edge_indices, 
					 refinement_edge_boundary_indices);
}


template <int dim, int spacedim>
inline
void
MGConstrainedDoFs::initialize(
  const MGDoFHandler<dim,spacedim>& dof,
  const typename FunctionMap<dim>::type& function_map,
  const std::vector<bool>& component_mask)
{
  const unsigned int nlevels = dof.get_tria().n_levels();
  boundary_indices.resize(nlevels);
  refinement_edge_indices.resize(nlevels);
  refinement_edge_boundary_indices.resize(nlevels);
  
  for(unsigned int l=0; l<nlevels; ++l)
    {
      boundary_indices[l].clear();
      refinement_edge_indices[l].resize(dof.n_dofs(l));
      refinement_edge_boundary_indices[l].resize(dof.n_dofs(l));
    }
  
  MGTools::make_boundary_list (dof, function_map, boundary_indices, component_mask);
  MGTools::extract_inner_interface_dofs (dof, refinement_edge_indices, 
					 refinement_edge_boundary_indices);
}


inline
void
MGConstrainedDoFs::clear() 
{
  for(unsigned int l=0; l<boundary_indices.size(); ++l)
    boundary_indices[l].clear();
  
  for(unsigned int l=0; l<refinement_edge_indices.size(); ++l)
    refinement_edge_indices[l].clear();
  
  for(unsigned int l=0; l<refinement_edge_boundary_indices.size(); ++l)
    refinement_edge_boundary_indices[l].clear();
}


inline
bool
MGConstrainedDoFs::is_boundary_index (const unsigned int level, 
				  const unsigned int index) const
{
  AssertIndexRange(level, boundary_indices.size());
  if(boundary_indices[level].find(index) != boundary_indices[level].end())
    return true;
  else
    return false;
}


inline
bool
MGConstrainedDoFs::at_refinement_edge (const unsigned int level,
				   const unsigned int index) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  AssertIndexRange(index, refinement_edge_indices[level].size());

  return refinement_edge_indices[level][index];
}


inline
bool
MGConstrainedDoFs::at_refinement_edge_boundary (const unsigned int level,
					    const unsigned int index) const
{
  AssertIndexRange(level, refinement_edge_boundary_indices.size());
  AssertIndexRange(index, refinement_edge_boundary_indices[level].size());

  return refinement_edge_boundary_indices[level][index];
}

inline
const std::vector<std::set<unsigned int> > &
MGConstrainedDoFs::get_boundary_indices () const
{
  return boundary_indices;
}

inline
const std::vector<std::vector<bool> > &
MGConstrainedDoFs::get_refinement_edge_indices () const
{
  return refinement_edge_indices;
}

inline
const std::vector<std::vector<bool> > &
MGConstrainedDoFs::get_refinement_edge_boundary_indices () const
{
  return refinement_edge_boundary_indices;
}

inline
bool 
MGConstrainedDoFs::set_boundary_values () const
{
  const bool boundary_values_need_to_be_set 
    = boundary_indices.size()!=0;
  return boundary_values_need_to_be_set;
}

inline
bool 
MGConstrainedDoFs::continuity_across_edges () const
{
  bool is_continuous = false;
  for(unsigned int l=0; l<refinement_edge_indices.size(); ++l)
    for(unsigned int i=0; i<refinement_edge_indices[l].size(); ++i)
      if(refinement_edge_indices[l][i])
      {
        is_continuous = true;
        return is_continuous;
      }
  return is_continuous;
}

DEAL_II_NAMESPACE_CLOSE

#endif
