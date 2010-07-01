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
 */
class MGConstraints : public Subscriptor
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
    bool is_boundary_index(unsigned int level, unsigned int index) const;

				     /**
				      * Determine whether a dof index
				      * is at the refinement edge.
				      */
    bool at_refinement_edge(unsigned int level, unsigned int index) const;

				     /**
				      * Determine whether a dof index
				      * is at the refinement edge and 
                                      * subject to a boundary
				      * constraint .
				      */
    bool at_refinement_edge_boundary(unsigned int level, unsigned int index) const;

				     /**
				      * The indices of boundary dofs
				      * for each level.
				      */
    std::vector<std::set<unsigned int> > boundary_indices;

  private:

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
void MGConstraints::initialize(const MGDoFHandler<dim,spacedim>& dof)
{
  MGTools::extract_inner_interface_dofs (dof, refinement_edge_indices, 
      refinement_edge_boundary_indices);
}

template <int dim, int spacedim>
void MGConstraints::initialize(const MGDoFHandler<dim,spacedim>& dof,
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

void
MGConstraints::clear() 
{
  for(unsigned int l=0; l<boundary_indices.size(); ++l)
    boundary_indices[l].clear();

  for(unsigned int l=0; l<refinement_edge_indices.size(); ++l)
    refinement_edge_indices[l].clear();

  for(unsigned int l=0; l<refinement_edge_boundary_indices.size(); ++l)
    refinement_edge_boundary_indices[l].clear();
}

bool
MGConstraints::is_boundary_index(unsigned int level, unsigned int index) const
{
  AssertIndexRange(level, boundary_indices.size());
  if(boundary_indices[level].find(index) != boundary_indices[level].end())
    return true;
  else
    return false;
}

bool
MGConstraints::at_refinement_edge(unsigned int level, unsigned int index) const
{
  AssertIndexRange(level, refinement_edge_indices.size());
  AssertIndexRange(index, refinement_edge_indices[level].size());

  return refinement_edge_indices[level][index];
}

bool
MGConstraints::at_refinement_edge_boundary(unsigned int level, unsigned int index) const
{
  AssertIndexRange(level, refinement_edge_boundary_indices.size());
  AssertIndexRange(index, refinement_edge_boundary_indices[level].size());

  return refinement_edge_boundary_indices[level][index];
}


DEAL_II_NAMESPACE_CLOSE

#endif
