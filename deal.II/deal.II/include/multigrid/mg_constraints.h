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

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

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
    bool at_boundary(unsigned int level, unsigned int index) const;
    
				     /**
				      * Determine whether a dof index
				      * is at the refinement edge.
				      */
    bool at_refinement_edge(unsigned int level, unsigned int index) const;
    
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


bool
MGConstraints::at_boundary(unsigned int level, unsigned int index)
{
  AssertIndexRange(level, boundary_indices.size());
  AssertIndexRange(level, boundary_indices.size());  
}


DEAL_II_NAMESPACE_CLOSE

#endif
