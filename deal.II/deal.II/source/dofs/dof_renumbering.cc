//----------------------------  dof_renumbering.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_renumbering.cc  ---------------------------


#include <lac/sparsity_pattern.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_handler.h>
#include <multigrid/mg_dof_tools.h>
#include <fe/fe.h>
#include <numerics/dof_renumbering.h>

#include <vector>
#include <map>
#include <algorithm>


template <int dim>
void DoFRenumbering::Cuthill_McKee (DoFHandler<dim>            &dof_handler,
				    const bool                  reversed_numbering,
				    const bool                  use_constraints,
				    const vector<unsigned int> &starting_indices)
{
				   // make the connection graph
  SparsityPattern sparsity (dof_handler.n_dofs(),
			    dof_handler.max_couplings_between_dofs());
  DoFTools::make_sparsity_pattern (dof_handler, sparsity);

  if (use_constraints) 
    {
      ConstraintMatrix constraints;
      DoFTools::make_hanging_node_constraints (dof_handler, constraints);
      constraints.close ();
      constraints.condense (sparsity);
    };
    
  const unsigned int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; invalid_dof_index means
				   // that no new number was chosen yet
  vector<unsigned int> new_number(sparsity.n_rows(), DoFHandler<dim>::invalid_dof_index);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<unsigned int> last_round_dofs (starting_indices);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]==DoFHandler<dim>::invalid_dof_index) ||
	(last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = DoFHandler<dim>::invalid_dof_index;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<unsigned int>(), DoFHandler<dim>::invalid_dof_index));
  
				   // now if no valid points remain:
				   // find dof with lowest coordination
				   // number
  
  if (last_round_dofs.size() == 0)
    {
      unsigned int starting_point   = DoFHandler<dim>::invalid_dof_index;
      unsigned int min_coordination = n_dofs;
      for (unsigned int row=0; row<n_dofs; ++row) 
	{
	  unsigned int j;

					   // loop until we hit the end
					   // of this row's entries
	  for (j=sparsity.get_rowstart_indices()[row];
	       j<sparsity.get_rowstart_indices()[row+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	      break;
					   // post-condition after loop:
					   // coordination, i.e. the number
					   // of entries in this row is now
					   // j-rowstart[row]
	  if (j-sparsity.get_rowstart_indices()[row] <  min_coordination)
	    {
	      min_coordination = j-sparsity.get_rowstart_indices()[row];
	      starting_point   = row;
	    };
	};
      
				       // now we still have to care for the
				       // case that no dof has a coordination
				       // number less than n_dofs. this rather
				       // exotic case only happens if we only
				       // have one cell, as far as I can see,
				       // but there may be others as well.
				       //
				       // if that should be the case, we can
				       // chose an arbitrary dof as starting
				       // point, e.g. the one with number zero
      if (starting_point == DoFHandler<dim>::invalid_dof_index)
	starting_point = 0;
      
				       // initialize the first dof
      last_round_dofs.push_back (starting_point);
    };


				   // store next free dof index
  unsigned int next_free_number = 0;

				   // enumerate the first round dofs
  for (unsigned int i=0; i!=last_round_dofs.size(); ++i)
    new_number[last_round_dofs[i]] = next_free_number++;

  bool all_dofs_renumbered = false;

				   // now do as many steps as needed to
				   // renumber all dofs
  while (!all_dofs_renumbered) 
    {
				       // store the indices of the dofs to be
				       // renumbered in the next round
      vector<unsigned int> next_round_dofs;

				       // find all neighbors of the
				       // dofs numbered in the last
				       // round
      for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	     j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	  if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
				       // sort dof numbers
      sort (next_round_dofs.begin(), next_round_dofs.end());

				       // delete multiple entries
      vector<unsigned int>::iterator end_sorted;
      end_sorted = unique (next_round_dofs.begin(), next_round_dofs.end());
      next_round_dofs.erase (end_sorted, next_round_dofs.end());

				       // eliminate dofs which are
				       // already numbered
      for (int s=next_round_dofs.size()-1; s>=0; --s)
	if (new_number[next_round_dofs[s]] != DoFHandler<dim>::invalid_dof_index)
	  next_round_dofs.erase (&next_round_dofs[s]);

				       // check whether there are any new
				       // dofs
      all_dofs_renumbered = (next_round_dofs.size() == 0);
      if (all_dofs_renumbered)
					 // end loop if possible
	continue;


				       // store for each coordination
				       // number the dofs with these
				       // coordination number
      multimap<unsigned int, int> dofs_by_coordination;
      
				       // find coordination number for
				       // each of these dofs
      for (vector<unsigned int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (unsigned int j=sparsity.get_rowstart_indices()[*s];
	       j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	      break;
	    else
	      ++coordination;

					   // insert this dof at its
					   // coordination number
	  const pair<const unsigned int, int> new_entry (coordination, *s);
	  dofs_by_coordination.insert (new_entry);
	};
      
				       // assign new DoF numbers to
				       // the elements of the present
				       // front:
      multimap<unsigned int, int>::iterator i;
      for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i) 
	new_number[i->second] = next_free_number++;

				       // after that: copy this round's
				       // dofs for the next round
      last_round_dofs = next_round_dofs;
    };

//TODO: Allow incomplete renumbering for non-discretization values

#ifdef DEBUG
				   // test for all indices
				   // numbered. this mostly tests
				   // whether the
				   // front-marching-algorithm (which
				   // Cuthill-McKee actually is) has
				   // reached all points. it should
				   // usually do so, but might not for
				   // two reasons:
				   //
				   // - The algorithm above has a bug, or
				   // - The domain is not connected and
				   // consists of separate parts.
				   //
				   // In any case, if not all DoFs
				   // have been reached, renumbering
				   // will not be possible
  if (find (new_number.begin(), new_number.end(), DoFHandler<dim>::invalid_dof_index)
      !=
      new_number.end())
    Assert (false, ExcRenumberingIncomplete());
  Assert (next_free_number == n_dofs,
	  ExcRenumberingIncomplete());
#endif

  if (reversed_numbering)
    for (vector<unsigned int>::iterator i=new_number.begin(); i!=new_number.end(); ++i)
      *i = n_dofs-*i;

				   // actually perform renumbering;
				   // this is dimension specific and
				   // thus needs an own function
  dof_handler.renumber_dofs (new_number);
};


template <int dim>
void DoFRenumbering::Cuthill_McKee (MGDoFHandler<dim>          &dof_handler,
				    const unsigned int          level,
				    const bool                  reversed_numbering,
				    const vector<unsigned int> &starting_indices) {
				   // make the connection graph
  SparsityPattern sparsity (dof_handler.n_dofs(level),
			    dof_handler.max_couplings_between_dofs());
  MGDoFTools::make_sparsity_pattern (dof_handler, sparsity, level);
    
  const unsigned int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; invalid_dof_index means
				   // that no new number was chosen yet
  vector<unsigned int> new_number(n_dofs, DoFHandler<dim>::invalid_dof_index);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<unsigned int> last_round_dofs (starting_indices);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]==DoFHandler<dim>::invalid_dof_index) ||
	(last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = DoFHandler<dim>::invalid_dof_index;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<unsigned int>(), DoFHandler<dim>::invalid_dof_index));
  
				   // now if no valid points remain:
				   // find dof with lowest coordination
				   // number
  
  if (last_round_dofs.size() == 0)
    {
      unsigned int starting_point   = DoFHandler<dim>::invalid_dof_index;
      unsigned int min_coordination = n_dofs;
      for (unsigned int row=0; row<n_dofs; ++row) 
	{
	  unsigned int j;
	  for (j=sparsity.get_rowstart_indices()[row];
	       j<sparsity.get_rowstart_indices()[row+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	      break;
					   // post-condition after loop:
					   // coordination is now
					   // j-rowstart[row]
	  if (j-sparsity.get_rowstart_indices()[row] <  min_coordination)
	    {
	      min_coordination = j-sparsity.get_rowstart_indices()[row];
	      starting_point   = row;
	    };
	};
				       // initialize the first dof
      last_round_dofs.push_back (starting_point);
    };


				   // store next free dof index
  unsigned int next_free_number = 0;

				   // enumerate the first round dofs
  for (unsigned int i=0; i!=last_round_dofs.size(); ++i)
    new_number[last_round_dofs[i]] = next_free_number++;

  bool all_dofs_renumbered = false;

				   // now do as many steps as needed to
				   // renumber all dofs
  while (!all_dofs_renumbered) 
    {
				       // store the indices of the dofs to be
				       // renumbered in the next round
      vector<unsigned int> next_round_dofs;

				       // find all neighbors of the
				       // dofs numbered in the last
				       // round
      for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	     j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	  if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
				       // sort dof numbers
      sort (next_round_dofs.begin(), next_round_dofs.end());

				       // delete multiple entries
      vector<unsigned int>::iterator end_sorted;
      end_sorted = unique (next_round_dofs.begin(), next_round_dofs.end());
      next_round_dofs.erase (end_sorted, next_round_dofs.end());

				       // eliminate dofs which are
				       // already numbered
      for (int s=next_round_dofs.size()-1; s>=0; --s)
	if (new_number[next_round_dofs[s]] != DoFHandler<dim>::invalid_dof_index)
	  next_round_dofs.erase (&next_round_dofs[s]);

				       // check whether there are any new
				       // dofs
      all_dofs_renumbered = (next_round_dofs.size() == 0);
      if (all_dofs_renumbered)
					 // end loop if possible
	continue;


				       // store for each coordination
				       // number the dofs with these
				       // coordination number
      multimap<unsigned int, int> dofs_by_coordination;
      
				       // find coordination number for
				       // each of these dofs
      for (vector<unsigned int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (unsigned int j=sparsity.get_rowstart_indices()[*s];
	       j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	      break;
	    else
	      ++coordination;

					   // insert this dof at its
					   // coordination number
	  const pair<const unsigned int, int> new_entry (coordination, *s);
	  dofs_by_coordination.insert (new_entry);
	};
      
				       ////
      multimap<unsigned int, int>::iterator i;
      for (i = dofs_by_coordination.begin(); i!=dofs_by_coordination.end(); ++i) 
	new_number[i->second] = next_free_number++;

				       // after that: copy this round's
				       // dofs for the next round
      last_round_dofs = next_round_dofs;
    };

#ifdef DEBUG
				   //  test for all indices numbered
  if (find (new_number.begin(), new_number.end(),
	    DoFHandler<dim>::invalid_dof_index)
      !=
      new_number.end())
    Assert (false, ExcRenumberingIncomplete());
  Assert (next_free_number == n_dofs,
	  ExcRenumberingIncomplete());
#endif

  if (reversed_numbering)
    for (vector<unsigned int>::iterator i=new_number.begin(); i!=new_number.end(); ++i)
      *i = n_dofs-*i;

				   // actually perform renumbering;
				   // this is dimension specific and
				   // thus needs an own function
  dof_handler.renumber_dofs (level, new_number);
};



template <int dim>
void DoFRenumbering::component_wise (DoFHandler<dim>            &dof_handler,
				     const vector<unsigned int> &component_order_arg)
{
  const unsigned int dofs_per_cell = dof_handler.get_fe().dofs_per_cell;

				   // do nothing if the FE has only
				   // one component
  if (dof_handler.get_fe().n_components() == 1)
    return;
  
  vector<unsigned int> component_order (component_order_arg);
  if (component_order.size() == 0)
    for (unsigned int i=0; i<dof_handler.get_fe().n_components(); ++i)
      component_order.push_back (i);

				   // check whether the component list has
				   // the right length and contains all
				   // component numbers
  Assert (component_order.size() == dof_handler.get_fe().n_components(),
	  ExcInvalidComponentOrder());
  for (unsigned int i=0; i<dof_handler.get_fe().n_components(); ++i)
    Assert (find (component_order.begin(), component_order.end(), i)
	    != component_order.end(),
	    ExcInvalidComponentOrder ());

				   // vector to hold the dof indices on
				   // the cell we visit at a time
  vector<unsigned int> local_dof_indices(dofs_per_cell);
				   // prebuilt list to which component
				   // a given dof on a cell belongs
  vector<unsigned int> component_list (dofs_per_cell);
  for (unsigned int i=0; i<component_list.size(); ++i)
    component_list[i] = dof_handler.get_fe().system_to_component_index(i).first;
  
				   // set up a map where for each component
				   // the respective degrees of freedom are
				   // collected.
				   //
				   // note that this map is sorted by component
				   // but that within each component it is NOT
				   // sorted by dof index. note also that some
				   // dof indices are entered multiply, so we
				   // will have to take care of that
  vector<vector<unsigned int> > component_to_dof_map (dof_handler.get_fe().n_components());
  for (typename DoFHandler<dim>::active_cell_iterator cell=dof_handler.begin_active();
       cell!=dof_handler.end(); ++cell)
    {
				       // on each cell: get dof indices
				       // and insert them into the global
				       // list using their component
      cell->get_dof_indices (local_dof_indices);
      for (unsigned int i=0; i<dofs_per_cell; ++i)
	component_to_dof_map[component_list[i]].push_back (local_dof_indices[i]);
    };
  
				   // now we've got all indices sorted into
				   // buckets labelled with their component
				   // number. we've only got to traverse this
				   // list and assign the new indices
				   //
				   // however, we first want to sort the
				   // indices entered into the buckets to
				   // preserve the order within each component
				   // and during this also remove duplicate
				   // entries
  for (unsigned int component=0; component<dof_handler.get_fe().n_components(); ++component)
    {
      sort (component_to_dof_map[component].begin(),
	    component_to_dof_map[component].end());
      component_to_dof_map[component].erase (unique (component_to_dof_map[component].begin(),
						     component_to_dof_map[component].end()),
					     component_to_dof_map[component].end());
    };
  
  unsigned int next_free_index = 0;
  vector<unsigned int> new_indices (dof_handler.n_dofs(), DoFHandler<dim>::invalid_dof_index);
  for (unsigned int component=0; component<dof_handler.get_fe().n_components(); ++component)
    {
      const typename vector<unsigned int>::const_iterator
	begin_of_component = component_to_dof_map[component].begin(),
	end_of_component   = component_to_dof_map[component].end();
            
      for (typename vector<unsigned int>::const_iterator dof_index = begin_of_component;
	   dof_index != end_of_component; ++dof_index)
	new_indices[*dof_index] = next_free_index++;
    };

  Assert (next_free_index == dof_handler.n_dofs(),
	  ExcInternalError());

  dof_handler.renumber_dofs (new_indices);
};



template <int dim>
void
DoFRenumbering::sort_selected_dofs_back (const vector<bool> &selected_dofs,
					 DoFHandler<dim>    &dof_handler)
{
  const unsigned int n_dofs = dof_handler.n_dofs();
  Assert (selected_dofs.size() == n_dofs,
	  ExcInvalidArraySize (selected_dofs.size(), n_dofs));

				   // re-sort the dofs according to
				   // their selection state
  vector<unsigned int> new_dof_indices (n_dofs);
  const unsigned int   n_selected_dofs = count (selected_dofs.begin(),
						selected_dofs.end(),
						false);
  
  unsigned int next_unselected = 0;
  unsigned int next_selected   = n_selected_dofs;
  for (unsigned int i=0; i<n_dofs; ++i)
    if (selected_dofs[i] == false)
      {
	new_dof_indices[i] = next_unselected;
	++next_unselected;
      }
    else
      {
	new_dof_indices[i] = next_selected;
	++next_selected;
      };
  Assert (next_unselected == n_selected_dofs, ExcInternalError());
  Assert (next_selected == n_dofs, ExcInternalError());

				   // now perform the renumbering
  dof_handler.renumber_dofs (new_dof_indices);
};





// explicit instantiations
template
void DoFRenumbering::Cuthill_McKee (DoFHandler<deal_II_dimension> &dof_handler,
				    const bool                     reversed_numbering,
				    const bool                     use_constraints,
				    const vector<unsigned int>    &starting_indices);

template
void DoFRenumbering::Cuthill_McKee (MGDoFHandler<deal_II_dimension> &dof_handler,
				    const unsigned int               level,
				    const bool                       reversed_numbering,
				    const vector<unsigned int>      &starting_indices);

template
void DoFRenumbering::component_wise (DoFHandler<deal_II_dimension> &dof_handler,
				     const vector<unsigned int>    &component_order_arg);


template
void DoFRenumbering::sort_selected_dofs_back (const vector<bool> &selected_dofs,
					      DoFHandler<deal_II_dimension> &dof_handler);
