/* $Id$ */

#include <lac/sparsity_pattern.h>
#include <dofs/dof_accessor.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_handler.h>
#include <dofs/mg_dof_handler.h>
#include <dofs/dof_constraints.h>
#include <dofs/dof_tools.h>
#include <dofs/mg_dof_tools.h>
#include <fe/fe.h>
#include <numerics/dof_renumbering.h>

#include <vector>
#include <map>
#include <algorithm>



template <int dim>
void DoFRenumbering::Cuthill_McKee (DoFHandler<dim>   &dof_handler,
				    const bool         reversed_numbering,
				    const bool         use_constraints,
				    const vector<int> &starting_indices) {
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
    
  const int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; -1 means
				   // that no new number was chosen yet
  vector<int> new_number(sparsity.n_rows(), -1);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<int> last_round_dofs (starting_indices);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]<0) || (last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = -1;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<int>(), -1));
  
				   // now if no valid points remain:
				   // find dof with lowest coordination
				   // number
  
  if (last_round_dofs.size() == 0)
    {
      int          starting_point   = -1;
      unsigned int min_coordination = n_dofs;
      for (int row=0; row<n_dofs; ++row) 
	{
	  unsigned int j;

					   // loop until we hit the end
					   // of this row's entries
	  for (j=sparsity.get_rowstart_indices()[row];
	       j<sparsity.get_rowstart_indices()[row+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
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
      if (starting_point == -1)
	starting_point = 0;
      
				       // initialize the first dof
      last_round_dofs.push_back (starting_point);
    };
  

				   // store next free dof index
  int         next_free_number = 0;

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
      vector<int> next_round_dofs;

				       // find all neighbors of the
				       // dofs numbered in the last
				       // round
      for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	     j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	  if (sparsity.get_column_numbers()[j] == -1)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
				       // sort dof numbers
      sort (next_round_dofs.begin(), next_round_dofs.end());

				       // delete multiple entries
      vector<int>::iterator end_sorted;
      end_sorted = unique (next_round_dofs.begin(), next_round_dofs.end());
      next_round_dofs.erase (end_sorted, next_round_dofs.end());

				       // eliminate dofs which are
				       // already numbered
      for (int s=next_round_dofs.size()-1; s>=0; --s)
	if (new_number[next_round_dofs[s]] != -1)
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
      for (vector<int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (unsigned int j=sparsity.get_rowstart_indices()[*s];
	       j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
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

//TODO: Allow incomplete renumbering for non-discretization values

#ifdef DEBUG
				   //  test for all indices numbered
  if (find (new_number.begin(), new_number.end(), -1) != new_number.end())
    Assert (false, ExcRenumberingIncomplete());
  Assert (next_free_number == n_dofs,
	  ExcRenumberingIncomplete());
#endif

  if (reversed_numbering)
    for (vector<int>::iterator i=new_number.begin(); i!=new_number.end(); ++i)
      *i = n_dofs-*i;

				   // actually perform renumbering;
				   // this is dimension specific and
				   // thus needs an own function
  dof_handler.renumber_dofs (new_number);
};




template <int dim>
void DoFRenumbering::Cuthill_McKee (MGDoFHandler<dim>      &dof_handler,
				    const unsigned int      level,
				    const bool              reversed_numbering,
				    const vector<int>      &starting_indices) {
				   // make the connection graph
  SparsityPattern sparsity (dof_handler.n_dofs(level),
			    dof_handler.max_couplings_between_dofs());
  MGDoFTools::make_sparsity_pattern (dof_handler, level, sparsity);
    
  const int n_dofs = sparsity.n_rows();
				   // store the new dof numbers; -1 means
				   // that no new number was chosen yet
  vector<int> new_number(n_dofs, -1);
  
				   // store the indices of the dofs renumbered
				   // in the last round. Default to starting
				   // points
  vector<int> last_round_dofs (starting_indices);
  
				   // delete disallowed elements
  for (unsigned int i=0; i<last_round_dofs.size(); ++i)
    if ((last_round_dofs[i]<0) || (last_round_dofs[i]>=n_dofs))
      last_round_dofs[i] = -1;
  
  remove_if (last_round_dofs.begin(), last_round_dofs.end(),
	     bind2nd(equal_to<int>(), -1));
  
				   // now if no valid points remain:
				   // find dof with lowest coordination
				   // number
  
  if (last_round_dofs.size() == 0)
    {
      int          starting_point   = -1;
      unsigned int min_coordination = n_dofs;
      for (int row=0; row<n_dofs; ++row) 
	{
	  unsigned int j;
	  for (j=sparsity.get_rowstart_indices()[row];
	       j<sparsity.get_rowstart_indices()[row+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
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
  int         next_free_number = 0;

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
      vector<int> next_round_dofs;

				       // find all neighbors of the
				       // dofs numbered in the last
				       // round
      for (unsigned int i=0; i<last_round_dofs.size(); ++i)
	for (unsigned int j=sparsity.get_rowstart_indices()[last_round_dofs[i]];
	     j<sparsity.get_rowstart_indices()[last_round_dofs[i]+1]; ++j)
	  if (sparsity.get_column_numbers()[j] == -1)
	    break;
	  else
	    next_round_dofs.push_back (sparsity.get_column_numbers()[j]);
      
				       // sort dof numbers
      sort (next_round_dofs.begin(), next_round_dofs.end());

				       // delete multiple entries
      vector<int>::iterator end_sorted;
      end_sorted = unique (next_round_dofs.begin(), next_round_dofs.end());
      next_round_dofs.erase (end_sorted, next_round_dofs.end());

				       // eliminate dofs which are
				       // already numbered
      for (int s=next_round_dofs.size()-1; s>=0; --s)
	if (new_number[next_round_dofs[s]] != -1)
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
      for (vector<int>::iterator s=next_round_dofs.begin();
	   s!=next_round_dofs.end(); ++s) 
	{
	  unsigned int coordination = 0;
	  for (unsigned int j=sparsity.get_rowstart_indices()[*s];
	       j<sparsity.get_rowstart_indices()[*s+1]; ++j)
	    if (sparsity.get_column_numbers()[j] == -1)
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
  if (find (new_number.begin(), new_number.end(), -1) != new_number.end())
    Assert (false, ExcRenumberingIncomplete());
  Assert (next_free_number == n_dofs,
	  ExcRenumberingIncomplete());
#endif

  if (reversed_numbering)
    for (vector<int>::iterator i=new_number.begin(); i!=new_number.end(); ++i)
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
  vector<int> local_dof_indices(dofs_per_cell);
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
  vector<vector<int> > component_to_dof_map (dof_handler.get_fe().n_components());
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
  
  int next_free_index = 0;
  vector<int> new_indices (dof_handler.n_dofs(), -1);
  for (unsigned int component=0; component<dof_handler.get_fe().n_components(); ++component)
    {
      const typename vector<int>::const_iterator
	begin_of_component = component_to_dof_map[component].begin(),
	end_of_component   = component_to_dof_map[component].end();
            
      for (typename vector<int>::const_iterator dof_index = begin_of_component;
	   dof_index != end_of_component; ++dof_index)
	new_indices[*dof_index] = next_free_index++;
    };

  Assert (next_free_index == static_cast<int>(dof_handler.n_dofs()),
	  ExcInternalError());

  dof_handler.renumber_dofs (new_indices);
};

  
  



// explicit instantiations
template
void DoFRenumbering::Cuthill_McKee (DoFHandler<deal_II_dimension> &dof_handler,
				    const bool                     reversed_numbering,
				    const bool                     use_constraints,
				    const vector<int>             &starting_indices);

template
void DoFRenumbering::Cuthill_McKee (MGDoFHandler<deal_II_dimension> &dof_handler,
				    const unsigned int               level,
				    const bool                       reversed_numbering,
				    const vector<int>               &starting_indices);

template
void DoFRenumbering::component_wise (DoFHandler<deal_II_dimension> &dof_handler,
				     const vector<unsigned int>    &component_order_arg);
