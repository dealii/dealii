/* $Id$ */

#include <grid/dof.h>
#include <grid/mg_dof.h>
#include <grid/dof_constraints.h>
#include <numerics/dof_renumbering.h>
#include <lac/sparsematrix.h>

#include <vector>
#include <map>
#include <algorithm>



template <int dim>
void DoFRenumbering::renumber_Cuthill_McKee (DoFHandler<dim>   &dof_handler,
					     const bool         reversed_numbering,
					     const bool         use_constraints,
					     const vector<int> &starting_indices) {
				   // make the connection graph
  SparseMatrixStruct sparsity (dof_handler.n_dofs(),
			       dof_handler.max_couplings_between_dofs());
  dof_handler.make_sparsity_pattern (sparsity);

  if (use_constraints) 
    {
      ConstraintMatrix constraints;
      dof_handler.make_hanging_node_constraints (constraints);
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
void DoFRenumbering::renumber_Cuthill_McKee (MGDoFHandler<dim>      &dof_handler,
					     const unsigned int      level,
					     const bool              reversed_numbering,
					     const vector<int>      &starting_indices) {
				   // make the connection graph
  SparseMatrixStruct sparsity (dof_handler.n_dofs(level),
			       dof_handler.max_couplings_between_dofs());
  dof_handler.make_sparsity_pattern (level, sparsity);
    
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






// explicit instantiations
template
void DoFRenumbering::renumber_Cuthill_McKee (DoFHandler<deal_II_dimension> &dof_handler,
					     const bool                     reversed_numbering,
					     const bool                     use_constraints,
					     const vector<int>             &starting_indices);

template
void DoFRenumbering::renumber_Cuthill_McKee (MGDoFHandler<deal_II_dimension> &dof_handler,
					     const unsigned int               level,
					     const bool                       reversed_numbering,
					     const vector<int>               &starting_indices);

