//----------------------------  dof_constraints.templates.h  ---------------------------
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
//----------------------------  dof_constraints.templates.h  ---------------------------
#ifndef __deal2__dof_constraints_templates_h
#define __deal2__dof_constraints_templates_h


#include <dofs/dof_constraints.h>




template<typename number>
void
ConstraintMatrix::condense (const SparseMatrix<number> &uncondensed,
			    SparseMatrix<number>       &condensed) const
{
  const SparsityPattern &uncondensed_struct = uncondensed.get_sparsity_pattern ();
  
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.is_compressed() == true, ExcMatrixNotClosed());
  Assert (condensed.get_sparsity_pattern().is_compressed() == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.n_rows() == uncondensed_struct.n_cols(),
	  ExcMatrixNotSquare());
  Assert (condensed.n() == condensed.m(),
	  ExcMatrixNotSquare());
  Assert (condensed.n()+n_constraints() == uncondensed.n(),
	  ExcWrongDimension());

				   // store for each line of the matrix
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  vector<int> new_line;

  new_line.reserve (uncondensed_struct.n_rows());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  const unsigned int n_rows = uncondensed_struct.n_rows();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      new_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  new_line.push_back (-1);
					   // note that @p{lines} is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		new_line.push_back (i-shift);
	      break;
	    };
	}
      else
	new_line.push_back (row-shift);


  next_constraint = lines.begin();
				   // note: in this loop we need not check
				   // whether @p{next_constraint} is a valid
				   // iterator, since @p{next_constraint} is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int row=0; row<uncondensed_struct.n_rows(); ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entries if column will not
				       // be condensed away, distribute
				       // otherwise
      for (unsigned int j=uncondensed_struct.get_rowstart_indices()[row];
	   j<uncondensed_struct.get_rowstart_indices()[row+1]; ++j)
	if (new_line[uncondensed_struct.get_column_numbers()[j]] != -1)
	  condensed.add (new_line[row], new_line[uncondensed_struct.get_column_numbers()[j]],
			 uncondensed.global_entry(j));
	else 
	  {
					     // let c point to the constraint
					     // of this column
	    vector<ConstraintLine>::const_iterator c = lines.begin();
	    while (c->line != uncondensed_struct.get_column_numbers()[j])
	      ++c;

	    for (unsigned int q=0; q!=c->entries.size(); ++q)
					       // distribute to rows with
					       // appropriate weight
	      condensed.add (new_line[row], new_line[c->entries[q].first],
			     uncondensed.global_entry(j) * c->entries[q].second);
	  }
    else
				       // line must be distributed
      {
	for (unsigned int j=uncondensed_struct.get_rowstart_indices()[row];
	     j<uncondensed_struct.get_rowstart_indices()[row+1]; ++j)
					   // for each column: distribute
	  if (new_line[uncondensed_struct.get_column_numbers()[j]] != -1)
					     // column is not constrained
	    for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
		condensed.add (new_line[next_constraint->entries[q].first],
			       new_line[uncondensed_struct.get_column_numbers()[j]],
			       uncondensed.global_entry(j) *
			       next_constraint->entries[q].second);
	
	  else
					     // not only this line but
					     // also this col is constrained
	    {
					       // let c point to the constraint
					       // of this column
	      vector<ConstraintLine>::const_iterator c = lines.begin();
	      while (c->line != uncondensed_struct.get_column_numbers()[j])
		++c;
	      
	      for (unsigned int p=0; p!=c->entries.size(); ++p)
		for (unsigned int q=0; q!=next_constraint->entries.size(); ++q)
		    condensed.add (new_line[next_constraint->entries[q].first],
				   new_line[c->entries[p].first],
				   uncondensed.global_entry(j) *
				   next_constraint->entries[q].second *
				   c->entries[p].second);
	    };
	
	++next_constraint;
      };
};



template<typename number>
void
ConstraintMatrix::condense (SparseMatrix<number> &uncondensed) const
{
  const SparsityPattern &sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
  vector<int> distribute (sparsity.n_rows(), -1);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = static_cast<signed int>(c);

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == -1)
					 // regular line. loop over cols
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
					   // end of row reached?
	  if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	    {
					       // this should not happen, since
					       // we only operate on compressed
					       // matrices!
	      Assert (false, ExcMatrixNotClosed());
	      break;
	    }
	  else
	    {
	      if (distribute[sparsity.get_column_numbers()[j]] != -1)
						 // distribute entry at regular
						 // row @p{row} and irregular column
						 // sparsity.get_column_numbers()[j]; set old
						 // entry to zero
		{
		  for (unsigned int q=0;
		       q!=lines[distribute[sparsity.get_column_numbers()[j]]]
				      .entries.size(); ++q)
		    uncondensed.add (row,
				     lines[distribute[sparsity.get_column_numbers()[j]]]
				     .entries[q].first,
				     uncondensed.global_entry(j) *
				     lines[distribute[sparsity.get_column_numbers()[j]]]
				     .entries[q].second);
		
		  uncondensed.global_entry(j) = 0.;
		};
	    }
      else
					 // row must be distributed
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
					   // end of row reached?
	  if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	    {
					       // this should not happen, since
					       // we only operate on compressed
					       // matrices!
	      Assert (false, ExcMatrixNotClosed());
	      break;
	    }
	  else
	    {
	      if (distribute[sparsity.get_column_numbers()[j]] == -1)
						 // distribute entry at irregular
						 // row @p{row} and regular column
						 // sparsity.get_column_numbers()[j]. set old
						 // entry to zero
		{
		  for (unsigned int q=0;
		       q!=lines[distribute[row]].entries.size(); ++q) 
		    uncondensed.add (lines[distribute[row]].entries[q].first,
				     sparsity.get_column_numbers()[j],
				     uncondensed.global_entry(j) *
				     lines[distribute[row]].entries[q].second);
		
		  uncondensed.global_entry(j) = 0.;
		}
	      else
						 // distribute entry at irregular
						 // row @p{row} and irregular column
						 // sparsity.get_column_numbers()[j]
						 // set old entry to one if on main
						 // diagonal, zero otherwise
		{
		  for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		    for (unsigned int q=0;
			 q!=lines[distribute[sparsity.get_column_numbers()[j]]]
					.entries.size(); ++q)
		      uncondensed.add (lines[distribute[row]].entries[p].first,
				       lines[distribute[sparsity.get_column_numbers()[j]]]
				       .entries[q].first,
				       uncondensed.global_entry(j) *
				       lines[distribute[row]].entries[p].second *
				       lines[distribute[sparsity.get_column_numbers()[j]]]
				       .entries[q].second);
		
		  uncondensed.global_entry(j) = (row == sparsity.get_column_numbers()[j] ?
						 1. : 0. );
		};
	    };
    };
};



template <typename number, int blocks>
void
ConstraintMatrix::condense (BlockSparseMatrix<number,blocks,blocks> &uncondensed) const
{
  const BlockSparsityPattern &
    sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcMatrixNotSquare());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcMatrixNotSquare());

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
  vector<int> distribute (sparsity.n_rows(), -1);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = static_cast<signed int>(c);

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // get index of this row
				       // within the blocks
      const pair<unsigned int,unsigned int>
	block_index = index_mapping.global_to_local(row);
      const unsigned int block_row = block_index.first;
      
      if (distribute[row] == -1)
					 // regular line. loop over
					 // all columns and see
					 // whether this column must
					 // be distributed
	{

					   // to loop over all entries
					   // in this row, we have to
					   // loop over all blocks in
					   // this blockrow and the
					   // corresponding row
					   // therein
	  for (unsigned int block_col=0; block_col<blocks; ++block_col)
	    {
	      const SparsityPattern &
		block_sparsity = sparsity.block(block_row, block_col);
	      
	      const unsigned int
		first = block_sparsity.get_rowstart_indices()[block_index.second],
		last  = block_sparsity.get_rowstart_indices()[block_index.second+1];
	      for (unsigned int j=first; j<last; ++j)
						 // end of row reached?
		if (block_sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
		  {
						     // this should not happen, since
						     // we only operate on compressed
						     // matrices!
		    Assert (false, ExcMatrixNotClosed());
		    break;
		  }
		else
		  {
		    const unsigned int global_col
		      = index_mapping.local_to_global(block_col,
						      block_sparsity.get_column_numbers()[j]);
		    
		    if (distribute[global_col] != -1)
						       // distribute entry at regular
						       // row @p{row} and irregular column
						       // global_col;
						       // set old entry to zero
		      {
			const double old_value =
			  uncondensed.block(block_row,block_col).global_entry(j);
			
			for (unsigned int q=0;
			     q!=lines[distribute[global_col]]
					    .entries.size(); ++q)
			  uncondensed.add (row,
					   lines[distribute[global_col]].entries[q].first,
					   old_value *
					   lines[distribute[global_col]].entries[q].second);
			
			uncondensed.block(block_row,block_col).global_entry(j)
			  = 0.;
		      };
		  };
	    };
	}
      else
	{
					   // row must be
					   // distributed. split the
					   // whole row into the
					   // chunks defined by the
					   // blocks
	  for (unsigned int block_col=0; block_col<blocks; ++block_col)
	    {
	      const SparsityPattern &
		block_sparsity = sparsity.block(block_row,block_col);
	      
	      const unsigned int
		first = block_sparsity.get_rowstart_indices()[block_index.second],
		last  = block_sparsity.get_rowstart_indices()[block_index.second+1];
      
	      for (unsigned int j=first; j<last; ++j)
						 // end of row reached?
		if (block_sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
		  {
						     // this should not happen, since
						     // we only operate on compressed
						     // matrices!
		    Assert (false, ExcMatrixNotClosed());
		    break;
		  }
		else
		  {
		    const unsigned int global_col
		      = index_mapping.local_to_global (block_col,
						       block_sparsity.get_column_numbers()[j]);
		    
		    if (distribute[global_col] == -1)
						       // distribute entry at irregular
						       // row @p{row} and regular column
						       // global_col. set old
						       // entry to zero
		      {
			const double old_value
			  = uncondensed.block(block_row,block_col).global_entry(j);
			  
			for (unsigned int q=0;
			     q!=lines[distribute[row]].entries.size(); ++q) 
			  uncondensed.add (lines[distribute[row]].entries[q].first,
					   global_col,
					   old_value *
					   lines[distribute[row]].entries[q].second);
			
			uncondensed.block(block_row,block_col).global_entry(j) = 0.;
		      }
		    else
						       // distribute entry at irregular
						       // row @p{row} and irregular column
						       // @p{global_col}
						       // set old entry to one if on main
						       // diagonal, zero otherwise
		      {
			const double old_value
			  = uncondensed.block(block_row,block_col).global_entry(j);
			  
			for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
			  for (unsigned int q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
			    uncondensed.add (lines[distribute[row]].entries[p].first,
					     lines[distribute[global_col]].entries[q].first,
					     old_value *
					     lines[distribute[row]].entries[p].second *
					     lines[distribute[global_col]].entries[q].second);
		
			uncondensed.block(block_row,block_col).global_entry(j)
			  = (row == global_col ? 1. : 0. );
		      };
		  };
	    };
	};
    };
};



template<class VectorType>
void
ConstraintMatrix::condense (const VectorType &uncondensed,
			    VectorType       &condensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (condensed.size()+n_constraints() == uncondensed.size(),
	  ExcWrongDimension());
  
				   // store for each line of the vector
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  vector<int> new_line;

  new_line.reserve (uncondensed.size());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  unsigned int n_rows = uncondensed.size();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      new_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  new_line.push_back (-1);
					   // note that @p{lines} is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		new_line.push_back (i-shift);
	      break;
	    };
	}
      else
	new_line.push_back (row-shift);


  next_constraint = lines.begin();
				   // note: in this loop we need not check
				   // whether @p{next_constraint} is a valid
				   // iterator, since @p{next_constraint} is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int row=0; row<uncondensed.size(); ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entry
      condensed(new_line[row]) += uncondensed(row);

    else
				       // line must be distributed
      {
	for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	  condensed(new_line[next_constraint->entries[q].first])
	    +=
	    uncondensed(row) * next_constraint->entries[q].second;

	++next_constraint;
      };
};



template <class VectorType>
void
ConstraintMatrix::condense (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  if (lines.size() == 0)
				     // nothing to do
    return;
  
  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (unsigned int row=0; row<vec.size(); ++row)
    if (row == next_constraint->line)
				       // line must be distributed
      {
	for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	  vec(next_constraint->entries[q].first)
	    +=
	    vec(row) * next_constraint->entries[q].second;
					 // set entry to zero
	vec(row) = 0.;
	
	++next_constraint;
	if (next_constraint == lines.end())
					   // nothing more to do
	  break;
      };
};



template <class VectorType>
void
ConstraintMatrix::set_zero (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  if (lines.size() == 0)
				     // nothing to do
    return;
  
  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (unsigned int row=0; row<vec.size(); ++row)
    if (row == next_constraint->line)
      {
					 // set entry to zero
	vec(row) = 0.;
	
	++next_constraint;
	if (next_constraint == lines.end())
					   // nothing more to do
	  break;
      };
};



template<class VectorType>
void
ConstraintMatrix::distribute (const VectorType &condensed,
			      VectorType       &uncondensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (condensed.size()+n_constraints() == uncondensed.size(),
	  ExcWrongDimension());

				   // store for each line of the new vector
				   // its old line number before
				   // distribution. If the shift is
				   // -1, this line was condensed away
  vector<int> old_line;

  old_line.reserve (uncondensed.size());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  unsigned int n_rows = uncondensed.size();

  if (next_constraint == lines.end()) 
				     // if no constraint is to be handled
    for (unsigned int row=0; row!=n_rows; ++row)
      old_line.push_back (row);
  else
    for (unsigned int row=0; row!=n_rows; ++row)
      if (row == next_constraint->line)
	{
					   // this line is constrained
	  old_line.push_back (-1);
					   // note that @p{lines} is ordered
	  ++shift;
	  ++next_constraint;
	  if (next_constraint == lines.end())
					     // nothing more to do; finish rest
					     // of loop
	    {
	      for (unsigned int i=row+1; i<n_rows; ++i)
		old_line.push_back (i-shift);
	      break;
	    };
	}
      else
	old_line.push_back (row-shift);


  next_constraint = lines.begin();
				   // note: in this loop we need not check
				   // whether @p{next_constraint} is a valid
				   // iterator, since @p{next_constraint} is
				   // only evaluated so often as there are
				   // entries in new_line[*] which tells us
				   // which constraints exist
  for (unsigned int line=0; line<uncondensed.size(); ++line) 
    if (old_line[line] != -1)
				       // line was not condensed away
      uncondensed(line) = condensed(old_line[line]);
    else
      {
					 // line was condensed away, create it newly
					 // first set it to zero
	uncondensed(line) = 0.;
					 // then add the different contributions
	for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
	  uncondensed(line) += (condensed(old_line[next_constraint->entries[i].first]) *
				next_constraint->entries[i].second);
	++next_constraint;
      };
};



template<class VectorType>
void
ConstraintMatrix::distribute (VectorType &vec) const
{
  Assert (sorted == true, ExcMatrixNotClosed());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (; next_constraint != lines.end(); ++next_constraint) 
    {
				       // make entry in line next_constraint.line
				       // first set it to zero
      vec(next_constraint->line) = 0.;
				       // then add the different contributions
      for (unsigned int i=0; i<next_constraint->entries.size(); ++i)
	vec(next_constraint->line) += (vec(next_constraint->entries[i].first) *
				       next_constraint->entries[i].second);
    };
};


#endif
