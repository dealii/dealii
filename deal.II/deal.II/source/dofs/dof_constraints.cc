//----------------------------  dof_constraints.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  dof_constraints.cc  ---------------------------


#include <dofs/dof_constraints.h>
#include <dofs/dof_constraints.templates.h>

#include <base/memory_consumption.h>
#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <set>



inline
bool
ConstraintMatrix::ConstraintLine::operator < (const ConstraintLine &a) const
{
  return line < a.line;
};



unsigned int
ConstraintMatrix::ConstraintLine::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (line) +
	  MemoryConsumption::memory_consumption (entries));
};



ConstraintMatrix::ConstraintMatrix () :
		lines(),
		sorted(false)
{};


void ConstraintMatrix::add_line (const unsigned int line)
{
  Assert (sorted==false, ExcMatrixIsClosed());

				   // check whether line already exists;
				   // it may, but then we need to quit
  for (unsigned int i=0; i!=lines.size(); ++i)
    if (lines[i].line == line)
      return;

				   // push a new line to the end of the
				   // list
  lines.push_back (ConstraintLine());
  lines.back().line = line;
};



void ConstraintMatrix::add_entry (const unsigned int line,
				  const unsigned int column,
				  const double       value)
{
  Assert (sorted==false, ExcMatrixIsClosed());

  std::vector<ConstraintLine>::iterator line_ptr;
  const std::vector<ConstraintLine>::const_iterator start=lines.begin();
				   // the usual case is that the line where
				   // a value is entered is the one we
				   // added last, so we search backward
  for (line_ptr=(lines.end()-1); line_ptr!=start; --line_ptr)
    if (line_ptr->line == line)
      break;

				   // if the loop didn't break, then
				   // line_ptr must be begin().
				   // we have an error if that doesn't
				   // point to 'line' then
  Assert (line_ptr->line==line, ExcLineInexistant(line));

				   // if in debug mode, check whether an
				   // entry for this column already
				   // exists and if its the same as
				   // the one entered at present
				   //
				   // in any case: exit the function if an
				   // entry for this column already exists,
				   // since we don't want to enter it twice
  for (std::vector<std::pair<unsigned int,double> >::const_iterator p=line_ptr->entries.begin();
       p != line_ptr->entries.end(); ++p)
    if (p->first == column)
      {
	Assert (p->second == value,
		ExcEntryAlreadyExists(line, column, p->second, value));
	return;
      };
  
  line_ptr->entries.push_back (std::make_pair(column,value));
};



void ConstraintMatrix::add_entries (const unsigned int                        line,
				    const std::vector<std::pair<unsigned int,double> > &col_val_pairs)
{
  Assert (sorted==false, ExcMatrixIsClosed());

  std::vector<ConstraintLine>::iterator line_ptr;
  const std::vector<ConstraintLine>::const_iterator start=lines.begin();
				   // the usual case is that the line where
				   // a value is entered is the one we
				   // added last, so we search backward
  for (line_ptr=(lines.end()-1); line_ptr!=start; --line_ptr)
    if (line_ptr->line == line)
      break;

				   // if the loop didn't break, then
				   // line_ptr must be begin().
				   // we have an error if that doesn't
				   // point to 'line' then
  Assert (line_ptr->line==line, ExcLineInexistant(line));

				   // if in debug mode, check whether an
				   // entry for this column already
				   // exists and if its the same as
				   // the one entered at present
				   //
				   // in any case: skip this entry if
				   // an entry for this column already
				   // exists, since we don't want to
				   // enter it twice
  for (std::vector<std::pair<unsigned int,double> >::const_iterator col_val_pair = col_val_pairs.begin();
       col_val_pair!=col_val_pairs.end(); ++col_val_pair)
    {
      for (std::vector<std::pair<unsigned int,double> >::const_iterator p=line_ptr->entries.begin();
	   p != line_ptr->entries.end(); ++p)
	if (p->first == col_val_pair->first)
	  {
					     // entry exists, break
					     // innermost loop
	    Assert (p->second == col_val_pair->second,
		    ExcEntryAlreadyExists(line, col_val_pair->first,
					  p->second, col_val_pair->second));
	    break;
	  };
      
      line_ptr->entries.push_back (*col_val_pair);
    };
};



void ConstraintMatrix::close ()
{
  Assert (sorted==false, ExcMatrixIsClosed());

				   // sort the entries in the different lines
				   // and strip zero entries
  std::vector<ConstraintLine>::iterator line = lines.begin(),
					endl = lines.end();
  for (; line!=endl; ++line)
    {
				       // first remove zero
				       // entries. that would mean
				       // that in the linear
				       // constraint for a node,
				       // x_i = ax_1 + bx_2 + ...,
				       // another node times 0
				       // appears. obviously,
				       // 0*something can be omitted
      line->entries.erase (remove_if (line->entries.begin(),
				      line->entries.end(),
				      std::compose1 (std::bind2nd (std::equal_to<double>(), 0),
						     std::select2nd<std::pair<unsigned int,double> >())),
                           line->entries.end());

				       // now sort the remainder
      sort (line->entries.begin(), line->entries.end());
    };
  
				   // sort the lines
  sort (lines.begin(), lines.end());

#ifdef DEBUG
				   // if in debug mode: check that no
				   // dof is constraint to another dof
				   // that is also constrained
  for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line)
    for (std::vector<std::pair<unsigned int,double> >::const_iterator entry=line->entries.begin();
	 entry!=line->entries.end(); ++entry)
      {
					 // make sure that
					 // entry->first is not the
					 // index of a line itself
	ConstraintLine test_line;
	test_line.line = entry->first;
	const std::vector<ConstraintLine>::const_iterator
	  test_line_position = lower_bound (lines.begin(),
					    lines.end(),
					    test_line);
	Assert ((test_line_position == lines.end())
		||
		(test_line_position->line != entry->first),
		ExcDoFConstrainedToConstrainedDoF(line->line, entry->first));
      };
#endif
  
  sorted = true;
};



void ConstraintMatrix::merge (const ConstraintMatrix &other_constraints)
{
				   // first check whether the
				   // constraints in the two objects
				   // are for different degrees of
				   // freedom
  if (true)
    {
				       // first insert all dofs in
				       // this object into a list...
      std::set<unsigned int> this_dofs;
      for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
	   line!=lines.end(); ++line)
	this_dofs.insert (line->line);

				       // ...then check whether it
				       // appears in the other object
				       // as well. note that we have
				       // to do this in a somewhat
				       // complicated style since the
				       // two objects may not be
				       // sorted
      for (std::vector<ConstraintLine>::const_iterator
	     line=other_constraints.lines.begin();
	   line!=other_constraints.lines.end(); ++line)
	AssertThrow (this_dofs.find (line->line) == this_dofs.end(),
		     ExcDoFIsConstrainedFromBothObjects (line->line));
    };

				   // store the previous state with
				   // respect to sorting
  const bool object_was_sorted = sorted;
  sorted = false;

				   // append new lines at the end
  lines.insert (lines.end(),
		other_constraints.lines.begin(),
		other_constraints.lines.end());

				   // if the object was sorted before,
				   // then make sure it is so
				   // afterwards as well. otherwise
				   // leave everything in the unsorted
				   // state
  if (object_was_sorted == true)
    close ();
};



void ConstraintMatrix::clear ()
{
  std::vector<ConstraintLine> tmp;
  lines.swap (tmp);
  sorted = false;
};



void ConstraintMatrix::condense (const SparsityPattern &uncondensed,
				 SparsityPattern       &condensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed.is_compressed() == true, ExcMatrixNotClosed());
  Assert (uncondensed.n_rows() == uncondensed.n_cols(),
	  ExcMatrixNotSquare());


				   // store for each line of the matrix
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  std::vector<int> new_line;

  new_line.reserve (uncondensed.n_rows());

  std::vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                                shift           = 0;
  unsigned int n_rows = uncondensed.n_rows();

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
  for (unsigned int row=0; row<uncondensed.n_rows(); ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entries if column will not
				       // be condensed away, distribute
				       // otherwise
      for (unsigned int j=uncondensed.get_rowstart_indices()[row];
	   j<uncondensed.get_rowstart_indices()[row+1]; ++j)
	if (new_line[uncondensed.get_column_numbers()[j]] != -1)
	  condensed.add (new_line[row], new_line[uncondensed.get_column_numbers()[j]]);
	else 
	  {
					     // let c point to the constraint
					     // of this column
	    std::vector<ConstraintLine>::const_iterator c = lines.begin();
	    while (c->line != uncondensed.get_column_numbers()[j])
	      ++c;

	    for (unsigned int q=0; q!=c->entries.size(); ++q) 
	      condensed.add (new_line[row], new_line[c->entries[q].first]);
	  }
    else
				       // line must be distributed
      {
	for (unsigned int j=uncondensed.get_rowstart_indices()[row];
	     j<uncondensed.get_rowstart_indices()[row+1]; ++j)
					   // for each entry: distribute
	  if (new_line[uncondensed.get_column_numbers()[j]] != -1)
					     // column is not constrained
	    for (unsigned int q=0; q!=next_constraint->entries.size(); ++q) 
	      condensed.add (new_line[next_constraint->entries[q].first],
			     new_line[uncondensed.get_column_numbers()[j]]);
	
	  else
					     // not only this line but
					     // also this col is constrained
	    {
					       // let c point to the constraint
					       // of this column
	      std::vector<ConstraintLine>::const_iterator c = lines.begin();
	      while (c->line != uncondensed.get_column_numbers()[j]) ++c;
	      
	      for (unsigned int p=0; p!=c->entries.size(); ++p)
		for (unsigned int q=0; q!=next_constraint->entries.size(); ++q)
		  condensed.add (new_line[next_constraint->entries[q].first],
				 new_line[c->entries[p].first]);
	    };
	
	++next_constraint;
      };

  condensed.compress();
};



void ConstraintMatrix::condense (SparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
  std::vector<int> distribute(sparsity.n_rows(), -1);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = static_cast<signed int>(c);

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == -1)
					 // regular line. loop over cols
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
	  {
	    const unsigned int column = sparsity.get_column_numbers()[j];
	    
					     // end of row reached?
	    if (column == SparsityPattern::invalid_entry)
	      break;
	    else
	      if (distribute[column] != -1)
		{
						   // distribute entry at regular
						   // row @p{row} and irregular column
						   // sparsity.colnums[j]
		  for (unsigned int q=0;
		       q!=lines[distribute[column]].entries.size();
		       ++q) 
		    sparsity.add (row,
				  lines[distribute[column]].entries[q].first);
		};
	  }
      else
					 // row must be distributed
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
					   // end of row reached?
	  if (sparsity.get_column_numbers()[j] == SparsityPattern::invalid_entry)
	    break;
	  else
	    {
	      if (distribute[sparsity.get_column_numbers()[j]] == -1)
						 // distribute entry at irregular
						 // row @p{row} and regular column
						 // sparsity.colnums[j]
		for (unsigned int q=0;
		     q!=lines[distribute[row]].entries.size(); ++q) 
		  sparsity.add (lines[distribute[row]].entries[q].first,
				sparsity.get_column_numbers()[j]);
	      else
						 // distribute entry at irregular
						 // row @p{row} and irregular column
						 // sparsity.get_column_numbers()[j]
		for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		  for (unsigned int q=0;
		       q!=lines[distribute[sparsity.get_column_numbers()[j]]]
				      .entries.size(); ++q)
		    sparsity.add (lines[distribute[row]].entries[p].first,
				  lines[distribute[sparsity.get_column_numbers()[j]]]
				  .entries[q].first);
	    };
    };
  
  sparsity.compress();
};



void ConstraintMatrix::condense (BlockSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcMatrixNotSquare());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcMatrixNotSquare());
  
  const BlockIndices &
    index_mapping = sparsity.get_column_indices();

  const unsigned int n_blocks = sparsity.n_block_rows();
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
  std::vector<int> distribute (sparsity.n_rows(), -1);
  
  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = static_cast<signed int>(c);

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // get index of this row
				       // within the blocks
      const std::pair<unsigned int,unsigned int>
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
	  for (unsigned int block_col=0; block_col<n_blocks; ++block_col)
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
						     // nothing more
						     // to do
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
						       // global_col
		      {
			for (unsigned int q=0;
			     q!=lines[distribute[global_col]]
					    .entries.size(); ++q)
			  sparsity.add (row,
					lines[distribute[global_col]].entries[q].first);
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
	  for (unsigned int block_col=0; block_col<n_blocks; ++block_col)
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
						     // nothing more to do
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
						       // global_col.
		      {
			for (unsigned int q=0; q!=lines[distribute[row]].entries.size(); ++q) 
			  sparsity.add (lines[distribute[row]].entries[q].first,
					global_col);
		      }
		    else
						       // distribute entry at irregular
						       // row @p{row} and irregular column
						       // @p{global_col}
		      {
			for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
			  for (unsigned int q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
			    sparsity.add (lines[distribute[row]].entries[p].first,
					  lines[distribute[global_col]].entries[q].first);
		      };
		  };
	    };
	};
    };
  
  sparsity.compress();
};



unsigned int ConstraintMatrix::n_constraints () const
{
  return lines.size();
};



bool ConstraintMatrix::is_constrained (const unsigned int index) const 
{
  if (sorted == true)
    {
      
      ConstraintLine index_comparison;
      index_comparison.line = index;
      
      return binary_search (lines.begin (),
			    lines.end (),
			    index_comparison);
    }
  else
    {
      for (std::vector<ConstraintLine>::const_iterator i=lines.begin();
	   i!=lines.end(); ++i)
	if (i->line == index)
	  return true;

      return false;
    };
};



unsigned int ConstraintMatrix::max_constraint_indirections () const 
{
  unsigned int return_value = 0;
  for (std::vector<ConstraintLine>::const_iterator i=lines.begin();
       i!=lines.end(); ++i)
				     // use static cast, since
				     // typeof(size)==size_t, which is
				     // != unsigned int on AIX
    return_value = std::max(return_value,
			    static_cast<unsigned int>(i->entries.size()));

  return return_value;
};

    

void ConstraintMatrix::print (std::ostream &out) const
{
  for (unsigned int i=0; i!=lines.size(); ++i)
    for (unsigned int j=0; j!=lines[i].entries.size(); ++j)
      out << "    " << lines[i].line
	  << " " << lines[i].entries[j].first
	  << ":  " << lines[i].entries[j].second << std::endl;

  AssertThrow (out, ExcIO());
};



unsigned int
ConstraintMatrix::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (lines) +
	  MemoryConsumption::memory_consumption (sorted));
};





// explicit instantiations
//
// define a list of functions for vectors and matrices, respectively,
// where the vector/matrix can be replaced using a preprocessor
// variable VectorType/MatrixType. note that we cannot do so by using
// a preprocessor function with one arg, since
// #vector_functions(BlockVector<double>)# is not recognized as one
// arg, and putting parentheses around the arg yields incorrect
// syntax...

#define vector_functions \
  template void ConstraintMatrix::condense(const VectorType &uncondensed,\
					   VectorType       &condensed) const;\
  template void ConstraintMatrix::condense(VectorType &vec) const;\
  template void ConstraintMatrix::set_zero(VectorType &vec) const;\
  template void ConstraintMatrix::distribute(const VectorType &condensed,\
					     VectorType       &uncondensed) const;\
  template void ConstraintMatrix::distribute(VectorType &vec) const;



#define matrix_functions_1 \
  template void ConstraintMatrix::condense(const MatrixType &uncondensed,\
					   MatrixType       &condensed) const;
#define matrix_functions_2 \
  template void ConstraintMatrix::condense(MatrixType &uncondensed) const;



#define VectorType Vector<float>
vector_functions;
#undef VectorType

#define VectorType Vector<double>
vector_functions;
#undef VectorType

#define VectorType BlockVector<double>
vector_functions;
#undef VectorType



#define MatrixType SparseMatrix<float>
matrix_functions_1;
matrix_functions_2;
#undef MatrixType

#define MatrixType SparseMatrix<double>
matrix_functions_1;
matrix_functions_2;
#undef MatrixType

// block sparse matrices are only implemented for one of the two matrix functions
#define MatrixType BlockSparseMatrix<double>
matrix_functions_2;
#undef MatrixType


