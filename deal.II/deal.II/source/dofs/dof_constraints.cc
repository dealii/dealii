//----------------------------  dof_constraints.cc  ---------------------------
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
//----------------------------  dof_constraints.cc  ---------------------------


#include <dofs/dof_constraints.h>
#include <dofs/dof_constraints.templates.h>

#include <lac/sparsity_pattern.h>
#include <lac/vector.h>
#include <lac/block_vector.h>
#include <lac/sparse_matrix.h>
#include <lac/block_sparse_matrix.h>
#include <iostream>
#include <algorithm>
#include <numeric>


inline
bool
ConstraintMatrix::ConstraintLine::operator < (const ConstraintLine &a) const
{
  return line < a.line;
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

  vector<ConstraintLine>::iterator line_ptr;
  const vector<ConstraintLine>::const_iterator start=lines.begin();
				   // the usual case is that the line where
				   // a value is entered is the one we
				   // added last, so we search backward
  for (line_ptr=&lines.back(); line_ptr!=start; --line_ptr)
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
  for (vector<pair<unsigned int,double> >::const_iterator p=line_ptr->entries.begin();
       p != line_ptr->entries.end(); ++p)
    if (p->first == column)
      {
	Assert (p->second == value,
		ExcEntryAlreadyExists(line, column, p->second, value));
	return;
      };
  
  line_ptr->entries.push_back (make_pair(column,value));
};



void ConstraintMatrix::close ()
{
  Assert (sorted==false, ExcMatrixIsClosed());

				   // sort the entries in the different lines
				   // and strip zero entries
  vector<ConstraintLine>::iterator line = lines.begin(),
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
				      compose1 (bind2nd (equal_to<double>(), 0),
						select2nd<pair<unsigned int,double> >())),
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
  for (vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line)
    for (vector<pair<unsigned int,double> >::const_iterator entry=line->entries.begin();
	 entry!=line->entries.end(); ++entry)
      {
					 // make sure that
					 // entry->first is not the
					 // index of a line itself
	const ConstraintLine test_line = { entry->first,
					   vector<pair<unsigned int,double> >() };
	const vector<ConstraintLine>::const_iterator
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



void ConstraintMatrix::clear ()
{
  lines = vector<ConstraintLine>();
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
  vector<int> new_line;

  new_line.reserve (uncondensed.n_rows());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
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
					   // note that #lines# is ordered	  
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
				   // whether #next_constraint# is a valid
				   // iterator, since #next_constraint# is
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
	    vector<ConstraintLine>::const_iterator c = lines.begin();
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
	      vector<ConstraintLine>::const_iterator c = lines.begin();
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
  vector<int> distribute(sparsity.n_rows(), -1);
  
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
						   // row #row# and irregular column
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
						 // row #row# and regular column
						 // sparsity.colnums[j]
		for (unsigned int q=0;
		     q!=lines[distribute[row]].entries.size(); ++q) 
		  sparsity.add (lines[distribute[row]].entries[q].first,
				sparsity.get_column_numbers()[j]);
	      else
						 // distribute entry at irregular
						 // row #row# and irregular column
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
      for (vector<ConstraintLine>::const_iterator i=lines.begin();
	   i!=lines.end(); ++i)
	if (i->line == index)
	  return true;

      return false;
    };
};



unsigned int ConstraintMatrix::max_constraint_indirections () const 
{
  unsigned int return_value = 0;
  for (vector<ConstraintLine>::const_iterator i=lines.begin();
       i!=lines.end(); ++i)
				     // use static cast, since
				     // typeof(size)==size_t, which is
				     // != unsigned int on AIX
    return_value = max(return_value,
		       static_cast<unsigned int>(i->entries.size()));

  return return_value;
};

    

void ConstraintMatrix::print (ostream &out) const {
  for (unsigned int i=0; i!=lines.size(); ++i)
    for (unsigned int j=0; j!=lines[i].entries.size(); ++j)
      out << "    " << lines[i].line
	  << " " << lines[i].entries[j].first
	  << ":  " << lines[i].entries[j].second << endl;

  AssertThrow (out, ExcIO());
};





// explicit instantiations
//
// define a list of functions for vectors and matrices, respectively,
// where the vector/matrix can be replaced using a preprocessor
// variable VectorType/MatrixType. note that we cannot do so by using
// a preprocessor function with one arg, since
// #vector_functions(BlockVector<2,double>)# is not recognized as one
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

#define VectorType BlockVector<2,double>
vector_functions;
#undef VectorType

#define VectorType BlockVector<3,double>
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
#define MatrixType BlockSparseMatrix<double,2,2>
matrix_functions_2;
#undef MatrixType

#define MatrixType BlockSparseMatrix<double,3,3>
matrix_functions_2;
#undef MatrixType


#define MatrixType BlockSparsityPattern<2,2>
matrix_functions_2;
#undef MatrixType

#define MatrixType BlockSparsityPattern<3,3>
matrix_functions_2;
#undef MatrixType

