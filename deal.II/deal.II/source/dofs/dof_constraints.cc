/* $Id$ */
/* Copyright W. Bangerth, University of Heidelberg, 1998 */



#include <grid/dof_constraints.h>
#include <lac/sparsematrix.h>
#include <lac/vector.h>
#include <iostream>
#include <algorithm>
#include <numeric>



bool ConstraintMatrix::ConstraintLine::operator < (const ConstraintMatrix::ConstraintLine &a) const {
  return line < a.line;
};




ConstraintMatrix::ConstraintMatrix () :
		lines(), sorted(false) {};



void ConstraintMatrix::add_line (const unsigned int line) {
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
				  const double       value) {
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




void ConstraintMatrix::close () {
  Assert (sorted==false, ExcMatrixIsClosed());

				   // sort the entries in the different lines
				   // and strip zero entries
  vector<ConstraintLine>::iterator line = lines.begin(),
				   endl = lines.end();
  for (; line!=endl; ++line)
    {
				       // first remove zero entries
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
  
  sorted = true;
};

	  

void ConstraintMatrix::clear () {
  lines = vector<ConstraintLine>();
  sorted = false;
};



void ConstraintMatrix::condense (const SparseMatrixStruct &uncondensed,
				 SparseMatrixStruct       &condensed) const {
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
	    while (c->line != (unsigned int)uncondensed.get_column_numbers()[j])
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
	      while (c->line != (unsigned int)uncondensed.get_column_numbers()[j]) ++c;
	      
	      for (unsigned int p=0; p!=c->entries.size(); ++p)
		for (unsigned int q=0; q!=next_constraint->entries.size(); ++q)
		    condensed.add (new_line[next_constraint->entries[q].first],
				   new_line[c->entries[p].first]);
	    };
	
	++next_constraint;
      };

  condensed.compress();
};



void ConstraintMatrix::condense (SparseMatrixStruct &sparsity) const {
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

  int n_rows = sparsity.n_rows();
  for (int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == -1)
					 // regular line. loop over cols
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
	  {
	    const int column = sparsity.get_column_numbers()[j];
	    
					     // end of row reached?
	    if (column == -1)
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
	  if (sparsity.get_column_numbers()[j] == -1)
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



void ConstraintMatrix::condense (const SparseMatrix<double> &uncondensed,
				 SparseMatrix<double>       &condensed) const {
  const SparseMatrixStruct &uncondensed_struct = uncondensed.get_sparsity_pattern ();
  
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
  unsigned int n_rows = uncondensed_struct.n_rows();

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
	    while (c->line != (unsigned int)uncondensed_struct.get_column_numbers()[j])
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
	      while (c->line != (unsigned int)uncondensed_struct.get_column_numbers()[j])
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



void ConstraintMatrix::condense (SparseMatrix<double> &uncondensed) const {
  const SparseMatrixStruct &sparsity = uncondensed.get_sparsity_pattern ();

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

  int n_rows = sparsity.n_rows();
  for (int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == -1)
					 // regular line. loop over cols
	for (unsigned int j=sparsity.get_rowstart_indices()[row];
	     j<sparsity.get_rowstart_indices()[row+1]; ++j)
					   // end of row reached?
	  if (sparsity.get_column_numbers()[j] == -1)
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
						 // row #row# and irregular column
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
	  if (sparsity.get_column_numbers()[j] == -1)
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
						 // row #row# and regular column
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
						 // row #row# and irregular column
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






void ConstraintMatrix::condense (const Vector<double> &uncondensed,
				 Vector<double>       &condensed) const {
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



void ConstraintMatrix::condense (Vector<double> &vec) const {
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
  



void ConstraintMatrix::distribute (const Vector<double> &condensed,
				   Vector<double>       &uncondensed) const {
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
					   // note that #lines# is ordered
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
				   // whether #next_constraint# is a valid
				   // iterator, since #next_constraint# is
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



void ConstraintMatrix::distribute (Vector<double> &vec) const {
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




unsigned int ConstraintMatrix::n_constraints () const {
  return lines.size();
};



void ConstraintMatrix::print (ostream &out) const {
  for (unsigned int i=0; i!=lines.size(); ++i)
    for (unsigned int j=0; j!=lines[i].entries.size(); ++j)
      out << "    " << lines[i].line
	  << " " << lines[i].entries[j].first
	  << ":  " << lines[i].entries[j].second << endl;

  AssertThrow (out, ExcIO());
};
