/* $Id$ */


#include <algo.h>
#include <grid/dof_constraints.h>
#include <lac/dsmatrix.h>
#include <lac/dvector.h>
#include <iostream.h>




bool ConstraintMatrix::ConstraintLine::operator < (const ConstraintMatrix::ConstraintLine &a) const {
  return line < a.line;
};




ConstraintMatrix::ConstraintMatrix () :
		lines(), sorted(false) {};



void ConstraintMatrix::add_line (const unsigned int line) {
				   // check whether line already exists
  Assert (sorted==false, ExcMatrixIsClosed());
#ifdef DEBUG
  for (unsigned int i=0; i!=lines.size(); ++i)
    Assert (lines[i].line != line, ExcLineExists(line));
#endif

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
  if (lines.back().line == line)
    line_ptr = &lines.back();
  else
    {
      for (line_ptr = &lines.back()-1; line_ptr>=lines.begin(); --line_ptr)
	if ((*line_ptr).line == line)
	  break;
      Assert (false, ExcLineInexistant(line));
    };

  (*line_ptr).entries.push_back (make_pair(column,value));
};




void ConstraintMatrix::close () {
  Assert (sorted==false, ExcMatrixIsClosed());

				   // sort the entries in the different lines
  vector<ConstraintLine>::iterator line = lines.begin(),
				   endl = lines.end();
  for (; line!=endl; ++line)
    sort ((*line).entries.begin(), (*line).entries.end());

				   // sort the lines
  sort (lines.begin(), lines.end());
  
  sorted = true;
};

	  

void ConstraintMatrix::clear () {
  lines = vector<ConstraintLine>();
  sorted = false;
};



void ConstraintMatrix::condense (const dSMatrixStruct &uncondensed,
				 dSMatrixStruct       &condensed) const {
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed.compressed == true, ExcMatrixNotClosed());
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
  unsigned int n_rows = (unsigned int)uncondensed.n_rows();
  for (unsigned int row=0; row!=n_rows; ++row)
    if (row == (*next_constraint).line)
      {
					 // this line is constrained
	new_line.push_back (-1);
					 // note that #lines# is ordered
	++next_constraint;
	++shift;
      }
    else
      new_line.push_back (row-shift);

 
  next_constraint = lines.begin();
  for (int row=0; row<uncondensed.rows; ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entries if column will not
				       // be condensed away, distribute
				       // otherwise
      for (int j=uncondensed.rowstart[row]; j<uncondensed.rowstart[row+1]; ++j)
	if (new_line[uncondensed.colnums[j]] != -1)
	  condensed.add (new_line[row], new_line[uncondensed.colnums[j]]);
	else 
	  {
					     // let c point to the constraint
					     // of this column
	    vector<ConstraintLine>::const_iterator c = lines.begin();
	    while ((*c).line != (unsigned int)uncondensed.colnums[j]) ++c;

	    for (unsigned int q=0; q!=(*c).entries.size(); ++q) 
	      condensed.add (new_line[row], new_line[(*c).entries[q].first]);
	  }
    else
				       // line must be distributed
      {
	for (int j=uncondensed.rowstart[row]; j<uncondensed.rowstart[row+1]; ++j)
					   // for each entry: distribute
	  if (new_line[uncondensed.colnums[j]] != -1)
					     // column is not constrained
	    for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q) 
		condensed.add (new_line[(*next_constraint).entries[q].first],
			       new_line[uncondensed.colnums[j]]);
	
	  else
					     // not only this line but
					     // also this col is constrained
	    {
					       // let c point to the constraint
					       // of this column
	      vector<ConstraintLine>::const_iterator c = lines.begin();
	      while ((*c).line != (unsigned int)uncondensed.colnums[j]) ++c;
	      
	      for (unsigned int p=0; p!=(*c).entries.size(); ++p)
		for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q)
		    condensed.add (new_line[(*next_constraint).entries[q].first],
				   new_line[(*c).entries[p].first]);
	    };
	
	++next_constraint;
      };

  condensed.compress();
};



void ConstraintMatrix::condense (dSMatrixStruct &sparsity) const {
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.compressed == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
				   //
				   // for gcc2.9: replace this by
				   // distribute(sparsity.n_rows(), -1)
  vector<int> distribute;
  distribute.reserve (sparsity.n_rows());
				   // replace the following line, since
				   // gcc2.8 chokes over it.
//  distribute.insert (distribute.end(), sparsity.n_rows(), -1);
  for (unsigned int i=0; i<(unsigned int)sparsity.n_rows(); ++i)
    distribute.push_back (-1);
  
  for (int c=0; c<(signed int)lines.size(); ++c)
    distribute[lines[c].line] = c;

  int n_rows = sparsity.n_rows();
  for (int row=0; row<n_rows; ++row)
    if (distribute[row] == -1)
				       // regular line. loop over cols
      for (int j=sparsity.rowstart[row]; j<sparsity.rowstart[row+1]; ++j)
					 // end of row reached?
	if (sparsity.colnums[j] == -1)
	  break;
	else
	  {
	    if (distribute[sparsity.colnums[j]] != -1)
					       // distribute entry at regular
					       // row #row# and irregular column
					       // sparsity.colnums[j]
	      for (unsigned int q=0;
		   q!=lines[distribute[sparsity.colnums[j]]].entries.size(); ++q) 
		sparsity.add (row,
			      lines[distribute[sparsity.colnums[j]]].entries[q].first);
	  }
    else
				       // row must be distributed
      for (int j=sparsity.rowstart[row]; j<sparsity.rowstart[row+1]; ++j)
					 // end of row reached?
	if (sparsity.colnums[j] == -1)
	  break;
	else
	  {
	    if (distribute[sparsity.colnums[j]] == -1)
					       // distribute entry at irregular
					       // row #row# and regular column
					       // sparsity.colnums[j]
	      for (unsigned int q=0;
		   q!=lines[distribute[row]].entries.size(); ++q) 
		sparsity.add (lines[distribute[row]].entries[q].first,
			      sparsity.colnums[j]);
	    else
					       // distribute entry at irregular
					       // row #row# and irregular column
					       // sparsity.colnums[j]
	      for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		for (unsigned int q=0;
		     q!=lines[distribute[sparsity.colnums[j]]].entries.size(); ++q)
		  sparsity.add (lines[distribute[row]].entries[p].first,
				lines[distribute[sparsity.colnums[j]]].entries[q].first);
	  };
  sparsity.compress();
};



void ConstraintMatrix::condense (const dSMatrix &uncondensed,
				 dSMatrix       &condensed) const {
  const dSMatrixStruct &uncondensed_struct = uncondensed.get_sparsity_pattern ();
  
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.compressed == true, ExcMatrixNotClosed());
  Assert (condensed.get_sparsity_pattern().compressed == true, ExcMatrixNotClosed());
  Assert (uncondensed_struct.n_rows() == uncondensed_struct.n_cols(),
	  ExcMatrixNotSquare());
  Assert (condensed.n() == condensed.m(),
	  ExcMatrixNotSquare());
  Assert ((unsigned int)condensed.n()+n_constraints() == (unsigned int)uncondensed.n(),
	  ExcWrongDimension());

				   // store for each line of the matrix
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  vector<int> new_line;

  new_line.reserve (uncondensed_struct.n_rows());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  unsigned int n_rows = (unsigned int)uncondensed_struct.n_rows();
  for (unsigned int row=0; row!=n_rows; ++row)
    if (row == (*next_constraint).line)
      {
					 // this line is constrained
	new_line.push_back (-1);
					 // note that #lines# is ordered
	++next_constraint;
	++shift;
      }
    else
      new_line.push_back (row-shift);

 
  next_constraint = lines.begin();
  for (int row=0; row<uncondensed_struct.rows; ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entries if column will not
				       // be condensed away, distribute
				       // otherwise
      for (int j=uncondensed_struct.rowstart[row]; j<uncondensed_struct.rowstart[row+1]; ++j)
	if (new_line[uncondensed_struct.colnums[j]] != -1)
	  condensed.add (new_line[row], new_line[uncondensed_struct.colnums[j]],
			 uncondensed.val[j]);
	else 
	  {
					     // let c point to the constraint
					     // of this column
	    vector<ConstraintLine>::const_iterator c = lines.begin();
	    while ((*c).line != (unsigned int)uncondensed_struct.colnums[j]) ++c;

	    for (unsigned int q=0; q!=(*c).entries.size(); ++q)
					       // distribute to rows with
					       // appropriate weight
	      condensed.add (new_line[row], new_line[(*c).entries[q].first],
			     uncondensed.val[j] * (*c).entries[q].second);
	  }
    else
				       // line must be distributed
      {
	for (int j=uncondensed_struct.rowstart[row]; j<uncondensed_struct.rowstart[row+1]; ++j)
					   // for each column: distribute
	  if (new_line[uncondensed_struct.colnums[j]] != -1)
					     // column is not constrained
	    for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q) 
		condensed.add (new_line[(*next_constraint).entries[q].first],
			       new_line[uncondensed_struct.colnums[j]],
			       uncondensed.val[j] * (*next_constraint).entries[q].second);
	
	  else
					     // not only this line but
					     // also this col is constrained
	    {
					       // let c point to the constraint
					       // of this column
	      vector<ConstraintLine>::const_iterator c = lines.begin();
	      while ((*c).line != (unsigned int)uncondensed_struct.colnums[j]) ++c;
	      
	      for (unsigned int p=0; p!=(*c).entries.size(); ++p)
		for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q)
		    condensed.add (new_line[(*next_constraint).entries[q].first],
				   new_line[(*c).entries[p].first],
				   uncondensed.val[j] *
				   (*next_constraint).entries[q].second *
				   (*c).entries[p].second);
	    };
	
	++next_constraint;
      };
};



void ConstraintMatrix::condense (dSMatrix &uncondensed) const {
  const dSMatrixStruct &sparsity = uncondensed.get_sparsity_pattern ();

  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.compressed == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcMatrixNotSquare());
  
				   // store for each index whether it
				   // must be distributed or not. If entry
				   // is -1, no distribution is necessary.
				   // otherwise, the number states which
				   // line in the constraint matrix handles
				   // this index
  				   //
				   // for gcc2.9: replace this by
				   // distribute(sparsity.n_rows(), -1)
  vector<int> distribute;
  distribute.reserve (sparsity.n_rows());
				   // replace the following line, since
				   // gcc2.8 chokes over it.
//  distribute.insert (distribute.end(), sparsity.n_rows(), -1);
  for (unsigned int i=0; i<(unsigned int)sparsity.n_rows(); ++i)
    distribute.push_back (-1);

  
  for (int c=0; c<(signed int)lines.size(); ++c)
    distribute[lines[c].line] = c;

  int n_rows = sparsity.n_rows();
  for (int row=0; row<n_rows; ++row) 
    {
      
    if (distribute[row] == -1)
				       // regular line. loop over cols
      for (int j=sparsity.rowstart[row]; j<sparsity.rowstart[row+1]; ++j)
					 // end of row reached?
	if (sparsity.colnums[j] == -1)
	  break;
	else
	  {
	    if (distribute[sparsity.colnums[j]] != -1)
					       // distribute entry at regular
					       // row #row# and irregular column
					       // sparsity.colnums[j]; set old
					       // entry to zero
	      {
		for (unsigned int q=0;
		     q!=lines[distribute[sparsity.colnums[j]]].entries.size(); ++q) 
		  uncondensed.add (row,
				   lines[distribute[sparsity.colnums[j]]].entries[q].first,
				   uncondensed.val[j] *
				   lines[distribute[sparsity.colnums[j]]].entries[q].second);
		uncondensed.val[j] = 0.;
	      };
	  }
    else
				       // row must be distributed
      for (int j=sparsity.rowstart[row]; j<sparsity.rowstart[row+1]; ++j)
					 // end of row reached?
	if (sparsity.colnums[j] == -1)
	  break;
	else
	  {
	    if (distribute[sparsity.colnums[j]] == -1)
					       // distribute entry at irregular
					       // row #row# and regular column
					       // sparsity.colnums[j]. set old
					       // entry to zero
	      {
		for (unsigned int q=0;
		     q!=lines[distribute[row]].entries.size(); ++q) 
		  uncondensed.add (lines[distribute[row]].entries[q].first,
				   sparsity.colnums[j],
				   uncondensed.val[j] *
				   lines[distribute[row]].entries[q].second);
		
		uncondensed.val[j] = 0.;
	      }
	    else
					       // distribute entry at irregular
					       // row #row# and irregular column
					       // sparsity.colnums[j]
					       // set old entry to one if on main
					       // diagonal, zero otherwise
	      {
		for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		  for (unsigned int q=0;
		       q!=lines[distribute[sparsity.colnums[j]]].entries.size(); ++q)
		    uncondensed.add (lines[distribute[row]].entries[p].first,
				     lines[distribute[sparsity.colnums[j]]].entries[q].first,
				     uncondensed.val[j] *
				     lines[distribute[row]].entries[p].second *
				     lines[distribute[sparsity.colnums[j]]].entries[q].second);
		
		uncondensed.val[j] = (row == sparsity.colnums[j] ?
				      1. : 0. );
	      };
	  };
    };
};






void ConstraintMatrix::condense (const dVector &uncondensed,
				 dVector       &condensed) const {
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert ((unsigned int)condensed.n()+n_constraints() == (unsigned int)uncondensed.n(),
	  ExcWrongDimension());
  
				   // store for each line of the vector
				   // its new line number
				   // after compression. If the shift is
				   // -1, this line will be condensed away
  vector<int> new_line;

  new_line.reserve (uncondensed.n());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  unsigned int n_rows = (unsigned int)uncondensed.n();
  for (unsigned int row=0; row!=n_rows; ++row)
    if (row == (*next_constraint).line)
      {
					 // this line is constrained
	new_line.push_back (-1);
					 // note that #lines# is ordered
	++next_constraint;
	++shift;
      }
    else
      new_line.push_back (row-shift);

 
  next_constraint = lines.begin();
  for (int row=0; row<uncondensed.n(); ++row)
    if (new_line[row] != -1)
				       // line not constrained
				       // copy entry
      condensed(new_line[row]) += uncondensed(row);

    else
				       // line must be distributed
      {
	for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q) 
	  condensed(new_line[(*next_constraint).entries[q].first])
	    +=
	    uncondensed(row) * (*next_constraint).entries[q].second;

	++next_constraint;
      };
};



void ConstraintMatrix::condense (dVector &vec) const {
  Assert (sorted == true, ExcMatrixNotClosed());
  
  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (unsigned int row=0; row<(unsigned int)vec.n(); ++row)
    if (row == (*next_constraint).line)
				       // line must be distributed
      {
	for (unsigned int q=0; q!=(*next_constraint).entries.size(); ++q) 
	  vec((*next_constraint).entries[q].first)
	    +=
	    vec(row) * (*next_constraint).entries[q].second;
					 // set entry to zero
	vec(row) = 0.;
	
	++next_constraint;
      };
};
  



void ConstraintMatrix::distribute (const dVector &condensed,
				   dVector       &uncondensed) const {
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert ((unsigned int)condensed.n()+n_constraints() == (unsigned int)uncondensed.n(),
	  ExcWrongDimension());
  
				   // store for each line of the new vector
				   // its old line number before
				   // distribution. If the shift is
				   // -1, this line was condensed away
  vector<int> old_line;

  old_line.reserve (uncondensed.n());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  unsigned int                           shift           = 0;
  unsigned int n_rows = (unsigned int)uncondensed.n();
  for (unsigned int row=0; row!=n_rows; ++row)
    if (row == (*next_constraint).line)
      {
					 // this line is constrained
	old_line.push_back (-1);
					 // note that #lines# is ordered
	++next_constraint;
	++shift;
      }
    else
      old_line.push_back (row-shift);

 
  next_constraint = lines.begin();
  for (unsigned int line=0; line<(unsigned int)uncondensed.n(); ++line) 
    if (old_line[line] != -1)
				       // line was not condensed away
      uncondensed(line) = condensed(old_line[line]);
    else
      {
					 // line was condensed away, create it newly
					 // first set it to zero
	uncondensed(line) = 0.;
					 // then add the different contributions
	for (unsigned int i=0; i<(*next_constraint).entries.size(); ++i)
	  uncondensed(line) += (condensed(old_line[(*next_constraint).entries[i].first]) *
				(*next_constraint).entries[i].second);
	++next_constraint;
      };
};



void ConstraintMatrix::distribute (dVector &vec) const {
  Assert (sorted == true, ExcMatrixNotClosed());

  vector<ConstraintLine>::const_iterator next_constraint = lines.begin();
  for (; next_constraint != lines.end(); ++next_constraint) 
    {
				       // make entry in line next_constraint.line
				       // first set it to zero
      vec((*next_constraint).line) = 0.;
				       // then add the different contributions
      for (unsigned int i=0; i<(*next_constraint).entries.size(); ++i)
	vec((*next_constraint).line) += (vec((*next_constraint).entries[i].first) *
					 (*next_constraint).entries[i].second);
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
};
