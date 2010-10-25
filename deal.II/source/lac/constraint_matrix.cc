//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/constraint_matrix.h>
#include <lac/constraint_matrix.templates.h>

#include <base/memory_consumption.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/block_vector.h>
#include <lac/block_sparse_matrix.h>
#include <lac/sparse_matrix_ez.h>
#include <lac/block_sparse_matrix_ez.h>
#include <lac/petsc_vector.h>
#include <lac/petsc_block_vector.h>
#include <lac/petsc_sparse_matrix.h>
#include <lac/petsc_block_sparse_matrix.h>
#include <lac/petsc_parallel_vector.h>
#include <lac/petsc_parallel_block_vector.h>
#include <lac/petsc_parallel_sparse_matrix.h>
#include <lac/petsc_parallel_block_sparse_matrix.h>
#include <lac/trilinos_vector.h>
#include <lac/trilinos_block_vector.h>
#include <lac/trilinos_sparse_matrix.h>
#include <lac/trilinos_block_sparse_matrix.h>
#include <lac/matrix_block.h>

#include <algorithm>
#include <numeric>
#include <set>

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif

DEAL_II_NAMESPACE_OPEN



				        // Static member variable
const Table<2,bool> ConstraintMatrix::default_empty_table = Table<2,bool>();



bool
ConstraintMatrix::check_zero_weight (const std::pair<unsigned int, double> &p)
{
  return (p.second == 0);
}



bool
ConstraintMatrix::ConstraintLine::operator < (const ConstraintLine &a) const
{
  return line < a.line;
}



bool
ConstraintMatrix::ConstraintLine::operator == (const ConstraintLine &a) const
{
  return line == a.line;
}



unsigned int
ConstraintMatrix::ConstraintLine::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (line) +
	  MemoryConsumption::memory_consumption (entries) +
	  MemoryConsumption::memory_consumption (inhomogeneity));
}



void
ConstraintMatrix::add_lines (const std::set<unsigned int> &lines)
{
  for (std::set<unsigned int>::const_iterator
	 i = lines.begin(); i != lines.end(); ++i)
    add_line (*i);
}



void
ConstraintMatrix::add_lines (const std::vector<bool> &lines)
{
  for (unsigned int i=0; i<lines.size(); ++i)
    if (lines[i] == true)
      add_line (i);
}



void
ConstraintMatrix::add_lines (const IndexSet &lines)
{
  for (unsigned int i=0; i<lines.n_elements(); ++i)
    add_line (lines.nth_index_in_set(i));
}



void
ConstraintMatrix::add_entries (const unsigned int                        line,
                               const std::vector<std::pair<unsigned int,double> > &col_val_pairs)
{
  Assert (sorted==false, ExcMatrixIsClosed());
  Assert (is_constrained(line), ExcLineInexistant(line));

  ConstraintLine * line_ptr = &lines[lines_cache[calculate_line_index(line)]];
  Assert (line_ptr->line == line, ExcInternalError());

				   // if the loop didn't break, then
				   // line_ptr must be begin().
				   // we have an error if that doesn't
				   // point to 'line' then

				   // if in debug mode, check whether an
				   // entry for this column already
				   // exists and if its the same as
				   // the one entered at present
				   //
				   // in any case: skip this entry if
				   // an entry for this column already
				   // exists, since we don't want to
				   // enter it twice
  for (std::vector<std::pair<unsigned int,double> >::const_iterator
         col_val_pair = col_val_pairs.begin();
       col_val_pair!=col_val_pairs.end(); ++col_val_pair)
    {
      Assert (line != col_val_pair->first,
	      ExcMessage ("Can't constrain a degree of freedom to itself"));

      for (std::vector<std::pair<unsigned int,double> >::const_iterator
             p=line_ptr->entries.begin();
	   p != line_ptr->entries.end(); ++p)
	if (p->first == col_val_pair->first)
	  {
					     // entry exists, break
					     // innermost loop
	    Assert (p->second == col_val_pair->second,
		    ExcEntryAlreadyExists(line, col_val_pair->first,
					  p->second, col_val_pair->second));
	    break;
	  }

      line_ptr->entries.push_back (*col_val_pair);
    }
}



void ConstraintMatrix::add_selected_constraints
  (const ConstraintMatrix &constraints,
   const IndexSet         &filter)
{
  Assert (filter.size() > constraints.lines.back().line,
	  ExcMessage ("Filter needs to be larger than constraint matrix size."));
  for (std::vector<ConstraintLine>::const_iterator line=constraints.lines.begin();
       line!=constraints.lines.end(); ++line)
    if (filter.is_element(line->line))
      {
	const unsigned int row = filter.index_within_set (line->line);
	add_line (row);
	set_inhomogeneity (row, line->inhomogeneity);
	for (unsigned int i=0; i<line->entries.size(); ++i)
	  if (filter.is_element(line->entries[i].first))
	    add_entry (row, filter.index_within_set (line->entries[i].first),
		       line->entries[i].second);
      }
}



void ConstraintMatrix::close ()
{
  if (sorted == true)
    return;

				   // sort the lines
  std::sort (lines.begin(), lines.end());

				   // update list of pointers and give the
				   // vector a sharp size since we won't
				   // modify the size any more after this
				   // point.
  {
    std::vector<unsigned int> new_lines (lines_cache.size(),
					 numbers::invalid_unsigned_int);
    unsigned int counter = 0;
    for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
	 line!=lines.end(); ++line, ++counter)
      new_lines[calculate_line_index(line->line)] = counter;
    std::swap (lines_cache, new_lines);
  }

				   // in debug mode: check whether we really
				   // set the pointers correctly.
  for (unsigned int i=0; i<lines_cache.size(); ++i)
    if (lines_cache[i] != numbers::invalid_unsigned_int)
      Assert (i == calculate_line_index(lines[lines_cache[i]].line),
	      ExcInternalError());

				   // first, strip zero entries, as we
				   // have to do that only once
  for (std::vector<ConstraintLine>::iterator line = lines.begin();
       line!=lines.end(); ++line)
				     // first remove zero
				     // entries. that would mean that
				     // in the linear constraint for a
				     // node, x_i = ax_1 + bx_2 + ...,
				     // another node times 0
				     // appears. obviously,
				     // 0*something can be omitted
    line->entries.erase (std::remove_if (line->entries.begin(),
					 line->entries.end(),
					 &check_zero_weight),
			 line->entries.end());

				   // replace references to dofs that
				   // are themselves constrained. note
				   // that because we may replace
				   // references to other dofs that
				   // may themselves be constrained to
				   // third ones, we have to iterate
				   // over all this until we replace
				   // no chains of constraints any
				   // more
                                   //
                                   // the iteration replaces
                                   // references to constrained
                                   // degrees of freedom by
                                   // second-order references. for
                                   // example if x3=x0/2+x2/2 and
                                   // x2=x0/2+x1/2, then the new list
                                   // will be x3=x0/2+x0/4+x1/4. note
                                   // that x0 appear twice. we will
                                   // throw this duplicate out in the
                                   // following step, where we sort
                                   // the list so that throwing out
                                   // duplicates becomes much more
                                   // efficient. also, we have to do
                                   // it only once, rather than in
                                   // each iteration
  unsigned int iteration = 0;
  while (true)
    {
      bool chained_constraint_replaced = false;

      for (std::vector<ConstraintLine>::iterator line = lines.begin();
	   line!=lines.end(); ++line)
	{
					   // loop over all entries of
					   // this line (including
					   // ones that we have
					   // appended in this go
					   // around) and see whether
					   // they are further
					   // constrained. ignore
					   // elements that we don't
					   // store on the current
					   // processor
	  unsigned int entry = 0;
	  while (entry < line->entries.size())
	    if (((local_lines.size() == 0)
		 ||
		 (local_lines.is_element(line->entries[entry].first)))
		&&
		is_constrained (line->entries[entry].first))
	      {
						 // ok, this entry is
						 // further
						 // constrained:
		chained_constraint_replaced = true;

						 // look up the chain
						 // of constraints for
						 // this entry
		const unsigned int dof_index = line->entries[entry].first;
		const double       weight = line->entries[entry].second;

		Assert (dof_index != line->line,
			ExcMessage ("Cycle in constraints detected!"));

		const ConstraintLine * constrained_line =
		  &lines[lines_cache[calculate_line_index(dof_index)]];
		Assert (constrained_line->line == dof_index,
			ExcInternalError());

						 // now we have to
						 // replace an entry
						 // by its
						 // expansion. we do
						 // that by
						 // overwriting the
						 // entry by the first
						 // entry of the
						 // expansion and
						 // adding the
						 // remaining ones to
						 // the end, where we
						 // will later process
						 // them once more
						 //
						 // we can of course
						 // only do that if
						 // the DoF that we
						 // are currently
						 // handle is
						 // constrained by a
						 // linear combination
						 // of other dofs:
		if (constrained_line->entries.size() > 0)
		  {
		    for (unsigned int i=0; i<constrained_line->entries.size(); ++i)
		      Assert (dof_index != constrained_line->entries[i].first,
			      ExcMessage ("Cycle in constraints detected!"));

						     // replace first
						     // entry, then tack
						     // the rest to the
						     // end of the list
		    line->entries[entry] =
		      std::make_pair (constrained_line->entries[0].first,
				      constrained_line->entries[0].second *
				      weight);

		    for (unsigned int i=1; i<constrained_line->entries.size(); ++i)
		      line->entries
			.push_back (std::make_pair (constrained_line->entries[i].first,
						    constrained_line->entries[i].second *
						    weight));
		  }
		else
						   // the DoF that we
						   // encountered is not
						   // constrained by a
						   // linear combination of
						   // other dofs but is
						   // equal to zero
						   // (i.e. its chain of
						   // entries is empty). in
						   // that case, we can't
						   // just overwrite the
						   // current entry, but we
						   // have to actually
						   // eliminate it
		  {
		    line->entries.erase (line->entries.begin()+entry);
		  }

		line->inhomogeneity += constrained_line->inhomogeneity *
		                       weight;

						 // now that we're here, do
						 // not increase index by
						 // one but rather make
						 // another pass for the
						 // present entry because we
						 // have replaced the
						 // present entry by another
						 // one, or because we have
						 // deleted it and shifted
						 // all following ones one
						 // forward
	      }
	    else
					       // entry not further
					       // constrained. just move
					       // ahead by one
	      ++entry;
	}

				       // if we didn't do anything in
				       // this round, then quit the
				       // loop
      if (chained_constraint_replaced == false)
	break;

				       // increase iteration count. note
				       // that we should not iterate more
				       // times than there are constraints,
				       // since this puts a natural upper
				       // bound on the length of constraint
				       // chains
      ++iteration;
      Assert (iteration <= lines.size(),
	      ExcInternalError());
    }

				   // finally sort the entries and re-scale
				   // them if necessary. in this step, we
				   // also throw out duplicates as mentioned
				   // above
  for (std::vector<ConstraintLine>::iterator line = lines.begin();
       line!=lines.end(); ++line)
    {
      std::sort (line->entries.begin(), line->entries.end());

                                       // loop over the now sorted list and
                                       // see whether any of the entries
                                       // references the same dofs more than
                                       // once
      for (unsigned int i=1; i<line->entries.size(); ++i)
        if (line->entries[i].first == line->entries[i-1].first)
          {
                                             // ok, we've found a
                                             // duplicate. go on to count
                                             // how many duplicates there
                                             // are so that we can allocate
                                             // the right amount of memory
            unsigned int duplicates = 1;
            for (unsigned int j=i+1; j<line->entries.size(); ++j)
              if (line->entries[j].first == line->entries[j-1].first)
                ++duplicates;

            std::vector<std::pair<unsigned int, double> > new_entries;
            new_entries.reserve (line->entries.size() - duplicates);

                                             // now copy all entries
                                             // with unique keys. copy
                                             // verbatim the entries
                                             // at the front for which
                                             // we have already
                                             // determined that they
                                             // have no duplicate
                                             // copies
            new_entries.insert (new_entries.begin(),
                                line->entries.begin(),
                                line->entries.begin()+i);
                                             // now for the rest
            for (unsigned int j=i; j<line->entries.size(); ++j)
              if (line->entries[j].first == line->entries[j-1].first)
                {
                  Assert (new_entries.back().first == line->entries[j].first,
                          ExcInternalError());
                  new_entries.back().second += line->entries[j].second;
                }
              else
                new_entries.push_back (line->entries[j]);

            Assert (new_entries.size() == line->entries.size() - duplicates,
                    ExcInternalError());

                                             // make sure there are
                                             // really no duplicates
                                             // left and that the list
                                             // is still sorted
            for (unsigned int j=1; j<new_entries.size(); ++j)
              {
                Assert (new_entries[j].first != new_entries[j-1].first,
                        ExcInternalError());
                Assert (new_entries[j].first > new_entries[j-1].first,
                        ExcInternalError());
              }


                                             // replace old list of
                                             // constraints for this
                                             // dof by the new one and
                                             // quit loop
            line->entries.swap (new_entries);
            break;
          }

				       // finally do the following
				       // check: if the sum of
				       // weights for the
				       // constraints is close to
				       // one, but not exactly
				       // one, then rescale all
				       // the weights so that they
				       // sum up to 1. this adds a
				       // little numerical
				       // stability and avoids all
				       // sorts of problems where
				       // the actual value is
				       // close to, but not quite
				       // what we expected
				       //
				       // the case where the
				       // weights don't quite sum
				       // up happens when we
				       // compute the
				       // interpolation weights
				       // "on the fly", i.e. not
				       // from precomputed
				       // tables. in this case,
				       // the interpolation
				       // weights are also subject
				       // to round-off
      double sum = 0;
      for (unsigned int i=0; i<line->entries.size(); ++i)
	sum += line->entries[i].second;
      if ((sum != 1.0) && (std::fabs (sum-1.) < 1.e-13))
	{
	  for (unsigned int i=0; i<line->entries.size(); ++i)
	    line->entries[i].second /= sum;
	  line->inhomogeneity /= sum;
	}
    }

#ifdef DEBUG
				   // if in debug mode: check that no
				   // dof is constraint to another dof
				   // that is also
				   // constrained. exclude dofs from
				   // this check whose constraint
				   // lines would not be store on the
				   // local processor
  for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line)
    for (std::vector<std::pair<unsigned int,double> >::const_iterator
	   entry=line->entries.begin();
	 entry!=line->entries.end(); ++entry)
      if ((local_lines.size() == 0)
	  ||
	  (local_lines.is_element(entry->first)))
	{
					   // make sure that
					   // entry->first is not the
					   // index of a line itself
	  const bool is_circle = is_constrained(entry->first);
	  Assert (is_circle == false,
		  ExcDoFConstrainedToConstrainedDoF(line->line, entry->first));
	}
#endif

  sorted = true;
}



void
ConstraintMatrix::merge (const ConstraintMatrix &other_constraints,
			 const MergeConflictBehavior merge_conflict_behavior)
{
				   //TODO: this doesn't work with IndexSets yet. [TH]
  AssertThrow(local_lines.size()==0, ExcNotImplemented());
  AssertThrow(other_constraints.local_lines.size()==0, ExcNotImplemented());

				   // store the previous state with
				   // respect to sorting
  const bool object_was_sorted = sorted;
  sorted = false;

  if (other_constraints.lines_cache.size() > lines_cache.size())
    lines_cache.resize(other_constraints.lines_cache.size(),
		       numbers::invalid_unsigned_int);

				   // first action is to fold into the present
				   // object possible constraints in the
				   // second object. we don't strictly need to
				   // do this any more since the
				   // ConstraintMatrix has learned to deal
				   // with chains of constraints in the
				   // close() function, but we have
				   // traditionally done this and it's not
				   // overly hard to do.
				   //
				   // for this, loop over all
				   // constraints and replace the
				   // constraint lines with a new one
				   // where constraints are replaced
				   // if necessary.
  std::vector<std::pair<unsigned int,double> > tmp;
  for (std::vector<ConstraintLine>::iterator line=lines.begin();
       line!=lines.end(); ++line)
    {
      tmp.clear ();
      for (unsigned int i=0; i<line->entries.size(); ++i)
	{
					   // if the present dof is not
					   // constrained, or if we won't take
					   // the constraint from the other
					   // object, then simply copy it over
	  if (!other_constraints.is_constrained(line->entries[i].first)
	      ||
	      ((merge_conflict_behavior != right_object_wins)
	       &&
	       other_constraints.is_constrained(line->entries[i].first)
	       &&
	       this->is_constrained(line->entries[i].first)))
	    tmp.push_back(line->entries[i]);
	  else
					     // otherwise resolve
					     // further constraints by
					     // replacing the old
					     // entry by a sequence of
					     // new entries taken from
					     // the other object, but
					     // with multiplied
					     // weights
	    {
	      const std::vector<std::pair<unsigned int,double> >*
		other_line
		= other_constraints.get_constraint_entries (line->entries[i].first);
	      Assert (other_line != 0,
		      ExcInternalError());
	      
	      const double weight = line->entries[i].second;
	      
	      for (std::vector<std::pair<unsigned int, double> >::const_iterator
		     j=other_line->begin();
		   j!=other_line->end(); ++j)
		tmp.push_back (std::make_pair(j->first, j->second*weight));

	      line->inhomogeneity += other_constraints.get_inhomogeneity(line->entries[i].first) *
				     line->entries[i].second;
	    }
	}
				       // finally exchange old and
				       // newly resolved line
      line->entries.swap (tmp);
    }



				   // next action: append those lines at the
				   // end that we want to add
  for (std::vector<ConstraintLine>::const_iterator
	 line=other_constraints.lines.begin();
       line!=other_constraints.lines.end(); ++line)
    if (!is_constrained(line->line))
      lines.push_back (*line);
    else
      {
					 // the constrained dof we want to
					 // copy from the other object is also
					 // constrained here. let's see what
					 // we should do with that
	switch (merge_conflict_behavior)
	  {
	    case no_conflicts_allowed:
		  AssertThrow (false,
			       ExcDoFIsConstrainedFromBothObjects (line->line));
		  break;
		  
	    case left_object_wins:
						   // ignore this constraint
		  break;
		  
	    case right_object_wins:
						   // we need to replace the
						   // existing constraint by
						   // the one from the other
						   // object
		  lines[lines_cache[calculate_line_index(line->line)]].entries
		    = line->entries;
		  lines[lines_cache[calculate_line_index(line->line)]].inhomogeneity
		    = line->inhomogeneity;
		  break;

	    default:
		  Assert (false, ExcNotImplemented());
	  }
      }
  
				   // if the object was sorted before,
				   // then make sure it is so
				   // afterwards as well. otherwise
				   // leave everything in the unsorted
				   // state
  unsigned int counter = 0;
  for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line, ++counter)
    lines_cache[line->line] = counter;
  if (object_was_sorted == true)
    close ();
}



void ConstraintMatrix::shift (const unsigned int offset)
{
				   //TODO: this doesn't work with IndexSets yet. [TH]
  AssertThrow(local_lines.size()==0, ExcNotImplemented());

  lines_cache.insert (lines_cache.begin(), offset,
		      numbers::invalid_unsigned_int);

  for (std::vector<ConstraintLine>::iterator i = lines.begin();
       i != lines.end(); ++i)
    {
      i->line += offset;
      for (std::vector<std::pair<unsigned int,double> >::iterator
	     j = i->entries.begin();
	   j != i->entries.end(); ++j)
	j->first += offset;
    }
}



void ConstraintMatrix::clear ()
{
  {
    std::vector<ConstraintLine> tmp;
    lines.swap (tmp);
  }

  {
    std::vector<unsigned int> tmp;
    lines_cache.swap (tmp);
  }

  sorted = false;
}



void ConstraintMatrix::reinit (const IndexSet & local_constraints)
{
  local_lines = local_constraints;
  clear();
}



void ConstraintMatrix::condense (const SparsityPattern &uncondensed,
				 SparsityPattern       &condensed) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (uncondensed.is_compressed() == true, ExcMatrixNotClosed());
  Assert (uncondensed.n_rows() == uncondensed.n_cols(),
	  ExcNotQuadratic());


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
}



void ConstraintMatrix::condense (SparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute(sparsity.n_rows(),
                                       numbers::invalid_unsigned_int);

  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_unsigned_int)
	{
					   // regular line. loop over cols all
					   // valid cols. note that this
					   // changes the line we are
					   // presently working on: we add
					   // additional entries. these are
					   // put to the end of the
					   // row. however, as constrained
					   // nodes cannot be constrained to
					   // other constrained nodes, nothing
					   // will happen if we run into these
					   // added nodes, as they can't be
					   // distributed further. we might
					   // store the position of the last
					   // old entry and stop work there,
					   // but since operating on the newly
					   // added ones only takes two
					   // comparisons (column index valid,
					   // distribute[column] necessarily
					   // ==numbers::invalid_unsigned_int),
					   // it is cheaper to not do so and
					   // run right until the end of the
					   // line
          for (SparsityPattern::iterator entry = sparsity.begin(row);
               ((entry != sparsity.end(row)) &&
                entry->is_valid_entry());
               ++entry)
	    {
	      const unsigned int column = entry->column();

              if (distribute[column] != numbers::invalid_unsigned_int)
                {
                                                   // distribute entry
                                                   // at regular row
                                                   // @p{row} and
                                                   // irregular column
                                                   // sparsity.colnums[j]
                  for (unsigned int q=0;
                       q!=lines[distribute[column]].entries.size();
                       ++q)
                    sparsity.add (row,
                                  lines[distribute[column]].entries[q].first);
                }
	    }
	}
      else
					 // row must be
					 // distributed. note that
					 // here the present row is
					 // not touched (unlike above)
	{
          for (SparsityPattern::iterator entry = sparsity.begin(row);
               (entry != sparsity.end(row)) && entry->is_valid_entry(); ++entry)
            {
              const unsigned int column = entry->column();
              if (distribute[column] == numbers::invalid_unsigned_int)
                                                 // distribute entry at irregular
                                                 // row @p{row} and regular column
                                                 // sparsity.colnums[j]
                for (unsigned int q=0;
                     q!=lines[distribute[row]].entries.size(); ++q)
                  sparsity.add (lines[distribute[row]].entries[q].first,
                                column);
              else
                                                 // distribute entry at irregular
                                                 // row @p{row} and irregular column
                                                 // sparsity.get_column_numbers()[j]
                for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
                  for (unsigned int q=0;
                       q!=lines[distribute[column]].entries.size(); ++q)
                    sparsity.add (lines[distribute[row]].entries[p].first,
                                  lines[distribute[column]].entries[q].first);
            }
	}
    }

  sparsity.compress();
}



void ConstraintMatrix::condense (CompressedSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute(sparsity.n_rows(),
                                       numbers::invalid_unsigned_int);

  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // cols. note that as we
					 // proceed to distribute
					 // cols, the loop may get
					 // longer
	for (unsigned int j=0; j<sparsity.row_length(row); ++j)
	  {
	    const unsigned int column = sparsity.column_number(row,j);

 	    if (distribute[column] != numbers::invalid_unsigned_int)
	      {
						 // distribute entry
						 // at regular row
						 // @p{row} and
						 // irregular column
						 // column. note that
						 // this changes the
						 // line we are
						 // presently working
						 // on: we add
						 // additional
						 // entries. if we add
						 // another entry at a
						 // column behind the
						 // present one, we
						 // will encounter it
						 // later on (but
						 // since it can't be
						 // further
						 // constrained, won't
						 // have to do
						 // anything about
						 // it). if we add it
						 // up front of the
						 // present column, we
						 // will find the
						 // present column
						 // later on again as
						 // it was shifted
						 // back (again
						 // nothing happens,
						 // in particular no
						 // endless loop, as
						 // when we encounter
						 // it the second time
						 // we won't be able
						 // to add more
						 // entries as they
						 // all already exist,
						 // but we do the same
						 // work more often
						 // than necessary,
						 // and the loop gets
						 // longer), so move
						 // the cursor one to
						 // the right in the
						 // case that we add
						 // an entry up front
						 // that did not exist
						 // before. check
						 // whether it existed
						 // before by tracking
						 // the length of this
						 // row
		unsigned int old_rowlength = sparsity.row_length(row);
		for (unsigned int q=0;
		     q!=lines[distribute[column]].entries.size();
		     ++q)
		  {
		    const unsigned int
		      new_col = lines[distribute[column]].entries[q].first;

		    sparsity.add (row, new_col);

		    const unsigned int new_rowlength = sparsity.row_length(row);
		    if ((new_col < column) && (old_rowlength != new_rowlength))
		      ++j;
		    old_rowlength = new_rowlength;
		  };
	      };
	  }
      else
					 // row must be distributed
	for (unsigned int j=0; j<sparsity.row_length(row); ++j)
	  {
	    const unsigned int column = sparsity.column_number(row,j);

	    if (distribute[column] == numbers::invalid_unsigned_int)
					       // distribute entry at irregular
					       // row @p{row} and regular column
					       // sparsity.colnums[j]
	      for (unsigned int q=0;
		   q!=lines[distribute[row]].entries.size(); ++q)
		sparsity.add (lines[distribute[row]].entries[q].first,
			      column);
	    else
					       // distribute entry at irregular
					       // row @p{row} and irregular column
					       // sparsity.get_column_numbers()[j]
	      for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		for (unsigned int q=0;
		     q!=lines[distribute[sparsity.column_number(row,j)]]
				    .entries.size(); ++q)
		  sparsity.add (lines[distribute[row]].entries[p].first,
				lines[distribute[sparsity.column_number(row,j)]]
				.entries[q].first);
	  };
    };
}



void ConstraintMatrix::condense (CompressedSetSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute(sparsity.n_rows(),
                                       numbers::invalid_unsigned_int);

  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_unsigned_int)
	{
					   // regular line. loop over
					   // cols. note that as we proceed to
					   // distribute cols, the loop may
					   // get longer
	  CompressedSetSparsityPattern::row_iterator col_num = sparsity.row_begin (row);

	  for (; col_num != sparsity.row_end (row); ++col_num)
	    {
	      const unsigned int column = *col_num;

	      if (distribute[column] != numbers::invalid_unsigned_int)
		{
		  // row
		  for (unsigned int q=0;
		       q!=lines[distribute[column]].entries.size();
		       ++q)
		    {
		      const unsigned int
			new_col = lines[distribute[column]].entries[q].first;

		      sparsity.add (row, new_col);
		    }
		}
	    }
	}
      else
	// row must be distributed
	{
	  CompressedSetSparsityPattern::row_iterator col_num = sparsity.row_begin (row);

	  for (; col_num != sparsity.row_end (row); ++col_num)
	    {
	      const unsigned int column = *col_num;

	      if (distribute[column] == numbers::invalid_unsigned_int)
		// distribute entry at irregular
		// row @p{row} and regular column
		// sparsity.colnums[j]
		for (unsigned int q=0;
		     q!=lines[distribute[row]].entries.size(); ++q)
		  sparsity.add (lines[distribute[row]].entries[q].first,
				column);
	      else
		// distribute entry at irregular
		// row @p{row} and irregular column
		// sparsity.get_column_numbers()[j]
		for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		  for (unsigned int q=0;
		       q!=lines[distribute[column]]
			 .entries.size(); ++q)
		    sparsity.add (lines[distribute[row]].entries[p].first,
				  lines[distribute[column]]
				  .entries[q].first);
	    };
	}
    };
}



void ConstraintMatrix::condense (CompressedSimpleSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute(sparsity.n_rows(),
                                       numbers::invalid_unsigned_int);

  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // cols. note that as we
					 // proceed to distribute
					 // cols, the loop may get
					 // longer
	for (unsigned int j=0; j<sparsity.row_length(row); ++j)
	  {
	    const unsigned int column = sparsity.column_number(row,j);

 	    if (distribute[column] != numbers::invalid_unsigned_int)
	      {
						 // distribute entry
						 // at regular row
						 // @p{row} and
						 // irregular column
						 // column. note that
						 // this changes the
						 // line we are
						 // presently working
						 // on: we add
						 // additional
						 // entries. if we add
						 // another entry at a
						 // column behind the
						 // present one, we
						 // will encounter it
						 // later on (but
						 // since it can't be
						 // further
						 // constrained, won't
						 // have to do
						 // anything about
						 // it). if we add it
						 // up front of the
						 // present column, we
						 // will find the
						 // present column
						 // later on again as
						 // it was shifted
						 // back (again
						 // nothing happens,
						 // in particular no
						 // endless loop, as
						 // when we encounter
						 // it the second time
						 // we won't be able
						 // to add more
						 // entries as they
						 // all already exist,
						 // but we do the same
						 // work more often
						 // than necessary,
						 // and the loop gets
						 // longer), so move
						 // the cursor one to
						 // the right in the
						 // case that we add
						 // an entry up front
						 // that did not exist
						 // before. check
						 // whether it existed
						 // before by tracking
						 // the length of this
						 // row
		unsigned int old_rowlength = sparsity.row_length(row);
		for (unsigned int q=0;
		     q!=lines[distribute[column]].entries.size();
		     ++q)
		  {
		    const unsigned int
		      new_col = lines[distribute[column]].entries[q].first;

		    sparsity.add (row, new_col);

		    const unsigned int new_rowlength = sparsity.row_length(row);
		    if ((new_col < column) && (old_rowlength != new_rowlength))
		      ++j;
		    old_rowlength = new_rowlength;
		  };
	      };
	  }
      else
					 // row must be distributed
	for (unsigned int j=0; j<sparsity.row_length(row); ++j)
	  {
	    const unsigned int column = sparsity.column_number(row,j);

	    if (distribute[column] == numbers::invalid_unsigned_int)
					       // distribute entry at irregular
					       // row @p{row} and regular column
					       // sparsity.colnums[j]
	      for (unsigned int q=0;
		   q!=lines[distribute[row]].entries.size(); ++q)
		sparsity.add (lines[distribute[row]].entries[q].first,
			      column);
	    else
					       // distribute entry at irregular
					       // row @p{row} and irregular column
					       // sparsity.get_column_numbers()[j]
	      for (unsigned int p=0; p!=lines[distribute[row]].entries.size(); ++p)
		for (unsigned int q=0;
		     q!=lines[distribute[sparsity.column_number(row,j)]]
				    .entries.size(); ++q)
		  sparsity.add (lines[distribute[row]].entries[p].first,
				lines[distribute[sparsity.column_number(row,j)]]
				.entries[q].first);
	  };
    };
}



void ConstraintMatrix::condense (BlockSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcNotQuadratic());

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();

  const unsigned int n_blocks = sparsity.n_block_rows();

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);

  for (unsigned int c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const unsigned int n_rows = sparsity.n_rows();
  for (unsigned int row=0; row<n_rows; ++row)
    {
				       // get index of this row
				       // within the blocks
      const std::pair<unsigned int,unsigned int>
	block_index = index_mapping.global_to_local(row);
      const unsigned int block_row = block_index.first;

      if (distribute[row] == numbers::invalid_unsigned_int)
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

              for (SparsityPattern::const_iterator
                     entry = block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                     entry->is_valid_entry();
                   ++entry)
                {
                  const unsigned int global_col
                    = index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] != numbers::invalid_unsigned_int)
                                                     // distribute entry at regular
                                                     // row @p{row} and irregular column
                                                     // global_col
                    {
                      for (unsigned int q=0;
                           q!=lines[distribute[global_col]].entries.size(); ++q)
                        sparsity.add (row,
                                      lines[distribute[global_col]].entries[q].first);
                    }
                }
	    }
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

              for (SparsityPattern::const_iterator
                     entry = block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                     entry->is_valid_entry();
                   ++entry)
                {
                  const unsigned int global_col
                    = index_mapping.local_to_global (block_col, entry->column());

                  if (distribute[global_col] == numbers::invalid_unsigned_int)
                                                     // distribute entry at irregular
                                                     // row @p{row} and regular column
                                                     // global_col.
                    {
                      for (unsigned int q=0; q!=lines[distribute[row]].entries.size(); ++q)
                        sparsity.add (lines[distribute[row]].entries[q].first, global_col);
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
                    }
                }
	    }
	}
    }

  sparsity.compress();
}



void ConstraintMatrix::condense (BlockCompressedSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcNotQuadratic());

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();

  const unsigned int n_blocks = sparsity.n_block_rows();

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);

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
      const unsigned int local_row = block_index.second;

      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // all columns and see
					 // whether this column must
					 // be distributed. note that
					 // as we proceed to
					 // distribute cols, the loop
					 // over cols may get longer.
					 //
					 // don't try to be clever
					 // here as in the algorithm
					 // for the
					 // CompressedSparsityPattern,
					 // as that would be much more
					 // complicated here. after
					 // all, we know that
					 // compressed patterns are
					 // inefficient...
	{

					   // to loop over all entries
					   // in this row, we have to
					   // loop over all blocks in
					   // this blockrow and the
					   // corresponding row
					   // therein
	  for (unsigned int block_col=0; block_col<n_blocks; ++block_col)
	    {
	      const CompressedSparsityPattern &
		block_sparsity = sparsity.block(block_row, block_col);

	      for (unsigned int j=0; j<block_sparsity.row_length(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global(block_col,
						    block_sparsity.column_number(local_row,j));

		  if (distribute[global_col] != numbers::invalid_unsigned_int)
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
	      const CompressedSparsityPattern &
		block_sparsity = sparsity.block(block_row,block_col);

	      for (unsigned int j=0; j<block_sparsity.row_length(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global (block_col,
						     block_sparsity.column_number(local_row,j));

		  if (distribute[global_col] == numbers::invalid_unsigned_int)
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
}



void ConstraintMatrix::condense (BlockCompressedSetSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcNotQuadratic());

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();

  const unsigned int n_blocks = sparsity.n_block_rows();

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);

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
      const unsigned int local_row = block_index.second;

      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // all columns and see
					 // whether this column must
					 // be distributed. note that
					 // as we proceed to
					 // distribute cols, the loop
					 // over cols may get longer.
					 //
					 // don't try to be clever
					 // here as in the algorithm
					 // for the
					 // CompressedSparsityPattern,
					 // as that would be much more
					 // complicated here. after
					 // all, we know that
					 // compressed patterns are
					 // inefficient...
	{

					   // to loop over all entries
					   // in this row, we have to
					   // loop over all blocks in
					   // this blockrow and the
					   // corresponding row
					   // therein
	  for (unsigned int block_col=0; block_col<n_blocks; ++block_col)
	    {
	      const CompressedSetSparsityPattern &
		block_sparsity = sparsity.block(block_row, block_col);

	      for (CompressedSetSparsityPattern::row_iterator
		     j = block_sparsity.row_begin(local_row);
		   j != block_sparsity.row_end(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global(block_col, *j);

		  if (distribute[global_col] != numbers::invalid_unsigned_int)
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
	      const CompressedSetSparsityPattern &
		block_sparsity = sparsity.block(block_row,block_col);

	      for (CompressedSetSparsityPattern::row_iterator
		     j = block_sparsity.row_begin(local_row);
		   j != block_sparsity.row_end(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global (block_col, *j);

		  if (distribute[global_col] == numbers::invalid_unsigned_int)
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
}



void ConstraintMatrix::condense (BlockCompressedSimpleSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.n_block_rows() == sparsity.n_block_cols(),
	  ExcNotQuadratic());
  Assert (sparsity.get_column_indices() == sparsity.get_row_indices(),
	  ExcNotQuadratic());

  const BlockIndices &
    index_mapping = sparsity.get_column_indices();

  const unsigned int n_blocks = sparsity.n_block_rows();

				   // store for each index whether it must be
				   // distributed or not. If entry is
				   // numbers::invalid_unsigned_int,
				   // no distribution is necessary.
				   // otherwise, the number states which line
				   // in the constraint matrix handles this
				   // index
  std::vector<unsigned int> distribute (sparsity.n_rows(),
                                        numbers::invalid_unsigned_int);

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
      const unsigned int local_row = block_index.second;

      if (distribute[row] == numbers::invalid_unsigned_int)
					 // regular line. loop over
					 // all columns and see
					 // whether this column must
					 // be distributed. note that
					 // as we proceed to
					 // distribute cols, the loop
					 // over cols may get longer.
					 //
					 // don't try to be clever
					 // here as in the algorithm
					 // for the
					 // CompressedSparsityPattern,
					 // as that would be much more
					 // complicated here. after
					 // all, we know that
					 // compressed patterns are
					 // inefficient...
	{

					   // to loop over all entries
					   // in this row, we have to
					   // loop over all blocks in
					   // this blockrow and the
					   // corresponding row
					   // therein
	  for (unsigned int block_col=0; block_col<n_blocks; ++block_col)
	    {
	      const CompressedSimpleSparsityPattern &
		block_sparsity = sparsity.block(block_row, block_col);

	      for (unsigned int j=0; j<block_sparsity.row_length(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global(block_col,
						    block_sparsity.column_number(local_row,j));

		  if (distribute[global_col] != numbers::invalid_unsigned_int)
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
	      const CompressedSimpleSparsityPattern &
		block_sparsity = sparsity.block(block_row,block_col);

	      for (unsigned int j=0; j<block_sparsity.row_length(local_row); ++j)
		{
		  const unsigned int global_col
		    = index_mapping.local_to_global (block_col,
						     block_sparsity.column_number(local_row,j));

		  if (distribute[global_col] == numbers::invalid_unsigned_int)
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
}



#ifdef DEAL_II_USE_TRILINOS

				   // this is a specialization for a
				   // parallel (non-block) Trilinos
				   // vector. The basic idea is to just work
				   // on the local range of the vector. But
				   // we need access to values that the
				   // local nodes are constrained to.

template<>
void
ConstraintMatrix::distribute (TrilinosWrappers::MPI::Vector &vec) const
{
				   //TODO: not implemented yet, we need to fix
				   //LocalRange() first to only include
				   //"owned" indices. For this we need to keep
				   //track of the owned indices, because
				   //Trilinos doesn't. Use same constructor
				   //interface as in PETSc with two IndexSets!
  AssertThrow (vec.vector_partitioner().IsOneToOne(),
	       ExcMessage ("Distribute does not work on vectors with overlapping parallel partitioning."));


  typedef std::vector<ConstraintLine>::const_iterator constraint_iterator;
  ConstraintLine index_comparison;
  index_comparison.line = vec.local_range().first;
  const constraint_iterator begin_my_constraints =
    std::lower_bound (lines.begin(),lines.end(),index_comparison);

  index_comparison.line = vec.local_range().second;
  const constraint_iterator end_my_constraints
    = std::lower_bound(lines.begin(),lines.end(),index_comparison);

				   // Here we search all the indices that we
				   // need to have read-access to - the
				   // local nodes and all the nodes that the
				   // constraints indicate.
  IndexSet my_indices (vec.size());
  {
    const std::pair<unsigned int, unsigned int>
      local_range = vec.local_range();

    my_indices.add_range (local_range.first, local_range.second);

    std::set<unsigned int> individual_indices;
    for (constraint_iterator it = begin_my_constraints;
	 it != end_my_constraints; ++it)
      for (unsigned int i=0; i<it->entries.size(); ++i)
	if ((it->entries[i].first < local_range.first)
	    ||
	    (it->entries[i].first >= local_range.second))
	  individual_indices.insert (it->entries[i].first);

    my_indices.add_indices (individual_indices.begin(),
			    individual_indices.end());
  }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  const Epetra_MpiComm *mpi_comm
    = dynamic_cast<const Epetra_MpiComm*>(&vec.trilinos_vector().Comm());

  Assert (mpi_comm != 0, ExcInternalError());

  TrilinosWrappers::MPI::Vector vec_distribute
    (my_indices.make_trilinos_map (mpi_comm->Comm(), true));
#else
  TrilinosWrappers::MPI::Vector vec_distribute
    (my_indices.make_trilinos_map (MPI_COMM_WORLD, true));
#endif

				   // here we import the data
  vec_distribute.reinit(vec,false,true);

  for (constraint_iterator it = begin_my_constraints;
       it != end_my_constraints; ++it)
    {
				       // fill entry in line
				       // next_constraint.line by adding the
				       // different contributions
      double new_value = it->inhomogeneity;
      for (unsigned int i=0; i<it->entries.size(); ++i)
	new_value += (vec_distribute(it->entries[i].first) *
                      it->entries[i].second);
      vec(it->line) = new_value;
    }

				   // some processes might not apply
				   // constraints, so we need to explicitly
				   // state, that the others are doing an
				   // insert here:
  vec.compress (Insert);
}



template<>
void
ConstraintMatrix::distribute (TrilinosWrappers::MPI::BlockVector &vec) const
{
  IndexSet my_indices (vec.size());
  for (unsigned int block=0; block<vec.n_blocks(); ++block)
    {
      typedef std::vector<ConstraintLine>::const_iterator constraint_iterator;
      ConstraintLine index_comparison;
      index_comparison.line = vec.block(block).local_range().first
	+vec.get_block_indices().block_start(block);
      const constraint_iterator begin_my_constraints =
	std::lower_bound (lines.begin(),lines.end(),index_comparison);

      index_comparison.line = vec.block(block).local_range().second
	+vec.get_block_indices().block_start(block);

      const constraint_iterator end_my_constraints
	= std::lower_bound(lines.begin(),lines.end(),index_comparison);

				   // Here we search all the indices that we
				   // need to have read-access to - the local
				   // nodes and all the nodes that the
				   // constraints indicate. No caching done
				   // yet. would need some more clever data
				   // structures for doing that.
      const std::pair<unsigned int, unsigned int>
	local_range = vec.block(block).local_range();

      my_indices.add_range (local_range.first, local_range.second);

      std::set<unsigned int> individual_indices;
      for (constraint_iterator it = begin_my_constraints;
	   it != end_my_constraints; ++it)
	for (unsigned int i=0; i<it->entries.size(); ++i)
	  if ((it->entries[i].first < local_range.first)
	      ||
	      (it->entries[i].first >= local_range.second))
	    individual_indices.insert (it->entries[i].first);

      my_indices.add_indices (individual_indices.begin(),
			      individual_indices.end());
    }

#ifdef DEAL_II_COMPILER_SUPPORTS_MPI
  const Epetra_MpiComm *mpi_comm
    = dynamic_cast<const Epetra_MpiComm*>(&vec.block(0).trilinos_vector().Comm());

  Assert (mpi_comm != 0, ExcInternalError());

  TrilinosWrappers::MPI::Vector vec_distribute
    (my_indices.make_trilinos_map (mpi_comm->Comm(), true));
#else
  TrilinosWrappers::MPI::Vector vec_distribute
    (my_indices.make_trilinos_map (MPI_COMM_WORLD, true));
#endif

				   // here we import the data
  vec_distribute.reinit(vec,true);

  for (unsigned int block=0; block<vec.n_blocks(); ++block)
    {
      typedef std::vector<ConstraintLine>::const_iterator constraint_iterator;
      ConstraintLine index_comparison;
      index_comparison.line = vec.block(block).local_range().first
	+vec.get_block_indices().block_start(block);
      const constraint_iterator begin_my_constraints =
	std::lower_bound (lines.begin(),lines.end(),index_comparison);

      index_comparison.line = vec.block(block).local_range().second
	+vec.get_block_indices().block_start(block);

      const constraint_iterator end_my_constraints
	= std::lower_bound(lines.begin(),lines.end(),index_comparison);

      for (constraint_iterator it = begin_my_constraints;
	   it != end_my_constraints; ++it)
	{
				       // fill entry in line
				       // next_constraint.line by adding the
				       // different contributions
	  double new_value = it->inhomogeneity;
	  for (unsigned int i=0; i<it->entries.size(); ++i)
	    new_value += (vec_distribute(it->entries[i].first) *
			  it->entries[i].second);
	  vec(it->line) = new_value;
	}
    }

				   // force every processor to write something
  unsigned int idx = vec.block(0).local_range().first;
  vec(idx) = vec(idx);
  vec.compress ();
}

#endif

#ifdef DEAL_II_USE_PETSC

				   // this is a specialization for a
				   // parallel (non-block) PETSc
				   // vector. The basic idea is to just work
				   // on the local range of the vector. But
				   // we need access to values that the
				   // local nodes are constrained to.

template<>
void
ConstraintMatrix::distribute (PETScWrappers::MPI::Vector &vec) const
{
  typedef std::vector<ConstraintLine>::const_iterator constraint_iterator;
  ConstraintLine index_comparison;
  index_comparison.line = vec.local_range().first;
  const constraint_iterator begin_my_constraints =
    std::lower_bound (lines.begin(),lines.end(),index_comparison);

  index_comparison.line = vec.local_range().second;
  const constraint_iterator end_my_constraints
    = std::lower_bound(lines.begin(),lines.end(),index_comparison);

				   // all indices we need to read from
  IndexSet my_indices (vec.size());

  const std::pair<unsigned int, unsigned int>
    local_range = vec.local_range();

  my_indices.add_range (local_range.first, local_range.second);

  std::set<unsigned int> individual_indices;
  for (constraint_iterator it = begin_my_constraints;
       it != end_my_constraints; ++it)
    for (unsigned int i=0; i<it->entries.size(); ++i)
      if ((it->entries[i].first < local_range.first)
	  ||
	  (it->entries[i].first >= local_range.second))
	individual_indices.insert (it->entries[i].first);

  my_indices.add_indices (individual_indices.begin(),
			  individual_indices.end());

  IndexSet local_range_is (vec.size());
  local_range_is.add_range(local_range.first, local_range.second);


				   // create a vector and import those indices
  PETScWrappers::MPI::Vector ghost_vec (vec.get_mpi_communicator(),
					local_range_is,
					my_indices);
  ghost_vec = vec;
  ghost_vec.update_ghost_values();

				   // finally do the distribution on own
				   // constraints
  for (constraint_iterator it = begin_my_constraints;
       it != end_my_constraints; ++it)
    {
				       // fill entry in line
				       // next_constraint.line by adding the
				       // different contributions
      PetscScalar new_value = it->inhomogeneity;
      for (unsigned int i=0; i<it->entries.size(); ++i)
	new_value += (PetscScalar(ghost_vec(it->entries[i].first)) *
                      it->entries[i].second);
      vec(it->line) = new_value;
    }

				   // force every processor to write something
  vec(local_range.first) = vec(local_range.first);

  vec.compress ();
}


template<>
void
ConstraintMatrix::distribute (PETScWrappers::MPI::BlockVector &/*vec*/) const
{
  AssertThrow (false, ExcNotImplemented());
}

#endif



unsigned int ConstraintMatrix::n_constraints () const
{
  return lines.size();
}



bool ConstraintMatrix::is_identity_constrained (const unsigned int index) const
{
  if (is_constrained(index) == false)
    return false;

  const ConstraintLine & p = lines[lines_cache[calculate_line_index(index)]];
  Assert (p.line == index, ExcInternalError());

				       // return if an entry for this
				       // line was found and if it has
				       // only one entry equal to 1.0
  return ((p.entries.size() == 1) &&
	  (p.entries[0].second == 1.0));
}



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
}



bool ConstraintMatrix::has_inhomogeneities () const
{
  for (std::vector<ConstraintLine>::const_iterator i=lines.begin();
       i!=lines.end(); ++i)
    if (i->inhomogeneity != 0.)
      return true;

  return false;
}


void ConstraintMatrix::print (std::ostream &out) const
{
  for (unsigned int i=0; i!=lines.size(); ++i)
    {
				       // output the list of
				       // constraints as pairs of dofs
				       // and their weights
      if (lines[i].entries.size() > 0)
	{
	  for (unsigned int j=0; j<lines[i].entries.size(); ++j)
	    out << "    " << lines[i].line
		<< " " << lines[i].entries[j].first
		<< ":  " << lines[i].entries[j].second << "\n";

				       // print out inhomogeneity.
	  if (lines[i].inhomogeneity != 0)
	    out << "    " << lines[i].line
		<< ": " << lines[i].inhomogeneity << "\n";
	}
      else
					 // but also output something
					 // if the constraint simply
					 // reads x[13]=0, i.e. where
					 // the right hand side is not
					 // a linear combination of
					 // other dofs
	{
	  if (lines[i].inhomogeneity != 0)
	    out << "    " << lines[i].line
		<< " = " << lines[i].inhomogeneity
		<< "\n";
	  else
	    out << "    " << lines[i].line << " = 0\n";
	}
    }

  AssertThrow (out, ExcIO());
}



void
ConstraintMatrix::write_dot (std::ostream &out) const
{
  out << "digraph constraints {"
      << std::endl;
  for (unsigned int i=0; i!=lines.size(); ++i)
    {
				       // same concept as in the
				       // previous function
      if (lines[i].entries.size() > 0)
	for (unsigned int j=0; j<lines[i].entries.size(); ++j)
	  out << "  " << lines[i].line << "->" << lines[i].entries[j].first
	      << "; // weight: "
	      << lines[i].entries[j].second
	      << "\n";
      else
	out << "  " << lines[i].line << "\n";
    }
  out << "}" << std::endl;
}



unsigned int
ConstraintMatrix::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (lines) +
	  MemoryConsumption::memory_consumption (lines_cache) +
	  MemoryConsumption::memory_consumption (sorted) +
	  MemoryConsumption::memory_consumption (local_lines));
}





// explicit instantiations
//
// define a list of functions for vectors and matrices, respectively, where
// the vector/matrix can be replaced using a preprocessor variable
// VectorType/MatrixType. note that we need a space between "VectorType" and
// ">" to disambiguate ">>" when VectorType trails in an angle bracket

// TODO: The way we define all the instantiations is probably not the very
// best one. Try to find a better description.

#define VECTOR_FUNCTIONS(VectorType) \
  template void ConstraintMatrix::condense<VectorType >(const VectorType &uncondensed,\
					               VectorType       &condensed) const;\
  template void ConstraintMatrix::condense<VectorType >(VectorType &vec) const;\
  template void ConstraintMatrix::condense<float,VectorType >(const SparseMatrix<float> &uncondensed, \
							      const VectorType &uncondensed_vector, \
							      SparseMatrix<float> &condensed, \
							      VectorType       &condensed_vector) const; \
  template void ConstraintMatrix::condense<double,VectorType >(const SparseMatrix<double> &uncondensed, \
							       const VectorType &uncondensed_vector, \
							       SparseMatrix<double> &condensed, \
							       VectorType       &condensed_vector) const; \
  template void ConstraintMatrix::set_zero<VectorType >(VectorType &vec) const;\
  template void ConstraintMatrix:: \
    distribute_local_to_global<VectorType > (const Vector<double>            &, \
                                             const std::vector<unsigned int> &, \
                                             VectorType                      &, \
					     const FullMatrix<double>        &) const; \
  template void ConstraintMatrix::distribute<VectorType >(const VectorType &condensed,\
					                 VectorType       &uncondensed) const;\
  template void ConstraintMatrix::distribute<VectorType >(VectorType &vec) const

#define PARALLEL_VECTOR_FUNCTIONS(VectorType) \
  template void ConstraintMatrix:: \
    distribute_local_to_global<VectorType > (const Vector<double>            &, \
                                             const std::vector<unsigned int> &, \
                                             VectorType                      &, \
					     const FullMatrix<double>        &) const



VECTOR_FUNCTIONS(Vector<float>);
VECTOR_FUNCTIONS(Vector<double>);
VECTOR_FUNCTIONS(BlockVector<double>);
VECTOR_FUNCTIONS(BlockVector<float>);

// TODO: Can PETSc really do all the operations required by the above
// condense/distribute function etc also on distributed vectors? Trilinos
// can't do that - we have to rewrite those functions by hand if we want to
// use them. The key is to use local ranges etc., which still needs to be
// implemented.
#ifdef DEAL_II_USE_PETSC
VECTOR_FUNCTIONS(PETScWrappers::Vector);
VECTOR_FUNCTIONS(PETScWrappers::BlockVector);
VECTOR_FUNCTIONS(PETScWrappers::MPI::Vector);
VECTOR_FUNCTIONS(PETScWrappers::MPI::BlockVector);
#endif

#ifdef DEAL_II_USE_TRILINOS
VECTOR_FUNCTIONS(TrilinosWrappers::Vector);
VECTOR_FUNCTIONS(TrilinosWrappers::BlockVector);
PARALLEL_VECTOR_FUNCTIONS(TrilinosWrappers::MPI::Vector);
PARALLEL_VECTOR_FUNCTIONS(TrilinosWrappers::MPI::BlockVector);
#endif

#define CONDENSE_FUNCTIONS(VectorType, number, MatrixType)		\
  template void ConstraintMatrix::condense<number,VectorType >(MatrixType &uncondensed, \
							       VectorType &vec) const \


CONDENSE_FUNCTIONS(Vector<float>,float,SparseMatrix<float>);
CONDENSE_FUNCTIONS(Vector<double>,float,SparseMatrix<float>);
CONDENSE_FUNCTIONS(Vector<double>,double,SparseMatrix<double>);
CONDENSE_FUNCTIONS(Vector<float>,double,SparseMatrix<double>);
CONDENSE_FUNCTIONS(BlockVector<double>,float,BlockSparseMatrix<float>);
CONDENSE_FUNCTIONS(BlockVector<float>,float,BlockSparseMatrix<float>);
CONDENSE_FUNCTIONS(BlockVector<double>,double,BlockSparseMatrix<double>);
CONDENSE_FUNCTIONS(BlockVector<float>,double,BlockSparseMatrix<double>);


template
void
ConstraintMatrix::condense<float>(const SparseMatrix<float> &uncondensed,
                                  SparseMatrix<float> &condensed) const;
template
void
ConstraintMatrix::condense<float>(SparseMatrix<float> &uncondensed) const;

template
void
ConstraintMatrix::condense<double>(const SparseMatrix<double> &uncondensed,
                                   SparseMatrix<double> &condensed) const;

template
void
ConstraintMatrix::condense<double>(SparseMatrix<double> &uncondensed) const;


// block sparse matrices are only implemented for one of the two matrix
// functions (the single-argument, in-place function)
template
void
ConstraintMatrix::condense<double>(BlockSparseMatrix<double> &uncondensed) const;

template
void
ConstraintMatrix::condense<float>(BlockSparseMatrix<float> &uncondensed) const;


#define MATRIX_VECTOR_FUNCTIONS(MatrixType, VectorType) \
template void ConstraintMatrix:: \
distribute_local_to_global<MatrixType,VectorType > (const FullMatrix<double>        &, \
                                                    const Vector<double>            &, \
                                                    const std::vector<unsigned int> &, \
                                                    MatrixType                      &, \
                                                    VectorType                      &, \
                                                    internal::bool2type<false>) const
#define MATRIX_FUNCTIONS(MatrixType) \
template void ConstraintMatrix:: \
distribute_local_to_global<MatrixType,Vector<double> > (const FullMatrix<double>        &, \
                                                        const Vector<double>            &, \
                                                        const std::vector<unsigned int> &, \
                                                        MatrixType                      &, \
                                                        Vector<double>                  &, \
                                                        internal::bool2type<false>) const
#define BLOCK_MATRIX_VECTOR_FUNCTIONS(MatrixType, VectorType)   \
template void ConstraintMatrix:: \
distribute_local_to_global<MatrixType,VectorType > (const FullMatrix<double>        &, \
                                                    const Vector<double>            &, \
                                                    const std::vector<unsigned int> &, \
                                                    MatrixType                      &, \
                                                    VectorType                      &, \
                                                    internal::bool2type<true>) const
#define BLOCK_MATRIX_FUNCTIONS(MatrixType)      \
template void ConstraintMatrix:: \
distribute_local_to_global<MatrixType,Vector<double> > (const FullMatrix<double>        &, \
                                                        const Vector<double>            &, \
                                                        const std::vector<unsigned int> &, \
                                                        MatrixType                      &, \
                                                        Vector<double>                  &, \
                                                        internal::bool2type<true>) const

MATRIX_FUNCTIONS(SparseMatrix<double>);
MATRIX_FUNCTIONS(SparseMatrix<float>);
MATRIX_VECTOR_FUNCTIONS(SparseMatrix<float>, Vector<float>);

BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrix<double>);
BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrix<float>);
BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrix<double>, BlockVector<double>);
BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrix<float>,  BlockVector<float>);
BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrix<float>,  BlockVector<double>);

MATRIX_FUNCTIONS(SparseMatrixEZ<double>);
MATRIX_FUNCTIONS(SparseMatrixEZ<float>);
MATRIX_VECTOR_FUNCTIONS(SparseMatrixEZ<float>,  Vector<float>);

// BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrixEZ<double>);
// BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrixEZ<float>,  Vector<float>);

#ifdef DEAL_II_USE_PETSC
MATRIX_FUNCTIONS(PETScWrappers::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(PETScWrappers::BlockSparseMatrix);
MATRIX_FUNCTIONS(PETScWrappers::MPI::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(PETScWrappers::MPI::BlockSparseMatrix);
MATRIX_VECTOR_FUNCTIONS(PETScWrappers::SparseMatrix, PETScWrappers::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(PETScWrappers::BlockSparseMatrix, PETScWrappers::BlockVector);
MATRIX_VECTOR_FUNCTIONS(PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(PETScWrappers::MPI::BlockSparseMatrix ,PETScWrappers::MPI::BlockVector);
#endif

#ifdef DEAL_II_USE_TRILINOS
MATRIX_FUNCTIONS(TrilinosWrappers::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix);
MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::SparseMatrix, TrilinosWrappers::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix, TrilinosWrappers::BlockVector);
MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::SparseMatrix, TrilinosWrappers::MPI::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix, TrilinosWrappers::MPI::BlockVector);
#endif


#define SPARSITY_FUNCTIONS(SparsityType) \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityType> (\
    const std::vector<unsigned int> &, \
    SparsityType &,                    \
    const bool,                        \
    const Table<2,bool> &, \
    internal::bool2type<false>) const; \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityType> (\
    const std::vector<unsigned int> &, \
    const std::vector<unsigned int> &, \
    SparsityType &,                    \
    const bool,                        \
    const Table<2,bool> &) const
#define BLOCK_SPARSITY_FUNCTIONS(SparsityType) \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityType> (\
    const std::vector<unsigned int> &, \
    SparsityType &,                    \
    const bool,                        \
    const Table<2,bool> &, \
    internal::bool2type<true>) const; \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityType> (\
    const std::vector<unsigned int> &, \
    const std::vector<unsigned int> &, \
    SparsityType &,                    \
    const bool,                        \
    const Table<2,bool> &) const

SPARSITY_FUNCTIONS(SparsityPattern);
SPARSITY_FUNCTIONS(CompressedSparsityPattern);
SPARSITY_FUNCTIONS(CompressedSetSparsityPattern);
SPARSITY_FUNCTIONS(CompressedSimpleSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockCompressedSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockCompressedSetSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockCompressedSimpleSparsityPattern);

#ifdef DEAL_II_USE_TRILINOS
SPARSITY_FUNCTIONS(TrilinosWrappers::SparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(TrilinosWrappers::BlockSparsityPattern);
#endif


#define ONLY_MATRIX_FUNCTIONS(MatrixType) \
  template void ConstraintMatrix::distribute_local_to_global<MatrixType > (\
  const FullMatrix<double>        &, \
  const std::vector<unsigned int> &, \
  const std::vector<unsigned int> &, \
  MatrixType                      &) const

ONLY_MATRIX_FUNCTIONS(SparseMatrix<float>);
ONLY_MATRIX_FUNCTIONS(SparseMatrix<double>);
ONLY_MATRIX_FUNCTIONS(MatrixBlock<SparseMatrix<float> >);
ONLY_MATRIX_FUNCTIONS(MatrixBlock<SparseMatrix<double> >);

DEAL_II_NAMESPACE_CLOSE
