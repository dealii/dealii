// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/constraint_matrix.templates.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/parallel_vector.h>
#include <deal.II/lac/parallel_block_vector.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_vector.h>
#include <deal.II/lac/petsc_parallel_block_vector.h>
#include <deal.II/lac/petsc_parallel_sparse_matrix.h>
#include <deal.II/lac/petsc_parallel_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/trilinos_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/matrix_block.h>

#include <algorithm>
#include <numeric>
#include <set>
#include <ostream>

DEAL_II_NAMESPACE_OPEN



// Static member variable
const Table<2,bool> ConstraintMatrix::default_empty_table = Table<2,bool>();



bool
ConstraintMatrix::check_zero_weight (const std::pair<size_type, double> &p)
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



std::size_t
ConstraintMatrix::ConstraintLine::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (line) +
          MemoryConsumption::memory_consumption (entries) +
          MemoryConsumption::memory_consumption (inhomogeneity));
}



void
ConstraintMatrix::add_lines (const std::set<size_type> &lines)
{
  for (std::set<size_type>::const_iterator
       i = lines.begin(); i != lines.end(); ++i)
    add_line (*i);
}



void
ConstraintMatrix::add_lines (const std::vector<bool> &lines)
{
  for (size_type i=0; i<lines.size(); ++i)
    if (lines[i] == true)
      add_line (i);
}



void
ConstraintMatrix::add_lines (const IndexSet &lines)
{
  for (size_type i=0; i<lines.n_elements(); ++i)
    add_line (lines.nth_index_in_set(i));
}



void
ConstraintMatrix::add_entries
(const size_type                                  line,
 const std::vector<std::pair<size_type,double> > &col_val_pairs)
{
  Assert (sorted==false, ExcMatrixIsClosed());
  Assert (is_constrained(line), ExcLineInexistant(line));

  ConstraintLine *line_ptr = &lines[lines_cache[calculate_line_index(line)]];
  Assert (line_ptr->line == line, ExcInternalError());

  // if in debug mode, check whether an entry for this column already
  // exists and if its the same as the one entered at present
  //
  // in any case: skip this entry if an entry for this column already
  // exists, since we don't want to enter it twice
  for (std::vector<std::pair<size_type,double> >::const_iterator
       col_val_pair = col_val_pairs.begin();
       col_val_pair!=col_val_pairs.end(); ++col_val_pair)
    {
      Assert (line != col_val_pair->first,
              ExcMessage ("Can't constrain a degree of freedom to itself"));

      for (ConstraintLine::Entries::const_iterator
           p=line_ptr->entries.begin();
           p != line_ptr->entries.end(); ++p)
        if (p->first == col_val_pair->first)
          {
            // entry exists, break innermost loop
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
  if (constraints.n_constraints() == 0)
    return;

  Assert (filter.size() > constraints.lines.back().line,
          ExcMessage ("Filter needs to be larger than constraint matrix size."));
  for (std::vector<ConstraintLine>::const_iterator line=constraints.lines.begin();
       line!=constraints.lines.end(); ++line)
    if (filter.is_element(line->line))
      {
        const size_type row = filter.index_within_set (line->line);
        add_line (row);
        set_inhomogeneity (row, line->inhomogeneity);
        for (size_type i=0; i<line->entries.size(); ++i)
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

  // update list of pointers and give the vector a sharp size since we
  // won't modify the size any more after this point.
  {
    std::vector<size_type> new_lines (lines_cache.size(),
                                      numbers::invalid_size_type);
    size_type counter = 0;
    for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
         line!=lines.end(); ++line, ++counter)
      new_lines[calculate_line_index(line->line)] = counter;
    std::swap (lines_cache, new_lines);
  }

  // in debug mode: check whether we really set the pointers correctly.
  for (size_type i=0; i<lines_cache.size(); ++i)
    if (lines_cache[i] != numbers::invalid_size_type)
      Assert (i == calculate_line_index(lines[lines_cache[i]].line),
              ExcInternalError());

  // first, strip zero entries, as we have to do that only once
  for (std::vector<ConstraintLine>::iterator line = lines.begin();
       line!=lines.end(); ++line)
    // first remove zero entries. that would mean that in the linear
    // constraint for a node, x_i = ax_1 + bx_2 + ..., another node times 0
    // appears. obviously, 0*something can be omitted
    line->entries.erase (std::remove_if (line->entries.begin(),
                                         line->entries.end(),
                                         &check_zero_weight),
                         line->entries.end());



#ifdef DEBUG
  // In debug mode we are computing an estimate for the maximum number
  // of constraints so that we can bail out if there is a cycle in the
  // constraints (which is easier than searching for cycles in the graph).
  //
  // Let us figure out the largest dof index. This is an upper bound for the
  // number of constraints because it is an approximation for the number of dofs
  // in our system.
  size_type largest_idx = 0;
  for (std::vector<ConstraintLine>::iterator line = lines.begin();
       line!=lines.end(); ++line)
    {
      for (ConstraintLine::Entries::iterator it = line->entries.begin(); it!=line->entries.end(); ++it)
        {
          largest_idx=std::max(largest_idx, it->first);
        }
    }
#endif

  // replace references to dofs that are themselves constrained. note that
  // because we may replace references to other dofs that may themselves be
  // constrained to third ones, we have to iterate over all this until we
  // replace no chains of constraints any more
  //
  // the iteration replaces references to constrained degrees of freedom by
  // second-order references. for example if x3=x0/2+x2/2 and x2=x0/2+x1/2,
  // then the new list will be x3=x0/2+x0/4+x1/4. note that x0 appear
  // twice. we will throw this duplicate out in the following step, where
  // we sort the list so that throwing out duplicates becomes much more
  // efficient. also, we have to do it only once, rather than in each
  // iteration
  size_type iteration = 0;
  while (true)
    {
      bool chained_constraint_replaced = false;

      for (std::vector<ConstraintLine>::iterator line = lines.begin();
           line!=lines.end(); ++line)
        {
#ifdef DEBUG
          // we need to keep track of how many replacements we do in this line, because we can
          // end up in a cycle A->B->C->A without the number of entries growing.
          size_type n_replacements = 0;
#endif

          // loop over all entries of this line (including ones that we
          // have appended in this go around) and see whether they are
          // further constrained. ignore elements that we don't store on
          // the current processor
          size_type entry = 0;
          while (entry < line->entries.size())
            if (((local_lines.size() == 0)
                 ||
                 (local_lines.is_element(line->entries[entry].first)))
                &&
                is_constrained (line->entries[entry].first))
              {
                // ok, this entry is further constrained:
                chained_constraint_replaced = true;

                // look up the chain of constraints for this entry
                const size_type  dof_index = line->entries[entry].first;
                const double     weight = line->entries[entry].second;

                Assert (dof_index != line->line,
                        ExcMessage ("Cycle in constraints detected!"));

                const ConstraintLine *constrained_line =
                  &lines[lines_cache[calculate_line_index(dof_index)]];
                Assert (constrained_line->line == dof_index,
                        ExcInternalError());

                // now we have to replace an entry by its expansion. we do
                // that by overwriting the entry by the first entry of the
                // expansion and adding the remaining ones to the end,
                // where we will later process them once more
                //
                // we can of course only do that if the DoF that we are
                // currently handle is constrained by a linear combination
                // of other dofs:
                if (constrained_line->entries.size() > 0)
                  {
                    for (size_type i=0; i<constrained_line->entries.size(); ++i)
                      Assert (dof_index != constrained_line->entries[i].first,
                              ExcMessage ("Cycle in constraints detected!"));

                    // replace first entry, then tack the rest to the end
                    // of the list
                    line->entries[entry] =
                      std::make_pair (constrained_line->entries[0].first,
                                      constrained_line->entries[0].second *
                                      weight);

                    for (size_type i=1; i<constrained_line->entries.size(); ++i)
                      line->entries
                      .push_back (std::make_pair (constrained_line->entries[i].first,
                                                  constrained_line->entries[i].second *
                                                  weight));

#ifdef DEBUG
                    // keep track of how many entries we replace in this
                    // line. If we do more than there are constraints or
                    // dofs in our system, we must have a cycle.
                    ++n_replacements;
                    Assert(n_replacements/2<largest_idx, ExcMessage("Cycle in constraints detected!"));
                    if (n_replacements/2>=largest_idx)
                      return; // this enables us to test for this Exception.
#endif
                  }
                else
                  // the DoF that we encountered is not constrained by a
                  // linear combination of other dofs but is equal to just
                  // the inhomogeneity (i.e. its chain of entries is
                  // empty). in that case, we can't just overwrite the
                  // current entry, but we have to actually eliminate it
                  {
                    line->entries.erase (line->entries.begin()+entry);
                  }

                line->inhomogeneity += constrained_line->inhomogeneity *
                                       weight;

                // now that we're here, do not increase index by one but
                // rather make another pass for the present entry because
                // we have replaced the present entry by another one, or
                // because we have deleted it and shifted all following
                // ones one forward
              }
            else
              // entry not further constrained. just move ahead by one
              ++entry;
        }

      // if we didn't do anything in this round, then quit the loop
      if (chained_constraint_replaced == false)
        break;

      // increase iteration count. note that we should not iterate more
      // times than there are constraints, since this puts a natural upper
      // bound on the length of constraint chains
      ++iteration;
      Assert (iteration <= lines.size(), ExcInternalError());
    }

  // finally sort the entries and re-scale them if necessary. in this step,
  // we also throw out duplicates as mentioned above. moreover, as some
  // entries might have had zero weights, we replace them by a vector with
  // sharp sizes.
  for (std::vector<ConstraintLine>::iterator line = lines.begin();
       line!=lines.end(); ++line)
    {
      std::sort (line->entries.begin(), line->entries.end());

      // loop over the now sorted list and see whether any of the entries
      // references the same dofs more than once in order to find how many
      // non-duplicate entries we have. This lets us allocate the correct
      // amount of memory for the constraint entries.
      size_type duplicates = 0;
      for (size_type i=1; i<line->entries.size(); ++i)
        if (line->entries[i].first == line->entries[i-1].first)
          duplicates++;

      if (duplicates > 0 || line->entries.size() < line->entries.capacity())
        {
          ConstraintLine::Entries new_entries;

          // if we have no duplicates, copy verbatim the entries. this way,
          // the final size is of the vector is correct.
          if (duplicates == 0)
            new_entries = line->entries;
          else
            {
              // otherwise, we need to go through the list by and and
              // resolve the duplicates
              new_entries.reserve (line->entries.size() - duplicates);
              new_entries.push_back(line->entries[0]);
              for (size_type j=1; j<line->entries.size(); ++j)
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

              // make sure there are really no duplicates left and that the
              // list is still sorted
              for (size_type j=1; j<new_entries.size(); ++j)
                {
                  Assert (new_entries[j].first != new_entries[j-1].first,
                          ExcInternalError());
                  Assert (new_entries[j].first > new_entries[j-1].first,
                          ExcInternalError());
                }
            }

          // replace old list of constraints for this dof by the new one
          line->entries.swap (new_entries);
        }

      // finally do the following check: if the sum of weights for the
      // constraints is close to one, but not exactly one, then rescale all
      // the weights so that they sum up to 1. this adds a little numerical
      // stability and avoids all sorts of problems where the actual value
      // is close to, but not quite what we expected
      //
      // the case where the weights don't quite sum up happens when we
      // compute the interpolation weights "on the fly", i.e. not from
      // precomputed tables. in this case, the interpolation weights are
      // also subject to round-off
      double sum = 0;
      for (size_type i=0; i<line->entries.size(); ++i)
        sum += line->entries[i].second;
      if ((sum != 1.0) && (std::fabs (sum-1.) < 1.e-13))
        {
          for (size_type i=0; i<line->entries.size(); ++i)
            line->entries[i].second /= sum;
          line->inhomogeneity /= sum;
        }
    } // end of loop over all constraint lines

#ifdef DEBUG
  // if in debug mode: check that no dof is constrained to another dof that
  // is also constrained. exclude dofs from this check whose constraint
  // lines are not stored on the local processor
  for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line)
    for (ConstraintLine::Entries::const_iterator
         entry=line->entries.begin();
         entry!=line->entries.end(); ++entry)
      if ((local_lines.size() == 0)
          ||
          (local_lines.is_element(entry->first)))
        {
          // make sure that entry->first is not the index of a line itself
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
  AssertThrow(local_lines == other_constraints.local_lines,
              ExcNotImplemented());

  // store the previous state with respect to sorting
  const bool object_was_sorted = sorted;
  sorted = false;

  if (other_constraints.lines_cache.size() > lines_cache.size())
    lines_cache.resize(other_constraints.lines_cache.size(),
                       numbers::invalid_size_type);

  // first action is to fold into the present object possible constraints
  // in the second object. we don't strictly need to do this any more since
  // the ConstraintMatrix has learned to deal with chains of constraints in
  // the close() function, but we have traditionally done this and it's not
  // overly hard to do.
  //
  // for this, loop over all constraints and replace the constraint lines
  // with a new one where constraints are replaced if necessary.
  ConstraintLine::Entries tmp;
  for (std::vector<ConstraintLine>::iterator line=lines.begin();
       line!=lines.end(); ++line)
    {
      tmp.clear ();
      for (size_type i=0; i<line->entries.size(); ++i)
        {
          // if the present dof is not constrained, or if we won't take the
          // constraint from the other object, then simply copy it over
          if (other_constraints.is_constrained(line->entries[i].first) == false
              ||
              ((merge_conflict_behavior != right_object_wins)
               &&
               other_constraints.is_constrained(line->entries[i].first)
               &&
               this->is_constrained(line->entries[i].first)))
            tmp.push_back(line->entries[i]);
          else
            // otherwise resolve further constraints by replacing the old
            // entry by a sequence of new entries taken from the other
            // object, but with multiplied weights
            {
              const ConstraintLine::Entries *other_line
                = other_constraints.get_constraint_entries (line->entries[i].first);
              Assert (other_line != 0,
                      ExcInternalError());

              const double weight = line->entries[i].second;

              for (ConstraintLine::Entries::const_iterator j=other_line->begin();
                   j!=other_line->end(); ++j)
                tmp.push_back (std::pair<size_type,double>(j->first,
                                                           j->second*weight));

              line->inhomogeneity += other_constraints.get_inhomogeneity(line->entries[i].first) *
                                     weight;
            }
        }
      // finally exchange old and newly resolved line
      line->entries.swap (tmp);
    }



  // next action: append those lines at the end that we want to add
  for (std::vector<ConstraintLine>::const_iterator
       line=other_constraints.lines.begin();
       line!=other_constraints.lines.end(); ++line)
    if (is_constrained(line->line) == false)
      lines.push_back (*line);
    else
      {
        // the constrained dof we want to copy from the other object is
        // also constrained here. let's see what we should do with that
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
            // we need to replace the existing constraint by the one from
            // the other object
            lines[lines_cache[calculate_line_index(line->line)]].entries
              = line->entries;
            lines[lines_cache[calculate_line_index(line->line)]].inhomogeneity
              = line->inhomogeneity;
            break;

          default:
            Assert (false, ExcNotImplemented());
          }
      }

  // update the lines cache
  size_type counter = 0;
  for (std::vector<ConstraintLine>::const_iterator line=lines.begin();
       line!=lines.end(); ++line, ++counter)
    lines_cache[calculate_line_index(line->line)] = counter;

  // if the object was sorted before, then make sure it is so afterward as
  // well. otherwise leave everything in the unsorted state
  if (object_was_sorted == true)
    close ();
}



void ConstraintMatrix::shift (const size_type offset)
{
  //TODO: this doesn't work with IndexSets yet. [TH]
  AssertThrow(local_lines.size()==0, ExcNotImplemented());

  lines_cache.insert (lines_cache.begin(), offset,
                      numbers::invalid_size_type);

  for (std::vector<ConstraintLine>::iterator i = lines.begin();
       i != lines.end(); ++i)
    {
      i->line += offset;
      for (ConstraintLine::Entries::iterator
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
    std::vector<size_type> tmp;
    lines_cache.swap (tmp);
  }

  sorted = false;
}



void ConstraintMatrix::reinit (const IndexSet &local_constraints)
{
  local_lines = local_constraints;

  // make sure the IndexSet is compressed. Otherwise this can lead to crashes
  // that are hard to find (only happen in release mode).
  // see tests/mpi/constraint_matrix_crash_01
  local_lines.compress();

  clear();
}



void ConstraintMatrix::condense (SparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraint matrix
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_size_type)
        {
          // regular line. loop over cols all valid cols. note that this
          // changes the line we are presently working on: we add additional
          // entries. these are put to the end of the row. however, as
          // constrained nodes cannot be constrained to other constrained
          // nodes, nothing will happen if we run into these added nodes, as
          // they can't be distributed further. we might store the position of
          // the last old entry and stop work there, but since operating on
          // the newly added ones only takes two comparisons (column index
          // valid, distribute[column] necessarily
          // ==numbers::invalid_size_type), it is cheaper to not do so and
          // run right until the end of the line
          for (SparsityPattern::iterator entry = sparsity.begin(row);
               ((entry != sparsity.end(row)) &&
                entry->is_valid_entry());
               ++entry)
            {
              const size_type column = entry->column();

              if (distribute[column] != numbers::invalid_size_type)
                {
                  // distribute entry at regular row @p{row} and irregular
                  // column sparsity.colnums[j]
                  for (size_type q=0;
                       q!=lines[distribute[column]].entries.size();
                       ++q)
                    sparsity.add (row,
                                  lines[distribute[column]].entries[q].first);
                }
            }
        }
      else
        // row must be distributed. note that here the present row is not
        // touched (unlike above)
        {
          for (SparsityPattern::iterator entry = sparsity.begin(row);
               (entry != sparsity.end(row)) && entry->is_valid_entry(); ++entry)
            {
              const size_type column = entry->column();
              if (distribute[column] == numbers::invalid_size_type)
                // distribute entry at irregular row @p{row} and regular
                // column sparsity.colnums[j]
                for (size_type q=0;
                     q!=lines[distribute[row]].entries.size(); ++q)
                  sparsity.add (lines[distribute[row]].entries[q].first,
                                column);
              else
                // distribute entry at irregular row @p{row} and irregular
                // column sparsity.get_column_numbers()[j]
                for (size_type p=0; p!=lines[distribute[row]].entries.size(); ++p)
                  for (size_type q=0;
                       q!=lines[distribute[column]].entries.size(); ++q)
                    sparsity.add (lines[distribute[row]].entries[p].first,
                                  lines[distribute[column]].entries[q].first);
            }
        }
    }

  sparsity.compress();
}




void ConstraintMatrix::condense (DynamicSparsityPattern &sparsity) const
{
  Assert (sorted == true, ExcMatrixNotClosed());
  Assert (sparsity.n_rows() == sparsity.n_cols(),
          ExcNotQuadratic());

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraint matrix
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over cols. note that as we proceed to
        // distribute cols, the loop may get longer
        for (size_type j=0; j<sparsity.row_length(row); ++j)
          {
            const size_type column = sparsity.column_number(row,j);

            if (distribute[column] != numbers::invalid_size_type)
              {
                // distribute entry at regular row @p{row} and irregular
                // column column. note that this changes the line we are
                // presently working on: we add additional entries. if we
                // add another entry at a column behind the present one, we
                // will encounter it later on (but since it can't be
                // further constrained, won't have to do anything about
                // it). if we add it up front of the present column, we
                // will find the present column later on again as it was
                // shifted back (again nothing happens, in particular no
                // endless loop, as when we encounter it the second time we
                // won't be able to add more entries as they all already
                // exist, but we do the same work more often than
                // necessary, and the loop gets longer), so move the cursor
                // one to the right in the case that we add an entry up
                // front that did not exist before. check whether it
                // existed before by tracking the length of this row
                size_type old_rowlength = sparsity.row_length(row);
                for (size_type q=0;
                     q!=lines[distribute[column]].entries.size();
                     ++q)
                  {
                    const size_type
                    new_col = lines[distribute[column]].entries[q].first;

                    sparsity.add (row, new_col);

                    const size_type new_rowlength = sparsity.row_length(row);
                    if ((new_col < column) && (old_rowlength != new_rowlength))
                      ++j;
                    old_rowlength = new_rowlength;
                  };
              };
          }
      else
        // row must be distributed
        for (size_type j=0; j<sparsity.row_length(row); ++j)
          {
            const size_type column = sparsity.column_number(row,j);

            if (distribute[column] == numbers::invalid_size_type)
              // distribute entry at irregular row @p{row} and regular
              // column sparsity.colnums[j]
              for (size_type q=0;
                   q!=lines[distribute[row]].entries.size(); ++q)
                sparsity.add (lines[distribute[row]].entries[q].first,
                              column);
            else
              // distribute entry at irregular row @p{row} and irregular
              // column sparsity.get_column_numbers()[j]
              for (size_type p=0; p!=lines[distribute[row]].entries.size(); ++p)
                for (size_type q=0;
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

  const size_type n_blocks = sparsity.n_block_rows();

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraint matrix
  // handles this index
  std::vector<size_type> distribute (sparsity.n_rows(),
                                     numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      // get index of this row within the blocks
      const std::pair<size_type,size_type>
      block_index = index_mapping.global_to_local(row);
      const size_type block_row = block_index.first;

      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over all columns and see whether this column
        // must be distributed
        {

          // to loop over all entries in this row, we have to loop over all
          // blocks in this blockrow and the corresponding row therein
          for (size_type block_col=0; block_col<n_blocks; ++block_col)
            {
              const SparsityPattern &
              block_sparsity = sparsity.block(block_row, block_col);

              for (SparsityPattern::const_iterator
                   entry = block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                   entry->is_valid_entry();
                   ++entry)
                {
                  const size_type global_col
                    = index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at regular row @p{row} and
                    // irregular column global_col
                    {
                      for (size_type q=0;
                           q!=lines[distribute[global_col]].entries.size(); ++q)
                        sparsity.add (row,
                                      lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
      else
        {
          // row must be distributed. split the whole row into the chunks
          // defined by the blocks
          for (size_type block_col=0; block_col<n_blocks; ++block_col)
            {
              const SparsityPattern &
              block_sparsity = sparsity.block(block_row,block_col);

              for (SparsityPattern::const_iterator
                   entry = block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                   entry->is_valid_entry();
                   ++entry)
                {
                  const size_type global_col
                    = index_mapping.local_to_global (block_col, entry->column());

                  if (distribute[global_col] == numbers::invalid_size_type)
                    // distribute entry at irregular row @p{row} and
                    // regular column global_col.
                    {
                      for (size_type q=0; q!=lines[distribute[row]].entries.size(); ++q)
                        sparsity.add (lines[distribute[row]].entries[q].first, global_col);
                    }
                  else
                    // distribute entry at irregular row @p{row} and
                    // irregular column @p{global_col}
                    {
                      for (size_type p=0; p!=lines[distribute[row]].entries.size(); ++p)
                        for (size_type q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
                          sparsity.add (lines[distribute[row]].entries[p].first,
                                        lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
    }

  sparsity.compress();
}




void ConstraintMatrix::condense (BlockDynamicSparsityPattern &sparsity) const
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

  const size_type n_blocks = sparsity.n_block_rows();

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraint matrix
  // handles this index
  std::vector<size_type> distribute (sparsity.n_rows(),
                                     numbers::invalid_size_type);

  for (size_type c=0; c<lines.size(); ++c)
    distribute[lines[c].line] = static_cast<signed int>(c);

  const size_type n_rows = sparsity.n_rows();
  for (size_type row=0; row<n_rows; ++row)
    {
      // get index of this row within the blocks
      const std::pair<size_type,size_type>
      block_index = index_mapping.global_to_local(row);
      const size_type block_row = block_index.first;
      const size_type local_row = block_index.second;

      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over all columns and see whether this column
        // must be distributed. note that as we proceed to distribute cols,
        // the loop over cols may get longer.
        //
        // don't try to be clever here as in the algorithm for the
        // DynamicSparsityPattern, as that would be much more
        // complicated here. after all, we know that compressed patterns
        // are inefficient...
        {

          // to loop over all entries in this row, we have to loop over all
          // blocks in this blockrow and the corresponding row therein
          for (size_type block_col=0; block_col<n_blocks; ++block_col)
            {
              const DynamicSparsityPattern &
              block_sparsity = sparsity.block(block_row, block_col);

              for (size_type j=0; j<block_sparsity.row_length(local_row); ++j)
                {
                  const size_type global_col
                    = index_mapping.local_to_global(block_col,
                                                    block_sparsity.column_number(local_row,j));

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at regular row @p{row} and
                    // irregular column global_col
                    {
                      for (size_type q=0;
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
          // row must be distributed. split the whole row into the chunks
          // defined by the blocks
          for (size_type block_col=0; block_col<n_blocks; ++block_col)
            {
              const DynamicSparsityPattern &
              block_sparsity = sparsity.block(block_row,block_col);

              for (size_type j=0; j<block_sparsity.row_length(local_row); ++j)
                {
                  const size_type global_col
                    = index_mapping.local_to_global (block_col,
                                                     block_sparsity.column_number(local_row,j));

                  if (distribute[global_col] == numbers::invalid_size_type)
                    // distribute entry at irregular row @p{row} and
                    // regular column global_col.
                    {
                      for (size_type q=0;
                           q!=lines[distribute[row]].entries.size(); ++q)
                        sparsity.add (lines[distribute[row]].entries[q].first,
                                      global_col);
                    }
                  else
                    // distribute entry at irregular row @p{row} and
                    // irregular column @p{global_col}
                    {
                      for (size_type p=0;
                           p!=lines[distribute[row]].entries.size(); ++p)
                        for (size_type q=0; q!=lines[distribute[global_col]].entries.size(); ++q)
                          sparsity.add (lines[distribute[row]].entries[p].first,
                                        lines[distribute[global_col]].entries[q].first);
                    };
                };
            };
        };
    };
}



bool ConstraintMatrix::is_identity_constrained (const size_type index) const
{
  if (is_constrained(index) == false)
    return false;

  const ConstraintLine &p = lines[lines_cache[calculate_line_index(index)]];
  Assert (p.line == index, ExcInternalError());

  // return if an entry for this line was found and if it has only one
  // entry equal to 1.0
  return ((p.entries.size() == 1) &&
          (p.entries[0].second == 1.0));
}


bool ConstraintMatrix::are_identity_constrained (const size_type index1,
                                                 const size_type index2) const
{
  if (is_constrained(index1) == true)
    {
      const ConstraintLine &p = lines[lines_cache[calculate_line_index(index1)]];
      Assert (p.line == index1, ExcInternalError());

      // return if an entry for this line was found and if it has only one
      // entry equal to 1.0 and that one is index2
      return ((p.entries.size() == 1) &&
              (p.entries[0].first == index2) &&
              (p.entries[0].second == 1.0));
    }
  else if (is_constrained(index2) == true)
    {
      const ConstraintLine &p = lines[lines_cache[calculate_line_index(index2)]];
      Assert (p.line == index2, ExcInternalError());

      // return if an entry for this line was found and if it has only one
      // entry equal to 1.0 and that one is index1
      return ((p.entries.size() == 1) &&
              (p.entries[0].first == index1) &&
              (p.entries[0].second == 1.0));
    }
  else
    return false;
}



ConstraintMatrix::size_type
ConstraintMatrix::max_constraint_indirections () const
{
  size_type return_value = 0;
  for (std::vector<ConstraintLine>::const_iterator i=lines.begin();
       i!=lines.end(); ++i)
    // use static cast, since typeof(size)==std::size_t, which is !=
    // size_type on AIX
    return_value = std::max(return_value,
                            static_cast<size_type>(i->entries.size()));

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
  for (size_type i=0; i!=lines.size(); ++i)
    {
      // output the list of constraints as pairs of dofs and their weights
      if (lines[i].entries.size() > 0)
        {
          for (size_type j=0; j<lines[i].entries.size(); ++j)
            out << "    " << lines[i].line
                << " " << lines[i].entries[j].first
                << ":  " << lines[i].entries[j].second << "\n";

          // print out inhomogeneity.
          if (lines[i].inhomogeneity != 0)
            out << "    " << lines[i].line
                << ": " << lines[i].inhomogeneity << "\n";
        }
      else
        // but also output something if the constraint simply reads
        // x[13]=0, i.e. where the right hand side is not a linear
        // combination of other dofs
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
  for (size_type i=0; i!=lines.size(); ++i)
    {
      // same concept as in the previous function
      if (lines[i].entries.size() > 0)
        for (size_type j=0; j<lines[i].entries.size(); ++j)
          out << "  " << lines[i].line << "->" << lines[i].entries[j].first
              << "; // weight: "
              << lines[i].entries[j].second
              << "\n";
      else
        out << "  " << lines[i].line << "\n";
    }
  out << "}" << std::endl;
}



std::size_t
ConstraintMatrix::memory_consumption () const
{
  return (MemoryConsumption::memory_consumption (lines) +
          MemoryConsumption::memory_consumption (lines_cache) +
          MemoryConsumption::memory_consumption (sorted) +
          MemoryConsumption::memory_consumption (local_lines));
}



void
ConstraintMatrix::resolve_indices (std::vector<types::global_dof_index> &indices) const
{
  const unsigned int indices_size = indices.size();
  const std::vector<std::pair<types::global_dof_index,double> > *line_ptr;
  for (unsigned int i=0; i<indices_size; ++i)
    {
      line_ptr = get_constraint_entries(indices[i]);
      // if the index is constraint, the constraints indices are added to the
      // indices vector
      if (line_ptr!=NULL)
        {
          const unsigned int line_size = line_ptr->size();
          for (unsigned int j=0; j<line_size; ++j)
            indices.push_back((*line_ptr)[j].first);
        }
    }

  // keep only the unique elements
  std::sort(indices.begin(),indices.end());
  std::vector<types::global_dof_index>::iterator it;
  it = std::unique(indices.begin(),indices.end());
  indices.resize(it-indices.begin());
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
  template void ConstraintMatrix:: \
  distribute_local_to_global<VectorType > (const Vector<VectorType::value_type>            &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           VectorType                      &, \
                                           const FullMatrix<VectorType::value_type>        &) const;\
  template void ConstraintMatrix:: \
  distribute_local_to_global<VectorType > (const Vector<VectorType::value_type>            &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           VectorType                      &, \
                                           const FullMatrix<VectorType::value_type> &, \
                                           bool) const

#define PARALLEL_VECTOR_FUNCTIONS(VectorType) \
  template void ConstraintMatrix:: \
  distribute_local_to_global<VectorType > (const Vector<VectorType::value_type>            &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           VectorType                      &, \
                                           const FullMatrix<VectorType::value_type>        &) const;\
  template void ConstraintMatrix:: \
  distribute_local_to_global<VectorType > (const Vector<VectorType::value_type>            &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           const std::vector<ConstraintMatrix::size_type>  &, \
                                           VectorType                      &, \
                                           const FullMatrix<VectorType::value_type> &, \
                                           bool) const

#ifdef DEAL_II_WITH_PETSC
VECTOR_FUNCTIONS(PETScWrappers::MPI::Vector);
VECTOR_FUNCTIONS(PETScWrappers::MPI::BlockVector);
#endif

#ifdef DEAL_II_WITH_TRILINOS
PARALLEL_VECTOR_FUNCTIONS(TrilinosWrappers::MPI::Vector);
PARALLEL_VECTOR_FUNCTIONS(TrilinosWrappers::MPI::BlockVector);
#endif

#define MATRIX_VECTOR_FUNCTIONS(MatrixType, VectorType) \
  template void ConstraintMatrix:: \
  distribute_local_to_global<MatrixType,VectorType > (const FullMatrix<MatrixType::value_type>        &, \
                                                      const Vector<VectorType::value_type>            &, \
                                                      const std::vector<ConstraintMatrix::size_type> &, \
                                                      MatrixType                      &, \
                                                      VectorType                      &, \
                                                      bool                             , \
                                                      internal::bool2type<false>) const
#define MATRIX_FUNCTIONS(MatrixType) \
  template void ConstraintMatrix:: \
  distribute_local_to_global<MatrixType,Vector<MatrixType::value_type> > (const FullMatrix<MatrixType::value_type>        &, \
      const Vector<MatrixType::value_type>            &, \
      const std::vector<ConstraintMatrix::size_type> &, \
      MatrixType                      &, \
      Vector<MatrixType::value_type>                  &, \
      bool                             , \
      internal::bool2type<false>) const
#define BLOCK_MATRIX_VECTOR_FUNCTIONS(MatrixType, VectorType)   \
  template void ConstraintMatrix:: \
  distribute_local_to_global<MatrixType,VectorType > (const FullMatrix<MatrixType::value_type>        &, \
                                                      const Vector<VectorType::value_type>            &, \
                                                      const std::vector<ConstraintMatrix::size_type> &, \
                                                      MatrixType                      &, \
                                                      VectorType                      &, \
                                                      bool                             , \
                                                      internal::bool2type<true>) const
#define BLOCK_MATRIX_FUNCTIONS(MatrixType)      \
  template void ConstraintMatrix:: \
  distribute_local_to_global<MatrixType,Vector<MatrixType::value_type> > (const FullMatrix<MatrixType::value_type>        &, \
      const Vector<MatrixType::value_type>            &, \
      const std::vector<ConstraintMatrix::size_type> &, \
      MatrixType                      &, \
      Vector<MatrixType::value_type>                  &, \
      bool                             , \
      internal::bool2type<true>) const

MATRIX_FUNCTIONS(SparseMatrix<long double>);
MATRIX_FUNCTIONS(SparseMatrix<double>);
MATRIX_FUNCTIONS(SparseMatrix<float>);
MATRIX_FUNCTIONS(FullMatrix<double>);
MATRIX_FUNCTIONS(FullMatrix<float>);
MATRIX_FUNCTIONS(FullMatrix<std::complex<double> >);
MATRIX_FUNCTIONS(SparseMatrix<std::complex<double> >);
MATRIX_FUNCTIONS(SparseMatrix<std::complex<long double> >);
MATRIX_FUNCTIONS(SparseMatrix<std::complex<float> >);

BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrix<double>);
BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrix<float>);
BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrix<double>, BlockVector<double>);
BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrix<float>,  BlockVector<float>);

MATRIX_FUNCTIONS(SparseMatrixEZ<double>);
MATRIX_FUNCTIONS(SparseMatrixEZ<float>);
MATRIX_FUNCTIONS(ChunkSparseMatrix<double>);
MATRIX_FUNCTIONS(ChunkSparseMatrix<float>);

// BLOCK_MATRIX_FUNCTIONS(BlockSparseMatrixEZ<double>);
// BLOCK_MATRIX_VECTOR_FUNCTIONS(BlockSparseMatrixEZ<float>,  Vector<float>);

#ifdef DEAL_II_WITH_PETSC
MATRIX_FUNCTIONS(PETScWrappers::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(PETScWrappers::BlockSparseMatrix);
MATRIX_FUNCTIONS(PETScWrappers::MPI::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(PETScWrappers::MPI::BlockSparseMatrix);
MATRIX_VECTOR_FUNCTIONS(PETScWrappers::SparseMatrix, PETScWrappers::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(PETScWrappers::BlockSparseMatrix, PETScWrappers::BlockVector);
MATRIX_VECTOR_FUNCTIONS(PETScWrappers::MPI::SparseMatrix, PETScWrappers::MPI::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(PETScWrappers::MPI::BlockSparseMatrix ,PETScWrappers::MPI::BlockVector);
#endif

#ifdef DEAL_II_WITH_TRILINOS
MATRIX_FUNCTIONS(TrilinosWrappers::SparseMatrix);
BLOCK_MATRIX_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix);
MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::SparseMatrix, TrilinosWrappers::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix, TrilinosWrappers::BlockVector);
MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::SparseMatrix, TrilinosWrappers::MPI::Vector);
BLOCK_MATRIX_VECTOR_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix, TrilinosWrappers::MPI::BlockVector);
#endif


#define SPARSITY_FUNCTIONS(SparsityPatternType)                                      \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityPatternType> ( \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      SparsityPatternType &,                                                         \
      const bool,                                                                    \
      const Table<2,bool> &,                                                         \
      internal::bool2type<false>) const;                                             \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityPatternType> ( \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      SparsityPatternType &,                                                         \
      const bool,                                                                    \
      const Table<2,bool> &) const
#define BLOCK_SPARSITY_FUNCTIONS(SparsityPatternType)                                \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityPatternType> ( \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      SparsityPatternType &,                                                         \
      const bool,                                                                    \
      const Table<2,bool> &,                                                         \
      internal::bool2type<true>) const;                                              \
  template void ConstraintMatrix::add_entries_local_to_global<SparsityPatternType> ( \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      const std::vector<ConstraintMatrix::size_type> &,                              \
      SparsityPatternType &,                                                         \
      const bool,                                                                    \
      const Table<2,bool> &) const

SPARSITY_FUNCTIONS(SparsityPattern);
SPARSITY_FUNCTIONS(DynamicSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockSparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(BlockDynamicSparsityPattern);

#ifdef DEAL_II_WITH_TRILINOS
SPARSITY_FUNCTIONS(TrilinosWrappers::SparsityPattern);
BLOCK_SPARSITY_FUNCTIONS(TrilinosWrappers::BlockSparsityPattern);
#endif


#define ONLY_MATRIX_FUNCTIONS(MatrixType) \
  template void ConstraintMatrix::distribute_local_to_global<MatrixType > (\
      const FullMatrix<MatrixType::value_type>        &, \
      const std::vector<ConstraintMatrix::size_type> &, \
      const std::vector<ConstraintMatrix::size_type> &, \
      MatrixType                      &) const

ONLY_MATRIX_FUNCTIONS(FullMatrix<float>);
ONLY_MATRIX_FUNCTIONS(FullMatrix<double>);
ONLY_MATRIX_FUNCTIONS(SparseMatrix<float>);
ONLY_MATRIX_FUNCTIONS(SparseMatrix<double>);
ONLY_MATRIX_FUNCTIONS(MatrixBlock<SparseMatrix<float> >);
ONLY_MATRIX_FUNCTIONS(MatrixBlock<SparseMatrix<double> >);
ONLY_MATRIX_FUNCTIONS(BlockSparseMatrix<float>);
ONLY_MATRIX_FUNCTIONS(BlockSparseMatrix<double>);

#ifdef DEAL_II_WITH_TRILINOS
ONLY_MATRIX_FUNCTIONS(TrilinosWrappers::SparseMatrix);
ONLY_MATRIX_FUNCTIONS(TrilinosWrappers::BlockSparseMatrix);
#endif

#ifdef DEAL_II_WITH_PETSC
ONLY_MATRIX_FUNCTIONS(PETScWrappers::SparseMatrix);
ONLY_MATRIX_FUNCTIONS(PETScWrappers::BlockSparseMatrix);
ONLY_MATRIX_FUNCTIONS(PETScWrappers::MPI::SparseMatrix);
ONLY_MATRIX_FUNCTIONS(PETScWrappers::MPI::BlockSparseMatrix);
#endif

#include "constraint_matrix.inst"

// allocate scratch data. Cannot use the generic template instantiation
// because we need to provide an initializer object of type
// internals::ConstraintMatrixData<Number> that can be passed to the
// constructor of scratch_data (it won't allow one to be constructed in place).
namespace internals
{
#define SCRATCH_INITIALIZER(Number,Name)                                \
  ConstraintMatrixData<Number>::ScratchData scratch_data_initializer_##Name; \
  template<> Threads::ThreadLocalStorage<ConstraintMatrixData<Number>::ScratchData> \
  ConstraintMatrixData<Number>::scratch_data(scratch_data_initializer_##Name)

  SCRATCH_INITIALIZER(double,double);
  SCRATCH_INITIALIZER(float,float);
  SCRATCH_INITIALIZER(long double,ldouble);
  SCRATCH_INITIALIZER(std::complex<double>,cdouble);
  SCRATCH_INITIALIZER(std::complex<float>,cfloat);
  SCRATCH_INITIALIZER(std::complex<long double>,cldouble);
#undef SCRATCH_INITIALIZER
}


DEAL_II_NAMESPACE_CLOSE
