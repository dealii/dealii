// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2008 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <set>

DEAL_II_NAMESPACE_OPEN



template <typename ForwardIterator>
void
DynamicSparsityPattern::Line::add_entries(ForwardIterator begin,
                                          ForwardIterator end,
                                          const bool      indices_are_sorted)
{
  const int n_elements = end - begin;
  if (n_elements <= 0)
    return;

  const size_type stop_size = entries.size() + n_elements;

  if (indices_are_sorted == true && n_elements > 3)
    {
      // in debug mode, check whether the
      // indices really are sorted.
      if constexpr (running_in_debug_mode())
        {
          {
            ForwardIterator test = begin, test1 = begin;
            ++test1;
            for (; test1 != end; ++test, ++test1)
              Assert(*test1 > *test, ExcInternalError());
          }
        }

      if (entries.empty() || entries.back() < *begin)
        {
          entries.insert(entries.end(), begin, end);
          return;
        }

      // find a possible insertion point for
      // the first entry. check whether the
      // first entry is a duplicate before
      // actually doing something.
      ForwardIterator                  my_it = begin;
      size_type                        col   = *my_it;
      std::vector<size_type>::iterator it =
        Utilities::lower_bound(entries.begin(), entries.end(), col);
      while (*it == col)
        {
          ++my_it;
          if (my_it == end)
            break;
          col = *my_it;
          // check the very next entry in the
          // current array
          ++it;
          if (it == entries.end())
            break;
          if (*it > col)
            break;
          if (*it == col)
            continue;
          // ok, it wasn't the very next one, do a
          // binary search to find the insert point
          it = Utilities::lower_bound(it, entries.end(), col);
          if (it == entries.end())
            break;
        }
      // all input entries were duplicates.
      if (my_it == end)
        return;

      // resize vector by just inserting the
      // list
      const size_type pos1 = it - entries.begin();
      Assert(pos1 <= entries.size(), ExcInternalError());
      entries.insert(it, my_it, end);
      it = entries.begin() + pos1;
      Assert(entries.size() >= static_cast<size_type>(it - entries.begin()),
             ExcInternalError());

      // now merge the two lists.
      std::vector<size_type>::iterator it2 = it + (end - my_it);

      // as long as there are indices both in
      // the end of the entries list and in the
      // input list
      while (my_it != end && it2 != entries.end())
        {
          if (*my_it < *it2)
            *it++ = *my_it++;
          else if (*my_it == *it2)
            {
              *it++ = *it2++;
              ++my_it;
            }
          else
            *it++ = *it2++;
        }
      // in case there are indices left in the
      // input list
      while (my_it != end)
        *it++ = *my_it++;

      // in case there are indices left in the
      // end of entries
      while (it2 != entries.end())
        *it++ = *it2++;

      // resize and return
      const size_type new_size = it - entries.begin();
      Assert(new_size <= stop_size, ExcInternalError());
      entries.resize(new_size);
      return;
    }

  // unsorted case or case with too few
  // elements
  ForwardIterator my_it = begin;

  // If necessary, increase the size of the
  // array.
  if (stop_size > entries.capacity())
    entries.reserve(stop_size);

  size_type                        col = *my_it;
  std::vector<size_type>::iterator it, it2;
  // insert the first element as for one
  // entry only first check the last
  // element (or if line is still empty)
  if ((entries.empty()) || (entries.back() < col))
    {
      entries.push_back(col);
      it = entries.end() - 1;
    }
  else
    {
      // do a binary search to find the place
      // where to insert:
      it2 = Utilities::lower_bound(entries.begin(), entries.end(), col);

      // If this entry is a duplicate, continue
      // immediately Insert at the right place
      // in the vector. Vector grows
      // automatically to fit elements. Always
      // doubles its size.
      if (*it2 != col)
        it = entries.insert(it2, col);
      else
        it = it2;
    }

  ++my_it;
  // Now try to be smart and insert with
  // bias in the direction we are
  // walking. This has the advantage that
  // for sorted lists, we always search in
  // the right direction, what should
  // decrease the work needed in here.
  for (; my_it != end; ++my_it)
    {
      col = *my_it;
      // need a special insertion command when
      // we're at the end of the list
      if (col > entries.back())
        {
          entries.push_back(col);
          it = entries.end() - 1;
        }
      // search to the right (preferred search
      // direction)
      else if (col > *it)
        {
          it2 = Utilities::lower_bound(it++, entries.end(), col);
          if (*it2 != col)
            it = entries.insert(it2, col);
        }
      // search to the left
      else if (col < *it)
        {
          it2 = Utilities::lower_bound(entries.begin(), it, col);
          if (*it2 != col)
            it = entries.insert(it2, col);
        }
      // if we're neither larger nor smaller,
      // then this was a duplicate and we can
      // just continue.
    }
}


DynamicSparsityPattern::size_type
DynamicSparsityPattern::Line::memory_consumption() const
{
  return entries.capacity() * sizeof(size_type) + sizeof(Line);
}


DynamicSparsityPattern::DynamicSparsityPattern()
  : SparsityPatternBase()
  , have_entries(false)
  , rowset(0)
{}



DynamicSparsityPattern::DynamicSparsityPattern(const DynamicSparsityPattern &s)
  : SparsityPatternBase()
  , have_entries(false)
  , rowset(0)
{
  Assert(s.rows == 0 && s.cols == 0,
         ExcMessage(
           "This constructor can only be called if the provided argument "
           "is the sparsity pattern for an empty matrix. This constructor can "
           "not be used to copy-construct a non-empty sparsity pattern."));
}



DynamicSparsityPattern::DynamicSparsityPattern(const size_type m,
                                               const size_type n,
                                               const IndexSet &rowset_)
  : SparsityPatternBase()
  , have_entries(false)
  , rowset(0)
{
  reinit(m, n, rowset_);
}


DynamicSparsityPattern::DynamicSparsityPattern(const IndexSet &rowset_)
  : DynamicSparsityPattern(rowset_.size(), rowset_.size(), rowset_)
{}


DynamicSparsityPattern::DynamicSparsityPattern(const size_type n)
  : SparsityPatternBase()
  , have_entries(false)
  , rowset(0)
{
  reinit(n, n);
}



DynamicSparsityPattern &
DynamicSparsityPattern::operator=(const DynamicSparsityPattern &s)
{
  Assert(s.n_rows() == 0 && s.n_cols() == 0,
         ExcMessage(
           "This operator can only be called if the provided argument "
           "is the sparsity pattern for an empty matrix. This operator can "
           "not be used to copy a non-empty sparsity pattern."));

  Assert(n_rows() == 0 && n_cols() == 0,
         ExcMessage("This operator can only be called if the current object is "
                    "empty."));

  return *this;
}



void
DynamicSparsityPattern::reinit(const size_type m,
                               const size_type n,
                               const IndexSet &rowset_)
{
  resize(m, n);
  have_entries = false;
  rowset       = rowset_;

  Assert(rowset.size() == 0 || rowset.size() == m,
         ExcMessage(
           "The IndexSet argument to this function needs to either "
           "be empty (indicating the complete set of rows), or have size "
           "equal to the desired number of rows as specified by the "
           "first argument to this function. (Of course, the number "
           "of indices in this IndexSet may be less than the number "
           "of rows, but the *size* of the IndexSet must be equal.)"));

  std::vector<Line> new_lines(rowset.size() == 0 ? n_rows() :
                                                   rowset.n_elements());
  lines.swap(new_lines);
}



void
DynamicSparsityPattern::compress()
{}



bool
DynamicSparsityPattern::empty() const
{
  return ((rows == 0) && (cols == 0));
}



DynamicSparsityPattern::size_type
DynamicSparsityPattern::max_entries_per_row() const
{
  if (!have_entries)
    return 0;

  size_type m = 0;
  for (const auto &line : lines)
    {
      m = std::max(m, static_cast<size_type>(line.entries.size()));
    }

  return m;
}



void
DynamicSparsityPattern::add_row_entries(
  const size_type                  &row,
  const ArrayView<const size_type> &columns,
  const bool                        indices_are_sorted)
{
  add_entries(row, columns.begin(), columns.end(), indices_are_sorted);
}



bool
DynamicSparsityPattern::exists(const size_type i, const size_type j) const
{
  AssertIndexRange(i, rows);
  AssertIndexRange(j, cols);
  Assert(
    rowset.size() == 0 || rowset.is_element(i),
    ExcMessage(
      "The row IndexSet does not contain the index i. This sparsity pattern "
      "object cannot know whether the entry (i, j) exists or not."));

  // Avoid a segmentation fault in below code if the row index happens to
  // not be present in the IndexSet rowset:
  if (!(rowset.size() == 0 || rowset.is_element(i)))
    return false;

  if (!have_entries)
    return false;

  const size_type rowindex =
    rowset.size() == 0 ? i : rowset.index_within_set(i);

  return std::binary_search(lines[rowindex].entries.begin(),
                            lines[rowindex].entries.end(),
                            j);
}



void
DynamicSparsityPattern::symmetrize()
{
  Assert(rows == cols, ExcNotQuadratic());

  // loop over all elements presently in the sparsity pattern and add the
  // transpose element. note:
  //
  // 1. that the sparsity pattern changes which we work on, but not the present
  // row
  //
  // 2. that the @p{add} function can be called on elements that already exist
  // without any harm
  for (size_type row = 0; row < lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size() == 0 ? row : rowset.nth_index_in_set(row);

      for (const size_type row_entry : lines[row].entries)
        // add the transpose entry if this is not the diagonal
        if (rowindex != row_entry)
          add(row_entry, rowindex);
    }
}



void
DynamicSparsityPattern::clear_row(const size_type row)
{
  AssertIndexRange(row, n_rows());
  if (!have_entries)
    return;

  if (rowset.size() > 0 && !rowset.is_element(row))
    return;

  const size_type rowindex =
    rowset.size() == 0 ? row : rowset.index_within_set(row);

  AssertIndexRange(rowindex, lines.size());
  lines[rowindex].entries = std::vector<size_type>();
}



DynamicSparsityPattern
DynamicSparsityPattern::get_view(const IndexSet &rows) const
{
  DynamicSparsityPattern view;
  view.reinit(rows.n_elements(), this->n_cols());
  AssertDimension(rows.size(), this->n_rows());

  const auto                        end      = rows.end();
  DynamicSparsityPattern::size_type view_row = 0;
  for (auto it = rows.begin(); it != end; ++it, ++view_row)
    {
      const size_type rowindex =
        rowset.size() == 0 ? *it : rowset.index_within_set(*it);

      view.lines[view_row].entries = lines[rowindex].entries;
      view.have_entries |= (lines[rowindex].entries.size() > 0);
    }
  return view;
}



template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
void
DynamicSparsityPattern::compute_Tmmult_pattern(
  const SparsityPatternTypeLeft  &sp_A,
  const SparsityPatternTypeRight &sp_B)
{
  Assert(sp_A.n_rows() == sp_B.n_rows(),
         ExcDimensionMismatch(sp_A.n_rows(), sp_B.n_rows()));

  this->reinit(sp_A.n_cols(), sp_B.n_cols());
  // we will go through all the
  // rows in the matrix A, and for each column in a row we add the whole
  // row of matrix B with that row number. This means that we will insert
  // a lot of entries to each row, which is best handled by the
  // DynamicSparsityPattern class.

  std::vector<size_type> new_cols;
  new_cols.reserve(sp_B.max_entries_per_row());

  // C_{kl} = A_{ik} B_{il}
  for (size_type i = 0; i < sp_A.n_rows(); ++i)
    {
      // get all column numbers from sp_B in a temporary vector:
      new_cols.resize(sp_B.row_length(i));
      {
        const auto last_il = sp_B.end(i);
        auto      *col_ptr = new_cols.data();
        for (auto il = sp_B.begin(i); il != last_il; ++il)
          *col_ptr++ = il->column();
      }
      std::sort(new_cols.begin(), new_cols.end());

      // now for each k, add new_cols to the target sparsity
      const auto last_ik = sp_A.end(i);
      for (auto ik = sp_A.begin(i); ik != last_ik; ++ik)
        this->add_entries(ik->column(), new_cols.begin(), new_cols.end(), true);
    }
}



template <typename SparsityPatternTypeLeft, typename SparsityPatternTypeRight>
void
DynamicSparsityPattern::compute_mmult_pattern(
  const SparsityPatternTypeLeft  &left,
  const SparsityPatternTypeRight &right)
{
  Assert(left.n_cols() == right.n_rows(),
         ExcDimensionMismatch(left.n_cols(), right.n_rows()));

  this->reinit(left.n_rows(), right.n_cols());

  typename SparsityPatternTypeLeft::iterator it_left  = left.begin(),
                                             end_left = left.end();
  for (; it_left != end_left; ++it_left)
    {
      const unsigned int j = it_left->column();

      // We are sitting on entry (i,j) of the left sparsity pattern. We then
      // need to add all entries (i,k) to the final sparsity pattern where (j,k)
      // exists in the right sparsity pattern -- i.e., we need to iterate over
      // row j.
      typename SparsityPatternTypeRight::iterator it_right  = right.begin(j),
                                                  end_right = right.end(j);
      for (; it_right != end_right; ++it_right)
        this->add(it_left->row(), it_right->column());
    }
}



void
DynamicSparsityPattern::print(std::ostream &out) const
{
  for (size_type row = 0; row < lines.size(); ++row)
    {
      out << '[' << (rowset.size() == 0 ? row : rowset.nth_index_in_set(row));

      for (const auto entry : lines[row].entries)
        out << ',' << entry;

      out << ']' << std::endl;
    }

  AssertThrow(out.fail() == false, ExcIO());
}



void
DynamicSparsityPattern::print_gnuplot(std::ostream &out) const
{
  for (size_type row = 0; row < lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size() == 0 ? row : rowset.nth_index_in_set(row);

      for (const auto entry : lines[row].entries)
        // while matrix entries are usually
        // written (i,j), with i vertical and
        // j horizontal, gnuplot output is
        // x-y, that is we have to exchange
        // the order of output
        out << entry << " " << -static_cast<signed int>(rowindex) << std::endl;
    }


  AssertThrow(out.fail() == false, ExcIO());
}



DynamicSparsityPattern::size_type
DynamicSparsityPattern::bandwidth() const
{
  size_type b = 0;
  for (size_type row = 0; row < lines.size(); ++row)
    {
      const size_type rowindex =
        rowset.size() == 0 ? row : rowset.nth_index_in_set(row);

      for (const auto entry : lines[row].entries)
        if (static_cast<size_type>(
              std::abs(static_cast<int>(rowindex - entry))) > b)
          b = std::abs(static_cast<signed int>(rowindex - entry));
    }

  return b;
}



DynamicSparsityPattern::size_type
DynamicSparsityPattern::n_nonzero_elements() const
{
  if (!have_entries)
    return 0;

  size_type n = 0;
  for (const auto &line : lines)
    {
      n += line.entries.size();
    }

  return n;
}



IndexSet
DynamicSparsityPattern::nonempty_cols() const
{
  std::set<types::global_dof_index> cols;
  for (const auto &line : lines)
    cols.insert(line.entries.begin(), line.entries.end());

  IndexSet res(this->n_cols());
  res.add_indices(cols.begin(), cols.end());
  return res;
}



IndexSet
DynamicSparsityPattern::nonempty_rows() const
{
  const IndexSet  all_rows            = complete_index_set(this->n_rows());
  const IndexSet &locally_stored_rows = rowset.size() == 0 ? all_rows : rowset;

  std::vector<types::global_dof_index> rows;
  auto                                 line = lines.begin();
  AssertDimension(locally_stored_rows.n_elements(), lines.size());
  for (const auto row : locally_stored_rows)
    {
      if (line->entries.size() > 0)
        rows.push_back(row);

      ++line;
    }

  IndexSet res(this->n_rows());
  res.add_indices(rows.begin(), rows.end());
  return res;
}



DynamicSparsityPattern::size_type
DynamicSparsityPattern::memory_consumption() const
{
  size_type mem = sizeof(DynamicSparsityPattern) +
                  MemoryConsumption::memory_consumption(rowset) -
                  sizeof(rowset);

  for (const auto &line : lines)
    mem += MemoryConsumption::memory_consumption(line);

  return mem;
}



types::global_dof_index
DynamicSparsityPattern::column_index(
  const DynamicSparsityPattern::size_type row,
  const DynamicSparsityPattern::size_type col) const
{
  AssertIndexRange(row, n_rows());
  AssertIndexRange(col, n_cols());
  Assert(rowset.size() == 0 || rowset.is_element(row), ExcInternalError());

  const DynamicSparsityPattern::size_type local_row =
    rowset.size() != 0u ? rowset.index_within_set(row) : row;

  // now we need to do a binary search. Note that col indices are assumed to
  // be sorted.
  const auto &cols = lines[local_row].entries;
  auto        it   = Utilities::lower_bound(cols.begin(), cols.end(), col);

  if ((it != cols.end()) && (*it == col))
    return (it - cols.begin());
  else
    return numbers::invalid_size_type;
}



// explicit instantiations
template void
DynamicSparsityPattern::Line::add_entries(size_type *, size_type *, const bool);
template void
DynamicSparsityPattern::Line::add_entries(const size_type *,
                                          const size_type *,
                                          const bool);
#ifndef DEAL_II_VECTOR_ITERATOR_IS_POINTER
template void
DynamicSparsityPattern::Line::add_entries(std::vector<size_type>::iterator,
                                          std::vector<size_type>::iterator,
                                          const bool);
template void
DynamicSparsityPattern::Line::add_entries(
  std::vector<size_type>::const_iterator,
  std::vector<size_type>::const_iterator,
  const bool);
#endif

template void
DynamicSparsityPattern::compute_mmult_pattern(const DynamicSparsityPattern &,
                                              const DynamicSparsityPattern &);
template void
DynamicSparsityPattern::compute_mmult_pattern(const DynamicSparsityPattern &,
                                              const SparsityPattern &);
template void
DynamicSparsityPattern::compute_mmult_pattern(const SparsityPattern &,
                                              const DynamicSparsityPattern &);
template void
DynamicSparsityPattern::compute_mmult_pattern(const SparsityPattern &,
                                              const SparsityPattern &);

template void
DynamicSparsityPattern::compute_Tmmult_pattern(const SparsityPattern &,
                                               const SparsityPattern &);
template void
DynamicSparsityPattern::compute_Tmmult_pattern(const DynamicSparsityPattern &,
                                               const SparsityPattern &);
template void
DynamicSparsityPattern::compute_Tmmult_pattern(const SparsityPattern &,
                                               const DynamicSparsityPattern &);
template void
DynamicSparsityPattern::compute_Tmmult_pattern(const DynamicSparsityPattern &,
                                               const DynamicSparsityPattern &);

DEAL_II_NAMESPACE_CLOSE
