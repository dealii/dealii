// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>

#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/trilinos_tpetra_types.h>

#include <boost/container/small_vector.hpp>

#include <vector>

#ifdef DEAL_II_WITH_TRILINOS
DEAL_II_DISABLE_EXTRA_DIAGNOSTICS
#  ifdef DEAL_II_WITH_MPI
#    include <Epetra_MpiComm.h>
#  endif
#  include <Epetra_Map.h>
#  include <Epetra_SerialComm.h>
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
#    include <Tpetra_Map.hpp>
#  endif
DEAL_II_ENABLE_EXTRA_DIAGNOSTICS
#endif

DEAL_II_NAMESPACE_OPEN



#ifdef DEAL_II_WITH_TRILINOS

#  ifdef DEAL_II_TRILINOS_WITH_TPETRA

template <typename NodeType>
IndexSet::IndexSet(
  const Teuchos::RCP<
    const Tpetra::Map<int, types::signed_global_dof_index, NodeType>> &map)
  : is_compressed(true)
  , index_space_size(1 + map->getMaxAllGlobalIndex())
  , largest_range(numbers::invalid_unsigned_int)
{
  Assert(map->getMinAllGlobalIndex() == 0,
         ExcMessage(
           "The Tpetra::Map does not contain the global index 0, "
           "which means some entries are not present on any processor."));

  // For a contiguous map, we do not need to go through the whole data...
  if (map->isContiguous())
    add_range(size_type(map->getMinGlobalIndex()),
              size_type(map->getMaxGlobalIndex() + 1));
  else
    {
#    if DEAL_II_TRILINOS_VERSION_GTE(13, 4, 0)
      const size_type n_indices = map->getLocalNumElements();
#    else
      const size_type n_indices = map->getNodeNumElements();
#    endif
      const types::signed_global_dof_index *indices =
        map->getMyGlobalIndices().data();
      add_indices(indices, indices + n_indices);
    }
  compress();
}

#  endif // DEAL_II_TRILINOS_WITH_TPETRA


// the 64-bit path uses a few different names, so put that into a separate
// implementation

#  ifdef DEAL_II_WITH_64BIT_INDICES

IndexSet::IndexSet(const Epetra_BlockMap &map)
  : is_compressed(true)
  , index_space_size(1 + map.MaxAllGID64())
  , largest_range(numbers::invalid_unsigned_int)
{
  Assert(map.MinAllGID64() == 0,
         ExcMessage(
           "The Epetra_BlockMap does not contain the global index 0, "
           "which means some entries are not present on any processor."));

  // For a contiguous map, we do not need to go through the whole data...
  if (map.LinearMap())
    add_range(size_type(map.MinMyGID64()), size_type(map.MaxMyGID64() + 1));
  else
    {
      const size_type n_indices = map.NumMyElements();
      size_type      *indices =
        reinterpret_cast<size_type *>(map.MyGlobalElements64());
      add_indices(indices, indices + n_indices);
    }
  compress();
}

#  else

// this is the standard 32-bit implementation

IndexSet::IndexSet(const Epetra_BlockMap &map)
  : is_compressed(true)
  , index_space_size(1 + map.MaxAllGID())
  , largest_range(numbers::invalid_unsigned_int)
{
  Assert(map.MinAllGID() == 0,
         ExcMessage(
           "The Epetra_BlockMap does not contain the global index 0, "
           "which means some entries are not present on any processor."));

  // For a contiguous map, we do not need to go through the whole data...
  if (map.LinearMap())
    add_range(size_type(map.MinMyGID()), size_type(map.MaxMyGID() + 1));
  else
    {
      const size_type n_indices = map.NumMyElements();
      unsigned int   *indices =
        reinterpret_cast<unsigned int *>(map.MyGlobalElements());
      add_indices(indices, indices + n_indices);
    }
  compress();
}

#  endif

#endif // ifdef DEAL_II_WITH_TRILINOS



void
IndexSet::do_compress() const
{
  {
    // we will, in the following, modify mutable variables. this can only
    // work in multithreaded applications if we lock the data structures
    // via a mutex, so that users can call 'const' functions from threads
    // in parallel (and these 'const' functions can then call compress()
    // which itself calls the current function)
    std::lock_guard<std::mutex> lock(compress_mutex);

    // see if any of the contiguous ranges can be merged. do not use
    // std::vector::erase in-place as it is quadratic in the number of
    // ranges. since the ranges are sorted by their first index, determining
    // overlap isn't all that hard
    std::vector<Range>::iterator store = ranges.begin();
    for (std::vector<Range>::iterator i = ranges.begin(); i != ranges.end();)
      {
        std::vector<Range>::iterator next = i;
        ++next;

        size_type first_index = i->begin;
        size_type last_index  = i->end;

        // see if we can merge any of the following ranges
        while (next != ranges.end() && (next->begin <= last_index))
          {
            last_index = std::max(last_index, next->end);
            ++next;
          }
        i = next;

        // store the new range in the slot we last occupied
        *store = Range(first_index, last_index);
        ++store;
      }
    // use a compact array with exactly the right amount of storage
    if (store != ranges.end())
      {
        std::vector<Range> new_ranges(ranges.begin(), store);
        ranges.swap(new_ranges);
      }

    // now compute indices within set and the range with most elements
    size_type next_index = 0, largest_range_size = 0;
    for (std::vector<Range>::iterator i = ranges.begin(); i != ranges.end();
         ++i)
      {
        Assert(i->begin < i->end, ExcInternalError());

        i->nth_index_in_set = next_index;
        next_index += (i->end - i->begin);
        if (i->end - i->begin > largest_range_size)
          {
            largest_range_size = i->end - i->begin;
            largest_range      = i - ranges.begin();
          }
      }
    is_compressed = true;

    // check that next_index is correct. needs to be after the previous
    // statement because we otherwise will get into an endless loop
    Assert(next_index == n_elements(), ExcInternalError());
  }

  if constexpr (running_in_debug_mode())
    {
      // A consistency check: We should only ever have added indices
      // that are within the range of the index set. Instead of doing
      // this in every one of the many functions that add indices,
      // do this in the current, central location
      for (const auto &range : ranges)
        Assert((range.begin < index_space_size) &&
                 (range.end <= index_space_size),
               ExcMessage(
                 "In the process of creating the current IndexSet "
                 "object, you added indices beyond the size of the index "
                 "space. Specifically, you added elements that form the "
                 "range [" +
                 std::to_string(range.begin) + "," + std::to_string(range.end) +
                 "), but the size of the index space is only " +
                 std::to_string(index_space_size) + "."));
    }
}



#ifndef DOXYGEN
IndexSet
IndexSet::operator&(const IndexSet &is) const
{
  Assert(size() == is.size(), ExcDimensionMismatch(size(), is.size()));

  compress();
  is.compress();

  std::vector<Range>::const_iterator r1 = ranges.begin(),
                                     r2 = is.ranges.begin();
  IndexSet result(size());

  while ((r1 != ranges.end()) && (r2 != is.ranges.end()))
    {
      // if r1 and r2 do not overlap at all, then move the pointer that sits
      // to the left of the other up by one
      if (r1->end <= r2->begin)
        ++r1;
      else if (r2->end <= r1->begin)
        ++r2;
      else
        {
          // the ranges must overlap somehow
          Assert(((r1->begin <= r2->begin) && (r1->end > r2->begin)) ||
                   ((r2->begin <= r1->begin) && (r2->end > r1->begin)),
                 ExcInternalError());

          // add the overlapping range to the result
          result.add_range(std::max(r1->begin, r2->begin),
                           std::min(r1->end, r2->end));

          // now move that iterator that ends earlier one up. note that it has
          // to be this one because a subsequent range may still have a chance
          // of overlapping with the range that ends later
          if (r1->end <= r2->end)
            ++r1;
          else
            ++r2;
        }
    }

  result.compress();
  return result;
}
#endif



IndexSet
IndexSet::get_view(const size_type begin, const size_type end) const
{
  Assert(begin <= end,
         ExcMessage("End index needs to be larger or equal to begin index!"));
  Assert(end <= size(),
         ExcMessage("You are asking for a view into an IndexSet object "
                    "that would cover the sub-range [" +
                    std::to_string(begin) + ',' + std::to_string(end) +
                    "). But this is not a subset of the range "
                    "of the current object, which is [0," +
                    std::to_string(size()) + ")."));

  IndexSet                           result(end - begin);
  std::vector<Range>::const_iterator r1 = ranges.begin();

  while (r1 != ranges.end())
    {
      if ((r1->end > begin) && (r1->begin < end))
        {
          result.add_range(std::max(r1->begin, begin) - begin,
                           std::min(r1->end, end) - begin);
        }
      else if (r1->begin >= end)
        break;

      ++r1;
    }

  result.compress();
  return result;
}



IndexSet
IndexSet::get_view(const IndexSet &mask) const
{
  Assert(size() == mask.size(),
         ExcMessage("The mask must have the same size index space "
                    "as the index set it is applied to."));

  // If 'other' is an empty set, then the view is also empty:
  if (mask == IndexSet())
    return {};

  // For everything, it is more efficient to work on compressed sets:
  compress();
  mask.compress();

  // If 'other' has a single range, then we can just defer to the
  // previous function
  if (mask.ranges.size() == 1)
    return get_view(mask.ranges[0].begin, mask.ranges[0].end);

  // For the general case where the mask is an arbitrary set,
  // the situation is slightly more complicated. We need to walk
  // the ranges of the two index sets in parallel and search for
  // overlaps, and then appropriately shift

  // we save all new ranges to our IndexSet in an temporary vector and
  // add all of them in one go at the end.
  std::vector<Range> new_ranges;

  std::vector<Range>::iterator own_it  = ranges.begin();
  std::vector<Range>::iterator mask_it = mask.ranges.begin();

  while ((own_it != ranges.end()) && (mask_it != mask.ranges.end()))
    {
      // If our own range lies completely ahead of the current
      // range in the mask, move forward and start the loop body
      // anew. If this was the last range, the 'while' loop above
      // will terminate, so we don't have to check for end iterators
      if (own_it->end <= mask_it->begin)
        {
          ++own_it;
          continue;
        }

      // Do the same if the current mask range lies completely ahead of
      // the current range of the this object:
      if (mask_it->end <= own_it->begin)
        {
          ++mask_it;
          continue;
        }

      // Now own_it and other_it overlap. Check that that is true by
      // enumerating the cases that can happen. This is
      // surprisingly tricky because the two intervals can intersect in
      // a number of different ways, but there really are only the four
      // following possibilities:

      // Case 1: our interval overlaps the left end of the other interval
      //
      // So we need to add the elements from the first element of the mask's
      // interval to the end of our own interval. But we need to shift the
      // indices so that they correspond to the how many'th element within the
      // mask this is; fortunately (because we compressed the mask), this
      // is recorded in the mask's ranges.
      if ((own_it->begin <= mask_it->begin) && (own_it->end <= mask_it->end))
        {
          new_ranges.emplace_back(mask_it->begin - mask_it->nth_index_in_set,
                                  own_it->end - mask_it->nth_index_in_set);
        }
      else
        // Case 2:our interval overlaps the tail end of the other interval
        if ((mask_it->begin <= own_it->begin) && (mask_it->end <= own_it->end))
          {
            const size_type offset_within_mask_interval =
              own_it->begin - mask_it->begin;
            new_ranges.emplace_back(mask_it->nth_index_in_set +
                                      offset_within_mask_interval,
                                    mask_it->nth_index_in_set +
                                      (mask_it->end - mask_it->begin));
          }
        else
          // Case 3: Our own interval completely encloses the other interval
          if ((own_it->begin <= mask_it->begin) &&
              (own_it->end >= mask_it->end))
            {
              new_ranges.emplace_back(mask_it->begin -
                                        mask_it->nth_index_in_set,
                                      mask_it->end - mask_it->nth_index_in_set);
            }
          else
            // Case 3: The other interval completely encloses our own interval
            if ((mask_it->begin <= own_it->begin) &&
                (mask_it->end >= own_it->end))
              {
                const size_type offset_within_mask_interval =
                  own_it->begin - mask_it->begin;
                new_ranges.emplace_back(mask_it->nth_index_in_set +
                                          offset_within_mask_interval,
                                        mask_it->nth_index_in_set +
                                          offset_within_mask_interval +
                                          (own_it->end - own_it->begin));
              }
            else
              DEAL_II_ASSERT_UNREACHABLE();

      // We considered the overlap of these two intervals. It may of course
      // be that one of them overlaps with another one, but that can only
      // be the case for the interval that extends further to the right. So
      // we can safely move on from the interval that terminates earlier:
      if (own_it->end < mask_it->end)
        ++own_it;
      else if (mask_it->end < own_it->end)
        ++mask_it;
      else
        {
          // The intervals ended at the same point. We can move on from both.
          // (The algorithm would also work if we only moved on from one,
          // but we can micro-optimize here without too much effort.)
          ++own_it;
          ++mask_it;
        }
    }

  // Now turn the ranges of overlap we have accumulated into an IndexSet in
  // its own right:
  IndexSet result(mask.n_elements());
  for (const auto &range : new_ranges)
    result.add_range(range.begin, range.end);
  result.compress();

  return result;
}



std::vector<IndexSet>
IndexSet::split_by_block(
  const std::vector<types::global_dof_index> &n_indices_per_block) const
{
  std::vector<IndexSet> partitioned;
  const unsigned int    n_blocks = n_indices_per_block.size();

  partitioned.reserve(n_blocks);
  types::global_dof_index start = 0;
  for (const auto n_block_indices : n_indices_per_block)
    {
      partitioned.push_back(this->get_view(start, start + n_block_indices));
      start += n_block_indices;
    }

  if constexpr (running_in_debug_mode())
    {
      types::global_dof_index sum = 0;
      for (const auto &partition : partitioned)
        {
          sum += partition.size();
        }
      AssertDimension(sum, this->size());
    }

  return partitioned;
}



void
IndexSet::subtract_set(const IndexSet &other)
{
  compress();
  other.compress();
  is_compressed = false;


  // we save all new ranges to our IndexSet in an temporary vector and
  // add all of them in one go at the end.
  std::vector<Range> new_ranges;

  std::vector<Range>::iterator own_it   = ranges.begin();
  std::vector<Range>::iterator other_it = other.ranges.begin();

  while (own_it != ranges.end() && other_it != other.ranges.end())
    {
      // advance own iterator until we get an overlap
      if (own_it->end <= other_it->begin)
        {
          new_ranges.push_back(*own_it);
          ++own_it;
          continue;
        }
      // we are done with other_it, so advance
      if (own_it->begin >= other_it->end)
        {
          ++other_it;
          continue;
        }

      // Now own_it and other_it overlap.  First save the part of own_it that
      // is before other_it (if not empty).
      if (own_it->begin < other_it->begin)
        {
          Range r(own_it->begin, other_it->begin);
          r.nth_index_in_set = 0; // fix warning of unused variable
          new_ranges.push_back(r);
        }
      // change own_it to the sub range behind other_it. Do not delete own_it
      // in any case. As removal would invalidate iterators, we just shrink
      // the range to an empty one.
      own_it->begin = other_it->end;
      if (own_it->begin > own_it->end)
        {
          own_it->begin = own_it->end;
          ++own_it;
        }

      // continue without advancing iterators, the right one will be advanced
      // next.
    }

  // make sure to take over the remaining ranges
  for (; own_it != ranges.end(); ++own_it)
    new_ranges.push_back(*own_it);

  ranges.clear();

  // done, now add the temporary ranges
  const std::vector<Range>::iterator end = new_ranges.end();
  for (std::vector<Range>::iterator it = new_ranges.begin(); it != end; ++it)
    add_range(it->begin, it->end);

  compress();
}



IndexSet
IndexSet::tensor_product(const IndexSet &other) const
{
  IndexSet set(this->size() * other.size());
  for (const auto el : *this)
    set.add_indices(other, el * other.size());
  set.compress();
  return set;
}



IndexSet::size_type
IndexSet::pop_back()
{
  Assert(is_empty() == false,
         ExcMessage(
           "pop_back() failed, because this IndexSet contains no entries."));

  const size_type index = ranges.back().end - 1;
  --ranges.back().end;

  if (ranges.back().begin == ranges.back().end)
    ranges.pop_back();

  return index;
}



IndexSet::size_type
IndexSet::pop_front()
{
  Assert(is_empty() == false,
         ExcMessage(
           "pop_front() failed, because this IndexSet contains no entries."));

  const size_type index = ranges.front().begin;
  ++ranges.front().begin;

  if (ranges.front().begin == ranges.front().end)
    ranges.erase(ranges.begin());

  // We have to set this in any case, because nth_index_in_set is no longer
  // up to date for all but the first range
  is_compressed = false;

  return index;
}



void
IndexSet::add_range_lower_bound(const Range &new_range)
{
  // if the inserted range is already within the range we find by lower_bound,
  // there is no need to do anything; we do not try to be clever here and
  // leave all other work to compress().
  const auto insert_position =
    Utilities::lower_bound(ranges.begin(), ranges.end(), new_range);
  if (insert_position == ranges.end() ||
      insert_position->begin > new_range.begin ||
      insert_position->end < new_range.end)
    ranges.insert(insert_position, new_range);
}



void
IndexSet::add_ranges_internal(
  boost::container::small_vector<std::pair<size_type, size_type>, 200>
            &tmp_ranges,
  const bool ranges_are_sorted)
{
  if (!ranges_are_sorted)
    std::sort(tmp_ranges.begin(), tmp_ranges.end());

  // if we have many ranges, we first construct a temporary index set (where
  // we add ranges in a consecutive way, so fast), otherwise, we work with
  // add_range(). the number 9 is chosen heuristically given the fact that
  // there are typically up to 8 independent ranges when adding the degrees of
  // freedom on a 3d cell or 9 when adding degrees of freedom of faces. if
  // doing cell-by-cell additions, we want to avoid repeated calls to
  // IndexSet::compress() which gets called upon merging two index sets, so we
  // want to be in the other branch then.
  if (tmp_ranges.size() > 9)
    {
      IndexSet tmp_set(size());
      tmp_set.ranges.reserve(tmp_ranges.size());
      for (const auto &i : tmp_ranges)
        tmp_set.add_range(i.first, i.second);

      // Case if we have zero or just one range: Add into the other set with
      // its indices, as that is cheaper
      if (this->ranges.size() <= 1)
        {
          if (this->ranges.size() == 1)
            tmp_set.add_range(ranges[0].begin, ranges[0].end);
          std::swap(*this, tmp_set);
        }
      else
        this->add_indices(tmp_set);
    }
  else
    for (const auto &i : tmp_ranges)
      add_range(i.first, i.second);
}



void
IndexSet::add_indices(const IndexSet &other, const size_type offset)
{
  if ((this == &other) && (offset == 0))
    return;

  if (other.ranges.size() != 0)
    {
      AssertIndexRange(other.ranges.back().end - 1, index_space_size);
    }

  compress();
  other.compress();

  std::vector<Range>::const_iterator r1 = ranges.begin(),
                                     r2 = other.ranges.begin();

  std::vector<Range> new_ranges;
  // just get the start and end of the ranges right in this method, everything
  // else will be done in compress()
  while (r1 != ranges.end() || r2 != other.ranges.end())
    {
      // the two ranges do not overlap or we are at the end of one of the
      // ranges
      if (r2 == other.ranges.end() ||
          (r1 != ranges.end() && r1->end < (r2->begin + offset)))
        {
          new_ranges.push_back(*r1);
          ++r1;
        }
      else if (r1 == ranges.end() || (r2->end + offset) < r1->begin)
        {
          new_ranges.emplace_back(r2->begin + offset, r2->end + offset);
          ++r2;
        }
      else
        {
          // ok, we do overlap, so just take the combination of the current
          // range (do not bother to merge with subsequent ranges)
          Range next(std::min(r1->begin, r2->begin + offset),
                     std::max(r1->end, r2->end + offset));
          new_ranges.push_back(next);
          ++r1;
          ++r2;
        }
    }
  ranges.swap(new_ranges);

  is_compressed = false;
  compress();
}



bool
IndexSet::is_subset_of(const IndexSet &other) const
{
  Assert(size() == other.size(),
         ExcMessage("One index set can only be a subset of another if they "
                    "describe index spaces of the same size. The ones in "
                    "question here have sizes " +
                    std::to_string(size()) + " and " +
                    std::to_string(other.size()) + "."));

  // See whether there are indices in the current set that are not in 'other'.
  // If so, then this is clearly not a subset of 'other'.
  IndexSet A_minus_B = *this;
  A_minus_B.subtract_set(other);
  if (A_minus_B.n_elements() > 0)
    return false;
  else
    // Else, every index in 'this' is also in 'other', since we ended up
    // with an empty set upon subtraction. This means that we have a subset:
    return true;
}



void
IndexSet::write(std::ostream &out) const
{
  compress();
  out << size() << " ";
  out << ranges.size() << std::endl;
  std::vector<Range>::const_iterator r = ranges.begin();
  for (; r != ranges.end(); ++r)
    {
      out << r->begin << " " << r->end << std::endl;
    }
}



void
IndexSet::read(std::istream &in)
{
  AssertThrow(in.fail() == false, ExcIO());

  size_type    s;
  unsigned int n_ranges;

  in >> s >> n_ranges;
  ranges.clear();
  set_size(s);
  for (unsigned int i = 0; i < n_ranges; ++i)
    {
      AssertThrow(in.fail() == false, ExcIO());

      size_type b, e;
      in >> b >> e;
      add_range(b, e);
    }
}


void
IndexSet::block_write(std::ostream &out) const
{
  AssertThrow(out.fail() == false, ExcIO());
  out.write(reinterpret_cast<const char *>(&index_space_size),
            sizeof(index_space_size));
  std::size_t n_ranges = ranges.size();
  out.write(reinterpret_cast<const char *>(&n_ranges), sizeof(n_ranges));
  if (ranges.empty() == false)
    out.write(reinterpret_cast<const char *>(&*ranges.begin()),
              ranges.size() * sizeof(Range));
  AssertThrow(out.fail() == false, ExcIO());
}

void
IndexSet::block_read(std::istream &in)
{
  size_type   size;
  std::size_t n_ranges;
  in.read(reinterpret_cast<char *>(&size), sizeof(size));
  in.read(reinterpret_cast<char *>(&n_ranges), sizeof(n_ranges));
  // we have to clear ranges first
  ranges.clear();
  set_size(size);
  ranges.resize(n_ranges, Range(0, 0));
  if (n_ranges != 0u)
    in.read(reinterpret_cast<char *>(&*ranges.begin()),
            ranges.size() * sizeof(Range));

  do_compress(); // needed so that largest_range can be recomputed
}



bool
IndexSet::is_element_binary_search(const size_type index) const
{
  // get the element after which we would have to insert a range that
  // consists of all elements from this element to the end of the index
  // range plus one. after this call we know that if p!=end() then
  // p->begin<=index unless there is no such range at all
  //
  // if the searched for element is an element of this range, then we're
  // done. otherwise, the element can't be in one of the following ranges
  // because otherwise p would be a different iterator
  //
  // since we already know the position relative to the largest range (we
  // called compress!), we can perform the binary search on ranges with
  // lower/higher number compared to the largest range
  std::vector<Range>::const_iterator p = std::upper_bound(
    ranges.begin() +
      (index < ranges[largest_range].begin ? 0 : largest_range + 1),
    index < ranges[largest_range].begin ? ranges.begin() + largest_range :
                                          ranges.end(),
    Range(index, size() + 1));

  if (p == ranges.begin())
    return ((index >= p->begin) && (index < p->end));

  Assert((p == ranges.end()) || (p->begin > index), ExcInternalError());

  // now move to that previous range
  --p;
  Assert(p->begin <= index, ExcInternalError());

  return (p->end > index);
}



IndexSet::size_type
IndexSet::nth_index_in_set_binary_search(const size_type n) const
{
  // find out which chunk the local index n belongs to by using a binary
  // search. the comparator is based on the end of the ranges.
  Range r(n, n + 1);
  r.nth_index_in_set = n;

  const std::vector<Range>::const_iterator p = Utilities::lower_bound(
    ranges.begin(), ranges.end(), r, Range::nth_index_compare);

  Assert(p != ranges.end(), ExcInternalError());
  return p->begin + (n - p->nth_index_in_set);
}



IndexSet::size_type
IndexSet::index_within_set_binary_search(const size_type n) const
{
  // we could try to use the main range for splitting up the search range, but
  // since we only come here when the largest range did not contain the index,
  // there is little gain from doing a first step manually.
  Range                              r(n, n);
  std::vector<Range>::const_iterator p =
    Utilities::lower_bound(ranges.begin(), ranges.end(), r, Range::end_compare);

  // if n is not in this set
  if (p == ranges.end() || p->end == n || p->begin > n)
    return numbers::invalid_dof_index;

  Assert(p != ranges.end(), ExcInternalError());
  Assert(p->begin <= n, ExcInternalError());
  Assert(n < p->end, ExcInternalError());
  return (n - p->begin) + p->nth_index_in_set;
}



IndexSet::ElementIterator
IndexSet::at(const size_type global_index) const
{
  compress();
  AssertIndexRange(global_index, size());

  if (ranges.empty())
    return end();

  std::vector<Range>::const_iterator main_range =
    ranges.begin() + largest_range;

  Range r(global_index, global_index + 1);
  // This optimization makes the bounds for lower_bound smaller by checking
  // the largest range first.
  std::vector<Range>::const_iterator range_begin, range_end;
  if (global_index < main_range->begin)
    {
      range_begin = ranges.begin();
      range_end   = main_range;
    }
  else
    {
      range_begin = main_range;
      range_end   = ranges.end();
    }

  // This will give us the first range p=[a,b[ with b>=global_index using
  // a binary search
  const std::vector<Range>::const_iterator p =
    Utilities::lower_bound(range_begin, range_end, r, Range::end_compare);

  // We couldn't find a range, which means we have no range that contains
  // global_index and also no range behind it, meaning we need to return end().
  if (p == ranges.end())
    return end();

  // Finally, we can have two cases: Either global_index is not in [a,b[,
  // which means we need to return an iterator to a because global_index, ...,
  // a-1 is not in the IndexSet (if branch). Alternatively, global_index is in
  // [a,b[ and we will return an iterator pointing directly at global_index
  // (else branch).
  if (global_index < p->begin)
    return {this, static_cast<size_type>(p - ranges.begin()), p->begin};
  else
    return {this, static_cast<size_type>(p - ranges.begin()), global_index};
}



std::vector<IndexSet::size_type>
IndexSet::get_index_vector() const
{
  compress();

  std::vector<size_type> indices;
  indices.reserve(n_elements());

  for (const auto &range : ranges)
    for (size_type entry = range.begin; entry < range.end; ++entry)
      indices.push_back(entry);

  Assert(indices.size() == n_elements(), ExcInternalError());

  return indices;
}



void
IndexSet::fill_index_vector(std::vector<size_type> &indices) const
{
  indices = get_index_vector();
}



#ifdef DEAL_II_WITH_TRILINOS
#  ifdef DEAL_II_TRILINOS_WITH_TPETRA

template <typename NodeType>
Tpetra::Map<int, types::signed_global_dof_index, NodeType>
IndexSet::make_tpetra_map(const MPI_Comm communicator,
                          const bool     overlapping) const
{
  return *make_tpetra_map_rcp<NodeType>(communicator, overlapping);
}



template <typename NodeType>
Teuchos::RCP<Tpetra::Map<int, types::signed_global_dof_index, NodeType>>
IndexSet::make_tpetra_map_rcp(const MPI_Comm communicator,
                              const bool     overlapping) const
{
  compress();
  (void)communicator;

  if constexpr (running_in_debug_mode())
    {
      if (!overlapping)
        {
          const size_type n_global_elements =
            Utilities::MPI::sum(n_elements(), communicator);
          Assert(n_global_elements == size(),
                 ExcMessage("You are trying to create an Tpetra::Map object "
                            "that partitions elements of an index set "
                            "between processors. However, the union of the "
                            "index sets on different processors does not "
                            "contain all indices exactly once: the sum of "
                            "the number of entries the various processors "
                            "want to store locally is " +
                            std::to_string(n_global_elements) +
                            " whereas the total size of the object to be "
                            "allocated is " +
                            std::to_string(size()) +
                            ". In other words, there are "
                            "either indices that are not spoken for "
                            "by any processor, or there are indices that are "
                            "claimed by multiple processors."));
        }
    }

  // Find out if the IndexSet is ascending and 1:1. This corresponds to a
  // linear Tpetra::Map. Overlapping IndexSets are never 1:1.
  const bool linear =
    overlapping ? false : is_ascending_and_one_to_one(communicator);
  if (linear)
    return Utilities::Trilinos::internal::make_rcp<
      Tpetra::Map<int, types::signed_global_dof_index, NodeType>>(
      size(),
      n_elements(),
      0,
#    ifdef DEAL_II_WITH_MPI
      Utilities::Trilinos::internal::make_rcp<Teuchos::MpiComm<int>>(
        communicator)
#    else
      Utilities::Trilinos::internal::make_rcp<Teuchos::Comm<int>>()
#    endif // DEAL_II_WITH_MPI
    );
  else
    {
      const std::vector<size_type>                indices = get_index_vector();
      std::vector<types::signed_global_dof_index> int_indices(indices.size());
      std::copy(indices.begin(), indices.end(), int_indices.begin());
      const Teuchos::ArrayView<types::signed_global_dof_index> arr_view(
        int_indices);

      return Utilities::Trilinos::internal::make_rcp<
        Tpetra::Map<int, types::signed_global_dof_index, NodeType>>(
        size(),
        arr_view,
        0,
#    ifdef DEAL_II_WITH_MPI
        Utilities::Trilinos::internal::make_rcp<Teuchos::MpiComm<int>>(
          communicator)
#    else
        Utilities::Trilinos::internal::make_rcp<Teuchos::Comm<int>>()
#    endif // DEAL_II_WITH_MPI
      );
    }
}
#  endif



Epetra_Map
IndexSet::make_trilinos_map(const MPI_Comm communicator,
                            const bool     overlapping) const
{
  compress();
  (void)communicator;

  if constexpr (running_in_debug_mode())
    {
      if (!overlapping)
        {
          const size_type n_global_elements =
            Utilities::MPI::sum(n_elements(), communicator);
          Assert(n_global_elements == size(),
                 ExcMessage("You are trying to create an Epetra_Map object "
                            "that partitions elements of an index set "
                            "between processors. However, the union of the "
                            "index sets on different processors does not "
                            "contain all indices exactly once: the sum of "
                            "the number of entries the various processors "
                            "want to store locally is " +
                            std::to_string(n_global_elements) +
                            " whereas the total size of the object to be "
                            "allocated is " +
                            std::to_string(size()) +
                            ". In other words, there are "
                            "either indices that are not spoken for "
                            "by any processor, or there are indices that are "
                            "claimed by multiple processors."));
        }
    }

  // Find out if the IndexSet is ascending and 1:1. This corresponds to a
  // linear EpetraMap. Overlapping IndexSets are never 1:1.
  const bool linear =
    overlapping ? false : is_ascending_and_one_to_one(communicator);

  if (linear)
    return Epetra_Map(TrilinosWrappers::types::int_type(size()),
                      TrilinosWrappers::types::int_type(n_elements()),
                      0,
#  ifdef DEAL_II_WITH_MPI
                      Epetra_MpiComm(communicator)
#  else
                      Epetra_SerialComm()
#  endif
    );
  else
    {
      const std::vector<size_type> indices = get_index_vector();
      return Epetra_Map(
        TrilinosWrappers::types::int_type(-1),
        TrilinosWrappers::types::int_type(n_elements()),
        (n_elements() > 0 ?
           reinterpret_cast<const TrilinosWrappers::types::int_type *>(
             indices.data()) :
           nullptr),
        0,
#  ifdef DEAL_II_WITH_MPI
        Epetra_MpiComm(communicator)
#  else
        Epetra_SerialComm()
#  endif
      );
    }
}
#endif


#ifdef DEAL_II_WITH_PETSC
IS
IndexSet::make_petsc_is(const MPI_Comm communicator) const
{
  std::vector<size_type> indices;
  fill_index_vector(indices);

  // If the size of the index set can be converted to a PetscInt then every
  // value can also be converted
  AssertThrowIntegerConversion(static_cast<PetscInt>(size()), size());
  const auto local_size = static_cast<PetscInt>(n_elements());
  AssertIntegerConversion(local_size, n_elements());

  size_type             i = 0;
  std::vector<PetscInt> petsc_indices(n_elements());
  for (const auto &index : *this)
    {
      const auto petsc_index = static_cast<PetscInt>(index);
      AssertIntegerConversion(petsc_index, index);
      petsc_indices[i] = petsc_index;
      ++i;
    }

  IS             is;
  PetscErrorCode ierr = ISCreateGeneral(
    communicator, local_size, petsc_indices.data(), PETSC_COPY_VALUES, &is);
  AssertThrow(ierr == 0, ExcPETScError(ierr));

  return is;
}
#endif



bool
IndexSet::is_ascending_and_one_to_one(const MPI_Comm communicator) const
{
  // If the sum of local elements does not add up to the total size,
  // the IndexSet can't be complete.
  const size_type n_global_elements =
    Utilities::MPI::sum(n_elements(), communicator);
  if (n_global_elements != size())
    return false;

  if (n_global_elements == 0)
    return true;

#ifdef DEAL_II_WITH_MPI
  // Non-contiguous IndexSets can't be linear.
  const bool all_contiguous =
    (Utilities::MPI::min(is_contiguous() ? 1 : 0, communicator) == 1);
  if (!all_contiguous)
    return false;

  bool is_globally_ascending = true;
  // we know that there is only one interval
  types::global_dof_index first_local_dof = (n_elements() > 0) ?
                                              *(begin_intervals()->begin()) :
                                              numbers::invalid_dof_index;

  const unsigned int my_rank = Utilities::MPI::this_mpi_process(communicator);
  const std::vector<types::global_dof_index> global_dofs =
    Utilities::MPI::gather(communicator, first_local_dof, 0);

  if (my_rank == 0)
    {
      // find out if the received std::vector is ascending
      types::global_dof_index index = 0;
      while (global_dofs[index] == numbers::invalid_dof_index)
        ++index;
      types::global_dof_index old_dof = global_dofs[index++];
      for (; index < global_dofs.size(); ++index)
        {
          const types::global_dof_index new_dof = global_dofs[index];
          if (new_dof != numbers::invalid_dof_index)
            {
              if (new_dof <= old_dof)
                {
                  is_globally_ascending = false;
                  break;
                }
              else
                old_dof = new_dof;
            }
        }
    }

  // now broadcast the result
  int is_ascending = is_globally_ascending ? 1 : 0;
  int ierr         = MPI_Bcast(&is_ascending, 1, MPI_INT, 0, communicator);
  AssertThrowMPI(ierr);

  return (is_ascending == 1);
#else
  return true;
#endif // DEAL_II_WITH_MPI
}



std::size_t
IndexSet::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(ranges) +
          MemoryConsumption::memory_consumption(is_compressed) +
          MemoryConsumption::memory_consumption(index_space_size) +
          sizeof(compress_mutex));
}

// explicit template instantiations

#ifndef DOXYGEN
#  ifdef DEAL_II_WITH_TRILINOS
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA

template IndexSet::IndexSet(
  const Teuchos::RCP<const Tpetra::Map<
    int,
    types::signed_global_dof_index,
    LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Host>>>
    &);

#      if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
        defined(KOKKOS_ENABLE_SYCL)
template IndexSet::IndexSet(
  const Teuchos::RCP<const Tpetra::Map<
    int,
    types::signed_global_dof_index,
    LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Default>>>
    &);
#      endif

template LinearAlgebra::TpetraWrappers::TpetraTypes::MapType<MemorySpace::Host>
dealii::IndexSet::make_tpetra_map<
  LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Host>>(
  const MPI_Comm,
  bool) const;

#      if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
        defined(KOKKOS_ENABLE_SYCL)
template LinearAlgebra::TpetraWrappers::TpetraTypes::MapType<
  MemorySpace::Default>
dealii::IndexSet::make_tpetra_map<
  LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Default>>(
  const MPI_Comm,
  bool) const;
#      endif

template Teuchos::RCP<
  LinearAlgebra::TpetraWrappers::TpetraTypes::MapType<MemorySpace::Host>>
dealii::IndexSet::make_tpetra_map_rcp<
  LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Host>>(
  const MPI_Comm,
  bool) const;

#      if defined(KOKKOS_ENABLE_CUDA) || defined(KOKKOS_ENABLE_HIP) || \
        defined(KOKKOS_ENABLE_SYCL)
template Teuchos::RCP<
  LinearAlgebra::TpetraWrappers::TpetraTypes::MapType<MemorySpace::Default>>
dealii::IndexSet::make_tpetra_map_rcp<
  LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<MemorySpace::Default>>(
  const MPI_Comm,
  bool) const;
#      endif

#    endif
#  endif
#endif

DEAL_II_NAMESPACE_CLOSE
