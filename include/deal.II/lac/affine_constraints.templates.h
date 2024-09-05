// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_affine_constraints_templates_h
#define dealii_affine_constraints_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi_compute_index_owner_internal.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_local_storage.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/trilinos_utilities.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/matrix_block.h>
#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_block_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_block_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_sparse_matrix.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>
#include <deal.II/lac/trilinos_vector.h>

#include <boost/serialization/complex.hpp>
#include <boost/serialization/utility.hpp>

#include <algorithm>
#include <complex>
#include <iomanip>
#include <numeric>
#include <ostream>
#include <set>

DEAL_II_NAMESPACE_OPEN



template <typename number>
std::size_t
AffineConstraints<number>::ConstraintLine::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(index) +
          MemoryConsumption::memory_consumption(entries) +
          MemoryConsumption::memory_consumption(inhomogeneity));
}



template <typename number>
void
AffineConstraints<number>::add_constraint(
  const size_type                                      constrained_dof,
  const ArrayView<const std::pair<size_type, number>> &dependencies,
  const number                                         inhomogeneity)
{
  Assert(sorted == false, ExcMatrixIsClosed());
  Assert(is_constrained(constrained_dof) == false,
         ExcMessage("You cannot add a constraint for a degree of freedom "
                    "that is already constrained."));

  // In debug mode, make sure that columns don't show up more than
  // once. Do this by sorting an array of the column indices. Try to
  // avoid memory allocation by using a vector type that allocates up
  // to 25 entries on the stack -- this is enough for the constraints
  // of Q4 elements in 3d, and so should cover the vast majority of
  // cases. If we have a constraint with more dependencies, then
  // that's just going to require a heap allocation.
#ifdef DEBUG
  {
    boost::container::small_vector<size_type, 25> column_indices;
    column_indices.reserve(dependencies.size());
    for (const auto &d : dependencies)
      column_indices.emplace_back(d.first);
    std::sort(column_indices.begin(), column_indices.end());
    Assert(std::adjacent_find(column_indices.begin(), column_indices.end()) ==
             column_indices.end(),
           ExcMessage(
             "You are trying to insert a constraint that lists the same "
             "degree of freedom more than once on the right hand side. This is "
             "not allowed."));
  }
#endif


  // The following can happen when we compute with distributed meshes and dof
  // handlers and we constrain a degree of freedom whose number we don't have
  // locally. if we don't abort here the program will try to allocate several
  // terabytes of memory to resize the various arrays below :-)
  Assert(constrained_dof != numbers::invalid_size_type, ExcInternalError());

  // if necessary enlarge vector of existing entries for cache
  const size_type line_index = calculate_line_index(constrained_dof);
  if (line_index >= lines_cache.size())
    lines_cache.resize(std::max(2 * static_cast<size_type>(lines_cache.size()),
                                line_index + 1),
                       numbers::invalid_size_type);

  // Push a new line to the end of the list and fill it with the
  // provided information:
  ConstraintLine &constraint = lines.emplace_back();
  constraint.index           = constrained_dof;
  constraint.entries.reserve(dependencies.size());
  for (const auto &[column, weight] : dependencies)
    {
      Assert(column != constrained_dof,
             ExcMessage("You cannot constrain a degree of freedom against "
                        "itself."));
      constraint.entries.emplace_back(column, weight);
    }
  constraint.inhomogeneity = inhomogeneity;

  // Record the new constraint in the cache:
  Assert(lines_cache[line_index] == numbers::invalid_size_type,
         ExcInternalError());
  lines_cache[line_index] = lines.size() - 1;
}



template <typename number>
void
AffineConstraints<number>::constrain_dof_to_zero(
  const size_type constrained_dof)
{
  Assert(sorted == false, ExcMatrixIsClosed());

  // The following can happen when we compute with distributed meshes and dof
  // handlers and we constrain a degree of freedom whose number we don't have
  // locally. if we don't abort here the program will try to allocate several
  // terabytes of memory to resize the various arrays below :-)
  Assert(constrained_dof != numbers::invalid_size_type, ExcInternalError());

  // if necessary enlarge vector of existing entries for cache
  const size_type line_index = calculate_line_index(constrained_dof);
  if (line_index >= lines_cache.size())
    lines_cache.resize(std::max(2 * static_cast<size_type>(lines_cache.size()),
                                line_index + 1),
                       numbers::invalid_size_type);

  // Let's check whether the DoF is already constrained. This is only allowed
  // if it had previously been constrained to zero, and only then.
  if (lines_cache[line_index] != numbers::invalid_size_type)
    {
      Assert(lines[lines_cache[line_index]].entries.empty() &&
               (lines[lines_cache[line_index]].inhomogeneity == number(0.)),
             ExcMessage("You cannot constrain a degree of freedom "
                        "to zero that is is already constrained to "
                        "something else."));
    }
  else
    {
      // Push a new line to the end of the list and fill it with the
      // provided information:
      ConstraintLine &constraint = lines.emplace_back();
      constraint.index           = constrained_dof;
      constraint.inhomogeneity   = 0.;

      // Record the new constraint in the cache:
      Assert(lines_cache[line_index] == numbers::invalid_size_type,
             ExcInternalError());
      lines_cache[line_index] = lines.size() - 1;
    }
}



template <typename number>
typename AffineConstraints<number>::LineRange
AffineConstraints<number>::get_lines() const
{
  return boost::make_iterator_range(lines.begin(), lines.end());
}



template <typename number>
bool
AffineConstraints<number>::is_consistent_in_parallel(
  const std::vector<IndexSet> &locally_owned_dofs,
  const IndexSet              &locally_active_dofs,
  const MPI_Comm               mpi_communicator,
  const bool                   verbose) const
{
  // Helper to return a ConstraintLine object that belongs to row @p row.
  // If @p row is not constrained or not stored locally, return an empty
  // constraint object that would correspond to a zero constraints
  auto get_line = [&](const size_type line_n) -> ConstraintLine {
    const size_type line_index = calculate_line_index(line_n);
    if (line_index >= lines_cache.size() ||
        lines_cache[line_index] == numbers::invalid_size_type)
      {
        const ConstraintLine empty = {line_n, {}, 0.0};
        return empty;
      }
    else
      return lines[lines_cache[line_index]];
  };

  // identify non-owned rows and send to owner:
  std::map<unsigned int, std::vector<ConstraintLine>> to_send;

  const unsigned int myid =
    dealii::Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int nproc =
    dealii::Utilities::MPI::n_mpi_processes(mpi_communicator);

  // We will send all locally active dofs that are not locally owned for
  // checking. Note that we allow constraints to differ on locally_relevant (and
  // not active) DoFs.
  // NOLINTNEXTLINE(performance-unnecessary-copy-initialization)
  IndexSet non_owned = locally_active_dofs;
  non_owned.subtract_set(locally_owned_dofs[myid]);
  for (unsigned int owner = 0; owner < nproc; ++owner)
    {
      // find all lines to send to @p owner
      const IndexSet indices_to_send = non_owned & locally_owned_dofs[owner];
      for (const auto row_idx : indices_to_send)
        {
          to_send[owner].emplace_back(get_line(row_idx));
        }
    }

  const std::map<unsigned int, std::vector<ConstraintLine>> received =
    Utilities::MPI::some_to_some(mpi_communicator, to_send);

  unsigned int inconsistent = 0;

  // from each processor:
  for (const auto &kv : received)
    {
      // for each incoming line:
      for (const auto &lineit : kv.second)
        {
          const ConstraintLine reference = get_line(lineit.index);

          if (lineit.inhomogeneity != reference.inhomogeneity)
            {
              ++inconsistent;

              if (verbose)
                std::cout << "Proc " << myid << " got line " << lineit.index
                          << " from " << kv.first << " inhomogeneity "
                          << lineit.inhomogeneity
                          << " != " << reference.inhomogeneity << std::endl;
            }
          else if (lineit.entries != reference.entries)
            {
              ++inconsistent;
              if (verbose)
                std::cout << "Proc " << myid << " got line " << lineit.index
                          << " from " << kv.first << " with wrong values!"
                          << std::endl;
            }
        }
    }

  const unsigned int total =
    Utilities::MPI::sum(inconsistent, mpi_communicator);
  if (verbose && total > 0 && myid == 0)
    std::cout << total << " inconsistent lines discovered!" << std::endl;
  return total == 0;
}


namespace internal
{
  template <typename number>
  std::vector<typename dealii::AffineConstraints<number>::ConstraintLine>
  compute_locally_relevant_constraints(
    const dealii::AffineConstraints<number> &constraints_in,
    const IndexSet                          &locally_owned_dofs,
    const IndexSet                          &locally_relevant_dofs,
    const MPI_Comm                           mpi_communicator)
  {
    // The result vector filled step by step.
    using ConstraintLine =
      typename dealii::AffineConstraints<number>::ConstraintLine;
    std::vector<ConstraintLine> locally_relevant_constraints;

#ifndef DEAL_II_WITH_MPI
    AssertThrow(false, ExcNotImplemented()); // one should not come here
    (void)constraints_in;
    (void)locally_owned_dofs;
    (void)locally_relevant_dofs;
    (void)mpi_communicator;
#else

    const unsigned int my_rank =
      Utilities::MPI::this_mpi_process(mpi_communicator);

    // First define a helper function that sorts and normalizes the constraints
    // provided through the function's argument:
    const auto sort_and_make_unique =
      [](std::vector<ConstraintLine> &constraints) {
        if (constraints.empty())
          return;

        // First sort the array of constraints by their index:
        std::sort(constraints.begin(),
                  constraints.end(),
                  [](const ConstraintLine &l1, const ConstraintLine &l2) {
                    return l1.index < l2.index;
                  });

        // It is possible that two processes have computed constraints
        // for the same DoF differently (for example, in parallel because
        // they have visited different faces of the same cell because
        // on different processes, different neighbor cells are artificial.
        //
        // This is ok, as long as after resolution of all chains of
        // constraints we end up with the same constraint. But at this
        // point, we just don't know yet. We deal with this by simply
        // dropping constraints for DoFs for which there is a previous
        // constraint in the list:
        constraints.erase(std::unique(constraints.begin(),
                                      constraints.end(),
                                      [](const ConstraintLine &a,
                                         const ConstraintLine &b) {
                                        return (a.index == b.index);
                                      }),
                          constraints.end());

        // For those constraints that survive, sort the right hand side arrays:
        for (ConstraintLine &entry : constraints)
          std::sort(entry.entries.begin(),
                    entry.entries.end(),
                    [](const auto &l1, const auto &l2) {
                      return l1.first < l2.first;
                    });
      };

    // step 0: Collect the indices of constrained DoFs we know of:
    IndexSet my_constraint_indices(locally_owned_dofs.size());
    {
      std::vector<types::global_dof_index> indices;
      indices.reserve(constraints_in.n_constraints());
      for (const ConstraintLine &line : constraints_in.get_lines())
        if (locally_owned_dofs.is_element(line.index) == false)
          indices.push_back(line.index);
      std::sort(indices.begin(), indices.end());
      my_constraint_indices.add_indices(indices.begin(), indices.end());
    }

    // step 1: Identify the owners of DoFs we know to be constrained but do
    //         not own.
    std::vector<unsigned int> owners_of_my_constraints(
      my_constraint_indices.n_elements());
    Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
      constrained_indices_process(locally_owned_dofs,
                                  my_constraint_indices,
                                  mpi_communicator,
                                  owners_of_my_constraints,
                                  true);

    Utilities::MPI::ConsensusAlgorithms::Selector<
      std::vector<std::pair<types::global_dof_index, types::global_dof_index>>,
      std::vector<unsigned int>>
      consensus_algorithm;
    consensus_algorithm.run(constrained_indices_process, mpi_communicator);

    // step 2: Collect all locally owned constraints into a data structure
    //         that we can send to other processes that want to know
    //         about them.
    const std::map<unsigned int, IndexSet> constrained_indices_by_ranks =
      constrained_indices_process.get_requesters();
    {
      std::map<unsigned int, std::vector<ConstraintLine>> send_data;

      // For each constraint we know of but owned by another process, create a
      // copy of the constraint and add it to the list of things to send to
      // other processes:
      for (unsigned int constraint_index = 0;
           constraint_index < owners_of_my_constraints.size();
           ++constraint_index)
        {
          ConstraintLine entry;

          const types::global_dof_index index =
            my_constraint_indices.nth_index_in_set(constraint_index);

          entry.index = index;
          if (const std::vector<std::pair<types::global_dof_index, number>>
                *constraints = constraints_in.get_constraint_entries(index))
            entry.entries = *constraints;
          entry.inhomogeneity = constraints_in.get_inhomogeneity(index);


          Assert(owners_of_my_constraints[constraint_index] != my_rank,
                 ExcInternalError());
          send_data[owners_of_my_constraints[constraint_index]].push_back(
            entry);
        }

      // Now exchange this data between processes:
      const std::map<unsigned int, std::vector<ConstraintLine>> received_data =
        Utilities::MPI::some_to_some(mpi_communicator, send_data);

      // Finally join things with the constraints we know about and own
      // ourselves, collate everything we received into one array, sort,
      // and make it unique:
      for (const ConstraintLine &line : constraints_in.get_lines())
        if (locally_owned_dofs.is_element(line.index))
          locally_relevant_constraints.push_back(line);
      for (const auto &[rank, constraints] : received_data)
        locally_relevant_constraints.insert(locally_relevant_constraints.end(),
                                            constraints.begin(),
                                            constraints.end());

      sort_and_make_unique(locally_relevant_constraints);
    }

    // step 3: communicate constraints so that each process knows how the
    // locally relevant dofs are constrained
    {
      // ... determine owners of locally relevant dofs
      IndexSet locally_relevant_dofs_non_local = locally_relevant_dofs;
      locally_relevant_dofs_non_local.subtract_set(locally_owned_dofs);

      std::vector<unsigned int> locally_relevant_dofs_owners(
        locally_relevant_dofs_non_local.n_elements());
      Utilities::MPI::internal::ComputeIndexOwner::ConsensusAlgorithmsPayload
        locally_relevant_dofs_process(locally_owned_dofs,
                                      locally_relevant_dofs_non_local,
                                      mpi_communicator,
                                      locally_relevant_dofs_owners,
                                      true);

      Utilities::MPI::ConsensusAlgorithms::Selector<
        std::vector<
          std::pair<types::global_dof_index, types::global_dof_index>>,
        std::vector<unsigned int>>
        consensus_algorithm;
      consensus_algorithm.run(locally_relevant_dofs_process, mpi_communicator);

      const auto locally_relevant_dofs_by_ranks =
        locally_relevant_dofs_process.get_requesters();

      std::map<unsigned int, std::vector<ConstraintLine>> send_data;

      for (const std::map<unsigned int, IndexSet>::value_type
             &rank_and_indices : locally_relevant_dofs_by_ranks)
        {
          Assert(rank_and_indices.first != my_rank, ExcInternalError());

          std::vector<ConstraintLine> data;

          for (const types::global_dof_index index : rank_and_indices.second)
            {
              // note: at this stage locally_relevant_constraints still
              // contains only locally owned constraints
              const typename std::vector<ConstraintLine>::iterator ptr =
                std::find_if(locally_relevant_constraints.begin(),
                             locally_relevant_constraints.end(),
                             [index](const auto &a) {
                               return a.index == index;
                             });
              if (ptr != locally_relevant_constraints.end())
                data.push_back(*ptr);
            }

          send_data[rank_and_indices.first] = std::move(data);
        }

      const std::map<unsigned int, std::vector<ConstraintLine>> received_data =
        Utilities::MPI::some_to_some(mpi_communicator, send_data);

      // Finally collate everything into one array, sort, and make it unique:
      for (const auto &[rank, constraints] : received_data)
        locally_relevant_constraints.insert(locally_relevant_constraints.end(),
                                            constraints.begin(),
                                            constraints.end());

      sort_and_make_unique(locally_relevant_constraints);
    }

#endif

    return locally_relevant_constraints;
  }
} // namespace internal



template <typename number>
void
AffineConstraints<number>::make_consistent_in_parallel(
  const IndexSet &locally_owned_dofs,
  const IndexSet &constraints_to_make_consistent_,
  const MPI_Comm  mpi_communicator)
{
  // We need to resolve chains of constraints in this function, and that
  // is exactly what close() does. So call this function at the beginning
  // to make sure we start working with constraints that are already
  // resolved to the degree possible before exchanging information.
  // We do this here so that it is also done if we choose the early
  // exit in the following statement.
  if (sorted == false)
    close();

  if (Utilities::MPI::n_mpi_processes(mpi_communicator) == 1)
    return; // Nothing to do, since serial.

  Assert(this->local_lines.size() > 0,
         ExcMessage(
           "This functionality requires that the AffineConstraints object "
           "knows for which degrees of freedom it can store constraints. "
           "Please initialize this object with the corresponding index sets."));

  // Container for indices that are constrained or that other indices are
  // constrained against. Generously reserve memory for it, assuming the worst
  // case that each constraint line involves two distinct DoF indices.
  std::vector<types::global_dof_index> constrained_indices;
  constrained_indices.reserve(2 * n_constraints());

  // This IndexSet keeps track of the locally stored constraints on this
  // AffineConstraints object: If we have received constrained indices that we
  // don't know of, we need to expand our playing field.
  IndexSet locally_stored_constraints;

  // This IndexSet contains those DoFs about which we want to know all
  // constraints. We will receive constraints for these DoFs against other DoFs,
  // which might be constrained themselves. To successfully resolve all chains
  // of constraints, we need to know the constraints of all these DoFs, and we
  // keep track of the relevant DoFs here.
  IndexSet constraints_to_make_consistent = constraints_to_make_consistent_;

  // This IndexSet stores DoFs from the previous iteration. If the old and new
  // index sets match, we have converged.
  IndexSet constraints_made_consistent;

  const unsigned int max_iterations  = 10;
  unsigned int       iteration_count = 0;
  for (; iteration_count < max_iterations; ++iteration_count)
    {
      // 1) Get all locally relevant constraints we need to know about:
      const std::vector<ConstraintLine> imported_constraints =
        internal::compute_locally_relevant_constraints(
          *this,
          locally_owned_dofs,
          constraints_to_make_consistent,
          mpi_communicator);

      // 2) Add untracked DoFs to the index sets.
      constrained_indices.clear();
      for (const auto &line : imported_constraints)
        {
          constrained_indices.push_back(line.index);
          for (const auto &entry : line.entries)
            constrained_indices.push_back(entry.first);
        }
      std::sort(constrained_indices.begin(), constrained_indices.end());

      locally_stored_constraints = this->local_lines;
      locally_stored_constraints.add_indices(constrained_indices.begin(),
                                             constrained_indices.end());

      constraints_made_consistent = constraints_to_make_consistent;
      constraints_to_make_consistent.add_indices(constrained_indices.begin(),
                                                 constrained_indices.end());

      // 3) Clear and refill this constraint matrix. Also resolve chains:
      this->reinit(locally_owned_dofs, locally_stored_constraints);
      for (const auto &line : imported_constraints)
        this->add_constraint(line.index, line.entries, line.inhomogeneity);
      this->close();

      // 4) Stop loop if converged.
      const bool constraints_converged =
        (Utilities::MPI::min(
           (constraints_to_make_consistent == constraints_made_consistent ? 1 :
                                                                            0),
           mpi_communicator) == 1);
      if (constraints_converged)
        break;
    }

  AssertThrow(iteration_count < max_iterations,
              ExcMessage(
                "make_consistent_in_parallel() did not converge after " +
                Utilities::to_string(max_iterations) + " iterations."));

  Assert(this->is_consistent_in_parallel(
           Utilities::MPI::all_gather(mpi_communicator, locally_owned_dofs),
           constraints_to_make_consistent,
           mpi_communicator),
         ExcInternalError());
}



template <typename number>
void
AffineConstraints<number>::add_lines(const std::set<size_type> &lines)
{
  for (const size_type &line : lines)
    add_line(line);
}



template <typename number>
void
AffineConstraints<number>::add_lines(const std::vector<bool> &lines)
{
  for (size_type i = 0; i < lines.size(); ++i)
    if (lines[i] == true)
      add_line(i);
}



template <typename number>
void
AffineConstraints<number>::add_lines(const IndexSet &lines)
{
  for (size_type i = 0; i < lines.n_elements(); ++i)
    add_line(lines.nth_index_in_set(i));
}



template <typename number>
void
AffineConstraints<number>::add_entries(
  const size_type                                  constrained_dof_index,
  const std::vector<std::pair<size_type, number>> &col_weight_pairs)
{
  Assert(sorted == false, ExcMatrixIsClosed());
  Assert(is_constrained(constrained_dof_index),
         ExcLineInexistent(constrained_dof_index));

  ConstraintLine &line =
    lines[lines_cache[calculate_line_index(constrained_dof_index)]];
  Assert(line.index == constrained_dof_index, ExcInternalError());

  // if in debug mode, check whether an entry for this column already
  // exists and if its the same as the one entered at present
  //
  // in any case: skip this entry if an entry for this column already
  // exists, since we don't want to enter it twice
  line.entries.reserve(line.entries.size() + col_weight_pairs.size());
  for (const std::pair<size_type, number> &col_weight_pair : col_weight_pairs)
    {
      Assert(constrained_dof_index != col_weight_pair.first,
             ExcMessage("Can't constrain a degree of freedom to itself"));
      bool entry_exists = false;
      for (const std::pair<size_type, number> &entry : line.entries)
        if (entry.first == col_weight_pair.first)
          {
            // entry exists, break innermost loop
            Assert(entry.second == col_weight_pair.second,
                   ExcEntryAlreadyExists(constrained_dof_index,
                                         col_weight_pair.first,
                                         entry.second,
                                         col_weight_pair.second));
            entry_exists = true;
            break;
          }

      if (entry_exists == false)
        line.entries.push_back(col_weight_pair);
    }
}



template <typename number>
void
AffineConstraints<number>::add_selected_constraints(
  const AffineConstraints &constraints,
  const IndexSet          &filter)
{
  if (constraints.n_constraints() == 0)
    return;

  Assert(filter.size() > constraints.lines.back().index,
         ExcMessage(
           "The filter must be larger than the given constraints object."));
  for (const ConstraintLine &line : constraints.lines)
    if (filter.is_element(line.index))
      {
        const size_type row = filter.index_within_set(line.index);
        add_line(row);
        set_inhomogeneity(row, line.inhomogeneity);
        for (const std::pair<size_type, number> &entry : line.entries)
          if (filter.is_element(entry.first))
            add_entry(row, filter.index_within_set(entry.first), entry.second);
      }
}



template <typename number>
void
AffineConstraints<number>::close()
{
  if (sorted == true)
    return;

  // sort the lines
  std::sort(lines.begin(),
            lines.end(),
            [](const ConstraintLine &l1, const ConstraintLine &l2) {
              return l1.index < l2.index;
            });

  // update list of pointers and give the vector a sharp size since we
  // won't modify the size any more after this point.
  {
    std::vector<size_type> new_lines(lines_cache.size(),
                                     numbers::invalid_size_type);
    size_type              counter = 0;
    for (const ConstraintLine &line : lines)
      {
        new_lines[calculate_line_index(line.index)] = counter;
        ++counter;
      }
    lines_cache = std::move(new_lines);
  }

  // in debug mode: check whether we really set the pointers correctly.
  for (size_type i = 0; i < lines_cache.size(); ++i)
    if (lines_cache[i] != numbers::invalid_size_type)
      Assert(i == calculate_line_index(lines[lines_cache[i]].index),
             ExcInternalError());

  // The second part is that we need to work on the individual lines.
  // Let us start by stripping zero entries. That would mean that in the linear
  // constraint for a node, x_i = ax_1 + bx_2 + ..., another node times 0
  // appears. obviously, 0*something can be omitted.
  //
  // This can be done in parallel:
  parallel::apply_to_subranges(
    lines.begin(),
    lines.end(),
    [](const typename std::vector<ConstraintLine>::iterator begin,
       const typename std::vector<ConstraintLine>::iterator end) {
      for (ConstraintLine &line : boost::iterator_range<
             typename std::vector<ConstraintLine>::iterator>(begin, end))
        line.entries.erase(
          std::remove_if(line.entries.begin(),
                         line.entries.end(),
                         [](const std::pair<size_type, number> &p) {
                           return p.second == number(0.);
                         }),
          line.entries.end());
    },
    /* grainsize = */ 100);



#ifdef DEBUG
  // In debug mode we are computing an estimate for the maximum number
  // of constraints so that we can bail out if there is a cycle in the
  // constraints (which is easier than searching for cycles in the graph).
  //
  // Let us figure out the largest dof index. This is an upper bound for the
  // number of constraints because it is an approximation for the number of dofs
  // in our system.
  size_type largest_idx = 0;
  for (const ConstraintLine &line : lines)
    for (const std::pair<size_type, number> &entry : line.entries)
      largest_idx = std::max(largest_idx, entry.first);
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
#ifdef DEBUG
  size_type iteration = 0;
#endif
  bool              chained_constraint_replaced = false;
  std::vector<bool> line_finalized(lines.size(), false);
  do
    {
      chained_constraint_replaced = false;
      for (unsigned int line_index = 0; line_index < lines.size(); ++line_index)
        if (line_finalized[line_index] == false)
          {
            ConstraintLine &line = lines[line_index];

#ifdef DEBUG
            // we need to keep track of how many replacements we do in this
            // line, because we can end up in a cycle A->B->C->A without the
            // number of entries growing.
            size_type n_replacements = 0;
#endif

            // loop over all entries of this line (including ones that we
            // have appended in this go around) and see whether they are
            // further constrained. ignore elements that we don't store on
            // the current processor.
            const unsigned int n_original_entries  = line.entries.size();
            bool               has_sub_constraints = false;
            const size_type    lines_cache_size    = lines_cache.size();
            for (unsigned int entry = 0; entry < n_original_entries; ++entry)
              {
                const size_type dof_index =
                  calculate_line_index(line.entries[entry].first);

                if (dof_index < lines_cache_size &&
                    lines_cache[dof_index] != numbers::invalid_size_type)
                  {
                    // Now we found an index that is further constrained, so
                    // we need to resolve the entries.  We may replace some of
                    // them, and we may add some, but we only walk through the
                    // original entries (and also don't touch the
                    // replacements) even though in principle we could also do
                    // the ones we have added or replaced. There is no harm in
                    // doing so: we will simply treat those in the next round
                    // around of the outer iteration.
                    has_sub_constraints = true;

                    const number weight = line.entries[entry].second;
                    Assert(line.entries[entry].first != line.index,
                           ExcMessage("Cycle in constraints detected!"));

                    const ConstraintLine &constrained_line =
                      lines[lines_cache[dof_index]];
                    Assert(constrained_line.index == line.entries[entry].first,
                           ExcInternalError());

                    // now we have to replace an entry by its expansion. we do
                    // that by overwriting the entry by the first entry of the
                    // expansion and adding the remaining ones to the end,
                    // where we will later process them once more
                    //
                    // we can of course only do that if the DoF that we are
                    // currently handling is constrained by a linear combination
                    // of other dofs:
                    if (constrained_line.entries.size() > 0)
                      {
                        for (size_type i = 0;
                             i < constrained_line.entries.size();
                             ++i)
                          Assert(line.entries[entry].first !=
                                   constrained_line.entries[i].first,
                                 ExcMessage("Cycle in constraints detected!"));

                        // replace first entry, then tack the rest to the end
                        // of the list
                        line.entries[entry] = std::pair<size_type, number>(
                          constrained_line.entries[0].first,
                          constrained_line.entries[0].second * weight);

                        for (size_type i = 1;
                             i < constrained_line.entries.size();
                             ++i)
                          line.entries.emplace_back(
                            constrained_line.entries[i].first,
                            constrained_line.entries[i].second * weight);

#ifdef DEBUG
                        // keep track of how many entries we replace in this
                        // line. If we do more than there are constraints or
                        // dofs in our system, we must have a cycle.
                        ++n_replacements;
                        Assert(n_replacements / 2 < largest_idx,
                               ExcMessage("Cycle in constraints detected!"));
#endif
                      }
                    else
                      // the DoF that we encountered is not constrained by a
                      // linear combination of other dofs but is equal to just
                      // the inhomogeneity (i.e. its chain of entries is
                      // empty). in that case, we can't just overwrite the
                      // current entry, but we have to actually eliminate
                      // it. we do not want to change the loop length above we
                      // do so by setting the 'first' entry to
                      // invalid_size_type here to finally remove entries in a
                      // second loop
                      {
                        line.entries[entry].first = numbers::invalid_size_type;
                      }

                    line.inhomogeneity +=
                      constrained_line.inhomogeneity * weight;
                  }
              }

            // If none of the entries in the current line refer to DoFs that are
            // themselves constrained, then we can move on:
            if (has_sub_constraints == false)
              line_finalized[line_index] = true;
            else
              {
                chained_constraint_replaced = true;

                // Now delete the elements we have marked for deletion.
                auto remaining_entries = line.entries.begin();
                for (const auto &entry : line.entries)
                  if (entry.first != numbers::invalid_size_type)
                    {
                      *remaining_entries = entry;
                      ++remaining_entries;
                    }
                line.entries.erase(remaining_entries, line.entries.end());
              }
          }

#ifdef DEBUG
      // increase iteration count. note that we should not iterate more
      // times than there are constraints, since this puts a natural upper
      // bound on the length of constraint chains
      ++iteration;
      Assert(iteration <= lines.size() + 1, ExcInternalError());
#endif
    }
  while (chained_constraint_replaced == true);

  // Finally sort the entries and re-scale them if necessary. in this step,
  // we also throw out duplicates as mentioned above. moreover, as some
  // entries might have had zero weights, we replace them by a vector with
  // sharp sizes.
  //
  // This is again an operation that works on each line separately. It can be
  // run in parallel:
  parallel::apply_to_subranges(
    lines.begin(),
    lines.end(),
    [](const typename std::vector<ConstraintLine>::iterator &begin,
       const typename std::vector<ConstraintLine>::iterator &end) {
      for (ConstraintLine &line : boost::iterator_range<
             typename std::vector<ConstraintLine>::iterator>(begin, end))
        {
          unsigned int duplicates                   = 0;
          bool         is_sorted_without_duplicates = true;
          for (unsigned int i = 1; i < line.entries.size(); ++i)
            if (!(line.entries[i - 1].first < line.entries[i].first))
              {
                is_sorted_without_duplicates = false;
                break;
              }
          if (is_sorted_without_duplicates == false)
            {
              std::sort(line.entries.begin(),
                        line.entries.end(),
                        [](const std::pair<unsigned int, number> &a,
                           const std::pair<unsigned int, number> &b) -> bool {
                          // Just look at the index, ignore the value.
                          return a.first < b.first;
                        });

              // loop over the now sorted list and see whether any of the
              // entries references the same dofs more than once in order to
              // find how many non-duplicate entries we have. This lets us
              // allocate the correct amount of memory for the constraint
              // entries.
              for (size_type i = 1; i < line.entries.size(); ++i)
                if (line.entries[i].first == line.entries[i - 1].first)
                  ++duplicates;
            }

          if (duplicates > 0 || (line.entries.size() < line.entries.capacity()))
            {
              typename ConstraintLine::Entries new_entries;

              // if we have no duplicates, copy verbatim the entries. this way,
              // the final size is of the vector is correct.
              if (duplicates == 0)
                new_entries = line.entries;
              else
                {
                  // otherwise, we need to go through the list and resolve the
                  // duplicates
                  new_entries.reserve(line.entries.size() - duplicates);
                  new_entries.push_back(line.entries[0]);
                  for (size_type j = 1; j < line.entries.size(); ++j)
                    if (line.entries[j].first == line.entries[j - 1].first)
                      {
                        Assert(new_entries.back().first ==
                                 line.entries[j].first,
                               ExcInternalError());
                        new_entries.back().second += line.entries[j].second;
                      }
                    else
                      new_entries.push_back(line.entries[j]);

                  Assert(new_entries.size() == line.entries.size() - duplicates,
                         ExcInternalError());

                  // make sure there are really no duplicates left and that the
                  // list is still sorted
                  for (size_type j = 1; j < new_entries.size(); ++j)
                    {
                      Assert(new_entries[j].first != new_entries[j - 1].first,
                             ExcInternalError());
                      Assert(new_entries[j].first > new_entries[j - 1].first,
                             ExcInternalError());
                    }
                }

              // replace old list of constraints for this dof by the new one
              line.entries.swap(new_entries);
            }

          // Finally do the following check: if the sum of weights for the
          // constraints is close to one, but not exactly one, then rescale all
          // the weights so that they sum up to 1. this adds a little numerical
          // stability and avoids all sorts of problems where the actual value
          // is close to, but not quite what we expected
          //
          // the case where the weights don't quite sum up happens when we
          // compute the interpolation weights "on the fly", i.e. not from
          // precomputed tables. in this case, the interpolation weights are
          // also subject to round-off
          number sum = 0.;
          for (const std::pair<size_type, number> &entry : line.entries)
            sum += entry.second;
          if (std::abs(sum - number(1.)) < 1.e-13 &&
              std::abs(sum - number(1.)) > 0.)
            {
              const number inverse_sum = number(1.) / sum;
              for (std::pair<size_type, number> &entry : line.entries)
                entry.second *= inverse_sum;
              line.inhomogeneity *= inverse_sum;
            }
        }
    },
    /* grainsize = */ 100);

  // if in debug mode: check that no dof is constrained to another dof that
  // is also constrained. exclude dofs from this check whose constraint
  // lines are not stored on the local processor
  Assert(std::none_of(lines.begin(),
                      lines.end(),
                      [this](const ConstraintLine &line) {
                        for (const std::pair<size_type, number> &entry :
                             line.entries)
                          if ((local_lines.size() == 0) ||
                              (local_lines.is_element(entry.first)))
                            {
                              // make sure that entry->first is not the index of
                              // a line itself
                              const bool is_cycle = is_constrained(entry.first);
                              if (is_cycle)
                                return true;
                            }
                        return false;
                      }),
         ExcMessage("The constraints represented by this object have a cycle. "
                    "This is not allowed."));

  // Now also compute which indices we need in distribute().
  //
  // This processor owns only part of the vectors we will work on. One may think
  // that every processor should be able to simply communicate those elements it
  // owns and for which it knows that they act as sources to constrained DoFs to
  // the owner of these DoFs. This would lead to a scheme where all we need to
  // do is to add some local elements to (possibly non-local) ones and then call
  // compress().
  //
  // Alas, this scheme does not work as evidenced by the disaster of
  // bug #51 (originally stored in the Google Code repository, but now
  // unfortunately no longer available because that platform has gone
  // away) and reversion of one attempt that implements this in svn
  // r29662. Rather, we need to get a vector that has all the
  // *sources* or constraints we own locally, possibly as ghost vector
  // elements, then read from them, and finally throw away the ghosted
  // vector. Implement this in the following.
  needed_elements_for_distribute = locally_owned_dofs;

  if (needed_elements_for_distribute !=
      complete_index_set(locally_owned_dofs.size()))
    {
      std::vector<types::global_dof_index> additional_elements;
      for (const ConstraintLine &line : lines)
        if (locally_owned_dofs.is_element(line.index))
          for (const std::pair<size_type, number> &entry : line.entries)
            if (!locally_owned_dofs.is_element(entry.first))
              additional_elements.emplace_back(entry.first);
      std::sort(additional_elements.begin(), additional_elements.end());
      needed_elements_for_distribute.add_indices(additional_elements.begin(),
                                                 additional_elements.end());
    }

  sorted = true;
}



template <typename number>
bool
AffineConstraints<number>::is_closed() const
{
  return sorted || (n_constraints() == 0);
}



template <typename number>
bool
AffineConstraints<number>::is_closed(const MPI_Comm comm) const
{
  return Utilities::MPI::min(static_cast<unsigned int>(is_closed()), comm) == 1;
}



template <typename number>
void
AffineConstraints<number>::shift(const size_type offset)
{
  if (local_lines.size() == 0)
    lines_cache.insert(lines_cache.begin(), offset, numbers::invalid_size_type);
  else
    {
      // shift local_lines
      IndexSet new_local_lines(local_lines.size() + offset);
      new_local_lines.add_indices(local_lines, offset);
      std::swap(local_lines, new_local_lines);
    }

  for (ConstraintLine &line : lines)
    {
      line.index += offset;
      for (std::pair<size_type, number> &entry : line.entries)
        entry.first += offset;
    }

#ifdef DEBUG
  // make sure that lines, lines_cache and local_lines
  // are still linked correctly
  for (size_type index = 0; index < lines_cache.size(); ++index)
    Assert(lines_cache[index] == numbers::invalid_size_type ||
             calculate_line_index(lines[lines_cache[index]].index) == index,
           ExcInternalError());
#endif
}



template <typename number>
AffineConstraints<number>
AffineConstraints<number>::get_view(const IndexSet &mask) const
{
  Assert(sorted == true, ExcMatrixNotClosed());

  AffineConstraints<number> view;

  // If index sets were associated with the current object, take views
  // of those as well:
  if (locally_owned_dofs != IndexSet())
    {
      Assert(locally_owned_dofs.size() == mask.size(),
             ExcMessage("This operation only makes sense if the index "
                        "space described by the mask is the same as "
                        "the index space associated with the "
                        "AffineConstraints object into which you take "
                        "a view."));
      Assert(local_lines.size() == mask.size(),
             ExcMessage("This operation only makes sense if the index "
                        "space described by the mask is the same as "
                        "the index space associated with the "
                        "AffineConstraints object into which you take "
                        "a view."));

      view.reinit(locally_owned_dofs.get_view(mask),
                  local_lines.get_view(mask));
    }

  for (const ConstraintLine &line : lines)
    if (mask.is_element(line.index))
      {
#ifdef DEBUG
        for (const std::pair<size_type, number> &entry : line.entries)
          {
            Assert(
              mask.is_element(entry.first),
              ExcMessage(
                "In creating a view of an AffineConstraints "
                "object, the constraint on degree of freedom " +
                std::to_string(line.index) + " (which corresponds to the " +
                std::to_string(mask.index_within_set(line.index)) +
                "th degree of freedom selected in the mask) "
                "is constrained against degree of freedom " +
                std::to_string(entry.first) +
                ", but this degree of freedom is not listed in the mask and "
                "consequently cannot be transcribed into the index space "
                "of the output object."));
          }
#endif

        std::vector<std::pair<size_type, number>> translated_entries =
          line.entries;
        for (auto &entry : translated_entries)
          entry.first = mask.index_within_set(entry.first);

        view.add_constraint(mask.index_within_set(line.index),
                            translated_entries,
                            line.inhomogeneity);
      }

  view.close();
  return view;
}



template <typename number>
void
AffineConstraints<number>::clear()
{
  {
    std::vector<ConstraintLine> tmp;
    lines.swap(tmp);
  }

  {
    std::vector<size_type> tmp;
    lines_cache.swap(tmp);
  }

  locally_owned_dofs             = {};
  local_lines                    = {};
  needed_elements_for_distribute = {};

  sorted = false;
}



template <typename number>
void
AffineConstraints<number>::reinit()
{
  reinit(IndexSet(), IndexSet());
}



template <typename number>
void
AffineConstraints<number>::reinit(const IndexSet &locally_stored_constraints)
{
  reinit(locally_stored_constraints, locally_stored_constraints);
}



template <typename number>
void
AffineConstraints<number>::reinit(const IndexSet &locally_owned_dofs,
                                  const IndexSet &locally_stored_constraints)
{
  // First clear previous content
  clear();

  // Then set the objects that describe the index sets of DoFs we care about:
  Assert(locally_owned_dofs.is_subset_of(locally_stored_constraints),
         ExcMessage("The set of locally stored constraints needs to be a "
                    "superset of the locally owned DoFs."));

  this->locally_owned_dofs = locally_owned_dofs;
  this->local_lines        = locally_stored_constraints;

  // make sure the IndexSet is compressed. Otherwise this can lead to crashes
  // that are hard to find (only happen in release mode).
  // see tests/mpi/affine_constraints_crash_01
  this->locally_owned_dofs.compress();
  this->local_lines.compress();
}



template <typename number>
bool
AffineConstraints<number>::is_identity_constrained(const size_type line_n) const
{
  if (is_constrained(line_n) == false)
    return false;

  const ConstraintLine &line = lines[lines_cache[calculate_line_index(line_n)]];
  Assert(line.index == line_n, ExcInternalError());

  // return if an entry for this line was found and if it has only one
  // entry equal to 1.0
  return (line.entries.size() == 1) && (line.entries[0].second == number(1.0));
}


template <typename number>
bool
AffineConstraints<number>::are_identity_constrained(
  const size_type line_n_1,
  const size_type line_n_2) const
{
  if (is_constrained(line_n_1) == true)
    {
      const ConstraintLine &line =
        lines[lines_cache[calculate_line_index(line_n_1)]];
      Assert(line.index == line_n_1, ExcInternalError());

      // return if an entry for this line was found and if it has only one
      // entry equal to 1.0 and that one is index2
      return ((line.entries.size() == 1) &&
              (line.entries[0].first == line_n_2) &&
              (line.entries[0].second == number(1.0)));
    }
  else if (is_constrained(line_n_2) == true)
    {
      const ConstraintLine &line =
        lines[lines_cache[calculate_line_index(line_n_2)]];
      Assert(line.index == line_n_2, ExcInternalError());

      // return if an entry for this line was found and if it has only one
      // entry equal to 1.0 and that one is line_n_1
      return ((line.entries.size() == 1) &&
              (line.entries[0].first == line_n_1) &&
              (line.entries[0].second == number(1.0)));
    }
  else
    return false;
}



template <typename number>
typename AffineConstraints<number>::size_type
AffineConstraints<number>::max_constraint_indirections() const
{
  size_type return_value = 0;
  for (const ConstraintLine &line : lines)
    // use static cast, since typeof(size)==std::size_t, which is !=
    // size_type on AIX
    return_value =
      std::max(return_value, static_cast<size_type>(line.entries.size()));

  return return_value;
}



template <typename number>
bool
AffineConstraints<number>::has_inhomogeneities() const
{
  for (const ConstraintLine &line : lines)
    if (line.inhomogeneity != number(0.))
      return true;

  return false;
}



template <typename number>
void
AffineConstraints<number>::print(std::ostream &out) const
{
  for (const ConstraintLine &line : lines)
    {
      // output the list of constraints as pairs of dofs and their weights
      if (line.entries.size() > 0)
        {
          for (const std::pair<size_type, number> &entry : line.entries)
            out << "    " << line.index << ' ' << entry.first << ":  "
                << entry.second << '\n';

          // print out inhomogeneity.
          if (line.inhomogeneity != number(0.))
            out << "    " << line.index << ": " << line.inhomogeneity << '\n';
        }
      else
        // but also output something if the constraint simply reads
        // x[13]=0, i.e. where the right hand side is not a linear
        // combination of other dofs
        {
          if (line.inhomogeneity != number(0.))
            out << "    " << line.index << " = " << line.inhomogeneity << '\n';
          else
            out << "    " << line.index << " = 0\n";
        }
    }

  AssertThrow(out.fail() == false, ExcIO());
}



template <typename number>
void
AffineConstraints<number>::write_dot(std::ostream &out) const
{
  out << "digraph constraints {" << std::endl;
  for (size_type i = 0; i != lines.size(); ++i)
    {
      // same concept as in the previous function
      if (lines[i].entries.size() > 0)
        for (size_type j = 0; j < lines[i].entries.size(); ++j)
          out << "  " << lines[i].index << "->" << lines[i].entries[j].first
              << "; // weight: " << lines[i].entries[j].second << '\n';
      else
        out << "  " << lines[i].index << '\n';
    }
  out << '}' << std::endl;
}



template <typename number>
std::size_t
AffineConstraints<number>::memory_consumption() const
{
  return (MemoryConsumption::memory_consumption(lines) +
          MemoryConsumption::memory_consumption(lines_cache) +
          MemoryConsumption::memory_consumption(sorted) +
          MemoryConsumption::memory_consumption(local_lines));
}



template <typename number>
void
AffineConstraints<number>::resolve_indices(
  std::vector<types::global_dof_index> &indices) const
{
  const unsigned int indices_size = indices.size();
  const std::vector<std::pair<types::global_dof_index, number>> *line_ptr;
  for (unsigned int i = 0; i < indices_size; ++i)
    {
      line_ptr = get_constraint_entries(indices[i]);
      // if the index is constraint, the constraints indices are added to the
      // indices vector
      if (line_ptr != nullptr)
        {
          const unsigned int line_size = line_ptr->size();
          for (unsigned int j = 0; j < line_size; ++j)
            indices.push_back((*line_ptr)[j].first);
        }
    }

  // keep only the unique elements
  std::sort(indices.begin(), indices.end());
  indices.erase(std::unique(indices.begin(), indices.end()), indices.end());
}



template <typename number>
void
AffineConstraints<number>::condense(SparsityPattern &sparsity) const
{
  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the AffineConstraints object
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
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
               ((entry != sparsity.end(row)) && entry->is_valid_entry());
               ++entry)
            {
              const size_type column = entry->column();

              if (distribute[column] != numbers::invalid_size_type)
                {
                  // distribute entry at regular row @p{row} and irregular
                  // column sparsity.colnums[j]
                  for (size_type q = 0;
                       q != lines[distribute[column]].entries.size();
                       ++q)
                    sparsity.add(row,
                                 lines[distribute[column]].entries[q].first);
                }
            }
        }
      else
        // row must be distributed. note that here the present row is not
        // touched (unlike above)
        {
          for (SparsityPattern::iterator entry = sparsity.begin(row);
               (entry != sparsity.end(row)) && entry->is_valid_entry();
               ++entry)
            {
              const size_type column = entry->column();
              if (distribute[column] == numbers::invalid_size_type)
                // distribute entry at irregular row @p{row} and regular
                // column sparsity.colnums[j]
                for (size_type q = 0;
                     q != lines[distribute[row]].entries.size();
                     ++q)
                  sparsity.add(lines[distribute[row]].entries[q].first, column);
              else
                // distribute entry at irregular row @p{row} and irregular
                // column sparsity.get_column_numbers()[j]
                for (size_type p = 0;
                     p != lines[distribute[row]].entries.size();
                     ++p)
                  for (size_type q = 0;
                       q != lines[distribute[column]].entries.size();
                       ++q)
                    sparsity.add(lines[distribute[row]].entries[p].first,
                                 lines[distribute[column]].entries[q].first);
            }
        }
    }

  sparsity.compress();
}



template <typename number>
void
AffineConstraints<number>::condense(BlockSparsityPattern &sparsity) const
{
  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.is_compressed() == false, ExcMatrixIsClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());
  Assert(sparsity.n_block_rows() == sparsity.n_block_cols(), ExcNotQuadratic());
  Assert(sparsity.get_column_indices() == sparsity.get_row_indices(),
         ExcNotQuadratic());

  const BlockIndices &index_mapping = sparsity.get_column_indices();

  const size_type n_blocks = sparsity.n_block_rows();

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the AffineConstraints object
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
    {
      // get index of this row within the blocks
      const std::pair<size_type, size_type> block_index =
        index_mapping.global_to_local(row);
      const size_type block_row = block_index.first;

      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over all columns and see whether this column
        // must be distributed
        {
          // to loop over all entries in this row, we have to loop over all
          // blocks in this blockrow and the corresponding row therein
          for (size_type block_col = 0; block_col < n_blocks; ++block_col)
            {
              const SparsityPattern &block_sparsity =
                sparsity.block(block_row, block_col);

              for (SparsityPattern::const_iterator entry =
                     block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                   entry->is_valid_entry();
                   ++entry)
                {
                  const size_type global_col =
                    index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at regular row @p{row} and
                    // irregular column global_col
                    {
                      for (size_type q = 0;
                           q != lines[distribute[global_col]].entries.size();
                           ++q)
                        sparsity.add(
                          row, lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
      else
        {
          // row must be distributed. split the whole row into the chunks
          // defined by the blocks
          for (size_type block_col = 0; block_col < n_blocks; ++block_col)
            {
              const SparsityPattern &block_sparsity =
                sparsity.block(block_row, block_col);

              for (SparsityPattern::const_iterator entry =
                     block_sparsity.begin(block_index.second);
                   (entry != block_sparsity.end(block_index.second)) &&
                   entry->is_valid_entry();
                   ++entry)
                {
                  const size_type global_col =
                    index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] == numbers::invalid_size_type)
                    // distribute entry at irregular row @p{row} and
                    // regular column global_col.
                    {
                      for (size_type q = 0;
                           q != lines[distribute[row]].entries.size();
                           ++q)
                        sparsity.add(lines[distribute[row]].entries[q].first,
                                     global_col);
                    }
                  else
                    // distribute entry at irregular row @p{row} and
                    // irregular column @p{global_col}
                    {
                      for (size_type p = 0;
                           p != lines[distribute[row]].entries.size();
                           ++p)
                        for (size_type q = 0;
                             q != lines[distribute[global_col]].entries.size();
                             ++q)
                          sparsity.add(
                            lines[distribute[row]].entries[p].first,
                            lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
    }

  sparsity.compress();
}



template <typename number>
void
AffineConstraints<number>::condense(DynamicSparsityPattern &sparsity) const
{
  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraints object
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over cols. note that as we proceed to
        // distribute cols, the loop may get longer
        for (size_type j = 0; j < sparsity.row_length(row); ++j)
          {
            const size_type column = sparsity.column_number(row, j);

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
                for (size_type q = 0;
                     q != lines[distribute[column]].entries.size();
                     ++q)
                  {
                    const size_type new_col =
                      lines[distribute[column]].entries[q].first;

                    sparsity.add(row, new_col);

                    const size_type new_rowlength = sparsity.row_length(row);
                    if ((new_col < column) && (old_rowlength != new_rowlength))
                      ++j;
                    old_rowlength = new_rowlength;
                  }
              }
          }
      else
        // row must be distributed
        for (size_type j = 0; j < sparsity.row_length(row); ++j)
          {
            const size_type column = sparsity.column_number(row, j);

            if (distribute[column] == numbers::invalid_size_type)
              // distribute entry at irregular row @p{row} and regular
              // column sparsity.colnums[j]
              for (size_type q = 0; q != lines[distribute[row]].entries.size();
                   ++q)
                sparsity.add(lines[distribute[row]].entries[q].first, column);
            else
              // distribute entry at irregular row @p{row} and irregular
              // column sparsity.get_column_numbers()[j]
              for (size_type p = 0; p != lines[distribute[row]].entries.size();
                   ++p)
                for (size_type q = 0;
                     q != lines[distribute[sparsity.column_number(row, j)]]
                            .entries.size();
                     ++q)
                  sparsity.add(lines[distribute[row]].entries[p].first,
                               lines[distribute[sparsity.column_number(row, j)]]
                                 .entries[q]
                                 .first);
          }
    }
}



template <typename number>
void
AffineConstraints<number>::condense(BlockDynamicSparsityPattern &sparsity) const
{
  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());
  Assert(sparsity.n_block_rows() == sparsity.n_block_cols(), ExcNotQuadratic());
  Assert(sparsity.get_column_indices() == sparsity.get_row_indices(),
         ExcNotQuadratic());

  const BlockIndices &index_mapping = sparsity.get_column_indices();

  const size_type n_blocks = sparsity.n_block_rows();

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_unsigned_int, no distribution is necessary.
  // otherwise, the number states which line in the constraints object
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = static_cast<signed int>(c);

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
    {
      // get index of this row within the blocks
      const std::pair<size_type, size_type> block_index =
        index_mapping.global_to_local(row);
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
          for (size_type block_col = 0; block_col < n_blocks; ++block_col)
            {
              const DynamicSparsityPattern &block_sparsity =
                sparsity.block(block_row, block_col);

              for (size_type j = 0; j < block_sparsity.row_length(local_row);
                   ++j)
                {
                  const size_type global_col = index_mapping.local_to_global(
                    block_col, block_sparsity.column_number(local_row, j));

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at regular row @p{row} and
                    // irregular column global_col
                    {
                      for (size_type q = 0;
                           q != lines[distribute[global_col]].entries.size();
                           ++q)
                        sparsity.add(
                          row, lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
      else
        {
          // row must be distributed. split the whole row into the chunks
          // defined by the blocks
          for (size_type block_col = 0; block_col < n_blocks; ++block_col)
            {
              const DynamicSparsityPattern &block_sparsity =
                sparsity.block(block_row, block_col);

              for (size_type j = 0; j < block_sparsity.row_length(local_row);
                   ++j)
                {
                  const size_type global_col = index_mapping.local_to_global(
                    block_col, block_sparsity.column_number(local_row, j));

                  if (distribute[global_col] == numbers::invalid_size_type)
                    // distribute entry at irregular row @p{row} and
                    // regular column global_col.
                    {
                      for (size_type q = 0;
                           q != lines[distribute[row]].entries.size();
                           ++q)
                        sparsity.add(lines[distribute[row]].entries[q].first,
                                     global_col);
                    }
                  else
                    // distribute entry at irregular row @p{row} and
                    // irregular column @p{global_col}
                    {
                      for (size_type p = 0;
                           p != lines[distribute[row]].entries.size();
                           ++p)
                        for (size_type q = 0;
                             q != lines[distribute[global_col]].entries.size();
                             ++q)
                          sparsity.add(
                            lines[distribute[row]].entries[p].first,
                            lines[distribute[global_col]].entries[q].first);
                    }
                }
            }
        }
    }
}



template <typename number>
void
AffineConstraints<number>::condense(SparseMatrix<number> &uncondensed) const
{
  Vector<number> dummy(0);
  condense(uncondensed, dummy);
}



template <typename number>
void
AffineConstraints<number>::condense(
  BlockSparseMatrix<number> &uncondensed) const
{
  BlockVector<number> dummy(0);
  condense(uncondensed, dummy);
}



template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::condense(const VectorType &vec_ghosted,
                                    VectorType       &vec) const
{
  Assert(sorted == true, ExcMatrixNotClosed());

  // if this is called with different arguments, we need to copy the data
  // over:
  if (&vec != &vec_ghosted)
    vec = vec_ghosted;

  // distribute all entries, and set them to zero. do so in two loops
  // because in the first one we need to add to elements and in the second
  // one we need to set elements to zero. for parallel vectors, this can
  // only work if we can put a compress() in between, but we don't want to
  // call compress() twice per entry
  for (const ConstraintLine &line : lines)
    {
      // in case the constraint is inhomogeneous, this function is not
      // appropriate. Throw an exception.
      Assert(line.inhomogeneity == number(0.),
             ExcMessage("Inhomogeneous constraint cannot be condensed "
                        "without any matrix specified."));

      const typename VectorType::value_type old_value = vec_ghosted(line.index);
      for (const std::pair<size_type, number> &entry : line.entries)
        if (vec.in_local_range(entry.first) == true)
          vec(entry.first) +=
            (static_cast<typename VectorType::value_type>(old_value) *
             entry.second);
    }

  vec.compress(VectorOperation::add);

  for (const ConstraintLine &line : lines)
    if (vec.in_local_range(line.index) == true)
      vec(line.index) = 0.;

  vec.compress(VectorOperation::insert);
}



template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::condense(VectorType &vec) const
{
  condense(vec, vec);
}



template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::condense(SparseMatrix<number> &uncondensed,
                                    VectorType           &vec) const
{
  // check whether we work on real vectors or we just used a dummy when
  // calling the other function above.
  const bool use_vectors = vec.size() == 0 ? false : true;

  const SparsityPattern &sparsity = uncondensed.get_sparsity_pattern();

  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());
  if (use_vectors == true)
    AssertDimension(vec.size(), sparsity.n_rows());

  number average_diagonal = 0.;
  for (size_type i = 0; i < uncondensed.m(); ++i)
    average_diagonal += std::abs(uncondensed.diag_element(i));
  average_diagonal /= uncondensed.m();

  // store for each index whether it must be distributed or not. If entry
  // is invalid_size_type, no distribution is necessary.  otherwise, the
  // number states which line in the constraints object handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
    {
      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over cols
        {
          for (typename SparseMatrix<number>::iterator entry =
                 uncondensed.begin(row);
               entry != uncondensed.end(row);
               ++entry)
            {
              const size_type column = entry->column();

              // end of row reached? this should not happen, since we only
              // operate on compressed matrices!
              Assert(column != SparsityPattern::invalid_entry,
                     ExcMatrixNotClosed());

              if (distribute[column] != numbers::invalid_size_type)
                // distribute entry at regular row @p row and irregular
                // column sparsity.get_column_numbers()[j]; set old entry
                // to zero
                {
                  for (size_type q = 0;
                       q != lines[distribute[column]].entries.size();
                       ++q)
                    {
                      // need a temporary variable to avoid errors like no
                      // known conversion from 'complex<typename
                      // ProductType<float, double>::type>' to 'const
                      // complex<float>' for 3rd argument
                      number v = static_cast<number>(entry->value());
                      v *= lines[distribute[column]].entries[q].second;
                      uncondensed.add(
                        row, lines[distribute[column]].entries[q].first, v);
                    }

                  // need to subtract this element from the vector. this
                  // corresponds to an explicit elimination in the
                  // respective row of the inhomogeneous constraint in the
                  // matrix with Gauss elimination
                  if (use_vectors == true)
                    vec(row) -= static_cast<number>(entry->value()) *
                                lines[distribute[column]].inhomogeneity;

                  // set old value to zero
                  entry->value() = 0.;
                }
            }
        }
      else
        // row must be distributed
        {
          for (typename SparseMatrix<number>::iterator entry =
                 uncondensed.begin(row);
               entry != uncondensed.end(row);
               ++entry)
            {
              const size_type column = entry->column();

              // end of row reached? this should not happen, since we only
              // operate on compressed matrices!
              Assert(column != SparsityPattern::invalid_entry,
                     ExcMatrixNotClosed());

              if (distribute[column] == numbers::invalid_size_type)
                // distribute entry at irregular row @p row and regular
                // column column. set old entry to zero
                {
                  for (size_type q = 0;
                       q != lines[distribute[row]].entries.size();
                       ++q)
                    {
                      // need a temporary variable to avoid errors like
                      // no known conversion from 'complex<typename
                      // ProductType<float, double>::type>' to 'const
                      // complex<float>' for 3rd argument
                      number v = static_cast<number>(entry->value());
                      v *= lines[distribute[row]].entries[q].second;
                      uncondensed.add(lines[distribute[row]].entries[q].first,
                                      column,
                                      v);
                    }

                  // set old entry to zero
                  entry->value() = 0.;
                }
              else
                // distribute entry at irregular row @p row and irregular
                // column @p column set old entry to one on main diagonal,
                // zero otherwise
                {
                  for (size_type p = 0;
                       p != lines[distribute[row]].entries.size();
                       ++p)
                    {
                      for (size_type q = 0;
                           q != lines[distribute[column]].entries.size();
                           ++q)
                        {
                          // need a temporary variable to avoid errors like
                          // no known conversion from 'complex<typename
                          // ProductType<float, double>::type>' to 'const
                          // complex<float>' for 3rd argument
                          number v = static_cast<number>(entry->value());
                          v *= lines[distribute[row]].entries[p].second *
                               lines[distribute[column]].entries[q].second;
                          uncondensed.add(
                            lines[distribute[row]].entries[p].first,
                            lines[distribute[column]].entries[q].first,
                            v);
                        }

                      if (use_vectors == true)
                        vec(lines[distribute[row]].entries[p].first) -=
                          static_cast<number>(entry->value()) *
                          lines[distribute[row]].entries[p].second *
                          lines[distribute[column]].inhomogeneity;
                    }

                  // set old entry to correct value
                  entry->value() = (row == column ? average_diagonal : 0.);
                }
            }

          // take care of vector
          if (use_vectors == true)
            {
              for (size_type q = 0; q != lines[distribute[row]].entries.size();
                   ++q)
                vec(lines[distribute[row]].entries[q].first) +=
                  (vec(row) * lines[distribute[row]].entries[q].second);

              vec(lines[distribute[row]].index) = 0.;
            }
        }
    }
}



template <typename number>
template <typename BlockVectorType>
void
AffineConstraints<number>::condense(BlockSparseMatrix<number> &uncondensed,
                                    BlockVectorType           &vec) const
{
  // check whether we work on real vectors or we just used a dummy when
  // calling the other function above.
  const bool use_vectors = vec.n_blocks() == 0 ? false : true;

  const size_type blocks = uncondensed.n_block_rows();

  const BlockSparsityPattern &sparsity = uncondensed.get_sparsity_pattern();

  Assert(sorted == true, ExcMatrixNotClosed());
  Assert(sparsity.is_compressed() == true, ExcMatrixNotClosed());
  Assert(sparsity.n_rows() == sparsity.n_cols(), ExcNotQuadratic());
  Assert(sparsity.n_block_rows() == sparsity.n_block_cols(), ExcNotQuadratic());
  Assert(sparsity.n_block_rows() == sparsity.n_block_cols(), ExcNotQuadratic());
  Assert(sparsity.get_column_indices() == sparsity.get_row_indices(),
         ExcNotQuadratic());

  if (use_vectors == true)
    {
      AssertDimension(vec.size(), sparsity.n_rows());
      AssertDimension(vec.n_blocks(), sparsity.n_block_rows());
    }

  number average_diagonal = 0.;
  for (size_type b = 0; b < uncondensed.n_block_rows(); ++b)
    for (size_type i = 0; i < uncondensed.block(b, b).m(); ++i)
      average_diagonal += std::abs(uncondensed.block(b, b).diag_element(i));
  average_diagonal /= uncondensed.m();

  const BlockIndices &index_mapping = sparsity.get_column_indices();

  // store for each index whether it must be distributed or not. If entry
  // is numbers::invalid_size_type, no distribution is necessary.
  // otherwise, the number states which line in the constraints object
  // handles this index
  std::vector<size_type> distribute(sparsity.n_rows(),
                                    numbers::invalid_size_type);

  for (size_type c = 0; c < lines.size(); ++c)
    distribute[lines[c].index] = c;

  const size_type n_rows = sparsity.n_rows();
  for (size_type row = 0; row < n_rows; ++row)
    {
      // get index of this row within the blocks
      const std::pair<size_type, size_type> block_index =
        index_mapping.global_to_local(row);
      const size_type block_row = block_index.first;

      if (distribute[row] == numbers::invalid_size_type)
        // regular line. loop over all columns and see whether this column
        // must be distributed
        {
          // to loop over all entries in this row, we have to loop over all
          // blocks in this blockrow and the corresponding row therein
          for (size_type block_col = 0; block_col < blocks; ++block_col)
            {
              for (typename SparseMatrix<number>::iterator entry =
                     uncondensed.block(block_row, block_col)
                       .begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col)
                              .end(block_index.second);
                   ++entry)
                {
                  const size_type global_col =
                    index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] != numbers::invalid_size_type)
                    // distribute entry at regular row @p row and irregular
                    // column global_col; set old entry to zero
                    {
                      const number old_value = entry->value();

                      for (size_type q = 0;
                           q != lines[distribute[global_col]].entries.size();
                           ++q)
                        uncondensed.add(
                          row,
                          lines[distribute[global_col]].entries[q].first,
                          old_value *
                            lines[distribute[global_col]].entries[q].second);

                      // need to subtract this element from the vector.
                      // this corresponds to an explicit elimination in the
                      // respective row of the inhomogeneous constraint in
                      // the matrix with Gauss elimination
                      if (use_vectors == true)
                        vec(row) -= static_cast<number>(entry->value()) *
                                    lines[distribute[global_col]].inhomogeneity;

                      entry->value() = 0.;
                    }
                }
            }
        }
      else
        {
          // row must be distributed. split the whole row into the chunks
          // defined by the blocks
          for (size_type block_col = 0; block_col < blocks; ++block_col)
            {
              for (typename SparseMatrix<number>::iterator entry =
                     uncondensed.block(block_row, block_col)
                       .begin(block_index.second);
                   entry != uncondensed.block(block_row, block_col)
                              .end(block_index.second);
                   ++entry)
                {
                  const size_type global_col =
                    index_mapping.local_to_global(block_col, entry->column());

                  if (distribute[global_col] == numbers::invalid_size_type)
                    // distribute entry at irregular row @p row and regular
                    // column global_col. set old entry to zero
                    {
                      const number old_value = entry->value();

                      for (size_type q = 0;
                           q != lines[distribute[row]].entries.size();
                           ++q)
                        uncondensed.add(
                          lines[distribute[row]].entries[q].first,
                          global_col,
                          old_value * lines[distribute[row]].entries[q].second);

                      entry->value() = 0.;
                    }
                  else
                    // distribute entry at irregular row @p row and
                    // irregular column @p global_col set old entry to one
                    // if on main diagonal, zero otherwise
                    {
                      const number old_value = entry->value();

                      for (size_type p = 0;
                           p != lines[distribute[row]].entries.size();
                           ++p)
                        {
                          for (size_type q = 0;
                               q !=
                               lines[distribute[global_col]].entries.size();
                               ++q)
                            uncondensed.add(
                              lines[distribute[row]].entries[p].first,
                              lines[distribute[global_col]].entries[q].first,
                              old_value *
                                lines[distribute[row]].entries[p].second *
                                lines[distribute[global_col]]
                                  .entries[q]
                                  .second);

                          if (use_vectors == true)
                            vec(lines[distribute[row]].entries[p].first) -=
                              old_value *
                              lines[distribute[row]].entries[p].second *
                              lines[distribute[global_col]].inhomogeneity;
                        }

                      entry->value() =
                        (row == global_col ? average_diagonal : 0.);
                    }
                }
            }

          // take care of vector
          if (use_vectors == true)
            {
              for (size_type q = 0; q != lines[distribute[row]].entries.size();
                   ++q)
                vec(lines[distribute[row]].entries[q].first) +=
                  (vec(row) * lines[distribute[row]].entries[q].second);

              vec(lines[distribute[row]].index) = 0.;
            }
        }
    }
}



// TODO: I'm sure the following could be made more elegant by using a bit of
// introspection using static member variables of the various vector
// classes to dispatch between the different functions, rather than using
// knowledge of the individual types

// number of functions to select the right implementation for set_zero().
namespace internal
{
  namespace AffineConstraintsImplementation
  {
    using size_type = types::global_dof_index;

    template <typename VectorType>
    void
    set_zero_parallel(const std::vector<size_type> &cm,
                      VectorType                   &vec,
                      size_type                     shift = 0)
    {
      Assert(!vec.has_ghost_elements(), ExcInternalError());
      IndexSet locally_owned = vec.locally_owned_elements();
      for (const size_type index : cm)
        {
          // If shift>0 then we are working on a part of a BlockVector
          // so vec(i) is actually the global entry i+shift.
          // We first make sure the line falls into the range of vec,
          // then check if is part of the local part of the vector, before
          // finally setting value to 0.
          if (index < shift)
            continue;
          const size_type idx = index - shift;
          if (idx < vec.size() && locally_owned.is_element(idx))
            internal::ElementAccess<VectorType>::set(0., idx, vec);
        }
    }

    template <typename number>
    void
    set_zero_parallel(const std::vector<size_type>               &cm,
                      LinearAlgebra::distributed::Vector<number> &vec,
                      size_type                                   shift = 0)
    {
      for (const size_type index : cm)
        {
          // If shift>0 then we are working on a part of a BlockVector
          // so vec(i) is actually the global entry i+shift.
          // We first make sure the line falls into the range of vec,
          // then check if is part of the local part of the vector, before
          // finally setting the value to 0.
          if (index < shift)
            continue;
          const size_type idx = index - shift;
          if (vec.in_local_range(idx))
            vec(idx) = 0.;
        }
      vec.zero_out_ghost_values();
    }

    template <typename number>
    void
    set_zero_parallel(
      const std::vector<size_type>                                     &cm,
      LinearAlgebra::distributed::Vector<number, MemorySpace::Default> &vec,
      size_type shift = 0)
    {
      Assert(shift == 0, ExcNotImplemented());
      (void)shift;
      std::vector<size_type> constrained_local_dofs_host;
      constrained_local_dofs_host.reserve(cm.size());

      for (const auto global_index : cm)
        if (vec.in_local_range(global_index))
          constrained_local_dofs_host.push_back(
            vec.get_partitioner()->global_to_local(global_index));

      const int n_constraints = constrained_local_dofs_host.size();
      Kokkos::View<size_type *, MemorySpace::Default::kokkos_space>
        constrained_local_dofs_device(
          Kokkos::view_alloc(Kokkos::WithoutInitializing,
                             "constrained_local_dofs_device"),
          n_constraints);
      Kokkos::deep_copy(constrained_local_dofs_device,
                        Kokkos::View<size_type *, Kokkos::HostSpace>(
                          constrained_local_dofs_host.data(),
                          constrained_local_dofs_host.size()));

      using ExecutionSpace =
        MemorySpace::Default::kokkos_space::execution_space;
      ExecutionSpace exec;
      auto          *local_values = vec.get_values();
      Kokkos::parallel_for(
        "dealii::set_zero_parallel",
        Kokkos::RangePolicy<ExecutionSpace>(exec, 0, n_constraints),
        KOKKOS_LAMBDA(int i) {
          local_values[constrained_local_dofs_device[i]] = 0;
        });

      vec.zero_out_ghost_values();
    }

    template <typename VectorType>
    void
    set_zero_in_parallel(const std::vector<size_type> &cm,
                         VectorType                   &vec,
                         std::bool_constant<false>)
    {
      set_zero_parallel(cm, vec, 0);
    }

    // in parallel for BlockVectors
    template <typename VectorType>
    void
    set_zero_in_parallel(const std::vector<size_type> &cm,
                         VectorType                   &vec,
                         std::bool_constant<true>)
    {
      size_type start_shift = 0;
      for (size_type j = 0; j < vec.n_blocks(); ++j)
        {
          set_zero_parallel(cm, vec.block(j), start_shift);
          start_shift += vec.block(j).size();
        }
    }

    template <typename VectorType>
    void
    set_zero_serial(const std::vector<size_type> &cm, VectorType &vec)
    {
      for (const size_type index : cm)
        vec(index) = 0.;
    }

    template <typename VectorType>
    void
    set_zero_all(const std::vector<size_type> &cm, VectorType &vec)
    {
      set_zero_in_parallel<VectorType>(
        cm, vec, std::bool_constant<IsBlockVector<VectorType>::value>());
      vec.compress(VectorOperation::insert);
    }

    template <class T>
    void
    set_zero_all(const std::vector<size_type> &cm, dealii::Vector<T> &vec)
    {
      set_zero_serial(cm, vec);
    }

    template <class T>
    void
    set_zero_all(const std::vector<size_type> &cm, dealii::BlockVector<T> &vec)
    {
      set_zero_serial(cm, vec);
    }
  } // namespace AffineConstraintsImplementation
} // namespace internal



template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::distribute_local_to_global(
  const Vector<number>         &local_vector,
  const std::vector<size_type> &local_dof_indices,
  VectorType                   &global_vector,
  const FullMatrix<number>     &local_matrix) const
{
  distribute_local_to_global(local_vector,
                             local_dof_indices,
                             local_dof_indices,
                             global_vector,
                             local_matrix,
                             true);
}



template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::distribute_local_to_global(
  const Vector<number>         &local_vector,
  const std::vector<size_type> &local_dof_indices_row,
  const std::vector<size_type> &local_dof_indices_col,
  VectorType                   &global_vector,
  const FullMatrix<number>     &local_matrix,
  bool                          diagonal) const
{
  Assert(sorted == true, ExcMatrixNotClosed());
  AssertDimension(local_vector.size(), local_dof_indices_row.size());
  AssertDimension(local_matrix.m(), local_dof_indices_row.size());
  AssertDimension(local_matrix.n(), local_dof_indices_col.size());
  Assert(global_vector.has_ghost_elements() == false, ExcGhostsPresent());

  // diagonal checks if we have only one index set (if both are equal
  // diagonal should be set to true).
  // If true we do both, assembly of the right hand side (next lines)
  // and (see further below) modifications of the right hand side
  // according to the inhomogeneous constraints.
  // Otherwise we only modify the right hand side according to
  // local_matrix and the inhomogeneous constraints, and omit the vector add.

  const size_type m_local_dofs = local_dof_indices_row.size();
  const size_type n_local_dofs = local_dof_indices_col.size();
  if (lines.empty())
    {
      if (diagonal)
        global_vector.add(local_dof_indices_row, local_vector);
    }
  else
    for (size_type i = 0; i < n_local_dofs; ++i)
      {
        // check whether the current index is
        // constrained. if not, just write the entry
        // into the vector. otherwise, need to resolve
        // the constraint
        if (is_constrained(local_dof_indices_col[i]) == false)
          {
            if (diagonal)
              global_vector(local_dof_indices_row[i]) += local_vector(i);
            continue;
          }

        // find the constraint line to the given
        // global dof index
        const size_type line_index =
          calculate_line_index(local_dof_indices_col[i]);
        AssertIndexRange(line_index, lines_cache.size());
        AssertIndexRange(lines_cache[line_index], lines.size());
        const ConstraintLine &position = lines[lines_cache[line_index]];

        // Gauss elimination of the matrix columns with the inhomogeneity.
        // Go through them one by one and again check whether they are
        // constrained. If so, distribute the constraint
        const auto val = position.inhomogeneity;
        if (val != number(0.))
          for (size_type j = 0; j < m_local_dofs; ++j)
            {
              if (is_constrained(local_dof_indices_row[j]) == false)
                {
                  global_vector(local_dof_indices_row[j]) -=
                    val * local_matrix(j, i);
                  continue;
                }

              const number matrix_entry = local_matrix(j, i);

              if (matrix_entry == number())
                continue;

              const ConstraintLine &position_j =
                lines[lines_cache[calculate_line_index(
                  local_dof_indices_row[j])]];

              for (size_type q = 0; q < position_j.entries.size(); ++q)
                {
                  Assert(!(!local_lines.size() ||
                           local_lines.is_element(
                             position_j.entries[q].first)) ||
                           is_constrained(position_j.entries[q].first) == false,
                         ExcMessage("Tried to distribute to a fixed dof."));
                  global_vector(position_j.entries[q].first) -=
                    val * position_j.entries[q].second * matrix_entry;
                }
            }

        // now distribute the constraint,
        // but make sure we don't touch
        // the entries of fixed dofs
        if (diagonal)
          {
            for (size_type j = 0; j < position.entries.size(); ++j)
              {
                Assert(!(!local_lines.size() ||
                         local_lines.is_element(position.entries[j].first)) ||
                         is_constrained(position.entries[j].first) == false,
                       ExcMessage("Tried to distribute to a fixed dof."));
                global_vector(position.entries[j].first) +=
                  local_vector(i) * position.entries[j].second;
              }
          }
      }
}

namespace internal
{
  // create an output vector that consists of the input vector's locally owned
  // elements plus some ghost elements that need to be imported from elsewhere
  //
  // this is an operation that is different for all vector types and so we
  // need a few overloads
#ifdef DEAL_II_WITH_TRILINOS
  inline void
  import_vector_with_ghost_elements(
    const TrilinosWrappers::MPI::Vector &vec,
    const IndexSet & /*locally_owned_elements*/,
    const IndexSet                &needed_elements,
    TrilinosWrappers::MPI::Vector &output,
    const std::bool_constant<false> /*is_block_vector*/)
  {
    Assert(!vec.has_ghost_elements(), ExcGhostsPresent());
#  ifdef DEAL_II_WITH_MPI
    const Epetra_MpiComm *mpi_comm =
      dynamic_cast<const Epetra_MpiComm *>(&vec.trilinos_vector().Comm());

    Assert(mpi_comm != nullptr, ExcInternalError());
    output.reinit(needed_elements, mpi_comm->GetMpiComm());
#  else
    output.reinit(needed_elements, MPI_COMM_SELF);
#  endif
    output = vec;
  }
#endif



#ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number, typename MemorySpace>
  inline void
  import_vector_with_ghost_elements(
    const LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &vec,
    const IndexSet &locally_owned_elements,
    const IndexSet &needed_elements,
    LinearAlgebra::TpetraWrappers::Vector<Number, MemorySpace> &output,
    const std::bool_constant<false> /*is_block_vector*/)
  {
    Assert(!vec.has_ghost_elements(), ExcGhostsPresent());
    IndexSet parallel_partitioner = locally_owned_elements;
    parallel_partitioner.add_indices(needed_elements);

    const MPI_Comm mpi_comm = Utilities::Trilinos::teuchos_comm_to_mpi_comm(
      vec.trilinos_vector().getMap()->getComm());

    output.reinit(locally_owned_elements, needed_elements, mpi_comm);

    output = vec;
  }
#endif // DEAL_II_TRILINOS_WITH_TPETRA



#ifdef DEAL_II_WITH_PETSC
  inline void
  import_vector_with_ghost_elements(
    const PETScWrappers::MPI::Vector &vec,
    const IndexSet                   &locally_owned_elements,
    const IndexSet                   &needed_elements,
    PETScWrappers::MPI::Vector       &output,
    const std::bool_constant<false> /*is_block_vector*/)
  {
    output.reinit(locally_owned_elements,
                  needed_elements,
                  vec.get_mpi_communicator());
    output = vec;
  }
#endif



  template <typename number>
  void
  import_vector_with_ghost_elements(
    const LinearAlgebra::distributed::Vector<number> &vec,
    const IndexSet                                   &locally_owned_elements,
    const IndexSet                                   &needed_elements,
    LinearAlgebra::distributed::Vector<number>       &output,
    const std::bool_constant<false> /*is_block_vector*/)
  {
    // TODO: the in vector might already have all elements. need to find a
    // way to efficiently avoid the copy then
    const_cast<LinearAlgebra::distributed::Vector<number> &>(vec)
      .zero_out_ghost_values();
    output.reinit(locally_owned_elements,
                  needed_elements,
                  vec.get_mpi_communicator());
    output = vec;
    output.update_ghost_values();
  }



  // all other vector non-block vector types are sequential and we should
  // not have this function called at all -- so throw an exception
  template <typename Vector>
  void
  import_vector_with_ghost_elements(
    const Vector & /*vec*/,
    const IndexSet & /*locally_owned_elements*/,
    const IndexSet & /*needed_elements*/,
    Vector & /*output*/,
    const std::bool_constant<false> /*is_block_vector*/)
  {
    Assert(false, ExcMessage("We shouldn't even get here!"));
  }

  // for block vectors, simply dispatch to the individual blocks
  template <typename VectorType>
  void
  import_vector_with_ghost_elements(
    const VectorType &vec,
    const IndexSet   &locally_owned_elements,
    const IndexSet   &needed_elements,
    VectorType       &output,
    const std::bool_constant<true> /*is_block_vector*/)
  {
    output.reinit(vec.n_blocks());

    types::global_dof_index block_start = 0;
    for (unsigned int b = 0; b < vec.n_blocks(); ++b)
      {
        import_vector_with_ghost_elements(
          vec.block(b),
          locally_owned_elements.get_view(block_start,
                                          block_start + vec.block(b).size()),
          needed_elements.get_view(block_start,
                                   block_start + vec.block(b).size()),
          output.block(b),
          std::bool_constant<false>());
        block_start += vec.block(b).size();
      }

    output.collect_sizes();
  }
} // namespace internal


namespace internal
{
  namespace AffineConstraints
  {
    template <typename VectorType,
              std::enable_if_t<internal::is_block_vector<VectorType> == false>
                * = nullptr>
    MPI_Comm
    get_mpi_communicator(const VectorType &vec)
    {
      return vec.get_mpi_communicator();
    }



    template <typename VectorType,
              std::enable_if_t<internal::is_block_vector<VectorType> == true>
                * = nullptr>
    MPI_Comm
    get_mpi_communicator(const VectorType &vec)
    {
      Assert(vec.n_blocks() > 0, ExcInternalError());
      return vec.block(0).get_mpi_communicator();
    }
  } // namespace AffineConstraints
} // namespace internal


template <typename number>
template <typename VectorType>
void
AffineConstraints<number>::distribute(VectorType &vec) const
{
  Assert(sorted == true, ExcMatrixNotClosed());

  // if the vector type supports parallel storage and if the vector actually
  // does store only part of the vector, distributing is slightly more
  // complicated. we might be able to skip the complicated part if one
  // processor owns everything and pretend that this is a sequential vector,
  // but it is difficult for the other processors to know whether they should
  // not do anything or if other processors will create a temporary vector,
  // exchange data (requiring communication, maybe even with the processors
  // that do not own anything because of that particular parallel model), and
  // call compress() finally. the first case here is for the complicated case,
  // the last else is for the simple case (sequential vector)
  const IndexSet vec_owned_elements = vec.locally_owned_elements();

  if constexpr (dealii::is_serial_vector<VectorType>::value == false)
    {
      // First check whether there are any constraints at all. If
      // there aren't, then there is nothing we need to do and we can
      // save the considerable effort below to communicate information
      // -- this is for sure the case on the first (not yet adaptively
      // refined) level of computations, for example:
      if (Utilities::MPI::sum(lines.size(),
                              internal::AffineConstraints::get_mpi_communicator(
                                vec)) > 0)
        {
          // We are looking at parallel vectors here. Of course, they could be
          // used in sequential contexts, but at least if we know that we are
          // in a parallel context with more than one process, we should check
          // that the user of this class provided index sets to the constructor
          // or the reinit() functions:
          if (Utilities::MPI::n_mpi_processes(
                internal::AffineConstraints::get_mpi_communicator(vec)) > 1)
            {
              Assert(local_lines != IndexSet(),
                     ExcMessage(
                       "You are using the AffineConstraints class in a "
                       "program that is running in parallel. In this case, "
                       "it is inefficient to store all constraints: This class "
                       "should really only store constraints for those degrees "
                       "of freedom that are locally relevant. It is considered "
                       "a bug not to do that. Please call either the "
                       "constructor of this class, or one of its reinit() "
                       "functions that provide this class with index sets "
                       "of the locally owned and the locally relevant "
                       "degrees of freedom."));

              Assert(vec_owned_elements.size() == locally_owned_dofs.size(),
                     ExcMessage("You have previously initialized this "
                                "AffineConstraints object with an index set "
                                "that stated that vectors have size " +
                                std::to_string(locally_owned_dofs.size()) +
                                " entries, but you are now calling "
                                "AffineConstraints::distribute() with a vector "
                                "of size " +
                                std::to_string(vec_owned_elements.size()) +
                                "."));
            }


            // Check that the set of indices we will import is a superset of
            // the locally-owned ones. This *should* be the case if, as one
            // would expect, the AffineConstraint object was initialized
            // with a locally-relevant index set that is indeed a superset
            // of the locally-owned indices. But you never know what people
            // pass as arguments...
#ifdef DEBUG
          if (needed_elements_for_distribute != IndexSet())
            {
              Assert(vec_owned_elements.size() ==
                       needed_elements_for_distribute.size(),
                     ExcMessage("You have previously initialized this "
                                "AffineConstraints object with an index set "
                                "that stated that vectors have size " +
                                std::to_string(locally_owned_dofs.size()) +
                                " entries, but you are now calling "
                                "AffineConstraints::distribute() with a vector "
                                "of size " +
                                std::to_string(vec_owned_elements.size()) +
                                "."));

              for (const auto i : vec_owned_elements)
                Assert(needed_elements_for_distribute.is_element(i),
                       ExcInternalError());
            }
#endif

          VectorType ghosted_vector;

          // It is possible that the user is using a parallel vector type,
          // but is running a non-parallel program (for example, step-31
          // does this). In this case, it is (perhaps?) not a bug to not
          // set an IndexSet for the local_lines and the
          // locally_owned_lines -- they are simply both empty sets in
          // that case. If that is so, we could just assign
          // 'ghosted_vector=vec;'. But this is dangerous. Not having set
          // index sets could also have been a bug, and we would get
          // downstream errors about accessing elements that are not
          // locally available. Rather, if no index sets were provided,
          // simply import the *entire* vector.
          if (local_lines != IndexSet())
            {
              Assert(vec_owned_elements.size() ==
                       needed_elements_for_distribute.size(),
                     ExcInternalError());

              Assert(needed_elements_for_distribute != IndexSet(),
                     ExcInternalError());
              internal::import_vector_with_ghost_elements(
                vec,
                vec_owned_elements,
                needed_elements_for_distribute,
                ghosted_vector,
                std::bool_constant<IsBlockVector<VectorType>::value>());
            }
          else
            {
              Assert(needed_elements_for_distribute == IndexSet(),
                     ExcInternalError());

              // TODO: We should really consider it a bug if this parallel
              // truly is distributed (and not just a parallel vector type
              // used for a sequential program). Assert that we are really
              // working in a sequential context.

              internal::import_vector_with_ghost_elements(
                vec,
                vec_owned_elements,
                complete_index_set(vec_owned_elements.size()),
                ghosted_vector,
                std::bool_constant<IsBlockVector<VectorType>::value>());
            }

          for (const ConstraintLine &line : lines)
            if (vec_owned_elements.is_element(line.index))
              {
                typename VectorType::value_type new_value = line.inhomogeneity;
                for (const std::pair<size_type, number> &entry : line.entries)
                  new_value +=
                    (static_cast<typename VectorType::value_type>(
                       internal::ElementAccess<VectorType>::get(ghosted_vector,
                                                                entry.first)) *
                     entry.second);
                AssertIsFinite(new_value);
                internal::ElementAccess<VectorType>::set(new_value,
                                                         line.index,
                                                         vec);
              }

          // now compress to communicate the entries that we added to
          // and that weren't to local processors to the owner
          //
          // this shouldn't be strictly necessary but it probably doesn't
          // hurt either
          vec.compress(VectorOperation::insert);
        }
    }
  else
    // purely sequential vector (either because the type doesn't
    // support anything else or because it's completely stored
    // locally)
    {
      for (const ConstraintLine &next_constraint : lines)
        {
          // fill entry in line
          // next_constraint.index by adding the
          // different contributions
          typename VectorType::value_type new_value =
            next_constraint.inhomogeneity;
          for (const std::pair<size_type, number> &entry :
               next_constraint.entries)
            new_value +=
              (static_cast<typename VectorType::value_type>(
                 internal::ElementAccess<VectorType>::get(vec, entry.first)) *
               entry.second);
          AssertIsFinite(new_value);
          internal::ElementAccess<VectorType>::set(new_value,
                                                   next_constraint.index,
                                                   vec);
        }
    }
}

// Some helper definitions for the local_to_global functions.
namespace internal
{
  namespace AffineConstraints
  {
    inline Distributing::Distributing(const size_type global_row,
                                      const size_type local_row)
      : global_row(global_row)
      , local_row(local_row)
      , constraint_position(numbers::invalid_size_type)
    {}



    inline Distributing::Distributing(const Distributing &in)
      : constraint_position(numbers::invalid_size_type)
    {
      *this = (in);
    }



    inline Distributing &
    Distributing::operator=(const Distributing &in)
    {
      global_row = in.global_row;
      local_row  = in.local_row;
      // the constraints pointer should not contain any data here.
      Assert(constraint_position == numbers::invalid_size_type,
             ExcInternalError());

      if (in.constraint_position != numbers::invalid_size_type)
        {
          constraint_position    = in.constraint_position;
          in.constraint_position = numbers::invalid_size_type;
        }
      return *this;
    }



    template <typename number>
    DataCache<number>::DataCache()
      : row_length(8)
    {}



    template <typename number>
    void
    DataCache<number>::reinit()
    {
      individual_size.resize(0);
      data.resize(0);
    }



    template <typename number>
    size_type
    DataCache<number>::insert_new_index(
      const std::pair<size_type, number> &pair)
    {
      Assert(row_length > 0, ExcInternalError());
      const unsigned int index = individual_size.size();
      individual_size.push_back(1);
      data.resize(individual_size.size() * row_length);
      data[index * row_length] = pair;
      individual_size[index]   = 1;
      return index;
    }



    template <typename number>
    void
    DataCache<number>::append_index(const size_type                     index,
                                    const std::pair<size_type, number> &pair)
    {
      AssertIndexRange(index, individual_size.size());
      const size_type my_length = individual_size[index];
      if (my_length == row_length)
        {
          AssertDimension(data.size(), individual_size.size() * row_length);
          // no space left in this row, need to double row_length and
          // rearrange the data items. Move all items to the right except the
          // first one, starting at the back. Since individual_size contains
          // at least one element when we get here, subtracting 1 works fine.
          data.resize(2 * data.size());
          for (size_type i = individual_size.size() - 1; i > 0; --i)
            {
              const auto ptr = data.data();
              std::move_backward(ptr + i * row_length,
                                 ptr + i * row_length + individual_size[i],
                                 ptr + i * 2 * row_length + individual_size[i]);
            }
          row_length *= 2;
        }
      data[index * row_length + my_length] = pair;
      individual_size[index]               = my_length + 1;
    }



    template <typename number>
    size_type
    DataCache<number>::get_size(const size_type index) const
    {
      return individual_size[index];
    }



    template <typename number>
    const std::pair<size_type, number> *
    DataCache<number>::get_entry(const size_type index) const
    {
      return &data[index * row_length];
    }



    template <typename number>
    GlobalRowsFromLocal<number>::GlobalRowsFromLocal()
      : n_active_rows(0)
      , n_inhomogeneous_rows(0)
    {}



    template <typename number>
    void
    GlobalRowsFromLocal<number>::reinit(const size_type n_local_rows)
    {
      total_row_indices.resize(n_local_rows);
      for (unsigned int i = 0; i < n_local_rows; ++i)
        total_row_indices[i].constraint_position = numbers::invalid_size_type;
      n_active_rows        = n_local_rows;
      n_inhomogeneous_rows = 0;
      data_cache.reinit();
    }



    template <typename number>
    void
    GlobalRowsFromLocal<number>::print(std::ostream &os)
    {
      os << "Active rows " << n_active_rows << std::endl
         << "Constr rows " << n_constraints() << std::endl
         << "Inhom  rows " << n_inhomogeneous_rows << std::endl
         << "Local: ";
      for (const auto &total_row_index : total_row_indices)
        os << ' ' << std::setw(4) << total_row_index.local_row;
      os << std::endl << "Global:";
      for (const auto &total_row_index : total_row_indices)
        os << ' ' << std::setw(4) << total_row_index.global_row;
      os << std::endl << "ConPos:";
      for (const auto &total_row_index : total_row_indices)
        os << ' ' << std::setw(4) << total_row_index.constraint_position;
      os << std::endl;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::size() const
    {
      return n_active_rows;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::size(const size_type counter_index) const
    {
      return (total_row_indices[counter_index].constraint_position ==
                  numbers::invalid_size_type ?
                0 :
                data_cache.get_size(
                  total_row_indices[counter_index].constraint_position));
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::global_row(const size_type counter_index) const
    {
      return total_row_indices[counter_index].global_row;
    }



    template <typename number>
    size_type &
    GlobalRowsFromLocal<number>::global_row(const size_type counter_index)
    {
      return total_row_indices[counter_index].global_row;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::local_row(const size_type counter_index) const
    {
      return total_row_indices[counter_index].local_row;
    }



    template <typename number>
    size_type &
    GlobalRowsFromLocal<number>::local_row(const size_type counter_index)
    {
      return total_row_indices[counter_index].local_row;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::local_row(
      const size_type counter_index,
      const size_type index_in_constraint) const
    {
      return (data_cache.get_entry(total_row_indices[counter_index]
                                     .constraint_position)[index_in_constraint])
        .first;
    }



    template <typename number>
    number
    GlobalRowsFromLocal<number>::constraint_value(
      const size_type counter_index,
      const size_type index_in_constraint) const
    {
      return (data_cache.get_entry(total_row_indices[counter_index]
                                     .constraint_position)[index_in_constraint])
        .second;
    }



    template <typename number>
    bool
    GlobalRowsFromLocal<number>::have_indirect_rows() const
    {
      return data_cache.individual_size.empty() == false;
    }



    template <typename number>
    void
    GlobalRowsFromLocal<number>::insert_constraint(
      const size_type constrained_local_dof)
    {
      --n_active_rows;
      total_row_indices[n_active_rows].local_row  = constrained_local_dof;
      total_row_indices[n_active_rows].global_row = numbers::invalid_size_type;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::n_constraints() const
    {
      return total_row_indices.size() - n_active_rows;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::n_inhomogeneities() const
    {
      return n_inhomogeneous_rows;
    }



    template <typename number>
    void
    GlobalRowsFromLocal<number>::set_ith_constraint_inhomogeneous(
      const size_type i)
    {
      Assert(i >= n_inhomogeneous_rows, ExcInternalError());
      std::swap(total_row_indices[n_active_rows + i],
                total_row_indices[n_active_rows + n_inhomogeneous_rows]);
      ++n_inhomogeneous_rows;
    }



    template <typename number>
    size_type
    GlobalRowsFromLocal<number>::constraint_origin(size_type i) const
    {
      return total_row_indices[n_active_rows + i].local_row;
    }


    // a function that appends an additional row to the list of values, or
    // appends a value to an already existing row. Similar functionality as for
    // std::map<size_type,Distributing>, but here done for a
    // std::vector<Distributing>, much faster for short lists as we have them
    // here
    template <typename number>
    inline void
    GlobalRowsFromLocal<number>::insert_index(const size_type global_row,
                                              const size_type local_row,
                                              const number    constraint_value)
    {
      using index_iterator = std::vector<Distributing>::iterator;
      index_iterator               pos, pos1;
      Distributing                 row_value(global_row);
      std::pair<size_type, number> constraint(local_row, constraint_value);

      // check whether the list was really sorted before entering here
      for (size_type i = 1; i < n_active_rows; ++i)
        Assert(total_row_indices[i - 1] < total_row_indices[i],
               ExcInternalError());

      pos = Utilities::lower_bound(total_row_indices.begin(),
                                   total_row_indices.begin() + n_active_rows,
                                   row_value);
      if (pos->global_row == global_row)
        pos1 = pos;
      else
        {
          pos1 = total_row_indices.insert(pos, row_value);
          ++n_active_rows;
        }

      if (pos1->constraint_position == numbers::invalid_size_type)
        pos1->constraint_position = data_cache.insert_new_index(constraint);
      else
        data_cache.append_index(pos1->constraint_position, constraint);
    }

    // this sort algorithm sorts std::vector<Distributing>, but does not take
    // the constraints into account. this means that in case that constraints
    // are already inserted, this function does not work as expected. Use
    // shellsort, which is very fast in case the indices are already sorted
    // (which is the usual case with DG elements), and not too slow in other
    // cases
    template <typename number>
    inline void
    GlobalRowsFromLocal<number>::sort()
    {
      size_type i, j, j2, temp, templ, istep;
      size_type step;

      // check whether the constraints are really empty.
      const size_type length = size();

      // make sure that we are in the range of the vector
      AssertIndexRange(length, total_row_indices.size() + 1);
      for (size_type i = 0; i < length; ++i)
        Assert(total_row_indices[i].constraint_position ==
                 numbers::invalid_size_type,
               ExcInternalError());

      step = length / 2;
      while (step > 0)
        {
          for (i = step; i < length; ++i)
            {
              istep = step;
              j     = i;
              j2    = j - istep;
              temp  = total_row_indices[i].global_row;
              templ = total_row_indices[i].local_row;
              if (total_row_indices[j2].global_row > temp)
                {
                  while ((j >= istep) &&
                         (total_row_indices[j2].global_row > temp))
                    {
                      total_row_indices[j].global_row =
                        total_row_indices[j2].global_row;
                      total_row_indices[j].local_row =
                        total_row_indices[j2].local_row;
                      j = j2;
                      j2 -= istep;
                    }
                  total_row_indices[j].global_row = temp;
                  total_row_indices[j].local_row  = templ;
                }
            }
          step = step >> 1;
        }
    }



    /**
     * This class is an accessor class to scratch data that is used
     * during calls to distribute_local_to_global and
     * add_entries_local_to_global. In order to avoid frequent memory
     * allocation, we keep the data alive from one call to the next in
     * a static variable. Since we want to allow for different number
     * types in matrices, this is a template.
     *
     * Since each thread gets its private version of scratch data out of the
     * ThreadLocalStorage, no conflicting access can occur. For this to be
     * valid, we need to make sure that no call within
     * distribute_local_to_global is made that by itself can spawn tasks.
     * Otherwise, we might end up in a situation where several threads fight for
     * the data.
     */
    template <typename number>
    class ScratchDataAccessor
    {
    public:
      /**
       * Constructor. Takes the scratch data object for the current
       * thread out of the provided object and marks it as used.
       */
      ScratchDataAccessor(
        Threads::ThreadLocalStorage<ScratchData<number>> &tls_scratch_data)
        : my_scratch_data(&tls_scratch_data.get())
      {
        Assert(my_scratch_data->in_use == false,
               ExcMessage(
                 "Access to thread-local scratch data tried, but it is already "
                 "in use"));
        my_scratch_data->in_use = true;
      }

      /**
       * Destructor. Mark scratch data as available again.
       */
      ~ScratchDataAccessor()
      {
        my_scratch_data->in_use = false;
      }

      /**
       * Dereferencing operator.
       */
      ScratchData<number> &
      operator*()
      {
        return *my_scratch_data;
      }

      /**
       * Dereferencing operator.
       */
      ScratchData<number> *
      operator->()
      {
        return my_scratch_data;
      }

    private:
      ScratchData<number> *my_scratch_data;
    };


    // function for block matrices: Find out where in the list of local dofs
    // (sorted according to global ids) the individual blocks start.
    // Transform the global indices to block-local indices in order to be
    // able to use functions like vector.block(1)(block_local_id), instead of
    // vector(global_id). This avoids transforming indices one-by-one later
    // on.
    template <typename number, typename BlockType>
    inline void
    make_block_starts(const BlockType             &block_object,
                      GlobalRowsFromLocal<number> &global_rows,
                      std::vector<size_type>      &block_starts)
    {
      AssertDimension(block_starts.size(), block_object.n_block_rows() + 1);

      using row_iterator         = std::vector<Distributing>::iterator;
      row_iterator block_indices = global_rows.total_row_indices.begin();

      const size_type num_blocks    = block_object.n_block_rows();
      const size_type n_active_rows = global_rows.size();

      // find end of rows.
      block_starts[0] = 0;
      for (size_type i = 1; i < num_blocks; ++i)
        {
          row_iterator first_block = Utilities::lower_bound(
            block_indices,
            global_rows.total_row_indices.begin() + n_active_rows,
            Distributing(block_object.get_row_indices().block_start(i)));
          block_starts[i] = first_block - global_rows.total_row_indices.begin();
          block_indices   = first_block;
        }
      block_starts[num_blocks] = n_active_rows;

      // transform row indices to block-local index space
      for (size_type i = block_starts[1]; i < n_active_rows; ++i)
        global_rows.global_row(i) =
          block_object.get_row_indices()
            .global_to_local(global_rows.global_row(i))
            .second;
    }

    // same as before, but for std::vector<uint> instead of
    // GlobalRowsFromLocal. Used in functions for sparsity patterns.
    template <typename BlockType>
    inline void
    make_block_starts(const BlockType        &block_object,
                      std::vector<size_type> &row_indices,
                      std::vector<size_type> &block_starts)
    {
      AssertDimension(block_starts.size(), block_object.n_block_rows() + 1);

      using row_iterator       = std::vector<size_type>::iterator;
      row_iterator col_indices = row_indices.begin();

      const size_type num_blocks = block_object.n_block_rows();

      // find end of rows.
      block_starts[0] = 0;
      for (size_type i = 1; i < num_blocks; ++i)
        {
          row_iterator first_block =
            Utilities::lower_bound(col_indices,
                                   row_indices.end(),
                                   block_object.get_row_indices().block_start(
                                     i));
          block_starts[i] = first_block - row_indices.begin();
          col_indices     = first_block;
        }
      block_starts[num_blocks] = row_indices.size();

      // transform row indices to local index space
      for (size_type i = block_starts[1]; i < row_indices.size(); ++i)
        row_indices[i] =
          block_object.get_row_indices().global_to_local(row_indices[i]).second;
    }

    // resolves constraints of one column at the innermost loop. goes through
    // the origin of each global entry and finds out which data we need to
    // collect.
    template <typename number>
    static inline number
    resolve_matrix_entry(const GlobalRowsFromLocal<number> &global_rows,
                         const GlobalRowsFromLocal<number> &global_cols,
                         const size_type                    i,
                         const size_type                    j,
                         const size_type                    loc_row,
                         const FullMatrix<number>          &local_matrix)
    {
      const size_type loc_col = global_cols.local_row(j);
      number          col_val;

      // case 1: row has direct contribution in local matrix. decide whether col
      // has a direct contribution. if not, set the value to zero.
      if (loc_row != numbers::invalid_size_type)
        {
          col_val = ((loc_col != numbers::invalid_size_type) ?
                       local_matrix(loc_row, loc_col) :
                       0);

          // account for indirect contributions by constraints in column
          for (size_type p = 0; p < global_cols.size(j); ++p)
            col_val += (local_matrix(loc_row, global_cols.local_row(j, p)) *
                        global_cols.constraint_value(j, p));
        }

      // case 2: row has no direct contribution in local matrix
      else
        col_val = 0;

      // account for indirect contributions by constraints in row, going through
      // the direct and indirect references in the given column.
      for (size_type q = 0; q < global_rows.size(i); ++q)
        {
          number add_this =
            (loc_col != numbers::invalid_size_type) ?
              local_matrix(global_rows.local_row(i, q), loc_col) :
              0;

          for (size_type p = 0; p < global_cols.size(j); ++p)
            add_this += (local_matrix(global_rows.local_row(i, q),
                                      global_cols.local_row(j, p)) *
                         global_cols.constraint_value(j, p));
          col_val += add_this * global_rows.constraint_value(i, q);
        }
      return col_val;
    }

    // computes all entries that need to be written into global_rows[i]. Lists
    // the resulting values in val_ptr, and the corresponding column indices in
    // col_ptr.
    template <typename number>
    inline void
    resolve_matrix_row(const GlobalRowsFromLocal<number> &global_rows,
                       const GlobalRowsFromLocal<number> &global_cols,
                       const size_type                    i,
                       const size_type                    column_start,
                       const size_type                    column_end,
                       const FullMatrix<number>          &local_matrix,
                       size_type                        *&col_ptr,
                       number                           *&val_ptr)
    {
      if (column_end == column_start)
        return;

      AssertIndexRange(column_end - 1, global_cols.size());
      const size_type loc_row = global_rows.local_row(i);

      // fast function if there are no indirect references to any of the local
      // rows at all on this set of dofs (saves a lot of checks). the only check
      // we actually need to perform is whether the matrix element is zero.
      if (global_rows.have_indirect_rows() == false &&
          global_cols.have_indirect_rows() == false)
        {
          AssertIndexRange(loc_row, local_matrix.m());
          const number *matrix_ptr = &local_matrix(loc_row, 0);

          for (size_type j = column_start; j < column_end; ++j)
            {
              const size_type loc_col = global_cols.local_row(j);
              AssertIndexRange(loc_col, local_matrix.n());
              const number col_val = matrix_ptr[loc_col];
              if (col_val != number())
                {
                  *val_ptr++ = static_cast<number>(col_val);
                  *col_ptr++ = global_cols.global_row(j);
                }
            }
        }

      // more difficult part when there are indirect references and when we need
      // to do some more checks.
      else
        {
          for (size_type j = column_start; j < column_end; ++j)
            {
              number col_val = resolve_matrix_entry(
                global_rows, global_cols, i, j, loc_row, local_matrix);

              // if we got some nontrivial value, append it to the array of
              // values.
              if (col_val != number())
                {
                  *val_ptr++ = static_cast<number>(col_val);
                  *col_ptr++ = global_cols.global_row(j);
                }
            }
        }
    }

    // specialized function that can write into the row of a
    // SparseMatrix<number>.
    namespace dealiiSparseMatrix
    {
      template <typename SparseMatrixIterator, typename LocalType>
      static inline void
      add_value(const LocalType       value,
                const size_type       row,
                const size_type       column,
                SparseMatrixIterator &matrix_values)
      {
        (void)row;
        if (value != LocalType())
          {
            while (matrix_values->column() < column)
              ++matrix_values;

            // Ensure that we haven't walked off the end of the row and
            // accidentally found the column in a later row. This should
            // not have happened since we believed we know that we enter
            // entries sorted by column, and that all columns are represented
            // in the current matrix row. But it's worth checking.
            Assert(matrix_values->row() == row, ExcInternalError());

            // Then also assert that the column index actually exists.
            Assert(matrix_values->column() == column,
                   typename SparseMatrix<
                     typename SparseMatrixIterator::MatrixType::value_type>::
                     ExcInvalidIndex(row, column));

            // Now so convinced, let us add the relevant value:
            matrix_values->value() += value;
          }
      }
    } // namespace dealiiSparseMatrix

    // similar as before, now with shortcut for deal.II sparse matrices. this
    // lets us avoid using extra arrays, and does all the operations just in
    // place, i.e., in the respective matrix row
    template <typename number>
    inline void
    resolve_matrix_row(const GlobalRowsFromLocal<number> &global_rows,
                       const size_type                    i,
                       const size_type                    column_start,
                       const size_type                    column_end,
                       const FullMatrix<number>          &local_matrix,
                       SparseMatrix<number>              *sparse_matrix)
    {
      if (column_end == column_start)
        return;

      AssertIndexRange(column_end - 1, global_rows.size());
      const SparsityPattern &sparsity = sparse_matrix->get_sparsity_pattern();

      if (sparsity.n_nonzero_elements() == 0)
        return;

      const size_type row     = global_rows.global_row(i);
      const size_type loc_row = global_rows.local_row(i);

      typename SparseMatrix<number>::iterator matrix_values =
        sparse_matrix->begin(row);
      const bool optimize_diagonal = sparsity.n_rows() == sparsity.n_cols();

      // distinguish three cases about what can happen for checking whether the
      // diagonal is the first element of the row. this avoids if statements at
      // the innermost loop positions

      if (!optimize_diagonal) // case 1: no diagonal optimization in matrix
        {
          if (global_rows.have_indirect_rows() == false)
            {
              AssertIndexRange(loc_row, local_matrix.m());
              const number *matrix_ptr = &local_matrix(loc_row, 0);

              for (size_type j = column_start; j < column_end; ++j)
                {
                  const size_type loc_col = global_rows.local_row(j);
                  const number    col_val = matrix_ptr[loc_col];
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
            }
          else
            {
              for (size_type j = column_start; j < column_end; ++j)
                {
                  number col_val = resolve_matrix_entry(
                    global_rows, global_rows, i, j, loc_row, local_matrix);
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
            }
        }
      else if (i >= column_start && i < column_end) // case 2: can split loop
        {
          ++matrix_values; // jump over diagonal element
          if (global_rows.have_indirect_rows() == false)
            {
              AssertIndexRange(loc_row, local_matrix.m());
              const number *matrix_ptr = &local_matrix(loc_row, 0);

              sparse_matrix->begin(row)->value() += matrix_ptr[loc_row];
              for (size_type j = column_start; j < i; ++j)
                {
                  const size_type loc_col = global_rows.local_row(j);
                  const number    col_val = matrix_ptr[loc_col];
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
              for (size_type j = i + 1; j < column_end; ++j)
                {
                  const size_type loc_col = global_rows.local_row(j);
                  const number    col_val = matrix_ptr[loc_col];
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
            }
          else
            {
              sparse_matrix->begin(row)->value() += resolve_matrix_entry(
                global_rows, global_rows, i, i, loc_row, local_matrix);
              for (size_type j = column_start; j < i; ++j)
                {
                  number col_val = resolve_matrix_entry(
                    global_rows, global_rows, i, j, loc_row, local_matrix);
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
              for (size_type j = i + 1; j < column_end; ++j)
                {
                  number col_val = resolve_matrix_entry(
                    global_rows, global_rows, i, j, loc_row, local_matrix);
                  dealiiSparseMatrix::add_value(col_val,
                                                row,
                                                global_rows.global_row(j),
                                                matrix_values);
                }
            }
        }
      // case 3: can't say - need to check inside the loop
      else if (global_rows.have_indirect_rows() == false)
        {
          ++matrix_values; // jump over diagonal element
          AssertIndexRange(loc_row, local_matrix.m());
          const number *matrix_ptr = &local_matrix(loc_row, 0);

          for (size_type j = column_start; j < column_end; ++j)
            {
              const size_type loc_col = global_rows.local_row(j);
              const number    col_val = matrix_ptr[loc_col];
              if (row == global_rows.global_row(j))
                sparse_matrix->begin(row)->value() += col_val;
              else
                dealiiSparseMatrix::add_value(col_val,
                                              row,
                                              global_rows.global_row(j),
                                              matrix_values);
            }
        }
      else
        {
          ++matrix_values; // jump over diagonal element
          for (size_type j = column_start; j < column_end; ++j)
            {
              number col_val = resolve_matrix_entry(
                global_rows, global_rows, i, j, loc_row, local_matrix);
              if (row == global_rows.global_row(j))
                sparse_matrix->begin(row)->value() += col_val;
              else
                dealiiSparseMatrix::add_value(col_val,
                                              row,
                                              global_rows.global_row(j),
                                              matrix_values);
            }
        }
    }

    // Same function to resolve all entries that will be added to the given
    // global row global_rows[i] as before, now for sparsity pattern
    template <typename number>
    inline void
    resolve_matrix_row(const GlobalRowsFromLocal<number> &global_rows,
                       const size_type                    i,
                       const size_type                    column_start,
                       const size_type                    column_end,
                       const Table<2, bool>              &dof_mask,
                       std::vector<size_type>::iterator  &col_ptr)
    {
      if (column_end == column_start)
        return;

      const size_type loc_row = global_rows.local_row(i);

      // fast function if there are no indirect references to any of the local
      // rows at all on this set of dofs
      if (global_rows.have_indirect_rows() == false)
        {
          Assert(loc_row < dof_mask.n_rows(), ExcInternalError());

          for (size_type j = column_start; j < column_end; ++j)
            {
              const size_type loc_col = global_rows.local_row(j);
              Assert(loc_col < dof_mask.n_cols(), ExcInternalError());

              if (dof_mask(loc_row, loc_col) == true)
                *col_ptr++ = global_rows.global_row(j);
            }
        }

      // slower functions when there are indirect references and when we need to
      // do some more checks.
      else
        {
          for (size_type j = column_start; j < column_end; ++j)
            {
              const size_type loc_col = global_rows.local_row(j);
              if (loc_row != numbers::invalid_size_type)
                {
                  Assert(loc_row < dof_mask.n_rows(), ExcInternalError());
                  if (loc_col != numbers::invalid_size_type)
                    {
                      Assert(loc_col < dof_mask.n_cols(), ExcInternalError());
                      if (dof_mask(loc_row, loc_col) == true)
                        goto add_this_index;
                    }

                  for (size_type p = 0; p < global_rows.size(j); ++p)
                    if (dof_mask(loc_row, global_rows.local_row(j, p)) == true)
                      goto add_this_index;
                }

              for (size_type q = 0; q < global_rows.size(i); ++q)
                {
                  if (loc_col != numbers::invalid_size_type)
                    {
                      Assert(loc_col < dof_mask.n_cols(), ExcInternalError());
                      if (dof_mask(global_rows.local_row(i, q), loc_col) ==
                          true)
                        goto add_this_index;
                    }

                  for (size_type p = 0; p < global_rows.size(j); ++p)
                    if (dof_mask(global_rows.local_row(i, q),
                                 global_rows.local_row(j, p)) == true)
                      goto add_this_index;
                }

              continue;
              // if we got some nontrivial value, append it to the array of
              // values.
            add_this_index:
              *col_ptr++ = global_rows.global_row(j);
            }
        }
    }

    // to make sure that the global matrix remains invertible, we need to do
    // something with the diagonal elements. Add the average of the
    // absolute values of the local matrix diagonals, so the resulting entry
    // will always be positive and furthermore be in the same order of magnitude
    // as the other elements of the matrix. If all local matrix diagonals are
    // zero, add the l1 norm of the local matrix divided by the matrix size
    // to the diagonal of the global matrix. If the entire local matrix is zero,
    // add 1 to the diagonal of the global matrix.
    //
    // note that this also captures the special case that a dof is both
    // constrained and fixed (this can happen for hanging nodes in 3d that also
    // happen to be on the boundary). in that case, following the program flow
    // in distribute_local_to_global, it is realized that when distributing the
    // row and column no elements of the matrix are actually touched if all the
    // degrees of freedom to which this dof is constrained are also constrained
    // (the usual case with hanging nodes in 3d). however, in the line below, we
    // do actually do something with this dof
    template <typename number, typename MatrixType, typename VectorType>
    inline void
    set_matrix_diagonals(
      const internal::AffineConstraints::GlobalRowsFromLocal<number>
                                              &global_rows,
      const std::vector<size_type>            &local_dof_indices,
      const FullMatrix<number>                &local_matrix,
      const dealii::AffineConstraints<number> &constraints,
      MatrixType                              &global_matrix,
      VectorType                              &global_vector,
      bool                                     use_inhomogeneities_for_rhs)
    {
      if (global_rows.n_constraints() > 0)
        {
          number average_diagonal = number();
          for (size_type i = 0; i < local_matrix.m(); ++i)
            average_diagonal += std::abs(local_matrix(i, i));
          average_diagonal /= static_cast<number>(local_matrix.m());

          // handle the case that all diagonal elements are zero
          if (average_diagonal == static_cast<number>(0.))
            {
              average_diagonal = static_cast<number>(local_matrix.l1_norm()) /
                                 static_cast<number>(local_matrix.m());
              // if the entire matrix is zero, use 1. for the diagonal
              if (average_diagonal == static_cast<number>(0.))
                average_diagonal = static_cast<number>(1.);
            }

          for (size_type i = 0; i < global_rows.n_constraints(); ++i)
            {
              const size_type local_row  = global_rows.constraint_origin(i);
              const size_type global_row = local_dof_indices[local_row];

              const number current_diagonal =
                local_matrix(local_row, local_row);
              if (std::abs(current_diagonal) != 0.)
                {
                  global_matrix.add(global_row,
                                    global_row,
                                    std::abs(current_diagonal));
                  // if the use_inhomogeneities_for_rhs flag is set to true, the
                  // inhomogeneities are used to create the global vector.
                  // instead of fill in a zero in the ith components with an
                  // inhomogeneity, we set those to:
                  // inhomogeneity(i)*global_matrix (i,i).
                  if (use_inhomogeneities_for_rhs == true)
                    global_vector(global_row) +=
                      current_diagonal *
                      constraints.get_inhomogeneity(global_row);
                }
              else
                {
                  global_matrix.add(global_row, global_row, average_diagonal);

                  if (use_inhomogeneities_for_rhs == true)
                    global_vector(global_row) +=
                      average_diagonal *
                      constraints.get_inhomogeneity(global_row);
                }
            }
        }
    }

    // similar function as the one above for setting matrix diagonals, but now
    // doing that for sparsity patterns when setting them up using
    // add_entries_local_to_global(). In case we keep constrained entries, add
    // all the rows and columns related to the constrained dof, otherwise just
    // add the diagonal
    template <typename number>
    inline void
    set_sparsity_diagonals(
      const internal::AffineConstraints::GlobalRowsFromLocal<number>
                                   &global_rows,
      const std::vector<size_type> &local_dof_indices,
      const Table<2, bool>         &dof_mask,
      const bool                    keep_constrained_entries,
      ScratchData<number>          &scratch_data,
      SparsityPatternBase          &sparsity_pattern)
    {
      // if we got constraints, need to add the diagonal element and, if the
      // user requested so, also the rest of the entries in rows and columns
      // that have been left out above
      if (global_rows.n_constraints() > 0)
        {
          std::vector<std::pair<size_type, size_type>> &cell_entries =
            scratch_data.new_entries;
          cell_entries.resize(0);
          cell_entries.reserve(local_dof_indices.size());
          for (size_type i = 0; i < global_rows.n_constraints(); ++i)
            {
              const size_type local_row  = global_rows.constraint_origin(i);
              const size_type global_row = local_dof_indices[local_row];
              if (keep_constrained_entries == true)
                {
                  for (size_type j = 0; j < local_dof_indices.size(); ++j)
                    {
                      if (dof_mask(local_row, j) == true)
                        cell_entries.emplace_back(global_row,
                                                  local_dof_indices[j]);
                      if (dof_mask(j, local_row) == true)
                        cell_entries.emplace_back(local_dof_indices[j],
                                                  global_row);
                    }
                }
              else
                // don't keep constrained entries - just add the diagonal.
                cell_entries.emplace_back(global_row, global_row);
            }
          sparsity_pattern.add_entries(make_array_view(cell_entries));
        }
    }

  } // end of namespace AffineConstraints
} // end of namespace internal



// Basic idea of setting up a list of
// all global dofs: first find all rows and columns
// that we are going to write touch,
// and then go through the
// lines and collect all the local rows that
// are related to it.
template <typename number>
void
AffineConstraints<number>::make_sorted_row_list(
  const std::vector<size_type>                             &local_dof_indices,
  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows) const
{
  const size_type n_local_dofs = local_dof_indices.size();
  AssertDimension(n_local_dofs, global_rows.size());

  // when distributing the local data to the global matrix, we can quite
  // cheaply sort the indices (obviously, this introduces the need for
  // allocating some memory on the way, but we need to do this only for rows,
  // whereas the distribution process itself goes over rows and columns). This
  // has the advantage that when writing into the global matrix, we can make
  // use of the sortedness.

  // so the first step is to create a sorted list of all row values that are
  // possible. these values are either the rows from unconstrained dofs, or
  // some indices introduced by dofs constrained to a combination of some
  // other dofs. regarding the data type, choose a <tt>std::vector</tt> of a
  // pair of unsigned ints (for global columns) and internal data (containing
  // local columns + possible jumps from constraints). Choosing
  // <tt>std::map</tt> or anything else M.K. knows of would be much more
  // expensive here!

  // cache whether we have to resolve any indirect rows generated from
  // resolving constrained dofs.
  size_type added_rows = 0;

  // first add the indices in an unsorted way and only keep track of the
  // constraints that appear. They are resolved in a second step.
  for (size_type i = 0; i < n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
        {
          global_rows.global_row(added_rows)  = local_dof_indices[i];
          global_rows.local_row(added_rows++) = i;
        }
      else
        global_rows.insert_constraint(i);
    }
  global_rows.sort();

  const size_type n_constrained_rows = n_local_dofs - added_rows;
  for (size_type i = 0; i < n_constrained_rows; ++i)
    {
      const size_type local_row = global_rows.constraint_origin(i);
      AssertIndexRange(local_row, n_local_dofs);
      const size_type global_row = local_dof_indices[local_row];
      Assert(is_constrained(global_row), ExcInternalError());
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(global_row)]];
      if (position.inhomogeneity != number(0.))
        global_rows.set_ith_constraint_inhomogeneous(i);
      for (size_type q = 0; q < position.entries.size(); ++q)
        global_rows.insert_index(position.entries[q].first,
                                 local_row,
                                 position.entries[q].second);
    }
}



// Same function as before, but now do only extract the global indices that
// come from the local ones without storing their origin. Used for sparsity
// pattern generation.
template <typename number>
inline void
AffineConstraints<number>::make_sorted_row_list(
  const std::vector<size_type> &local_dof_indices,
  std::vector<size_type>       &active_dofs) const
{
  const size_type n_local_dofs = local_dof_indices.size();
  size_type       added_rows   = 0;
  for (size_type i = 0; i < n_local_dofs; ++i)
    {
      if (is_constrained(local_dof_indices[i]) == false)
        {
          active_dofs[added_rows++] = local_dof_indices[i];
          continue;
        }

      active_dofs[n_local_dofs - i + added_rows - 1] = i;
    }
  std::sort(active_dofs.begin(), active_dofs.begin() + added_rows);

  const size_type n_constrained_dofs = n_local_dofs - added_rows;
  for (size_type i = n_constrained_dofs; i > 0; --i)
    {
      const size_type local_row = active_dofs.back();

      // remove constrained entry since we are going to resolve it in place
      active_dofs.pop_back();
      const size_type       global_row = local_dof_indices[local_row];
      const ConstraintLine &position =
        lines[lines_cache[calculate_line_index(global_row)]];
      for (size_type q = 0; q < position.entries.size(); ++q)
        {
          const size_type new_index = position.entries[q].first;
          // in case all dofs are constrained, we might insert at
          // active_dofs.begin(), but we should never insert before that
          AssertIndexRange(i - 1, active_dofs.size() + 1);
          if (i > active_dofs.size() ||
              active_dofs[active_dofs.size() - i] < new_index)
            active_dofs.insert(active_dofs.end() - i + 1, new_index);

          // make binary search to find where to put the new index in order to
          // keep the list sorted
          else
            {
              std::vector<size_type>::iterator it =
                Utilities::lower_bound(active_dofs.begin(),
                                       active_dofs.end() - i + 1,
                                       new_index);
              if (*it != new_index)
                active_dofs.insert(it, new_index);
            }
        }
    }
}



// Resolve the constraints from the vector and apply inhomogeneities.
template <typename number>
template <typename MatrixScalar, typename VectorScalar>
inline typename ProductType<VectorScalar, MatrixScalar>::type
AffineConstraints<number>::resolve_vector_entry(
  const size_type                                                 i,
  const internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows,
  const Vector<VectorScalar>                                     &local_vector,
  const std::vector<size_type>   &local_dof_indices,
  const FullMatrix<MatrixScalar> &local_matrix) const
{
  const size_type loc_row              = global_rows.local_row(i);
  const size_type n_inhomogeneous_rows = global_rows.n_inhomogeneities();
  typename ProductType<VectorScalar, MatrixScalar>::type val = 0;
  // has a direct contribution from some local entry. If we have inhomogeneous
  // constraints, compute the contribution of the inhomogeneity in the current
  // row.
  if (loc_row != numbers::invalid_size_type)
    {
      val = local_vector(loc_row);
      for (size_type i = 0; i < n_inhomogeneous_rows; ++i)
        val -= (local_matrix(loc_row, global_rows.constraint_origin(i)) *
                lines[lines_cache[calculate_line_index(
                        local_dof_indices[global_rows.constraint_origin(i)])]]
                  .inhomogeneity);
    }

  // go through the indirect contributions
  for (size_type q = 0; q < global_rows.size(i); ++q)
    {
      const size_type loc_row_q = global_rows.local_row(i, q);
      typename ProductType<VectorScalar, MatrixScalar>::type add_this =
        local_vector(loc_row_q);
      for (size_type k = 0; k < n_inhomogeneous_rows; ++k)
        add_this -=
          (local_matrix(loc_row_q, global_rows.constraint_origin(k)) *
           lines[lines_cache[calculate_line_index(
                   local_dof_indices[global_rows.constraint_origin(k)])]]
             .inhomogeneity);
      val += add_this * global_rows.constraint_value(i, q);
    }
  return val;
}



// internal implementation for distribute_local_to_global for standard
// (non-block) matrices
template <typename number>
template <typename MatrixType, typename VectorType>
void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number>     &local_matrix,
  const Vector<number>         &local_vector,
  const std::vector<size_type> &local_dof_indices,
  MatrixType                   &global_matrix,
  VectorType                   &global_vector,
  const bool                    use_inhomogeneities_for_rhs,
  const std::bool_constant<false>) const
{
  // FIXME: static_assert MatrixType::value_type == number

  // check whether we work on real vectors or we just used a dummy when
  // calling the other function above.
  const bool use_vectors =
    (local_vector.size() == 0 && global_vector.size() == 0) ? false : true;
  const bool use_dealii_matrix =
    std::is_same_v<MatrixType, SparseMatrix<number>>;

  AssertDimension(local_matrix.n(), local_dof_indices.size());
  AssertDimension(local_matrix.m(), local_dof_indices.size());
  Assert(global_vector.has_ghost_elements() == false, ExcGhostsPresent());
  Assert(global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  if (use_vectors == true)
    {
      AssertDimension(local_matrix.m(), local_vector.size());
      AssertDimension(global_matrix.m(), global_vector.size());
    }
  Assert(lines.empty() || sorted == true, ExcMatrixNotClosed());

  const size_type n_local_dofs = local_dof_indices.size();

  typename internal::AffineConstraints::ScratchDataAccessor<number>
    scratch_data(this->scratch_data);

  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows =
    scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);
  make_sorted_row_list(local_dof_indices, global_rows);

  const size_type n_actual_dofs = global_rows.size();

  // create arrays for the column data (indices and values) that will then be
  // written into the matrix. Shortcut for deal.II sparse matrix. We can use
  // the scratch data if we have a double matrix. Otherwise, we need to create
  // an array in any case since we cannot know about the actual data type in
  // the AffineConstraints class (unless we do cast). This involves a little
  // bit of logic to determine the type of the matrix value.
  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>    &vals = scratch_data->values;
  // create arrays for writing into the vector as well
  std::vector<size_type> &vector_indices = scratch_data->vector_indices;
  std::vector<typename VectorType::value_type> &vector_values =
    scratch_data->vector_values;
  vector_indices.resize(n_actual_dofs);
  vector_values.resize(n_actual_dofs);
  SparseMatrix<number> *sparse_matrix =
    dynamic_cast<SparseMatrix<number> *>(&global_matrix);
  if (use_dealii_matrix == false)
    {
      cols.resize(n_actual_dofs);
      vals.resize(n_actual_dofs);
    }
  else
    Assert(sparse_matrix != nullptr, ExcInternalError());

  // now do the actual job. go through all the global rows that we will touch
  // and call resolve_matrix_row for each of those.
  size_type local_row_n = 0;
  for (size_type i = 0; i < n_actual_dofs; ++i)
    {
      const size_type row = global_rows.global_row(i);

      // calculate all the data that will be written into the matrix row.
      if (use_dealii_matrix == false)
        {
          size_type *col_ptr = cols.data();
          // cast is uncritical here and only used to avoid compiler
          // warnings. We never access a non-double array
          number *val_ptr = vals.data();
          internal::AffineConstraints::resolve_matrix_row(global_rows,
                                                          global_rows,
                                                          i,
                                                          0,
                                                          n_actual_dofs,
                                                          local_matrix,
                                                          col_ptr,
                                                          val_ptr);
          const size_type n_values = col_ptr - cols.data();
          if (n_values > 0)
            global_matrix.add(row,
                              n_values,
                              cols.data(),
                              vals.data(),
                              /* elide zero additions */ false,
                              /* sorted by column index */ true);
        }
      else
        internal::AffineConstraints::resolve_matrix_row(
          global_rows, i, 0, n_actual_dofs, local_matrix, sparse_matrix);

      // now to the vectors. besides doing the same job as we did above (i.e.,
      // distribute the content of the local vector into the global one), need
      // to account for inhomogeneities here: this corresponds to eliminating
      // the respective column in the local matrix with value on the right
      // hand side.
      if (use_vectors == true)
        {
          const typename VectorType::value_type val = resolve_vector_entry(
            i, global_rows, local_vector, local_dof_indices, local_matrix);
          AssertIsFinite(val);

          if (val != typename VectorType::value_type())
            {
              vector_indices[local_row_n] = row;
              vector_values[local_row_n]  = val;
              ++local_row_n;
            }
        }
    }
  // Drop the elements of vector_indices and vector_values that we do not use
  // (we may always elide writing zero values to vectors)
  const size_type n_local_rows = local_row_n;
  vector_indices.resize(n_local_rows);
  vector_values.resize(n_local_rows);

  // While the standard case is that these types are equal, they need not be, so
  // only do a bulk update if they are. Note that the types in the arguments to
  // add must be equal if we have a Trilinos or PETSc vector but do not have to
  // be if we have a deal.II native vector: one could further optimize this for
  // Vector, LinearAlgebra::distributed::vector, etc.
  if constexpr (std::is_same_v<typename VectorType::value_type, number>)
    {
      global_vector.add(vector_indices, vector_values);
    }
  else
    {
      for (size_type row_n = 0; row_n < n_local_rows; ++row_n)
        {
          global_vector(vector_indices[row_n]) +=
            static_cast<typename VectorType::value_type>(vector_values[row_n]);
        }
    }

  internal::AffineConstraints::set_matrix_diagonals(
    global_rows,
    local_dof_indices,
    local_matrix,
    *this,
    global_matrix,
    global_vector,
    use_inhomogeneities_for_rhs);
}



// similar function as above, but now specialized for block matrices. See the
// other function for additional comments.
template <typename number>
template <typename MatrixType, typename VectorType>
void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number>     &local_matrix,
  const Vector<number>         &local_vector,
  const std::vector<size_type> &local_dof_indices,
  MatrixType                   &global_matrix,
  VectorType                   &global_vector,
  const bool                    use_inhomogeneities_for_rhs,
  const std::bool_constant<true>) const
{
  const bool use_vectors =
    (local_vector.size() == 0 && global_vector.size() == 0) ? false : true;
  const bool use_dealii_matrix =
    std::is_same_v<MatrixType, BlockSparseMatrix<number>>;

  AssertDimension(local_matrix.n(), local_dof_indices.size());
  AssertDimension(local_matrix.m(), local_dof_indices.size());
  Assert(global_vector.has_ghost_elements() == false, ExcGhostsPresent());
  Assert(global_matrix.m() == global_matrix.n(), ExcNotQuadratic());
  Assert(global_matrix.n_block_rows() == global_matrix.n_block_cols(),
         ExcNotQuadratic());
  if (use_vectors == true)
    {
      AssertDimension(local_matrix.m(), local_vector.size());
      AssertDimension(global_matrix.m(), global_vector.size());
    }
  Assert(sorted == true, ExcMatrixNotClosed());

  typename internal::AffineConstraints::ScratchDataAccessor<number>
    scratch_data(this->scratch_data);

  const size_type n_local_dofs = local_dof_indices.size();
  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows =
    scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);

  make_sorted_row_list(local_dof_indices, global_rows);
  const size_type n_actual_dofs = global_rows.size();

  std::vector<size_type> &global_indices = scratch_data->vector_indices;
  if (use_vectors == true)
    {
      global_indices.resize(n_actual_dofs);
      for (size_type i = 0; i < n_actual_dofs; ++i)
        global_indices[i] = global_rows.global_row(i);
    }

  // additional construct that also takes care of block indices.
  const size_type         num_blocks   = global_matrix.n_block_rows();
  std::vector<size_type> &block_starts = scratch_data->block_starts;
  block_starts.resize(num_blocks + 1);
  internal::AffineConstraints::make_block_starts(global_matrix,
                                                 global_rows,
                                                 block_starts);

  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>    &vals = scratch_data->values;
  if (use_dealii_matrix == false)
    {
      cols.resize(n_actual_dofs);
      vals.resize(n_actual_dofs);
    }

  // the basic difference to the non-block variant from now onwards is that we
  // go through the blocks of the matrix separately, which allows us to set
  // the block entries individually
  for (size_type block = 0; block < num_blocks; ++block)
    {
      const size_type next_block = block_starts[block + 1];
      for (size_type i = block_starts[block]; i < next_block; ++i)
        {
          const size_type row = global_rows.global_row(i);

          for (size_type block_col = 0; block_col < num_blocks; ++block_col)
            {
              const size_type start_block = block_starts[block_col],
                              end_block   = block_starts[block_col + 1];
              if (use_dealii_matrix == false)
                {
                  size_type *col_ptr = cols.data();
                  number    *val_ptr = vals.data();
                  internal::AffineConstraints::resolve_matrix_row(global_rows,
                                                                  global_rows,
                                                                  i,
                                                                  start_block,
                                                                  end_block,
                                                                  local_matrix,
                                                                  col_ptr,
                                                                  val_ptr);
                  const size_type n_values = col_ptr - cols.data();
                  if (n_values > 0)
                    global_matrix.block(block, block_col)
                      .add(
                        row, n_values, cols.data(), vals.data(), false, true);
                }
              else
                {
                  SparseMatrix<number> *sparse_matrix =
                    dynamic_cast<SparseMatrix<number> *>(
                      &global_matrix.block(block, block_col));
                  Assert(sparse_matrix != nullptr, ExcInternalError());
                  internal::AffineConstraints::resolve_matrix_row(
                    global_rows,
                    i,
                    start_block,
                    end_block,
                    local_matrix,
                    sparse_matrix);
                }
            }

          if (use_vectors == true)
            {
              const number val = resolve_vector_entry(
                i, global_rows, local_vector, local_dof_indices, local_matrix);

              if (val != number())
                global_vector(global_indices[i]) +=
                  static_cast<typename VectorType::value_type>(val);
            }
        }
    }

  internal::AffineConstraints::set_matrix_diagonals(
    global_rows,
    local_dof_indices,
    local_matrix,
    *this,
    global_matrix,
    global_vector,
    use_inhomogeneities_for_rhs);
}



template <typename number>
template <typename MatrixType>
void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number>     &local_matrix,
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices,
  MatrixType                   &global_matrix) const
{
  distribute_local_to_global(
    local_matrix, row_indices, *this, col_indices, global_matrix);
}



template <typename number>
template <typename MatrixType>
void
AffineConstraints<number>::distribute_local_to_global(
  const FullMatrix<number>        &local_matrix,
  const std::vector<size_type>    &row_indices,
  const AffineConstraints<number> &col_constraint_matrix,
  const std::vector<size_type>    &col_indices,
  MatrixType                      &global_matrix) const
{
  AssertDimension(local_matrix.m(), row_indices.size());
  AssertDimension(local_matrix.n(), col_indices.size());

  const size_type n_local_row_dofs = row_indices.size();
  const size_type n_local_col_dofs = col_indices.size();

  typename internal::AffineConstraints::ScratchDataAccessor<
    typename MatrixType::value_type>
    scratch_data(this->scratch_data);

  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows =
    scratch_data->global_rows;
  global_rows.reinit(n_local_row_dofs);

  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_cols =
    scratch_data->global_columns;
  global_cols.reinit(n_local_col_dofs);

  make_sorted_row_list(row_indices, global_rows);
  col_constraint_matrix.make_sorted_row_list(col_indices, global_cols);

  const size_type n_actual_row_dofs = global_rows.size();
  const size_type n_actual_col_dofs = global_cols.size();

  // create arrays for the column data (indices and values) that will then be
  // written into the matrix. Shortcut for deal.II sparse matrix
  std::vector<size_type> &cols = scratch_data->columns;
  std::vector<number>    &vals = scratch_data->values;
  cols.resize(n_actual_col_dofs);
  vals.resize(n_actual_col_dofs);

  // now do the actual job.
  for (size_type i = 0; i < n_actual_row_dofs; ++i)
    {
      const size_type row = global_rows.global_row(i);

      // calculate all the data that will be written into the matrix row.
      size_type *col_ptr = cols.data();
      number    *val_ptr = vals.data();
      internal::AffineConstraints::resolve_matrix_row(global_rows,
                                                      global_cols,
                                                      i,
                                                      0,
                                                      n_actual_col_dofs,
                                                      local_matrix,
                                                      col_ptr,
                                                      val_ptr);
      const size_type n_values = col_ptr - cols.data();
      if (n_values > 0)
        global_matrix.add(row, n_values, cols.data(), vals.data(), false, true);
    }
}



template <typename number>
void
AffineConstraints<number>::add_entries_local_to_global(
  const std::vector<size_type> &local_dof_indices,
  SparsityPatternBase          &sparsity_pattern,
  const bool                    keep_constrained_entries,
  const Table<2, bool>         &dof_mask) const
{
  Assert(sparsity_pattern.n_rows() == sparsity_pattern.n_cols(),
         ExcNotQuadratic());

  const size_type n_local_dofs = local_dof_indices.size();
  typename internal::AffineConstraints::ScratchDataAccessor<number>
    scratch_data(this->scratch_data);

  const bool dof_mask_is_active = (dof_mask.n_rows() == n_local_dofs);
  if (dof_mask_is_active == true)
    {
      AssertDimension(dof_mask.n_cols(), n_local_dofs);
    }
  else
    {
      // if the dof mask is not active, all we have to do is to add some indices
      // in a matrix format. To do this, we first create an array of all the
      // indices that are to be added. these indices are the local dof indices
      // plus some indices that come from constraints.
      std::vector<size_type> &actual_dof_indices = scratch_data->columns;
      actual_dof_indices.resize(n_local_dofs);
      make_sorted_row_list(local_dof_indices, actual_dof_indices);
      const size_type n_actual_dofs = actual_dof_indices.size();

      // now add the indices we collected above to the sparsity pattern. Very
      // easy here - just add the same array to all the rows...
      for (size_type i = 0; i < n_actual_dofs; ++i)
        sparsity_pattern.add_row_entries(actual_dof_indices[i],
                                         make_array_view(actual_dof_indices),
                                         true);

      // need to add the whole row and column structure in case we keep
      // constrained entries.
      std::vector<std::pair<size_type, size_type>> &cell_entries =
        scratch_data->new_entries;
      cell_entries.resize(0);
      cell_entries.reserve(n_local_dofs);
      for (size_type i = 0; i < n_local_dofs; ++i)
        if (is_constrained(local_dof_indices[i]))
          {
            if (keep_constrained_entries == true)
              for (size_type j = 0; j < n_local_dofs; ++j)
                {
                  cell_entries.emplace_back(local_dof_indices[i],
                                            local_dof_indices[j]);
                  cell_entries.emplace_back(local_dof_indices[j],
                                            local_dof_indices[i]);
                }
            else
              {
                cell_entries.emplace_back(local_dof_indices[i],
                                          local_dof_indices[i]);
              }
          }
      sparsity_pattern.add_entries(make_array_view(cell_entries));

      return;
    }

  // complicated case: we need to filter out some indices. then the function
  // gets similar to the function for distributing matrix entries, see there
  // for additional comments.
  internal::AffineConstraints::GlobalRowsFromLocal<number> &global_rows =
    scratch_data->global_rows;
  global_rows.reinit(n_local_dofs);
  make_sorted_row_list(local_dof_indices, global_rows);
  const size_type n_actual_dofs = global_rows.size();

  // create arrays for the column indices that will then be written into the
  // sparsity pattern.
  std::vector<size_type> &cols = scratch_data->columns;
  cols.resize(n_actual_dofs);

  for (size_type i = 0; i < n_actual_dofs; ++i)
    {
      std::vector<size_type>::iterator col_ptr = cols.begin();
      const size_type                  row     = global_rows.global_row(i);
      internal::AffineConstraints::resolve_matrix_row(
        global_rows, i, 0, n_actual_dofs, dof_mask, col_ptr);

      // finally, write all the information that accumulated under the given
      // process into the global matrix row and into the vector
      if (col_ptr != cols.begin())
        sparsity_pattern.add_row_entries(row,
                                         make_array_view(cols.begin(), col_ptr),
                                         true);
    }
  internal::AffineConstraints::set_sparsity_diagonals(global_rows,
                                                      local_dof_indices,
                                                      dof_mask,
                                                      keep_constrained_entries,
                                                      *scratch_data,
                                                      sparsity_pattern);
}



template <typename number>
void
AffineConstraints<number>::add_entries_local_to_global(
  const std::vector<size_type>    &row_indices,
  const AffineConstraints<number> &col_constraints,
  const std::vector<size_type>    &col_indices,
  SparsityPatternBase             &sparsity_pattern,
  const bool                       keep_constrained_entries,
  const Table<2, bool>            &dof_mask) const
{
  const size_type n_local_rows = row_indices.size();
  const size_type n_local_cols = col_indices.size();

  // Early return if the length of row and column indices is zero, relevant
  // for the usage with FE_Nothing.
  if (n_local_cols == 0 && n_local_rows == 0)
    return;

  typename internal::AffineConstraints::ScratchDataAccessor<number>
    scratch_data(this->scratch_data);
  std::vector<std::pair<size_type, size_type>> &cell_entries =
    scratch_data->new_entries;
  cell_entries.resize(0);
  cell_entries.reserve(n_local_rows * n_local_cols);

  // if constrained entries should be kept, need to add rows and columns of
  // those to the sparsity pattern
  if (keep_constrained_entries == true)
    {
      for (const size_type row_index : row_indices)
        if (is_constrained(row_index))
          for (const size_type col_index : col_indices)
            cell_entries.emplace_back(row_index, col_index);
      for (const size_type col_index : col_indices)
        if (col_constraints.is_constrained(col_index))
          for (const size_type row_index : row_indices)
            cell_entries.emplace_back(row_index, col_index);
      sparsity_pattern.add_entries(make_array_view(cell_entries));
    }

  // if the dof mask is not active, all we have to do is to add some indices
  // in a matrix format. To do this, we first create an array of all the
  // indices that are to be added. these indices are the local dof indices
  // plus some indices that come from constraints.
  const bool dof_mask_is_active =
    dof_mask.n_rows() == n_local_rows && dof_mask.n_cols() == n_local_cols;
  if (dof_mask_is_active == false)
    {
      std::vector<size_type> &rows = scratch_data->rows;
      std::vector<size_type> &cols = scratch_data->columns;
      rows.resize(n_local_rows);
      cols.resize(n_local_cols);
      // TODO these fills may not be necessary: previously we assumed all zeros
      // which seems incorrect. At least this way things will crash
      std::fill(rows.begin(),
                rows.end(),
                std::numeric_limits<size_type>::max());
      std::fill(cols.begin(),
                cols.end(),
                std::numeric_limits<size_type>::max());
      make_sorted_row_list(row_indices, rows);
      col_constraints.make_sorted_row_list(col_indices, cols);
      const size_type n_actual_rows = rows.size();

      // now add the indices we collected above to the sparsity pattern. Very
      // easy here - just add the same array to all the rows...
      for (size_type i = 0; i < n_actual_rows; ++i)
        sparsity_pattern.add_row_entries(rows[i], make_array_view(cols), true);
      return;
    }

  // TODO: implement this
  DEAL_II_NOT_IMPLEMENTED();
}



template <typename number>
void
AffineConstraints<number>::add_entries_local_to_global(
  const std::vector<size_type> &row_indices,
  const std::vector<size_type> &col_indices,
  SparsityPatternBase          &sparsity_pattern,
  const bool                    keep_constrained_entries,
  const Table<2, bool>         &dof_mask) const
{
  // Call the function with the same name that takes a column constraint as well
  add_entries_local_to_global(row_indices,
                              *this,
                              col_indices,
                              sparsity_pattern,
                              keep_constrained_entries,
                              dof_mask);
}

DEAL_II_NAMESPACE_CLOSE

#endif
