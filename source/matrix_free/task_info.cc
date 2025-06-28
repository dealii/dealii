// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <deal.II/matrix_free/task_info.h>
#include <deal.II/matrix_free/util.h>


#ifdef DEAL_II_WITH_TBB
#  include <tbb/blocked_range.h>
#  include <tbb/parallel_for.h>
#  include <tbb/task.h>
#  ifndef DEAL_II_TBB_WITH_ONEAPI
#    include <tbb/task_scheduler_init.h>
#  endif
#endif

#include <iostream>
#include <set>

//
// TBB with oneAPI API has deprecated and removed the
// <code>tbb::tasks</code> backend. With this it is no longer possible to
// compile the following code that builds a directed acyclic graph (DAG) of
// (thread parallel) tasks without a major porting effort. It turned out
// that such a dynamic handling of dependencies and structures is not as
// competitive as initially assumed. Consequently, this part of the matrix
// free infrastructure has seen less attention than the rest over the last
// years and is (presumably) not used that often.
//
// In case of detected oneAPI backend we simply disable threading in the
// matrix free backend for now.
//
// Matthias Maier, Martin Kronbichler, 2021
//

DEAL_II_NAMESPACE_OPEN



/*-------------------- Implementation of the matrix-free loop --------------*/
namespace internal
{
  namespace MatrixFreeFunctions
  {
#if defined(DEAL_II_WITH_TBB) && !defined(DEAL_II_TBB_WITH_ONEAPI)

    // This defines the TBB data structures that are needed to schedule the
    // partition-partition variant

    namespace partition
    {
      class ActualCellWork
      {
      public:
        ActualCellWork(MFWorkerInterface **worker_pointer,
                       const unsigned int  partition,
                       const TaskInfo     &task_info)
          : worker(nullptr)
          , worker_pointer(worker_pointer)
          , partition(partition)
          , task_info(task_info)
        {}

        ActualCellWork(MFWorkerInterface &worker,
                       const unsigned int partition,
                       const TaskInfo    &task_info)
          : worker(&worker)
          , worker_pointer(nullptr)
          , partition(partition)
          , task_info(task_info)
        {}

        void
        operator()() const
        {
          MFWorkerInterface *used_worker =
            worker != nullptr ? worker : *worker_pointer;
          Assert(used_worker != nullptr, ExcInternalError());
          used_worker->cell(partition);

          if (task_info.face_partition_data.empty() == false)
            {
              used_worker->face(partition);
              used_worker->boundary(partition);
            }
        }

      private:
        MFWorkerInterface  *worker;
        MFWorkerInterface **worker_pointer;
        const unsigned int  partition;
        const TaskInfo     &task_info;
      };

      class CellWork : public tbb::task
      {
      public:
        CellWork(MFWorkerInterface &worker,
                 const unsigned int partition,
                 const TaskInfo    &task_info,
                 const bool         is_blocked)
          : dummy(nullptr)
          , work(worker, partition, task_info)
          , is_blocked(is_blocked)
        {}

        tbb::task *
        execute() override
        {
          work();

          if (is_blocked == true)
            tbb::empty_task::spawn(*dummy);
          return nullptr;
        }

        tbb::empty_task *dummy;

      private:
        ActualCellWork work;
        const bool     is_blocked;
      };



      class PartitionWork : public tbb::task
      {
      public:
        PartitionWork(MFWorkerInterface &function_in,
                      const unsigned int partition_in,
                      const TaskInfo    &task_info_in,
                      const bool         is_blocked_in = false)
          : dummy(nullptr)
          , function(function_in)
          , partition(partition_in)
          , task_info(task_info_in)
          , is_blocked(is_blocked_in)
        {}

        tbb::task *
        execute() override
        {
          tbb::empty_task *root =
            new (tbb::task::allocate_root()) tbb::empty_task;
          const unsigned int evens = task_info.partition_evens[partition];
          const unsigned int odds  = task_info.partition_odds[partition];
          const unsigned int n_blocked_workers =
            task_info.partition_n_blocked_workers[partition];
          const unsigned int n_workers =
            task_info.partition_n_workers[partition];
          std::vector<CellWork *> worker(n_workers);
          std::vector<CellWork *> blocked_worker(n_blocked_workers);

          root->set_ref_count(evens + 1);
          for (unsigned int j = 0; j < evens; ++j)
            {
              worker[j] = new (root->allocate_child())
                CellWork(function,
                         task_info.partition_row_index[partition] + 2 * j,
                         task_info,
                         false);
              if (j > 0)
                {
                  worker[j]->set_ref_count(2);
                  blocked_worker[j - 1]->dummy =
                    new (worker[j]->allocate_child()) tbb::empty_task;
                  tbb::task::spawn(*blocked_worker[j - 1]);
                }
              else
                worker[j]->set_ref_count(1);
              if (j < evens - 1)
                {
                  blocked_worker[j] = new (worker[j]->allocate_child())
                    CellWork(function,
                             task_info.partition_row_index[partition] + 2 * j +
                               1,
                             task_info,
                             true);
                }
              else
                {
                  if (odds == evens)
                    {
                      worker[evens] = new (worker[j]->allocate_child())
                        CellWork(function,
                                 task_info.partition_row_index[partition] +
                                   2 * j + 1,
                                 task_info,
                                 false);
                      tbb::task::spawn(*worker[evens]);
                    }
                  else
                    {
                      tbb::empty_task *child =
                        new (worker[j]->allocate_child()) tbb::empty_task();
                      tbb::task::spawn(*child);
                    }
                }
            }

          root->wait_for_all();
          root->destroy(*root);
          if (is_blocked == true)
            tbb::empty_task::spawn(*dummy);
          return nullptr;
        }

        tbb::empty_task *dummy;

      private:
        MFWorkerInterface &function;
        const unsigned int partition;
        const TaskInfo    &task_info;
        const bool         is_blocked;
      };

    } // end of namespace partition



    namespace color
    {
      class CellWork
      {
      public:
        CellWork(MFWorkerInterface &worker_in,
                 const TaskInfo    &task_info_in,
                 const unsigned int partition_in)
          : worker(worker_in)
          , task_info(task_info_in)
          , partition(partition_in)
        {}

        void
        operator()(const tbb::blocked_range<unsigned int> &r) const
        {
          const unsigned int start_index =
            task_info.cell_partition_data[partition] +
            task_info.block_size * r.begin();
          const unsigned int end_index =
            std::min(start_index + task_info.block_size * (r.end() - r.begin()),
                     task_info.cell_partition_data[partition + 1]);
          worker.cell(std::make_pair(start_index, end_index));

          if (task_info.face_partition_data.empty() == false)
            {
              AssertThrow(false, ExcNotImplemented());
            }
        }

      private:
        MFWorkerInterface &worker;
        const TaskInfo    &task_info;
        const unsigned int partition;
      };



      class PartitionWork : public tbb::task
      {
      public:
        PartitionWork(MFWorkerInterface &worker_in,
                      const unsigned int partition_in,
                      const TaskInfo    &task_info_in,
                      const bool         is_blocked_in)
          : dummy(nullptr)
          , worker(worker_in)
          , partition(partition_in)
          , task_info(task_info_in)
          , is_blocked(is_blocked_in)
        {}

        tbb::task *
        execute() override
        {
          const unsigned int n_chunks =
            (task_info.cell_partition_data[partition + 1] -
             task_info.cell_partition_data[partition] + task_info.block_size -
             1) /
            task_info.block_size;
          parallel_for(tbb::blocked_range<unsigned int>(0, n_chunks, 1),
                       CellWork(worker, task_info, partition));
          if (is_blocked == true)
            tbb::empty_task::spawn(*dummy);
          return nullptr;
        }

        tbb::empty_task *dummy;

      private:
        MFWorkerInterface &worker;
        const unsigned int partition;
        const TaskInfo    &task_info;
        const bool         is_blocked;
      };

    } // end of namespace color



    class MPICommunication : public tbb::task
    {
    public:
      MPICommunication(MFWorkerInterface &worker_in, const bool do_compress)
        : worker(worker_in)
        , do_compress(do_compress)
      {}

      tbb::task *
      execute() override
      {
        if (do_compress == false)
          worker.vector_update_ghosts_finish();
        else
          worker.vector_compress_start();
        return nullptr;
      }

    private:
      MFWorkerInterface &worker;
      const bool         do_compress;
    };

#endif // DEAL_II_WITH_TBB



    void
    TaskInfo::loop(MFWorkerInterface &funct) const
    {
      // If we use thread parallelism, we do not currently support to schedule
      // pieces of updates within the loop, so this index will collect all
      // calls in that case and work like a single complete loop over all
      // cells
      if (scheme != none)
        funct.cell_loop_pre_range(numbers::invalid_unsigned_int);
      else
        funct.cell_loop_pre_range(
          partition_row_index[partition_row_index.size() - 2]);

      funct.vector_update_ghosts_start();

#if defined(DEAL_II_WITH_TBB) && !defined(DEAL_II_TBB_WITH_ONEAPI)

      if (scheme != none)
        {
          funct.zero_dst_vector_range(numbers::invalid_unsigned_int);
          if (scheme == partition_partition && evens > 0)
            {
              tbb::empty_task *root =
                new (tbb::task::allocate_root()) tbb::empty_task;
              root->set_ref_count(evens + 1);
              std::vector<partition::PartitionWork *> worker(n_workers);
              std::vector<partition::PartitionWork *> blocked_worker(
                n_blocked_workers);
              MPICommunication *worker_compr =
                new (root->allocate_child()) MPICommunication(funct, true);
              worker_compr->set_ref_count(1);
              for (unsigned int j = 0; j < evens; ++j)
                {
                  if (j > 0)
                    {
                      worker[j] = new (root->allocate_child())
                        partition::PartitionWork(funct, 2 * j, *this, false);
                      worker[j]->set_ref_count(2);
                      blocked_worker[j - 1]->dummy =
                        new (worker[j]->allocate_child()) tbb::empty_task;
                      tbb::task::spawn(*blocked_worker[j - 1]);
                    }
                  else
                    {
                      worker[j] = new (worker_compr->allocate_child())
                        partition::PartitionWork(funct, 2 * j, *this, false);
                      worker[j]->set_ref_count(2);
                      MPICommunication *worker_dist =
                        new (worker[j]->allocate_child())
                          MPICommunication(funct, false);
                      tbb::task::spawn(*worker_dist);
                    }
                  if (j < evens - 1)
                    {
                      blocked_worker[j] = new (worker[j]->allocate_child())
                        partition::PartitionWork(funct, 2 * j + 1, *this, true);
                    }
                  else
                    {
                      if (odds == evens)
                        {
                          worker[evens] = new (worker[j]->allocate_child())
                            partition::PartitionWork(funct,
                                                     2 * j + 1,
                                                     *this,
                                                     false);
                          tbb::task::spawn(*worker[evens]);
                        }
                      else
                        {
                          tbb::empty_task *child =
                            new (worker[j]->allocate_child()) tbb::empty_task();
                          tbb::task::spawn(*child);
                        }
                    }
                }

              root->wait_for_all();
              root->destroy(*root);
            }
          else if (scheme == partition_partition)
            {
              // catch the case of empty partition list: we still need to call
              // the vector communication routines to clean up and initiate
              // things
              funct.vector_update_ghosts_finish();
              funct.vector_compress_start();
            }
          else // end of partition-partition, start of partition-color
            {
              // check whether there is only one partition. if not, build up the
              // tree of partitions
              if (odds > 0)
                {
                  tbb::empty_task *root =
                    new (tbb::task::allocate_root()) tbb::empty_task;
                  root->set_ref_count(evens + 1);
                  const unsigned int n_blocked_workers =
                    odds - (odds + evens + 1) % 2;
                  const unsigned int n_workers =
                    cell_partition_data.size() - 1 - n_blocked_workers;
                  std::vector<color::PartitionWork *> worker(n_workers);
                  std::vector<color::PartitionWork *> blocked_worker(
                    n_blocked_workers);
                  unsigned int      worker_index = 0, slice_index = 0;
                  int               spawn_index_child = -2;
                  MPICommunication *worker_compr =
                    new (root->allocate_child()) MPICommunication(funct, true);
                  worker_compr->set_ref_count(1);
                  for (unsigned int part = 0;
                       part < partition_row_index.size() - 1;
                       part++)
                    {
                      if (part == 0)
                        worker[worker_index] =
                          new (worker_compr->allocate_child())
                            color::PartitionWork(funct,
                                                 slice_index,
                                                 *this,
                                                 false);
                      else
                        worker[worker_index] = new (root->allocate_child())
                          color::PartitionWork(funct,
                                               slice_index,
                                               *this,
                                               false);
                      ++slice_index;
                      for (; slice_index < partition_row_index[part + 1];
                           slice_index++)
                        {
                          worker[worker_index]->set_ref_count(1);
                          ++worker_index;
                          worker[worker_index] =
                            new (worker[worker_index - 1]->allocate_child())
                              color::PartitionWork(funct,
                                                   slice_index,
                                                   *this,
                                                   false);
                        }
                      worker[worker_index]->set_ref_count(2);
                      if (part > 0)
                        {
                          blocked_worker[(part - 1) / 2]->dummy =
                            new (worker[worker_index]->allocate_child())
                              tbb::empty_task;
                          ++worker_index;
                          if (spawn_index_child == -1)
                            tbb::task::spawn(*blocked_worker[(part - 1) / 2]);
                          else
                            {
                              Assert(spawn_index_child >= 0,
                                     ExcInternalError());
                              tbb::task::spawn(*worker[spawn_index_child]);
                            }
                          spawn_index_child = -2;
                        }
                      else
                        {
                          MPICommunication *worker_dist =
                            new (worker[worker_index]->allocate_child())
                              MPICommunication(funct, false);
                          tbb::task::spawn(*worker_dist);
                          ++worker_index;
                        }
                      part += 1;
                      if (part < partition_row_index.size() - 1)
                        {
                          if (part < partition_row_index.size() - 2)
                            {
                              blocked_worker[part / 2] =
                                new (worker[worker_index - 1]->allocate_child())
                                  color::PartitionWork(funct,
                                                       slice_index,
                                                       *this,
                                                       true);
                              ++slice_index;
                              if (slice_index < partition_row_index[part + 1])
                                {
                                  blocked_worker[part / 2]->set_ref_count(1);
                                  worker[worker_index] = new (
                                    blocked_worker[part / 2]->allocate_child())
                                    color::PartitionWork(funct,
                                                         slice_index,
                                                         *this,
                                                         false);
                                  ++slice_index;
                                }
                              else
                                {
                                  spawn_index_child = -1;
                                  continue;
                                }
                            }
                          for (; slice_index < partition_row_index[part + 1];
                               slice_index++)
                            {
                              if (slice_index > partition_row_index[part])
                                {
                                  worker[worker_index]->set_ref_count(1);
                                  ++worker_index;
                                }
                              worker[worker_index] =
                                new (worker[worker_index - 1]->allocate_child())
                                  color::PartitionWork(funct,
                                                       slice_index,
                                                       *this,
                                                       false);
                            }
                          spawn_index_child = worker_index;
                          ++worker_index;
                        }
                      else
                        {
                          tbb::empty_task *final =
                            new (worker[worker_index - 1]->allocate_child())
                              tbb::empty_task;
                          tbb::task::spawn(*final);
                          spawn_index_child = worker_index - 1;
                        }
                    }
                  if (evens == odds)
                    {
                      Assert(spawn_index_child >= 0, ExcInternalError());
                      tbb::task::spawn(*worker[spawn_index_child]);
                    }
                  root->wait_for_all();
                  root->destroy(*root);
                }
              // case when we only have one partition: this is the usual
              // coloring scheme, and we just schedule a parallel for loop for
              // each color
              else
                {
                  Assert(evens <= 1, ExcInternalError());
                  funct.vector_update_ghosts_finish();

                  for (unsigned int color = 0; color < partition_row_index[1];
                       ++color)
                    {
                      tbb::empty_task *root =
                        new (tbb::task::allocate_root()) tbb::empty_task;
                      root->set_ref_count(2);
                      color::PartitionWork *worker =
                        new (root->allocate_child())
                          color::PartitionWork(funct, color, *this, false);
                      tbb::empty_task::spawn(*worker);
                      root->wait_for_all();
                      root->destroy(*root);
                    }

                  funct.vector_compress_start();
                }
            }
        }
      else
#endif
        // serial loop, go through up to three times and do the MPI transfer at
        // the beginning/end of the second part
        {
          for (unsigned int part = 0; part < partition_row_index.size() - 2;
               ++part)
            {
              if (part == 1)
                funct.vector_update_ghosts_finish();

              for (unsigned int i = partition_row_index[part];
                   i < partition_row_index[part + 1];
                   ++i)
                {
                  funct.cell_loop_pre_range(i);
                  funct.zero_dst_vector_range(i);
                  AssertIndexRange(i + 1, cell_partition_data.size());
                  if (cell_partition_data[i + 1] > cell_partition_data[i])
                    {
                      funct.cell(i);
                    }

                  if (face_partition_data.empty() == false)
                    {
                      if (face_partition_data[i + 1] > face_partition_data[i])
                        funct.face(i);
                      if (boundary_partition_data[i + 1] >
                          boundary_partition_data[i])
                        funct.boundary(i);
                    }
                  funct.cell_loop_post_range(i);
                }

              if (part == 1)
                funct.vector_compress_start();
            }
        }
      funct.vector_compress_finish();

      if (scheme != none)
        funct.cell_loop_post_range(numbers::invalid_unsigned_int);
      else
        funct.cell_loop_post_range(
          partition_row_index[partition_row_index.size() - 2]);
    }



    TaskInfo::TaskInfo()
    {
      clear();
    }



    void
    TaskInfo::clear()
    {
      n_active_cells       = 0;
      n_ghost_cells        = 0;
      vectorization_length = 1;
      block_size           = 0;
      n_blocks             = 0;
      scheme               = none;
      partition_row_index.clear();
      partition_row_index.resize(2);
      cell_partition_data.clear();
      face_partition_data.clear();
      boundary_partition_data.clear();
      evens             = 0;
      odds              = 0;
      n_blocked_workers = 0;
      n_workers         = 0;
      partition_evens.clear();
      partition_odds.clear();
      partition_n_blocked_workers.clear();
      partition_n_workers.clear();
      communicator = MPI_COMM_SELF;
      my_pid       = 0;
      n_procs      = 1;
    }



    template <typename StreamType>
    void
    TaskInfo::print_memory_statistics(StreamType       &out,
                                      const std::size_t data_length) const
    {
      Utilities::MPI::MinMaxAvg memory_c =
        Utilities::MPI::min_max_avg(1e-6 * data_length, communicator);
      if (n_procs < 2)
        out << memory_c.min;
      else
        out << memory_c.min << "/" << memory_c.avg << "/" << memory_c.max;
      out << " MB" << std::endl;
    }



    std::size_t
    TaskInfo::memory_consumption() const
    {
      return (
        sizeof(*this) +
        MemoryConsumption::memory_consumption(partition_row_index) +
        MemoryConsumption::memory_consumption(cell_partition_data) +
        MemoryConsumption::memory_consumption(face_partition_data) +
        MemoryConsumption::memory_consumption(boundary_partition_data) +
        MemoryConsumption::memory_consumption(partition_evens) +
        MemoryConsumption::memory_consumption(partition_odds) +
        MemoryConsumption::memory_consumption(partition_n_blocked_workers) +
        MemoryConsumption::memory_consumption(partition_n_workers));
    }



    void
    TaskInfo::make_boundary_cells_divisible(
      std::vector<unsigned int> &boundary_cells)
    {
      // try to make the number of boundary cells divisible by the number of
      // vectors in vectorization
      unsigned int fillup_needed =
        (vectorization_length - boundary_cells.size() % vectorization_length) %
        vectorization_length;
      if (fillup_needed > 0 && boundary_cells.size() < n_active_cells)
        {
          // fill additional cells into the list of boundary cells to get a
          // balanced number. Go through the indices successively until we
          // found enough indices
          std::vector<unsigned int> new_boundary_cells;
          new_boundary_cells.reserve(boundary_cells.size());

          unsigned int next_free_slot = 0, bound_index = 0;
          while (fillup_needed > 0 && bound_index < boundary_cells.size())
            {
              if (next_free_slot < boundary_cells[bound_index])
                {
                  // check if there are enough cells to fill with in the
                  // current slot
                  if (next_free_slot + fillup_needed <=
                      boundary_cells[bound_index])
                    {
                      for (unsigned int j =
                             boundary_cells[bound_index] - fillup_needed;
                           j < boundary_cells[bound_index];
                           ++j)
                        new_boundary_cells.push_back(j);
                      fillup_needed = 0;
                    }
                  // ok, not enough indices, so just take them all up to the
                  // next boundary cell
                  else
                    {
                      for (unsigned int j = next_free_slot;
                           j < boundary_cells[bound_index];
                           ++j)
                        new_boundary_cells.push_back(j);
                      fillup_needed -=
                        boundary_cells[bound_index] - next_free_slot;
                    }
                }
              new_boundary_cells.push_back(boundary_cells[bound_index]);
              next_free_slot = boundary_cells[bound_index] + 1;
              ++bound_index;
            }
          while (fillup_needed > 0 &&
                 (new_boundary_cells.empty() ||
                  new_boundary_cells.back() < n_active_cells - 1))
            new_boundary_cells.push_back(new_boundary_cells.back() + 1);
          while (bound_index < boundary_cells.size())
            new_boundary_cells.push_back(boundary_cells[bound_index++]);

          boundary_cells.swap(new_boundary_cells);
        }

      // set the number of cells
      std::sort(boundary_cells.begin(), boundary_cells.end());

      // check that number of boundary cells is divisible by
      // vectorization_length or that it contains all cells
      Assert(boundary_cells.size() % vectorization_length == 0 ||
               boundary_cells.size() == n_active_cells,
             ExcInternalError());
    }



    void
    TaskInfo::create_blocks_serial(
      const std::vector<unsigned int> &cells_with_comm,
      const unsigned int               dofs_per_cell,
      const bool                       categories_are_hp,
      const std::vector<unsigned int> &cell_vectorization_categories,
      const bool                       cell_vectorization_categories_strict,
      const std::vector<unsigned int> &parent_relation,
      std::vector<unsigned int>       &renumbering,
      std::vector<unsigned char>      &incompletely_filled_vectorization)
    {
      Assert(dofs_per_cell > 0, ExcInternalError());
      // This function is decomposed into several steps to determine a good
      // ordering that satisfies the following constraints:
      // a. Only cells belonging to the same category (or next higher if the
      // cell_vectorization_categories_strict is false) can be grouped into
      // the same SIMD batch
      // b. hp-adaptive computations must form contiguous ranges for the same
      // degree (category) in cell_partition_data
      // c. We want to group the cells with the same parent in the same SIMD
      // lane if possible
      // d. The cell order should be similar to the initial one
      // e. Form sets without MPI communication and those with to overlap
      // communication with computation
      //
      // These constraints are satisfied by first grouping by the categories
      // and, within the groups, to distinguish between cells with a parent
      // and those without. All of this is set up with batches of cells (with
      // padding if the size does not match). Then we define a vector of
      // arrays where we define sorting criteria for the cell batches to
      // satisfy the items b and d together, split by different parts to
      // satisfy item e.

      // Give the compiler a chance to detect that vectorization_length is a
      // power of two, which allows it to replace integer divisions by shifts
      const unsigned int n_lanes = indicate_power_of_two(vectorization_length);

      // Step 1: find tight map of categories for not taking exceeding amounts
      // of memory below. Sort the new categories by the numbers in the
      // old one to ensure we respect the given rules
      unsigned int              n_categories = 1;
      std::vector<unsigned int> tight_category_map(n_active_cells, 0);
      if (cell_vectorization_categories.empty() == false)
        {
          AssertDimension(cell_vectorization_categories.size(),
                          n_active_cells + n_ghost_cells);

          std::set<unsigned int> used_categories;
          for (unsigned int i = 0; i < n_active_cells; ++i)
            used_categories.insert(cell_vectorization_categories[i]);
          std::vector<unsigned int> used_categories_vector(
            used_categories.size());
          n_categories = 0;
          for (const auto &it : used_categories)
            used_categories_vector[n_categories++] = it;
          for (unsigned int i = 0; i < n_active_cells; ++i)
            {
              const unsigned int index =
                std::lower_bound(used_categories_vector.begin(),
                                 used_categories_vector.end(),
                                 cell_vectorization_categories[i]) -
                used_categories_vector.begin();
              AssertIndexRange(index, used_categories_vector.size());
              tight_category_map[i] = index;
            }
        }

      // Step 2: Sort the cells by the category. If we want to fill up the
      // ranges in vectorization, promote some of the cells to a higher
      // category
      std::vector<std::vector<unsigned int>> renumbering_category(n_categories);
      for (unsigned int i = 0; i < n_active_cells; ++i)
        renumbering_category[tight_category_map[i]].push_back(i);

      if (cell_vectorization_categories_strict == false && n_categories > 1)
        for (unsigned int j = n_categories - 1; j > 0; --j)
          {
            unsigned int lower_index = j - 1;
            while ((renumbering_category[j].size() % n_lanes) != 0u)
              {
                while (((renumbering_category[j].size() % n_lanes) != 0u) &&
                       !renumbering_category[lower_index].empty())
                  {
                    renumbering_category[j].push_back(
                      renumbering_category[lower_index].back());
                    renumbering_category[lower_index].pop_back();
                  }
                if (lower_index == 0)
                  break;
                else
                  --lower_index;
              }
          }

      // Step 3: Use the parent relation to find a good grouping of cells. To
      // do this, we first put cells of each category defined above into two
      // bins, those which we know can be grouped together by the given parent
      // relation and those which cannot
      std::vector<unsigned int> temporary_numbering;
      temporary_numbering.reserve(n_active_cells +
                                  (n_lanes - 1) * n_categories);
      const unsigned int n_cells_per_parent =
        std::count(parent_relation.begin(), parent_relation.end(), 0);
      std::vector<unsigned int> category_size;
      for (unsigned int j = 0; j < n_categories; ++j)
        {
          std::vector<std::pair<unsigned int, unsigned int>> grouped_cells;
          std::vector<unsigned int>                          other_cells;
          for (const unsigned int cell : renumbering_category[j])
            if (parent_relation.empty() ||
                parent_relation[cell] == numbers::invalid_unsigned_int)
              other_cells.push_back(cell);
            else
              grouped_cells.emplace_back(parent_relation[cell], cell);

          // Compute the number of cells per group
          std::sort(grouped_cells.begin(), grouped_cells.end());
          std::vector<unsigned int> n_cells_per_group;
          unsigned int              length = 0;
          for (unsigned int i = 0; i < grouped_cells.size(); ++i, ++length)
            if (i > 0 && grouped_cells[i].first != grouped_cells[i - 1].first)
              {
                n_cells_per_group.push_back(length);
                length = 0;
              }
          if (length > 0)
            n_cells_per_group.push_back(length);

          // Move groups that do not have the complete size (due to
          // categories) to the 'other_cells'. The cells with correct group
          // size are immediately appended to the temporary cell numbering
          auto group_it = grouped_cells.begin();
          for (const unsigned int length : n_cells_per_group)
            if (length < n_cells_per_parent)
              for (unsigned int j = 0; j < length; ++j)
                other_cells.push_back((group_it++)->second);
            else
              {
                // we should not have more cells in a group than in the first
                // check we did above
                AssertDimension(length, n_cells_per_parent);
                for (unsigned int j = 0; j < length; ++j)
                  temporary_numbering.push_back((group_it++)->second);
              }

          // Sort the remaining cells and append them as well
          std::sort(other_cells.begin(), other_cells.end());
          temporary_numbering.insert(temporary_numbering.end(),
                                     other_cells.begin(),
                                     other_cells.end());

          while (temporary_numbering.size() % n_lanes != 0)
            temporary_numbering.push_back(numbers::invalid_unsigned_int);

          category_size.push_back(temporary_numbering.size());
        }

      // Step 4: Identify the batches with cells marked as "comm"
      std::vector<bool> batch_with_comm(temporary_numbering.size() / n_lanes,
                                        false);
      std::vector<unsigned int> temporary_numbering_inverse(n_active_cells);
      for (unsigned int i = 0; i < temporary_numbering.size(); ++i)
        if (temporary_numbering[i] != numbers::invalid_unsigned_int)
          temporary_numbering_inverse[temporary_numbering[i]] = i;
      for (const unsigned int cell : cells_with_comm)
        batch_with_comm[temporary_numbering_inverse[cell] / n_lanes] = true;

      // Step 5: Sort the batches of cells by their last cell index to get
      // good locality, assuming that the initial cell order is of good
      // locality. In case we have hp-calculations with categories, we need to
      // sort also by the category.
      std::vector<std::array<unsigned int, 3>> batch_order;
      std::vector<std::array<unsigned int, 3>> batch_order_comm;
      for (unsigned int i = 0; i < temporary_numbering.size(); i += n_lanes)
        {
          unsigned int max_index = 0;
          for (unsigned int j = 0; j < n_lanes; ++j)
            if (temporary_numbering[i + j] < numbers::invalid_unsigned_int)
              max_index = std::max(temporary_numbering[i + j], max_index);
          const unsigned int category_hp =
            categories_are_hp ?
              std::upper_bound(category_size.begin(), category_size.end(), i) -
                category_size.begin() :
              0;
          const std::array<unsigned int, 3> next{{category_hp, max_index, i}};
          if (batch_with_comm[i / n_lanes])
            batch_order_comm.emplace_back(next);
          else
            batch_order.emplace_back(next);
        }

      std::sort(batch_order.begin(), batch_order.end());
      std::sort(batch_order_comm.begin(), batch_order_comm.end());

      // Step 6: Put the cells with communication in the middle of the cell
      // range. For the MPI case, we need three groups to enable overlap for
      // communication and computation (part before comm, part with comm, part
      // after comm), whereas we need one for the other case. And in each
      // case, we allow for a slot of "ghosted" cells.
      std::vector<unsigned int> blocks;
      if (n_procs == 1)
        {
          if (batch_order.empty())
            std::swap(batch_order_comm, batch_order);
          Assert(batch_order_comm.empty(), ExcInternalError());
          partition_row_index.resize(3);
          blocks = {0, static_cast<unsigned int>(batch_order.size())};
        }
      else
        {
          partition_row_index.resize(5);
          const unsigned int comm_begin = batch_order.size() / 2;
          batch_order.insert(batch_order.begin() + comm_begin,
                             batch_order_comm.begin(),
                             batch_order_comm.end());
          const unsigned int comm_end = comm_begin + batch_order_comm.size();
          const unsigned int end      = batch_order.size();
          blocks                      = {0, comm_begin, comm_end, end};
        }

      // Step 7: sort ghost cells according to the category
      std::vector<std::array<unsigned int, 2>> tight_category_map_ghost;

      if (cell_vectorization_categories.empty() == false)
        {
          tight_category_map_ghost.reserve(n_ghost_cells);

          std::set<unsigned int> used_categories;
          for (unsigned int i = 0; i < n_ghost_cells; ++i)
            used_categories.insert(
              cell_vectorization_categories[i + n_active_cells]);

          std::vector<unsigned int> used_categories_vector(
            used_categories.size());
          n_categories = 0;
          for (const auto &it : used_categories)
            used_categories_vector[n_categories++] = it;

          std::vector<unsigned int> counters(n_categories, 0);

          for (unsigned int i = 0; i < n_ghost_cells; ++i)
            {
              const unsigned int index =
                std::lower_bound(
                  used_categories_vector.begin(),
                  used_categories_vector.end(),
                  cell_vectorization_categories[i + n_active_cells]) -
                used_categories_vector.begin();
              AssertIndexRange(index, used_categories_vector.size());
              tight_category_map_ghost.emplace_back(
                std::array<unsigned int, 2>{{index, i}});

              // account for padding in the hp and strict case
              if (categories_are_hp || cell_vectorization_categories_strict)
                counters[index]++;
            }

          // insert padding
          for (unsigned int i = 0; i < counters.size(); ++i)
            if (counters[i] % n_lanes != 0)
              for (unsigned int j = counters[i] % n_lanes; j < n_lanes; ++j)
                tight_category_map_ghost.emplace_back(
                  std::array<unsigned int, 2>{
                    {i, numbers::invalid_unsigned_int}});

          std::sort(tight_category_map_ghost.begin(),
                    tight_category_map_ghost.end());
        }

      // Step 8: Fill in the data by batches for the locally owned cells.
      const unsigned int n_cell_batches = batch_order.size();
      const unsigned int n_ghost_batches =
        ((tight_category_map_ghost.empty() ? n_ghost_cells :
                                             tight_category_map_ghost.size()) +
         n_lanes - 1) /
        n_lanes;
      incompletely_filled_vectorization.resize(n_cell_batches +
                                               n_ghost_batches);

      cell_partition_data.clear();
      cell_partition_data.resize(1, 0);

      renumbering.clear();
      renumbering.resize(n_active_cells + n_ghost_cells,
                         numbers::invalid_unsigned_int);

      unsigned int counter = 0;
      for (unsigned int block = 0; block < blocks.size() - 1; ++block)
        {
          const unsigned int grain_size =
            std::max((2048U / dofs_per_cell) / 8 * 4, 2U);
          for (unsigned int k = blocks[block]; k < blocks[block + 1];
               k += grain_size)
            cell_partition_data.push_back(
              std::min(k + grain_size, blocks[block + 1]));
          partition_row_index[block + 1] = cell_partition_data.size() - 1;

          // Set the numbering according to the reordered temporary one
          for (unsigned int k = blocks[block]; k < blocks[block + 1]; ++k)
            {
              const unsigned int pos = batch_order[k][2];
              unsigned int       j   = 0;
              for (; j < n_lanes && temporary_numbering[pos + j] !=
                                      numbers::invalid_unsigned_int;
                   ++j)
                renumbering[counter++] = temporary_numbering[pos + j];
              if (j < n_lanes)
                incompletely_filled_vectorization[k] = j;
            }
        }
      AssertDimension(counter, n_active_cells);

      // Step 9: Treat the ghost cells
      if (tight_category_map_ghost.empty())
        {
          for (unsigned int cell = 0; cell < n_ghost_cells; ++cell)
            renumbering[n_active_cells + cell] = n_active_cells + cell;

          if ((n_ghost_cells % n_lanes) != 0u)
            incompletely_filled_vectorization.back() = n_ghost_cells % n_lanes;
        }
      else
        {
          for (unsigned int k = 0, ptr = 0; k < n_ghost_batches;
               ++k, ptr += n_lanes)
            {
              unsigned int j = 0;

              for (;
                   j < n_lanes && (ptr + j < tight_category_map_ghost.size()) &&
                   (tight_category_map_ghost[ptr + j][1] !=
                    numbers::invalid_unsigned_int);
                   ++j)
                renumbering[counter++] =
                  n_active_cells + tight_category_map_ghost[ptr + j][1];

              if (j < n_lanes)
                incompletely_filled_vectorization[n_cell_batches + k] = j;
            }

          AssertDimension(counter, n_active_cells + n_ghost_cells);
        }

      cell_partition_data.push_back(n_cell_batches + n_ghost_batches);
      partition_row_index.back() = cell_partition_data.size() - 1;

      if constexpr (running_in_debug_mode())
        {
          std::vector<unsigned int> renumber_cpy(renumbering);
          std::sort(renumber_cpy.begin(), renumber_cpy.end());
          for (unsigned int i = 0; i < renumber_cpy.size(); ++i)
            AssertDimension(i, renumber_cpy[i]);
        }
    }



    void
    TaskInfo::initial_setup_blocks_tasks(
      const std::vector<unsigned int> &boundary_cells,
      std::vector<unsigned int>       &renumbering,
      std::vector<unsigned char>      &incompletely_filled_vectorization)
    {
      const unsigned int n_cell_batches =
        (n_active_cells + vectorization_length - 1) / vectorization_length;
      const unsigned int n_ghost_slots =
        (n_ghost_cells + vectorization_length - 1) / vectorization_length;
      incompletely_filled_vectorization.resize(n_cell_batches + n_ghost_slots);
      if (n_cell_batches * vectorization_length > n_active_cells)
        incompletely_filled_vectorization[n_cell_batches - 1] =
          vectorization_length -
          (n_cell_batches * vectorization_length - n_active_cells);
      if (n_ghost_slots * vectorization_length > n_ghost_cells)
        incompletely_filled_vectorization[n_cell_batches + n_ghost_slots - 1] =
          vectorization_length -
          (n_ghost_slots * vectorization_length - n_ghost_cells);

      std::vector<unsigned int> reverse_numbering(
        n_active_cells, numbers::invalid_unsigned_int);
      for (unsigned int j = 0; j < boundary_cells.size(); ++j)
        reverse_numbering[boundary_cells[j]] = j;
      unsigned int counter = boundary_cells.size();
      for (unsigned int j = 0; j < n_active_cells; ++j)
        if (reverse_numbering[j] == numbers::invalid_unsigned_int)
          reverse_numbering[j] = counter++;

      AssertDimension(counter, n_active_cells);
      renumbering = Utilities::invert_permutation(reverse_numbering);

      for (unsigned int j = n_active_cells; j < n_active_cells + n_ghost_cells;
           ++j)
        renumbering.push_back(j);

      // TODO: might be able to simplify this code by not relying on the cell
      // partition data while computing the thread graph
      cell_partition_data.clear();
      cell_partition_data.push_back(0);
      if (n_procs > 1)
        {
          const unsigned int n_macro_boundary_cells =
            (boundary_cells.size() + vectorization_length - 1) /
            vectorization_length;
          cell_partition_data.push_back(
            (n_cell_batches - n_macro_boundary_cells) / 2);
          cell_partition_data.push_back(cell_partition_data[1] +
                                        n_macro_boundary_cells);
        }
      else
        AssertDimension(boundary_cells.size(), 0);
      cell_partition_data.push_back(n_cell_batches);
      cell_partition_data.push_back(cell_partition_data.back() + n_ghost_slots);
      partition_row_index.resize(n_procs > 1 ? 4 : 2);
      partition_row_index[0] = 0;
      partition_row_index[1] = 1;
      if (n_procs > 1)
        {
          partition_row_index[2] = 2;
          partition_row_index[3] = 3;
        }
    }



    void
    TaskInfo::guess_block_size(const unsigned int dofs_per_cell)
    {
      // user did not say a positive number, so we have to guess
      if (block_size == 0)
        {
          // we would like to have enough work to do, so as first guess, try
          // to get 16 times as many chunks as we have threads on the system.
          block_size = n_active_cells / (MultithreadInfo::n_threads() * 16 *
                                         vectorization_length);

          // if there are too few degrees of freedom per cell, need to
          // increase the block size
          const unsigned int minimum_parallel_grain_size = 200;
          if (dofs_per_cell * block_size < minimum_parallel_grain_size)
            block_size = (minimum_parallel_grain_size / dofs_per_cell + 1);
          if (dofs_per_cell * block_size > 10000)
            block_size /= 4;

          block_size =
            1 << static_cast<unsigned int>(std::log2(block_size + 1));
        }
      if (block_size > n_active_cells)
        block_size = std::max(1U, n_active_cells);
    }



    void
    TaskInfo::make_thread_graph_partition_color(
      DynamicSparsityPattern     &connectivity_large,
      std::vector<unsigned int>  &renumbering,
      std::vector<unsigned char> &irregular_cells,
      const bool)
    {
      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      if (n_cell_batches == 0)
        return;

      Assert(vectorization_length > 0, ExcInternalError());

      unsigned int partition = 0, counter = 0;

      // Create connectivity graph for blocks based on connectivity graph for
      // cells.
      DynamicSparsityPattern connectivity(n_blocks, n_blocks);
      make_connectivity_cells_to_blocks(irregular_cells,
                                        connectivity_large,
                                        connectivity);

      // Create cell-block  partitioning.

      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to -1 to mark them as not yet assigned
      // a partition.
      std::vector<unsigned int> cell_partition(n_blocks,
                                               numbers::invalid_unsigned_int);

      // In element j of this variable, one puts the old number of the block
      // that should be the jth block in the new numeration.
      std::vector<unsigned int> partition_list(n_blocks, 0);
      std::vector<unsigned int> partition_color_list(n_blocks, 0);

      // This vector points to the start of each partition.
      std::vector<unsigned int> partition_size(2, 0);

      // blocking_connectivity = true;

      // The cluster_size in make_partitioning defines that the no. of cells
      // in each partition should be a multiple of cluster_size.
      unsigned int cluster_size = 1;

      // Make the partitioning of the first layer of the blocks of cells.
      make_partitioning(connectivity,
                        cluster_size,
                        cell_partition,
                        partition_list,
                        partition_size,
                        partition);

      // Color the cells within each partition
      make_coloring_within_partitions_pre_blocked(connectivity,
                                                  partition,
                                                  cell_partition,
                                                  partition_list,
                                                  partition_size,
                                                  partition_color_list);

      partition_list = renumbering;

      if constexpr (running_in_debug_mode())
        {
          // in debug mode, check that the partition color list is one-to-one
          {
            std::vector<unsigned int> sorted_pc_list(partition_color_list);
            std::sort(sorted_pc_list.begin(), sorted_pc_list.end());
            for (unsigned int i = 0; i < sorted_pc_list.size(); ++i)
              Assert(sorted_pc_list[i] == i, ExcInternalError());
          }
        }

      // set the start list for each block and compute the renumbering of
      // cells
      std::vector<unsigned int>  block_start(n_cell_batches + 1);
      std::vector<unsigned char> irregular(n_cell_batches);

      unsigned int mcell_start = 0;
      block_start[0]           = 0;
      for (unsigned int block = 0; block < n_blocks; ++block)
        {
          block_start[block + 1] = block_start[block];
          for (unsigned int mcell = mcell_start;
               mcell < std::min(mcell_start + block_size, n_cell_batches);
               ++mcell)
            {
              unsigned int n_comp = (irregular_cells[mcell] > 0) ?
                                      irregular_cells[mcell] :
                                      vectorization_length;
              block_start[block + 1] += n_comp;
              ++counter;
            }
          mcell_start += block_size;
        }
      counter                    = 0;
      unsigned int counter_macro = 0;
      unsigned int block_size_last =
        n_cell_batches - block_size * (n_blocks - 1);
      if (block_size_last == 0)
        block_size_last = block_size;

      unsigned int tick = 0;
      for (unsigned int block = 0; block < n_blocks; ++block)
        {
          unsigned int present_block = partition_color_list[block];
          for (unsigned int cell = block_start[present_block];
               cell < block_start[present_block + 1];
               ++cell)
            renumbering[counter++] = partition_list[cell];
          unsigned int this_block_size =
            (present_block == n_blocks - 1) ? block_size_last : block_size;

          // Also re-compute the content of cell_partition_data to
          // contain the numbers of cells, not blocks
          if (cell_partition_data[tick] == block)
            cell_partition_data[tick++] = counter_macro;

          for (unsigned int j = 0; j < this_block_size; ++j)
            irregular[counter_macro++] =
              irregular_cells[present_block * block_size + j];
        }
      AssertDimension(tick + 1, cell_partition_data.size());
      cell_partition_data.back() = counter_macro;

      irregular_cells.swap(irregular);
      AssertDimension(counter, n_active_cells);
      AssertDimension(counter_macro, n_cell_batches);

      // check that the renumbering is one-to-one
      if constexpr (running_in_debug_mode())
        {
          {
            std::vector<unsigned int> sorted_renumbering(renumbering);
            std::sort(sorted_renumbering.begin(), sorted_renumbering.end());
            for (unsigned int i = 0; i < sorted_renumbering.size(); ++i)
              Assert(sorted_renumbering[i] == i, ExcInternalError());
          }
        }


      update_task_info(
        partition); // Actually sets too much for partition color case

      AssertDimension(cell_partition_data.back(), n_cell_batches);
    }



    void
    TaskInfo::make_thread_graph(
      const std::vector<unsigned int> &cell_active_fe_index,
      DynamicSparsityPattern          &connectivity,
      std::vector<unsigned int>       &renumbering,
      std::vector<unsigned char>      &irregular_cells,
      const bool                       hp_bool)
    {
      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      if (n_cell_batches == 0)
        return;

      Assert(vectorization_length > 0, ExcInternalError());

      // if we want to block before partitioning, create connectivity graph
      // for blocks based on connectivity graph for cells.
      DynamicSparsityPattern connectivity_blocks(n_blocks, n_blocks);
      make_connectivity_cells_to_blocks(irregular_cells,
                                        connectivity,
                                        connectivity_blocks);

      unsigned int n_blocks = 0;
      if (scheme == partition_color ||
          scheme == color) // blocking_connectivity == true
        n_blocks = this->n_blocks;
      else
        n_blocks = n_active_cells;

      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to -1 to mark them as not yet assigned
      // a partition.
      std::vector<unsigned int> cell_partition(n_blocks,
                                               numbers::invalid_unsigned_int);

      // In element j of this variable, one puts the old number (but after
      // renumbering according to the input renumbering) of the block that
      // should be the jth block in the new numeration.
      std::vector<unsigned int> partition_list(n_blocks, 0);
      std::vector<unsigned int> partition_2layers_list(n_blocks, 0);

      // This vector points to the start of each partition.
      std::vector<unsigned int> partition_size(2, 0);

      unsigned int partition = 0;

      // Within the partitions we want to be able to block for the case that
      // we do not block already in the connectivity. The cluster_size in
      // make_partitioning defines that the no. of cells in each partition
      // should be a multiple of cluster_size.
      unsigned int cluster_size = 1;
      if (scheme == partition_partition)
        cluster_size = block_size * vectorization_length;

      // Make the partitioning of the first layer of the blocks of cells.
      if (scheme == partition_color || scheme == color)
        make_partitioning(connectivity_blocks,
                          cluster_size,
                          cell_partition,
                          partition_list,
                          partition_size,
                          partition);
      else
        make_partitioning(connectivity,
                          cluster_size,
                          cell_partition,
                          partition_list,
                          partition_size,
                          partition);

      // Partition or color second layer
      if (scheme == partition_partition)

        {
          // Partition within partitions.
          make_partitioning_within_partitions_post_blocked(
            connectivity,
            cell_active_fe_index,
            partition,
            cluster_size,
            hp_bool,
            cell_partition,
            partition_list,
            partition_size,
            partition_2layers_list,
            irregular_cells);
        }
      else if (scheme == partition_color || scheme == color)
        {
          make_coloring_within_partitions_pre_blocked(connectivity_blocks,
                                                      partition,
                                                      cell_partition,
                                                      partition_list,
                                                      partition_size,
                                                      partition_2layers_list);
        }

      // in debug mode, check that the partition_2layers_list is one-to-one
      if constexpr (running_in_debug_mode())
        {
          {
            std::vector<unsigned int> sorted_pc_list(partition_2layers_list);
            std::sort(sorted_pc_list.begin(), sorted_pc_list.end());
            for (unsigned int i = 0; i < sorted_pc_list.size(); ++i)
              Assert(sorted_pc_list[i] == i, ExcInternalError());
          }
        }

      // Set the new renumbering
      std::vector<unsigned int> renumbering_in(n_active_cells, 0);
      renumbering_in.swap(renumbering);
      if (scheme == partition_partition) // blocking_connectivity == false
        {
          // This is the simple case. The renumbering is just a combination of
          // the renumbering that we were given as an input and the
          // renumbering of partition/coloring given in partition_2layers_list
          for (unsigned int j = 0; j < renumbering.size(); ++j)
            renumbering[j] = renumbering_in[partition_2layers_list[j]];
          // Account for the ghost cells, finally.
          for (unsigned int i = 0; i < n_ghost_cells; ++i)
            renumbering.push_back(i + n_active_cells);
        }
      else
        {
          // set the start list for each block and compute the renumbering of
          // cells
          std::vector<unsigned int>  block_start(n_cell_batches + 1);
          std::vector<unsigned char> irregular(n_cell_batches);

          unsigned int counter     = 0;
          unsigned int mcell_start = 0;
          block_start[0]           = 0;
          for (unsigned int block = 0; block < n_blocks; ++block)
            {
              block_start[block + 1] = block_start[block];
              for (unsigned int mcell = mcell_start;
                   mcell < std::min(mcell_start + block_size, n_cell_batches);
                   ++mcell)
                {
                  unsigned int n_comp = (irregular_cells[mcell] > 0) ?
                                          irregular_cells[mcell] :
                                          vectorization_length;
                  block_start[block + 1] += n_comp;
                  ++counter;
                }
              mcell_start += block_size;
            }
          counter                    = 0;
          unsigned int counter_macro = 0;
          unsigned int block_size_last =
            n_cell_batches - block_size * (n_blocks - 1);
          if (block_size_last == 0)
            block_size_last = block_size;

          unsigned int tick = 0;
          for (unsigned int block = 0; block < n_blocks; ++block)
            {
              unsigned int present_block = partition_2layers_list[block];
              for (unsigned int cell = block_start[present_block];
                   cell < block_start[present_block + 1];
                   ++cell)
                renumbering[counter++] = renumbering_in[cell];
              unsigned int this_block_size =
                (present_block == n_blocks - 1) ? block_size_last : block_size;

              // Also re-compute the content of cell_partition_data to
              // contain the numbers of cells, not blocks
              if (cell_partition_data[tick] == block)
                cell_partition_data[tick++] = counter_macro;

              for (unsigned int j = 0; j < this_block_size; ++j)
                irregular[counter_macro++] =
                  irregular_cells[present_block * block_size + j];
            }
          AssertDimension(tick + 1, cell_partition_data.size());
          cell_partition_data.back() = counter_macro;

          irregular_cells.swap(irregular);
          AssertDimension(counter, n_active_cells);
          AssertDimension(counter_macro, n_cell_batches);
          // check that the renumbering is one-to-one
          if constexpr (running_in_debug_mode())
            {
              {
                std::vector<unsigned int> sorted_renumbering(renumbering);
                std::sort(sorted_renumbering.begin(), sorted_renumbering.end());
                for (unsigned int i = 0; i < sorted_renumbering.size(); ++i)
                  Assert(sorted_renumbering[i] == i, ExcInternalError());
              }
            }
        }

      // Update the task_info with the more information for the thread graph.
      update_task_info(partition);
    }



    void
    TaskInfo::make_thread_graph_partition_partition(
      const std::vector<unsigned int> &cell_active_fe_index,
      DynamicSparsityPattern          &connectivity,
      std::vector<unsigned int>       &renumbering,
      std::vector<unsigned char>      &irregular_cells,
      const bool                       hp_bool)
    {
      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      if (n_cell_batches == 0)
        return;

      const unsigned int cluster_size = block_size * vectorization_length;

      // Create cell-block  partitioning.

      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to n_cell_batches to mark them as not
      // yet assigned a partition.
      std::vector<unsigned int> cell_partition(n_active_cells,
                                               numbers::invalid_unsigned_int);


      // In element j of this variable, one puts the old number of the block
      // that should be the jth block in the new numeration.
      std::vector<unsigned int> partition_list(n_active_cells, 0);
      std::vector<unsigned int> partition_partition_list(n_active_cells, 0);

      // This vector points to the start of each partition.
      std::vector<unsigned int> partition_size(2, 0);

      unsigned int partition = 0;
      // Here, we do not block inside the connectivity graph
      // blocking_connectivity = false;

      // Make the partitioning of the first layer of the blocks of cells.
      make_partitioning(connectivity,
                        cluster_size,
                        cell_partition,
                        partition_list,
                        partition_size,
                        partition);

      // Partition within partitions.
      make_partitioning_within_partitions_post_blocked(connectivity,
                                                       cell_active_fe_index,
                                                       partition,
                                                       cluster_size,
                                                       hp_bool,
                                                       cell_partition,
                                                       partition_list,
                                                       partition_size,
                                                       partition_partition_list,
                                                       irregular_cells);

      partition_list.swap(renumbering);

      for (unsigned int j = 0; j < renumbering.size(); ++j)
        renumbering[j] = partition_list[partition_partition_list[j]];

      for (unsigned int i = 0; i < n_ghost_cells; ++i)
        renumbering.push_back(i + n_active_cells);

      update_task_info(partition);
    }



    void
    TaskInfo::make_connectivity_cells_to_blocks(
      const std::vector<unsigned char> &irregular_cells,
      const DynamicSparsityPattern     &connectivity_cells,
      DynamicSparsityPattern           &connectivity_blocks) const
    {
      std::vector<std::vector<unsigned int>> cell_blocks(n_blocks);
      std::vector<unsigned int>              touched_cells(n_active_cells);
      unsigned int                           cell = 0;
      for (unsigned int i = 0, mcell = 0; i < n_blocks; ++i)
        {
          for (unsigned int c = 0;
               c < block_size && mcell < *(cell_partition_data.end() - 2);
               ++c, ++mcell)
            {
              unsigned int ncomp = (irregular_cells[mcell] > 0) ?
                                     irregular_cells[mcell] :
                                     vectorization_length;
              for (unsigned int c = 0; c < ncomp; ++c, ++cell)
                {
                  cell_blocks[i].push_back(cell);
                  touched_cells[cell] = i;
                }
            }
        }
      AssertDimension(cell, n_active_cells);
      for (unsigned int i = 0; i < cell_blocks.size(); ++i)
        for (unsigned int col = 0; col < cell_blocks[i].size(); ++col)
          {
            for (DynamicSparsityPattern::iterator it =
                   connectivity_cells.begin(cell_blocks[i][col]);
                 it != connectivity_cells.end(cell_blocks[i][col]);
                 ++it)
              {
                if (touched_cells[it->column()] != i)
                  connectivity_blocks.add(i, touched_cells[it->column()]);
              }
          }
    }



    // Function to create partitioning on the second layer within each
    // partition. Version without preblocking.
    void
    TaskInfo::make_partitioning_within_partitions_post_blocked(
      const DynamicSparsityPattern    &connectivity,
      const std::vector<unsigned int> &cell_active_fe_index,
      const unsigned int               partition,
      const unsigned int               cluster_size,
      const bool                       hp_bool,
      const std::vector<unsigned int> &cell_partition,
      const std::vector<unsigned int> &partition_list,
      const std::vector<unsigned int> &partition_size,
      std::vector<unsigned int>       &partition_partition_list,
      std::vector<unsigned char>      &irregular_cells)
    {
      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      const unsigned int n_ghost_slots =
        *(cell_partition_data.end() - 1) - n_cell_batches;

      // List of cells in previous partition
      std::vector<unsigned int> neighbor_list;
      // List of cells in current partition for use as neighbors in next
      // partition
      std::vector<unsigned int> neighbor_neighbor_list;

      std::vector<unsigned int> renumbering(n_active_cells);

      irregular_cells.back() = 0;
      irregular_cells.resize(n_active_cells + n_ghost_slots);

      unsigned int max_fe_index = 0;
      for (const unsigned int fe_index : cell_active_fe_index)
        max_fe_index = std::max(fe_index, max_fe_index);

      Assert(!hp_bool || cell_active_fe_index.size() == n_active_cells,
             ExcInternalError());

      {
        unsigned int n_cell_batches_before = 0;
        // Create partitioning within partitions.

        // For each block of cells, this variable saves to which partitions
        // the block belongs. Initialize all to n_cell_batches to mark them as
        // not yet assigned a partition.
        std::vector<unsigned int> cell_partition_l2(
          n_active_cells, numbers::invalid_unsigned_int);
        partition_row_index.clear();
        partition_row_index.resize(partition + 1, 0);
        cell_partition_data.resize(1, 0);

        unsigned int counter = 0;
        unsigned int missing_macros;
        for (unsigned int part = 0; part < partition; ++part)
          {
            neighbor_neighbor_list.resize(0);
            neighbor_list.resize(0);
            bool         work              = true;
            unsigned int partition_l2      = 0;
            unsigned int start_up          = partition_size[part];
            unsigned int partition_counter = 0;
            while (work)
              {
                if (neighbor_list.empty())
                  {
                    work              = false;
                    partition_counter = 0;
                    for (unsigned int j = start_up;
                         j < partition_size[part + 1];
                         ++j)
                      if (cell_partition[partition_list[j]] == part &&
                          cell_partition_l2[partition_list[j]] ==
                            numbers::invalid_unsigned_int)
                        {
                          start_up          = j;
                          work              = true;
                          partition_counter = 1;
                          // To start up, set the start_up cell to partition
                          // and list all its neighbors.
                          AssertIndexRange(start_up, partition_size[part + 1]);
                          cell_partition_l2[partition_list[start_up]] =
                            partition_l2;
                          neighbor_neighbor_list.push_back(
                            partition_list[start_up]);
                          partition_partition_list[counter++] =
                            partition_list[start_up];
                          ++start_up;
                          break;
                        }
                  }
                else
                  {
                    partition_counter = 0;
                    for (const unsigned int neighbor : neighbor_list)
                      {
                        Assert(cell_partition[neighbor] == part,
                               ExcInternalError());
                        Assert(cell_partition_l2[neighbor] == partition_l2 - 1,
                               ExcInternalError());
                        auto       neighbor_it = connectivity.begin(neighbor);
                        const auto end_it      = connectivity.end(neighbor);
                        for (; neighbor_it != end_it; ++neighbor_it)
                          {
                            if (cell_partition[neighbor_it->column()] == part &&
                                cell_partition_l2[neighbor_it->column()] ==
                                  numbers::invalid_unsigned_int)
                              {
                                cell_partition_l2[neighbor_it->column()] =
                                  partition_l2;
                                neighbor_neighbor_list.push_back(
                                  neighbor_it->column());
                                partition_partition_list[counter++] =
                                  neighbor_it->column();
                                ++partition_counter;
                              }
                          }
                      }
                  }
                if (partition_counter > 0)
                  {
                    int index_before = neighbor_neighbor_list.size(),
                        index        = index_before;
                    {
                      // put the cells into separate lists for each FE index
                      // within one partition-partition
                      missing_macros = 0;
                      std::vector<unsigned int> remaining_per_cell_batch(
                        max_fe_index + 1);
                      std::vector<std::vector<unsigned int>>
                                   renumbering_fe_index;
                      unsigned int cell;
                      bool         filled = true;
                      if (hp_bool == true)
                        {
                          renumbering_fe_index.resize(max_fe_index + 1);
                          for (cell = counter - partition_counter;
                               cell < counter;
                               ++cell)
                            {
                              renumbering_fe_index
                                [cell_active_fe_index.empty() ?
                                   0 :
                                   cell_active_fe_index
                                     [partition_partition_list[cell]]]
                                  .push_back(partition_partition_list[cell]);
                            }
                          // check how many more cells are needed in the lists
                          for (unsigned int j = 0; j < max_fe_index + 1; ++j)
                            {
                              remaining_per_cell_batch[j] =
                                renumbering_fe_index[j].size() %
                                vectorization_length;
                              if (remaining_per_cell_batch[j] != 0)
                                filled = false;
                              missing_macros +=
                                ((renumbering_fe_index[j].size() +
                                  vectorization_length - 1) /
                                 vectorization_length);
                            }
                        }
                      else
                        {
                          remaining_per_cell_batch.resize(1);
                          remaining_per_cell_batch[0] =
                            partition_counter % vectorization_length;
                          missing_macros =
                            partition_counter / vectorization_length;
                          if (remaining_per_cell_batch[0] != 0)
                            {
                              filled = false;
                              ++missing_macros;
                            }
                        }
                      missing_macros =
                        cluster_size - (missing_macros % cluster_size);

                      // now we realized that there are some cells missing.
                      while (missing_macros > 0 || filled == false)
                        {
                          if (index == 0)
                            {
                              index = neighbor_neighbor_list.size();
                              if (index == index_before)
                                {
                                  if (missing_macros != 0)
                                    {
                                      neighbor_neighbor_list.resize(0);
                                    }
                                  start_up--;
                                  break; // not connected - start again
                                }
                              index_before = index;
                            }
                          index--;
                          unsigned int additional =
                            neighbor_neighbor_list[index];

                          // go through the neighbors of the last cell in the
                          // current partition and check if we find some to
                          // fill up with.
                          DynamicSparsityPattern::iterator neighbor =
                                                             connectivity.begin(
                                                               additional),
                                                           end =
                                                             connectivity.end(
                                                               additional);
                          for (; neighbor != end; ++neighbor)
                            {
                              if (cell_partition[neighbor->column()] == part &&
                                  cell_partition_l2[neighbor->column()] ==
                                    numbers::invalid_unsigned_int)
                                {
                                  unsigned int this_index = 0;
                                  if (hp_bool == true)
                                    this_index =
                                      cell_active_fe_index.empty() ?
                                        0 :
                                        cell_active_fe_index[neighbor
                                                               ->column()];

                                  // Only add this cell if we need more macro
                                  // cells in the current block or if there is
                                  // a macro cell with the FE index that is
                                  // not yet fully populated
                                  if (missing_macros > 0 ||
                                      remaining_per_cell_batch[this_index] > 0)
                                    {
                                      cell_partition_l2[neighbor->column()] =
                                        partition_l2;
                                      neighbor_neighbor_list.push_back(
                                        neighbor->column());
                                      if (hp_bool == true)
                                        renumbering_fe_index[this_index]
                                          .push_back(neighbor->column());
                                      partition_partition_list[counter] =
                                        neighbor->column();
                                      ++counter;
                                      ++partition_counter;
                                      if (remaining_per_cell_batch
                                              [this_index] == 0 &&
                                          missing_macros > 0)
                                        missing_macros--;
                                      remaining_per_cell_batch[this_index]++;
                                      if (remaining_per_cell_batch
                                            [this_index] ==
                                          vectorization_length)
                                        {
                                          remaining_per_cell_batch[this_index] =
                                            0;
                                        }
                                      if (missing_macros == 0)
                                        {
                                          filled = true;
                                          for (unsigned int fe_ind = 0;
                                               fe_ind < max_fe_index + 1;
                                               ++fe_ind)
                                            if (remaining_per_cell_batch
                                                  [fe_ind] != 0)
                                              filled = false;
                                        }
                                      if (filled == true)
                                        break;
                                    }
                                }
                            }
                        }
                      if (hp_bool == true)
                        {
                          // set the renumbering according to their active FE
                          // index within one partition-partition which was
                          // implicitly assumed above
                          cell = counter - partition_counter;
                          for (unsigned int j = 0; j < max_fe_index + 1; ++j)
                            {
                              for (const unsigned int jj :
                                   renumbering_fe_index[j])
                                renumbering[cell++] = jj;
                              if (renumbering_fe_index[j].size() %
                                    vectorization_length !=
                                  0)
                                irregular_cells[renumbering_fe_index[j].size() /
                                                  vectorization_length +
                                                n_cell_batches_before] =
                                  renumbering_fe_index[j].size() %
                                  vectorization_length;
                              n_cell_batches_before +=
                                (renumbering_fe_index[j].size() +
                                 vectorization_length - 1) /
                                vectorization_length;
                              renumbering_fe_index[j].resize(0);
                            }
                        }
                      else
                        {
                          n_cell_batches_before +=
                            partition_counter / vectorization_length;
                          if (partition_counter % vectorization_length != 0)
                            {
                              irregular_cells[n_cell_batches_before] =
                                partition_counter % vectorization_length;
                              ++n_cell_batches_before;
                            }
                        }
                    }
                    cell_partition_data.push_back(n_cell_batches_before);
                    partition_l2++;
                  }
                neighbor_list = neighbor_neighbor_list;
                neighbor_neighbor_list.resize(0);
              }
            partition_row_index[part + 1] =
              partition_row_index[part] + partition_l2;
          }
      }
      if (hp_bool == true)
        {
          partition_partition_list.swap(renumbering);
        }
    }



    // Function to create coloring on the second layer within each partition.
    // Version assumes preblocking.
    void
    TaskInfo::make_coloring_within_partitions_pre_blocked(
      const DynamicSparsityPattern    &connectivity,
      const unsigned int               partition,
      const std::vector<unsigned int> &cell_partition,
      const std::vector<unsigned int> &partition_list,
      const std::vector<unsigned int> &partition_size,
      std::vector<unsigned int>       &partition_color_list)
    {
      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      std::vector<unsigned int> cell_color(n_blocks, n_cell_batches);
      std::vector<bool>         color_finder;

      partition_row_index.resize(partition + 1);
      cell_partition_data.clear();
      unsigned int color_counter = 0, index_counter = 0;
      for (unsigned int part = 0; part < partition; ++part)
        {
          partition_row_index[part] = index_counter;
          unsigned int max_color    = 0;
          for (unsigned int k = partition_size[part];
               k < partition_size[part + 1];
               k++)
            {
              unsigned int cell        = partition_list[k];
              unsigned int n_neighbors = connectivity.row_length(cell);

              // In the worst case, each neighbor has a different color. So we
              // find at least one available color between 0 and n_neighbors.
              color_finder.resize(n_neighbors + 1);
              for (unsigned int j = 0; j <= n_neighbors; ++j)
                color_finder[j] = true;
              DynamicSparsityPattern::iterator neighbor =
                                                 connectivity.begin(cell),
                                               end = connectivity.end(cell);
              for (; neighbor != end; ++neighbor)
                {
                  // Mark the color that a neighbor within the partition has
                  // as taken
                  if (cell_partition[neighbor->column()] == part &&
                      cell_color[neighbor->column()] <= n_neighbors)
                    color_finder[cell_color[neighbor->column()]] = false;
                }
              // Choose the smallest color that is not taken for the block
              cell_color[cell] = 0;
              while (color_finder[cell_color[cell]] == false)
                cell_color[cell]++;
              if (cell_color[cell] > max_color)
                max_color = cell_color[cell];
            }
          // Reorder within partition: First, all blocks that belong the 0 and
          // then so on until those with color max (Note that the smaller the
          // number the larger the partition)
          for (unsigned int color = 0; color <= max_color; ++color)
            {
              cell_partition_data.push_back(color_counter);
              ++index_counter;
              for (unsigned int k = partition_size[part];
                   k < partition_size[part + 1];
                   k++)
                {
                  unsigned int cell = partition_list[k];
                  if (cell_color[cell] == color)
                    {
                      partition_color_list[color_counter++] = cell;
                    }
                }
            }
        }
      cell_partition_data.push_back(n_blocks);
      partition_row_index[partition] = index_counter;
      AssertDimension(color_counter, n_blocks);
    }


    // Function to create partitioning on the first layer.
    void
    TaskInfo::make_partitioning(const DynamicSparsityPattern &connectivity,
                                const unsigned int            cluster_size,
                                std::vector<unsigned int>    &cell_partition,
                                std::vector<unsigned int>    &partition_list,
                                std::vector<unsigned int>    &partition_size,
                                unsigned int                 &partition) const

    {
      // For each block of cells, this variable saves to which partitions the
      // block belongs. Initialize all to n_cell_batches to mark them as not
      // yet assigned a partition.
      // std::vector<unsigned int> cell_partition (n_active_cells,
      //                                          numbers::invalid_unsigned_int);
      // List of cells in previous partition
      std::vector<unsigned int> neighbor_list;
      // List of cells in current partition for use as neighbors in next
      // partition
      std::vector<unsigned int> neighbor_neighbor_list;

      // In element j of this variable, one puts the old number of the block
      // that should be the jth block in the new numeration.
      // std::vector<unsigned int> partition_list(n_active_cells,0);

      // This vector points to the start of each partition.
      // std::vector<unsigned int> partition_size(2,0);

      partition            = 0;
      unsigned int counter = 0;
      unsigned int start_nonboundary =
        cell_partition_data.size() == 5 ?
          vectorization_length *
            (cell_partition_data[2] - cell_partition_data[1]) :
          0;

      const unsigned int n_cell_batches = *(cell_partition_data.end() - 2);
      if (n_cell_batches == 0)
        return;
      if (scheme == color)
        start_nonboundary = n_cell_batches;
      if (scheme == partition_color ||
          scheme == color) // blocking_connectivity == true
        start_nonboundary = ((start_nonboundary + block_size - 1) / block_size);
      unsigned int n_blocks;
      if (scheme == partition_color ||
          scheme == color) // blocking_connectivity == true
        n_blocks = this->n_blocks;
      else
        n_blocks = n_active_cells;

      if (start_nonboundary > n_blocks)
        start_nonboundary = n_blocks;


      unsigned int start_up  = 0;
      bool         work      = true;
      unsigned int remainder = cluster_size;

      // this performs a classical breath-first search in the connectivity
      // graph of the cells under the restriction that the size of the
      // partitions should be a multiple of the given block size
      while (work)
        {
          // put the cells with neighbors on remote MPI processes up front
          if (start_nonboundary > 0)
            {
              for (unsigned int cell = 0; cell < start_nonboundary; ++cell)
                {
                  const unsigned int cell_nn = cell;
                  cell_partition[cell_nn]    = partition;
                  neighbor_list.push_back(cell_nn);
                  partition_list[counter++] = cell_nn;
                  partition_size.back()++;
                }
              start_nonboundary = 0;
              remainder -= (start_nonboundary % cluster_size);
              if (remainder == cluster_size)
                remainder = 0;
            }
          else
            {
              // To start up, set the start_up cell to partition and list all
              // its neighbors.
              cell_partition[start_up] = partition;
              neighbor_list.push_back(start_up);
              partition_list[counter++] = start_up;
              partition_size.back()++;
              ++start_up;
              remainder--;
              if (remainder == cluster_size)
                remainder = 0;
            }
          int index_before = neighbor_list.size(), index = index_before,
              index_stop = 0;
          while (remainder > 0)
            {
              if (index == index_stop)
                {
                  index = neighbor_list.size();
                  if (index == index_before)
                    {
                      neighbor_list.resize(0);
                      goto not_connect;
                    }
                  index_stop   = index_before;
                  index_before = index;
                }
              index--;
              unsigned int additional = neighbor_list[index];
              DynamicSparsityPattern::iterator neighbor =
                                                 connectivity.begin(additional),
                                               end =
                                                 connectivity.end(additional);
              for (; neighbor != end; ++neighbor)
                {
                  if (cell_partition[neighbor->column()] ==
                      numbers::invalid_unsigned_int)
                    {
                      partition_size.back()++;
                      cell_partition[neighbor->column()] = partition;
                      neighbor_list.push_back(neighbor->column());
                      partition_list[counter++] = neighbor->column();
                      remainder--;
                      if (remainder == 0)
                        break;
                    }
                }
            }

          while (neighbor_list.size() > 0)
            {
              ++partition;

              // counter for number of cells so far in current partition
              unsigned int partition_counter = 0;

              // Mark the start of the new partition
              partition_size.push_back(partition_size.back());

              // Loop through the list of cells in previous partition and put
              // all their neighbors in current partition
              for (const unsigned int cell : neighbor_list)
                {
                  Assert(cell_partition[cell] == partition - 1,
                         ExcInternalError());
                  auto       neighbor = connectivity.begin(cell);
                  const auto end      = connectivity.end(cell);
                  for (; neighbor != end; ++neighbor)
                    {
                      if (cell_partition[neighbor->column()] ==
                          numbers::invalid_unsigned_int)
                        {
                          partition_size.back()++;
                          cell_partition[neighbor->column()] = partition;

                          // collect the cells of the current partition for
                          // use as neighbors in next partition
                          neighbor_neighbor_list.push_back(neighbor->column());
                          partition_list[counter++] = neighbor->column();
                          ++partition_counter;
                        }
                    }
                }
              remainder = cluster_size - (partition_counter % cluster_size);
              if (remainder == cluster_size)
                remainder = 0;
              int index_stop   = 0;
              int index_before = neighbor_neighbor_list.size(),
                  index        = index_before;
              while (remainder > 0)
                {
                  if (index == index_stop)
                    {
                      index = neighbor_neighbor_list.size();
                      if (index == index_before)
                        {
                          neighbor_neighbor_list.resize(0);
                          break;
                        }
                      index_stop   = index_before;
                      index_before = index;
                    }
                  index--;
                  unsigned int additional = neighbor_neighbor_list[index];
                  DynamicSparsityPattern::iterator neighbor =
                                                     connectivity.begin(
                                                       additional),
                                                   end = connectivity.end(
                                                     additional);
                  for (; neighbor != end; ++neighbor)
                    {
                      if (cell_partition[neighbor->column()] ==
                          numbers::invalid_unsigned_int)
                        {
                          partition_size.back()++;
                          cell_partition[neighbor->column()] = partition;
                          neighbor_neighbor_list.push_back(neighbor->column());
                          partition_list[counter++] = neighbor->column();
                          remainder--;
                          if (remainder == 0)
                            break;
                        }
                    }
                }

              neighbor_list = neighbor_neighbor_list;
              neighbor_neighbor_list.resize(0);
            }
        not_connect:
          // One has to check if the graph is not connected so we have to find
          // another partition.
          work = false;
          for (unsigned int j = start_up; j < n_blocks; ++j)
            if (cell_partition[j] == numbers::invalid_unsigned_int)
              {
                start_up = j;
                work     = true;
                if (remainder == 0)
                  remainder = cluster_size;
                break;
              }
        }
      if (remainder != 0)
        ++partition;

      AssertDimension(partition_size[partition], n_blocks);
    }


    void
    TaskInfo::update_task_info(const unsigned int partition)
    {
      evens             = (partition + 1) / 2;
      odds              = partition / 2;
      n_blocked_workers = odds - (odds + evens + 1) % 2;
      n_workers         = evens + odds - n_blocked_workers;
      // From here only used for partition partition option.
      partition_evens.resize(partition);
      partition_odds.resize(partition);
      partition_n_blocked_workers.resize(partition);
      partition_n_workers.resize(partition);
      for (unsigned int part = 0; part < partition; ++part)
        {
          partition_evens[part] =
            (partition_row_index[part + 1] - partition_row_index[part] + 1) / 2;
          partition_odds[part] =
            (partition_row_index[part + 1] - partition_row_index[part]) / 2;
          partition_n_blocked_workers[part] =
            partition_odds[part] -
            (partition_odds[part] + partition_evens[part] + 1) % 2;
          partition_n_workers[part] = partition_evens[part] +
                                      partition_odds[part] -
                                      partition_n_blocked_workers[part];
        }
    }
  } // namespace MatrixFreeFunctions
} // namespace internal



// explicit instantiations of template functions
template void
internal::MatrixFreeFunctions::TaskInfo::print_memory_statistics<std::ostream>(
  std::ostream &,
  const std::size_t) const;
template void
internal::MatrixFreeFunctions::TaskInfo::print_memory_statistics<
  ConditionalOStream>(ConditionalOStream &, const std::size_t) const;


DEAL_II_NAMESPACE_CLOSE
