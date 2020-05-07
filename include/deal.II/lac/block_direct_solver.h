/*
 * block_direct_solver.h
 *
 *  Created on: Mar 18, 2020
 *      Author: mwichro
 */

#ifndef INCLUDE_BLOCK_DIRECT_SOLVER_H_
#define INCLUDE_BLOCK_DIRECT_SOLVER_H_

#include <deal.II/base/index_set.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

#include <deal.II/lac/trilinos_parallel_block_vector.h>
#include <deal.II/lac/trilinos_precondition.h>
#include <deal.II/lac/trilinos_solver.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/multigrid/mg_coarse.h>

#include <deal.II/numerics/vector_tools.h>

using namespace dealii;

template <typename LevelNumber>
class BlockDirectSolver
  : public MGCoarseGridBase<
      LinearAlgebra::distributed::BlockVector<LevelNumber>>
{
public:
  friend class BlockCoarseTest;
  typedef LevelNumber                                          number;
  typedef LinearAlgebra::distributed::BlockVector<LevelNumber> BlockVectorType;


  // Constructor. Requires owned_partitioning[block][MPI process],
  // relevant_partitioning[block] and MPI communicator.
  // The owned partitioning for each block may be obtained by using
  // DoFHandler::compute_locally_owned_dofs_per_processor() or in case of MG
  // DoFHandler::compute_locally_owned_mg_dofs_per_processor(level)


  BlockDirectSolver(
    const std::vector<std::vector<IndexSet>> &owned_partitioning_,
    const std::vector<IndexSet> &             relevant_dofs_,
    const MPI_Comm &                          communicator);

  // Computes TrilinosWrappers::SparseMatrix using internal DoF renumbering,
  // initializes direct solver (TrilinosWrappers::SolverDirect)
  // Level will be checked if  operator() is used
  template <class Operator>
  void
  initialize(const Operator &vmult_operator, const unsigned int &lvl = 0);
  //
  //	template<class Operator, int dim>
  //		void initialize(const Operator & vmult_operator,
  //				const std::vector<DoFHandler<dim>*> &dof_handlers,
  //				const unsigned int& lvl=0);
  //


  void
  vmult(BlockVectorType &dst, const BlockVectorType &src) const;
  // Interface to act as a coarse grid solver.
  // Applies inverse of operator used in initialize.
  virtual void
  operator()(const unsigned int     lvl,
             BlockVectorType &      dst,
             const BlockVectorType &src) const;

private:
  typedef TrilinosWrappers::MPI::Vector InternalVectorType;

  const MPI_Comm     mpi_comm;
  const unsigned int this_mpi_process;
  const unsigned int n_mpi_process;

  int                level;
  const unsigned int n_blocks;

  types::global_dof_index n_local_dofs;

  std::vector<IndexSet> owned_dofs;
  //[block][process]


  const std::vector<std::vector<IndexSet>> owned_partitioning;
  const std::vector<IndexSet>              relevant_dofs;

  // Searching for owning processor is in order specified by following list
  mutable std::list<unsigned int> searching_order;

  // Shitfs of dof indices, used for internal renumbering
  std::vector<std::vector<types::global_dof_index>> dof_shifts;

  std::vector<types::global_dof_index> first_block_index;
  std::vector<types::global_dof_index> last_block_index;

  // IndexSet of reordered owned dofs.
  IndexSet shifted_owned;

  std::unique_ptr<TrilinosWrappers::SparseMatrix> matrix;
  SolverControl                                   solver_control;
  mutable TrilinosWrappers::SolverDirect          inverse;


  // Computes internal DoF number from external (dof, block) numbering.
  // First checks if dof is locally owned, then searches other index sets, in
  // order
  // specified by searching_order.
  types::global_dof_index
  ext2int(const types::global_dof_index &external_dof,
          const unsigned int             block) const;


  // Not implemented, not needed.
//  std::pair<types::global_dof_index, unsigned int>
//  int2ext(types::global_dof_index &internal_dof) const;
};

template <typename LevelNumber>
BlockDirectSolver<LevelNumber>::BlockDirectSolver(
  const std::vector<std::vector<IndexSet>> &owned_partitioning_,
  const std::vector<IndexSet> &             relevant_dofs_,
  const MPI_Comm &                          communicator)
  : mpi_comm(communicator)
  , this_mpi_process(Utilities::MPI::this_mpi_process(communicator))
  , n_mpi_process(Utilities::MPI::n_mpi_processes(mpi_comm))
  , level(-1)
  , n_blocks(relevant_dofs_.size())
  , n_local_dofs(0)
  , owned_partitioning(owned_partitioning_)
  , relevant_dofs(relevant_dofs_)
  , solver_control(10, 1e-14)
  , inverse(solver_control)
{
  AssertDimension(owned_partitioning.size(), n_blocks);
  owned_dofs.resize(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      AssertDimension(owned_partitioning[b].size(), n_mpi_process);
      owned_dofs[b] = owned_partitioning[b][this_mpi_process];

      Assert(owned_dofs[b].is_ascending_and_one_to_one(mpi_comm),
             ExcMessage("Owned sets have to ve  ascending and one to one!"));
    }

  searching_order.clear();
  for (unsigned int i = 0; i < n_mpi_process; ++i)
    {
      for (unsigned int b = 0; b < n_blocks; ++b)
        {
          Assert(owned_partitioning[b][i].is_contiguous(),
                 ExcMessage("Owned index set not contiguous"));
        }
      searching_order.push_front(i);
    }

  AssertDimension(searching_order.size(), n_mpi_process);

  Assert(0 != n_blocks, ExcEmptyObject());
  dof_shifts.resize(n_mpi_process);
  for (unsigned int p = 0; p < n_mpi_process; ++p)
    dof_shifts[p].resize(n_blocks);


  std::vector<types::global_dof_index> starts(n_mpi_process, 0);
  types::global_dof_index              n_dofs = 0;



  for (unsigned int i = 1; i < n_mpi_process; ++i)
    {
      starts[i] = starts[i - 1];

      for (unsigned int b = 0; b < n_blocks; ++b)
        {
          starts[i] += owned_partitioning[b][i - 1].n_elements();
        }
    }


  for (unsigned int b = 0; b < n_blocks; ++b)
    n_dofs += owned_partitioning[b][this_mpi_process].size();


  for (unsigned int p = 0; p < n_mpi_process; ++p)
    {
      dof_shifts[p][0] =
        starts[p] - owned_partitioning[0][p].nth_index_in_set(0);
      for (unsigned int b = 1; b < n_blocks; ++b)
        dof_shifts[p][b] = owned_partitioning[b - 1][p].nth_index_in_set(0) +
                           dof_shifts[p][b - 1] +
                           owned_partitioning[b - 1][p].n_elements() -
                           owned_partitioning[b][p].nth_index_in_set(0);
    }

  shifted_owned.set_size(n_dofs);

  const unsigned int begin_range =
    dof_shifts[this_mpi_process][0] + owned_dofs[0].nth_index_in_set(0);
  const unsigned int end_range = owned_dofs[n_blocks - 1].nth_index_in_set(0) +
                                 dof_shifts[this_mpi_process][n_blocks - 1] +
                                 owned_dofs[n_blocks - 1].n_elements();

  shifted_owned.add_range(begin_range, end_range);


  n_local_dofs = 0;
  for (unsigned int b = 0; b < n_blocks; ++b)
    n_local_dofs += owned_dofs[b].n_elements();

  AssertDimension(n_local_dofs, shifted_owned.n_elements());

  Assert(shifted_owned.is_ascending_and_one_to_one(mpi_comm),
         ExcMessage(
           "The shifted dofs are not valid index set (internal error)"));


  first_block_index.resize(n_blocks);
  last_block_index.resize(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      first_block_index[b] =
        dof_shifts[this_mpi_process][b] + owned_dofs[b].nth_index_in_set(0);
      last_block_index[b] =
        first_block_index[b] + owned_dofs[b].n_elements() - 1;
    }

  AssertDimension(last_block_index[n_blocks - 1],
                  first_block_index[0] + n_local_dofs - 1);
}



template <typename LevelNumber>
template <class Operator>
void
BlockDirectSolver<LevelNumber>::initialize(const Operator &    vmult_operator,
                                           const unsigned int &lvl)
{
  matrix = std::make_unique<TrilinosWrappers::SparseMatrix>(shifted_owned,
                                                            shifted_owned,
                                                            mpi_comm,
                                                            100);
  level  = lvl;

  BlockVectorType dst;
  BlockVectorType src;

  src.reinit(n_blocks);
  dst.reinit(n_blocks);
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      src.block(b).reinit(owned_dofs[b], relevant_dofs[b], mpi_comm);
      dst.block(b).reinit(owned_dofs[b], relevant_dofs[b], mpi_comm);
    }
  src.collect_sizes();
  dst.collect_sizes();

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (types::global_dof_index i = 0; i < owned_dofs[b].size(); ++i)
        {
          src = 0;
          dst = 0;
          if (owned_dofs[b].is_element(i))
            {
              src.block(b)(i) = 1;
            }

          src.compress(VectorOperation::insert);
          vmult_operator.vmult(dst, src);

          dst.update_ghost_values();


          for (unsigned int k = 0; k < n_blocks; ++k)
            {
              for (IndexSet::ElementIterator index_iter = owned_dofs[k].begin();
                   index_iter != owned_dofs[k].end();
                   ++index_iter)
                {
                  if (dst.block(k)(*index_iter) != 0)
                    matrix->set(ext2int(*index_iter, k),
                                ext2int(i, b),
                                dst.block(k)(*index_iter));
                }
            }
        }
    }
  matrix->compress(VectorOperation::unknown);
  inverse.initialize(*matrix);
}


template <typename LevelNumber>
types::global_dof_index
BlockDirectSolver<LevelNumber>::ext2int(
  const types::global_dof_index &external_dof,
  const unsigned int             block) const
{
  // the provided dof in in our local range
  if (owned_dofs[block].is_element(external_dof))
    {
      const types::global_dof_index result =
        external_dof + dof_shifts[this_mpi_process][block];
      Assert(shifted_owned.is_element(result),
             ExcMessage(
               "The computed internal DoF index is not in owned range"));
      return result;
    }
  // the provided dof is in "ghost range", search for the right process
  // we are going to access dofs owned by only few neigbour MPI processes
  // We keep the search orders of IndexSets as a list
  // and bring to front index set that we will find.
  // assuming that the ghost range of vector is shared with low number of
  // processors
  // this is expected to be more efficient that binary search

  AssertDimension(n_mpi_process, searching_order.size());
  std::list<unsigned int>::iterator iter;
  for (iter = searching_order.begin(); iter != searching_order.end(); ++iter)
    if (owned_partitioning[block][*iter].is_element(external_dof))
      {
        // compute internal index and return
        const unsigned int            proc = *iter;
        const types::global_dof_index result =
          external_dof + dof_shifts[proc][block];

        searching_order.erase(iter);
        searching_order.push_front(proc);

        return result;
      }


  Assert(false, ExcInternalError());
  return -1;
}



//template <typename LevelNumber>
//std::pair<types::global_dof_index, unsigned int>
//BlockDirectSolver<LevelNumber>::int2ext(
//  types::global_dof_index &internal_dof) const
//{
//  Assert(shifted_owned.is_element(internal_dof),
//         ExcMessage("The computed internal DoF index is not in owned range"));
//  Assert(false, ExcNotImplemented());
//}



template <typename LevelNumber>
void
BlockDirectSolver<LevelNumber>::vmult(BlockVectorType &      dst,
                                      const BlockVectorType &src) const
{
  InternalVectorType src_internal(shifted_owned, mpi_comm);
  InternalVectorType dst_internal(shifted_owned, mpi_comm);


  // Reorder src vector:
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (IndexSet::ElementIterator index_iter = owned_dofs[b].begin();
           index_iter != owned_dofs[b].end();
           ++index_iter)
        {
          src_internal(ext2int(*index_iter, b)) = src.block(b)(*index_iter);
        }
    }
  src_internal.compress(VectorOperation::insert);

  // apply direct solver
  inverse.solve(dst_internal, src_internal);

  // fill the dst vector
  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (IndexSet::ElementIterator index_iter = owned_dofs[b].begin();
           index_iter != owned_dofs[b].end();
           ++index_iter)
        {
          dst.block(b)(*index_iter) = dst_internal(ext2int(*index_iter, b));
        }
    }
  dst.compress(VectorOperation::insert);
}

template <typename LevelNumber>
void
BlockDirectSolver<LevelNumber>::operator()(const unsigned int     lvl,
                                           BlockVectorType &      dst,
                                           const BlockVectorType &src) const
{
  AssertDimension(lvl, level);
  (void) lvl;
  this->vmult(dst, src);
}



#endif /* INCLUDE_BLOCK_DIRECT_SOLVER_H_ */
