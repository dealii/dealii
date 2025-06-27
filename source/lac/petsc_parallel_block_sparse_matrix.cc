// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2004 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/mpi_stub.h>

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_compatibility.h>


#ifdef DEAL_II_WITH_PETSC

#  include <petscmat.h>


DEAL_II_NAMESPACE_OPEN

namespace
{
  // A dummy utility routine to create an empty matrix in case we import
  // a MATNEST with NULL blocks
  static Mat
  create_dummy_mat(MPI_Comm comm,
                   PetscInt lr,
                   PetscInt gr,
                   PetscInt lc,
                   PetscInt gc)
  {
    Mat            dummy;
    PetscErrorCode ierr;

    ierr = MatCreate(comm, &dummy);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatSetSizes(dummy, lr, lc, gr, gc);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatSetType(dummy, MATAIJ);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatSeqAIJSetPreallocation(dummy, 0, nullptr);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatMPIAIJSetPreallocation(dummy, 0, nullptr, 0, nullptr);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatSetUp(dummy);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatSetOption(dummy, MAT_NO_OFF_PROC_ENTRIES, PETSC_TRUE);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatAssemblyBegin(dummy, MAT_FINAL_ASSEMBLY);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    ierr = MatAssemblyEnd(dummy, MAT_FINAL_ASSEMBLY);
    AssertThrow(ierr == 0, dealii::ExcPETScError(ierr));
    return dummy;
  }
} // namespace


namespace PETScWrappers
{
  namespace MPI
  {
    BlockSparseMatrix &
    BlockSparseMatrix::operator=(const BlockSparseMatrix &m)
    {
      BaseClass::operator=(m);

      return *this;
    }



    BlockSparseMatrix::~BlockSparseMatrix()
    {
      PetscErrorCode ierr = MatDestroy(&petsc_nest_matrix);
      AssertNothrow(ierr == 0, ExcPETScError(ierr));
    }



#  ifndef DOXYGEN
    void
    BlockSparseMatrix::reinit(const size_type n_block_rows,
                              const size_type n_block_columns)
    {
      // first delete previous content of
      // the subobjects array
      clear();

      // then resize. set sizes of blocks to
      // zero. user will later have to call
      // collect_sizes for this
      this->sub_objects.reinit(n_block_rows, n_block_columns);
      this->row_block_indices.reinit(n_block_rows, 0);
      this->column_block_indices.reinit(n_block_columns, 0);

      // and reinitialize the blocks
      for (size_type r = 0; r < this->n_block_rows(); ++r)
        for (size_type c = 0; c < this->n_block_cols(); ++c)
          {
            BlockType *p            = new BlockType();
            this->sub_objects[r][c] = p;
          }
    }
#  endif



    void
    BlockSparseMatrix::reinit(const std::vector<IndexSet>       &rows,
                              const std::vector<IndexSet>       &cols,
                              const BlockDynamicSparsityPattern &bdsp,
                              const MPI_Comm                     com)
    {
      Assert(rows.size() == bdsp.n_block_rows(), ExcMessage("invalid size"));
      Assert(cols.size() == bdsp.n_block_cols(), ExcMessage("invalid size"));


      clear();
      this->sub_objects.reinit(bdsp.n_block_rows(), bdsp.n_block_cols());

      std::vector<types::global_dof_index> row_sizes;
      for (unsigned int r = 0; r < bdsp.n_block_rows(); ++r)
        row_sizes.push_back(bdsp.block(r, 0).n_rows());
      this->row_block_indices.reinit(row_sizes);

      std::vector<types::global_dof_index> col_sizes;
      for (unsigned int c = 0; c < bdsp.n_block_cols(); ++c)
        col_sizes.push_back(bdsp.block(0, c).n_cols());
      this->column_block_indices.reinit(col_sizes);

      for (unsigned int r = 0; r < this->n_block_rows(); ++r)
        for (unsigned int c = 0; c < this->n_block_cols(); ++c)
          {
            Assert(rows[r].size() == bdsp.block(r, c).n_rows(),
                   ExcMessage("invalid size"));
            Assert(cols[c].size() == bdsp.block(r, c).n_cols(),
                   ExcMessage("invalid size"));

            BlockType *p = new BlockType();
            p->reinit(rows[r], cols[c], bdsp.block(r, c), com);
            this->sub_objects[r][c] = p;
          }

      this->collect_sizes();
    }

    void
    BlockSparseMatrix::reinit(const std::vector<IndexSet>       &sizes,
                              const BlockDynamicSparsityPattern &bdsp,
                              const MPI_Comm                     com)
    {
      reinit(sizes, sizes, bdsp, com);
    }



    void
    BlockSparseMatrix::create_empty_matrices_if_needed()
    {
      auto           m = this->n_block_rows();
      auto           n = this->n_block_cols();
      PetscErrorCode ierr;

      // Create empty matrices if needed
      // This is needed by the base class
      // not by MATNEST
      std::vector<size_type> row_sizes(m, size_type(-1));
      std::vector<size_type> col_sizes(n, size_type(-1));
      std::vector<size_type> row_local_sizes(m, size_type(-1));
      std::vector<size_type> col_local_sizes(n, size_type(-1));
      MPI_Comm               comm = MPI_COMM_NULL;
      for (size_type r = 0; r < m; r++)
        {
          for (size_type c = 0; c < n; c++)
            {
              if (this->sub_objects[r][c])
                {
                  comm = this->sub_objects[r][c]->get_mpi_communicator();
                  row_sizes[r]       = this->sub_objects[r][c]->m();
                  col_sizes[c]       = this->sub_objects[r][c]->n();
                  row_local_sizes[r] = this->sub_objects[r][c]->local_size();
                  col_local_sizes[c] =
                    this->sub_objects[r][c]->local_domain_size();
                }
            }
        }
      for (size_type r = 0; r < m; r++)
        {
          for (size_type c = 0; c < n; c++)
            {
              if (!this->sub_objects[r][c])
                {
                  Assert(
                    row_sizes[r] != size_type(-1),
                    ExcMessage(
                      "When passing empty sub-blocks of a block matrix, you need to make "
                      "sure that at least one block in each block row and block column is "
                      "non-empty. However, block row " +
                      std::to_string(r) +
                      " is completely empty "
                      "and so it is not possible to determine how many rows it should have."));
                  Assert(
                    col_sizes[c] != size_type(-1),
                    ExcMessage(
                      "When passing empty sub-blocks of a block matrix, you need to make "
                      "sure that at least one block in each block row and block column is "
                      "non-empty. However, block column " +
                      std::to_string(c) +
                      " is completely empty "
                      "and so it is not possible to determine how many columns it should have."));
                  Mat dummy =
                    create_dummy_mat(comm,
                                     static_cast<PetscInt>(row_local_sizes[r]),
                                     static_cast<PetscInt>(row_sizes[r]),
                                     static_cast<PetscInt>(col_local_sizes[c]),
                                     static_cast<PetscInt>(col_sizes[c]));
                  this->sub_objects[r][c] = new BlockType(dummy);

                  // the new object got a reference on dummy, we can safely
                  // call destroy here
                  ierr = MatDestroy(&dummy);
                  AssertThrow(ierr == 0, ExcPETScError(ierr));
                }
            }
        }
    }


    void
    BlockSparseMatrix::collect_sizes()
    {
      this->create_empty_matrices_if_needed();
      BaseClass::collect_sizes();
      this->setup_nest_mat();
    }

    void
    BlockSparseMatrix::setup_nest_mat()
    {
      auto           m = this->n_block_rows();
      auto           n = this->n_block_cols();
      PetscErrorCode ierr;

      MPI_Comm comm = PETSC_COMM_SELF;

      ierr = MatDestroy(&petsc_nest_matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      std::vector<Mat> psub_objects(m * n);
      for (unsigned int r = 0; r < m; r++)
        for (unsigned int c = 0; c < n; c++)
          {
            comm = this->sub_objects[r][c]->get_mpi_communicator();
            psub_objects[r * n + c] = this->sub_objects[r][c]->petsc_matrix();
          }
      ierr = MatCreateNest(
        comm, m, nullptr, n, nullptr, psub_objects.data(), &petsc_nest_matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));

      ierr = MatNestSetVecType(petsc_nest_matrix, VECNEST);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
    }



    void
    BlockSparseMatrix::compress(VectorOperation::values operation)
    {
      BaseClass::compress(operation);
      petsc_increment_state_counter(petsc_nest_matrix);
    }



    std::vector<IndexSet>
    BlockSparseMatrix::locally_owned_domain_indices() const
    {
      std::vector<IndexSet> index_sets;

      for (unsigned int i = 0; i < this->n_block_cols(); ++i)
        index_sets.push_back(this->block(0, i).locally_owned_domain_indices());

      return index_sets;
    }



    std::vector<IndexSet>
    BlockSparseMatrix::locally_owned_range_indices() const
    {
      std::vector<IndexSet> index_sets;

      for (unsigned int i = 0; i < this->n_block_rows(); ++i)
        index_sets.push_back(this->block(i, 0).locally_owned_range_indices());

      return index_sets;
    }



    std::uint64_t
    BlockSparseMatrix::n_nonzero_elements() const
    {
      std::uint64_t n_nonzero = 0;
      for (size_type rows = 0; rows < this->n_block_rows(); ++rows)
        for (size_type cols = 0; cols < this->n_block_cols(); ++cols)
          n_nonzero += this->block(rows, cols).n_nonzero_elements();

      return n_nonzero;
    }



    MPI_Comm
    BlockSparseMatrix::get_mpi_communicator() const
    {
      return PetscObjectComm(reinterpret_cast<PetscObject>(petsc_nest_matrix));
    }

    BlockSparseMatrix::operator const Mat &() const
    {
      return petsc_nest_matrix;
    }



    Mat &
    BlockSparseMatrix::petsc_matrix()
    {
      return petsc_nest_matrix;
    }

    void
    BlockSparseMatrix::reinit(Mat A)
    {
      clear();

      PetscBool isnest;
      PetscInt  nr = 1, nc = 1;

      PetscErrorCode ierr =
        PetscObjectTypeCompare(reinterpret_cast<PetscObject>(A),
                               MATNEST,
                               &isnest);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      std::vector<Mat> mats;
      bool             need_empty_matrices = false;
      if (isnest)
        {
          ierr = MatNestGetSize(A, &nr, &nc);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          for (PetscInt i = 0; i < nr; ++i)
            {
              for (PetscInt j = 0; j < nc; ++j)
                {
                  Mat sA;
                  ierr = MatNestGetSubMat(A, i, j, &sA);
                  mats.push_back(sA);
                  if (!sA)
                    need_empty_matrices = true;
                }
            }
        }
      else
        {
          mats.push_back(A);
        }

      std::vector<size_type> r_block_sizes(nr, 0);
      std::vector<size_type> c_block_sizes(nc, 0);
      this->row_block_indices.reinit(r_block_sizes);
      this->column_block_indices.reinit(c_block_sizes);
      this->sub_objects.reinit(nr, nc);
      for (PetscInt i = 0; i < nr; ++i)
        {
          for (PetscInt j = 0; j < nc; ++j)
            {
              if (mats[i * nc + j])
                this->sub_objects[i][j] = new BlockType(mats[i * nc + j]);
              else
                this->sub_objects[i][j] = nullptr;
            }
        }
      if (need_empty_matrices)
        this->create_empty_matrices_if_needed();

      BaseClass::collect_sizes();
      if (need_empty_matrices || !isnest)
        {
          setup_nest_mat();
        }
      else
        {
          ierr = PetscObjectReference(reinterpret_cast<PetscObject>(A));
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          PetscErrorCode ierr = MatDestroy(&petsc_nest_matrix);
          AssertThrow(ierr == 0, ExcPETScError(ierr));
          petsc_nest_matrix = A;
        }
    }

  } // namespace MPI
} // namespace PETScWrappers



DEAL_II_NAMESPACE_CLOSE

#endif
