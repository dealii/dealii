// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#include <deal.II/lac/petsc_block_sparse_matrix.h>
#include <deal.II/lac/petsc_compatibility.h>

#ifdef DEAL_II_WITH_PETSC

DEAL_II_NAMESPACE_OPEN

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
      PetscErrorCode ierr = destroy_matrix(petsc_nest_matrix);
      (void)ierr;
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
    BlockSparseMatrix::reinit(const std::vector<IndexSet> &      rows,
                              const std::vector<IndexSet> &      cols,
                              const BlockDynamicSparsityPattern &bdsp,
                              const MPI_Comm &                   com)
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

      collect_sizes();
    }

    void
    BlockSparseMatrix::reinit(const std::vector<IndexSet> &      sizes,
                              const BlockDynamicSparsityPattern &bdsp,
                              const MPI_Comm &                   com)
    {
      reinit(sizes, sizes, bdsp, com);
    }

    void
    BlockSparseMatrix::collect_sizes()
    {
      BaseClass::collect_sizes();

      auto m = this->n_block_cols();
      auto n = this->n_block_cols();

      PetscErrorCode ierr = destroy_matrix(petsc_nest_matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
      std::vector<Mat> psub_objects(m * n);
      for (unsigned int r = 0; r < m; r++)
        for (unsigned int c = 0; c < n; c++)
          psub_objects[r * n + c] = this->sub_objects[r][c]->petsc_matrix();

      MPI_Comm comm =
        psub_objects.size() > 0 ?
          PetscObjectComm(reinterpret_cast<PetscObject>(psub_objects[0])) :
          PETSC_COMM_SELF;

      ierr = MatCreateNest(
        comm, m, nullptr, n, nullptr, psub_objects.data(), &petsc_nest_matrix);
      AssertThrow(ierr == 0, ExcPETScError(ierr));
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

    const MPI_Comm &
    BlockSparseMatrix::get_mpi_communicator() const
    {
      this->returncomm =
        PetscObjectComm(reinterpret_cast<PetscObject>(petsc_nest_matrix));
      return this->returncomm;
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
              // TODO: MATNEST supports NULL blocks
              this->sub_objects[i][j] = new BlockType(mats[i * nc + j]);
            }
        }

      collect_sizes();
    }

  } // namespace MPI
} // namespace PETScWrappers



DEAL_II_NAMESPACE_CLOSE

#endif
