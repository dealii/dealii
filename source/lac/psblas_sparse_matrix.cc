// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include "deal.II/base/exception_macros.h"
#include "deal.II/base/exceptions.h"
#include "deal.II/base/index_set.h"
#include "deal.II/base/mpi.h"
#include "deal.II/base/tensor.h"
#include "deal.II/base/types.h"

#include "deal.II/lac/dynamic_sparsity_pattern.h"
#include "deal.II/lac/exceptions.h"
#include "deal.II/lac/psblas_common.h"
#include <deal.II/lac/psblas_sparse_matrix.h>

#include <muParserDef.h>
#include <psb_types.h>

#include <utility>

#ifdef DEAL_II_WITH_PSBLAS

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  SparseMatrix::SparseMatrix()
    : psblas_sparse_matrix(nullptr)
    , psblas_descriptor(nullptr)
    , state(internal::State::Default)
  {}

  SparseMatrix::SparseMatrix(const SparsityPattern &psblas_sparsity_pattern,
                             const MPI_Comm         communicator)
  {
    Assert(psblas_descriptor.get() == nullptr,
           ExcMessage("PSBLAS matrix descriptor must not be initialized."));

    Assert(psblas_sparsity_pattern.psblas_descriptor.get() != nullptr,
           ExcMessage("The given SparsityPattern is not valid."));

    this->communicator = communicator;
    Assert(communicator != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to SparseMatrix::reinit()."));

    psblas_descriptor = psblas_sparsity_pattern.psblas_descriptor;

    // Create a new PSBLAS sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    // Initialize the sparse matrix with the descriptor
    int err =
      psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcAllocationPSBLASMatrix(err));
  }



  SparseMatrix::~SparseMatrix()
  {
    if (psblas_sparse_matrix != nullptr && psblas_descriptor.get() != nullptr)
      {
        // We clear the underlying PSBLAS sparse matrix
        int err = psb_c_dspfree(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcCallingPSBLASFunction(err, "psb_c_dspfree"));
      }
  }



  void
  SparseMatrix::copy_from(const SparseMatrix &other)
  {
    if (this == &other)
      return;

    DEAL_II_NOT_IMPLEMENTED();
  }


  void
  SparseMatrix::reinit(const IndexSet &index_set, const MPI_Comm comm)
  {
    Assert(index_set.n_elements() > 0,
           ExcMessage("An empty IndexSet has been given."));

    Assert(comm != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to SparseMatrix::reinit()."));
    communicator = comm;

    // Free old resources before reinitializing
    if (psblas_sparse_matrix != nullptr && psblas_descriptor.get() != nullptr)
      {
        int err = psb_c_dspfree(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcCallingPSBLASFunction(err, "psb_c_dspfree"));
        psblas_sparse_matrix = nullptr;
      }

    // Create a new PSBLAS descriptor
    psblas_descriptor.reset(psb_c_new_descriptor(),
                            PSCToolkitWrappers::internal::DescriptorDeleter());

    // Use get_index_vector() from IndexSet to get the indexes
    const std::vector<types::global_dof_index> &indexes =
      index_set.get_index_vector();

    psb_i_t number_of_local_indexes = indexes.size(); // Number of local indexes

    // Copy the indexes into an array called vl
    std::vector<psb_l_t> vl(number_of_local_indexes);
    for (psb_i_t i = 0; i < number_of_local_indexes; ++i)
      {
        const auto psblas_index = static_cast<psb_l_t>(indexes[i]);
        AssertIntegerConversion(psblas_index, indexes[i]);
        vl[i] = psblas_index;
      }

    // Insert the indexes into the descriptor
    psb_c_cdall_vl(number_of_local_indexes,
                   vl.data(),
                   *InitFinalize::get_psblas_context(),
                   psblas_descriptor.get());

    // Create a new PSBLAS sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    // Initialize the sparse matrix with the descriptor
    int err =
      psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcAllocationPSBLASMatrix(err));

    state = internal::State::Build;
  }



  void
  SparseMatrix::reinit(const SparsityPattern &psblas_sparsity_pattern,
                       const MPI_Comm         communicator)
  {
    Assert(psblas_sparsity_pattern.psblas_descriptor.get() != nullptr,
           ExcMessage("The given SparsityPattern is not valid."));

    this->communicator = communicator;
    Assert(communicator != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to SparseMatrix::reinit()."));

    // Free old resources before reinitializing
    if (psblas_sparse_matrix != nullptr && psblas_descriptor.get() != nullptr)
      {
        int err = psb_c_dspfree(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcCallingPSBLASFunction(err, "psb_c_dspfree"));
        psblas_sparse_matrix = nullptr;
      }

    psblas_descriptor = psblas_sparsity_pattern.psblas_descriptor;

    // Create a new PSBLAS sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    // Initialize the sparse matrix with the descriptor
    int err =
      psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcAllocationPSBLASMatrix(err));

    state = internal::State::Build;
  }



  void
  SparseMatrix::reinit(const IndexSet               &local_rows,
                       const DynamicSparsityPattern &sparsity_pattern,
                       const MPI_Comm                communicator)
  {
    Assert(sparsity_pattern.n_rows() == sparsity_pattern.n_cols(),
           ExcNotQuadratic());

    // Check dimensions match IndexSet
    Assert(sparsity_pattern.n_rows() == local_rows.size(),
           ExcMessage(
             "SparsityPattern and IndexSet have different number of rows"));

    Assert(local_rows.is_ascending_and_one_to_one(communicator),
           ExcNotImplemented());

    // Free old resources before reinitializing
    if (psblas_sparse_matrix != nullptr && psblas_descriptor.get() != nullptr)
      {
        int err = psb_c_dspfree(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcCallingPSBLASFunction(err, "psb_c_dspfree"));
        psblas_sparse_matrix = nullptr;
      }


    if constexpr (running_in_debug_mode())
      {
        types::global_dof_index row_owners =
          Utilities::MPI::sum(local_rows.n_elements(), communicator);
        Assert(row_owners == sparsity_pattern.n_rows(),
               ExcMessage(
                 std::string(
                   "Each row has to be owned by exactly one owner (n_rows()=") +
                 std::to_string(sparsity_pattern.n_rows()) +
                 " but sum(local_rows.n_elements())=" +
                 std::to_string(row_owners) + ")"));
      }

    this->communicator = communicator;

    // Set up the PSBLAS descriptor from the local IndexSet
    psblas_descriptor.reset(psb_c_new_descriptor(),
                            PSCToolkitWrappers::internal::DescriptorDeleter());

    {
      const std::vector<types::global_dof_index> indexes =
        local_rows.get_index_vector();
      const psb_i_t n_local = static_cast<psb_i_t>(indexes.size());

      std::vector<psb_l_t> vl(n_local);
      for (psb_i_t i = 0; i < n_local; ++i)
        {
          const auto idx = static_cast<psb_l_t>(indexes[i]);
          AssertIntegerConversion(idx, indexes[i]);
          vl[i] = idx;
        }

      psb_c_cdall_vl(n_local,
                     vl.data(),
                     *InitFinalize::get_psblas_context(),
                     psblas_descriptor.get());
    }

    // Create and initialize the sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    int err;

    if (local_rows.n_elements() > 0)
      {
        const psb_l_t local_row_start =
          static_cast<psb_l_t>(local_rows.nth_index_in_set(0));
        const psb_l_t local_row_end =
          local_row_start + static_cast<psb_l_t>(local_rows.n_elements());

        // Insert entries from the sparsity pattern row by row
        for (psb_l_t i = local_row_start; i < local_row_end; ++i)
          {
            const auto row_length =
              static_cast<psb_i_t>(sparsity_pattern.row_length(i));
            if (row_length == 0)
              continue;

            std::vector<psb_l_t> irw(row_length, i);
            std::vector<psb_l_t> icl(row_length);

            psb_i_t k = 0;
            for (typename DynamicSparsityPattern::iterator p =
                   sparsity_pattern.begin(i);
                 p != sparsity_pattern.end(i);
                 ++p, ++k)
              {
                const auto col = static_cast<psb_l_t>(p->column());
                AssertIntegerConversion(col, p->column());
                icl[k] = col;
              }

            err = psb_c_cdins(row_length,
                              irw.data(),
                              icl.data(),
                              psblas_descriptor.get());
            Assert(err == 0, ExcInsertionInPSBLASMatrix(err));
          }

        // Initialize the sparse matrix with the descriptor
        err =
          psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcAllocationPSBLASMatrix(err));
      }
    else
      {
        // empty local partition: just initialize
        err =
          psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcAllocationPSBLASMatrix(err));
      }

    state = internal::State::Build;
  }



  SparseMatrix::size_type
  SparseMatrix::local_size() const
  {
    return psb_c_cd_get_local_rows(psblas_descriptor.get());
  }


  std::pair<SparseMatrix::size_type, SparseMatrix::size_type>
  SparseMatrix::local_range() const
  {
    std::vector<psb_l_t> local_indices(local_size());
    int                  err = psb_c_cd_get_global_indices(local_indices.data(),
                                          local_size(),
                                          true /*only owned indices*/,
                                          psblas_descriptor.get());
    Assert(err == 0,
           ExcCallingPSBLASFunction(err, "psb_c_cd_get_global_indices"));

    return {local_indices[0], local_indices.back() + 1};
  }



  bool
  SparseMatrix::in_local_range(const size_type index) const
  {
    std::pair<size_type, size_type> local_range = this->local_range();
    return index >= local_range.first && index < local_range.second;
  }



  SparseMatrix::size_type
  SparseMatrix::m() const
  {
    return psb_c_cd_get_global_rows(psblas_descriptor.get());
  }



  SparseMatrix::size_type
  SparseMatrix::n() const
  {
    return psb_c_cd_get_global_cols(psblas_descriptor.get());
  }



  SparseMatrix::size_type
  SparseMatrix::n_nonzero_elements() const
  {
    return psb_c_dnnz(psblas_sparse_matrix, psblas_descriptor.get());
  }



  void
  SparseMatrix::set(const SparseMatrix::size_type  i,
                    const SparseMatrix::size_type  j,
                    const SparseMatrix::value_type value)
  {
    // Insert a value into the sparse matrix
    psb_l_t irw = i;
    psb_l_t icl = j;
    psb_d_t val = value;

    int err = psb_c_dspins(
      1, &irw, &icl, &val, psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcInsertionInPSBLASMatrix(err));
  }



  SparseMatrix::value_type
  SparseMatrix::el(const SparseMatrix::size_type i,
                   const SparseMatrix::size_type j) const
  {
    const auto psblas_index_i = static_cast<psb_l_t>(i);
    AssertIntegerConversion(psblas_index_i, i);
    const auto psblas_index_j = static_cast<psb_l_t>(j);
    AssertIntegerConversion(psblas_index_j, j);

    return psb_c_dmatgetelem(psblas_sparse_matrix,
                             psblas_index_i,
                             psblas_index_j,
                             psblas_descriptor.get());
  }



  SparseMatrix::value_type
  SparseMatrix::diag_element(const SparseMatrix::size_type i) const
  {
    Assert(m() == n(), ExcNotQuadratic());
    return el(i, i);
  }



  SparseMatrix::value_type
  SparseMatrix::operator()(const SparseMatrix::size_type i,
                           const SparseMatrix::size_type j) const
  {
    Assert(m() == n(), ExcNotQuadratic());
    return el(i, j);
  }



  psb_c_dspmat *
  SparseMatrix::get_psblas_matrix() const
  {
    return psblas_sparse_matrix;
  }



  psb_c_descriptor *
  SparseMatrix::get_psblas_descriptor() const
  {
    return psblas_descriptor.get();
  }



  MPI_Comm
  SparseMatrix::get_mpi_communicator() const
  {
    return communicator;
  }



  void
  SparseMatrix::set(const std::vector<SparseMatrix::size_type> &indices,
                    const FullMatrix<double>                   &matrix)
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(matrix.m() == indices.size(),
           ExcDimensionMismatch(matrix.m(), indices.size()));
    Assert(matrix.n() == indices.size(),
           ExcDimensionMismatch(matrix.n(), indices.size()));

    // Get the number of indices.
    const unsigned int n_indices = indices.size();
    psb_i_t            nz = n_indices * n_indices; // Number of non-zero entries

    // Fill the arrays with row indices, column indices, and values
    std::vector<psb_l_t> irw(nz);
    std::vector<psb_l_t> icl(nz);
    std::vector<psb_d_t> val(nz);
    for (unsigned int i = 0; i < n_indices; ++i)
      {
        for (unsigned int j = 0; j < n_indices; ++j)
          {
            const auto psblas_row_index = static_cast<psb_l_t>(indices[i]);
            const auto psblas_col_index = static_cast<psb_l_t>(indices[j]);
            AssertIntegerConversion(psblas_row_index, indices[i]);
            AssertIntegerConversion(psblas_col_index, indices[j]);
            irw[i * n_indices + j] = psblas_row_index;
            icl[i * n_indices + j] = psblas_col_index;
            val[i * n_indices + j] = matrix(i, j);
          }
      }

    // Insert the values into the sparse matrix
    int err = psb_c_dspins(nz,
                           irw.data(),
                           icl.data(),
                           val.data(),
                           psblas_sparse_matrix,
                           psblas_descriptor.get());

    Assert(err == 0, ExcInsertionInPSBLASMatrix(err));
  }



  void
  SparseMatrix::add(const SparseMatrix::size_type  i,
                    const SparseMatrix::size_type  j,
                    const SparseMatrix::value_type value)
  {
    AssertIsFinite(value);

    psb_l_t irw  = i;
    psb_l_t icl  = j;
    int     info = psb_c_dspins(1 /*nz*/,
                            &irw,
                            &icl,
                            &value,
                            psblas_sparse_matrix,
                            psblas_descriptor.get());
    Assert(info == 0, ExcInsertionInPSBLASMatrix(info));
  }



  void
  SparseMatrix::add(const SparseMatrix::size_type               row,
                    const SparseMatrix::size_type               ncols,
                    const std::vector<SparseMatrix::size_type> &col_indices,
                    const SparseMatrix::value_type             *values,
                    const bool,
                    const bool)
  {
    Assert(col_indices.size() == ncols,
           ExcDimensionMismatch(col_indices.size(), ncols));
    // Call the function below taking a raw pointer from the vector
    add(row, ncols, col_indices.data(), values, false, false);
  }



  void
  SparseMatrix::add(const SparseMatrix::size_type                row,
                    const SparseMatrix::size_type                ncols,
                    const std::vector<SparseMatrix::size_type>  &col_indices,
                    const std::vector<SparseMatrix::value_type> &values,
                    const bool,
                    const bool)
  {
    Assert(col_indices.size() == ncols,
           ExcDimensionMismatch(col_indices.size(), ncols));
    Assert(values.size() == ncols, ExcDimensionMismatch(ncols, values.size()));
    // Call the function below taking a raw pointer from the vector
    add(row, ncols, col_indices.data(), values.data(), false, false);
  }



  void
  SparseMatrix::add(const SparseMatrix::size_type   row,
                    const SparseMatrix::size_type   ncols,
                    const SparseMatrix::size_type  *col_indices,
                    const SparseMatrix::value_type *values,
                    const bool,
                    const bool)
  {
    std::vector<psb_l_t> irw(ncols);
    std::vector<psb_l_t> icl(ncols);
    for (SparseMatrix::size_type i = 0; i < ncols; ++i)
      {
        irw[i] = static_cast<psb_l_t>(row);
        icl[i] = static_cast<psb_l_t>(col_indices[i]);
      }

    int info = psb_c_dspins(ncols /*nz*/,
                            irw.data(),
                            icl.data(),
                            values,
                            psblas_sparse_matrix,
                            psblas_descriptor.get());
    Assert(info == 0, ExcInsertionInPSBLASMatrix(info));
  }



  void
  SparseMatrix::compress()
  {
    // We start by checking if the descriptor has already been assembled
    // elsewhere
    int err = -1;
    if (!psb_c_cd_is_asb(psblas_descriptor.get()))
      {
        err = psb_c_cdasb(psblas_descriptor.get());
        Assert(err == 0, ExcAssemblePSBLASDescriptor(err));
      }

    // Check if the sparse matrix is not already assembled
    if (!psb_c_dis_matasb(psblas_sparse_matrix, psblas_descriptor.get()))
      {
        // ... and the sparse matrix
        err = psb_c_dspasb(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0, ExcAssemblePSBLASMatrix(err));
        state = internal::State::Assembled;
      }
  }



  void
  SparseMatrix::vmult(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(state == internal::State::Assembled,
           ExcMessage("PSBLAS matrix has not been assembled."));

    int err = psb_c_dspmm(value_type(1.0), // alpha
                          psblas_sparse_matrix,
                          src.psblas_vector,
                          value_type(0.0), // beta
                          dst.psblas_vector,
                          psblas_descriptor.get());
    Assert(err == 0, ExcMatVecPSBLAS(err));
  }



  void
  SparseMatrix::vmult_add(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(state == internal::State::Assembled,
           ExcMessage("PSBLAS matrix has not been assembled."));

    int err = psb_c_dspmm(value_type(1.0), // alpha
                          psblas_sparse_matrix,
                          src.psblas_vector,
                          value_type(1.0), // beta
                          dst.psblas_vector,
                          psblas_descriptor.get());
    Assert(err == 0, ExcMatVecPSBLAS(err));
  }



  void
  SparseMatrix::Tvmult(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(state == internal::State::Assembled,
           ExcMessage("PSBLAS matrix has not been assembled."));

    char transpose = 'T';
    int  err       = psb_c_dspmm_opt(value_type(1.0), // alpha
                              psblas_sparse_matrix,
                              src.psblas_vector,
                              value_type(0.0), // beta
                              dst.psblas_vector,
                              psblas_descriptor.get(),
                              &transpose,
                              true);
    Assert(err == 0, ExcMatVecPSBLAS(err));
  }



  void
  SparseMatrix::Tvmult_add(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(&src != &dst, ExcSourceEqualsDestination());
    Assert(state == internal::State::Assembled,
           ExcMessage("PSBLAS matrix has not been assembled."));

    char transpose = 'T';
    int  err       = psb_c_dspmm_opt(value_type(1.0), // alpha
                              psblas_sparse_matrix,
                              src.psblas_vector,
                              value_type(1.0), // beta
                              dst.psblas_vector,
                              psblas_descriptor.get(),
                              &transpose,
                              true);
    Assert(err == 0, ExcMatVecPSBLAS(err));
  }



  SparseMatrix::value_type
  SparseMatrix::l1_norm() const
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  SparseMatrix::value_type
  SparseMatrix::linfty_norm() const
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  SparseMatrix::value_type
  SparseMatrix::frobenius_norm() const
  {
    DEAL_II_NOT_IMPLEMENTED();
  }



  SparseMatrix::value_type
  SparseMatrix::trace() const
  {
    Assert(m() == n(), ExcNotQuadratic());
    value_type local_diagonal_sum = 0;
    const auto local_indices      = local_range();
    for (types::global_dof_index idx = local_indices.first;
         idx < local_indices.second;
         ++idx)
      local_diagonal_sum += diag_element(idx);

    return Utilities::MPI::sum(local_diagonal_sum, communicator);
  }



} // namespace PSCToolkitWrappers

DEAL_II_NAMESPACE_CLOSE
#endif // DEAL_II_WITH_PSBLAS
