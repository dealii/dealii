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

#include <deal.II/lac/psblas_sparse_matrix.h>

#include <psb_types.h>

#ifdef DEAL_II_WITH_PSBLAS

DEAL_II_NAMESPACE_OPEN

namespace PSCToolkitWrappers
{

  SparseMatrix::SparseMatrix()
    : psblas_sparse_matrix(nullptr)
    , psblas_descriptor(nullptr)
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
    Assert(err == 0, ExcMessage("Error initializing PSBLAS sparse matrix."));
  }



  SparseMatrix::~SparseMatrix()
  {
    Assert((psblas_sparse_matrix != nullptr &&
            psblas_descriptor.get() != nullptr),
           ExcMessage("PSBLAS sparse matrix or descriptor is null."));

    // We first clear the sparse matrix
    int err = psb_c_dspfree(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcMessage("Error freeing PSBLAS sparse matrix."));

    // ... and then the descriptor
    err = psb_c_cdfree(psblas_descriptor.get());
    Assert(err == 0, ExcMessage("Error freeing PSBLAS descriptor."));
  }



  void
  SparseMatrix::reinit(const IndexSet &index_set, const MPI_Comm comm)
  {
    Assert(index_set.n_elements() > 0,
           ExcMessage("An empty IndexSet has been given."));

    Assert(comm != MPI_COMM_NULL,
           ExcMessage("MPI_COMM_NULL passed to SparseMatrix::reinit()."));
    communicator = comm;

    // Create a new PSBLAS descriptor
    psblas_descriptor.reset(psb_c_new_descriptor());

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
    psblas_context = InitFinalize::get_psblas_context();
    psb_c_cdall_vl(number_of_local_indexes,
                   vl.data(),
                   *psblas_context,
                   psblas_descriptor.get());

    // Create a new PSBLAS sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    // Initialize the sparse matrix with the descriptor
    int err =
      psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcMessage("Error initializing PSBLAS sparse matrix."));
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

    psblas_descriptor = psblas_sparsity_pattern.psblas_descriptor;

    // Create a new PSBLAS sparse matrix
    psblas_sparse_matrix = psb_c_new_dspmat();

    // Initialize the sparse matrix with the descriptor
    int err =
      psb_c_dspall_remote(psblas_sparse_matrix, psblas_descriptor.get());
    Assert(err == 0, ExcMessage("Error initializing PSBLAS sparse matrix."));
  }



  SparseMatrix::size_type
  SparseMatrix::local_size() const
  {
    return psb_c_cd_get_local_rows(psblas_descriptor.get());
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
    Assert(err == 0, ExcMessage("Failed insertion into PSBLAS sparse matrix."));
  }



  SparseMatrix::value_type
  SparseMatrix::el(const SparseMatrix::size_type,
                   const SparseMatrix::size_type) const
  {
    // Not implemented yet from PSBLAS.
    AssertThrow(false, ExcNotImplemented());
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

    Assert(err == 0, ExcMessage("Failed insertion into PSBLAS sparse matrix."));
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
    Assert(info == 0,
           ExcMessage("Error inserting values into PSBLAS sparse matrix."));
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
    ExcMessage("Error inserting values into PSBLAS sparse matrix.");
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
    ExcMessage("Error inserting values into PSBLAS sparse matrix.");
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
    Assert(info == 0,
           ExcMessage("Error inserting values into PSBLAS sparse matrix."));
  }



  void
  SparseMatrix::compress()
  {
    // We start by checking if the vector has already been assembled elsewhere
    int err = -1;
    if (!psb_c_cd_is_asb(psblas_descriptor.get()))
      {
        err = psb_c_cdasb(psblas_descriptor.get());
        Assert(err == 0, ExcMessage("Error while finalizing descriptor."));
      }

    // Check if the sparse matrix is not already assembled
    if (!psb_c_dis_matasb(psblas_sparse_matrix, psblas_descriptor.get()))
      {
        // ... and the sparse matrix
        err = psb_c_dspasb(psblas_sparse_matrix, psblas_descriptor.get());
        Assert(err == 0,
               ExcMessage("Error while assembling the PSBLAS sparse matrix."));
      }
  }



  void
  SparseMatrix::vmult(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(dst.psblas_vector != nullptr,
           ExcMessage("Destination PSBLAS vector has not been initialized."));

    int err = psb_c_dspmm(1.0, // alpha
                          psblas_sparse_matrix,
                          src.psblas_vector,
                          0.0, // beta
                          dst.psblas_vector,
                          psblas_descriptor.get());
    Assert(err == 0,
           ExcMessage("Error in PSBLAS matrix-vector multiplication."));
  }



  void
  SparseMatrix::Tvmult(Vector &dst, const Vector &src) const
  {
    Assert(psblas_sparse_matrix != nullptr,
           ExcMessage("PSBLAS matrix has not been initialized."));
    Assert(src.psblas_vector != nullptr,
           ExcMessage("Source PSBLAS vector has not been initialized."));
    Assert(dst.psblas_vector != nullptr,
           ExcMessage("Destination PSBLAS vector has not been initialized."));

    char transpose = 'T';
    int  err       = psb_c_dspmm_opt(1.0, // alpha
                              psblas_sparse_matrix,
                              src.psblas_vector,
                              0.0, // beta
                              dst.psblas_vector,
                              psblas_descriptor.get(),
                              &transpose,
                              true);
    Assert(err == 0,
           ExcMessage(
             "Error in PSBLAS transposed matrix-vector multiplication."));
  }

} // namespace PSCToolkitWrappers

DEAL_II_NAMESPACE_CLOSE
#endif // DEAL_II_WITH_PSBLAS