// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_meshworker_copy_data_h
#define dealii_meshworker_copy_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/ndarray.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <algorithm>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * Helper copy data struct.
   *
   * This class is a good default drop in CopyData object for the
   * WorkStream::run() and MeshWorker::mesh_loop() functions.
   *
   * It arrays of (local) full matrices, vectors, and local degrees of freedom
   * index vectors, with size determined by the corresponding template argument.
   *
   * In particular, you can specify the following template arguments
   *
   * - @tparam n_matrices: Size of the array of matrices
   * - @tparam n_vectors: size of the array of vectors
   * - @tparam n_dof_indices: size of the array of local dof indices
   * - @tparam ScalarType: The data type stored by the vectors and matrices.
   *           This must be a real or complex floating point number.
   */
  template <int n_matrices      = 1,
            int n_vectors       = n_matrices,
            int n_dof_indices   = n_matrices,
            typename ScalarType = double>
  struct CopyData
  {
    /**
     * Default constructor. All of the @p matrices, @p vectors and
     * @p local_dof_indices are empty, and should be initialized using
     * one of the reinit() functions.
     */
    CopyData() = default;

    /**
     * Initialize everything with the same @p size. This is usually the number
     * of local degrees of freedom.
     */
    explicit CopyData(const unsigned int size);

    /**
     * For every object, specify the size they should have.
     */
    explicit CopyData(
      const ndarray<unsigned int, n_matrices, 2>    &matrix_sizes,
      const std::array<unsigned int, n_vectors>     &vector_sizes,
      const std::array<unsigned int, n_dof_indices> &dof_indices_sizes);

    /**
     * Copy constructor.
     */
    CopyData(const CopyData &other) = default;

    /**
     * Copy operator.
     */
    CopyData &
    operator=(const CopyData &other) = default;

    /**
     * Reinitialize everything the same @p size. This is usually the number of
     * local degrees of freedom.
     */
    void
    reinit(const unsigned int size);

    /**
     * Reinitialize the @p index'th matrix, vector and local DoF index vector
     * with the same @p size. This is usually the number of local degrees of
     * freedom.
     *
     * @note In the event that there are a different number of @p matrices,
     * @p vectors and @p local_dof_indices, this function will not throw an
     * error when the index is over the size of one or two of these data
     * structures. However, if the index exceeds the size of all three of them
     * then an error will be thrown.
     */
    void
    reinit(const unsigned int index, const unsigned int size);

    /**
     * Reinitialize the @p index'th matrix with @p size_rows `x` @p size_columns,
     * and the vector and local DoF index vector with @p size_columns. These
     * sizes usually correspond to some local degrees of freedom (e.g., the
     * number of cell DoFs for @p size_rows and number of neighbor cell DoFs
     * for @p size_columns).
     *
     * @note In the event that there are a different number of @p matrices,
     * @p vectors and @p local_dof_indices, this function will not throw an
     * error when the index is over the size of one or two of these data
     * structures. However, if the index exceeds the size of all three of them
     * then an error will be thrown.
     */
    void
    reinit(const unsigned int index,
           const unsigned int size_rows,
           const unsigned int size_columns);

    /**
     * An array of local matrices.
     */
    std::array<FullMatrix<ScalarType>, n_matrices> matrices;

    /**
     * An array of local vectors.
     */
    std::array<Vector<ScalarType>, n_vectors> vectors;

    /**
     * An array of local degrees of freedom indices.
     */
    std::array<std::vector<types::global_dof_index>, n_dof_indices>
      local_dof_indices;
  };


#ifndef DOXYGEN
  //
  // Template definitions
  //
  template <int n_matrices,
            int n_vectors,
            int n_dof_indices,
            typename ScalarType>
  CopyData<n_matrices, n_vectors, n_dof_indices, ScalarType>::CopyData(
    const unsigned int size)
  {
    reinit(size);
  }



  template <int n_matrices,
            int n_vectors,
            int n_dof_indices,
            typename ScalarType>
  CopyData<n_matrices, n_vectors, n_dof_indices, ScalarType>::CopyData(
    const ndarray<unsigned int, n_matrices, 2>    &matrix_sizes,
    const std::array<unsigned int, n_vectors>     &vector_sizes,
    const std::array<unsigned int, n_dof_indices> &dof_indices_sizes)
  {
    for (unsigned int i = 0; i < n_matrices; ++i)
      matrices[i].reinit(matrix_sizes[i++]);

    for (unsigned int i = 0; i < n_vectors; ++i)
      vectors[i].reinit(vector_sizes[i++]);

    for (unsigned int i = 0; i < n_dof_indices; ++i)
      local_dof_indices[i].resize(dof_indices_sizes[i++]);
  }



  template <int n_matrices,
            int n_vectors,
            int n_dof_indices,
            typename ScalarType>
  void
  CopyData<n_matrices, n_vectors, n_dof_indices, ScalarType>::reinit(
    const unsigned int size)
  {
    for (auto &m : matrices)
      m.reinit({size, size});
    for (auto &v : vectors)
      v.reinit(size);
    for (auto &d : local_dof_indices)
      d.resize(size);
  }



  template <int n_matrices,
            int n_vectors,
            int n_dof_indices,
            typename ScalarType>
  void
  CopyData<n_matrices, n_vectors, n_dof_indices, ScalarType>::reinit(
    const unsigned int index,
    const unsigned int size)
  {
    reinit(index, size, size);
  }



  template <int n_matrices,
            int n_vectors,
            int n_dof_indices,
            typename ScalarType>
  void
  CopyData<n_matrices, n_vectors, n_dof_indices, ScalarType>::reinit(
    const unsigned int index,
    const unsigned int size_rows,
    const unsigned int size_columns)
  {
    // We permit different numbers of matrices, vectors and DoF index vectors.
    // So we have to be a bit permissive here.
    constexpr int max_index = std::max({n_matrices, n_vectors, n_dof_indices});
    Assert(index < max_index, ExcIndexRange(index, 0, max_index));

    if (index < n_matrices)
      matrices[index].reinit({size_rows, size_columns});

    if (index < n_vectors)
      vectors[index].reinit(size_columns);

    if (index < n_dof_indices)
      local_dof_indices[index].resize(size_columns);
  }

#endif // DOXYGEN
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif
