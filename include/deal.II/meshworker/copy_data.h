// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2020 by the deal.II authors
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

#ifndef dealii_meshworker_copy_data_h
#define dealii_meshworker_copy_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/types.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

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
   *
   * @author Luca Heltai, 2019.
   */
  template <int n_matrices    = 1,
            int n_vectors     = n_matrices,
            int n_dof_indices = n_matrices>
  struct CopyData
  {
    /**
     * Initialize everything with the same @p size. This is usually the number
     * of local degrees of freedom.
     */
    explicit CopyData(const unsigned int size);

    /**
     * For every object, specify the size they should have.
     */
    explicit CopyData(
      const std::array<std::array<unsigned int, 2>, n_matrices> &matrix_sizes,
      const std::array<unsigned int, n_vectors> &                vector_sizes,
      const std::array<unsigned int, n_dof_indices> &dof_indices_sizes);

    /**
     * Deep copy constructor.
     */
    CopyData(const CopyData<n_matrices, n_vectors, n_dof_indices> &other) =
      default;

    /**
     * An array of local matrices.
     */
    std::array<FullMatrix<double>, n_matrices> matrices;

    /**
     * An array of local vectors.
     */
    std::array<Vector<double>, n_vectors> vectors;

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
  template <int n_matrices, int n_vectors, int n_dof_indices>
  CopyData<n_matrices, n_vectors, n_dof_indices>::CopyData(
    const unsigned int size)
  {
    for (auto &m : matrices)
      m.reinit({size, size});
    for (auto &v : vectors)
      v.reinit(size);
    for (auto &d : local_dof_indices)
      d.resize(size);
  }



  template <int n_matrices, int n_vectors, int n_dof_indices>
  CopyData<n_matrices, n_vectors, n_dof_indices>::CopyData(
    const std::array<std::array<unsigned int, 2>, n_matrices> &matrix_sizes,
    const std::array<unsigned int, n_vectors> &                vector_sizes,
    const std::array<unsigned int, n_dof_indices> &dof_indices_sizes)
  {
    for (unsigned int i = 0; i < n_matrices; ++i)
      matrices[i].reinit(matrix_sizes[i++]);

    for (unsigned int i = 0; i < n_vectors; ++i)
      vectors[i].reinit(vector_sizes[i++]);

    for (unsigned int i = 0; i < n_dof_indices; ++i)
      local_dof_indices[i].resize(dof_indices_sizes[i++]);
  }

#endif // DOXYGEN
} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif
