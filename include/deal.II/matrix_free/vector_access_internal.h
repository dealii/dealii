// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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


#ifndef dealii_matrix_free_vector_access_internal_h
#define dealii_matrix_free_vector_access_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/type_traits.h>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  // below we use type-traits from matrix-free/type_traits.h

  // access to generic const vectors that have operator ().
  // FIXME: this is wrong for Trilinos/Petsc MPI vectors
  // where we should first do Partitioner::local_to_global()
  template <typename VectorType,
            typename std::enable_if<!has_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline typename VectorType::value_type
  vector_access(const VectorType &vec, const unsigned int entry)
  {
    return vec(entry);
  }



  // access to generic non-const vectors that have operator ().
  // FIXME: this is wrong for Trilinos/Petsc MPI vectors
  // where we should first do Partitioner::local_to_global()
  template <typename VectorType,
            typename std::enable_if<!has_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline typename VectorType::value_type &
  vector_access(VectorType &vec, const unsigned int entry)
  {
    return vec(entry);
  }



  // access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <typename VectorType,
            typename std::enable_if<has_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline typename VectorType::value_type &
  vector_access(VectorType &vec, const unsigned int entry)
  {
    return vec.local_element(entry);
  }



  // same for const access
  template <typename VectorType,
            typename std::enable_if<has_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline typename VectorType::value_type
  vector_access(const VectorType &vec, const unsigned int entry)
  {
    return vec.local_element(entry);
  }



  template <typename VectorType,
            typename std::enable_if<has_add_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_add(VectorType &                           vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vec.add_local_element(entry, val);
  }



  template <typename VectorType,
            typename std::enable_if<!has_add_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_add(VectorType &                           vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vector_access(vec, entry) += val;
  }



  template <typename VectorType,
            typename std::enable_if<has_add_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_add_global(VectorType &                           vec,
                           const types::global_dof_index          entry,
                           const typename VectorType::value_type &val)
  {
    vec.add(entry, val);
  }



  template <typename VectorType,
            typename std::enable_if<!has_add_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_add_global(VectorType &                           vec,
                           const types::global_dof_index          entry,
                           const typename VectorType::value_type &val)
  {
    vec(entry) += val;
  }



  template <typename VectorType,
            typename std::enable_if<has_set_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_set(VectorType &                           vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vec.set_local_element(entry, val);
  }



  template <typename VectorType,
            typename std::enable_if<!has_set_local_element<VectorType>::value,
                                    VectorType>::type * = nullptr>
  inline void
  vector_access_set(VectorType &                           vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vector_access(vec, entry) = val;
  }



  // this is to make sure that the parallel partitioning in VectorType
  // is really the same as stored in MatrixFree.
  // version below is when has_partitioners_are_compatible == false
  // FIXME: this is incorrect for PETSc/Trilinos MPI vectors
  template <
    typename VectorType,
    typename std::enable_if<!has_partitioners_are_compatible<VectorType>::value,
                            VectorType>::type * = nullptr>
  inline void
  check_vector_compatibility(
    const VectorType &                            vec,
    const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    (void)vec;
    (void)dof_info;

    AssertDimension(vec.size(), dof_info.vector_partitioner->size());
  }



  // same as above for has_partitioners_are_compatible == true
  template <
    typename VectorType,
    typename std::enable_if<has_partitioners_are_compatible<VectorType>::value,
                            VectorType>::type * = nullptr>
  inline void
  check_vector_compatibility(
    const VectorType &                            vec,
    const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    (void)vec;
    (void)dof_info;
    Assert(vec.partitioners_are_compatible(*dof_info.vector_partitioner),
           ExcMessage(
             "The parallel layout of the given vector is not "
             "compatible with the parallel partitioning in MatrixFree. "
             "Use MatrixFree::initialize_dof_vector to get a "
             "compatible vector."));
  }



  // Below, three classes (VectorReader, VectorSetter,
  // VectorDistributorLocalToGlobal) implement the same interface and can be
  // used to to read from vector, set elements of a vector and add to elements
  // of the vector.

  // 1. A class to read data from vector
  template <typename Number, typename VectorizedArrayType>
  struct VectorReader
  {
    template <typename VectorType>
    void
    process_dof(const unsigned int index,
                const VectorType & vec,
                Number &           res) const
    {
      res = vector_access(vec, index);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType &         vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, true>) const
    {
      const Number *vec_ptr = vec.begin() + dof_index;
      for (unsigned int i = 0; i < dofs_per_cell;
           ++i, vec_ptr += VectorizedArrayType::size())
        dof_values[i].load(vec_ptr);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            const VectorType &   vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          dof_values[i][v] =
            vector_access(vec, dof_index + v + i * VectorizedArrayType::size());
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      VectorType &         vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, true>) const
    {
      dealii::vectorized_load_and_transpose(dofs_per_cell,
                                            vec.begin(),
                                            dof_indices,
                                            dof_values);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      const VectorType &   vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, false>) const
    {
      for (unsigned int d = 0; d < dofs_per_cell; ++d)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          dof_values[d][v] = vector_access(vec, dof_indices[v] + d);
    }



    // variant where VectorType::value_type is the same as Number -> can call
    // gather
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       VectorType &         vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, true>) const
    {
      res.gather(vec.begin() + constant_offset, indices);
    }



    // variant where VectorType::value_type is not the same as Number -> must
    // manually load the data
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       const VectorType &   vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        res[v] = vector_access(vec, indices[v] + constant_offset);
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       const VectorType &            vec,
                       Number &                      res) const
    {
      res = vec(index);
    }



    void
    pre_constraints(const Number &, Number &res) const
    {
      res = Number();
    }



    template <typename VectorType>
    void
    process_constraint(const unsigned int index,
                       const Number       weight,
                       const VectorType & vec,
                       Number &           res) const
    {
      res += weight * vector_access(vec, index);
    }



    void
    post_constraints(const Number &sum, Number &write_pos) const
    {
      write_pos = sum;
    }



    void
    process_empty(VectorizedArrayType &res) const
    {
      res = VectorizedArrayType();
    }
  };



  // 2. A class to add values to the vector during
  // FEEvaluation::distribute_local_to_global() call
  template <typename Number, typename VectorizedArrayType>
  struct VectorDistributorLocalToGlobal
  {
    template <typename VectorType>
    void
    process_dof(const unsigned int index, VectorType &vec, Number &res) const
    {
      vector_access_add(vec, index, res);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType &         vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, true>) const
    {
      Number *vec_ptr = vec.begin() + dof_index;
      for (unsigned int i = 0; i < dofs_per_cell;
           ++i, vec_ptr += VectorizedArrayType::size())
        {
          VectorizedArrayType tmp;
          tmp.load(vec_ptr);
          tmp += dof_values[i];
          tmp.store(vec_ptr);
        }
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType &         vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access_add(vec,
                            dof_index + v + i * VectorizedArrayType::size(),
                            dof_values[i][v]);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      VectorType &         vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, true>) const
    {
      vectorized_transpose_and_store(
        true, dofs_per_cell, dof_values, dof_indices, vec.begin());
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      VectorType &         vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, false>) const
    {
      for (unsigned int d = 0; d < dofs_per_cell; ++d)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access_add(vec, dof_indices[v] + d, dof_values[d][v]);
    }



    // variant where VectorType::value_type is the same as Number -> can call
    // scatter
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       VectorType &         vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, true>) const
    {
#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS < 512
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access(vec, indices[v] + constant_offset) += res[v];
#else
      // only use gather in case there is also scatter.
      VectorizedArrayType tmp;
      tmp.gather(vec.begin() + constant_offset, indices);
      tmp += res;
      tmp.scatter(indices, vec.begin() + constant_offset);
#endif
    }



    // variant where VectorType::value_type is not the same as Number -> must
    // manually append all data
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       VectorType &         vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access_add(vec, indices[v] + constant_offset, res[v]);
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       VectorType &                  vec,
                       Number &                      res) const
    {
      vector_access_add_global(vec, index, res);
    }



    void
    pre_constraints(const Number &input, Number &res) const
    {
      res = input;
    }



    template <typename VectorType>
    void
    process_constraint(const unsigned int index,
                       const Number       weight,
                       VectorType &       vec,
                       Number &           res) const
    {
      vector_access_add(vec, index, weight * res);
    }



    void
    post_constraints(const Number &, Number &) const
    {}



    void
    process_empty(VectorizedArrayType &) const
    {}
  };



  // 3. A class to set elements of the vector
  template <typename Number, typename VectorizedArrayType>
  struct VectorSetter
  {
    template <typename VectorType>
    void
    process_dof(const unsigned int index, VectorType &vec, Number &res) const
    {
      vector_access(vec, index) = res;
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType &         vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, true>) const
    {
      Number *vec_ptr = vec.begin() + dof_index;
      for (unsigned int i = 0; i < dofs_per_cell;
           ++i, vec_ptr += VectorizedArrayType::size())
        dof_values[i].store(vec_ptr);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType &         vec,
                            VectorizedArrayType *dof_values,
                            std::integral_constant<bool, false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access(vec, dof_index + v + i * VectorizedArrayType::size()) =
            dof_values[i][v];
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      VectorType &         vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, true>) const
    {
      vectorized_transpose_and_store(
        false, dofs_per_cell, dof_values, dof_indices, vec.begin());
    }



    template <typename VectorType, bool booltype>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int * dof_indices,
                                      VectorType &         vec,
                                      VectorizedArrayType *dof_values,
                                      std::integral_constant<bool, false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access(vec, dof_indices[v] + i) = dof_values[i][v];
    }



    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       VectorType &         vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, true>) const
    {
      res.scatter(indices, vec.begin() + constant_offset);
    }



    template <typename VectorType>
    void
    process_dof_gather(const unsigned int * indices,
                       VectorType &         vec,
                       const unsigned int   constant_offset,
                       VectorizedArrayType &res,
                       std::integral_constant<bool, false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access(vec, indices[v] + constant_offset) = res[v];
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       VectorType &                  vec,
                       Number &                      res) const
    {
      vec(index) = res;
    }



    void
    pre_constraints(const Number &, Number &) const
    {}



    template <typename VectorType>
    void
    process_constraint(const unsigned int,
                       const Number,
                       VectorType &,
                       Number &) const
    {}



    void
    post_constraints(const Number &, Number &) const
    {}



    void
    process_empty(VectorizedArrayType &) const
    {}
  };
} // namespace internal


DEAL_II_NAMESPACE_CLOSE

#endif
