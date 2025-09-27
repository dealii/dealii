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


#ifndef dealii_matrix_free_vector_access_internal_h
#define dealii_matrix_free_vector_access_internal_h

#include <deal.II/base/config.h>

#include <deal.II/base/vectorization.h>

#include <deal.II/matrix_free/dof_info.h>
#include <deal.II/matrix_free/matrix_free.h>
#include <deal.II/matrix_free/type_traits.h>

#include <boost/algorithm/string/join.hpp>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  // below we use type-traits from matrix-free/type_traits.h



  // access to serial const vectors that have operator[].
  template <typename VectorType,
            std::enable_if_t<is_serial_vector_or_array<VectorType>::value,
                             VectorType> * = nullptr>
  inline typename VectorType::value_type
  vector_access(const VectorType &vec, const unsigned int entry)
  {
    return vec[entry];
  }



  // access to serial non-const vectors that have operator[].
  template <typename VectorType,
            std::enable_if_t<is_serial_vector_or_array<VectorType>::value,
                             VectorType> * = nullptr>
  inline typename VectorType::value_type &
  vector_access(VectorType &vec, const unsigned int entry)
  {
    return vec[entry];
  }



  // access to distributed MPI vectors that have a local_element(uint)
  // method to access data in local index space, which is what we use in
  // DoFInfo and hence in read_dof_values etc.
  template <
    typename VectorType,
    std::enable_if_t<has_local_element<VectorType>, VectorType> * = nullptr>
  inline typename VectorType::value_type &
  vector_access(VectorType &vec, const unsigned int entry)
  {
    return vec.local_element(entry);
  }



  // same for const access
  template <
    typename VectorType,
    std::enable_if_t<has_local_element<VectorType>, VectorType> * = nullptr>
  inline typename VectorType::value_type
  vector_access(const VectorType &vec, const unsigned int entry)
  {
    return vec.local_element(entry);
  }



  template <
    typename VectorType,
    std::enable_if_t<has_add_local_element<VectorType>, VectorType> * = nullptr>
  inline void
  vector_access_add(VectorType                            &vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vec.add_local_element(entry, val);
  }



  template <typename VectorType,
            std::enable_if_t<!has_add_local_element<VectorType>, VectorType> * =
              nullptr>
  inline void
  vector_access_add(VectorType                            &vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vector_access(vec, entry) += val;
  }



  template <
    typename VectorType,
    std::enable_if_t<has_add_local_element<VectorType>, VectorType> * = nullptr>
  inline void
  vector_access_add_global(VectorType                            &vec,
                           const types::global_dof_index          entry,
                           const typename VectorType::value_type &val)
  {
    vec.add(entry, val);
  }



  template <typename VectorType,
            std::enable_if_t<!has_add_local_element<VectorType>, VectorType> * =
              nullptr>
  inline void
  vector_access_add_global(VectorType                            &vec,
                           const types::global_dof_index          entry,
                           const typename VectorType::value_type &val)
  {
    vec[entry] += val;
  }



  template <
    typename VectorType,
    std::enable_if_t<has_set_local_element<VectorType>, VectorType> * = nullptr>
  inline void
  vector_access_set(VectorType                            &vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vec.set_local_element(entry, val);
  }



  template <typename VectorType,
            std::enable_if_t<!has_set_local_element<VectorType>, VectorType> * =
              nullptr>
  inline void
  vector_access_set(VectorType                            &vec,
                    const unsigned int                     entry,
                    const typename VectorType::value_type &val)
  {
    vector_access(vec, entry) = val;
  }



  // this is to make sure that the parallel partitioning in VectorType
  // is really the same as stored in MatrixFree.
  // version below is when has_partitioners_are_compatible == false
  // FIXME: this is incorrect for PETSc/Trilinos MPI vectors
  template <int dim,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            std::enable_if_t<!has_partitioners_are_compatible<VectorType>,
                             VectorType> * = nullptr>
  inline void
  check_vector_compatibility(
    const VectorType &vec,
    const MatrixFree<dim, Number, VectorizedArrayType> & /*matrix_free*/,
    const internal::MatrixFreeFunctions::DoFInfo &dof_info)
  {
    AssertDimension(vec.size(), dof_info.vector_partitioner->size());
  }



  // same as above for has_partitioners_are_compatible == true
  template <int dim,
            typename Number,
            typename VectorizedArrayType,
            typename VectorType,
            std::enable_if_t<has_partitioners_are_compatible<VectorType>,
                             VectorType> * = nullptr>
  inline void
  check_vector_compatibility(
    const VectorType                                   &vec,
    const MatrixFree<dim, Number, VectorizedArrayType> &matrix_free,
    const internal::MatrixFreeFunctions::DoFInfo       &dof_info)
  {
    (void)vec;
    (void)matrix_free;
    (void)dof_info;

    if constexpr (running_in_debug_mode())
      {
        if (vec.partitioners_are_compatible(*dof_info.vector_partitioner) ==
            false)
          {
            unsigned int dof_index = numbers::invalid_unsigned_int;

            for (unsigned int i = 0; i < matrix_free.n_components(); ++i)
              if (&matrix_free.get_dof_info(i) == &dof_info)
                {
                  dof_index = i;
                  break;
                }

            Assert(dof_index != numbers::invalid_unsigned_int,
                   ExcInternalError());

            std::vector<std::string> dof_indices_with_compatible_partitioners;

            for (unsigned int i = 0; i < matrix_free.n_components(); ++i)
              if (vec.partitioners_are_compatible(
                    *matrix_free.get_dof_info(i).vector_partitioner))
                dof_indices_with_compatible_partitioners.push_back(
                  std::to_string(i));

            if (dof_indices_with_compatible_partitioners.empty())
              {
                Assert(false,
                       ExcMessage(
                         "The parallel layout of the given vector is "
                         "compatible neither with the Partitioner of the "
                         "current FEEvaluation with dof_handler_index=" +
                         std::to_string(dof_index) +
                         " nor with any Partitioner in MatrixFree. A "
                         "potential reason is that you did not use "
                         "MatrixFree::initialize_dof_vector() to get a "
                         "compatible vector."));
              }
            else
              {
                Assert(
                  false,
                  ExcMessage(
                    "The parallel layout of the given vector is "
                    "not compatible with the Partitioner of the "
                    "current FEEvaluation with dof_handler_index=" +
                    std::to_string(dof_index) +
                    ". However, the underlying "
                    "MatrixFree contains Partitioner objects that are compatible. "
                    "They have the following dof_handler_index values: " +
                    boost::algorithm::join(
                      dof_indices_with_compatible_partitioners, ", ") +
                    ". Did you want to pass any of these values to the "
                    "constructor of the current FEEvaluation object or "
                    "did you not use MatrixFree::initialize_dof_vector() "
                    "with dof_handler_index=" +
                    std::to_string(dof_index) +
                    " to get a "
                    "compatible vector?"));
              }
          }
      }
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
                const VectorType  &vec,
                Number            &res) const
    {
      res = vector_access(vec, index);
    }



    template <typename VectorNumberType>
    void
    process_dof(const VectorNumberType &global, Number &local) const
    {
      local = global;
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType          &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<true>) const
    {
      if constexpr (running_in_debug_mode())
        {
          // in debug mode, run non-vectorized version because this path
          // has additional checks (e.g., regarding ghosting)
          process_dofs_vectorized(dofs_per_cell,
                                  dof_index,
                                  vec,
                                  dof_values,
                                  std::bool_constant<false>());
        }
      else
        {
          const Number *vec_ptr = vec.begin() + dof_index;
          for (unsigned int i = 0; i < dofs_per_cell;
               ++i, vec_ptr += VectorizedArrayType::size())
            dof_values[i].load(vec_ptr);
        }
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            const VectorType    &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          dof_values[i][v] =
            vector_access(vec, dof_index + v + i * VectorizedArrayType::size());
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int  *dof_indices,
                                      VectorType          &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<true>) const
    {
      dealii::vectorized_load_and_transpose(dofs_per_cell,
                                            vec.begin() + constant_offset,
                                            dof_indices,
                                            dof_values);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int  *dof_indices,
                                      const VectorType    &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<false>) const
    {
      for (unsigned int d = 0; d < dofs_per_cell; ++d)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          dof_values[d][v] =
            vector_access(vec, dof_indices[v] + constant_offset + d);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int       dofs_per_cell,
                                      const unsigned int      *dof_indices,
                                      VectorType              &vec,
                                      VectorizedArrayType     *dof_values,
                                      std::bool_constant<true> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int        dofs_per_cell,
                                      const unsigned int       *dof_indices,
                                      const VectorType         &vec,
                                      VectorizedArrayType      *dof_values,
                                      std::bool_constant<false> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int                                        dofs_per_cell,
      const std::array<Number2 *, VectorizedArrayType::size()> &global_ptr,
      VectorizedArrayType                                      *dof_values,
      std::bool_constant<true>) const
    {
      dealii::vectorized_load_and_transpose(dofs_per_cell,
                                            global_ptr,
                                            dof_values);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int,
      const std::array<Number2 *, VectorizedArrayType::size()> &,
      VectorizedArrayType *,
      std::bool_constant<false>) const
    {
      DEAL_II_NOT_IMPLEMENTED();
    }



    // variant where VectorType::value_type is the same as Number -> can call
    // gather
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int              *indices,
                       VectorType                      &vec,
                       const unsigned int               constant_offset,
                       typename VectorType::value_type *vec_ptr,
                       VectorizedArrayType             &res,
                       std::bool_constant<true>) const
    {
      (void)constant_offset;
      (void)vec;

      if constexpr (running_in_debug_mode())
        {
          // in debug mode, run non-vectorized version because this path
          // has additional checks (e.g., regarding ghosting)
          Assert(vec_ptr == vec.begin() + constant_offset, ExcInternalError());
          process_dof_gather(indices,
                             vec,
                             constant_offset,
                             vec_ptr,
                             res,
                             std::bool_constant<false>());
        }
      else
        {
          res.gather(vec_ptr, indices);
        }
    }



    // variant where VectorType::value_type is not the same as Number -> must
    // manually load the data
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int *indices,
                       const VectorType   &vec,
                       const unsigned int  constant_offset,
                       typename VectorType::value_type *,
                       VectorizedArrayType &res,
                       std::bool_constant<false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        res[v] = vector_access(vec, indices[v] + constant_offset);
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       const VectorType             &vec,
                       Number                       &res) const
    {
      res = vec[index];
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
                       const VectorType  &vec,
                       Number            &res) const
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


    template <typename VectorNumberType>
    void
    process_dof(VectorNumberType &global, Number &local) const
    {
      global += local;
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType          &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<true>) const
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
                            VectorType          &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<false>) const
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
                                      const unsigned int  *dof_indices,
                                      VectorType          &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<true>) const
    {
      vectorized_transpose_and_store(true,
                                     dofs_per_cell,
                                     dof_values,
                                     dof_indices,
                                     vec.begin() + constant_offset);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int  *dof_indices,
                                      VectorType          &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<false>) const
    {
      for (unsigned int d = 0; d < dofs_per_cell; ++d)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access_add(vec,
                            dof_indices[v] + constant_offset + d,
                            dof_values[d][v]);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int       dofs_per_cell,
                                      const unsigned int      *dof_indices,
                                      VectorType              &vec,
                                      VectorizedArrayType     *dof_values,
                                      std::bool_constant<true> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int        dofs_per_cell,
                                      const unsigned int       *dof_indices,
                                      VectorType               &vec,
                                      VectorizedArrayType      *dof_values,
                                      std::bool_constant<false> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int                                  dofs_per_cell,
      std::array<Number2 *, VectorizedArrayType::size()> &global_ptr,
      VectorizedArrayType                                *dof_values,
      std::bool_constant<true>) const
    {
      vectorized_transpose_and_store(true,
                                     dofs_per_cell,
                                     dof_values,
                                     global_ptr);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int,
      std::array<Number2 *, VectorizedArrayType::size()> &,
      VectorizedArrayType *,
      std::bool_constant<false>) const
    {
      DEAL_II_NOT_IMPLEMENTED();
    }



    // variant where VectorType::value_type is the same as Number -> can call
    // scatter
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int              *indices,
                       VectorType                      &vec,
                       const unsigned int               constant_offset,
                       typename VectorType::value_type *vec_ptr,
                       VectorizedArrayType             &res,
                       std::bool_constant<true>) const
    {
      (void)constant_offset;
      (void)vec_ptr;
      (void)vec;

#if DEAL_II_VECTORIZATION_WIDTH_IN_BITS < 512
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access(vec, indices[v] + constant_offset) += res[v];
#else
      // only use gather in case there is also scatter.
      VectorizedArrayType tmp;
      tmp.gather(vec_ptr, indices);
      tmp += res;
      tmp.scatter(indices, vec_ptr);
#endif
    }



    // variant where VectorType::value_type is not the same as Number -> must
    // manually append all data
    template <typename VectorType>
    void
    process_dof_gather(const unsigned int *indices,
                       VectorType         &vec,
                       const unsigned int  constant_offset,
                       typename VectorType::value_type *,
                       VectorizedArrayType &res,
                       std::bool_constant<false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access_add(vec, indices[v] + constant_offset, res[v]);
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       VectorType                   &vec,
                       Number                       &res) const
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
                       VectorType        &vec,
                       Number            &res) const
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



    template <typename VectorNumberType>
    void
    process_dof(VectorNumberType &global, Number &local) const
    {
      global = local;
    }



    template <typename VectorType>
    void
    process_dofs_vectorized(const unsigned int   dofs_per_cell,
                            const unsigned int   dof_index,
                            VectorType          &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<true>) const
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
                            VectorType          &vec,
                            VectorizedArrayType *dof_values,
                            std::bool_constant<false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access(vec, dof_index + v + i * VectorizedArrayType::size()) =
            dof_values[i][v];
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int  *dof_indices,
                                      VectorType          &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<true>) const
    {
      vectorized_transpose_and_store(false,
                                     dofs_per_cell,
                                     dof_values,
                                     dof_indices,
                                     vec.begin() + constant_offset);
    }



    template <typename VectorType, bool booltype>
    void
    process_dofs_vectorized_transpose(const unsigned int   dofs_per_cell,
                                      const unsigned int  *dof_indices,
                                      VectorType          &vec,
                                      const unsigned int   constant_offset,
                                      VectorizedArrayType *dof_values,
                                      std::bool_constant<false>) const
    {
      for (unsigned int i = 0; i < dofs_per_cell; ++i)
        for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
          vector_access(vec, constant_offset + dof_indices[v] + i) =
            dof_values[i][v];
    }



    template <typename VectorType>
    void
    process_dofs_vectorized_transpose(const unsigned int       dofs_per_cell,
                                      const unsigned int      *dof_indices,
                                      VectorType              &vec,
                                      VectorizedArrayType     *dof_values,
                                      std::bool_constant<true> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename VectorType, bool booltype>
    void
    process_dofs_vectorized_transpose(const unsigned int        dofs_per_cell,
                                      const unsigned int       *dof_indices,
                                      VectorType               &vec,
                                      VectorizedArrayType      *dof_values,
                                      std::bool_constant<false> type) const
    {
      process_dofs_vectorized_transpose(
        dofs_per_cell, dof_indices, vec, 0, dof_values, type);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int                                  dofs_per_cell,
      std::array<Number2 *, VectorizedArrayType::size()> &global_ptr,
      VectorizedArrayType                                *dof_values,
      std::bool_constant<true>) const
    {
      vectorized_transpose_and_store(false,
                                     dofs_per_cell,
                                     dof_values,
                                     global_ptr);
    }



    template <typename Number2>
    void
    process_dofs_vectorized_transpose(
      const unsigned int,
      std::array<Number2 *, VectorizedArrayType::size()> &,
      VectorizedArrayType *,
      std::bool_constant<false>) const
    {
      DEAL_II_NOT_IMPLEMENTED();
    }



    template <typename VectorType>
    void
    process_dof_gather(const unsigned int              *indices,
                       VectorType                      &vec,
                       const unsigned int               constant_offset,
                       typename VectorType::value_type *vec_ptr,
                       VectorizedArrayType             &res,
                       std::bool_constant<true>) const
    {
      Assert(vec_ptr == vec.begin() + constant_offset, ExcInternalError());
      res.scatter(indices, vec_ptr);
    }



    template <typename VectorType>
    void
    process_dof_gather(const unsigned int *indices,
                       VectorType         &vec,
                       const unsigned int  constant_offset,
                       typename VectorType::value_type *,
                       VectorizedArrayType &res,
                       std::bool_constant<false>) const
    {
      for (unsigned int v = 0; v < VectorizedArrayType::size(); ++v)
        vector_access(vec, indices[v] + constant_offset) = res[v];
    }



    template <typename VectorType>
    void
    process_dof_global(const types::global_dof_index index,
                       VectorType                   &vec,
                       Number                       &res) const
    {
      vec[index] = res;
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
