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


#ifndef dealii_matrix_free_type_traits_h
#define dealii_matrix_free_type_traits_h

// various type-traits used exclusively within the matrix-free framework

#include <deal.II/base/config.h>

#include <deal.II/base/partitioner.h>

#include <deal.II/lac/vector_type_traits.h>


DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
namespace internal
{
  //
  // type traits for FEEvaluation
  //

  // a helper type-trait that leverage SFINAE to figure out if type T has
  // ... T::local_element() const
  template <typename T>
  using local_element_t = decltype(std::declval<T const>().local_element(0));

  template <typename T>
  using has_local_element = is_detected<local_element_t, T>;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::add_local_element(const uint, const typename T::value_type)
  template <typename T>
  using add_local_element_t =
    decltype(std::declval<T>().add_local_element(0, typename T::value_type()));

  template <typename T>
  using has_add_local_element = is_detected<add_local_element_t, T>;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::set_local_element(const uint, const typename T::value_type)
  template <typename T>
  using set_local_element_t =
    decltype(std::declval<T>().set_local_element(0, typename T::value_type()));

  template <typename T>
  using has_set_local_element = is_detected<set_local_element_t, T>;



  // same as above to check
  // bool T::partitioners_are_compatible(const Utilities::MPI::Partitioner &)
  // const
  template <typename T>
  using partitioners_are_compatible_t =
    decltype(std::declval<T const>().partitioners_are_compatible(
      std::declval<Utilities::MPI::Partitioner>()));

  template <typename T>
  using has_partitioners_are_compatible =
    is_detected<partitioners_are_compatible_t, T>;



  // same as above to check
  // ... T::begin() const
  template <typename T>
  using begin_t = decltype(std::declval<T const>().begin());

  template <typename T>
  using has_begin = is_detected<begin_t, T>;



  // same as above to check
  // ... T::shared_vector_data() const
  template <typename T>
  using shared_vector_data_t =
    decltype(std::declval<T const>().shared_vector_data());

  template <typename T>
  using has_shared_vector_data = is_detected<shared_vector_data_t, T>;



  // type trait for vector T and Number to see if
  // we can do vectorized load/save.
  // for VectorReader and VectorDistributorLocalToGlobal we assume that
  // if both begin() and local_element()
  // exist, then begin() + offset == local_element(offset)
  template <typename T, typename Number>
  struct is_vectorizable
  {
    static const bool value =
      has_begin<T>::value &&
      (has_local_element<T>::value ||
       is_serial_vector<typename std::remove_const<T>::type>::value) &&
      std::is_same<typename T::value_type, Number>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T, typename Number>
  const bool is_vectorizable<T, Number>::value;


  //
  // type-traits for Matrix-Free
  //
  // similar to type traits in FEEvaluation, below we add type-traits
  // to distinguish between vectors that provide different interface for
  // operations like update_ghost_values(), compress(), etc.
  // see internal::has_local_element in fe_evaluation.h that documents
  // how those type traits work.

  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::update_ghost_values_start(const uint) const
  template <typename T>
  using update_ghost_values_start_t =
    decltype(std::declval<T const>().update_ghost_values_start(0));

  template <typename T>
  using has_update_ghost_values_start =
    is_detected<update_ghost_values_start_t, T>;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::	compress_start(const uint, VectorOperation::values)
  template <typename T>
  using compress_start_t =
    decltype(std::declval<T>().compress_start(0, VectorOperation::add));

  template <typename T>
  using has_compress_start = is_detected<compress_start_t, T>;



  // type trait for vector T to see if
  // we do a custom data exchange route.
  // We assume that if both begin() and local_element()
  // exist, then begin() + offset == local_element(offset)
  template <typename T>
  struct has_exchange_on_subset
  {
    static const bool value = has_begin<T>::value &&
                              has_local_element<T>::value &&
                              has_partitioners_are_compatible<T>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_exchange_on_subset<T>::value;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // T::communication_block_size
  template <typename T>
  using communication_block_size_t = decltype(T::communication_block_size);

  template <typename T>
  using has_communication_block_size =
    is_detected<communication_block_size_t, T>;



  // type trait for vector T to see if
  // we need to do any data exchange for this vector type at all.
  // is_serial_vector<> would have been enough, but in some circumstances
  // (like calculation of diagonals for matrix-free operators)
  // a dummy InVector == unsigned int is provided.
  // Thus we have to treat this case as well.
  template <class T, class IsSerialVectorNotSpecialized = void>
  using not_parallel_vector_t =
    std::integral_constant<bool, is_serial_vector<T>::value>;

  template <class T>
  using is_not_parallel_vector =
    detected_or_t<std::true_type, not_parallel_vector_t, T>;
} // namespace internal
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
