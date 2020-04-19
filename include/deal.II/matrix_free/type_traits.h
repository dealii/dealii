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
  struct has_local_element
  {
  private:
    // this will work always.
    // we let it be void as we know T::local_element() (if exists) should
    // certainly return something
    static void
    detect(...);

    // this detecter will work only if we have "... T::local_element() const"
    // and its return type will be the same as local_element(),
    // that we expect to be T::value_type
    template <typename U>
    static decltype(std::declval<U const>().local_element(0))
    detect(const U &);

  public:
    // finally here we check if our detector has non-void return type
    // T::value_type. This will happen if compiler can use second detector,
    // otherwise SFINAE let it work with the more general first one that is void
    static const bool value =
      !std::is_same<void, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_local_element<T>::value;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::add_local_element(const uint, const typename T::value_type)
  template <typename T>
  struct has_add_local_element
  {
  private:
    static int
    detect(...);

    template <typename U>
    static decltype(
      std::declval<U>().add_local_element(0, typename T::value_type()))
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<int, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_add_local_element<T>::value;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::set_local_element(const uint, const typename T::value_type)
  template <typename T>
  struct has_set_local_element
  {
  private:
    static int
    detect(...);

    template <typename U>
    static decltype(
      std::declval<U>().set_local_element(0, typename T::value_type()))
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<int, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_set_local_element<T>::value;



  // same as above to check
  // bool T::partitioners_are_compatible(const Utilities::MPI::Partitioner &)
  // const
  template <typename T>
  struct has_partitioners_are_compatible
  {
  private:
    static void
    detect(...);

    template <typename U>
    static decltype(std::declval<U const>().partitioners_are_compatible(
      std::declval<Utilities::MPI::Partitioner>()))
    detect(const U &);

  public:
    static const bool value =
      std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_partitioners_are_compatible<T>::value;


  // same as above to check
  // ... T::begin() const
  template <typename T>
  struct has_begin
  {
  private:
    static void
    detect(...);

    template <typename U>
    static decltype(std::declval<U const>().begin())
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<void, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_begin<T>::value;


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
      (has_local_element<T>::value || is_serial_vector<T>::value) &&
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
  struct has_update_ghost_values_start
  {
  private:
    static bool
    detect(...);

    template <typename U>
    static decltype(std::declval<U const>().update_ghost_values_start(0))
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_update_ghost_values_start<T>::value;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // void T::	compress_start(const uint, VectorOperation::values)
  template <typename T>
  struct has_compress_start
  {
  private:
    static bool
    detect(...);

    template <typename U>
    static decltype(std::declval<U>().compress_start(0, VectorOperation::add))
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<bool, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_compress_start<T>::value;



  // type trait for vector T to see if
  // we do a custom data exchange route.
  // We assume that if both begin() and local_element()
  // exist, then begin() + offset == local_element(offset)
  template <typename T>
  struct has_exchange_on_subset
  {
    static const bool value =
      has_begin<T>::value && has_local_element<T>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_exchange_on_subset<T>::value;



  // a helper type-trait that leverage SFINAE to figure out if type T has
  // T::communication_block_size
  template <typename T>
  struct has_communication_block_size
  {
  private:
    static void
    detect(...);

    template <typename U>
    static decltype(U::communication_block_size)
    detect(const U &);

  public:
    static const bool value =
      !std::is_same<void, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool has_communication_block_size<T>::value;



  // type trait for vector T to see if
  // we need to do any data exchange for this vector type at all.
  // is_serial_vector<> would have been enough, but in some circumstances
  // (like calculation of diagonals for matrix-free operators)
  // a dummy InVector == unsigned int is provided.
  // Thus we have to treat this case as well.
  template <typename T>
  struct is_serial_or_dummy
  {
  private:
    // catches all cases including unsigned int
    static void
    detect(...);

    // catches serial vectors
    template <
      typename U,
      typename std::enable_if<is_serial_vector<U>::value, U>::type * = nullptr>
    static void
    detect(const U &);

    // catches parallel vectors
    template <
      typename U,
      typename std::enable_if<!is_serial_vector<U>::value, U>::type * = nullptr>
    static bool
    detect(const U &);

  public:
    static const bool value =
      std::is_same<void, decltype(detect(std::declval<T>()))>::value;
  };

  // We need to have a separate declaration for static const members
  template <typename T>
  const bool is_serial_or_dummy<T>::value;


} // namespace internal
#endif

DEAL_II_NAMESPACE_CLOSE

#endif
