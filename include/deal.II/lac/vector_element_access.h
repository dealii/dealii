// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_element_access_h
#define dealii_vector_element_access_h


#include <deal.II/base/config.h>

#include <deal.II/lac/trilinos_epetra_vector.h>
#include <deal.II/lac/trilinos_tpetra_block_vector.h>
#include <deal.II/lac/trilinos_tpetra_vector.h>

DEAL_II_NAMESPACE_OPEN

namespace internal
{
  template <typename VectorType>
  struct ElementAccess
  {
  public:
    static void
    add(const typename VectorType::value_type value,
        const types::global_dof_index         i,
        VectorType                           &V);

    static void
    set(typename VectorType::value_type value,
        const types::global_dof_index   i,
        VectorType                     &V);

    static typename VectorType::value_type
    get(const VectorType &V, const types::global_dof_index i);
  };



  template <typename VectorType>
  inline void
  ElementAccess<VectorType>::add(const typename VectorType::value_type value,
                                 const types::global_dof_index         i,
                                 VectorType                           &V)
  {
    V(i) += value;
  }



  template <typename VectorType>
  inline void
  ElementAccess<VectorType>::set(const typename VectorType::value_type value,
                                 const types::global_dof_index         i,
                                 VectorType                           &V)
  {
    V(i) = value;
  }



  template <typename VectorType>
  inline typename VectorType::value_type
  ElementAccess<VectorType>::get(const VectorType             &V,
                                 const types::global_dof_index i)
  {
    return V(i);
  }



#ifdef DEAL_II_WITH_TRILINOS
  template <>
  inline void
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::add(
    const double                           value,
    const types::global_dof_index          i,
    LinearAlgebra::EpetraWrappers::Vector &V)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i =
      vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    vector[0][trilinos_i] += value;
  }



  template <>
  inline void
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::set(
    const double                           value,
    const types::global_dof_index          i,
    LinearAlgebra::EpetraWrappers::Vector &V)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i =
      vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    vector[0][trilinos_i] = value;
  }

  template <>
  inline double
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::get(
    const LinearAlgebra::EpetraWrappers::Vector &V,
    const types::global_dof_index                i)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i =
      vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    return vector[0][trilinos_i];
  }



#  ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename NumberType, typename MemorySpace>
  struct ElementAccess<
    LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace>>
  {
  public:
    using VectorType =
      LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace>;
    static void
    add(const typename VectorType::value_type value,
        const types::global_dof_index         i,
        VectorType                           &V);

    static void
    set(typename VectorType::value_type value,
        const types::global_dof_index   i,
        VectorType                     &V);

    static typename VectorType::value_type
    get(const VectorType &V, const types::global_dof_index i);
  };



  template <typename NumberType, typename MemorySpace>
  inline void
  ElementAccess<
    LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace>>::
    add(const typename VectorType::value_type                           value,
        const types::global_dof_index                                   i,
        LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace> &V)
  {
    // Extract local indices in the vector.
    auto                              vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i =
      vector.getMap()->getLocalElement(
        static_cast<TrilinosWrappers::types::int_type>(i));

#    if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    auto vector_2d = vector.template getLocalView<Kokkos::HostSpace>(
      Tpetra::Access::ReadWrite);
#    else
    vector.template sync<Kokkos::HostSpace>();
    auto vector_2d = vector.template getLocalView<Kokkos::HostSpace>();
#    endif
    auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
#    if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    // We're going to modify the data on host.
    vector.template modify<Kokkos::HostSpace>();
#    endif
    vector_1d(trilinos_i) += value;
#    if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    vector.template sync<
      typename Tpetra::Vector<NumberType, int, types::signed_global_dof_index>::
        device_type::memory_space>();
#    endif
  }



  template <typename NumberType, typename MemorySpace>
  inline void
  ElementAccess<
    LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace>>::
    set(const typename VectorType::value_type                           value,
        const types::global_dof_index                                   i,
        LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace> &V)
  {
    // Extract local indices in the vector.
    auto                              vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i =
      vector.getMap()->getLocalElement(
        static_cast<TrilinosWrappers::types::int_type>(i));

#    if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    auto vector_2d = vector.template getLocalView<Kokkos::HostSpace>(
      Tpetra::Access::ReadWrite);
#    else
    vector.template sync<Kokkos::HostSpace>();
    auto vector_2d = vector.template getLocalView<Kokkos::HostSpace>();
#    endif
    auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
    // We're going to modify the data on host.
#    if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    vector.template modify<Kokkos::HostSpace>();
#    endif
    vector_1d(trilinos_i) = value;
#    if !DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    vector.template sync<
      typename Tpetra::Vector<NumberType, int, types::signed_global_dof_index>::
        device_type::memory_space>();
#    endif
  }



  template <typename NumberType, typename MemorySpace>
  inline typename LinearAlgebra::TpetraWrappers::Vector<NumberType,
                                                        MemorySpace>::value_type
  ElementAccess<
    LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace>>::
    get(const LinearAlgebra::TpetraWrappers::Vector<NumberType, MemorySpace> &V,
        const types::global_dof_index                                         i)
  {
    // Extract local indices in the vector.
#    if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
    const auto &vector = V.trilinos_vector();
    auto        vector_2d =
      vector.template getLocalView<Kokkos::HostSpace>(Tpetra::Access::ReadOnly);
#    else
    auto vector    = V.trilinos_vector();
    vector.template sync<Kokkos::HostSpace>();
    auto vector_2d = vector.template getLocalView<Kokkos::HostSpace>();
#    endif
    auto vector_1d = Kokkos::subview(vector_2d, Kokkos::ALL(), 0);
    TrilinosWrappers::types::int_type trilinos_i =
      vector.getMap()->getLocalElement(
        static_cast<TrilinosWrappers::types::int_type>(i));
    return vector_1d(trilinos_i);
  }
#  endif
#endif
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
