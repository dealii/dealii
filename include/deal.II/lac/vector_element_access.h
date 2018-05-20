// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef dealii_vector_element_access_h
#define dealii_vector_element_access_h

#include <deal.II/lac/trilinos_epetra_vector.h>

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
        VectorType&                           V);

    static void
    set(typename VectorType::value_type value,
        const types::global_dof_index   i,
        VectorType&                     V);

    static typename VectorType::value_type
    get(const VectorType& V, const types::global_dof_index i);
  };

  template <typename VectorType>
  inline void
  ElementAccess<VectorType>::add(const typename VectorType::value_type value,
                                 const types::global_dof_index         i,
                                 VectorType&                           V)
  {
    V(i) += value;
  }

  template <typename VectorType>
  inline void
  ElementAccess<VectorType>::set(const typename VectorType::value_type value,
                                 const types::global_dof_index         i,
                                 VectorType&                           V)
  {
    V(i) = value;
  }

  template <typename VectorType>
  inline typename VectorType::value_type
  ElementAccess<VectorType>::get(const VectorType&             V,
                                 const types::global_dof_index i)
  {
    return V(i);
  }

#if defined(DEAL_II_WITH_TRILINOS) && defined(DEAL_II_WITH_MPI)
  template <>
  inline void
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::add(
    const double                           value,
    const types::global_dof_index          i,
    LinearAlgebra::EpetraWrappers::Vector& V)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i
      = vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    vector[0][trilinos_i] += value;
  }

  template <>
  inline void
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::set(
    const double                           value,
    const types::global_dof_index          i,
    LinearAlgebra::EpetraWrappers::Vector& V)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i
      = vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    vector[0][trilinos_i] = value;
  }

  template <>
  inline double
  ElementAccess<LinearAlgebra::EpetraWrappers::Vector>::get(
    const LinearAlgebra::EpetraWrappers::Vector& V,
    const types::global_dof_index                i)
  {
    // Extract local indices in the vector.
    Epetra_FEVector                   vector = V.trilinos_vector();
    TrilinosWrappers::types::int_type trilinos_i
      = vector.Map().LID(static_cast<TrilinosWrappers::types::int_type>(i));

    return vector[0][trilinos_i];
  }
#endif
} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
