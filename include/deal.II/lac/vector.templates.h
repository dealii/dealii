// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_vector_templates_h
#define dealii_vector_templates_h


#include <deal.II/base/config.h>

#include <deal.II/base/numbers.h>
#include <deal.II/base/template_constraints.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/vector_operations_internal.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_vector_base.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_vector.h>
#endif

#ifdef DEAL_II_TRILINOS_WITH_TPETRA
#  include <deal.II/lac/trilinos_tpetra_vector.h>
#endif


#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>

DEAL_II_NAMESPACE_OPEN

template <typename Number>
Vector<Number>::Vector(const Vector<Number> &v)
{
  *this = v;
}



template <typename Number>
void
Vector<Number>::apply_givens_rotation(const std::array<Number, 3> &csr,
                                      const size_type              i,
                                      const size_type              k)
{
  auto        &V = *this;
  const Number t = V(i);
  V(i)           = csr[0] * V(i) + csr[1] * V(k);
  V(k)           = -csr[1] * t + csr[0] * V(k);
}



template <typename Number>
template <typename OtherNumber>
Vector<Number>::Vector(const Vector<OtherNumber> &v)
{
  *this = v;
}



#ifdef DEAL_II_WITH_PETSC
namespace internal
{
  template <typename Number>
  void
  copy_petsc_vector(const PETScWrappers::VectorBase &v,
                    ::dealii::Vector<Number>        &out)
  {
    if (v.size() == 0)
      {
        out.reinit(0);
        return;
      }
    // Create a sequential PETSc vector and then copy over the entries into
    // the deal.II vector.
    Vec        sequential_vector;
    VecScatter scatter_context;

    PetscErrorCode ierr =
      VecScatterCreateToAll(v, &scatter_context, &sequential_vector);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = VecScatterBegin(
      scatter_context, v, sequential_vector, INSERT_VALUES, SCATTER_FORWARD);
    AssertThrow(ierr == 0, ExcPETScError(ierr));
    ierr = VecScatterEnd(
      scatter_context, v, sequential_vector, INSERT_VALUES, SCATTER_FORWARD);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const PetscScalar *start_ptr;
    ierr = VecGetArrayRead(sequential_vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    const PETScWrappers::VectorBase::size_type v_size = v.size();
    if (out.size() != v_size)
      out.reinit(v_size, true);

    internal::VectorOperations::copy(start_ptr,
                                     start_ptr + out.size(),
                                     out.begin());
    ierr = VecRestoreArrayRead(sequential_vector, &start_ptr);
    AssertThrow(ierr == 0, ExcPETScError(ierr));

    ierr = VecScatterDestroy(&scatter_context);
    AssertNothrow(ierr == 0, ExcPETScError(ierr));
    ierr = VecDestroy(&sequential_vector);
    AssertNothrow(ierr == 0, ExcPETScError(ierr));
  }
} // namespace internal



template <typename Number>
Vector<Number>::Vector(const PETScWrappers::VectorBase &v)
{
  internal::copy_petsc_vector(v, *this);
}
#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number>::Vector(const TrilinosWrappers::MPI::Vector &v)
  : values(v.size())
{
  if (size() != 0)
    {
      // Copy the distributed vector to
      // a local one at all processors
      // that know about the original vector.
      // TODO: There could
      // be a better solution than
      // this, but it has not yet been
      // found.
      TrilinosWrappers::MPI::Vector localized_vector;
      localized_vector.reinit(complete_index_set(size()),
                              v.get_mpi_communicator());
      localized_vector.reinit(v, false, true);

      Assert(localized_vector.size() == size(),
             ExcDimensionMismatch(localized_vector.size(), size()));

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = localized_vector.trilinos_vector().ExtractView(&start_ptr);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      std::copy(start_ptr[0], start_ptr[0] + size(), begin());

      maybe_reset_thread_partitioner();
    }
}

#endif


#ifdef DEAL_II_TRILINOS_WITH_TPETRA

template <typename Number>
template <typename OtherNumber, typename MemorySpace>
Vector<Number>::Vector(
  const LinearAlgebra::TpetraWrappers::Vector<OtherNumber, MemorySpace> &v)
  : values(v.size())
{
  static_assert(
    std::is_same<Number, OtherNumber>::value,
    "TpetraWrappers::Vector and dealii::Vector must use the same number type here.");

  if (size() != 0)
    {
      // Copy the distributed vector to
      // a local one at all processors
      // that know about the original vector.
      LinearAlgebra::TpetraWrappers::TpetraTypes::VectorType<OtherNumber,
                                                             MemorySpace>
        localized_vector(
          complete_index_set(size())
            .template make_tpetra_map_rcp<
              LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<
                MemorySpace>>(),
          v.get_mpi_communicator());

      Teuchos::RCP<const LinearAlgebra::TpetraWrappers::TpetraTypes::ImportType<
        MemorySpace>>
        importer = Tpetra::createImport(v.trilinos_vector().getMap(),
                                        localized_vector.getMap());

      localized_vector.doImport(v.trilinos_vector(), *importer, Tpetra::INSERT);

      // get a kokkos view from the localized_vector
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto localized_vector_2d =
        localized_vector.template getLocalView<Kokkos::HostSpace>(
          Tpetra::Access::ReadOnly);
#  else
      localized_vector.template sync<Kokkos::HostSpace>();
      auto localized_vector_2d =
        localized_vector.template getLocalView<Kokkos::HostSpace>();
#  endif
      auto localized_vector_1d =
        Kokkos::subview(localized_vector_2d, Kokkos::ALL(), 0);
      const size_t local_length = localized_vector.getLocalLength();

      Kokkos::DefaultHostExecutionSpace         exec;
      Kokkos::View<Number *, Kokkos::HostSpace> values_view(values.data(),
                                                            local_length);
      Kokkos::deep_copy(exec, values_view, localized_vector_1d);
      exec.fence();
    }
}

#endif


template <typename Number>
inline Vector<Number> &
Vector<Number>::operator=(const Vector<Number> &v)
{
  if (PointerComparison::equal(this, &v))
    return *this;

  if (size() != v.size())
    reinit(v, true);

  if (0 < size())
    {
      dealii::internal::VectorOperations::Vector_copy<Number, Number> copier(
        v.begin(), begin());
      internal::VectorOperations::parallel_for(copier,
                                               0,
                                               size(),
                                               thread_loop_partitioner);
    }

  return *this;
}



template <typename Number>
template <typename Number2>
inline Vector<Number> &
Vector<Number>::operator=(const Vector<Number2> &v)
{
  if (size() != v.size())
    reinit(v, true);

  dealii::internal::VectorOperations::Vector_copy<Number, Number2> copier(
    v.begin(), begin());
  internal::VectorOperations::parallel_for(copier,
                                           0,
                                           size(),
                                           thread_loop_partitioner);

  return *this;
}



template <typename Number>
inline void
Vector<Number>::reinit(const size_type n, const bool omit_zeroing_entries)
{
  do_reinit(n, omit_zeroing_entries, true);
}



template <typename Number>
inline void
Vector<Number>::grow_or_shrink(const size_type n)
{
  values.resize(n);

  maybe_reset_thread_partitioner();
}



template <typename Number>
bool
Vector<Number>::all_zero() const
{
  Assert(size() != 0, ExcEmptyObject());

  for (size_type i = 0; i < size(); ++i)
    if (values[i] != Number())
      return false;
  return true;
}



template <typename Number>
bool
Vector<Number>::is_non_negative() const
{
  Assert(size() != 0, ExcEmptyObject());

  for (size_type i = 0; i < size(); ++i)
    if (!internal::VectorOperations::is_non_negative(values[i]))
      return false;

  return true;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator=(const Number s)
{
  AssertIsFinite(s);
  if (s != Number())
    Assert(size() != 0, ExcEmptyObject());

  if (size() > 0)
    {
      internal::VectorOperations::Vector_set<Number> setter(s, values.begin());
      internal::VectorOperations::parallel_for(setter,
                                               0,
                                               size(),
                                               thread_loop_partitioner);
    }

  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator*=(const Number factor)
{
  AssertIsFinite(factor);

  Assert(size() != 0, ExcEmptyObject());

  internal::VectorOperations::Vectorization_multiply_factor<Number>
    vector_multiply(values.begin(), factor);

  internal::VectorOperations::parallel_for(vector_multiply,
                                           0,
                                           size(),
                                           thread_loop_partitioner);

  return *this;
}



template <typename Number>
void
Vector<Number>::add(const Number a, const Vector<Number> &v)
{
  AssertIsFinite(a);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  internal::VectorOperations::Vectorization_add_av<Number> vector_add_av(
    values.begin(), v.values.begin(), a);
  internal::VectorOperations::parallel_for(vector_add_av,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::sadd(const Number x, const Number a, const Vector<Number> &v)
{
  AssertIsFinite(x);
  AssertIsFinite(a);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  internal::VectorOperations::Vectorization_sadd_xav<Number> vector_sadd_xav(
    values.begin(), v.values.begin(), a, x);
  internal::VectorOperations::parallel_for(vector_sadd_xav,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
Number
Vector<Number>::operator*(const Vector<Number2> &v) const
{
  Assert(size() != 0, ExcEmptyObject());

  if (PointerComparison::equal(this, &v))
    return norm_sqr();

  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  Number                                           sum;
  internal::VectorOperations::Dot<Number, Number2> dot(values.begin(),
                                                       v.values.begin());
  internal::VectorOperations::parallel_reduce(
    dot, 0, size(), sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::norm_sqr() const
{
  Assert(size() != 0, ExcEmptyObject());

  real_type                                            sum;
  internal::VectorOperations::Norm2<Number, real_type> norm2(values.begin());
  internal::VectorOperations::parallel_reduce(
    norm2, 0, size(), sum, thread_loop_partitioner);

  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
Number
Vector<Number>::mean_value() const
{
  Assert(size() != 0, ExcEmptyObject());

  Number                                        sum;
  internal::VectorOperations::MeanValue<Number> mean(values.begin());
  internal::VectorOperations::parallel_reduce(
    mean, 0, size(), sum, thread_loop_partitioner);

  return sum / real_type(size());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l1_norm() const
{
  Assert(size() != 0, ExcEmptyObject());

  real_type                                            sum;
  internal::VectorOperations::Norm1<Number, real_type> norm1(values.begin());
  internal::VectorOperations::parallel_reduce(
    norm1, 0, size(), sum, thread_loop_partitioner);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l2_norm() const
{
  // if l2_norm()^2 is finite and non-zero, the answer is computed as
  // std::sqrt(norm_sqr()). If norm_sqr() is infinite or zero, the l2 norm
  // might still be finite. In that case, recompute it (this is a rare case,
  // so working on the vector twice is uncritical and paid off by the extended
  // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
  Assert(size() != 0, ExcEmptyObject());

  real_type                                            norm_square;
  internal::VectorOperations::Norm2<Number, real_type> norm2(values.begin());
  internal::VectorOperations::parallel_reduce(
    norm2, 0, size(), norm_square, thread_loop_partitioner);
  if (numbers::is_finite(norm_square) &&
      norm_square >= std::numeric_limits<real_type>::min())
    return static_cast<typename Vector<Number>::real_type>(
      std::sqrt(norm_square));
  else
    {
      real_type scale = 0.;
      real_type sum   = 1.;
      for (size_type i = 0; i < size(); ++i)
        {
          if (values[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(values[i]);
              if (scale < abs_x)
                {
                  sum   = 1 + sum * (scale / abs_x) * (scale / abs_x);
                  scale = abs_x;
                }
              else
                sum += (abs_x / scale) * (abs_x / scale);
            }
        }
      AssertIsFinite(scale * std::sqrt(sum));
      return static_cast<typename Vector<Number>::real_type>(scale *
                                                             std::sqrt(sum));
    }
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::lp_norm(const real_type p) const
{
  Assert(size() != 0, ExcEmptyObject());

  if (p == 1.)
    return l1_norm();
  else if (p == 2.)
    return l2_norm();

  real_type                                            sum;
  internal::VectorOperations::NormP<Number, real_type> normp(values.begin(), p);
  internal::VectorOperations::parallel_reduce(
    normp, 0, size(), sum, thread_loop_partitioner);

  if (numbers::is_finite(sum) && sum >= std::numeric_limits<real_type>::min())
    return std::pow(sum, static_cast<real_type>(1. / p));
  else
    {
      real_type scale = 0.;
      real_type sum   = 1.;
      for (size_type i = 0; i < size(); ++i)
        {
          if (values[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(values[i]);
              if (scale < abs_x)
                {
                  sum   = 1. + sum * std::pow(scale / abs_x, p);
                  scale = abs_x;
                }
              else
                sum += std::pow(abs_x / scale, p);
            }
        }
      return scale * std::pow(sum, static_cast<real_type>(1. / p));
    }
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::linfty_norm() const
{
  Assert(size() != 0, ExcEmptyObject());

  real_type max = 0.;

  for (size_type i = 0; i < size(); ++i)
    max = std::max(numbers::NumberTraits<Number>::abs(values[i]), max);

  return max;
}



template <typename Number>
Number
Vector<Number>::add_and_dot(const Number          a,
                            const Vector<Number> &V,
                            const Vector<Number> &W)
{
  Assert(size() != 0, ExcEmptyObject());
  AssertDimension(size(), V.size());
  AssertDimension(size(), W.size());

  Number                                        sum;
  internal::VectorOperations::AddAndDot<Number> adder(values.begin(),
                                                      V.values.begin(),
                                                      W.values.begin(),
                                                      a);
  internal::VectorOperations::parallel_reduce(
    adder, 0, size(), sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
void
Vector<Number>::extract_subvector_to(
  const ArrayView<const types::global_dof_index> &indices,
  const ArrayView<Number>                        &elements) const
{
  AssertDimension(indices.size(), elements.size());
  for (unsigned int i = 0; i < indices.size(); ++i)
    {
      AssertIndexRange(indices[i], size());
      elements[i] = (*this)[indices[i]];
    }
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator+=(const Vector<Number> &v)
{
  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  internal::VectorOperations::Vectorization_add_v<Number> vector_add(
    values.begin(), v.values.begin());
  internal::VectorOperations::parallel_for(vector_add,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator-=(const Vector<Number> &v)
{
  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  internal::VectorOperations::Vectorization_subtract_v<Number> vector_subtract(
    values.begin(), v.values.begin());
  internal::VectorOperations::parallel_for(vector_subtract,
                                           0,
                                           size(),
                                           thread_loop_partitioner);

  return *this;
}



template <typename Number>
void
Vector<Number>::add(const Number v)
{
  Assert(size() != 0, ExcEmptyObject());

  internal::VectorOperations::Vectorization_add_factor<Number> vector_add(
    values.begin(), v);
  internal::VectorOperations::parallel_for(vector_add,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::add(const Number          a,
                    const Vector<Number> &v,
                    const Number          b,
                    const Vector<Number> &w)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));
  Assert(size() == w.size(), ExcDimensionMismatch(size(), w.size()));

  internal::VectorOperations::Vectorization_add_avpbw<Number> vector_add(
    values.begin(), v.values.begin(), w.values.begin(), a, b);
  internal::VectorOperations::parallel_for(vector_add,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::sadd(const Number x, const Vector<Number> &v)
{
  AssertIsFinite(x);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  internal::VectorOperations::Vectorization_sadd_xv<Number> vector_sadd(
    values.begin(), v.values.begin(), x);
  internal::VectorOperations::parallel_for(vector_sadd,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::scale(const Vector<Number> &s)
{
  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == s.size(), ExcDimensionMismatch(size(), s.size()));

  internal::VectorOperations::Vectorization_scale<Number> vector_scale(
    values.begin(), s.values.begin());
  internal::VectorOperations::parallel_for(vector_scale,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void
Vector<Number>::scale(const Vector<Number2> &s)
{
  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == s.size(), ExcDimensionMismatch(size(), s.size()));

  for (size_type i = 0; i < size(); ++i)
    values[i] *= Number(s.values[i]);
}



template <typename Number>
void
Vector<Number>::equ(const Number a, const Vector<Number> &u)
{
  AssertIsFinite(a);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == u.size(), ExcDimensionMismatch(size(), u.size()));

  internal::VectorOperations::Vectorization_equ_au<Number> vector_equ(
    values.begin(), u.values.begin(), a);
  internal::VectorOperations::parallel_for(vector_equ,
                                           0,
                                           size(),
                                           thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void
Vector<Number>::equ(const Number a, const Vector<Number2> &u)
{
  AssertIsFinite(a);

  Assert(size() != 0, ExcEmptyObject());
  Assert(size() == u.size(), ExcDimensionMismatch(size(), u.size()));

  // set the result vector to a*u. we have to
  // convert the elements of u to the type of
  // the result vector. this is necessary
  // because
  // operator*(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i = 0; i < size(); ++i)
    values[i] = a * Number(u.values[i]);
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator=(const BlockVector<Number> &v)
{
  if (v.size() != size())
    reinit(v.size(), true);

  size_type this_index = 0;
  for (size_type b = 0; b < v.n_blocks(); ++b)
    for (size_type i = 0; i < v.block(b).size(); ++i, ++this_index)
      values[this_index] = v.block(b)(i);

  return *this;
}



#ifdef DEAL_II_WITH_PETSC
template <typename Number>
Vector<Number> &
Vector<Number>::operator=(const PETScWrappers::VectorBase &v)
{
  internal::copy_petsc_vector(v, *this);
  return *this;
}
#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number> &
Vector<Number>::operator=(const TrilinosWrappers::MPI::Vector &v)
{
  if (v.size() != size())
    reinit(v.size(), true);
  if (size() != 0)
    {
      // Copy the distributed vector to
      // a local one at all processors
      // that know about the original vector.
      // TODO: There could
      // be a better solution than
      // this, but it has not yet been
      // found.
      TrilinosWrappers::MPI::Vector localized_vector;
      localized_vector.reinit(complete_index_set(size()),
                              v.get_mpi_communicator());
      localized_vector.reinit(v, false, true);

      Assert(localized_vector.size() == size(),
             ExcDimensionMismatch(localized_vector.size(), size()));

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = localized_vector.trilinos_vector().ExtractView(&start_ptr);
      AssertThrow(ierr == 0, ExcTrilinosError(ierr));

      std::copy(start_ptr[0], start_ptr[0] + size(), begin());
    }

  return *this;
}

#endif



#ifdef DEAL_II_TRILINOS_WITH_TPETRA

template <typename Number>
template <typename OtherNumber, typename MemorySpace>
Vector<Number> &
Vector<Number>::operator=(
  const LinearAlgebra::TpetraWrappers::Vector<OtherNumber, MemorySpace> &v)
{
  static_assert(
    std::is_same<Number, OtherNumber>::value,
    "TpetraWrappers::Vector and dealii::Vector must use the same number type here.");

  if (v.size() != size())
    reinit(v.size(), true);

  if (size() != 0)
    {
      // Copy the distributed vector to
      // a local one at all processors
      // that know about the original vector.
      LinearAlgebra::TpetraWrappers::TpetraTypes::VectorType<OtherNumber,
                                                             MemorySpace>
        localized_vector(
          complete_index_set(size())
            .template make_tpetra_map_rcp<
              LinearAlgebra::TpetraWrappers::TpetraTypes::NodeType<
                MemorySpace>>(),
          v.get_mpi_communicator());

      Teuchos::RCP<const LinearAlgebra::TpetraWrappers::TpetraTypes::ImportType<
        MemorySpace>>
        importer = Tpetra::createImport(v.trilinos_vector().getMap(),
                                        localized_vector.getMap());

      localized_vector.doImport(v.trilinos_vector(), *importer, Tpetra::INSERT);

      // get a kokkos view from the localized_vector
#  if DEAL_II_TRILINOS_VERSION_GTE(13, 2, 0)
      auto localized_vector_2d =
        localized_vector.template getLocalView<Kokkos::HostSpace>(
          Tpetra::Access::ReadOnly);
#  else
      localized_vector.template sync<Kokkos::HostSpace>();
      auto localized_vector_2d =
        localized_vector.template getLocalView<Kokkos::HostSpace>();
#  endif
      auto localized_vector_1d =
        Kokkos::subview(localized_vector_2d, Kokkos::ALL(), 0);
      const size_t local_length = localized_vector.getLocalLength();

      Kokkos::DefaultHostExecutionSpace         exec;
      Kokkos::View<Number *, Kokkos::HostSpace> values_view(values.data(),
                                                            local_length);
      Kokkos::deep_copy(exec, values_view, localized_vector_1d);
      exec.fence();
    }

  return *this;
}

#endif



template <typename Number>
template <typename Number2>
bool
Vector<Number>::operator==(const Vector<Number2> &v) const
{
  Assert(size() == v.size(), ExcDimensionMismatch(size(), v.size()));

  // compare the two vector. we have to
  // convert the elements of v to the type of
  // the result vector. this is necessary
  // because
  // operator==(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i = 0; i < size(); ++i)
    if (values[i] != Number(v.values[i]))
      return false;

  return true;
}



template <typename Number>
void
Vector<Number>::print(std::ostream      &out,
                      const unsigned int precision,
                      const bool         scientific,
                      const bool         across) const
{
  Assert(size() != 0, ExcEmptyObject());
  AssertThrow(out.fail() == false, ExcIO());

  std::ios::fmtflags old_flags     = out.flags();
  unsigned int       old_precision = out.precision(precision);

  out.precision(precision);
  if (scientific)
    out.setf(std::ios::scientific, std::ios::floatfield);
  else
    out.setf(std::ios::fixed, std::ios::floatfield);

  if (across)
    for (size_type i = 0; i < size(); ++i)
      out << values[i] << ' ';
  else
    for (size_type i = 0; i < size(); ++i)
      out << values[i] << std::endl;
  out << std::endl;

  AssertThrow(out.fail() == false, ExcIO());
  // reset output format
  out.flags(old_flags);
  out.precision(old_precision);
}



template <typename Number>
void
Vector<Number>::block_write(std::ostream &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  out << std::to_string(size()) << "\n[";
  out.write(reinterpret_cast<const char *>(begin()),
            reinterpret_cast<const char *>(end()) -
              reinterpret_cast<const char *>(begin()));
  out << ']';

  AssertThrow(out.fail() == false, ExcIO());
}



template <typename Number>
void
Vector<Number>::block_read(std::istream &in)
{
  AssertThrow(in.fail() == false, ExcIO());

  size_type sz;

  char buf[16];


  in.getline(buf, 16, '\n');
  sz = std::atoi(buf);

  // fast initialization, since the
  // data elements are overwritten anyway
  reinit(sz, true);

  char c;
  in.read(&c, 1);
  AssertThrow(c == '[', ExcIO());

  in.read(reinterpret_cast<char *>(begin()),
          reinterpret_cast<const char *>(end()) -
            reinterpret_cast<const char *>(begin()));

  in.read(&c, 1);
  AssertThrow(c == ']', ExcIO());
}



template <typename Number>
IndexSet
Vector<Number>::locally_owned_elements() const
{
  return complete_index_set(size());
}



template <typename Number>
std::size_t
Vector<Number>::memory_consumption() const
{
  return sizeof(*this) + values.memory_consumption() - sizeof(values);
}



template <typename Number>
void
Vector<Number>::maybe_reset_thread_partitioner()
{
  if (size() >= 4 * internal::VectorImplementation::minimum_parallel_grain_size)
    {
      if (thread_loop_partitioner == nullptr)
        thread_loop_partitioner =
          std::make_shared<parallel::internal::TBBPartitioner>();
    }
  else
    thread_loop_partitioner.reset();
}



template <typename Number>
void
Vector<Number>::do_reinit(const size_type new_size,
                          const bool      omit_zeroing_entries,
                          const bool      reset_partitioner)
{
  values.resize_fast(new_size);
  if (!omit_zeroing_entries)
    values.fill();

  if (reset_partitioner)
    maybe_reset_thread_partitioner();
}


DEAL_II_NAMESPACE_CLOSE

#endif
