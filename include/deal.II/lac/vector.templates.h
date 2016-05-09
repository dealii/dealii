// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#ifndef dealii__vector_templates_h
#define dealii__vector_templates_h


#include <deal.II/base/template_constraints.h>
#include <deal.II/base/numbers.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector_operations_internal.h>

#ifdef DEAL_II_WITH_PETSC
#  include <deal.II/lac/petsc_vector.h>
#  include <deal.II/lac/petsc_parallel_vector.h>
#endif

#ifdef DEAL_II_WITH_TRILINOS
#  include <deal.II/lac/trilinos_vector.h>
#endif


#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN

template <typename Number>
Vector<Number>::Vector (const Vector<Number> &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();
      *this = v;
    }
}



#ifdef DEAL_II_WITH_CXX11
template <typename Number>
Vector<Number>::Vector (Vector<Number> &&v)
  :
  Subscriptor(std::move(v)),
  vec_size(v.vec_size),
  max_vec_size(v.max_vec_size),
  val(v.val),
  thread_loop_partitioner(std::move(v.thread_loop_partitioner))
{
  v.vec_size = 0;
  v.max_vec_size = 0;
  v.val = nullptr;
}
#endif



#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG
template <typename Number>
template <typename OtherNumber>
Vector<Number>::Vector (const Vector<OtherNumber> &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();
      *this = v;
    }
}
#endif


#ifdef DEAL_II_WITH_PETSC

template <typename Number>
Vector<Number>::Vector (const PETScWrappers::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // get a representation of the vector
      // and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

      // restore the representation of the
      // vector
      ierr = VecRestoreArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }
}



template <typename Number>
Vector<Number>::Vector (const PETScWrappers::MPI::Vector &v)
  :
  Subscriptor(),
  vec_size(0),
  max_vec_size(0),
  val(0)
{
  if (v.size() != 0)
    {
      // do this in a two-stage process:
      // first convert to a sequential petsc
      // vector, then copy that
      PETScWrappers::Vector seq (v);
      *this = seq;
    }
}

#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number>::Vector (const TrilinosWrappers::MPI::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // Copy the distributed vector to
      // a local one at all
      // processors. TODO: There could
      // be a better solution than
      // this, but it has not yet been
      // found.
      TrilinosWrappers::Vector localized_vector (v);

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = localized_vector.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }
}



template <typename Number>
Vector<Number>::Vector (const TrilinosWrappers::Vector &v)
  :
  Subscriptor(),
  vec_size(v.size()),
  max_vec_size(v.size()),
  val(0)
{
  if (vec_size != 0)
    {
      allocate();

      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;

      int ierr = v.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }
}

#endif


template <typename Number>
inline
Vector<Number> &
Vector<Number>::operator= (const Vector<Number> &v)
{
  if (PointerComparison::equal(this, &v))
    return *this;

  thread_loop_partitioner = v.thread_loop_partitioner;
  if (vec_size != v.vec_size)
    reinit (v, true);

  dealii::internal::Vector_copy<Number,Number> copier(v.val, val);
  internal::parallel_for(copier,vec_size,thread_loop_partitioner);

  return *this;
}



#ifdef DEAL_II_WITH_CXX11
template <typename Number>
inline
Vector<Number> &
Vector<Number>::operator= (Vector<Number> &&v)
{
  Subscriptor::operator=(std::move(v));

  if (val) deallocate();

  vec_size = v.vec_size;
  max_vec_size = v.max_vec_size;
  val = v.val;
  thread_loop_partitioner = std::move(v.thread_loop_partitioner);

  v.vec_size = 0;
  v.max_vec_size = 0;
  v.val = nullptr;

  return *this;
}
#endif



template <typename Number>
template <typename Number2>
inline
Vector<Number> &
Vector<Number>::operator= (const Vector<Number2> &v)
{
  thread_loop_partitioner = v.thread_loop_partitioner;
  if (vec_size != v.vec_size)
    reinit (v, true);

  dealii::internal::Vector_copy<Number,Number2> copier(v.val, val);
  internal::parallel_for(copier,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
inline
void Vector<Number>::reinit (const size_type n,
                             const bool omit_zeroing_entries)
{
  if (n==0)
    {
      if (val) deallocate();
      val = 0;
      max_vec_size = vec_size = 0;
      thread_loop_partitioner.reset(new parallel::internal::TBBPartitioner());
      return;
    };

  if (n>max_vec_size)
    {
      if (val) deallocate();
      max_vec_size = n;
      allocate();
    };

  if (vec_size != n)
    {
      vec_size = n;

      // only reset the partitioner if we actually expect a significant vector
      // size
      if (vec_size >= 4*internal::Vector::minimum_parallel_grain_size)
        thread_loop_partitioner.reset(new parallel::internal::TBBPartitioner());
    }

  if (omit_zeroing_entries == false)
    *this = static_cast<Number>(0);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::reinit (const Vector<Number2> &v,
                             const bool omit_zeroing_entries)
{
  thread_loop_partitioner = v.thread_loop_partitioner;

  if (v.vec_size==0)
    {
      if (val) deallocate();
      val = 0;
      max_vec_size = vec_size = 0;
      return;
    };

  if (v.vec_size>max_vec_size)
    {
      if (val) deallocate();
      max_vec_size = v.vec_size;
      allocate();
    };
  vec_size = v.vec_size;
  if (omit_zeroing_entries == false)
    *this = static_cast<Number>(0);
}



template <typename Number>
bool
Vector<Number>::all_zero () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type i=0; i<vec_size; ++i)
    if (val[i] != Number(0))
      return false;
  return true;
}



template <typename Number>
bool
Vector<Number>::is_non_negative () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type i=0; i<vec_size; ++i)
    if ( ! internal::is_non_negative (val[i]))
      return false;

  return true;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const Number s)
{
  AssertIsFinite(s);
  if (s != Number())
    Assert (vec_size!=0, ExcEmptyObject());

  internal::Vector_set<Number> setter(s, val);

  internal::parallel_for(setter,vec_size,thread_loop_partitioner);

  return *this;
}



#ifdef DEAL_II_BOOST_BIND_COMPILER_BUG
template <>
Vector<std::complex<float> > &
Vector<std::complex<float> >::operator= (const std::complex<float> s)
{
  AssertIsFinite(s);
  if (s != std::complex<float>())
    Assert (vec_size!=0, ExcEmptyObject());
  if (vec_size!=0)
    std::fill (begin(), end(), s);

  return *this;
}
#endif



template <typename Number>
Vector<Number> &Vector<Number>::operator *= (const Number factor)
{
  AssertIsFinite(factor);

  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_multiply_factor<Number> vector_multiply(val, factor);

  internal::parallel_for(vector_multiply,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
void
Vector<Number>::add (const Number a,
                     const Vector<Number> &v)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_av<Number> vector_add_av(val, v.val, a);
  internal::parallel_for(vector_add_av,vec_size,thread_loop_partitioner);
}



template <typename Number>
void
Vector<Number>::sadd (const Number x,
                      const Number a,
                      const Vector<Number> &v)
{
  AssertIsFinite(x);
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xav<Number> vector_sadd_xav(val, v.val, a, x);
  internal::parallel_for(vector_sadd_xav,vec_size,thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
Number Vector<Number>::operator * (const Vector<Number2> &v) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (PointerComparison::equal (this, &v))
    return norm_sqr();

  Assert (vec_size == v.size(),
          ExcDimensionMismatch(vec_size, v.size()));

  Number sum;
  internal::Dot<Number,Number2> dot(val, v.val);
  internal::parallel_reduce (dot, vec_size, sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::norm_sqr () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Norm2<Number,real_type> norm2(val);
  internal::parallel_reduce (norm2, vec_size, sum, thread_loop_partitioner);

  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  Number sum;
  internal::MeanValue<Number> mean(val);
  internal::parallel_reduce (mean, vec_size, sum, thread_loop_partitioner);

  return sum / real_type(size());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l1_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Norm1<Number, real_type> norm1(val);
  internal::parallel_reduce (norm1, vec_size, sum, thread_loop_partitioner);

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l2_norm () const
{
  // if l2_norm()^2 is finite and non-zero, the answer is computed as
  // std::sqrt(norm_sqr()). If norm_sqr() is infinite or zero, the l2 norm
  // might still be finite. In that case, recompute it (this is a rare case,
  // so working on the vector twice is uncritical and paid off by the extended
  // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
  Assert (vec_size!=0, ExcEmptyObject());

  real_type norm_square;
  internal::Norm2<Number, real_type> norm2(val);
  internal::parallel_reduce (norm2, vec_size, norm_square,
                             thread_loop_partitioner);
  if (numbers::is_finite(norm_square) &&
      norm_square >= std::numeric_limits<real_type>::min())
    return std::sqrt(norm_square);
  else
    {
      real_type scale = 0.;
      real_type sum = 1.;
      for (size_type i=0; i<vec_size; ++i)
        {
          if (val[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(val[i]);
              if (scale < abs_x)
                {
                  sum = 1. + sum * (scale/abs_x) * (scale/abs_x);
                  scale = abs_x;
                }
              else
                sum += (abs_x/scale) * (abs_x/scale);
            }
        }
      AssertIsFinite(scale*std::sqrt(sum));
      return scale * std::sqrt(sum);
    }
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::lp_norm (const real_type p) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (p == 1.)
    return l1_norm();
  else if (p == 2.)
    return l2_norm();

  real_type sum;
  internal::NormP<Number, real_type> normp(val, p);
  internal::parallel_reduce (normp, vec_size, sum, thread_loop_partitioner);

  if (numbers::is_finite(sum) && sum >= std::numeric_limits<real_type>::min())
    return std::pow(sum, static_cast<real_type>(1./p));
  else
    {
      real_type scale = 0.;
      real_type sum = 1.;
      for (size_type i=0; i<vec_size; ++i)
        {
          if (val[i] != Number())
            {
              const real_type abs_x =
                numbers::NumberTraits<Number>::abs(val[i]);
              if (scale < abs_x)
                {
                  sum = 1. + sum * std::pow(scale/abs_x, p);
                  scale = abs_x;
                }
              else
                sum += std::pow(abs_x/scale, p);
            }
        }
      return scale * std::pow(sum, static_cast<real_type>(1./p));
    }
}



template <>
Vector<int>::real_type
Vector<int>::lp_norm (const real_type) const
{
  Assert(false, ExcMessage("No lp norm for integer vectors"));
  return -1;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::linfty_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type max = 0.;

  for (size_type i=0; i<vec_size; ++i)
    max = std::max (numbers::NumberTraits<Number>::abs(val[i]), max);

  return max;
}



template <typename Number>
Number
Vector<Number>::add_and_dot (const Number          a,
                             const Vector<Number> &V,
                             const Vector<Number> &W)
{
  Assert (vec_size!=0, ExcEmptyObject());
  AssertDimension (vec_size, V.size());
  AssertDimension (vec_size, W.size());

  Number sum;
  internal::AddAndDot<Number> adder(this->val, V.val, W.val, a);
  internal::parallel_reduce (adder, vec_size, sum, thread_loop_partitioner);
  AssertIsFinite(sum);

  return sum;
}



template <typename Number>
Vector<Number> &Vector<Number>::operator += (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  add (v);
  return *this;
}



template <typename Number>
Vector<Number> &Vector<Number>::operator -= (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_subtract_v<Number> vector_subtract(val, v.val);
  internal::parallel_for(vector_subtract,vec_size,thread_loop_partitioner);

  return *this;
}



template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_add_factor<Number> vector_add(val, v);
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::add (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_v<Number> vector_add(val, v.val);
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number> &v,
                          const Number b, const Vector<Number> &w)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_add_avpbw<Number> vector_add(val, v.val, w.val, a, b);
  internal::parallel_for(vector_add,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::sadd (const Number x,
                           const Vector<Number> &v)
{
  AssertIsFinite(x);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xv<Number> vector_sadd(val, v.val, x);
  internal::parallel_for(vector_sadd,vec_size,thread_loop_partitioner);
}



template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
                           const Vector<Number> &v, const Number b,
                           const Vector<Number> &w)
{
  AssertIsFinite(x);
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_sadd_xavbw<Number> vector_sadd(val, v.val, w.val, x,
                                                         a, b);
  internal::parallel_for(vector_sadd,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
                           const Vector<Number> &v, const Number b,
                           const Vector<Number> &w, const Number c,
                           const Vector<Number> &y)
{
  sadd (x, a, v, b, w);
  add (c, y);
}



template <typename Number>
void Vector<Number>::scale (const Vector<Number> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  internal::Vectorization_scale<Number> vector_scale(val, s.val);
  internal::parallel_for(vector_scale,vec_size,thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::scale (const Vector<Number2> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  for (size_type i=0; i<vec_size; ++i)
    val[i] *= Number(s.val[i]);
}



template <typename Number>
void Vector<Number>::equ (const Number a,
                          const Vector<Number> &u)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  internal::Vectorization_equ_au<Number> vector_equ(val, u.val, a);
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::equ (const Number a,
                          const Vector<Number2> &u)
{
  AssertIsFinite(a);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  // set the result vector to a*u. we have to
  // convert the elements of u to the type of
  // the result vector. this is necessary
  // because
  // operator*(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i=0; i<vec_size; ++i)
    val[i] = a * Number(u.val[i]);
}



template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                          const Number b, const Vector<Number> &v)
{
  AssertIsFinite(a);
  AssertIsFinite(b);

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_equ_aubv<Number> vector_equ(val, u.val, v.val, a, b);
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number> &u,
                          const Number b, const Vector<Number> &v,
                          const Number c, const Vector<Number> &w)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_equ_aubvcw<Number> vector_equ(val, u.val, v.val, w.val,
                                                        a, b, c);
  internal::parallel_for(vector_equ,vec_size,thread_loop_partitioner);
}


template <typename Number>
void Vector<Number>::ratio (const Vector<Number> &a,
                            const Vector<Number> &b)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (a.vec_size == b.vec_size,
          ExcDimensionMismatch (a.vec_size, b.vec_size));

  // no need to reinit with zeros, since
  // we overwrite them anyway
  reinit (a.size(), true);

  internal::Vectorization_ratio<Number> vector_ratio(val, a.val, b.val);
  internal::parallel_for(vector_ratio,vec_size,thread_loop_partitioner);
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const BlockVector<Number> &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);

  size_type this_index = 0;
  for (size_type b=0; b<v.n_blocks(); ++b)
    for (size_type i=0; i<v.block(b).size(); ++i, ++this_index)
      val[this_index] = v.block(b)(i);

  return *this;
}



#ifdef DEAL_II_WITH_PETSC

template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const PETScWrappers::Vector &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);
  if (vec_size != 0)
    {
      // get a representation of the vector
      // and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

      // restore the representation of the
      // vector
      ierr = VecRestoreArray (static_cast<const Vec &>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));
    }

  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const PETScWrappers::MPI::Vector &v)
{
  // do this in a two-stage process:
  // first convert to a sequential petsc
  // vector, then copy that
  PETScWrappers::Vector seq (v);
  *this = seq;

  return *this;
}

#endif


#ifdef DEAL_II_WITH_TRILINOS

template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const TrilinosWrappers::MPI::Vector &v)
{
  // Generate a localized version
  // of the Trilinos vectors and
  // then call the other =
  // operator.
  TrilinosWrappers::Vector localized_vector (v);
  *this = localized_vector;
  return *this;
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator= (const TrilinosWrappers::Vector &v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);
  if (vec_size != 0)
    {
      // get a representation of the vector
      // and copy it
      TrilinosScalar **start_ptr;
      int ierr = v.trilinos_vector().ExtractView (&start_ptr);
      AssertThrow (ierr == 0, ExcTrilinosError(ierr));

      std::copy (start_ptr[0], start_ptr[0]+vec_size, begin());
    }

  return *this;
}

#endif

template <typename Number>
template <typename Number2>
bool
Vector<Number>::operator== (const Vector<Number2> &v) const
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.size(), ExcDimensionMismatch(vec_size, v.size()));

  // compare the two vector. we have to
  // convert the elements of v to the type of
  // the result vector. this is necessary
  // because
  // operator==(complex<float>,complex<double>)
  // is not defined by default
  for (size_type i=0; i<vec_size; ++i)
    if (val[i] != Number(v.val[i]))
      return false;

  return true;
}



template <typename Number>
void Vector<Number>::print (const char *format) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (size_type j=0; j<size(); ++j)
    internal::print (val[j], format);
  std::printf ("\n");
}



template <typename Number>
void Vector<Number>::print (std::ostream      &out,
                            const unsigned int precision,
                            const bool         scientific,
                            const bool         across) const
{
  Assert (vec_size!=0, ExcEmptyObject());
  AssertThrow (out, ExcIO());

  std::ios::fmtflags old_flags = out.flags();
  unsigned int old_precision = out.precision (precision);

  out.precision (precision);
  if (scientific)
    out.setf (std::ios::scientific, std::ios::floatfield);
  else
    out.setf (std::ios::fixed, std::ios::floatfield);

  if (across)
    for (size_type i=0; i<size(); ++i)
      out << val[i] << ' ';
  else
    for (size_type i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;

  AssertThrow (out, ExcIO());
  // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
}



template <typename Number>
void
Vector<Number>::print (LogStream &out, const unsigned int width, const bool across) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (across)
    for (size_type i=0; i<size(); ++i)
      out << std::setw(width) << val[i] << ' ';
  else
    for (size_type i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;
}


template <typename Number>
void Vector<Number>::block_write (std::ostream &out) const
{
  AssertThrow (out, ExcIO());

  // other version of the following
  //  out << size() << std::endl << '[';
  // reason: operator<< seems to use
  // some resources that lead to
  // problems in a multithreaded
  // environment
  const size_type sz = size();
  char buf[16];

#ifdef DEAL_II_WITH_64BIT_INDICES
  std::sprintf(buf, "%llu", sz);
#else
  std::sprintf(buf, "%u", sz);
#endif
  std::strcat(buf, "\n[");

  out.write(buf, std::strlen(buf));
  out.write (reinterpret_cast<const char *>(begin()),
             reinterpret_cast<const char *>(end())
             - reinterpret_cast<const char *>(begin()));

  // out << ']';
  const char outro = ']';
  out.write (&outro, 1);

  AssertThrow (out, ExcIO());
}



template <typename Number>
void Vector<Number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  size_type sz;

  char buf[16];


  in.getline(buf,16,'\n');
  sz=std::atoi(buf);

  // fast initialization, since the
  // data elements are overwritten anyway
  reinit (sz, true);

  char c;
  //  in >> c;
  in.read (&c, 1);
  AssertThrow (c=='[', ExcIO());

  in.read (reinterpret_cast<char *>(begin()),
           reinterpret_cast<const char *>(end())
           - reinterpret_cast<const char *>(begin()));

  //  in >> c;
  in.read (&c, 1);
  AssertThrow (c==']', ExcIO());
}



template <typename Number>
IndexSet
Vector<Number>::locally_owned_elements() const
{
  return complete_index_set(size());
}



template <typename Number>
std::size_t
Vector<Number>::memory_consumption () const
{
  return sizeof(*this) + (max_vec_size * sizeof(Number));
}



template <typename Number>
void
Vector<Number>::allocate()
{
  // make sure that we don't create a memory leak
  Assert (val == 0, ExcInternalError());

  // then allocate memory with the proper alignment requirements of 64 bytes
  Utilities::System::posix_memalign ((void **)&val, 64, sizeof(Number)*max_vec_size);
}



template <typename Number>
void
Vector<Number>::deallocate()
{
  free(val);
  val = 0;
}

DEAL_II_NAMESPACE_CLOSE

#endif
