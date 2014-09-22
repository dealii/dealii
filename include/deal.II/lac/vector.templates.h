// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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

#ifndef __deal2__vector_templates_h
#define __deal2__vector_templates_h


#include <deal.II/base/template_constraints.h>
#include <deal.II/base/numbers.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/multithread_info.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/block_vector.h>

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
#include <mm_malloc.h>

DEAL_II_NAMESPACE_OPEN


namespace internal
{
  typedef types::global_dof_index size_type;

  template <typename T>
  bool is_non_negative (const T &t)
  {
    return t >= 0;
  }


  template <typename T>
  bool is_non_negative (const std::complex<T> &)
  {
    Assert (false,
            ExcMessage ("Complex numbers do not have an ordering."));

    return false;
  }


  template <typename T>
  void print (const T    &t,
              const char *format)
  {
    if (format != 0)
      std::printf (format, t);
    else
      std::printf (" %5.2f", double(t));
  }



  template <typename T>
  void print (const std::complex<T> &t,
              const char            *format)
  {
    if (format != 0)
      std::printf (format, t.real(), t.imag());
    else
      std::printf (" %5.2f+%5.2fi",
                   double(t.real()), double(t.imag()));
  }

  // call std::copy, except for in
  // the case where we want to copy
  // from std::complex to a
  // non-complex type
  template <typename T, typename U>
  void copy (const T *begin,
             const T *end,
             U       *dest)
  {
    std::copy (begin, end, dest);
  }

  template <typename T, typename U>
  void copy (const std::complex<T> *begin,
             const std::complex<T> *end,
             std::complex<U>       *dest)
  {
    std::copy (begin, end, dest);
  }

  template <typename T, typename U>
  void copy (const std::complex<T> *,
             const std::complex<T> *,
             U *)
  {
    Assert (false, ExcMessage ("Can't convert a vector of complex numbers "
                               "into a vector of reals/doubles"));
  }

  template <typename Functor>
  void vectorized_transform(Functor &functor,
                            size_type vec_size)
  {
#ifndef DEAL_II_WITH_THREADS
    functor(0,vec_size);
#else
    if (vec_size>internal::Vector::minimum_parallel_grain_size)
      {
        tbb::parallel_for (tbb::blocked_range<size_type> (0,
                                                          vec_size,
                                                          internal::Vector::minimum_parallel_grain_size),
                           functor,
                           tbb::auto_partitioner());
      }
    else if (vec_size > 0)
      functor(0,vec_size);
#endif
  }


  // Define the functors neccessary to use SIMD with TBB.
  template <typename Number>
  struct Vectorization_multiply_factor
  {
    Number *val;
    Number factor;

#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] *= factor;
    }
  };

  template <typename Number>
  struct Vectorization_add_av
  {
    Number *val;
    Number *v_val;
    Number factor;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] += factor*v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xav
  {
    Number *val;
    Number *v_val;
    Number a;
    Number x;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = x*val[i] + a*v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_subtract_v
  {
    Number *val;
    Number *v_val;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] -= v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_add_factor
  {
    Number *val;
    Number factor;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] += factor;
    }
  };

  template <typename Number>
  struct Vectorization_add_v
  {
    Number *val;
    Number *v_val;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] += v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_add_avpbw
  {
    Number *val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = val[i] + a*v_val[i] + b*w_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xv
  {
    Number *val;
    Number *v_val;
    Number x;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = x*val[i] + v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_sadd_xavbw
  {
    Number *val;
    Number *v_val;
    Number *w_val;
    Number x;
    Number a;
    Number b;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = x*val[i] + a*v_val[i] + b*w_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_scale
  {
    Number *val;
    Number *v_val;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] *= v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_equ_au
  {
    Number *val;
    Number *u_val;
    Number a;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = a*u_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_equ_aubv
  {
    Number *val;
    Number *u_val;
    Number *v_val;
    Number a;
    Number b;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = a*u_val[i] + b*v_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_equ_aubvcw
  {
    Number *val;
    Number *u_val;
    Number *v_val;
    Number *w_val;
    Number a;
    Number b;
    Number c;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = a*u_val[i] + b*v_val[i] + c*w_val[i];
    }
  };

  template <typename Number>
  struct Vectorization_ratio
  {
    Number *val;
    Number *a_val;
    Number *b_val;
#ifdef DEAL_II_WITH_THREADS
    void operator() (const tbb::blocked_range<size_type> &range) const
    {
      operator()(range.begin(),range.end());
    }
#endif

    void operator() (const size_type begin, const size_type end) const
    {
      DEAL_II_OPENMP_SIMD_PRAGMA
      for (size_type i=begin; i<end; ++i)
        val[i] = a_val[i]/b_val[i];
    }
  };

}




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
      allocate(max_vec_size);
      Assert (val != 0, ExcOutOfMemory());
      *this = v;
    }
}


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
      allocate(max_vec_size);
      Assert (val != 0, ExcOutOfMemory());
      std::copy (v.begin(), v.end(), begin());
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
      allocate(max_vec_size);
      Assert (val != 0, ExcOutOfMemory());

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
      allocate(max_vec_size);
      Assert (val != 0, ExcOutOfMemory());

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
      allocate(max_vec_size);
      Assert (val != 0, ExcOutOfMemory());

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
template <typename Number2>
void Vector<Number>::reinit (const Vector<Number2> &v,
                             const bool fast)
{
  reinit (v.size(), fast);
}

// Moved to vector.h as an inline function by Luca Heltai on
// 2009/04/12 to prevent strange compiling errors, after making swap
// virtual.
// template <typename Number>
// void
// Vector<Number>::swap (Vector<Number> &v)
// {
//   std::swap (vec_size,     v.vec_size);
//   std::swap (max_vec_size, v.max_vec_size);
//   std::swap (val,          v.val);
// }



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



namespace internal
{
  namespace Vector
  {
    template <typename T>
    void set_subrange (const T            s,
                       const typename dealii::Vector<T>::size_type begin,
                       const typename dealii::Vector<T>::size_type end,
                       dealii::Vector<T> &dst)
    {
      if (s == T())
        std::memset ((dst.begin()+begin),0,(end-begin)*sizeof(T));
      else
        std::fill (&*(dst.begin()+begin), &*(dst.begin()+end), s);
    }


    template <typename T>
    void copy_subrange (const typename dealii::Vector<T>::size_type         begin,
                        const typename dealii::Vector<T>::size_type         end,
                        const dealii::Vector<T> &src,
                        dealii::Vector<T>       &dst)
    {
      memcpy(&*(dst.begin()+begin), &*(src.begin()+begin),
             (end-begin)*sizeof(T));
    }


    template <typename T, typename U>
    void copy_subrange (const typename dealii::Vector<T>::size_type         begin,
                        const typename dealii::Vector<T>::size_type         end,
                        const dealii::Vector<T> &src,
                        dealii::Vector<U>       &dst)
    {
      const T *q = src.begin()+begin;
      const T *const end_q = src.begin()+end;
      U *p = dst.begin()+begin;
      for (; q!=end_q; ++q, ++p)
        *p = *q;
    }


    template <typename T, typename U>
    void copy_vector (const dealii::Vector<T> &src,
                      dealii::Vector<U>       &dst)
    {
      if (PointerComparison::equal(&src, &dst))
        return;

      const typename dealii::Vector<T>::size_type vec_size = src.size();
      const typename dealii::Vector<U>::size_type dst_size = dst.size();
      if (dst_size != vec_size)
        dst.reinit (vec_size, true);
      if (vec_size>internal::Vector::minimum_parallel_grain_size)
        parallel::apply_to_subranges (0U, vec_size,
                                      std_cxx11::bind(&internal::Vector::template
                                                      copy_subrange <T,U>,
                                                      std_cxx11::_1,
                                                      std_cxx11::_2,
                                                      std_cxx11::cref(src),
                                                      std_cxx11::ref(dst)),
                                      internal::Vector::minimum_parallel_grain_size);
      else if (vec_size > 0)
        copy_subrange (0U, vec_size, src, dst);
    }
  }
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator = (const Number s)
{
  Assert (numbers::is_finite(s), ExcNumberNotFinite());
  if (s != Number())
    Assert (vec_size!=0, ExcEmptyObject());
  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::apply_to_subranges (0U, vec_size,
                                  std_cxx11::bind(&internal::Vector::template
                                                  set_subrange<Number>,
                                                  s, std_cxx11::_1, std_cxx11::_2, std_cxx11::ref(*this)),
                                  internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    internal::Vector::set_subrange<Number>(s, 0U, vec_size, *this);

  return *this;
}



#ifdef DEAL_II_BOOST_BIND_COMPILER_BUG
template <>
Vector<std::complex<float> > &
Vector<std::complex<float> >::operator = (const std::complex<float> s)
{
  Assert (numbers::is_finite(s), ExcNumberNotFinite());
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
  Assert (numbers::is_finite(factor),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_multiply_factor<Number> vector_multiply;
  vector_multiply.val = val;
  vector_multiply.factor = factor;

  internal::vectorized_transform(vector_multiply,vec_size);

  return *this;
}



template <typename Number>
void
Vector<Number>::add (const Number a,
                     const Vector<Number> &v)
{
  Assert (numbers::is_finite(a),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_av<Number> vector_add_av;
  vector_add_av.val = val;
  vector_add_av.v_val = v.val;
  vector_add_av.factor = a;
  internal::vectorized_transform(vector_add_av,vec_size);
}



template <typename Number>
void
Vector<Number>::sadd (const Number x,
                      const Number a,
                      const Vector<Number> &v)
{
  Assert (numbers::is_finite(x),ExcNumberNotFinite());
  Assert (numbers::is_finite(a),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xav<Number> vector_sadd_xav;
  vector_sadd_xav.val = val;
  vector_sadd_xav.v_val = v.val;
  vector_sadd_xav.a = a;
  vector_sadd_xav.x = x;
  internal::vectorized_transform(vector_sadd_xav,vec_size);
}



namespace internal
{
  namespace Vector
  {
    // All sums over all the vector entries (l2-norm, inner product, etc.) are
    // performed with the same code, using a templated operation defined here
    template <typename Number, typename Number2>
    struct InnerProd
    {
      Number
      operator() (const Number *&X, const Number2 *&Y, const Number &) const
      {
        return *X++ * Number(numbers::NumberTraits<Number2>::conjugate(*Y++));
      }
    };

    template <typename Number, typename RealType>
    struct Norm2
    {
      RealType
      operator() (const Number  *&X, const Number  *&, const RealType &) const
      {
        return numbers::NumberTraits<Number>::abs_square(*X++);
      }
    };

    template <typename Number, typename RealType>
    struct Norm1
    {
      RealType
      operator() (const Number  *&X, const Number  *&, const RealType &) const
      {
        return numbers::NumberTraits<Number>::abs(*X++);
      }
    };

    template <typename Number, typename RealType>
    struct NormP
    {
      RealType
      operator() (const Number  *&X, const Number  *&, const RealType &p) const
      {
        return std::pow(numbers::NumberTraits<Number>::abs(*X++), p);
      }
    };

    template <typename Number>
    struct MeanValue
    {
      Number
      operator() (const Number  *&X, const Number  *&, const Number &) const
      {
        return *X++;
      }
    };

    // this is the main working loop for all vector sums using the templated
    // operation above. it accumulates the sums using a block-wise summation
    // algorithm with post-update. this blocked algorithm has been proposed in
    // a similar form by Castaldo, Whaley and Chronopoulos (SIAM
    // J. Sci. Comput. 31, 1156-1174, 2008) and we use the smallest possible
    // block size, 2. Sometimes it is referred to as pairwise summation. The
    // worst case error made by this algorithm is on the order O(eps *
    // log2(vec_size)), whereas a naive summation is O(eps * vec_size). Even
    // though the Kahan summation is even more accurate with an error O(eps)
    // by carrying along remainders not captured by the main sum, that involves
    // additional costs which are not worthwhile. See the Wikipedia article on
    // the Kahan summation algorithm.

    // The algorithm implemented here has the additional benefit that it is
    // easily parallelized without changing the order of how the elements are
    // added (floating point addition is not associative). For the same vector
    // size and minimum_parallel_grainsize, the blocks are always the
    // same and added pairwise. At the innermost level, eight values are added
    // consecutively in order to better balance multiplications and additions.

    // The code returns the result as the last argument in order to make
    // spawning tasks simpler and use automatic template deduction.
    template <typename Operation, typename Number, typename Number2,
              typename ResultType, typename size_type>
    void accumulate (const Operation   &op,
                     const Number      *X,
                     const Number2     *Y,
                     const ResultType   power,
                     const size_type    vec_size,
                     ResultType        &result,
                     const int          depth = -1)
    {
      if (vec_size <= 4096)
        {
          // the vector is short enough so we perform the summation. first
          // work on the regular part. The innermost 32 values are expanded in
          // order to obtain known loop bounds for most of the work.
          const Number *X_original = X;
          ResultType outer_results [128];
          size_type n_chunks = vec_size / 32;
          const size_type remainder = vec_size % 32;
          Assert (remainder == 0 || n_chunks < 128, ExcInternalError());

          for (size_type i=0; i<n_chunks; ++i)
            {
              ResultType r0 = op(X, Y, power);
              for (size_type j=1; j<8; ++j)
                r0 += op(X, Y, power);
              ResultType r1 = op(X, Y, power);
              for (size_type j=1; j<8; ++j)
                r1 += op(X, Y, power);
              r0 += r1;
              r1 = op(X, Y, power);
              for (size_type j=1; j<8; ++j)
                r1 += op(X, Y, power);
              ResultType r2 = op(X, Y, power);
              for (size_type j=1; j<8; ++j)
                r2 += op(X, Y, power);
              r1 += r2;
              r0 += r1;
              outer_results[i] = r0;
            }

          // now work on the remainder, i.e., the last
          // up to 32 values. Use switch statement with
          // fall-through to work on these values.
          if (remainder > 0)
            {
              const size_type inner_chunks = remainder / 8;
              Assert (inner_chunks <= 3, ExcInternalError());
              const size_type remainder_inner = remainder % 8;
              ResultType r0 = ResultType(), r1 = ResultType(),
                         r2 = ResultType();
              switch (inner_chunks)
                {
                case 3:
                  r2 = op(X, Y, power);
                  for (size_type j=1; j<8; ++j)
                    r2 += op(X, Y, power);
                // no break
                case 2:
                  r1 = op(X, Y, power);
                  for (size_type j=1; j<8; ++j)
                    r1 += op(X, Y, power);
                  r1 += r2;
                // no break
                case 1:
                  r2 = op(X, Y, power);
                  for (size_type j=1; j<8; ++j)
                    r2 += op(X, Y, power);
                // no break
                default:
                  for (size_type j=0; j<remainder_inner; ++j)
                    r0 += op(X, Y, power);
                  r0 += r2;
                  r0 += r1;
                  outer_results[n_chunks] = r0;
                  break;
                }
              n_chunks++;
            }
          AssertDimension(static_cast<size_type> (X - X_original), vec_size);

          // now sum the results from the chunks
          // recursively
          while (n_chunks > 1)
            {
              if (n_chunks % 2 == 1)
                outer_results[n_chunks++] = ResultType();
              for (size_type i=0; i<n_chunks; i+=2)
                outer_results[i/2] = outer_results[i] + outer_results[i+1];
              n_chunks /= 2;
            }
          result = outer_results[0];
        }
#ifdef DEAL_II_WITH_THREADS
      else if (multithread_info.n_threads() > 1 &&
               vec_size > 4 * internal::Vector::minimum_parallel_grain_size &&
               depth != 0)
        {
          // split the vector into smaller pieces to be worked on recursively
          // and create tasks for them. Make pieces divisible by 1024.
          const size_type new_size = (vec_size / 4096) * 1024;
          ResultType r0, r1, r2, r3;

          // find out how many recursions we should make (avoid too deep
          // hierarchies of tasks on large vectors), max use 8 *
          // multithread_info.n_threads()
          int next_depth = depth;
          if (depth == -1)
            next_depth = 8 * multithread_info.n_threads();
          next_depth /= 4;

          Threads::TaskGroup<> task_group;
          task_group += Threads::new_task(&accumulate<Operation,Number,Number2,
                                          ResultType,size_type>,
                                          op, X, Y, power, new_size, r0, next_depth);
          task_group += Threads::new_task(&accumulate<Operation,Number,Number2,
                                          ResultType,size_type>,
                                          op, X+new_size, Y+new_size, power,
                                          new_size, r1, next_depth);
          task_group += Threads::new_task(&accumulate<Operation,Number,Number2,
                                          ResultType,size_type>,
                                          op, X+2*new_size, Y+2*new_size, power,
                                          new_size, r2, next_depth);
          task_group += Threads::new_task(&accumulate<Operation,Number,Number2,
                                          ResultType,size_type>,
                                          op, X+3*new_size, Y+3*new_size, power,
                                          vec_size-3*new_size, r3, next_depth);
          task_group.join_all();
          r0 += r1;
          r2 += r3;
          result = r0 + r2;
        }
#endif
      else
        {
          // split vector into four pieces and work on
          // the pieces recursively. Make pieces (except last)
          // divisible by 1024.
          const size_type new_size = (vec_size / 4096) * 1024;
          ResultType r0, r1, r2, r3;
          accumulate (op, X, Y, power, new_size, r0);
          accumulate (op, X+new_size, Y+new_size, power, new_size, r1);
          accumulate (op, X+2*new_size, Y+2*new_size, power, new_size, r2);
          accumulate (op, X+3*new_size, Y+3*new_size, power, vec_size-3*new_size, r3);
          r0 += r1;
          r2 += r3;
          result = r0 + r2;
        }
    }
  }
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
  internal::Vector::accumulate (internal::Vector::InnerProd<Number,Number2>(),
                                val, v.val, Number(), vec_size, sum);
  Assert(numbers::is_finite(sum), ExcNumberNotFinite());

  return sum;
}


template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::norm_sqr () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Vector::accumulate (internal::Vector::Norm2<Number,real_type>(),
                                val, val, real_type(), vec_size, sum);

  Assert(numbers::is_finite(sum), ExcNumberNotFinite());

  return sum;
}


template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  Number sum;
  internal::Vector::accumulate (internal::Vector::MeanValue<Number>(),
                                val, val, Number(), vec_size, sum);

  return sum / real_type(size());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l1_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type sum;
  internal::Vector::accumulate (internal::Vector::Norm1<Number,real_type>(),
                                val, val, real_type(), vec_size, sum);

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
  internal::Vector::accumulate (internal::Vector::Norm2<Number,real_type>(),
                                val, val, real_type(), vec_size, norm_square);
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
      Assert(numbers::is_finite(scale)*std::sqrt(sum), ExcNumberNotFinite());
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
  internal::Vector::accumulate (internal::Vector::NormP<Number,real_type>(),
                                val, val, p, vec_size, sum);

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

  internal::Vectorization_subtract_v<Number> vector_subtract;
  vector_subtract.val = val;
  vector_subtract.v_val = v.val;
  internal::vectorized_transform(vector_subtract,vec_size);

  return *this;
}


template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  internal::Vectorization_add_factor<Number> vector_add;
  vector_add.val = val;
  vector_add.factor = v;
  internal::vectorized_transform(vector_add,vec_size);
}


template <typename Number>
void Vector<Number>::add (const Vector<Number> &v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_add_v<Number> vector_add;
  vector_add.val = val;
  vector_add.v_val = v.val;
  internal::vectorized_transform(vector_add,vec_size);
}


template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number> &v,
                          const Number b, const Vector<Number> &w)
{
  Assert (numbers::is_finite(a),ExcNumberNotFinite());
  Assert (numbers::is_finite(b),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_add_avpbw<Number> vector_add;
  vector_add.val = val;
  vector_add.v_val = v.val;
  vector_add.w_val = w.val;
  vector_add.a = a;
  vector_add.b = b;
  internal::vectorized_transform(vector_add,vec_size);
}



template <typename Number>
void Vector<Number>::sadd (const Number x,
                           const Vector<Number> &v)
{
  Assert (numbers::is_finite(x),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_sadd_xv<Number> vector_sadd;
  vector_sadd.val = val;
  vector_sadd.v_val = v.val;
  vector_sadd.x = x;
  internal::vectorized_transform(vector_sadd,vec_size);
}



template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
                           const Vector<Number> &v, const Number b,
                           const Vector<Number> &w)
{
  Assert (numbers::is_finite(x),ExcNumberNotFinite());
  Assert (numbers::is_finite(a),ExcNumberNotFinite());
  Assert (numbers::is_finite(b),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  internal::Vectorization_sadd_xavbw<Number> vector_sadd;
  vector_sadd.val = val;
  vector_sadd.v_val = v.val;
  vector_sadd.w_val = w.val;
  vector_sadd.x = x;
  vector_sadd.a = a;
  vector_sadd.b = b;
  internal::vectorized_transform(vector_sadd,vec_size);
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

  internal::Vectorization_scale<Number> vector_scale;
  vector_scale.val = val;
  vector_scale.v_val = s.val;
  internal::vectorized_transform(vector_scale,vec_size);
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
  Assert (numbers::is_finite(a), ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  internal::Vectorization_equ_au<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.a = a;
  internal::vectorized_transform(vector_equ,vec_size);
}



template <typename Number>
template <typename Number2>
void Vector<Number>::equ (const Number a,
                          const Vector<Number2> &u)
{
  Assert (numbers::is_finite(a), ExcNumberNotFinite());

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
  Assert (numbers::is_finite(a),ExcNumberNotFinite());
  Assert (numbers::is_finite(b),ExcNumberNotFinite());

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  internal::Vectorization_equ_aubv<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.v_val = v.val;
  vector_equ.a = a;
  vector_equ.b = b;
  internal::vectorized_transform(vector_equ,vec_size);
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

  internal::Vectorization_equ_aubvcw<Number> vector_equ;
  vector_equ.val = val;
  vector_equ.u_val = u.val;
  vector_equ.v_val = v.val;
  vector_equ.w_val = w.val;
  vector_equ.a = a;
  vector_equ.b = b;
  vector_equ.c = c;
  internal::vectorized_transform(vector_equ,vec_size);
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

  internal::Vectorization_ratio<Number> vector_ratio;
  vector_ratio.val = val;
  vector_ratio.a_val = a.val;
  vector_ratio.b_val = b.val;
  internal::vectorized_transform(vector_ratio,vec_size);
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator = (const BlockVector<Number> &v)
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
Vector<Number>::operator = (const PETScWrappers::Vector &v)
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
Vector<Number>::operator = (const PETScWrappers::MPI::Vector &v)
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
Vector<Number>::operator = (const TrilinosWrappers::MPI::Vector &v)
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
Vector<Number>::operator = (const TrilinosWrappers::Vector &v)
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
Vector<Number>::operator == (const Vector<Number2> &v) const
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
Vector<Number>::allocate(const size_type size)
{
  val = static_cast<Number *>(_mm_malloc (sizeof(Number)*size, 64));
}



template <typename Number>
void
Vector<Number>::deallocate()
{
  _mm_free(val);
}

DEAL_II_NAMESPACE_CLOSE

#endif
