//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__vector_templates_h
#define __deal2__vector_templates_h


#include <base/template_constraints.h>
#include <base/numbers.h>
#include <lac/vector.h>
#include <lac/block_vector.h>

#ifdef DEAL_II_USE_PETSC
#  include <lac/petsc_vector.h>
#  include <lac/petsc_parallel_vector.h>
#endif

#ifdef DEAL_II_USE_TRILINOS
#  include <lac/trilinos_vector.h>
#endif

#include <cmath>
#include <cstring>
#include <algorithm>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

#define BLOCK_LEVEL 6

namespace internal
{
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
	     U                     *)
  {
    Assert (false, ExcMessage ("Can't convert a vector of complex numbers "
			       "into a vector of reals/doubles"));
  }
}




template <typename Number>
Vector<Number>::Vector (const Vector<Number>& v)
                :
		Subscriptor(),
		vec_size(v.size()),
		max_vec_size(v.size()),
		val(0)
{
  if (vec_size != 0)
    {
      val = new Number[max_vec_size];
      Assert (val != 0, ExcOutOfMemory());
      std::copy (v.begin(), v.end(), begin());
    }
}


#ifndef DEAL_II_EXPLICIT_CONSTRUCTOR_BUG

template <typename Number>
template <typename OtherNumber>
Vector<Number>::Vector (const Vector<OtherNumber>& v)
                :
		Subscriptor(),
		vec_size(v.size()),
		max_vec_size(v.size()),
		val(0)
{
  if (vec_size != 0)
    {
      val = new Number[max_vec_size];
      Assert (val != 0, ExcOutOfMemory());
      std::copy (v.begin(), v.end(), begin());
    }
}

#endif

#ifdef DEAL_II_USE_PETSC


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
      val = new Number[max_vec_size];
      Assert (val != 0, ExcOutOfMemory());

                                       // get a representation of the vector
                                       // and copy it
      PetscScalar *start_ptr;
      int ierr = VecGetArray (static_cast<const Vec&>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

                                       // restore the representation of the
                                       // vector
      ierr = VecRestoreArray (static_cast<const Vec&>(v), &start_ptr);
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

#ifdef DEAL_II_USE_TRILINOS

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
      val = new Number[max_vec_size];
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
      val = new Number[max_vec_size];
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
void Vector<Number>::reinit (const Vector<Number2>& v,
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

  for (unsigned int i=0; i<vec_size; ++i)
    if (val[i] != Number(0))
      return false;
  return true;
}



template <typename Number>
bool
Vector<Number>::is_non_negative () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (unsigned int i=0; i<vec_size; ++i)
    if ( ! internal::is_non_negative (val[i]))
      return false;

  return true;
}



template <typename Number>
template <typename Number2>
Number Vector<Number>::operator * (const Vector<Number2>& v) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (PointerComparison::equal (this, &v))
    return norm_sqr();

  Assert (vec_size == v.size(),
	  ExcDimensionMismatch(vec_size, v.size()));

  const unsigned int blocking = 1<<BLOCK_LEVEL;
  register Number sum1, sum2, sum3, sum = Number();
  const Number * X = val, *X_end = X + vec_size, 
    *X_end3 = X + ((vec_size>>(BLOCK_LEVEL))<<(BLOCK_LEVEL)),
    *X_end2 = X + ((vec_size>>(2*BLOCK_LEVEL))<<(2*BLOCK_LEVEL)),
    *X_end1 = X + ((vec_size>>(3*BLOCK_LEVEL))<<(3*BLOCK_LEVEL));
  const Number2  * Y = v.val;

				   // multiply the two vectors. we have to
				   // convert the elements of u to the type of
				   // the result vector. this is necessary
				   // because
				   // operator*(complex<float>,complex<double>)
				   // is not defined by default. do the 
				   // operations block-wise with post-update.
				   // use three nested loops, which should make 
				   // the roundoff error very small up to 
				   // about 20m entries. in the
				   // end do extra loops with the remainders. 
				   // this blocked algorithm has been proposed
				   // by Castaldo, Whaley and Chronopoulos
				   // (SIAM J. Sci. Comput. 31, 1156-1174, 
				   // 2008)
  while (X != X_end1)
    {
      sum1 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum2 = 0;
	  for (unsigned int k=0; k<blocking; ++k)
	    {
	      sum3 = 0;
	      for (unsigned int i=0; i<blocking; ++i)
		sum3 += *X++ * Number(numbers::NumberTraits<Number2>::conjugate(*Y++));
	      sum2 += sum3;
	    }
	  sum1 += sum2;
	}
      sum += sum1;
    }
  while (X != X_end2)
    {
      sum2 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum3 = 0;
	  for (unsigned int i=0; i<blocking; ++i)
	    sum3 += *X++ * Number(numbers::NumberTraits<Number2>::conjugate(*Y++));
	  sum2 += sum3;
	}
      sum += sum2;
    }
  while (X != X_end3)
    {
      sum3 = 0;
      for (unsigned int i=0; i<blocking; ++i)
	sum3 += *X++ * Number(numbers::NumberTraits<Number2>::conjugate(*Y++));
      sum += sum3;
    }

  sum3 = 0;
  while (X != X_end) 
    sum3 += *X++ * Number(numbers::NumberTraits<Number2>::conjugate(*Y++));
  sum += sum3;

  return sum;
}


template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::norm_sqr () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  const unsigned int blocking = 1<<BLOCK_LEVEL;
  register real_type sum1, sum2, sum3, sum = 0.;
  const Number * X = val, *X_end = X + vec_size, 
    *X_end3 = X + ((vec_size>>(BLOCK_LEVEL))<<(BLOCK_LEVEL)),
    *X_end2 = X + ((vec_size>>(2*BLOCK_LEVEL))<<(2*BLOCK_LEVEL)),
    *X_end1 = X + ((vec_size>>(3*BLOCK_LEVEL))<<(3*BLOCK_LEVEL));

  while (X != X_end1)
    {
      sum1 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum2 = 0;
	  for (unsigned int k=0; k<blocking; ++k)
	    {
	      sum3 = 0;
	      for (unsigned int i=0; i<blocking; ++i)
		sum3 += numbers::NumberTraits<Number>::abs_square(*X++);
	      sum2 += sum3;
	    }
	  sum1 += sum2;
	}
      sum += sum1;
    }
  while (X != X_end2)
    {
      sum2 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum3 = 0;
	  for (unsigned int i=0; i<blocking; ++i)
	    sum3 += numbers::NumberTraits<Number>::abs_square(*X++);
	  sum2 += sum3;
	}
      sum += sum2;
    }
  while (X != X_end3)
    {
      sum3 = 0;
      for (unsigned int i=0; i<blocking; ++i)
	sum3 += numbers::NumberTraits<Number>::abs_square(*X++);
      sum += sum3;
    }

  sum3 = 0;
  while (X != X_end) 
    sum3 += numbers::NumberTraits<Number>::abs_square(*X++);
  sum += sum3;

  return sum;
}


template <typename Number>
Number Vector<Number>::mean_value () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  const unsigned int blocking = 1<<BLOCK_LEVEL;
  register Number sum1, sum2, sum3, sum = 0.;
  const Number * X = val, *X_end = X + vec_size, 
    *X_end3 = X + ((vec_size>>(BLOCK_LEVEL))<<(BLOCK_LEVEL)),
    *X_end2 = X + ((vec_size>>(2*BLOCK_LEVEL))<<(2*BLOCK_LEVEL)),
    *X_end1 = X + ((vec_size>>(3*BLOCK_LEVEL))<<(3*BLOCK_LEVEL));

  while (X != X_end1)
    {
      sum1 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum2 = 0;
	  for (unsigned int k=0; k<blocking; ++k)
	    {
	      sum3 = 0;
	      for (unsigned int i=0; i<blocking; ++i)
		sum3 += *X++;
	      sum2 += sum3;
	    }
	  sum1 += sum2;
	}
      sum += sum1;
    }
  while (X != X_end2)
    {
      sum2 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum3 = 0;
	  for (unsigned int i=0; i<blocking; ++i)
	    sum3 += *X++;
	  sum2 += sum3;
	}
      sum += sum2;
    }
  while (X != X_end3)
    {
      sum3 = 0;
      for (unsigned int i=0; i<blocking; ++i)
	sum3 += *X++;
      sum += sum3;
    }

  sum3 = 0;
  while (X != X_end) 
    sum3 += *X++;
  sum += sum3;

  return sum / real_type(size());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l1_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());
 
  const unsigned int blocking = 1<<BLOCK_LEVEL;
  register real_type sum1, sum2, sum3, sum = 0.;
  const Number * X = val, *X_end = X + vec_size, 
    *X_end3 = X + ((vec_size>>(BLOCK_LEVEL))<<(BLOCK_LEVEL)),
    *X_end2 = X + ((vec_size>>(2*BLOCK_LEVEL))<<(2*BLOCK_LEVEL)),
    *X_end1 = X + ((vec_size>>(3*BLOCK_LEVEL))<<(3*BLOCK_LEVEL));

  while (X != X_end1)
    {
      sum1 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum2 = 0;
	  for (unsigned int k=0; k<blocking; ++k)
	    {
	      sum3 = 0;
	      for (unsigned int i=0; i<blocking; ++i)
		sum3 += numbers::NumberTraits<Number>::abs(*X++);
	      sum2 += sum3;
	    }
	  sum1 += sum2;
	}
      sum += sum1;
    }
  while (X != X_end2)
    {
      sum2 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum3 = 0;
	  for (unsigned int i=0; i<blocking; ++i)
	    sum3 += numbers::NumberTraits<Number>::abs(*X++);
	  sum2 += sum3;
	}
      sum += sum2;
    }
  while (X != X_end3)
    {
      sum3 = 0;
      for (unsigned int i=0; i<blocking; ++i)
	sum3 += numbers::NumberTraits<Number>::abs(*X++);
      sum += sum3;
    }

  sum3 = 0;
  while (X != X_end) 
    sum3 += numbers::NumberTraits<Number>::abs(*X++);
  sum += sum3;

  return sum;
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::l2_norm () const
{
  return std::sqrt(norm_sqr());
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::lp_norm (const real_type p) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (p == 1.)
    return l1_norm();
  else if (p == 2.)
    return std::sqrt(norm_sqr());

  const unsigned int blocking = 1<<BLOCK_LEVEL;
  register real_type sum1, sum2, sum3, sum = 0.;
  const Number * X = val, *X_end = X + vec_size, 
    *X_end3 = X + ((vec_size>>(BLOCK_LEVEL))<<(BLOCK_LEVEL)),
    *X_end2 = X + ((vec_size>>(2*BLOCK_LEVEL))<<(2*BLOCK_LEVEL)),
    *X_end1 = X + ((vec_size>>(3*BLOCK_LEVEL))<<(3*BLOCK_LEVEL));

  while (X != X_end1)
    {
      sum1 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum2 = 0;
	  for (unsigned int k=0; k<blocking; ++k)
	    {
	      sum3 = 0;
	      for (unsigned int i=0; i<blocking; ++i)
		sum3 += std::pow(numbers::NumberTraits<Number>::abs(*X++),p);
	      sum2 += sum3;
	    }
	  sum1 += sum2;
	}
      sum += sum1;
    }
  while (X != X_end2)
    {
      sum2 = 0.;
      for (unsigned int j=0; j<blocking; ++j)
	{
	  sum3 = 0;
	  for (unsigned int i=0; i<blocking; ++i)
	    sum3 += std::pow(numbers::NumberTraits<Number>::abs(*X++),p);
	  sum2 += sum3;
	}
      sum += sum2;
    }
  while (X != X_end3)
    {
      sum3 = 0;
      for (unsigned int i=0; i<blocking; ++i)
	sum3 += std::pow(numbers::NumberTraits<Number>::abs(*X++),p);
      sum += sum3;
    }

  sum3 = 0;
  while (X != X_end) 
    sum3 += std::pow(numbers::NumberTraits<Number>::abs(*X++),p);
  sum += sum3;

  return std::pow(sum, static_cast<real_type>(1./p));
}



template <typename Number>
typename Vector<Number>::real_type
Vector<Number>::linfty_norm () const
{
  Assert (vec_size!=0, ExcEmptyObject());

  real_type max = 0;

  for (unsigned int i=0; i<vec_size; ++i)
    max = std::max (numbers::NumberTraits<Number>::abs(val[i]), max);

  return max;
}


template <typename Number>
Vector<Number>& Vector<Number>::operator += (const Vector<Number>& v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  add (v);
  return *this;
}


template <typename Number>
Vector<Number>& Vector<Number>::operator -= (const Vector<Number>& v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 v.val,
			 val,
			 (boost::lambda::_1 - boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] -= v.val[i];

  return *this;
}


template <typename Number>
void Vector<Number>::add (const Number v)
{
  Assert (vec_size!=0, ExcEmptyObject());

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 val,
			 (boost::lambda::_1 + v),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] += v;
}


template <typename Number>
void Vector<Number>::add (const Vector<Number>& v)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 v.val,
			 val,
			 (boost::lambda::_1 + boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] += v.val[i];
}


template <typename Number>
void Vector<Number>::add (const Number a, const Vector<Number>& v,
			  const Number b, const Vector<Number>& w)
{
  Assert (numbers::is_finite(a),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
  Assert (numbers::is_finite(b),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 v.val,
			 w.val,
			 val,
			 (boost::lambda::_1 + a*boost::lambda::_2 + b*boost::lambda::_3),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] += a * v.val[i] + b * w.val[i];
}


template <typename Number>
void Vector<Number>::sadd (const Number x,
			   const Vector<Number>& v)
{
  Assert (numbers::is_finite(x),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 v.val,
			 val,
			 (x*boost::lambda::_1 + boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = x * val[i] + v.val[i];
}



template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
			   const Vector<Number>& v, const Number b,
                           const Vector<Number>& w)
{
  Assert (numbers::is_finite(x),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
  Assert (numbers::is_finite(a),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
  Assert (numbers::is_finite(b),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
  parallel::transform (val,
		       val+vec_size,
		       v.val,
		       w.val,
		       val,
		       (x*boost::lambda::_1 + a*boost::lambda::_2 + b*boost::lambda::_3),
		       internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = x*val[i] + a * v.val[i] + b * w.val[i];
}


template <typename Number>
void Vector<Number>::sadd (const Number x, const Number a,
			   const Vector<Number>& v, const Number b,
			   const Vector<Number>& w, const Number c,
                           const Vector<Number>& y)
{
  sadd (x, a, v, b, w);
  add (c, y);
}



template <typename Number>
void Vector<Number>::scale (const Vector<Number> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (val,
			 val+vec_size,
			 s.val,
			 val,
			 (boost::lambda::_1*boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] *= s.val[i];
}



template <typename Number>
template <typename Number2>
void Vector<Number>::scale (const Vector<Number2> &s)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == s.vec_size, ExcDimensionMismatch(vec_size, s.vec_size));

  for (unsigned int i=0; i<vec_size; ++i)
    val[i] *= s.val[i];
}



template <typename Number>
void Vector<Number>::equ (const Number a,
			  const Vector<Number>& u)
{
  Assert (numbers::is_finite(a),
	  ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (u.val,
			 u.val+u.vec_size,
			 val,
			 (a*boost::lambda::_1),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = a * u.val[i];
}



template <typename Number>
template <typename Number2>
void Vector<Number>::equ (const Number a,
			  const Vector<Number2>& u)
{
  Assert (numbers::is_finite(a),
	  ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));

				   // set the result vector to a*u. we have to
				   // convert the elements of u to the type of
				   // the result vector. this is necessary
				   // because
				   // operator*(complex<float>,complex<double>)
				   // is not defined by default
  for (unsigned int i=0; i<vec_size; ++i)
    val[i] = a * Number(u.val[i]);
}



template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number>& u,
			  const Number b, const Vector<Number>& v)
{
  Assert (numbers::is_finite(a),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));
  Assert (numbers::is_finite(b),
          ExcMessage("The given value is not finite but either infinite or Not A Number (NaN)"));

  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (u.val,
			 u.val+u.vec_size,
			 v.val,
			 val,
			 (a*boost::lambda::_1 + b*boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = a * u.val[i] + b * v.val[i];
}


template <typename Number>
void Vector<Number>::equ (const Number a, const Vector<Number>& u,
			  const Number b, const Vector<Number>& v,
			  const Number c, const Vector<Number>& w)
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == u.vec_size, ExcDimensionMismatch(vec_size, u.vec_size));
  Assert (vec_size == v.vec_size, ExcDimensionMismatch(vec_size, v.vec_size));
  Assert (vec_size == w.vec_size, ExcDimensionMismatch(vec_size, w.vec_size));

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (u.val,
			 u.val+u.vec_size,
			 v.val,
			 w.val,
			 val,
			 (a*boost::lambda::_1 + b*boost::lambda::_2 + c*boost::lambda::_3),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = a * u.val[i] + b * v.val[i] + c * w.val[i];
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

  if (vec_size>internal::Vector::minimum_parallel_grain_size)
    parallel::transform (a.val,
			 a.val+a.vec_size,
			 b.val,
			 val,
			 (boost::lambda::_1 / boost::lambda::_2),
			 internal::Vector::minimum_parallel_grain_size);
  else if (vec_size > 0)
    for (unsigned int i=0; i<vec_size; ++i)
      val[i] = a.val[i]/b.val[i];
}



template <typename Number>
Vector<Number> &
Vector<Number>::operator = (const BlockVector<Number>& v)
{
  if (v.size() != vec_size)
    reinit (v.size(), true);

  unsigned int this_index = 0;
  for (unsigned int b=0; b<v.n_blocks(); ++b)
    for (unsigned int i=0; i<v.block(b).size(); ++i, ++this_index)
      val[this_index] = v.block(b)(i);

  return *this;
}



#ifdef DEAL_II_USE_PETSC

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
      int ierr = VecGetArray (static_cast<const Vec&>(v), &start_ptr);
      AssertThrow (ierr == 0, ExcPETScError(ierr));

      internal::copy (start_ptr, start_ptr+vec_size, begin());

                                       // restore the representation of the
                                       // vector
      ierr = VecRestoreArray (static_cast<const Vec&>(v), &start_ptr);
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


#ifdef DEAL_II_USE_TRILINOS

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
Vector<Number>::operator == (const Vector<Number2>& v) const
{
  Assert (vec_size!=0, ExcEmptyObject());
  Assert (vec_size == v.size(), ExcDimensionMismatch(vec_size, v.size()));

				   // compare the two vector. we have to
				   // convert the elements of v to the type of
				   // the result vector. this is necessary
				   // because
				   // operator==(complex<float>,complex<double>)
				   // is not defined by default
  for (unsigned int i=0; i<vec_size; ++i)
    if (val[i] != Number(v.val[i]))
      return false;

  return true;
}



template <typename Number>
void Vector<Number>::print (const char *format) const
{
  Assert (vec_size!=0, ExcEmptyObject());

  for (unsigned int j=0; j<size(); ++j)
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
    for (unsigned int i=0; i<size(); ++i)
      out << val[i] << ' ';
  else
    for (unsigned int i=0; i<size(); ++i)
      out << val[i] << std::endl;
  out << std::endl;

  AssertThrow (out, ExcIO());
                                   // reset output format
  out.flags (old_flags);
  out.precision(old_precision);
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
  const unsigned int sz = size();
  char buf[16];

  std::sprintf(buf, "%d", sz);
  std::strcat(buf, "\n[");

  out.write(buf, std::strlen(buf));
  out.write (reinterpret_cast<const char*>(begin()),
	     reinterpret_cast<const char*>(end())
	     - reinterpret_cast<const char*>(begin()));

				   // out << ']';
  const char outro = ']';
  out.write (&outro, 1);

  AssertThrow (out, ExcIO());
}



template <typename Number>
void Vector<Number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  unsigned int sz;

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

  in.read (reinterpret_cast<char*>(begin()),
	   reinterpret_cast<const char*>(end())
	   - reinterpret_cast<const char*>(begin()));

				   //  in >> c;
  in.read (&c, 1);
  AssertThrow (c==']', ExcIO());
}



template <typename Number>
unsigned int
Vector<Number>::memory_consumption () const
{
  return sizeof(*this) + (max_vec_size * sizeof(Number));
}


DEAL_II_NAMESPACE_CLOSE

#endif
