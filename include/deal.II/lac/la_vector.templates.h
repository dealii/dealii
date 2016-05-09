// ---------------------------------------------------------------------
//
// Copyright (C)  2015 by the deal.II authors
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

#ifndef dealii__la_vector_templates_h
#define dealii__la_vector_templates_h

#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/vector_operations_internal.h>
#include <iostream>
#include <iomanip>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  template <typename Number>
  Vector<Number> &Vector<Number>::operator= (const Vector<Number> &in_vector)
  {
    if (PointerComparison::equal(this, &in_vector))
      return *this;

    this->thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (this->size() != in_vector.size())
      this->reinit(in_vector, true);

    dealii::internal::Vector_copy<Number, Number> copier(in_vector.val, this->val);
    internal::parallel_for(copier, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  template <typename Number2>
  Vector<Number> &Vector<Number>::operator= (const Vector<Number2> &in_vector)
  {
    this->thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (this->size() != in_vector.size())
      this->reinit(in_vector, true);

    dealii::internal::Vector_copy<Number, Number2> copier(in_vector.val, this->val);
    internal::parallel_for(copier, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &Vector<Number>::operator= (const Number s)
  {
    Assert(s==static_cast<Number>(0), ExcMessage("Only 0 can be assigned to a vector."));
    (void) s;

    internal::Vector_set<Number> setter(Number(), this->val);
    internal::parallel_for(setter, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &Vector<Number>::operator*= (const Number factor)
  {
    AssertIsFinite(factor);

    internal::Vectorization_multiply_factor<Number> vector_multiply(this->val, factor);
    internal::parallel_for(vector_multiply, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &Vector<Number>::operator/= (const Number factor)
  {
    AssertIsFinite(factor);
    Assert(factor!=Number(0.), ExcZero());
    this->operator *= (Number(1.)/factor);

    return *this;
  }



  template <typename Number>
  Vector<Number> &Vector<Number>::operator+= (const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL, ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));

    internal::Vectorization_add_v<Number> vector_add(this->val, down_V.val);
    internal::parallel_for(vector_add, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &Vector<Number>::operator-= (const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL, ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot subtract two vectors with different numbers of elements"));
    internal::Vectorization_subtract_v<Number> vector_subtract(this->val, down_V.val);
    internal::parallel_for(vector_subtract, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Number Vector<Number>::operator* (const VectorSpaceVector<Number> &V) const
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot compute the scalar product "
                      "of two vectors with different numbers of elements"));
    Number sum;
    internal::Dot<Number, Number> dot(this->val, down_V.val);
    internal::parallel_reduce(dot, this->size(), sum, this->thread_loop_partitioner);

    return sum;
  }



  template <typename Number>
  void Vector<Number>::import(const ReadWriteVector<Number> &,
                              VectorOperation::values        ,
                              std_cxx11::shared_ptr<const CommunicationPatternBase>)
  {
    AssertThrow(false, ExcMessage("This function is not implemented."));
  }



  template <typename Number>
  inline
  void Vector<Number>::add(const Number a)
  {
    AssertIsFinite(a);

    internal::Vectorization_add_factor<Number> vector_add(this->val, a);
    internal::parallel_for(vector_add, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    AssertIsFinite(a);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));

    internal::Vectorization_add_av<Number> vector_add_av(this->val, down_V.val, a);
    internal::parallel_for(vector_add_av, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V,
                           const Number b, const VectorSpaceVector<Number> &W)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&W)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    // Downcast W. If fails, throws an exception.
    const Vector<Number> &down_W = dynamic_cast<const Vector<Number>&>(W);
    AssertIsFinite(a);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));
    AssertIsFinite(b);
    Assert(down_W.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));

    internal::Vectorization_add_avpbw<Number> vector_add(this->val, down_V.val,
                                                         down_W.val, a, b);
    internal::parallel_for(vector_add, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Number a,
                            const VectorSpaceVector<Number> &V)
  {
    AssertIsFinite(s);
    AssertIsFinite(a);

    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. It fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    internal::Vectorization_sadd_xav<Number> vector_sadd_xav(this->val, down_V.val,
                                                             a, s);
    internal::parallel_for(vector_sadd_xav, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&scaling_factors)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast scaling_factors. If fails, throws an exception.
    const Vector<Number> &down_scaling_factors =
      dynamic_cast<const Vector<Number>&>(scaling_factors);
    Assert(down_scaling_factors.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));

    internal::Vectorization_scale<Number> vector_scale(this->val, down_scaling_factors.val);
    internal::parallel_for(vector_scale, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
  {
    AssertIsFinite(a);

    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    internal::Vectorization_equ_au<Number> vector_equ(this->val, down_V.val, a);
    internal::parallel_for(vector_equ, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l1_norm() const
  {
    Assert (this->size(), ExcEmptyObject());

    typedef typename VectorSpaceVector<Number>::real_type real_type;
    real_type sum;
    internal::Norm1<Number, real_type> norm1(this->val);
    internal::parallel_reduce(norm1, this->size(), sum, this->thread_loop_partitioner);

    return sum;
  }




  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l2_norm() const
  {
    Assert (this->size(), ExcEmptyObject());

    // if l2_norm()^2 is finite and non-zero, the answer is computed as
    // std::sqrt(norm_sqr()). If norm_sqrt() is infinite or zero, the l2 norm
    // might still be finite. In that case, recompute ig (this is a rare case,
    // so working on the vector twice is uncritical and paid off by the extended
    // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
    typedef typename VectorSpaceVector<Number>::real_type real_type;
    real_type norm_square;
    internal::Norm2<Number, real_type> norm2(this->val);
    internal::parallel_reduce(norm2, this->size(), norm_square,
                              this->thread_loop_partitioner);
    if (numbers::is_finite(norm_square) &&
        norm_square>=std::numeric_limits<real_type>::min())
      return std::sqrt(norm_square);
    else
      {
        real_type scale = 0.;
        real_type sum = 1.;
        const size_type size = this->size();
        for (size_type i=0; i<size; ++i)
          {
            if (this->val[i] != Number())
              {
                const real_type abs_x =
                  numbers::NumberTraits<Number>::abs(this->val[i]);
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
        return scale*std::sqrt(sum);
      }
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::linfty_norm() const
  {
    typename ReadWriteVector<Number>::real_type norm = 0.;
    const size_type size = this->size();
    for (size_type i=0; i<size; ++i)
      norm = std::max(std::abs(this->val[i]),norm);

    return norm;
  }



  template <typename Number>
  Number Vector<Number>::add_and_dot(const Number a,
                                     const VectorSpaceVector<Number> &V,
                                     const VectorSpaceVector<Number> &W)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&V)!=NULL,
           ExcVectorTypeNotCompatible());
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number>*>(&W)!=NULL,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    // Downcast W. If fails, throws an exception.
    const Vector<Number> &down_W = dynamic_cast<const Vector<Number>&>(W);
    AssertIsFinite(a);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));
    Assert(down_W.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));

    Number sum;
    internal::AddAndDot<Number> adder(this->val, down_V.val, down_W.val, a);
    internal::parallel_reduce(adder, this->size(), sum, this->thread_loop_partitioner);
    AssertIsFinite(sum);

    return sum;
  }



  template <typename Number>
  void Vector<Number>::block_write(std::ostream &out) const
  {
    AssertThrow(out, ExcIO());

    // Other version of the following
    //  out << size() << std::endl << '[';
    // Reason: operator<< seems to use some resources  that lead to problems in
    // a multithreaded environment.
    const size_type sz = this->size();
    char buf[16];
#ifdef DEAL_II_WITH_64BIT_INDICES
    std::sprintf(buf, "%llu", sz);
#else
    std::sprintf(buf, "%u", sz);
#endif
    std::strcat(buf, "\n[");

    out.write(buf, std::strlen(buf));
    out.write(reinterpret_cast<const char *>(this->begin()),
              reinterpret_cast<const char *>(this->end())
              - reinterpret_cast<const char *>(this->begin()));

    // out << ']';
    const char outro = ']';
    out.write(&outro, 1);

    AssertThrow(out, ExcIO());
  }



  template <typename Number>
  void Vector<Number>::block_read(std::istream &in)
  {
    AssertThrow(in, ExcIO());

    size_type sz;

    char buf[16];

    in.getline(buf,16,'\n');
    sz = std::atoi(buf);
    this->reinit(sz,true);

    char c;
    // in >> c;
    in.read(&c,1);
    AssertThrow(c=='[', ExcIO());

    in.read(reinterpret_cast<char *>(this->begin()),
            reinterpret_cast<const char *>(this->end())
            - reinterpret_cast<const char *>(this->begin()));

    // in >> c;
    in.read(&c,1);
    AssertThrow(c==']', ExcIO());
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
