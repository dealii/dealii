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

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  template <typename Number>
  VectorSpaceVector<Number> &Vector<Number>::operator*= (const Number factor)
  {
    AssertIsFinite(factor);
    for (unsigned int i=0; i<this->size(); ++i)
      this->val[i] *= factor;

    return *this;
  }



  template <typename Number>
  VectorSpaceVector<Number> &Vector<Number>::operator/= (const Number factor)
  {
    AssertIsFinite(factor);
    Assert(factor!=Number(0.), ExcZero());
    this->operator *= (Number(1.)/factor);

    return *this;
  }



  template <typename Number>
  VectorSpaceVector<Number> &Vector<Number>::operator+= (const VectorSpaceVector<Number> &V)
  {
    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_V[i]);
        this->val[i] += down_V[i];
      }

    return *this;
  }



  template <typename Number>
  VectorSpaceVector<Number> &Vector<Number>::operator-= (const VectorSpaceVector<Number> &V)
  {
    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot subtract two vectors with different numbers of elements"));
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_V[i]);
        this->val[i] -= down_V[i];
      }

    return *this;
  }



  template <typename Number>
  Number Vector<Number>::operator* (const VectorSpaceVector<Number> &V)
  {
    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot compute the scalar product "
                      "of two vectors with different numbers of elements"));
    Number value = 0.;
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_V[i]);
        value += this->val[i]*down_V[i];
      }

    return value;
  }



  template <typename Number>
  void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
  {
    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    AssertIsFinite(a);
    Assert(down_V.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_V[i]);
        this->val[i] += a*down_V[i];
      }
  }



  template <typename Number>
  void Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V,
                           const Number b, const VectorSpaceVector<Number> &W)
  {
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
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_V[i]);
        AssertIsFinite(down_W[i]);
        this->val[i] += a*down_V[i]+b*down_W[i];
      }
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const VectorSpaceVector<Number> &V)
  {
    *this *= s;
    *this += V;
  }



  template <typename Number>
  void Vector<Number>::sadd(const Number s, const Number a,
                            const VectorSpaceVector<Number> &V)
  {
    *this *= s;
    // Downcast V. It fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    Vector<Number> tmp(down_V);
    tmp *= a;
    *this += tmp;
  }



  template <typename Number>
  void Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
  {
    // Downcast scaling_factors. If fails, throws an exception.
    const Vector<Number> &down_scaling_factors =
      dynamic_cast<const Vector<Number>&>(scaling_factors);
    Assert(down_scaling_factors.size()==this->size(),
           ExcMessage("Cannot add two vectors with different numbers of elements"));
    for (unsigned int i=0; i<this->size(); ++i)
      {
        AssertIsFinite(down_scaling_factors[i]);
        this->val[i] *= down_scaling_factors[i];
      }
  }



  template <typename Number>
  void Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
  {
    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number>&>(V);
    *this = down_V;
    *this *= a;
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l1_norm()
  {
    Assert (this->size(), ExcEmptyObject());

    typename ReadWriteVector<Number>::real_type norm = 0.;
    norm = l1_norm_recursive(0,this->size()-1);

    return norm;
  }




  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l2_norm()
  {
    Assert (this->size(), ExcEmptyObject());

    // if l2_norm()^2 is finite and non-zero, the answer is computed as
    // std::sqrt(norm_sqr()). If norm_sqrt() is infinite or zero, the l2 norm
    // might still be finite. In that case, recompute ig (this is a rare case,
    // so working on the vector twice is uncritical and paid off by the extended
    // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
    typedef typename ReadWriteVector<Number>::real_type real_type;
    real_type norm_square = 0.;
    norm_square = l2_norm_recursive(0,this->size()-1);
    if (numbers::is_finite(norm_square) &&
        norm_square>=std::numeric_limits<real_type>::min())
      return std::sqrt(norm_square);
    else
      {
        real_type scale = 0.;
        real_type sum = 1.;
        for (unsigned int i=0; i<this->size(); ++i)
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
  typename VectorSpaceVector<Number>::real_type Vector<Number>::linfty_norm()
  {
    typename ReadWriteVector<Number>::real_type norm = 0.;
    for (unsigned int i=0; i<this->size(); ++i)
      if (std::abs(this->val[i])>norm)
        norm = std::abs(this->val[i]);

    return norm;
  }



  template <typename Number>
  Number Vector<Number>::add_and_dot(const Number a,
                                     const VectorSpaceVector<Number> &V,
                                     const VectorSpaceVector<Number> &W)
  {
    this->add(a,V);

    return *this*W;
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l1_norm_recursive(unsigned int i,
      unsigned int j)
  {
    Assert(j>=i, ExcInternalError());
    typename ReadWriteVector<Number>::real_type norm = 0.;

    if ((j-i)!=0)
      {
        norm += l1_norm_recursive(i,(i+j)/2);
        norm += l1_norm_recursive((i+j)/2+1,j);
      }
    else
      norm += std::abs(this->val[i]);

    return norm;
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type Vector<Number>::l2_norm_recursive(unsigned int i,
      unsigned int j)
  {
    Assert(j>=i, ExcInternalError());

    typename ReadWriteVector<Number>::real_type norm = 0.;
    if ((j-i)!=0)
      {
        norm += l2_norm_recursive(i,(i+j)/2);
        norm += l2_norm_recursive((i+j)/2+1,j);
      }
    else
      norm += std::pow(std::abs(this->val[i]),2);

    return norm;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
