// ---------------------------------------------------------------------
//
// Copyright (C) 2015 - 2020 by the deal.II authors
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

#ifndef dealii_la_vector_templates_h
#define dealii_la_vector_templates_h

#include <deal.II/base/config.h>

#include <deal.II/lac/la_vector.h>
#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_operations_internal.h>

#include <boost/io/ios_state.hpp>

#include <iomanip>
#include <iostream>

DEAL_II_NAMESPACE_OPEN

namespace LinearAlgebra
{
  template <typename Number>
  void
  Vector<Number>::reinit(const size_type size, const bool omit_zeroing_entries)
  {
    ReadWriteVector<Number>::reinit(size, omit_zeroing_entries);
  }



  template <typename Number>
  template <typename Number2>
  void
  Vector<Number>::reinit(const ReadWriteVector<Number2> &in_vector,
                         const bool                      omit_zeroing_entries)
  {
    ReadWriteVector<Number>::reinit(in_vector, omit_zeroing_entries);
  }



  template <typename Number>
  void
  Vector<Number>::reinit(const IndexSet &locally_stored_indices,
                         const bool      omit_zeroing_entries)
  {
    ReadWriteVector<Number>::reinit(locally_stored_indices,
                                    omit_zeroing_entries);
  }



  template <typename Number>
  void
  Vector<Number>::reinit(const VectorSpaceVector<Number> &V,
                         const bool                       omit_zeroing_entries)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    ReadWriteVector<Number>::reinit(down_V, omit_zeroing_entries);
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator=(const Vector<Number> &in_vector)
  {
    if (PointerComparison::equal(this, &in_vector))
      return *this;

    this->thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (this->size() != in_vector.size())
      this->reinit(in_vector, true);

    dealii::internal::VectorOperations::Vector_copy<Number, Number> copier(
      in_vector.values.get(), this->values.get());
    dealii::internal::VectorOperations::parallel_for(
      copier,
      static_cast<size_type>(0),
      this->size(),
      this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  template <typename Number2>
  Vector<Number> &
  Vector<Number>::operator=(const Vector<Number2> &in_vector)
  {
    this->thread_loop_partitioner = in_vector.thread_loop_partitioner;
    if (this->size() != in_vector.size())
      ReadWriteVector<Number>::reinit(in_vector, true);

    dealii::internal::VectorOperations::Vector_copy<Number, Number2> copier(
      in_vector.values.get(), this->values.get());
    dealii::internal::VectorOperations::parallel_for(
      copier,
      static_cast<size_type>(0),
      this->size(),
      this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator=(const Number s)
  {
    Assert(s == static_cast<Number>(0),
           ExcMessage("Only 0 can be assigned to a vector."));
    (void)s;

    dealii::internal::VectorOperations::Vector_set<Number> setter(
      Number(), this->values.get());
    dealii::internal::VectorOperations::parallel_for(
      setter, 0, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator*=(const Number factor)
  {
    AssertIsFinite(factor);

    dealii::internal::VectorOperations::Vectorization_multiply_factor<Number>
      vector_multiply(this->values.get(), factor);
    dealii::internal::VectorOperations::parallel_for(
      vector_multiply,
      static_cast<size_type>(0),
      this->size(),
      this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator/=(const Number factor)
  {
    AssertIsFinite(factor);
    Assert(factor != Number(0.), ExcZero());
    this->operator*=(Number(1.) / factor);

    return *this;
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator+=(const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    dealii::internal::VectorOperations::Vectorization_add_v<Number> vector_add(
      this->values.get(), down_V.values.get());
    dealii::internal::VectorOperations::parallel_for(
      vector_add, 0, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Vector<Number> &
  Vector<Number>::operator-=(const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot subtract two vectors with different numbers of elements"));
    dealii::internal::VectorOperations::Vectorization_subtract_v<Number>
      vector_subtract(this->values.get(), down_V.values.get());
    dealii::internal::VectorOperations::parallel_for(
      vector_subtract, 0, this->size(), this->thread_loop_partitioner);

    return *this;
  }



  template <typename Number>
  Number Vector<Number>::operator*(const VectorSpaceVector<Number> &V) const
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    Assert(down_V.size() == this->size(),
           ExcMessage("Cannot compute the scalar product "
                      "of two vectors with different numbers of elements"));
    Number                                                  sum;
    dealii::internal::VectorOperations::Dot<Number, Number> dot(
      this->values.get(), down_V.values.get());
    dealii::internal::VectorOperations::parallel_reduce(
      dot, 0, this->size(), sum, this->thread_loop_partitioner);

    return sum;
  }



  template <typename Number>
  void
  Vector<Number>::import(const ReadWriteVector<Number> &,
                         VectorOperation::values,
                         std::shared_ptr<const CommunicationPatternBase>)
  {
    AssertThrow(false, ExcMessage("This function is not implemented."));
  }



  template <typename Number>
  inline void
  Vector<Number>::add(const Number a)
  {
    AssertIsFinite(a);

    dealii::internal::VectorOperations::Vectorization_add_factor<Number>
      vector_add(this->values.get(), a);
    dealii::internal::VectorOperations::parallel_for(
      vector_add, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void
  Vector<Number>::add(const Number a, const VectorSpaceVector<Number> &V)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    AssertIsFinite(a);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    dealii::internal::VectorOperations::Vectorization_add_av<Number>
      vector_add_av(this->values.get(), down_V.values.get(), a);
    dealii::internal::VectorOperations::parallel_for(
      vector_add_av, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void
  Vector<Number>::add(const Number                     a,
                      const VectorSpaceVector<Number> &V,
                      const Number                     b,
                      const VectorSpaceVector<Number> &W)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&W) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    // Downcast W. If fails, throws an exception.
    const Vector<Number> &down_W = dynamic_cast<const Vector<Number> &>(W);
    AssertIsFinite(a);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));
    AssertIsFinite(b);
    Assert(down_W.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    dealii::internal::VectorOperations::Vectorization_add_avpbw<Number>
      vector_add(
        this->values.get(), down_V.values.get(), down_W.values.get(), a, b);
    dealii::internal::VectorOperations::parallel_for(
      vector_add, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void
  Vector<Number>::sadd(const Number                     s,
                       const Number                     a,
                       const VectorSpaceVector<Number> &V)
  {
    AssertIsFinite(s);
    AssertIsFinite(a);

    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. It fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    dealii::internal::VectorOperations::Vectorization_sadd_xav<Number>
      vector_sadd_xav(this->values.get(), down_V.values.get(), a, s);
    dealii::internal::VectorOperations::parallel_for(
      vector_sadd_xav, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void
  Vector<Number>::scale(const VectorSpaceVector<Number> &scaling_factors)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&scaling_factors) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast scaling_factors. If fails, throws an exception.
    const Vector<Number> &down_scaling_factors =
      dynamic_cast<const Vector<Number> &>(scaling_factors);
    Assert(down_scaling_factors.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    dealii::internal::VectorOperations::Vectorization_scale<Number>
      vector_scale(this->values.get(), down_scaling_factors.values.get());
    dealii::internal::VectorOperations::parallel_for(
      vector_scale, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  void
  Vector<Number>::equ(const Number a, const VectorSpaceVector<Number> &V)
  {
    AssertIsFinite(a);

    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    dealii::internal::VectorOperations::Vectorization_equ_au<Number> vector_equ(
      this->values.get(), down_V.values.get(), a);
    dealii::internal::VectorOperations::parallel_for(
      vector_equ, 0, this->size(), this->thread_loop_partitioner);
  }



  template <typename Number>
  bool
  Vector<Number>::all_zero() const
  {
    Assert(this->size(), ExcEmptyObject());

    const size_type size = this->size();
    for (size_type i = 0; i < size; ++i)
      if (this->values[i] != Number())
        return false;

    return true;
  }



  template <typename Number>
  typename Vector<Number>::value_type
  Vector<Number>::mean_value() const
  {
    Assert(this->size(), ExcEmptyObject());

    using real_type = typename VectorSpaceVector<Number>::real_type;
    value_type                                            sum;
    dealii::internal::VectorOperations::MeanValue<Number> mean_value(
      this->values.get());
    dealii::internal::VectorOperations::parallel_reduce(
      mean_value, 0, this->size(), sum, this->thread_loop_partitioner);

    return sum / static_cast<real_type>(this->size());
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type
  Vector<Number>::l1_norm() const
  {
    Assert(this->size(), ExcEmptyObject());

    using real_type = typename VectorSpaceVector<Number>::real_type;
    real_type                                                    sum;
    dealii::internal::VectorOperations::Norm1<Number, real_type> norm1(
      this->values.get());
    dealii::internal::VectorOperations::parallel_reduce(
      norm1, 0, this->size(), sum, this->thread_loop_partitioner);

    return sum;
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type
  Vector<Number>::l2_norm() const
  {
    Assert(this->size(), ExcEmptyObject());

    // if l2_norm()^2 is finite and non-zero, the answer is computed as
    // std::sqrt(norm_sqr()). If norm_sqrt() is infinite or zero, the l2 norm
    // might still be finite. In that case, recompute ig (this is a rare case,
    // so working on the vector twice is uncritical and paid off by the extended
    // precision) using the BLAS approach with a weight, see e.g. dnrm2.f.
    using real_type = typename VectorSpaceVector<Number>::real_type;
    real_type                                                    norm_square;
    dealii::internal::VectorOperations::Norm2<Number, real_type> norm2(
      this->values.get());
    dealii::internal::VectorOperations::parallel_reduce(
      norm2, 0, this->size(), norm_square, this->thread_loop_partitioner);
    if (numbers::is_finite(norm_square) &&
        norm_square >= std::numeric_limits<real_type>::min())
      return std::sqrt(norm_square);
    else
      {
        real_type       scale = 0.;
        real_type       sum   = 1.;
        const size_type size  = this->size();
        for (size_type i = 0; i < size; ++i)
          {
            if (this->values[i] != Number())
              {
                const real_type abs_x =
                  numbers::NumberTraits<Number>::abs(this->values[i]);
                if (scale < abs_x)
                  {
                    sum   = 1. + sum * (scale / abs_x) * (scale / abs_x);
                    scale = abs_x;
                  }
                else
                  sum += (abs_x / scale) * (abs_x / scale);
              }
          }
        AssertIsFinite(scale * std::sqrt(sum));
        return scale * std::sqrt(sum);
      }
  }



  template <typename Number>
  typename VectorSpaceVector<Number>::real_type
  Vector<Number>::linfty_norm() const
  {
    typename ReadWriteVector<Number>::real_type norm = 0.;
    const size_type                             size = this->size();
    for (size_type i = 0; i < size; ++i)
      norm = std::max(std::abs(this->values[i]), norm);

    return norm;
  }



  template <typename Number>
  Number
  Vector<Number>::add_and_dot(const Number                     a,
                              const VectorSpaceVector<Number> &V,
                              const VectorSpaceVector<Number> &W)
  {
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&V) != nullptr,
           ExcVectorTypeNotCompatible());
    // Check that casting will work.
    Assert(dynamic_cast<const Vector<Number> *>(&W) != nullptr,
           ExcVectorTypeNotCompatible());

    // Downcast V. If fails, throws an exception.
    const Vector<Number> &down_V = dynamic_cast<const Vector<Number> &>(V);
    // Downcast W. If fails, throws an exception.
    const Vector<Number> &down_W = dynamic_cast<const Vector<Number> &>(W);
    AssertIsFinite(a);
    Assert(down_V.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));
    Assert(down_W.size() == this->size(),
           ExcMessage(
             "Cannot add two vectors with different numbers of elements"));

    Number                                                sum;
    dealii::internal::VectorOperations::AddAndDot<Number> adder(
      this->values.get(), down_V.values.get(), down_W.values.get(), a);
    dealii::internal::VectorOperations::parallel_reduce(
      adder, 0, this->size(), sum, this->thread_loop_partitioner);
    AssertIsFinite(sum);

    return sum;
  }



  template <typename Number>
  void
  Vector<Number>::print_as_numpy_array(std::ostream &     out,
                                       const unsigned int precision) const
  {
    AssertThrow(out, ExcIO());
    boost::io::ios_flags_saver restore_flags(out);

    out.precision(precision);

    const unsigned int n_elements = this->n_elements();
    for (unsigned int i = 0; i < n_elements; ++i)
      out << this->values[i] << ' ';
    out << '\n' << std::flush;

    AssertThrow(out, ExcIO());
  }



  template <typename Number>
  void
  Vector<Number>::block_write(std::ostream &out) const
  {
    AssertThrow(out, ExcIO());

    // Other version of the following
    //  out << size() << std::endl << '[';
    // Reason: operator<< seems to use some resources that lead to problems in
    // a multithreaded environment.  We convert the size index to
    // unsigned long long int that is at least 64 bits to be able to output it
    // on all platforms, since std::uint64_t is not in C.
    const unsigned long long int sz = this->size();
    char                         buf[16];

    std::sprintf(buf, "%llu", sz);
    std::strcat(buf, "\n[");

    out.write(buf, std::strlen(buf));
    out.write(reinterpret_cast<const char *>(this->begin()),
              reinterpret_cast<const char *>(this->end()) -
                reinterpret_cast<const char *>(this->begin()));

    // out << ']';
    const char outro = ']';
    out.write(&outro, 1);

    AssertThrow(out, ExcIO());
  }



  template <typename Number>
  void
  Vector<Number>::block_read(std::istream &in)
  {
    AssertThrow(in, ExcIO());

    size_type sz;

    char buf[16];

    in.getline(buf, 16, '\n');
    sz = std::atoi(buf);
    ReadWriteVector<Number>::reinit(sz, true);

    char c;
    // in >> c;
    in.read(&c, 1);
    AssertThrow(c == '[', ExcIO());

    in.read(reinterpret_cast<char *>(this->begin()),
            reinterpret_cast<const char *>(this->end()) -
              reinterpret_cast<const char *>(this->begin()));

    // in >> c;
    in.read(&c, 1);
    AssertThrow(c == ']', ExcIO());
  }
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif
