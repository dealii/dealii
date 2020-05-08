// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2019 by the deal.II authors
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

#ifndef dealii_cuda_vector_h
#define dealii_cuda_vector_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/index_set.h>

#include <deal.II/lac/vector_operation.h>
#include <deal.II/lac/vector_space_vector.h>

#ifdef DEAL_II_WITH_CUDA

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#  ifndef DOXYGEN
class CommunicationPatternBase;
template <typename Number>
class ReadWriteVector;
#  endif

namespace LinearAlgebra
{
  /**
   * A Namespace for the CUDA vectors.
   */
  namespace CUDAWrappers
  {
    /**
     * This class implements a vector using CUDA for use on Nvidia GPUs. This
     * class is derived from the LinearAlgebra::VectorSpaceVector class.
     *
     * @note Only float and double are supported.
     *
     * @see CUDAWrappers
     * @ingroup Vectors
     * @author Karl Ljungkvist, Bruno Turcksin, Daniel Arndt, 2016, 2018
     */
    template <typename Number>
    class Vector : public VectorSpaceVector<Number>
    {
    public:
      using value_type = typename VectorSpaceVector<Number>::value_type;
      using size_type  = typename VectorSpaceVector<Number>::size_type;
      using real_type  = typename VectorSpaceVector<Number>::real_type;

      /**
       * Constructor. Create a vector of dimension zero.
       */
      Vector();

      /**
       * Copy constructor.
       */
      Vector(const Vector<Number> &V);

      /**
       * Move constructor.
       */
      Vector(Vector<Number> &&) noexcept = default;

      /**
       * Constructor. Set dimension to @p n and initialize all elements with
       * zero.
       *
       * The constructor is made explicit to avoid accident like this:
       * <tt>v=0;</tt>. Presumably, the user wants to set every elements of
       * the vector to zero, but instead, what happens is this call:
       * <tt>v=Vector@<Number@>(0);</tt>, i.e. the vector is replaced by one
       * of length zero.
       */
      explicit Vector(const size_type n);

      /**
       * Copy assignment operator.
       */
      Vector &
      operator=(const Vector<Number> &v);

      /**
       * Move assignment operator.
       */
      Vector &
      operator=(Vector<Number> &&v) noexcept = default;

      /**
       * Swap the contents of this vector and the other vector @p v. One could do
       * this operation with a temporary variable and copying over the data
       * elements, but this function is significantly more efficient since it
       * only swaps the pointers to the data of the two vectors and therefore
       * does not need to allocate temporary storage and move data around.
       *
       * This function is analogous to the @p swap function of all C++
       * standard containers. Also, there is a global function
       * <tt>swap(u,v)</tt> that simply calls <tt>u.swap(v)</tt>, again in
       * analogy to standard functions.
       *
       * This function is virtual in order to allow for derived classes to
       * handle memory separately.
       */
      virtual void
      swap(Vector<Number> &v);

      /**
       * Reinit functionality. The flag <tt>omit_zeroing_entries</tt>
       * determines whether the vector should be filled with zero (false) or
       * left untouched (true).
       */
      void
      reinit(const size_type n, const bool omit_zeroing_entries = false);

      /**
       * Change the dimension to that of the vector V. The elements of V are not
       * copied.
       */
      virtual void
      reinit(const VectorSpaceVector<Number> &V,
             const bool omit_zeroing_entries = false) override;

      /**
       * Import all the element from the input vector @p V.
       * VectorOperation::values @p operation is used to decide if the
       * elements int @p V should be added to the current vector or replace
       * the current elements. The last parameter is not used. It is only used
       * for distributed vectors. This is the function that should be used to
       * copy a vector to the GPU.
       */
      virtual void
      import(
        const ReadWriteVector<Number> &                 V,
        VectorOperation::values                         operation,
        std::shared_ptr<const CommunicationPatternBase> communication_pattern =
          std::shared_ptr<const CommunicationPatternBase>()) override;

      /**
       * Sets all elements of the vector to the scalar @p s. This operation is
       * only allowed if @p s is equal to zero.
       */
      virtual Vector<Number> &
      operator=(const Number s) override;

      /**
       * Multiply the entive vector by a fixed factor.
       */
      virtual Vector<Number> &
      operator*=(const Number factor) override;

      /**
       * Divide the entire vector by a fixed factor.
       */
      virtual Vector<Number> &
      operator/=(const Number factor) override;

      /**
       * Add the vector @p V to the present one.
       */
      virtual Vector<Number> &
      operator+=(const VectorSpaceVector<Number> &V) override;

      /**
       * Subtract the vector @p V from the present one.
       */
      virtual Vector<Number> &
      operator-=(const VectorSpaceVector<Number> &V) override;

      /**
       * Return the scalar product of two vectors.
       */
      virtual Number
      operator*(const VectorSpaceVector<Number> &V) const override;

      /**
       * Add @p to all components. Note that @p a is a scalar not a vector.
       */
      virtual void
      add(const Number a) override;

      /**
       * Simple addition of a multiple of a vector, i.e. <tt>*this += a*V</tt>.
       */
      virtual void
      add(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Multiple additions of scaled vectors, i.e. <tt>*this += a*V+b*W</tt>.
       */
      virtual void
      add(const Number                     a,
          const VectorSpaceVector<Number> &V,
          const Number                     b,
          const VectorSpaceVector<Number> &W) override;

      /**
       * Scaling and simple addition of a multiple of a vector, i.e. <tt>*this
       * = s*(*this)+a*V</tt>
       */
      virtual void
      sadd(const Number                     s,
           const Number                     a,
           const VectorSpaceVector<Number> &V) override;

      /**
       * Scale each element of this vector by the corresponding element in the
       * argument. This function is mostly meant to simulate multiplication
       * (and immediate re-assignment) by a diagonal scaling matrix.
       */
      virtual void
      scale(const VectorSpaceVector<Number> &scaling_factors) override;

      /**
       * Assignment <tt>*this = a*V</tt>.
       */
      virtual void
      equ(const Number a, const VectorSpaceVector<Number> &V) override;

      /**
       * Return whether the vector contains only elements with value zero.
       */
      virtual bool
      all_zero() const override;

      /**
       * Return the mean value of all the entries of this vector.
       */
      virtual value_type
      mean_value() const override;

      /**
       * Return the l<sub>1</sub> norm of the vector (i.e., the sum of the
       * absolute values of all entries among all processors).
       */
      virtual real_type
      l1_norm() const override;

      /**
       * Return the l<sub>2</sub> norm of the vector (i.e., the square root of
       * the sum of the square of all entries among all processors).
       */
      virtual real_type
      l2_norm() const override;

      /**
       * Return the square of the $l_2$-norm.
       */
      real_type
      norm_sqr() const;

      /**
       * Return the maximum norm of the vector (i.e., the maximum absolute
       * value among all entries and among all processors).
       */
      virtual real_type
      linfty_norm() const override;

      /**
       * Perform a combined operation of a vector addition and a subsequent
       * inner product, returning the value of the inner product. In other
       * words, the result of this function is the same as if the user called
       * @code
       * this->add(a, V);
       * return_value = *this * W;
       * @endcode
       *
       * The reason this function exists is that this operation involves less
       * memory transfer than calling the two functions separately. This
       * method only needs to load three vectors, @p this, @p V, @p W, whereas
       * calling separate methods means to load the calling vector @p this
       * twice. Since most vector operations are memory transfer limited, this
       * reduces the time by 25\% (or 50\% if @p W equals @p this).
       *
       * For complex-valued vectors, the scalar product in the second step is
       * implemented as $\left<v,w\right>=\sum_i v_i \bar{w_i}$.
       */
      virtual Number
      add_and_dot(const Number                     a,
                  const VectorSpaceVector<Number> &V,
                  const VectorSpaceVector<Number> &W) override;

      /**
       * Return the pointer to the underlying array. Ownership still resides
       * with this class.
       */
      Number *
      get_values() const;

      /**
       * Return the size of the vector.
       */
      virtual size_type
      size() const override;

      /**
       * Return an index set that describe which elements of this vector are
       * owned by the current processor, i.e. [0, size).
       */
      virtual dealii::IndexSet
      locally_owned_elements() const override;

      /**
       * Print the vector to the output stream @p out.
       */
      virtual void
      print(std::ostream &     out,
            const unsigned int precision  = 2,
            const bool         scientific = true,
            const bool         across     = true) const override;

      /**
       * Return the memory consumption of this class in bytes.
       */
      virtual std::size_t
      memory_consumption() const override;

      /**
       * Attempt to perform an operation between two incompatible vector types.
       *
       * @ingroup Exceptions
       */
      DeclException0(ExcVectorTypeNotCompatible);

    private:
      /**
       * Pointer to the array of elements of this vector.
       */
      std::unique_ptr<Number[], void (*)(Number *)> val;

      /**
       * Number of elements in the vector.
       */
      size_type n_elements;
    };
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

// ---------------------------- Inline functions --------------------------

/**
 * Global function @p swap which overloads the default implementation of the
 * C++ standard library which uses a temporary object. The function simply
 * exchanges the data of the two vectors.
 *
 * @relatesalso Vector
 * @author Daniel Arndt, 2018
 */
template <typename Number>
inline void
swap(LinearAlgebra::CUDAWrappers::Vector<Number> &u,
     LinearAlgebra::CUDAWrappers::Vector<Number> &v)
{
  u.swap(v);
}

namespace LinearAlgebra
{
  namespace CUDAWrappers
  {
    template <typename Number>
    inline Number *
    Vector<Number>::get_values() const
    {
      return val.get();
    }



    template <typename Number>
    inline typename Vector<Number>::size_type
    Vector<Number>::size() const
    {
      return n_elements;
    }


    template <typename Number>
    inline IndexSet
    Vector<Number>::locally_owned_elements() const
    {
      return complete_index_set(n_elements);
    }



    template <typename Number>
    inline void
    Vector<Number>::swap(Vector<Number> &v)
    {
      std::swap(val, v.val);
      std::swap(n_elements, v.n_elements);
    }
  } // namespace CUDAWrappers
} // namespace LinearAlgebra

DEAL_II_NAMESPACE_CLOSE

#endif

#endif
