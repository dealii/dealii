// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_signaling_nan_h
#define dealii_signaling_nan_h

#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/point.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <limits>


DEAL_II_NAMESPACE_OPEN

namespace numbers
{
  namespace internal
  {
    /**
     * A namespace for the implementation of functions that create signaling
     * NaN objects. This is where the numbers::signaling_nan() function
     * calls into.
     */
    namespace SignalingNaN
    {
      /**
       * A general template for classes that know how to initialize objects of
       * type @p T with signaling NaNs to denote invalid values.
       *
       * The real implementation of this class happens in (partial)
       * specializations for particular values of the template argument @p T.
       */
      template <typename T>
      struct NaNInitializer;


      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a @p float value equal to the invalid signaling
       * NaN.
       */
      template <>
      struct NaNInitializer<float>
      {
        static float
        invalid_element()
        {
          return std::numeric_limits<float>::signaling_NaN();
        }
      };


      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a @p double value equal to the invalid
       * signaling NaN.
       */
      template <>
      struct NaNInitializer<double>
      {
        static double
        invalid_element()
        {
          return std::numeric_limits<double>::signaling_NaN();
        }
      };


      template <typename T, std::size_t width>
      struct NaNInitializer<VectorizedArray<T, width>>
      {
        static VectorizedArray<T, width>
        invalid_element()
        {
          return VectorizedArray<T, width>(
            NaNInitializer<T>::invalid_element());
        }
      };


      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a Tensor<1,dim> value whose components are
       * invalid signaling NaN values.
       */
      template <int dim, typename T>
      struct NaNInitializer<Tensor<1, dim, T>>
      {
        static Tensor<1, dim, T>
        invalid_element()
        {
          Tensor<1, dim, T> nan_tensor;

          for (unsigned int i = 0; i < dim; ++i)
            nan_tensor[i] = NaNInitializer<T>::invalid_element();

          return nan_tensor;
        }
      };



      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a Tensor<rank,dim> value whose components are
       * invalid signaling NaN values.
       */
      template <int rank, int dim, typename T>
      struct NaNInitializer<Tensor<rank, dim, T>>
      {
        static Tensor<rank, dim, T>
        invalid_element()
        {
          Tensor<rank, dim, T> nan_tensor;

          // recursively initialize sub-tensors with invalid elements
          for (unsigned int i = 0; i < dim; ++i)
            nan_tensor[i] =
              NaNInitializer<Tensor<rank - 1, dim, T>>::invalid_element();

          return nan_tensor;
        }
      };



      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a Tensor<rank,dim> value whose components are
       * invalid signaling NaN values.
       */
      template <int dim, typename T>
      struct NaNInitializer<Point<dim, T>>
      {
        static Point<dim, T>
        invalid_element()
        {
          Point<dim, T> nan_point;

          for (unsigned int i = 0; i < dim; ++i)
            nan_point[i] = NaNInitializer<T>::invalid_element();

          return nan_point;
        }
      };



      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a SymmetricTensor<rank,dim> value whose
       * components are invalid signaling NaN values.
       */
      template <int rank, int dim, typename T>
      struct NaNInitializer<SymmetricTensor<rank, dim, T>>
      {
        static SymmetricTensor<rank, dim, T>
        invalid_element()
        {
          // initialize symmetric tensors via the unrolled list of elements
          T initializers
            [SymmetricTensor<rank, dim, T>::n_independent_components];
          for (unsigned int i = 0;
               i < SymmetricTensor<rank, dim, T>::n_independent_components;
               ++i)
            initializers[i] = NaNInitializer<T>::invalid_element();

          return SymmetricTensor<rank, dim, T>(initializers);
        }
      };



      /**
       * A specialization of the general NaNInitializer class that provides a
       * function that returns a DerivativeForm<order,dim,spacedim> value
       * whose components are invalid signaling NaN values.
       */
      template <int order, int dim, int spacedim, typename T>
      struct NaNInitializer<DerivativeForm<order, dim, spacedim, T>>
      {
        static DerivativeForm<order, dim, spacedim, T>
        invalid_element()
        {
          DerivativeForm<order, dim, spacedim, T> form;

          // recursively initialize sub-tensors with invalid elements
          for (unsigned int i = 0; i < spacedim; ++i)
            form[i] = NaNInitializer<Tensor<order, dim, T>>::invalid_element();

          return form;
        }
      };
    } // namespace SignalingNaN
  }   // namespace internal



  /**
   * Provide an object of type @p T filled with a signaling NaN that will
   * cause an exception when used in a computation. The content of these
   * objects is a "signaling NaN" ("NaN" stands for "not a number", and
   * "signaling" implies that at least on platforms where this is supported,
   * any arithmetic operation using them terminates the program). The purpose
   * of such objects is to use them as markers for invalid objects and
   * arrays that are required to be initialized to valid values at some point,
   * and to trigger an error when this later initialization does not happen
   * before the first use. An example is code such as this:
   * @code
   *   double x = numbers::signaling_nan<double>();
   *   if (some condition)
   *   {
   *     ...much code computing a,b,c...
   *     x = f(a,b,c);
   *   }
   *   else
   *   {
   *     ...more code...
   *     // bug: we forgot to assign a value to 'x'.
   *   }
   *
   *   return std::sin(x);
   * @endcode
   * The bug is that the `else` branch forgot to write a value into the `x`
   * variable. If your platform supports signaling NaNs, then this mistake
   * will become apparent in the last line above because the program is
   * going to be terminated by a floating point exception: The processor
   * realizes that the code is attempting to do an operation on the signaling
   * NaN still stored in `x` and aborts the program, thereby facilitating
   * an easy way to find what the problem is. This would not have been an easy
   * bug to find if one had just initialized `x` to zero in the first line
   * (or just left it uninitialized altogether): In that case, the call to
   * `std::sin` in the last line would have simply computed the sine of
   * "something" if `some condition == false`, but this invalid results may
   * not have been obvious to the calling site and would have required
   * a substantial amount of debugging to uncover because downstream
   * computations would simply have been wrong, without any indication of
   * *why* they are wrong.
   *
   * @tparam T The type of the returned invalid object. This type can either
   * be a scalar, or of type Tensor, SymmetricTensor, or DerivativeForm. Other
   * types may be supported if there is a corresponding specialization of the
   * internal::SignalingNaN::NaNInitializer class for this type.
   *
   * @note Because the type @p T is not used as a function argument, the
   * compiler cannot deduce it from the type of arguments. Consequently, you
   * have to provide it explicitly. For example, the line
   *   @code
   *     Tensor<1,dim> tensor = Utilities::signaling_nan<Tensor<1,dim> >();
   *   @endcode
   * initializes a tensor with invalid values.
   */
  template <class T>
  T
  signaling_nan()
  {
    // dispatch to the classes in the internal namespace because there
    // we can do partial specializations, which is not possible for
    // template functions such as the current one
    return internal::SignalingNaN::NaNInitializer<T>::invalid_element();
  }
} // namespace numbers


DEAL_II_NAMESPACE_CLOSE

#endif
