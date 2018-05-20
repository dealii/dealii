// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2017 by the deal.II authors
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
     * NaN objects. This is where the Utilities::signaling_nan() function
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

          for(unsigned int i = 0; i < dim; ++i)
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
          for(unsigned int i = 0; i < dim; ++i)
            nan_tensor[i]
              = NaNInitializer<Tensor<rank - 1, dim, T>>::invalid_element();

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

          for(unsigned int i = 0; i < dim; ++i)
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
          for(unsigned int i = 0;
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
          for(unsigned int i = 0; i < spacedim; ++i)
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
   * of such objects is to use them as markers for uninitialized objects and
   * arrays that are required to be filled in other places, and to trigger an
   * error when this later initialization does not happen before the first
   * use.
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
