// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2021 by the deal.II authors
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

#ifndef dealii_differentiation_sd_symengine_product_types_h
#define dealii_differentiation_sd_symengine_product_types_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE


#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/template_constraints.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>

#  include <boost/type_traits.hpp>

#  include <type_traits>


DEAL_II_NAMESPACE_OPEN


template <>
struct EnableIfScalar<Differentiation::SD::Expression>
{
  using type = Differentiation::SD::Expression;
};


namespace internal
{
  namespace SD
  {
    /**
     * A more general implementation of product types.
     * There are so many permutation of admissible operations
     * that getting the compiler to determine the valid
     * combinations using template metaprogramming makes
     * more sense than manually maintaining the list by
     * hand.
     *
     * This class is a workaround for issue of non-deduction
     * of types in template partial specializations that
     * would otherwise occur if trying to directly implement
     * these as specializations of the ProductTypeImpl class
     * itself.
     */
    template <typename T, typename U, typename V = void>
    struct GeneralProductTypeImpl;

    template <typename T>
    struct GeneralProductTypeImpl<
      T,
      Differentiation::SD::Expression,
      typename std::enable_if<std::is_arithmetic<T>::value>::type>
    {
      using type = Differentiation::SD::Expression;
    };

    template <typename T>
    struct GeneralProductTypeImpl<
      T,
      Differentiation::SD::Expression,
      typename std::enable_if<
        boost::is_complex<T>::value &&
        std::is_arithmetic<typename T::value_type>::value>::type>
    {
      using type = Differentiation::SD::Expression;
    };

    template <int rank, int dim, typename T>
    struct GeneralProductTypeImpl<Tensor<rank, dim, T>,
                                  Differentiation::SD::Expression>
    {
      using type =
        Tensor<rank,
               dim,
               typename ProductType<T, Differentiation::SD::Expression>::type>;
    };

    template <int rank, int dim, typename T>
    struct GeneralProductTypeImpl<SymmetricTensor<rank, dim, T>,
                                  Differentiation::SD::Expression>
    {
      using type = SymmetricTensor<
        rank,
        dim,
        typename ProductType<T, Differentiation::SD::Expression>::type>;
    };

  } // namespace SD


  template <>
  struct ProductTypeImpl<Differentiation::SD::Expression,
                         Differentiation::SD::Expression>
  {
    using type = Differentiation::SD::Expression;
  };


  template <typename T>
  struct ProductTypeImpl<T, Differentiation::SD::Expression>
  {
    using type = typename SD::
      GeneralProductTypeImpl<T, Differentiation::SD::Expression>::type;
  };

  template <typename T>
  struct ProductTypeImpl<Differentiation::SD::Expression, T>
  {
    using type = typename SD::
      GeneralProductTypeImpl<T, Differentiation::SD::Expression>::type;
  };

} // namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
