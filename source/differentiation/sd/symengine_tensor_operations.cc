// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>
#  include <deal.II/base/utilities.h>

#  include <deal.II/differentiation/sd/symengine_product_types.h>
#  include <deal.II/differentiation/sd/symengine_scalar_operations.h>
#  include <deal.II/differentiation/sd/symengine_tensor_operations.h>
#  include <deal.II/differentiation/sd/symengine_types.h>
#  include <deal.II/differentiation/sd/symengine_utilities.h>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    // --- Symbolic variable creation ---

    namespace internal
    {
      namespace
      {
        template <int rank>
        std::string
        make_index_string(const TableIndices<rank> &indices)
        {
          std::string out;
          for (unsigned int i = 0; i < rank; ++i)
            out += std::to_string(indices[i]);
          return out;
        }


        template <int rank, int dim>
        struct Symbol_Tensor;


        template <int rank, int dim>
        struct Symbol_Tensor
        {
          static Tensor<rank, dim, Expression>
          create(const std::string &sym)
          {
            Tensor<rank, dim, Expression> out;
            for (unsigned int i = 0; i < out.n_independent_components; ++i)
              {
                const TableIndices<rank> indices(
                  out.unrolled_to_component_indices(i));
                out[indices] =
                  Expression(sym + "_" + internal::make_index_string(indices));
              }
            return out;
          }
        };


        template <int dim>
        struct Symbol_Tensor<0, dim>
        {
          static Tensor<0, dim, Expression>
          create(const std::string &sym)
          {
            return make_symbol(sym);
          }
        };


        template <int rank, int dim>
        struct Symbol_SymmetricTensor;


        template <int rank, int dim>
        struct Symbol_SymmetricTensor
        {
          static SymmetricTensor<2, dim, Expression>
          create(const std::string &sym)
          {
            SymmetricTensor<rank, dim, Expression> out;
            for (unsigned int i = 0; i < out.n_independent_components; ++i)
              {
                const TableIndices<rank> indices(
                  out.unrolled_to_component_indices(i));
                out[indices] =
                  Expression(sym + "_" + internal::make_index_string(indices));
              }
            return out;
          }
        };


        template <int dim>
        struct Symbol_SymmetricTensor<4, dim>
        {
          static SymmetricTensor<4, dim, Expression>
          create(const std::string &sym)
          {
            SymmetricTensor<4, dim, Expression> out;
            for (unsigned int i = 0;
                 i < SymmetricTensor<2, dim>::n_independent_components;
                 ++i)
              for (unsigned int j = 0;
                   j < SymmetricTensor<2, dim>::n_independent_components;
                   ++j)
                {
                  const TableIndices<4> indices =
                    make_rank_4_tensor_indices<dim>(i, j);
                  out[indices] = Expression(
                    sym + "_" + internal::make_index_string(indices));
                }
            return out;
          }
        };


        template <int rank, int dim>
        struct Symbol_Function_Tensor;


        template <int rank, int dim>
        struct Symbol_Function_Tensor
        {
          static Tensor<rank, dim, Expression>
          create(const std::string             &sym,
                 const types::substitution_map &arguments)
          {
            Tensor<rank, dim, Expression> out;
            const types::symbol_vector    args =
              Utilities::extract_symbols(arguments);
            for (unsigned int i = 0; i < out.n_independent_components; ++i)
              {
                const TableIndices<rank> indices(
                  out.unrolled_to_component_indices(i));
                out[indices] =
                  Expression(sym + "_" + internal::make_index_string(indices),
                             args);
              }
            return out;
          }
        };


        template <int dim>
        struct Symbol_Function_Tensor<0, dim>
        {
          static Tensor<0, dim, Expression>
          create(const std::string             &sym,
                 const types::substitution_map &arguments)
          {
            return make_symbolic_function(sym, arguments);
          }
        };


        template <int rank, int dim>
        struct Symbol_Function_SymmetricTensor;


        template <int rank, int dim>
        struct Symbol_Function_SymmetricTensor
        {
          static SymmetricTensor<2, dim, Expression>
          create(const std::string             &sym,
                 const types::substitution_map &arguments)
          {
            SymmetricTensor<rank, dim, Expression> out;
            const types::symbol_vector             args =
              Utilities::extract_symbols(arguments);
            for (unsigned int i = 0; i < out.n_independent_components; ++i)
              {
                const TableIndices<rank> indices(
                  out.unrolled_to_component_indices(i));
                out[indices] =
                  Expression(sym + "_" + internal::make_index_string(indices),
                             args);
              }
            return out;
          }
        };


        template <int dim>
        struct Symbol_Function_SymmetricTensor<4, dim>
        {
          static SymmetricTensor<4, dim, Expression>
          create(const std::string             &sym,
                 const types::substitution_map &arguments)
          {
            SymmetricTensor<4, dim, Expression> out;
            const types::symbol_vector          args =
              Utilities::extract_symbols(arguments);
            for (unsigned int i = 0;
                 i < SymmetricTensor<2, dim>::n_independent_components;
                 ++i)
              for (unsigned int j = 0;
                   j < SymmetricTensor<2, dim>::n_independent_components;
                   ++j)
                {
                  const TableIndices<4> indices =
                    make_rank_4_tensor_indices<dim>(i, j);
                  out[indices] =
                    Expression(sym + "_" + internal::make_index_string(indices),
                               args);
                }
            return out;
          }
        };
      } // namespace
    }   // namespace internal


    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbols(const std::string &sym)
    {
      return internal::Symbol_Tensor<1, dim>::create(sym);
    }


#  ifndef DOXYGEN
    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbolic_functions(const std::string             &sym,
                                      const types::substitution_map &arguments)
    {
      return internal::Symbol_Function_Tensor<1, dim>::create(sym, arguments);
    }
#  endif


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbols(const std::string &sym)
    {
      return internal::Symbol_Tensor<rank, dim>::create(sym);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbols(const std::string &sym)
    {
      return internal::Symbol_SymmetricTensor<rank, dim>::create(sym);
    }


#  ifndef DOXYGEN
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbolic_functions(const std::string             &sym,
                                      const types::substitution_map &arguments)
    {
      return internal::Symbol_Function_Tensor<rank, dim>::create(sym,
                                                                 arguments);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbolic_functions(
      const std::string             &sym,
      const types::substitution_map &arguments)
    {
      return internal::Symbol_Function_SymmetricTensor<rank, dim>::create(
        sym, arguments);
    }
#  endif // DOXYGEN

  } // namespace SD
} // namespace Differentiation

/* --- Explicit instantiations --- */

#  include "differentiation/sd/symengine_tensor_operations.inst"

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE
