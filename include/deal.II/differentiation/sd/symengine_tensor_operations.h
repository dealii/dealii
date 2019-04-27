// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

#ifndef dealii_differentiation_sd_symengine_tensor_operations_h
#define dealii_differentiation_sd_symengine_tensor_operations_h

#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_SYMENGINE

#  include <deal.II/base/symmetric_tensor.h>
#  include <deal.II/base/tensor.h>

#  include <deal.II/differentiation/sd/symengine_number_types.h>
#  include <deal.II/differentiation/sd/symengine_types.h>

#  include <utility>
#  include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Differentiation
{
  namespace SD
  {
    /**
     * @name Symbolic variable creation
     */
    //@{

    /**
     * Return a vector of Expressions representing a vectorial symbolic
     * variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"v"` then the vectorial
     * symbolic variable that is returned represents the vector $v$. Each
     * component of $v$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component index.
     *
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the vector of returned
     *            symbolic variables.
     * @return A vector (a rank-1 tensor) of symbolic variables with the name of
     *         each individual component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbols(const std::string &symbol);

    /**
     * Return a tensor of Expressions representing a tensorial symbolic
     * variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"T"` then the tensorial
     * symbolic variable that is returned represents the vector $T$. Each
     * component of $T$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component indices.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic variables.
     * @return A tensor of symbolic variables with the name of each individual
     *         component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbols(const std::string &symbol);

    /**
     * Return a symmetric tensor of Expressions representing a tensorial
     * symbolic variable with the identifier specified by @p symbol.
     *
     * For example, if the @p symbol is the string `"S"` then the tensorial
     * symbolic variable that is returned represents the vector $S$. Each
     * component of $S$ is prefixed by the given @p symbol, and has a suffix that
     * indicates its component indices.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic variables.
     * @return A tensor of symbolic variables with the name of each individual
     *         component prefixed by @p symbol.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbols(const std::string &symbol);

    /**
     * Return a vector of Expression representing a vectorial symbolic
     * function with the identifier specified by @p symbol. The functions'
     * symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the vector of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A vector (a rank-1 tensor) of generic symbolic functions with the
     *         name of each individual component prefixed by @p symbol, a suffix
     *         that indicates its component index, and the number of input
     *         arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int dim>
    Tensor<1, dim, Expression>
    make_vector_of_symbolic_functions(const std::string &            symbol,
                                      const types::substitution_map &arguments);

    /**
     * Return a tensor of Expression representing a tensorial symbolic
     * function with the identifier specified by @p symbol. The functions'
     * symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A tensor of generic symbolic functions with the name of each
     *         individual component prefixed by @p symbol, a suffix that
     *         indicates its component indeices, and the number of input
     *         arguments equal to  the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    make_tensor_of_symbolic_functions(const std::string &            symbol,
                                      const types::substitution_map &arguments);

    /**
     * Return a symmetric tensor of Expression representing a tensorial
     * symbolic function with the identifier specified by @p symbol. The
     * functions' symbolic dependencies are specified by the keys to the input
     * @p arguments map; the values stored in the map are ignored.
     *
     * @tparam rank The rank of the returned tensor.
     * @tparam dim The dimension of the returned tensor.
     * @param[in] symbol An indentifier (or name) for the tensor of returned
     *            symbolic functions.
     * @param[in] arguments A map of input arguments to the returned
     *            symbolic functions.
     * @return A symmetric tensor of generic symbolic functions with the name of
     *         each individual component prefixed by @p symbol, a suffix that
     *         indicates its component indeices, and the number of input
     *         arguments equal to the length of @p arguments.
     *
     * @warning It is up to the user to ensure that there is no ambiguity in the
     * symbols used within a section of code.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    make_symmetric_tensor_of_symbolic_functions(
      const std::string &            symbol,
      const types::substitution_map &arguments);

    //@}

    /**
     * @name Symbolic differentiation
     */
    //@{

    /**
     * Return the symbolic result of computing the partial derivative of the
     * scalar @p f with respect to the tensor @p T.
     *
     * @param[in] f A scalar symbolic function or (dependent) expression.
     * @param[in] T A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{T}}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Expression &f, const Tensor<rank, dim, Expression> &T);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * scalar @p f with respect to the symmetric tensor @p S.
     *
     * @param[in] f A scalar symbolic function or (dependent) expression.
     * @param[in] S A symmetric tensor of symbolic (independent) variables.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{S}}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Expression &                            f,
                  const SymmetricTensor<rank, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * rank-0 tensor (or scalar) @p f with respect to the tensor @p T.
     *
     * @param[in] f A rank-0 tensor symbolic function or (dependent) expression.
     * @param[in] T A tensor of symbolic (independent) variables.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{T}}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &   f,
                  const Tensor<rank, dim, Expression> &T);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * rank-0 tensor (or scalar) @p f with respect to the symmetric tensor @p S.
     *
     * @param[in] f A rank-0 tensor symbolic function or (dependent) expression.
     * @param[in] S A symmetric tensor of symbolic (independent) variables.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial f}{\partial \mathbf{S}}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &            f,
                  const SymmetricTensor<rank, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the scalar @p x.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{T}}{\partial x}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &T, const Expression &x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the scalar @p x.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] x A scalar symbolic (independent) variable.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{S}}{\partial x}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &S,
                  const Expression &                            x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the rank-0 tensor @p x.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] x A rank-0 tensor symbolic symbolic (independent) variable.
     * @return The tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{T}}{\partial x}$.
     */
    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &T,
                  const Tensor<0, dim, Expression> &   x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the rank-0 tensor @p x.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] x A rank-0 tensor symbolic symbolic (independent) variable.
     * @return The symmetric tensor of symbolic functions or expressions representing
     *         the result $\frac{\partial \mathbf{S}}{\partial x}$.
     */
    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &S,
                  const Tensor<0, dim, Expression> &            x);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T1 with respect to the tensor @p T2.
     *
     * @param[in] T1 A tensor of symbolic functions or (dependent) expressions.
     * @param[in] T2 A tensor symbolic symbolic (independent) variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{T}_{1}}{\partial
     * \mathbf{T}_{2}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &T1,
                  const Tensor<rank_2, dim, Expression> &T2);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S1 with respect to the symmetric tensor @p S2.
     *
     * @param[in] S1 A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] S2 A symmetric tensor symbolic symbolic (independent)
     *            variables.
     * @return The symmetric tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{S}_{1}}{\partial
     * \mathbf{S}_{2}}$.
     */
    template <int rank_1, int rank_2, int dim>
    SymmetricTensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &S1,
                  const SymmetricTensor<rank_2, dim, Expression> &S2);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * tensor @p T with respect to the symmetric tensor @p S.
     *
     * @param[in] T A tensor of symbolic functions or (dependent) expressions.
     * @param[in] S A symmetric tensor symbolic symbolic (independent)
     *            variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{T}}{\partial \mathbf{S}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &         T,
                  const SymmetricTensor<rank_2, dim, Expression> &S);

    /**
     * Return the symbolic result of computing the partial derivative of the
     * symmetric tensor @p S with respect to the tensor @p T.
     *
     * @param[in] S A symmetric tensor of symbolic functions or (dependent)
     *            expressions.
     * @param[in] T A tensor symbolic symbolic (independent) variables.
     * @return The tensor of symbolic functions or variables representing
     *         the result $\frac{\partial \mathbf{S}}{\partial \mathbf{T}}$.
     */
    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &S,
                  const Tensor<rank_2, dim, Expression> &         T);

    //@}

  } // namespace SD
} // namespace Differentiation


/* -------------------- inline and template functions ------------------ */


#  ifndef DOXYGEN

namespace Differentiation
{
  namespace SD
  {
    /* ---------------- Symbolic differentiation --------------*/


    namespace internal
    {
      template <int rank_1, int rank_2>
      TableIndices<rank_1 + rank_2>
      concatenate_indices(const TableIndices<rank_1> &indices_1,
                          const TableIndices<rank_2> &indices_2)
      {
        TableIndices<rank_1 + rank_2> indices_out;
        for (unsigned int i = 0; i < rank_1; ++i)
          indices_out[i] = indices_1[i];
        for (unsigned int j = 0; j < rank_2; ++j)
          indices_out[rank_1 + j] = indices_2[j];
        return indices_out;
      }


      template <int rank>
      TableIndices<rank>
      transpose_indices(const TableIndices<rank> &indices)
      {
        return indices;
      }


      template <>
      inline TableIndices<2>
      transpose_indices(const TableIndices<2> &indices)
      {
        return TableIndices<2>(indices[1], indices[0]);
      }


      template <int rank, int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<rank> &,
                             const Tensor<rank, dim, ValueType> &)
      {
        return false;
      }


      template <int rank, int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<rank> &,
                             const SymmetricTensor<rank, dim, ValueType> &)
      {
        static_assert(
          rank == 0 || rank == 2,
          "Querying symmetric component for non rank-2 symmetric tensor index is not allowed.");
        return false;
      }


      template <int dim, typename ValueType>
      bool
      is_symmetric_component(const TableIndices<2> &table_indices,
                             const SymmetricTensor<2, dim, ValueType> &)
      {
        return table_indices[0] != table_indices[1];
      }


      template <int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<0, dim, ValueType>
      scalar_diff_tensor(const ValueType &                    func,
                         const TensorType<0, dim, ValueType> &op)
      {
        return differentiate(func, op);
      }


      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      scalar_diff_tensor(const ValueType &                       func,
                         const TensorType<rank, dim, ValueType> &op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(func, op[indices]);

            if (is_symmetric_component(indices, op))
              out[indices] *= 0.5;
          }
        return out;
      }


      // Specialization for rank-0 tensor
      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_tensor(const TensorType<0, dim, ValueType> &   func,
                         const TensorType<rank, dim, ValueType> &op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(func, op[indices]);

            if (is_symmetric_component(indices, op))
              out[indices] *= 0.5;
          }
        return out;
      }


      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_scalar(const TensorType<rank, dim, ValueType> &funcs,
                         const ValueType &                       op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(funcs[indices], op);
          }
        return out;
      }


      // Specialization for rank-0 tensor
      template <int rank,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank, dim, ValueType>
      tensor_diff_tensor(const TensorType<rank, dim, ValueType> &funcs,
                         const TensorType<0, dim, ValueType> &   op)
      {
        TensorType<rank, dim, ValueType> out;
        for (unsigned int i = 0; i < out.n_independent_components; ++i)
          {
            const TableIndices<rank> indices(
              out.unrolled_to_component_indices(i));
            out[indices] = differentiate(funcs[indices], op);
          }
        return out;
      }


      // For either symmetric or normal tensors
      template <int rank_1,
                int rank_2,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType>
      TensorType<rank_1 + rank_2, dim, ValueType>
      tensor_diff_tensor(const TensorType<rank_1, dim, ValueType> &funcs,
                         const TensorType<rank_2, dim, ValueType> &op)
      {
        TensorType<rank_1 + rank_2, dim, ValueType> out;
        for (unsigned int i = 0; i < funcs.n_independent_components; ++i)
          {
            const TableIndices<rank_1> indices_i(
              funcs.unrolled_to_component_indices(i));
            for (unsigned int j = 0; j < op.n_independent_components; ++j)
              {
                const TableIndices<rank_2> indices_j(
                  op.unrolled_to_component_indices(j));
                const TableIndices<rank_1 + rank_2> indices_out =
                  concatenate_indices(indices_i, indices_j);

                out[indices_out] =
                  differentiate(funcs[indices_i], op[indices_j]);

                if (is_symmetric_component(indices_j, op))
                  out[indices_out] *= 0.5;
              }
          }
        return out;
      }


      // For mixed symmetric/standard tensors
      // The return type is always a standard tensor, since we cannot be sure
      // that any symmetries exist in either the function tensor or the
      // differential operator.
      template <int rank_1,
                int rank_2,
                int dim,
                typename ValueType = Expression,
                template <int, int, typename> class TensorType_1,
                template <int, int, typename> class TensorType_2>
      Tensor<rank_1 + rank_2, dim, ValueType>
      tensor_diff_tensor(const TensorType_1<rank_1, dim, ValueType> &funcs,
                         const TensorType_2<rank_2, dim, ValueType> &op)
      {
        Tensor<rank_1 + rank_2, dim, ValueType> out;
        for (unsigned int i = 0; i < funcs.n_independent_components; ++i)
          {
            const TableIndices<rank_1> indices_i(
              funcs.unrolled_to_component_indices(i));
            for (unsigned int j = 0; j < op.n_independent_components; ++j)
              {
                const TableIndices<rank_2> indices_j(
                  op.unrolled_to_component_indices(j));
                const TableIndices<rank_1 + rank_2> indices_out =
                  concatenate_indices(indices_i, indices_j);

                out[indices_out] =
                  differentiate(funcs[indices_i], op[indices_j]);

                if (is_symmetric_component(indices_j, op))
                  out[indices_out] *= 0.5;

                // TODO: Implement for SymmetricTensor<4,dim,...>
                if (std::is_same<TensorType_1<rank_1, dim, ValueType>,
                                 SymmetricTensor<2, dim, ValueType>>::
                      value) // Symmetric function
                  {
                    const TableIndices<rank_1 + rank_2> indices_out_t =
                      concatenate_indices(transpose_indices(indices_i),
                                          indices_j);
                    out[indices_out_t] = out[indices_out];
                  }
                else if (std::is_same<TensorType_2<rank_2, dim, ValueType>,
                                      SymmetricTensor<2, dim, ValueType>>::
                           value) // Symmetric operator
                  {
                    const TableIndices<rank_1 + rank_2> indices_out_t =
                      concatenate_indices(indices_i,
                                          transpose_indices(indices_j));
                    out[indices_out_t] = out[indices_out];
                  }
                else
                  {
                    Assert(
                      false,
                      ExcMessage(
                        "Expect mixed tensor differentiation to have at least "
                        "one component stemming from a symmetric tensor."));
                  }
              }
          }
        return out;
      }

    } // namespace internal


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Expression &                   func,
                  const Tensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Expression &                            func,
                  const SymmetricTensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &   func,
                  const Tensor<rank, dim, Expression> &op)
    {
      return internal::scalar_diff_tensor(func, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const Tensor<0, dim, Expression> &            func,
                  const SymmetricTensor<rank, dim, Expression> &op)
    {
      // Ensure that this returns a symmetric tensor by
      // invoking the scalar value function
      const Expression tmp = func;
      return internal::scalar_diff_tensor(tmp, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &symbol_tensor,
                  const Expression &                   op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    Tensor<rank, dim, Expression>
    differentiate(const Tensor<rank, dim, Expression> &symbol_tensor,
                  const Tensor<0, dim, Expression> &   op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &symbol_tensor,
                  const Expression &                            op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank, int dim>
    SymmetricTensor<rank, dim, Expression>
    differentiate(const SymmetricTensor<rank, dim, Expression> &symbol_tensor,
                  const Tensor<0, dim, Expression> &            op)
    {
      return internal::tensor_diff_scalar(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &symbol_tensor,
                  const Tensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    SymmetricTensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &symbol_tensor,
                  const SymmetricTensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const Tensor<rank_1, dim, Expression> &         symbol_tensor,
                  const SymmetricTensor<rank_2, dim, Expression> &op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }


    template <int rank_1, int rank_2, int dim>
    Tensor<rank_1 + rank_2, dim, Expression>
    differentiate(const SymmetricTensor<rank_1, dim, Expression> &symbol_tensor,
                  const Tensor<rank_2, dim, Expression> &         op)
    {
      return internal::tensor_diff_tensor(symbol_tensor, op);
    }

  } // namespace SD
} // namespace Differentiation

#  endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_SYMENGINE

#endif
