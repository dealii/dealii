// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2019 by the deal.II authors
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

#ifndef dealii_tensor_accessors_h
#define dealii_tensor_accessors_h

#include <deal.II/base/config.h>

#include <deal.II/base/table_indices.h>
#include <deal.II/base/template_constraints.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This namespace is a collection of algorithms working on generic tensorial
 * objects (of arbitrary rank).
 *
 * The rationale to implement such functionality in a generic fashion in a
 * separate namespace is
 *  - to easy code reusability and therefore avoid code duplication.
 *  - to have a well-defined interface that allows to exchange the low
 * level implementation.
 *
 *
 * A tensorial object has the notion of a rank and allows a rank-times
 * recursive application of the index operator, e.g., if <code>t</code> is a
 * tensorial object of rank 4, the following access is valid:
 * @code
 *   t[1][2][1][4]
 * @endcode
 *
 * deal.II has its own implementation for tensorial objects such as
 * dealii::Tensor<rank, dim, Number> and dealii::SymmetricTensor<rank, dim,
 * Number>
 *
 * The methods and algorithms implemented in this namespace, however, are
 * fully generic. More precisely, it can operate on nested c-style arrays, or
 * on class types <code>T</code> with a minimal interface that provides a
 * local alias <code>value_type</code> and an index operator
 * <code>operator[](unsigned int)</code> that returns a (const or non-const)
 * reference of <code>value_type</code>:
 * @code
 *   template <...>
 *   class T
 *   {
 *     using value_type = ...;
 *     value_type & operator[](unsigned int);
 *     const value_type & operator[](unsigned int) const;
 *   };
 * @endcode
 *
 * This namespace provides primitives for access, reordering and contraction
 * of such objects.
 *
 * @ingroup geomprimitives
 */
namespace TensorAccessors
{
  // forward declarations
  namespace internal
  {
    template <int index, int rank, typename T>
    class ReorderedIndexView;
    template <int position, int rank>
    struct ExtractHelper;
    template <int no_contr, int rank_1, int rank_2, int dim>
    class Contract;
    template <int rank_1, int rank_2, int dim>
    class Contract3;
  } // namespace internal


  /**
   * This class provides a local alias @p value_type denoting the resulting
   * type of an access with operator[](unsigned int). More precisely, @p
   * value_type will be
   *  - <code>T::value_type</code> if T is a tensorial class providing an
   * alias <code>value_type</code> and does not have a const qualifier.
   *  - <code>const T::value_type</code> if T is a tensorial class
   * providing an alias <code>value_type</code> and does have a const
   * qualifier.
   *  - <code>const T::value_type</code> if T is a tensorial class
   * providing an alias <code>value_type</code> and does have a const
   * qualifier.
   *  - <code>A</code> if T is of array type <code>A[...]</code>
   *  - <code>const A</code> if T is of array type <code>A[...]</code> and
   * does have a const qualifier.
   */
  template <typename T>
  struct ValueType
  {
    using value_type = typename T::value_type;
  };

  template <typename T>
  struct ValueType<const T>
  {
    using value_type = const typename T::value_type;
  };

  template <typename T, std::size_t N>
  struct ValueType<T[N]>
  {
    using value_type = T;
  };

  template <typename T, std::size_t N>
  struct ValueType<const T[N]>
  {
    using value_type = const T;
  };


  /**
   * This class provides a local alias @p value_type that is equal to the
   * alias <code>value_type</code> after @p deref_steps recursive
   * dereferences via ```operator[](unsigned int)```. Further, constness is
   * preserved via the ValueType type trait, i.e., if T is const,
   * ReturnType<rank, T>::value_type will also be const.
   */
  template <int deref_steps, typename T>
  struct ReturnType
  {
    using value_type =
      typename ReturnType<deref_steps - 1,
                          typename ValueType<T>::value_type>::value_type;
  };

  template <typename T>
  struct ReturnType<0, T>
  {
    using value_type = T;
  };


  /**
   * Provide a "tensorial view" to a reference @p t of a tensor object of rank
   * @p rank in which the index @p index is shifted to the end. As an example
   * consider a tensor of 5th order in dim=5 space dimensions that can be
   * accessed through 5 recursive <code>operator[]()</code> invocations:
   * @code
   *   Tensor<5, dim> tensor;
   *   tensor[0][1][2][3][4] = 42.;
   * @endcode
   * Index 1 (the 2nd index, count starts at 0) can now be shifted to the end
   * via
   * @code
   *   auto tensor_view = reordered_index_view<1, 5>(tensor);
   *   tensor_view[0][2][3][4][1] == 42.; // is true
   * @endcode
   * The usage of the dealii::Tensor type was solely for the sake of an
   * example. The mechanism implemented by this function is available for
   * fairly general tensorial types @p T.
   *
   * The purpose of this reordering facility is to be able to contract over an
   * arbitrary index of two (or more) tensors:
   *  - reorder the indices in mind to the end of the tensors
   *  - use the contract function below that contracts the _last_ elements of
   * tensors.
   *
   * @note This function returns an internal class object consisting of an
   * array subscript operator <code>operator[](unsigned int)</code> and an
   * alias <code>value_type</code> describing its return value.
   *
   * @tparam index The index to be shifted to the end. Indices are counted
   * from 0, thus the valid range is $0\le\text{index}<\text{rank}$.
   * @tparam rank Rank of the tensorial object @p t
   * @tparam T A tensorial object of rank @p rank. @p T must provide a local
   * alias <code>value_type</code> and an index operator
   * <code>operator[]()</code> that returns a (const or non-const) reference
   * of <code>value_type</code>.
   */
  template <int index, int rank, typename T>
  constexpr DEAL_II_ALWAYS_INLINE internal::ReorderedIndexView<index, rank, T>
                                  reordered_index_view(T &t)
  {
    static_assert(0 <= index && index < rank,
                  "The specified index must lie within the range [0,rank)");

    return internal::ReorderedIndexView<index, rank, T>(t);
  }


  /**
   * Return a reference (const or non-const) to a subobject of a tensorial
   * object @p t of type @p T, as described by an array type @p ArrayType
   * object @p indices. For example:
   * @code
   *   Tensor<5, dim> tensor;
   *   TableIndices<5> indices (0, 1, 2, 3, 4);
   *   TensorAccessors::extract(tensor, indices) = 42;
   * @endcode
   * This is equivalent to <code>tensor[0][1][2][3][4] = 42.</code>.
   *
   * @tparam T A tensorial object of rank @p rank. @p T must provide a local
   * alias <code>value_type</code> and an index operator
   * <code>operator[]()</code> that returns a (const or non-const) reference
   * of <code>value_type</code>. Further, its tensorial rank must be equal or
   * greater than @p rank.
   *
   * @tparam ArrayType An array like object, such as std::array, or
   * dealii::TableIndices  that stores at least @p rank indices that can be
   * accessed via operator[]().
   */
  template <int rank, typename T, typename ArrayType>
  constexpr DEAL_II_ALWAYS_INLINE typename ReturnType<rank, T>::value_type &
  extract(T &t, const ArrayType &indices)
  {
    return internal::ExtractHelper<0, rank>::template extract<T, ArrayType>(
      t, indices);
  }


  /**
   * This function contracts two tensorial objects @p left and @p right and
   * stores the result in @p result. The contraction is done over the _last_
   * @p no_contr indices of both tensorial objects:
   *
   * @f[
   *   \text{result}_{i_1,..,i_{r1},j_1,..,j_{r2}}
   *   = \sum_{k_1,..,k_{\mathrm{no\_contr}}}
   *     \mathrm{left}_{i_1,..,i_{r1},k_1,..,k_{\mathrm{no\_contr}}}
   *     \mathrm{right}_{j_1,..,j_{r2},k_1,..,k_{\mathrm{no\_contr}}}
   * @f]
   *
   * Calling this function is equivalent of writing the following low level
   * code:
   * @code
   *   for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
   *     ...
   *       for(unsigned int i_ = 0; i_ < dim; ++i_)
   *         for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
   *           ...
   *             for(unsigned int j_ = 0; j_ < dim; ++j_)
   *               {
   *                 result[i_0]..[i_][j_0]..[j_] = 0.;
   *                 for(unsigned int k_0 = 0; k_0 < dim; ++k_0)
   *                   ...
   *                     for(unsigned int k_ = 0; k_ < dim; ++k_)
   *                       result[i_0]..[i_][j_0]..[j_] +=
   *                         left[i_0]..[i_][k_0]..[k_]
   *                           * right[j_0]..[j_][k_0]..[k_];
   *               }
   * @endcode
   * with r = rank_1 + rank_2 - 2 * no_contr, l = rank_1 - no_contr, l1 =
   * rank_1, and c = no_contr.
   *
   * @note The Types @p T1, @p T2, and @p T3 must have rank rank_1 + rank_2 -
   * 2 * no_contr, rank_1, or rank_2, respectively. Obviously, no_contr must
   * be less or equal than rank_1 and rank_2.
   */
  template <int no_contr,
            int rank_1,
            int rank_2,
            int dim,
            typename T1,
            typename T2,
            typename T3>
  constexpr inline DEAL_II_ALWAYS_INLINE void
  contract(T1 &result, const T2 &left, const T3 &right)
  {
    static_assert(rank_1 >= no_contr,
                  "The rank of the left tensor must be "
                  "equal or greater than the number of "
                  "contractions");
    static_assert(rank_2 >= no_contr,
                  "The rank of the right tensor must be "
                  "equal or greater than the number of "
                  "contractions");

    internal::Contract<no_contr, rank_1, rank_2, dim>::
      template contract<T1, T2, T3>(result, left, right);
  }


  /**
   * Full contraction of three tensorial objects:
   *
   * @f[
   *   \sum_{i_1,..,i_{r1},j_1,..,j_{r2}}
   *   \text{left}_{i_1,..,i_{r1}}
   *   \text{middle}_{i_1,..,i_{r1},j_1,..,j_{r2}}
   *   \text{right}_{j_1,..,j_{r2}}
   * @f]
   *
   * Calling this function is equivalent of writing the following low level
   * code:
   * @code
   *   T1 result = T1();
   *   for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
   *     ...
   *       for(unsigned int i_ = 0; i_ < dim; ++i_)
   *         for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
   *           ...
   *             for(unsigned int j_ = 0; j_ < dim; ++j_)
   *               result += left[i_0]..[i_]
   *                           * middle[i_0]..[i_][j_0]..[j_]
   *                           * right[j_0]..[j_];
   * @endcode
   *
   * @note The Types @p T2, @p T3, and @p T4 must have rank rank_1, rank_1 +
   * rank_2, and rank_3, respectively. @p T1 must be a scalar type.
   */
  template <int rank_1,
            int rank_2,
            int dim,
            typename T1,
            typename T2,
            typename T3,
            typename T4>
  constexpr T1
  contract3(const T2 &left, const T3 &middle, const T4 &right)
  {
    return internal::Contract3<rank_1, rank_2, dim>::
      template contract3<T1, T2, T3, T4>(left, middle, right);
  }


  namespace internal
  {
    // -------------------------------------------------------------------------
    // Forward declarations and type traits
    // -------------------------------------------------------------------------

    template <int rank, typename S>
    class StoreIndex;
    template <typename T>
    class Identity;
    template <int no_contr, int dim>
    class Contract2;

    /**
     * An internally used type trait to allow nested application of the
     * function reordered_index_view(T &t).
     *
     * The problem is that when working with the actual tensorial types, we
     * have to return subtensors by reference - but sometimes, especially for
     * StoreIndex and ReorderedIndexView that return rvalues, we have to
     * return by value.
     */
    template <typename T>
    struct ReferenceType
    {
      using type = T &;
    };

    template <int rank, typename S>
    struct ReferenceType<StoreIndex<rank, S>>
    {
      using type = StoreIndex<rank, S>;
    };

    template <int index, int rank, typename T>
    struct ReferenceType<ReorderedIndexView<index, rank, T>>
    {
      using type = ReorderedIndexView<index, rank, T>;
    };


    // TODO: Is there a possibility to just have the following block of
    // explanation on an internal page in doxygen? If, yes. Doxygen
    // wizards, your call!

    // -------------------------------------------------------------------------
    // Implementation of helper classes for reordered_index_view
    // -------------------------------------------------------------------------

    // OK. This is utterly brutal template magic. Therefore, we will not
    // comment on the individual internal helper classes, because this is
    // of not much value, but explain the general recursion procedure.
    //
    // (In order of appearance)
    //
    // Our task is to reorder access to a tensor object where a specified
    // index is moved to the end. Thus we want to construct an object
    // <code>reordered</code> out of a <code>tensor</code> where the
    // following access patterns are equivalent:
    // @code
    //   tensor    [i_0]...[i_index-1][i_index][i_index+1]...[i_n]
    //   reordered [i_0]...[i_index_1][i_index+1]...[i_n][i_index]
    // @endcode
    //
    // The first task is to get rid of the application of
    // [i_0]...[i_index-1]. This is a classical recursion pattern - relay
    // the task from <index, rank> to <index-1, rank-1> by accessing the
    // subtensor object:

    template <int index, int rank, typename T>
    class ReorderedIndexView
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type = ReorderedIndexView<index - 1,
                                            rank - 1,
                                            typename ValueType<T>::value_type>;

      // Recurse by applying index j directly:
      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(t_[j]);
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // At some point we hit the condition index == 0 and rank > 1, i.e.,
    // the first index should be reordered to the end.
    //
    // At this point we cannot be lazy any more and have to start storing
    // indices because we get them in the wrong order. The user supplies
    //   [i_0][i_1]...[i_{rank - 1}]
    // but we have to call the subtensor object with
    //   [i_{rank - 1}[i_0][i_1]...[i_{rank-2}]
    //
    // So give up and relay the task to the StoreIndex class:

    template <int rank, typename T>
    class ReorderedIndexView<0, rank, T>
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type = StoreIndex<rank - 1, internal::Identity<T>>;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(Identity<T>(t_), j);
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // Sometimes, we're lucky and don't have to do anything. In this case
    // just return the original tensor.

    template <typename T>
    class ReorderedIndexView<0, 1, T>
    {
    public:
      constexpr ReorderedIndexView(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using value_type =
        typename ReferenceType<typename ValueType<T>::value_type>::type;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return t_[j];
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // Here, Identity is a helper class to ground the recursion in
    // StoreIndex. Its implementation is easy - we haven't stored any
    // indices yet. So, we just provide a function apply that returns the
    // application of an index j to the stored tensor t_:

    template <typename T>
    class Identity
    {
    public:
      constexpr Identity(typename ReferenceType<T>::type t)
        : t_(t)
      {}

      using return_type = typename ValueType<T>::value_type;

      constexpr DEAL_II_ALWAYS_INLINE typename ReferenceType<return_type>::type
      apply(unsigned int j) const
      {
        return t_[j];
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // StoreIndex is a class that stores an index recursively with every
    // invocation of operator[](unsigned int j): We do this by recursively
    // creating a new StoreIndex class of lower rank that stores the
    // supplied index j and holds a copy of the current class (with all
    // other stored indices). Again, we provide an apply member function
    // that knows how to apply an index on the highest rank and all
    // subsequently stored indices:

    template <int rank, typename S>
    class StoreIndex
    {
    public:
      constexpr StoreIndex(S s, int i)
        : s_(s)
        , i_(i)
      {}

      using value_type = StoreIndex<rank - 1, StoreIndex<rank, S>>;

      constexpr DEAL_II_ALWAYS_INLINE value_type
                                      operator[](unsigned int j) const
      {
        return value_type(*this, j);
      }

      using return_type =
        typename ValueType<typename S::return_type>::value_type;

      constexpr typename ReferenceType<return_type>::type
      apply(unsigned int j) const
      {
        return s_.apply(j)[i_];
      }

    private:
      const S   s_;
      const int i_;
    };

    // We have to store indices until we hit rank == 1. Then, upon the next
    // invocation of operator[](unsigned int j) we have all necessary
    // information available to return the actual object.

    template <typename S>
    class StoreIndex<1, S>
    {
    public:
      constexpr StoreIndex(S s, int i)
        : s_(s)
        , i_(i)
      {}

      using return_type =
        typename ValueType<typename S::return_type>::value_type;
      using value_type = return_type;

      constexpr DEAL_II_ALWAYS_INLINE return_type &
                                      operator[](unsigned int j) const
      {
        return s_.apply(j)[i_];
      }

    private:
      const S   s_;
      const int i_;
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for extract
    // -------------------------------------------------------------------------

    // Straightforward recursion implemented by specializing ExtractHelper
    // for position == rank. We use the type trait ReturnType<rank, T> to
    // have an idea what the final type will be.
    template <int position, int rank>
    struct ExtractHelper
    {
      template <typename T, typename ArrayType>
      constexpr static typename ReturnType<rank - position, T>::value_type &
      extract(T &t, const ArrayType &indices)
      {
        return ExtractHelper<position + 1, rank>::template extract<
          typename ValueType<T>::value_type,
          ArrayType>(t[indices[position]], indices);
      }
    };

    // For position == rank there is nothing to extract, just return the
    // object.
    template <int rank>
    struct ExtractHelper<rank, rank>
    {
      template <typename T, typename ArrayType>
      constexpr static T &
      extract(T &t, const ArrayType &)
      {
        return t;
      }
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for contract
    // -------------------------------------------------------------------------

    // Straightforward recursive pattern:
    //
    // As long as rank_1 > no_contr, assign indices from the left tensor to
    // result. This builds up the first part of the nested outer loops:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       [...]
    //         result[i_0]..[i_] ... left[i_0]..[i_] ...

    template <int no_contr, int rank_1, int rank_2, int dim>
    class Contract
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        for (unsigned int i = 0; i < dim; ++i)
          Contract<no_contr, rank_1 - 1, rank_2, dim>::contract(result[i],
                                                                left[i],
                                                                right);
      }
    };

    // If rank_1 == no_contr leave out the remaining no_contr indices for
    // the contraction and assign indices from the right tensor to the
    // result. This builds up the second part of the nested loops:
    //
    //  for(unsigned int i_0 = 0; i_0 < dim; ++i_0)
    //    ...
    //      for(unsigned int i_ = 0; i_ < dim; ++i_)
    //        for(unsigned int j_0 = 0; j_0 < dim; ++j_0)
    //          ...
    //            for(unsigned int j_ = 0; j_ < dim; ++j_)
    //             [...]
    //               result[i_0]..[i_][j_0]..[j_] ... left[i_0]..[i_] ...
    //               right[j_0]..[j_]
    //

    template <int no_contr, int rank_2, int dim>
    class Contract<no_contr, no_contr, rank_2, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        for (unsigned int i = 0; i < dim; ++i)
          Contract<no_contr, no_contr, rank_2 - 1, dim>::contract(result[i],
                                                                  left,
                                                                  right[i]);
      }
    };

    // If rank_1 == rank_2 == no_contr we have built up all of the outer
    // loop. Now, it is time to do the actual contraction:
    //
    // [...]
    //   {
    //     result[i_0]..[i_][j_0]..[j_] = 0.;
    //     for(unsigned int k_0 = 0; k_0 < dim; ++k_0)
    //       ...
    //         for(unsigned int k_ = 0; k_ < dim; ++k_)
    //           result[i_0]..[i_][j_0]..[j_] += left[i_0]..[i_][k_0]..[k_] *
    //           right[j_0]..[j_][k_0]..[k_];
    //   }
    //
    //  Relay this summation to another helper class.

    template <int no_contr, int dim>
    class Contract<no_contr, no_contr, no_contr, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static void
      contract(T1 &result, const T2 &left, const T3 &right)
      {
        result = Contract2<no_contr, dim>::template contract2<T1>(left, right);
      }
    };

    // Straightforward recursion:
    //
    // Contract leftmost index and recurse one down.

    template <int no_contr, int dim>
    class Contract2
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr inline DEAL_II_ALWAYS_INLINE static T1
      contract2(const T2 &left, const T3 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        if (dim == 0)
          {
            T1 result = dealii::internal::NumberType<T1>::value(0.0);
            return result;
          }
        else
          {
            T1 result =
              Contract2<no_contr - 1, dim>::template contract2<T1>(left[0],
                                                                   right[0]);
            for (unsigned int i = 1; i < dim; ++i)
              result +=
                Contract2<no_contr - 1, dim>::template contract2<T1>(left[i],
                                                                     right[i]);
            return result;
          }
      }
    };

    // A contraction of two objects of order 0 is just a scalar
    // multiplication:

    template <int dim>
    class Contract2<0, dim>
    {
    public:
      template <typename T1, typename T2, typename T3>
      constexpr DEAL_II_ALWAYS_INLINE static T1
      contract2(const T2 &left, const T3 &right)
      {
        return left * right;
      }
    };


    // -------------------------------------------------------------------------
    // Implementation of helper classes for contract3
    // -------------------------------------------------------------------------

    // Fully contract three tensorial objects
    //
    // As long as rank_1 > 0, recurse over left and middle:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       [...]
    //         left[i_0]..[i_] ... middle[i_0]..[i_] ... right

    template <int rank_1, int rank_2, int dim>
    class Contract3
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static inline T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        T1 result = dealii::internal::NumberType<T1>::value(0.0);
        for (unsigned int i = 0; i < dim; ++i)
          result += Contract3<rank_1 - 1, rank_2, dim>::template contract3<T1>(
            left[i], middle[i], right);
        return result;
      }
    };

    // If rank_1 ==0, continue to recurse over middle and right:
    //
    // for(unsigned int i_0; i_0 < dim; ++i_0)
    //   ...
    //     for(i_; i_ < dim; ++i_)
    //       for(unsigned int j_0; j_0 < dim; ++j_0)
    //         ...
    //           for(j_; j_ < dim; ++j_)
    //             [...]
    //               left[i_0]..[i_] ... middle[i_0]..[i_][j_0]..[j_] ...
    //               right[j_0]..[j_]

    template <int rank_2, int dim>
    class Contract3<0, rank_2, dim>
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static inline T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        // Some auto-differentiable numbers need explicit
        // zero initialization.
        T1 result = dealii::internal::NumberType<T1>::value(0.0);
        for (unsigned int i = 0; i < dim; ++i)
          result +=
            Contract3<0, rank_2 - 1, dim>::template contract3<T1>(left,
                                                                  middle[i],
                                                                  right[i]);
        return result;
      }
    };

    // Contraction of three tensorial objects of rank 0 is just a scalar
    // multiplication.

    template <int dim>
    class Contract3<0, 0, dim>
    {
    public:
      template <typename T1, typename T2, typename T3, typename T4>
      constexpr static T1
      contract3(const T2 &left, const T3 &middle, const T4 &right)
      {
        return left * middle * right;
      }
    };

    // -------------------------------------------------------------------------

  } /* namespace internal */
} /* namespace TensorAccessors */

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii_tensor_accessors_h */
