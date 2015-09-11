// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2015 by the deal.II authors
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

#ifndef dealii__tensor_accessors_h
#define dealii__tensor_accessors_h

#include <deal.II/base/config.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/table_indices.h>


DEAL_II_NAMESPACE_OPEN

/**
 * This namespace is a collection of algorithms working on generic
 * tensorial objects (of arbitrary rank).
 *
 * The rationale to implement such functionality in a generic fashion in a
 * separate namespace is
 *  - to easy code reusability and therefore avoid code duplication.
 *  - to have a well-defined interface that allows to exchange the low
 *    level implementation.
 *
 *
 * A tensorial object has the notion of a rank and allows a rank-times
 * recursive application of the index operator, e.g., if <code>t</code> is
 * a tensorial object of rank 4, the following access is valid:
 * @code
 *   t[1][2][1][4]
 * @endcode
 *
 * deal.II has its own implementation for tensorial objects such as
 * dealii::Tensor<rank, dim, Number> and
 * dealii::SymmetricTensor<rank, dim, Number>
 *
 * The methods and algorithms implemented in this namespace, however, are
 * fully generic. More precisely, it can operate on nested c-style arrays,
 * or on class types <code>T</code> with a minimal interface that provides
 * a local typedef <code>value_type</code> and an index operator
 * <code>operator[](unsigned int)</code> that returns a (const or
 * non-const) reference of <code>value_type</code>:
 * @code
 *   template<...>
 *   class T
 *   {
 *     typedef ... value_type;
 *     value_type & operator[](unsigned int);
 *     const value_type & operator[](unsigned int) const;
 *   };
 * @endcode
 *
 * This namespace provides primitves for access, reordering and contraction
 * of such objects.
 *
 * @ingroup geomprimitives
 */
namespace TensorAccessors
{
  // forward declarations
  namespace internal
  {
    template <int index, int rank, typename T> class ReorderedIndexView;
    template <int position, int rank> struct ExtractHelper;
  }


  /**
   * This class provides a local typedef @p value_type denoting the
   * resulting type of an access with operator[](unsigned int). More
   * precisely, @p value_type will be
   *  - <code>T::value_type</code> if T is a tensorial class providing a
   *    typedef <code>value_type</code> and does not have a const qualifier.
   *  - <code>const T::value_type</code> if T is a tensorial class
   *    providing a typedef <code>value_type</code> and does have a const
   *    qualifier.
   *  - <code>const T::value_type</code> if T is a tensorial class
   *    providing a typedef <code>value_type</code> and does have a const
   *    qualifier.
   *  - <code>A</code> if T is of array type <code>A[...]</code>
   *  - <code>const A</code> if T is of array type <code>A[...]</code> and
   *    does have a const qualifier.
   */
  template <typename T>
  struct ValueType
  {
    typedef typename T::value_type value_type;
  };

  template <typename T>
  struct ValueType<const T>
  {
    typedef const typename T::value_type value_type;
  };

  template <typename T, std::size_t N>
  struct ValueType<T[N]>
  {
    typedef T value_type;
  };

  template <typename T, std::size_t N>
  struct ValueType<const T[N]>
  {
    typedef const T value_type;
  };


  /**
   * This class provides a local typedef @p value_type that is equal to
   * the typedef <code>value_type</code> after @p deref_steps
   * recursive dereferenciations via ```operator[](unsigned int)```.
   * Further, constness is preserved via the ValueType
   * type trait, i.e., if T is const, ReturnType<rank, T>::value_type
   * will also be const.
   */
  template <int deref_steps, typename T>
  struct ReturnType
  {
    typedef typename ReturnType<deref_steps - 1, typename ValueType<T>::value_type>::value_type value_type;
  };

  template <typename T>
  struct ReturnType<0, T>
  {
    typedef T value_type;
  };


  /**
   * Provide a "tensorial view" to a reference @p t of a tensor object of
   * rank @p rank in which the index @p index is shifted to the
   * end. As an example consider a tensor of 5th order in dim=5 space
   * dimensions that can be accessed through 5 recursive
   * <code>operator[]()</code> invocations:
   * @code
   *   Tensor<5, dim> tensor;
   *   tensor[0][1][2][3][4] = 42.;
   * @endcode
   * Index 1 (the 2nd index, count starts at 0) can now be shifted to the
   * end via
   * @code
   *   auto tensor_view = reordered_index_view<1, 5>(tensor);
   *   tensor_view[0][2][3][4][1] == 42.; // is true
   * @endcode
   * The usage of the dealii::Tensor type was solely for the sake of an
   * example. The mechanism implemented by this function is available for
   * fairly general tensorial types @p T.
   *
   * The purpose of this reordering facility is to be able to contract over
   * an arbitrary index of two (ore more) tensors:
   *  - reorder the indices in mind to the end of the tensors
   *  - use belows contract function that contracts the _last_ elements of
   *    tensors.
   *
   * @note This function returns an internal class object consisting of an
   * array subscript operator <code>operator[](unsigned int)</code> and a
   * typedef <code>value_type</code> describing its return value.
   *
   * @tparam index The index to be shifted to the end. Indices are counted
   * from 0, thus the valid range is $0\le\text{index}<\text{rank}$.
   * @tparam rank Rank of the tensorial object @param t
   * @tparam T A tensorial object of rank @p rank. @p T must
   * provide a local typedef <code>value_type</code> and an index operator
   * <code>operator[]()</code> that returns a (const or non-const)
   * reference of <code>value_type</code>:
   * @code
   *   class T
   *   {
   *     typedef ... value_type
   *     value_type & operator[](unsigned int);
   *     const value_type & operator[](unsigned int) const;
   *   };
   * @endcode
   *
   * @relates ReorderedIndexView
   */
  template <int index, int rank, typename T>
  internal::ReorderedIndexView<index, rank, T>
  reordered_index_view(T &t)
  {
#ifdef DEAL_II_WITH_CXX11
    static_assert(0 <= index && index < rank,
                  "The specified index must lie within the range [0,rank)");
#endif

    return internal::ReorderedIndexView<index, rank, T>(t);
  }


  /**
   * Return a reference (const or non-const) to a subobject of a tensorial
   * object @p t of type @p T, as described by an array type @p ArrayType
   * object @p indices. For example: @code
   *   Tensor<5, dim> tensor;
   *   TableIndices<5> indices (0, 1, 2, 3, 4);
   *   TensorAccessors::extract(tensor, indices) = 42;
   * @endcode
   * This is equivalent to <code>tensor[0][1][2][3][4] = 42.</code>.
   *
   * @tparam T A tensorial object of rank @p rank. @p T must provide a
   * local typedef <code>value_type</code> and an index operator
   * <code>operator[]()</code> that returns a (const or non-const)
   * reference of <code>value_type</code>. Further, its tensorial rank must
   * be equal or greater than @p rank
   *
   * @tparam ArrayType An array like object, such as std::array, or
   * dealii::TableIndices  that stores at least @p rank indices that can be
   * accessed via operator[]().
   */
  template<int rank, typename T, typename ArrayType> typename
  ReturnType<rank, T>::value_type &
  extract(T &t, const ArrayType &indices)
  {
    return internal::ExtractHelper<0, rank>::template extract<T, ArrayType>(t, indices);
  }


  namespace internal
  {
    // -------------------------------------------------------------------------
    // Forward declarations and type traits
    // -------------------------------------------------------------------------

    template <int rank, typename S> class StoreIndex;
    template <typename T> class Identity;
    template <int no_contr, int dim> class Contract2;

    /**
     * An internally used type trait to allow nested application of the
     * function reordered_index_view(T &t).
     *
     * The problem is that when working with the actual tensorial types, we
     * have to return subtensors by reference - but sometimes, especially
     * for StoreIndex and ReorderedIndexView that return rvalues, we have
     * to return by value.
     */
    template<typename T>
    struct ReferenceType
    {
      typedef T &type;
    };

    template <int rank, typename S>
    struct ReferenceType<StoreIndex<rank, S> >
    {
      typedef StoreIndex<rank, S> type;
    };

    template <int index, int rank, typename T>
    struct ReferenceType<ReorderedIndexView<index, rank, T> >
    {
      typedef ReorderedIndexView<index, rank, T> type;
    };


    /**
     * An internally used type trait that strips StoreIndex<0, S> down to
     * its actual return value. This is needed to end the recursion in
     * StoreIndex as well as, to return something meaningful in case of a
     * nested application of the index reordering classes.
     */
    template<typename T>
    struct StripStoreIndex
    {
      typedef T type;
    };

    template <typename S>
    struct StripStoreIndex<StoreIndex<0, S> >
    {
      typedef typename StripStoreIndex<typename S::return_type>::type type;
    };


    // TODO: Is there a possibility ot just have the following block of
    // explanation on an internal page in doxygen? If, yes. Doxygen
    // wizards, your call!

    // -------------------------------------------------------------------------
    // Implemenation of helper classes for reordered_index_view
    // -------------------------------------------------------------------------

    // OK. This is utterly brutal template magic. Therefore, we will not
    // comment on the individual internal helper classes, because this is
    // of not much value, but explain the general recursion procedure.
    //
    // (In order of appearance)
    //
    // Our task is to reorder access to a tensor object where a specified
    // index is moved to the end. Thus we want to construct an object
    // <code>reorderd</code> out of a <code>tensor</code> where the
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
      ReorderedIndexView(typename ReferenceType<T>::type t) : t_(t) {}

      typedef ReorderedIndexView<index - 1, rank - 1, typename ValueType<T>::value_type>
      value_type;

      // Recurse by applying index j directly:
      inline
      value_type operator[](unsigned int j) const
      {
        return value_type(t_[j]);
      }

    private:
      typename ReferenceType<T>::type t_;
    };

    // At some point we hit the condition index == 0, i.e., the first index
    // should be reordered to the end.
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
      ReorderedIndexView(typename ReferenceType<T>::type t) : t_(t) {}

      typedef internal::StoreIndex<rank - 1, internal::Identity<T> > value_type;

      inline
      value_type operator[](unsigned int j) const
      {
        return value_type(internal::Identity<T>(t_), j);
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
      Identity(typename ReferenceType<T>::type t) : t_(t) {}

      typedef typename ValueType<T>::value_type return_type;

      inline
      typename ReferenceType<return_type>::type apply(unsigned int j) const
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
    class StoreIndex : private S
    {
    public:
      StoreIndex(S s, int i) : S(s), i_(i) {}

      typedef StoreIndex<rank - 1, StoreIndex<rank, S> > value_type;

      inline
      value_type operator[](unsigned int j) const
      {
        return value_type(*this, j);
      }

      typedef typename ValueType<typename S::return_type>::value_type return_type;

      inline
      typename ReferenceType<return_type>::type apply(unsigned int j) const
      {
        return S::apply(j)[i_];
      }

    private:
      const int i_;
    };

    // We can store indices until we hit rank == 0. Then, we have all
    // necessary indices and it is time to ground the recursion. For this,
    // StoreIndex is specialized for rank == 0. The specialization contains
    // to conversion operators to reference type to access the underlying
    // object. Just call apply(i_) on the StoreIndex object of rank 1:
    template <typename S>
    class StoreIndex<0, S> : private S
    {
    public:
      StoreIndex(S s, int i) : S(s), i_(i) {}

      // Strip nested StoreIndex<0, S2> objects and cast to tensor's value_type
      typedef typename StripStoreIndex<typename S::return_type>::type value_type;

      inline operator value_type &() const
      {
        return S::apply(i_);
      }

    private:
      const int i_;
    };


    // -------------------------------------------------------------------------
    // Implemenation of helper classes for extract
    // -------------------------------------------------------------------------

    // Straightforward recursion implemented by specializing ExtractHelper
    // for position == rank. We use the type trait ReturnType<rank, T> to
    // have an idea what the final type will be.
    template<int position, int rank>
    struct ExtractHelper
    {
      template<typename T, typename ArrayType>
      inline static
      typename ReturnType<rank - position, T>::value_type &
      extract(T &t, const ArrayType &indices)
      {
        return ExtractHelper<position + 1, rank>::
               template extract<typename ValueType<T>::value_type, ArrayType>
        (t[indices[position]], indices);
      }
    };

    // For dimension == rank there is nothing to extract, just return the
    // object.
    template<int rank>
    struct ExtractHelper<rank, rank>
    {
      template<typename T, typename ArrayType>
      inline static
      T &extract(T &t, const ArrayType &indices)
      {
        return t;
      }
    };

    // -------------------------------------------------------------------------

  } /* namespace internal */
} /* namespace TensorAccessors */

DEAL_II_NAMESPACE_CLOSE

#endif /* dealii__tensor_accessors_h */
