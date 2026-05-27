// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2026 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#ifndef dealii_cxx26_inplace_vector_h
#define dealii_cxx26_inplace_vector_h

#include <deal.II/base/config.h>


#ifdef DEAL_II_HAVE_CXX26
#  include <inplace_vector>
#else
#  include <deal.II/base/exception_macros.h>
#  include <deal.II/base/exceptions.h>

// Header optimization: we can skip including utility and iterator and instead
// rely on the implicit includes (which are smaller, e.g., <array> is less than
// half the preprocessed size of <iterator>) from <array>
#  include <algorithm>
#  include <array>
#  include <cstdint>
#  include <initializer_list>
#  include <limits>
#  include <type_traits>
#endif

#include <boost/serialization/array_wrapper.hpp>
#include <boost/serialization/collection_size_type.hpp>
#include <boost/serialization/split_free.hpp>

DEAL_II_NAMESPACE_OPEN

namespace std_cxx26
{
#ifndef DEAL_II_HAVE_CXX26
  DeclExceptionMsg(ExcCapacityExceeded,
                   "The current operation requires more capacity than the "
                   "container can provide.");

  /**
   * C++17-implementation of a subset of
   * [std::inplace_vector](https://en.cppreference.com/cpp/container/inplace_vector).
   *
   * Differences:
   *
   * 1. The number of elements is limited to the maximum value of an
   *    `unsigned short`, which is typically 65535.
   * 2. This class lacks optimized storage for the N = 0 case.
   * 3. Since std::lexicographical_compare() and std::equal() are not
   *    constexpr prior to C++20,  the comparison operators are not constexpr.
   * 4. Since operator<=>() is not available prior to C++20, this class
   *    implements the full set of comparison operators.
   * 5. The range constructor is not present sinces ranges were not available
   *    prior to C++20.
   * 6. Since std::optional<T&> was not available prior to C++26, this class
   *    does not support try_emplace_back() etc.
   * 7. Since placement new() is not constexpr prior to C++26, the majority of
   *    the constructors and assignment operators are not constexpr.
   *
   * Furthermore, copy and move constructors, copy and move assignment
   * operators, and the destructor are only trivial if both T is trivial for
   * that function and C++20 is available, since declaring this property
   * for a user-defined class requires concepts.
   *
   * See https://www.open-std.org/jtc1/sc22/wg21/docs/papers/2024/p2747r2.html
   * for more information.
   */
  template <typename T, std::size_t N>
  class inplace_vector
  {
  public:
    /**
     * Types.
     */
    /** @{ */
    using value_type      = T;
    using size_type       = std::size_t;
    using difference_type = std::ptrdiff_t;
    using reference       = value_type &;
    using const_reference = const value_type &;
    using pointer         = value_type *;
    using const_pointer   = const value_type *;
    using iterator        = pointer;
    using const_iterator  = const_pointer;

    using reverse_iterator       = std::reverse_iterator<iterator>;
    using const_reverse_iterator = std::reverse_iterator<const_iterator>;
    /** @} */

    /**
     * Constructors and destructor.
     */
    /** @{ */

    /**
     * Default constructor.
     */
    constexpr inplace_vector() noexcept
      : n_elements(0)
    {}

    /**
     * Constructor which default initializes @p n objects.
     */
    explicit inplace_vector(size_type n)
      : n_elements(0)
    {
      internal_resize<true>(n);
    }

    /**
     * Constructor which creates @p n objects by copying @p value.
     */
    inplace_vector(size_type n, const T &value)
      : n_elements(0)
    {
      if (n > N)
        throw std::bad_alloc();
      internal_append(SingleInputIterator(std::addressof(value), 0),
                      SingleInputIterator(std::addressof(value), n));
    }

    /**
     * Constructor which copies an iterator range.
     */
    template <class InputIterator>
    inplace_vector(InputIterator first, InputIterator last)
      : n_elements(0)
    {
      // If we have two integers (e.g., inplace_vector(42, 42)) then we should
      // convert to the correct types and call the other constructor
      if constexpr (std::is_integral_v<InputIterator>)
        *this = inplace_vector(static_cast<size_type>(first),
                               static_cast<const T &>(last));
      else
        internal_append<InputIterator, true, false>(first, last);
    }

#  ifdef DEAL_II_HAVE_CXX20
    /**
     * Trivial copy constructor. Requires C++20.
     */
    inplace_vector(const inplace_vector &)
      requires(std::is_trivially_copy_constructible_v<T>)
    = default;
#  endif

    /**
     * Copy constructor.
     */
    inplace_vector(const inplace_vector &other) noexcept(
      std::is_nothrow_copy_constructible_v<T>)
      : n_elements(0)
    {
      internal_append(other.begin(), other.end());
    }

#  ifdef DEAL_II_HAVE_CXX20
    /**
     * Trivial move constructor. Requires C++20.
     */
    inplace_vector(inplace_vector &&)
      requires(std::is_trivially_move_constructible_v<T>)
    = default;
#  endif

    /**
     * Move constructor.
     */
    inplace_vector(inplace_vector &&other) noexcept(
      std::is_nothrow_move_constructible_v<T>)
      : n_elements(0)
    {
      internal_append<iterator, false, true>(other.begin(), other.end());
    }

    /**
     * Copy constructor from an initializer list.
     */
    constexpr inplace_vector(std::initializer_list<T> init)
      : n_elements(0)
    {
      if (init.size() > N)
        throw std::bad_alloc();
      internal_append(init.begin(), init.end());
    }

#  ifdef DEAL_II_HAVE_CXX20
    /**
     * Trivial destructor. Requires C++20.
     */
    ~inplace_vector()
      requires(std::is_trivially_destructible_v<T>)
    = default;
#  endif

#  ifdef DEAL_II_HAVE_CXX20
    /**
     * Destructor.
     */
    constexpr
#  endif
      ~inplace_vector()
    {
      clear();
    }
    /** @} */

    /**
     * Assignment operators.
     */
    /** @{ */

    /**
     * Trivial copy assignment operator. Requires C++20.
     */
#  ifdef DEAL_II_HAVE_CXX20
    inplace_vector &
    operator=(const inplace_vector &)
      requires((std::is_trivially_copy_assignable_v<T> &&
                std::is_trivially_copy_constructible_v<T> &&
                std::is_trivially_destructible_v<T>))
    = default;
#  endif

    /**
     * Copy assignment operator.
     */
    inplace_vector &
    operator=(const inplace_vector &other)
    {
      if (std::addressof(other) != this)
        internal_assign(other.begin(), other.end());

      return *this;
    }

#  ifdef DEAL_II_HAVE_CXX20
    /**
     * Trivial move assignment operator. Requires C++20.
     */
    inplace_vector &
    operator=(inplace_vector &&)
      requires(std::is_trivially_move_assignable_v<T> &&
               std::is_trivially_move_constructible_v<T> &&
               std::is_trivially_destructible_v<T>)
    = default;
#  endif

    /**
     * Move assignment operator.
     */
    inplace_vector &
    operator=(inplace_vector &&other) noexcept(
      std::is_nothrow_move_assignable_v<T>
        &&std::is_nothrow_move_constructible_v<T>)
    {
      if (std::addressof(other) != this)
        internal_assign<iterator, false, true>(other.begin(), other.end());

      return *this;
    }

    /**
     * Copy assignment operator (from an initializer list).
     */
    inplace_vector &
    operator=(std::initializer_list<T> init)
    {
      if (init.size() > N)
        throw std::bad_alloc();

      internal_assign(init.begin(), init.end());

      return *this;
    }
    /** @} */

    /**
     * Assignment functions.
     */
    /** @{ */

    /**
     * Assign the contents of the vector to values in the input iterator range.
     */
    template <class InputIterator>
    void
    assign(InputIterator first, InputIterator last)
    {
      if constexpr (std::is_integral_v<InputIterator>)
        assign(static_cast<size_type>(first), static_cast<const T &>(last));
      else
        internal_assign(first, last);
    }

    /**
     * Assign the contents of the vector to @p n copies of @p value.
     */
    void
    assign(size_type n, const T &value)
    {
      if (n > N)
        throw std::bad_alloc();

      internal_assign(SingleInputIterator(std::addressof(value), 0),
                      SingleInputIterator(std::addressof(value), n));
    }

    /**
     * Assign the contents of the vector to values in the initializer list.
     */
    void
    assign(std::initializer_list<T> init)
    {
      if (init.size() > N)
        throw std::bad_alloc();

      internal_assign(init.begin(), init.end());
    }
    /** @} */

    /**
     * Iterators.
     */
    /** @{ */
    /**
     * Begin iterator.
     */
    constexpr iterator
    begin() noexcept
    {
      return std::addressof(elements[0]);
    }

    /**
     * End iterator.
     */
    constexpr iterator
    end() noexcept
    {
      return std::addressof(elements[0]) + size();
    }

    /**
     * Constant begin iterator.
     */
    constexpr const_iterator
    begin() const noexcept
    {
      return std::addressof(elements[0]);
    }

    /**
     * Constant end iterator.
     */
    constexpr const_iterator
    end() const noexcept
    {
      return std::addressof(elements[0]) + size();
    }

    /**
     * Reverse begin iterator.
     */
    constexpr reverse_iterator
    rbegin() noexcept
    {
      return reverse_iterator(end());
    }

    /**
     * Reverse end iterator.
     */
    constexpr reverse_iterator
    rend() noexcept
    {
      return reverse_iterator(begin());
    }

    /**
     * Constant reverse begin iterator.
     */
    constexpr const_reverse_iterator
    rbegin() const noexcept
    {
      return const_reverse_iterator(end());
    }

    /**
     * Constant reverse end iterator.
     */
    constexpr const_reverse_iterator
    rend() const noexcept
    {
      return const_reverse_iterator(begin());
    }

    /**
     * Constant begin iterator.
     */
    constexpr const_iterator
    cbegin() const noexcept
    {
      return std::addressof(elements[0]);
    }

    /**
     * Constant end iterator.
     */
    constexpr const_iterator
    cend() const noexcept
    {
      return std::addressof(elements[0]) + size();
    }

    /**
     * Constant reverse begin iterator.
     */
    constexpr const_reverse_iterator
    crbegin() const noexcept
    {
      return const_reverse_iterator(cend());
    }

    /**
     * Constant reverse end iterator.
     */
    constexpr const_reverse_iterator
    crend() const noexcept
    {
      return const_reverse_iterator(cbegin());
    }
    /** @} */

    /**
     * Capacity.
     */
    /** @{ */
    /**
     * Return whether or not the vector contains any objects.
     */
    constexpr bool
    empty() const noexcept
    {
      return size() == 0;
    }

    /**
     * Return the number of objects in the vector.
     */
    constexpr size_type
    size() const noexcept
    {
      Assert(n_elements <= N, ExcInternalError());
      return n_elements;
    }

    /**
     * Maximum possible size of the vector.
     */
    static constexpr size_type
    max_size() noexcept
    {
      return N;
    }

    /**
     * Maximum possible capacity of the vector.
     */
    static constexpr size_type
    capacity() noexcept
    {
      return N;
    }

    /**
     * Resize the vector to have @p n objects.
     */
    void
    resize(size_type n)
    {
      internal_resize<true>(n);
    }

    /**
     * Resize the vector to have @p n objects by constructing new objects with
     * copies of @p value.
     */
    void
    resize(size_type n, const T &value)
    {
      internal_resize<true>(n, value);
    }

    /**
     * If $n > N$ then throw an exception. This function exists for
     * compatibility with other containers but, since this vector has a fixed
     * internal storage size, calling this function does not affect the object
     * in any way.
     */
    static constexpr void
    reserve(size_type n)
    {
      if (n > N)
        throw std::bad_alloc();
    }

    /**
     * This function does nothing for inplace_vector and exists for
     * compatibility with other containers.
     */
    static constexpr void
    shrink_to_fit() noexcept
    {}
    /** @} */

    /**
     * Element access.
     */
    /** @{ */

    /**
     * Return a reference to the @p n-th element.
     */
    reference
    operator[](size_type n)
    {
      AssertIndexRange(n, size());
      return elements[n];
    }

    /**
     * Return a constant reference to the @p n-th element.
     */
    const_reference
    operator[](size_type n) const
    {
      AssertIndexRange(n, size());
      return elements[n];
    }

    /**
     * If $n$ is less than size() then return a reference to the @p n-th
     * element: otherwise throw an exception.
     */
    reference
    at(size_type n)
    {
      if (!(n < size()))
        throw std::out_of_range("inplace_vector::at(): index " +
                                std::to_string(n) + " out of range");
      return elements[n];
    }

    /**
     * If $n$ is less than size() then return a constant reference to the @p
     * n-th element: otherwise throw an exception.
     */
    const_reference
    at(size_type n) const
    {
      if (!(n < size()))
        throw std::out_of_range("inplace_vector::at(): index " +
                                std::to_string(n) + " out of range");
      return elements[n];
    }

    /**
     * Return a reference to the first object in the vector.
     */
    reference
    front()
    {
      Assert(!empty(), ExcEmptyObject());
      return elements[0];
    }

    /**
     * Return a constant reference to the first object in the vector.
     */
    const_reference
    front() const
    {
      Assert(!empty(), ExcEmptyObject());
      return elements[0];
    }

    /**
     * Return a reference to the last object in the vector.
     */
    reference
    back()
    {
      Assert(!empty(), ExcEmptyObject());
      return elements[size() - 1];
    }

    /**
     * Return a constant reference to the last object in the vector.
     */
    const_reference
    back() const
    {
      Assert(!empty(), ExcEmptyObject());
      return elements[size() - 1];
    }
    /** @} */

    /**
     * Data access.
     */
    /** @{ */

    /**
     * Return a pointer to all objects in the vector.
     *
     * @note This does not return a nullptr to signify no objects.
     */
    T *
    data() noexcept
    {
      return begin();
    }

    /**
     * Return a constant pointer to all objects in the vector.
     *
     * @note This does not return a nullptr to signify no objects.
     */
    const T *
    data() const noexcept
    {
      return cbegin();
    }
    /** @} */

    /**
     * Modifiers.
     */
    /** @{ */

    /**
     * Construct a new object at the end of the vector by forwarding all
     * arguments to its constructor.
     *
     * @note If the vector is full before this function is called then an
     * exception is thrown.
     */
    template <class... Args>
    reference
    emplace_back(Args &&...args)
    {
      internal_resize<true>(size() + 1, std::forward<Args>(args)...);
      return back();
    }

    /**
     * Construct a new object at the end of the vector by copying @p value.
     *
     * @note If the vector is full before this function is called then an
     * exception is thrown.
     */
    reference
    push_back(const T &value)
    {
      internal_resize<true>(size() + 1, value);
      return back();
    }

    /**
     * Construct a new object at the end of the vector by moving @p value.
     *
     * @note If the vector is full before this function is called then an
     * exception is thrown.
     */
    reference
    push_back(T &&value)
    {
      internal_resize<true>(size() + 1, std::forward<T>(value));
      return back();
    }

    /**
     * Destruct the last object in the vector and reduce the size of this object
     * by one.
     */
    void
    pop_back()
    {
      Assert(!empty(), ExcEmptyObject());
      internal_resize<true>(size() - 1);
    }

    /**
     * Construct a new object at the end of the vector by forwarding all objects
     * to its constructor.
     *
     * @note Unlike the checked versions, this function will not throw an out of
     * memory exception if the vector's capacity is exceeded.
     */
    template <class... Args>
    reference
    unchecked_emplace_back(Args &&...args)
    {
      Assert(size() < capacity(), ExcCapacityExceeded());
      internal_resize<false>(size() + 1, std::forward<Args>(args)...);
      return back();
    }

    /**
     * Construct a new object at the end of the vector by copying @p value.
     *
     * @note Unlike the checked versions, this function will not throw an out of
     * memory exception if the vector's capacity is exceeded.
     */
    reference
    unchecked_push_back(const T &value)
    {
      Assert(size() < capacity(), ExcCapacityExceeded());
      internal_resize<false>(size() + 1, value);
      return back();
    }

    /**
     * Construct a new object at the end of the vector by moving @p value.
     *
     * @note Unlike the checked versions, this function will not throw an out of
     * memory exception if the vector's capacity is exceeded.
     */
    reference
    unchecked_push_back(T &&value)
    {
      Assert(size() < capacity(), ExcCapacityExceeded());
      internal_resize<false>(size() + 1, std::forward<T>(value));
      return back();
    }

    /**
     * Construct a new object at the location of @p position by forwarding all
     * objects (@p args etc.) to its constructor.
     */
    template <class... Args>
    iterator
    emplace(const_iterator position, Args &&...args)
    {
      const auto index = position - cbegin();
      AssertIndexRange(index, size());
      // Since args may reference *this we have to construct the new object
      // first: do that in the buffer and then rotate so it is in the correct
      // place
      internal_resize<true>(size() + 1, std::forward<Args>(args)...);

      Assert(begin() < end(), ExcInternalError());
      Assert(begin() + index <= end() - 1, ExcInternalError());
      std::rotate(begin() + index, end() - 1, end());

      return begin() + index;
    }

    /**
     * Construct a new object at the location of @p position by copying @p value.
     */
    iterator
    insert(const_iterator position, const T &value)
    {
      return emplace(position, value);
    }

    /**
     * Construct a new object at the location of @p position by moving @p value.
     */
    iterator
    insert(const_iterator position, T &&value)
    {
      return emplace(position, std::move(value));
    }

    /**
     * Construct $n$ new objects at the location of @p position by copying @p value.
     */
    iterator
    insert(const_iterator position, size_type n, const T &value)
    {
      return insert(position,
                    SingleInputIterator(std::addressof(value), 0),
                    SingleInputIterator(std::addressof(value), n));
    }

    /**
     * Construct new objects at the location of @p position from the iterator
     * range defined by @p first and @p last.
     */
    template <typename InputIterator>
    iterator
    insert(const_iterator position, InputIterator first, InputIterator last)
    {
      if constexpr (std::is_integral_v<InputIterator>)
        {
          return insert(position,
                        static_cast<size_type>(first),
                        static_cast<const T &>(last));
        }
      else
        {
          const auto original_size = size();
          // TODO we need a nicer way to check for valid iterators
          const auto index = position - cbegin();
          Assert(position == cend() || index < size(),
                 ExcMessage("out of range"));
          const auto n_new_elements =
            internal_append<InputIterator>(first, last);
          Assert(position + n_new_elements <= end(), ExcInternalError());
          std::rotate(begin() + index, begin() + original_size, end());
          return begin() + index;
        }
    }

    /**
     * Insert values from an initializer list at a position.
     */
    iterator
    insert(const_iterator position, std::initializer_list<T> init)
    {
      return insert(position, init.begin(), init.end());
    }

    /**
     * Erase a value.
     */
    iterator
    erase(const_iterator position)
    {
      const auto index = position - cbegin();
      AssertIndexRange(index, size());
      Assert(begin() + index + 1 <= end(), ExcInternalError());
      std::move(begin() + index + 1, end(), begin() + index);
      pop_back();

      return begin() + index;
    }

    /**
     * Erase a range of values.
     */
    iterator
    erase(const_iterator first, const_iterator last)
    {
      const auto original_size = size();
      const auto distance      = std::distance(first, last);
      const auto first_index   = first - begin();
      const auto last_index    = last - cbegin();
      Assert(first_index <= last_index,
             ExcMessage("The given range is not valid."));
      std::rotate(begin() + first_index, begin() + last_index, end());
      internal_resize<true>(size() - distance);
      Assert(size() + distance == original_size, ExcInternalError());

      return begin() + first_index;
    }

    /**
     * Swap the contents of this vector with another.
     */
    void
    swap(inplace_vector &other) noexcept(
      std::is_nothrow_swappable_v<T> &&std::is_nothrow_move_constructible_v<T>)
    {
      if (std::addressof(other) == this)
        return;

      auto      &smaller      = size() < other.size() ? *this : other;
      const auto smaller_size = smaller.size();
      auto      &larger       = size() < other.size() ? other : *this;

      // Do the same thing as GCC: copy or move elements from the larger to the
      // smaller and then swap elements with the same indices.
      for (std::size_t i = smaller.size(); i < larger.size(); ++i)
        {
          if constexpr (std::is_nothrow_move_constructible_v<T>)
            smaller.push_back(std::move(larger[i]));
          else
            smaller.push_back(larger[i]);
        }
      larger.resize(smaller_size);
      using std::swap;
      for (std::size_t i = 0; i < smaller_size; ++i)
        swap(smaller[i], larger[i]);
    }

    /**
     * Friend swap function.
     */
    friend void
    swap(inplace_vector &x, inplace_vector &y) noexcept(
      std::is_nothrow_swappable_v<T> &&std::is_nothrow_move_constructible_v<T>)
    {
      x.swap(y);
    }

    /**
     * Clear the contents of the vector.
     */
    constexpr void
    clear() noexcept
    {
      internal_resize(0);
    }

    /** @} */
  private:
    /**
     * Wrapper class for an iterator-like interface to a single value.
     *
     * Most constructors and assign() overloads work with ranges of iterators.
     * SingleInputIterator allows for calling those functions (or common
     * utilities) from the functions which takes counts and values as well.
     *
     * @warning This class stores a pointer to the exemplar and it should be
     * used with care, e.g,. do not create pointers to objects currently inside
     * the inplace_vector.
     */
    class SingleInputIterator
    {
    public:
      using value_type        = const T;
      using reference         = const T &;
      using iterator_category = std::input_iterator_tag;

      SingleInputIterator(const T *value, size_type index)
        : ptr(value)
        , index(index)
      {}

      reference
      operator*() const
      {
        return *ptr;
      }

      bool
      operator==(const SingleInputIterator &other) const
      {
        return ptr == other.ptr && index == other.index;
      }

      bool
      operator!=(const SingleInputIterator &other) const
      {
        return !(*this == other);
      }

      SingleInputIterator &
      operator++()
      {
        ++index;
        return *this;
      }

    private:
      const T *ptr;

      std::size_t index;
    };

    /**
     * Common function for resizing the vector.
     */
    template <bool check = false, typename... Args>
    constexpr void
    internal_resize(const size_type n, Args &&...args) noexcept(!check)
    {
      Assert(n <= capacity(), ExcCapacityExceeded());
      static_assert(std::is_constructible_v<T, Args...>);
      if constexpr (check)
        if (n > N)
          throw std::bad_alloc();

      // Here and elsewhere we avoid std::uninitialized_copy() etc. so that we
      // maintain the invariant that n_elements is always correct
      if (n < size())
        for (size_type i = size(); i > n; --i)
          {
            Assert(size() > 0, ExcInternalError());
            back().~T();
            --n_elements;
          }
      else
        for (size_type i = size(); i < n; ++i)
          {
            AssertIndexRange(size(), capacity());
            new (end()) T(std::forward<Args>(args)...);
            ++n_elements;
          }
    }

    /**
     * Common function for assigning values to the vector from a range of input
     * iterators, first by copying and then by placement new.
     */
    template <typename InputIterator, bool check = false, bool move = false>
    void
    internal_assign(InputIterator first, InputIterator last) noexcept(
      !check && (move ? std::is_nothrow_move_assignable_v<T> :
                        std::is_nothrow_copy_assignable_v<T>))
    {
      static_assert(std::is_convertible_v<decltype(*first), T>);
      size_type i = 0;
      while (i < size() && first != last)
        {
          if constexpr (move)
            elements[i] = std::move(*first);
          else
            elements[i] = *first;
          ++i;
          ++first;
        }
      if (first == last)
        resize(i);
      else
        internal_append<InputIterator, check, move>(first, last);
    }

    /**
     * Common function for creating or moving new values at the end of the
     * vector.
     */
    template <typename InputIterator, bool check = false, bool move = false>
    constexpr size_type
    internal_append(InputIterator first, InputIterator last) noexcept(
      !check && (move ? std::is_nothrow_move_constructible_v<T> :
                        std::is_nothrow_copy_constructible_v<T>))
    {
      static_assert(std::is_convertible_v<decltype(*first), T>);

      // TODO: for constexpr-ification with C++17 we need to find a way to copy
      // trivial values over without using new() (since that is not constexpr
      // until C++26)

      size_type count = 0;
      while (first != last)
        {
          AssertIndexRange(size(), capacity());
          if constexpr (check)
            if (size() == N)
              throw std::bad_alloc();
          if constexpr (move)
            new (end()) T(std::move(*first));
          else
            new (end()) T(*first);
          ++n_elements;
          ++first;
          ++count;
        }
      return count;
    }

    /**
     * Prevent initialization of unused elements by placing the array in a
     * union.
     */
    union
    {
      T elements[N > 0 ? N : 1];
    };

    using buffer_size_type =
      std::conditional_t<N <= std::numeric_limits<std::uint8_t>::max(),
                         std::uint8_t,
                         std::uint16_t>;
    static_assert(
      N <= std::numeric_limits<buffer_size_type>::max(),
      "This class only supports objects of size <= the maximum size of a "
      "std::uint16_t (65535).");

    /**
     * Present number of elements in the buffer.
     */
    buffer_size_type n_elements;
  };

  /**
   * Comparison operators.
   */
  /** @{ */
  /**
   * Equality operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator==(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return (a.size() == b.size()) && std::equal(a.begin(), a.end(), b.begin());
  }

  /**
   * Inequality operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator!=(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return !(a == b);
  }

  /**
   * Less than operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator<(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
  }

  /**
   * Greater than operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator>(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return b < a;
  }

  /**
   * Less than or equal to operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator<=(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return !(b < a);
  }

  /**
   * Greater than or equal to operator for inplace_vector.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N>
  bool
  operator>=(const inplace_vector<T, N> &a, const inplace_vector<T, N> &b)
  {
    return !(a < b);
  }
  /** @} */

  /**
   * Erase all values equal to @p value in @p vec.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N, typename U = T>
  std::size_t
  erase(std_cxx26::inplace_vector<T, N> &vec, const U &value)
  {
    auto       it    = std::remove(vec.begin(), vec.end(), value);
    const auto count = std::distance(it, vec.end());
    Assert(vec.size() >= count, ExcInternalError());
    vec.resize(vec.size() - count);
    return count;
  }

  /**
   * Erase all values which satisfy the predicate @p pred in @p vec.
   *
   * @relates inplace_vector
   */
  template <typename T, std::size_t N, typename Predicate>
  std::size_t
  erase_if(std_cxx26::inplace_vector<T, N> &vec, Predicate pred)
  {
    auto       it    = std::remove_if(vec.begin(), vec.end(), pred);
    const auto count = std::distance(it, vec.end());
    Assert(vec.size() >= count, ExcInternalError());
    vec.resize(vec.size() - count);
    return count;
  }
#else
  using std::erase;
  using std::erase_if;
  using std::inplace_vector;
#endif
} // namespace std_cxx26

DEAL_II_NAMESPACE_CLOSE

namespace boost
{
  namespace serialization
  {
    /**
     * Write the data of this object to a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive, typename T, std::size_t N>
    inline void
    serialize(Archive                                 &ar,
              dealii::std_cxx26::inplace_vector<T, N> &t,
              const unsigned int                       file_version)
    {
      boost::serialization::split_free(ar, t, file_version);
    }

    /**
     * Write the data of this object to a stream for the purpose of
     * serialization using the [BOOST serialization
     * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
     */
    template <class Archive, typename T, std::size_t N>
    inline void
    save(Archive                                       &ar,
         const dealii::std_cxx26::inplace_vector<T, N> &vec,
         const unsigned int /*version*/)
    {
      const boost::serialization::collection_size_type count(vec.size());
      ar << count;
      if (!vec.empty())
        ar << boost::serialization::make_array(vec.data(), count);
    }

    template <class Archive, typename T, std::size_t N>
    inline void
    load(Archive                                 &ar,
         dealii::std_cxx26::inplace_vector<T, N> &vec,
         const unsigned int /*version*/)
    {
      boost::serialization::collection_size_type count(vec.size());
      ar >> count;
      vec.resize(count);
      if (!vec.empty())
        ar >> boost::serialization::make_array(vec.data(), count);
    }
  } // namespace serialization
} // namespace boost

#endif
