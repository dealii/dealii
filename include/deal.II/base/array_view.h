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

#ifndef dealii_array_view_h
#define dealii_array_view_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/memory_space.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/tensor.h>

#include <boost/container/small_vector.hpp>

#include <array>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declaration
template <class T>
class AlignedVector;

template <int N, typename T>
class Table;

template <typename number>
class LAPACKFullMatrix;


/**
 * A class that represents a window of memory locations of type @p ElementType
 * and presents it as if it was an array. In essence, this class is nothing more
 * than just a pointer to the first location and an integer that represents the
 * length of the array in elements. The memory remains owned by whoever
 * allocated it, as this class does not take over ownership.
 *
 * The advantage of using this class is that you don't have to pass around
 * pairs of pointers and that <code>operator[]</code> checks for the validity
 * of the index with which you subscript this array view. Note that accessing
 * elements is only allowed if the underlying data is stored in CPU memory.
 *
 * This class can handle views to both non-constant and constant memory
 * locations. If you want to represent a view of a constant array, then the
 * template argument type of this class needs to be @p const as well. The
 * following code snippet gives an example:
 * @code
 * std::vector<int> array = get_data(); // a writable array
 * ArrayView<int> view (&array[5], 5); // a view of elements 5..9 (inclusive)
 * view[2] = 42; // array[7] is set to 42
 * ArrayView<const int> const_view (&array[5], 5); // same view, but read-only
 * int element_7 = const_view[2]; // set element_7 to 42
 * const_view[2] = 42; // this line won't compile; can't write into this view
 * @endcode
 * In either case, accessing an element of a view does not change the
 * ArrayView object itself, and consequently ArrayView::operator[] is a @p
 * const function. This corresponds to the notion that a view simply
 * represents a, well, "view" of memory that is owned by someone else. Thus,
 * accessing elements of the view changes the memory managed by some other
 * object, but not the view itself, allowing us to make ArrayView::operator[]
 * a @p const member function. This is in contrast to, say, std::vector, which
 * manages the memory it points to and changing an element of the std::vector
 * therefore changes the std::vector object itself -- consequently, the
 * std::vector::operator[] is non-@p const.
 *
 * @note This class is similar to
 *   [`std::span`](https://en.cppreference.com/w/cpp/container/span), but the
 *   latter is only available starting in C++20.
 *
 * @ingroup data
 */
template <typename ElementType, typename MemorySpaceType = MemorySpace::Host>
class ArrayView
{
public:
  /**
   * An alias that denotes the "value_type" of this container-like class,
   * i.e., the type of the element it "stores" or points to.
   */
  using value_type = ElementType;

  /**
   * An alias for iterators pointing into the array.
   */
  using iterator = value_type *;

  /**
   * An alias for const iterators pointing into the array.
   */
  using const_iterator = const ElementType *;

  /**
   * Default constructor.
   */
  ArrayView();

  /**
   * Constructor.
   *
   * @param[in] starting_element A pointer to the first element of the array
   * this object should represent. The value of this argument is only evaluated
   * if `n_elements` is larger than zero. Otherwise, the value of this
   * argument is ignored as if the ArrayView object used a `nullptr`
   * to point to the first element of the array.
   * @param[in] n_elements The length (in elements) of the chunk of memory
   * this object should represent.
   *
   * @note The object that is constructed from these arguments has no
   * knowledge how large the object into which it points really is. As a
   * consequence, whenever you call ArrayView::operator[], the array view can
   * check that the given index is within the range of the view, but it can't
   * check that the view is indeed a subset of the valid range of elements of
   * the underlying object that allocated that range. In other words, you need
   * to ensure that the range of the view specified by the two arguments to
   * this constructor is in fact a subset of the elements of the array into
   * which it points. The appropriate way to do this is to use the
   * make_array_view() functions.
   */
  ArrayView(value_type *starting_element, const std::size_t n_elements);

  /**
   * Copy constructor from array views that point to non-@p const elements. If
   * the current object will point to non-@p const elements, then this is a
   * straight forward copy constructor. On the other hand, if the current
   * type's @p ElementType template argument is a @p const qualified type,
   * then the current constructor is a conversion constructor that converts a
   * non-@p const view to a @p const view, akin to converting a non-@p const
   * pointer to a @p const pointer.
   */
  ArrayView(
    const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType> &view);

  /**
   * A constructor that automatically creates a view from a single value_type
   * object. The view so created then has length one.
   */
  explicit ArrayView(value_type &element);

#ifdef DEAL_II_HAVE_CXX20
  /**
   * A constructor that automatically creates a view from a container
   * object that requires a contiguous array of elements (such as
   * `std::vector`, `std::array`, `boost::container::small_vector`,
   * and the like). The view encompasses all elements of the given
   * container.
   *
   * This implicit conversion constructor is particularly useful when calling
   * a function that takes an ArrayView object as argument, and passing in
   * a container.
   *
   * @note This constructor takes a reference to a @p const container as argument.
   *   It can only be used to initialize ArrayView objects that point to
   *   @p const memory locations, such as <code>ArrayView@<const double@></code>.
   *   You cannot initialize ArrayView objects to non-@p const memory with
   *   such arguments, such as <code>ArrayView@<double@></code>.
   */
  template <typename ContiguousContainer>
  ArrayView(const ContiguousContainer &container)
    requires(std::is_same_v<
               std::remove_cv_t<ElementType>,
               std::remove_cv_t<typename ContiguousContainer::value_type>> &&
             std::is_const_v<ElementType> &&
             concepts::is_contiguous_container<ContiguousContainer>);

  /**
   * A constructor that automatically creates a view from a container
   * object that requires a contiguous array of elements (such as
   * `std::vector`, `std::array`, `boost::container::small_vector`,
   * and the like). The view encompasses all elements of the given
   * container.
   *
   * This implicit conversion constructor is particularly useful when calling
   * a function that takes an ArrayView object as argument, and passing in
   * a container.
   *
   * @note This constructor takes a reference to a non-@p const container as
   *   argument. It can be used to initialize ArrayView objects that point to
   *   either @p const memory locations, such as
   *   <code>ArrayView@<const double@></code>, or to non-@p const memory,
   *   such as <code>ArrayView@<double@></code>.
   */
  template <typename ContiguousContainer>
  ArrayView(ContiguousContainer &container)
    requires(std::is_same_v<
               std::remove_cv_t<ElementType>,
               std::remove_cv_t<typename ContiguousContainer::value_type>> &&
             concepts::is_contiguous_container<ContiguousContainer>);

#else


  /**
   * A constructor that automatically creates a view from a container
   * object that requires a contiguous array of elements (such as
   * `std::vector`, `std::array`, `boost::container::small_vector`,
   * and the like). The view encompasses all elements of the given
   * container.
   *
   * This implicit conversion constructor is particularly useful when calling
   * a function that takes an ArrayView object as argument, and passing in
   * a container.
   *
   * @note This constructor takes a reference to a @p const container as argument.
   *   It can only be used to initialize ArrayView objects that point to
   *   @p const memory locations, such as <code>ArrayView@<const double@></code>.
   *   You cannot initialize ArrayView objects to non-@p const memory with
   *   such arguments, such as <code>ArrayView@<double@></code>.
   */
  template <typename ContiguousContainer,
            typename = decltype(std::data(std::declval<ContiguousContainer>())),
            typename = decltype(std::size(std::declval<ContiguousContainer>())),
            typename = std::enable_if_t<
              std::is_same_v<
                std::remove_cv_t<ElementType>,
                std::remove_cv_t<typename ContiguousContainer::value_type>> &&
              std::is_const_v<ElementType>>>
  ArrayView(const ContiguousContainer &container);

  /**
   * A constructor that automatically creates a view from a container
   * object that requires a contiguous array of elements (such as
   * `std::vector`, `std::array`, `boost::container::small_vector`,
   * and the like). The view encompasses all elements of the given
   * container.
   *
   * This implicit conversion constructor is particularly useful when calling
   * a function that takes an ArrayView object as argument, and passing in
   * a container.
   *
   * @note This constructor takes a reference to a non-@p const container as
   *   argument. It can be used to initialize ArrayView objects that point to
   *   either @p const memory locations, such as
   *   <code>ArrayView@<const double@></code>, or to non-@p const memory,
   *   such as <code>ArrayView@<double@></code>.
   */
  template <typename ContiguousContainer,
            typename = decltype(std::data(std::declval<ContiguousContainer>())),
            typename = decltype(std::size(std::declval<ContiguousContainer>())),
            typename = std::enable_if_t<std::is_same_v<
              std::remove_cv_t<ElementType>,
              std::remove_cv_t<typename ContiguousContainer::value_type>>>>
  ArrayView(ContiguousContainer &container);

#endif

  /**
   * A constructor that automatically creates a view for a given C-style array.
   * This constructor can be used as follows:
   * @code
   *   ArrayView<const int>
   *   get_data_table ()
   *   {
   *     static const int my_data[7] = { 1, 1, 2, 3, 5, 8, 13 };
   *     return {my_data};
   *   }
   * @endcode
   * The object so returned is then a view of the array, with the size 7
   * correctly deduced.
   */
  template <std::size_t N>
  ArrayView(value_type (&array)[N]);

  /**
   * A constructor that creates a view of the array that underlies a
   * [std::initializer_list](https://en.cppreference.com/w/cpp/utility/initializer_list).
   * This constructor allows for cases such as where one has a function
   * that takes an ArrayView object:
   * @code
   *   void f(const ArrayView<const int> &a);
   * @endcode
   * and then to call this function with a list of integers:
   * @code
   *   f({1,2,3});
   * @endcode
   * This also works with an empty list:
   * @code
   *   f({});
   * @encode
   *
   * @note This constructor only works if the template type is `const`
   * qualified. That is, you can initialize an `ArrayView<const int>`
   * object from a `std::initializer_list<int>`, but not an
   * `ArrayView<int>` object. This is because the elements of initializer
   * lists are `const`.
   *
   * @note `std::initializer_list` objects are temporary. They are constructed
   * where the compiler finds a brace-enclosed list, and so they only live
   * for at most the time it takes to execute the current statement. As a
   * consequence, creating an ArrayView object of such a `std::initializer_list`
   * also results in a view object that points to valid memory only for as long
   * as the current statement is executed. You shouldn't expect that the
   * resulting ArrayView can be used to point to useful memory content past
   * that point. In other words, while this code...
   * @code
   *   std::vector<int> v(10);
   *   ArrayView<int> a(v);
   *   f(a);
   * @endcode
   * ...works because the array `v` pointed to exists until after the call to
   * `f()`, the following code will not likely work as expected:
   * @code
   *   ArrayView<int> a({1,2,3});
   *   f(a);
   * @endcode
   */
  ArrayView(
    const std::initializer_list<std::remove_cv_t<value_type>> &initializer_list)
    DEAL_II_CXX20_REQUIRES(std::is_const_v<ElementType>);

  /**
   * Reinitialize a view.
   *
   * @param[in] starting_element A pointer to the first element of the array
   * this object should represent.
   * @param[in] n_elements The length (in elements) of the chunk of memory
   * this object should represent.
   *
   * @note The object that is constructed from these arguments has no
   * knowledge how large the object into which it points really is. As a
   * consequence, whenever you call ArrayView::operator[], the array view can
   * check that the given index is within the range of the view, but it can't
   * check that the view is indeed a subset of the valid range of elements of
   * the underlying object that allocated that range. In other words, you need
   * to ensure that the range of the view specified by the two arguments to
   * this constructor is in fact a subset of the elements of the array into
   * which it points. The appropriate way to do this is to use the
   * make_array_view() functions.
   */
  void
  reinit(value_type *starting_element, const std::size_t n_elements);

  /**
   * Compare two ArrayView objects of the same type. Two objects are considered
   * equal if they have the same size and the same starting pointer.
   *
   * Note that this means that the operation tests that the *views* are the
   * same. If they are, then of course the elements represented by the view
   * are also the same. But the converse is not true: Two ArrayView objects
   * may point to different parts of the memory space and in that case the
   * comparison for equality will return `false` even if the *elements* the
   * views point to are the same.
   *
   * This version always compares with the const value_type.
   */
  bool
  operator==(
    const ArrayView<const value_type, MemorySpaceType> &other_view) const;

  /**
   * Compare two ArrayView objects of the same type. Two objects are considered
   * equal if they have the same size and the same starting pointer.
   *
   * Note that this means that the operation tests that the *views* are the
   * same. If they are, then of course the elements represented by the view
   * are also the same. But the converse is not true: Two ArrayView objects
   * may point to different parts of the memory space and in that case the
   * comparison for equality will return `false` even if the *elements* the
   * views point to are the same.
   *
   * This version always compares with the non-const value_type.
   */
  bool
  operator==(const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType>
               &other_view) const;

  /**
   * Compare two ArrayView objects of the same type. Two objects are considered
   * equal if they have the same size and the same starting pointer, and the
   * current operation therefore returns `true` if the two views being compared
   * point to different memory locations, or if they point to the same memory
   * location but represent different sizes.
   *
   * Note that this means that the operation tests that the *views* are the
   * not the same. But this does not mean that the elements pointed to by
   * the view are not equal: Two ArrayView objects
   * may point to different parts of the memory space and in that case the
   * comparison for inequality will return `true` even if the *elements* the
   * views point to are the same.
   *
   * This version always compares with the const value_type.
   */
  bool
  operator!=(
    const ArrayView<const value_type, MemorySpaceType> &other_view) const;

  /**
   * Compare two ArrayView objects of the same type. Two objects are considered
   * equal if they have the same size and the same starting pointer.
   *
   * Note that this means that the operation tests that the *views* are the
   * not the same. But this does not mean that the elements pointed to by
   * the view are not equal: Two ArrayView objects
   * may point to different parts of the memory space and in that case the
   * comparison for inequality will return `true` even if the *elements* the
   * views point to are the same.
   *
   * This version always compares with the non-const value_type.
   */
  bool
  operator!=(const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType>
               &other_view) const;

  /**
   * Return the size (in elements) of the view of memory this object
   * represents.
   */
  std::size_t
  size() const;

  /**
   * Return a bool whether the array view is empty.
   */
  bool
  empty() const;

  /**
   * Return a pointer to the underlying array serving as element storage.
   * In case the container is empty a nullptr is returned.
   */
  DEAL_II_HOST_DEVICE value_type *
  data() const noexcept;

  /**
   * Return an iterator pointing to the beginning of the array view.
   */
  iterator
  begin() const;

  /**
   * Return an iterator pointing to one past the end of the array view.
   */
  iterator
  end() const;

  /**
   * Return a constant iterator pointing to the beginning of the array view.
   */
  const_iterator
  cbegin() const;

  /**
   * Return a constant iterator pointing to one past the end of the array view.
   */
  const_iterator
  cend() const;

  /**
   * Return a reference to the $i$th element of the range represented by the
   * current object.
   *
   * This function is marked as @p const because it does not change the
   * <em>view object</em>. It may however return a reference to a non-@p const
   * memory location depending on whether the template type of the class is @p
   * const or not.
   *
   * This function is only allowed to be called if the underlying data is indeed
   * stored in CPU memory.
   */
  value_type &
  operator[](const std::size_t i) const;

private:
  /**
   * A pointer to the first element of the range of locations in memory that
   * this object represents.
   */
  value_type *starting_element;

  /**
   * The length of the array this object represents.
   */
  std::size_t n_elements;

  friend class ArrayView<const ElementType, MemorySpaceType>;
};



//---------------------------------------------------------------------------


template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView()
  : starting_element(nullptr)
  , n_elements(0)
{}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  value_type       *starting_element,
  const std::size_t n_elements)
  : // In debug mode, make sure that n_elements>0 and if it is not, set
    // the pointer to a nullptr to trigger segfaults if anyone ever wanted
    // to access elements of the array. In release mode, just take the
    // pointer as given.
  starting_element((library_build_mode == LibraryBuildMode::release) ||
                       (n_elements > 0) ?
                     starting_element :
                     nullptr)
  , n_elements(n_elements)
{}



template <typename ElementType, typename MemorySpaceType>
inline void
ArrayView<ElementType, MemorySpaceType>::reinit(value_type *starting_element,
                                                const std::size_t n_elements)
{
  if constexpr (running_in_debug_mode())
    {
      if (n_elements > 0)
        this->starting_element = starting_element;
      else
        this->starting_element = nullptr;
    }
  else
    {
      this->starting_element = starting_element;
    }
  this->n_elements = n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(ElementType &element)
  : starting_element(&element)
  , n_elements(1)
{}



template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType> &view)
  : starting_element(view.starting_element)
  , n_elements(view.n_elements)
{}


#ifdef DEAL_II_HAVE_CXX20

template <typename ElementType, typename MemorySpaceType>
template <typename ContiguousContainer>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const ContiguousContainer &container)
  requires(std::is_same_v<
             std::remove_cv_t<ElementType>,
             std::remove_cv_t<typename ContiguousContainer::value_type>> &&
           std::is_const_v<ElementType> &&
           concepts::is_contiguous_container<ContiguousContainer>)
  : // use delegating constructor
  ArrayView(container.data(), container.size())
{}



template <typename ElementType, typename MemorySpaceType>
template <typename ContiguousContainer>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  ContiguousContainer &container)
  requires(std::is_same_v<
             std::remove_cv_t<ElementType>,
             std::remove_cv_t<typename ContiguousContainer::value_type>> &&
           concepts::is_contiguous_container<ContiguousContainer>)
  : // use delegating constructor
  ArrayView(std::data(container), std::size(container))
{}


#else

template <typename ElementType, typename MemorySpaceType>
template <typename ContiguousContainer, typename, typename, typename>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const ContiguousContainer &container)
  : // use delegating constructor
  ArrayView(container.data(), container.size())
{}



template <typename ElementType, typename MemorySpaceType>
template <typename ContiguousContainer, typename, typename, typename>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  ContiguousContainer &container)
  : // use delegating constructor
  ArrayView(std::data(container), std::size(container))
{}

#endif



template <typename ElementType, typename MemorySpaceType>
template <std::size_t N>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  ElementType (&array)[N])
  : ArrayView(&array[0], N)
{}


template <typename ElementType, typename MemorySpaceType>
inline ArrayView<ElementType, MemorySpaceType>::ArrayView(
  const std::initializer_list<std::remove_cv_t<value_type>> &initializer)
  DEAL_II_CXX20_REQUIRES(std::is_const_v<ElementType>)
  : // use delegating constructor
  ArrayView(initializer.begin(), initializer.size())
{}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::operator==(
  const ArrayView<const value_type, MemorySpaceType> &other_view) const
{
  return (other_view.data() == starting_element) &&
         (other_view.size() == n_elements);
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::operator==(
  const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType> &other_view)
  const
{
  return (other_view.data() == starting_element) &&
         (other_view.size() == n_elements);
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::operator!=(
  const ArrayView<const value_type, MemorySpaceType> &other_view) const
{
  return !(*this == other_view);
}



template <typename ElementType, typename MemorySpaceType>
inline DEAL_II_HOST_DEVICE
  typename ArrayView<ElementType, MemorySpaceType>::value_type *
  ArrayView<ElementType, MemorySpaceType>::data() const noexcept
{
  if (n_elements == 0)
    return nullptr;
  else
    return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::operator!=(
  const ArrayView<std::remove_cv_t<value_type>, MemorySpaceType> &other_view)
  const
{
  return !(*this == other_view);
}



template <typename ElementType, typename MemorySpaceType>
inline std::size_t
ArrayView<ElementType, MemorySpaceType>::size() const
{
  return n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline bool
ArrayView<ElementType, MemorySpaceType>::empty() const
{
  return n_elements == 0;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::iterator
ArrayView<ElementType, MemorySpaceType>::begin() const
{
  return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::iterator
ArrayView<ElementType, MemorySpaceType>::end() const
{
  return starting_element + n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::const_iterator
ArrayView<ElementType, MemorySpaceType>::cbegin() const
{
  return starting_element;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::const_iterator
ArrayView<ElementType, MemorySpaceType>::cend() const
{
  return starting_element + n_elements;
}



template <typename ElementType, typename MemorySpaceType>
inline typename ArrayView<ElementType, MemorySpaceType>::value_type &
ArrayView<ElementType, MemorySpaceType>::operator[](const std::size_t i) const
{
  AssertIndexRange(i, n_elements);

  return *(starting_element + i);
}



/**
 * A variation of @p ArrayView which allows strided access into the view.
 * This is particularly useful when you want to access only one lane of a
 * VectorizedArray.
 */
template <typename ElementType, std::size_t stride = 1>
class StridedArrayView
{
public:
  /**
   * An alias that denotes the "value_type" of this container-like class,
   * i.e., the type of the element it "stores" or points to.
   */
  using value_type = ElementType;

  /**
   * Constructor.
   *
   * @param[in] starting_element A pointer to the first element of the array
   * this object should represent.
   * @param[in] n_elements The length (in elements) of the chunk of memory
   * this object should represent.
   *
   * @note The object that is constructed from these arguments has no
   * knowledge how large the object into which it points really is. As a
   * consequence, whenever you call ArrayView::operator[], the array view can
   * check that the given index is within the range of the view, but it can't
   * check that the view is indeed a subset of the valid range of elements of
   * the underlying object that allocated that range. In other words, you need
   * to ensure that the range of the view specified by the two arguments to
   * this constructor is in fact a subset of the elements of the array into
   * which it points. The appropriate way to do this is to use the
   * make_array_view() functions.
   */
  StridedArrayView(value_type *starting_element, const std::size_t n_elements);

  /**
   * Return the size (in elements) of the view of memory this object
   * represents.
   */
  std::size_t
  size() const;

  /**
   * Return a bool whether the array view is empty.
   */
  bool
  empty() const;

  /**
   * Return a pointer to the underlying array serving as element storage.
   * In case the container is empty a nullptr is returned.
   */
  value_type *
  data() const noexcept;

  /**
   * Return a reference to the $i$th element of the range represented by the
   * current object.
   *
   * This function is marked as @p const because it does not change the
   * <em>view object</em>. It may however return a reference to a non-@p const
   * memory location depending on whether the template type of the class is @p
   * const or not.
   *
   * This function is only allowed to be called if the underlying data is indeed
   * stored in CPU memory.
   */
  value_type &
  operator[](const std::size_t i) const;

protected:
  /**
   * A pointer to the first element of the range of locations in memory that
   * this object represents.
   */
  value_type *starting_element;

  /**
   * The length of the array this object represents.
   */
  std::size_t n_elements;
};



template <typename ElementType, std::size_t stride>
typename StridedArrayView<ElementType, stride>::value_type &
StridedArrayView<ElementType, stride>::operator[](const std::size_t i) const
{
  AssertIndexRange(i, this->n_elements);

  return *(this->starting_element + stride * i);
}



template <typename ElementType, std::size_t stride>
typename StridedArrayView<ElementType, stride>::value_type *
StridedArrayView<ElementType, stride>::data() const noexcept
{
  if (this->n_elements == 0)
    return nullptr;
  else
    return this->starting_element;
}



template <typename ElementType, std::size_t stride>
bool
StridedArrayView<ElementType, stride>::empty() const
{
  return this->n_elements == 0;
}



template <typename ElementType, std::size_t stride>
std::size_t
StridedArrayView<ElementType, stride>::size() const
{
  return this->n_elements;
}



template <typename ElementType, std::size_t stride>
StridedArrayView<ElementType, stride>::StridedArrayView(
  value_type       *starting_element,
  const std::size_t n_elements)
  : starting_element(starting_element)
  , n_elements(n_elements)
{}



#ifndef DOXYGEN
namespace internal
{
  namespace ArrayViewHelper
  {
    /**
     * Return whether the objects one gets by dereferencing the
     * iterators within the given iterator range form a contiguous
     * range in memory.
     */
    template <class Iterator>
    bool
    is_contiguous(const Iterator &first, const Iterator &last)
    {
      const auto n = std::distance(first, last);
      for (std::decay_t<decltype(n)> i = 0; i < n; ++i)
        if (std::addressof(*(std::next(first, i))) !=
            std::next(std::addressof(*first), i))
          return false;
      return true;
    }


    /**
     * Return whether the objects one gets by dereferencing the
     * iterators within the given iterator range form a contiguous
     * range in memory.
     *
     * This specialization for (@p const or non-@p const) pointers
     * returns @p true unconditionally since the fact that objects
     * pointed to by pointers are contiguous is embedded in the memory
     * model of C++.
     */
    template <class T>
    constexpr bool
    is_contiguous(T *, T *)
    {
      return true;
    }
  } // namespace ArrayViewHelper
} // namespace internal
#endif



/**
 * Create an ArrayView that takes a pair of iterators as arguments. The type
 * of the ArrayView is inferred from the value type of the iterator (e.g., the
 * view created from two const iterators will have a const type).
 *
 * @warning The iterators @p begin and @p end must bound (in the usual half-open
 * way) a contiguous in memory range of values. This function is intended for
 * use with iterators into containers like
 * <code>boost::container::small_vector</code> or <code>std::vector</code> and
 * will not work correctly with, e.g.,
 * <code>boost::container::stable_vector</code> or <code>std::deque</code>.
 * In debug mode, we check that the provided iterators represent contiguous
 * memory indeed.
 *
 * @relatesalso ArrayView
 */
template <typename Iterator, typename MemorySpaceType = MemorySpace::Host>
ArrayView<
  std::remove_reference_t<typename std::iterator_traits<Iterator>::reference>,
  MemorySpaceType>
make_array_view(const Iterator begin, const Iterator end)
{
  static_assert(
    std::is_same_v<typename std::iterator_traits<Iterator>::iterator_category,
                   typename std::random_access_iterator_tag>
#ifdef DEAL_II_HAVE_CXX20
      ||
      std::is_same_v<typename std::iterator_traits<Iterator>::iterator_category,
                     typename std::contiguous_iterator_tag>
#endif
    ,
    "The provided iterator needs to be a random access iterator.");
  Assert(begin <= end,
         ExcMessage(
           "The beginning of the array view needs to be before the end."));
  Assert(internal::ArrayViewHelper::is_contiguous(begin, end),
         ExcMessage("The provided range isn't contiguous in memory!"));
  // the reference type, not the value type, knows the constness of the iterator
  return ArrayView<
    std::remove_reference_t<typename std::iterator_traits<Iterator>::reference>,
    MemorySpaceType>(std::addressof(*begin), end - begin);
}



/**
 * Create a view from a pair of pointers. <code>ElementType</code> may be
 * const-qualified.
 *
 * @warning The pointers @p begin and @p end must bound (in the usual
 * half-open way) a contiguous in memory range of values.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType, typename MemorySpaceType = MemorySpace::Host>
ArrayView<ElementType, MemorySpaceType>
make_array_view(ElementType *const begin, ElementType *const end)
{
  Assert(begin <= end,
         ExcMessage(
           "The beginning of the array view needs to be before the end."));
  return ArrayView<ElementType, MemorySpaceType>(begin, end - begin);
}



/**
 * Create a view from an ArrayView itself.
 *
 * This function is used for @p const references to objects of ArrayView type.
 * It only exists for compatibility purposes.
 *
 * @param[in] array_view The ArrayView that we wish to make a copy of.
 *
 * @relatesalso ArrayView
 */
template <typename Number, typename MemorySpaceType>
inline ArrayView<const Number, MemorySpaceType>
make_array_view(const ArrayView<Number, MemorySpaceType> &array_view)
{
  return make_array_view(array_view.cbegin(), array_view.cend());
}



/**
 * Create a view from an ArrayView itself.
 *
 * This function is used for non-@p const references to objects of ArrayView
 * type. It only exists for compatibility purposes.
 *
 * @param[in] array_view The ArrayView that we wish to make a copy of.
 *
 * @relatesalso ArrayView
 */
template <typename Number, typename MemorySpaceType>
inline ArrayView<Number, MemorySpaceType>
make_array_view(ArrayView<Number, MemorySpaceType> &array_view)
{
  return make_array_view(array_view.begin(), array_view.end());
}



/**
 * Create a view to an entire C-style array. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * Whether the resulting ArrayView is writable or not depends on the
 * ElementType being a const type or not.
 *
 * @param[in] array The C-style array for which we want to have an ArrayView
 * object. The ArrayView corresponds to the <em>entire</em> vector.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType, int N>
inline ArrayView<ElementType>
make_array_view(ElementType (&array)[N])
{
  return ArrayView<ElementType>(array, N);
}



/**
 * Create a view to an entire Vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for non-@p const references to objects of Vector
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] vector The Vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> Vector.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(Vector<ElementType> &vector)
{
  return ArrayView<ElementType>(vector.begin(), vector.size());
}



/**
 * Create a view to an entire Vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for @p const references to objects of Vector type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] vector The Vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> Vector.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Vector<ElementType> &vector)
{
  return ArrayView<const ElementType>(vector.begin(), vector.size());
}



/**
 * Create a view to an entire std::vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for non-@p const references to objects of vector
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> vector.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(std::vector<ElementType> &vector)
{
  return ArrayView<ElementType>(vector.data(), vector.size());
}



/**
 * Create a view to an entire std::vector object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for @p const references to objects of vector type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object. The array view corresponds to the <em>entire</em> vector.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const std::vector<ElementType> &vector)
{
  return ArrayView<const ElementType>(vector.data(), vector.size());
}



/**
 * Create a view to a part of a std::vector object. This is equivalent to
 * initializing the ArrayView object with a pointer to the
 * @p starting_index-th element and the @p size_of_view as the length of the view.
 *
 * This function is used for non-@p const references to objects of vector
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object.
 * @param[in] starting_index The index of the first element of the vector that
 * will be part of this view.
 * @param[in] size_of_view Number of elements in the new ArrayView.
 *
 * @pre <code>starting_index + size_of_view <= vector.size()</code>
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(std::vector<ElementType> &vector,
                const std::size_t         starting_index,
                const std::size_t         size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<ElementType>(&vector[starting_index], size_of_view);
}



/**
 * Create a view to a part of a std::vector object. This is equivalent to
 * initializing the ArrayView object with a pointer to the @p starting_index-th
 * element and the @p size_of_view as the length of the view.
 *
 * This function is used for @p const references to objects of vector type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] vector The vector for which we want to have an array view
 * object.
 * @param[in] starting_index The index of the first element of the vector that
 * will be part of this view.
 * @param[in] size_of_view Number of elements in the new ArrayView.
 *
 * @pre <code>starting_index + size_of_view <= vector.size()</code>
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const std::vector<ElementType> &vector,
                const std::size_t               starting_index,
                const std::size_t               size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<const ElementType>(&vector[starting_index], size_of_view);
}



/**
 * Create a writable view to an entire AlignedVector object. See the
 * documentation of the corresponding overload for std::vector for more
 * information.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(AlignedVector<ElementType> &vector)
{
  return ArrayView<ElementType>(vector.data(), vector.size());
}



/**
 * Create a read-only view to an entire AlignedVector object. See the
 * documentation of the corresponding overload for std::vector for more
 * information.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const AlignedVector<ElementType> &vector)
{
  return ArrayView<const ElementType>(vector.data(), vector.size());
}



/**
 * Create a writable view to a part of an AlignedVector object. See the
 * documentation of the corresponding overload for std::vector for more
 * information.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(AlignedVector<ElementType> &vector,
                const std::size_t           starting_index,
                const std::size_t           size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<ElementType>(&vector[starting_index], size_of_view);
}



/**
 * Create a read-only view to a part of an AlignedVector object. See the
 * documentation of the corresponding overload for std::vector for more
 * information.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const AlignedVector<ElementType> &vector,
                const std::size_t                 starting_index,
                const std::size_t                 size_of_view)
{
  Assert(starting_index + size_of_view <= vector.size(),
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of the given vector."));
  return ArrayView<const ElementType>(&vector[starting_index], size_of_view);
}



/**
 * Create a view to an entire std::array object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for non-@p const references to objects of array
 * type. Such objects contain elements that can be written to. Consequently,
 * the return type of this function is a view to a set of writable objects.
 *
 * @param[in] array The std::array object for which we want to have an array
 * view object. The array view corresponds to the <em>entire</em> array.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType, std::size_t N>
inline ArrayView<ElementType>
make_array_view(std::array<ElementType, N> &array)
{
  return ArrayView<ElementType>(array);
}



/**
 * Create a view to an entire std::array object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element and
 * the size of the given argument.
 *
 * This function is used for @p const references to objects of array type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] array The std::array object for which we want to have an array
 * view object. The array view corresponds to the <em>entire</em> array.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType, std::size_t N>
inline ArrayView<const ElementType>
make_array_view(const std::array<ElementType, N> &array)
{
  return ArrayView<const ElementType>(array);
}



/**
 * Create a view to an entire row of a Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given row, and the length of the row as the length of the view.
 *
 * This function is used for non-@p const references to objects of Table type.
 * Such objects contain elements that can be written to. Consequently, the
 * return type of this function is a view to a set of writable objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(Table<2, ElementType>                          &table,
                const typename Table<2, ElementType>::size_type row)
{
  AssertIndexRange(row, table.size()[0]);
  return ArrayView<ElementType>(&table[row][0], table.size()[1]);
}



/**
 * Create a view to an entire Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given table, and the number of table entries as the length of the view.
 *
 * This function is used for non-@p const references to objects of Table type.
 * Such objects contain elements that can be written to. Consequently, the
 * return type of this function is a view to a set of writable objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to the <em>entire</em> table but the order in
 * which the entries are presented in the array is an implementation detail
 * and should not be relied upon.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(Table<2, ElementType> &table)
{
  return ArrayView<ElementType>(&table[0][0], table.n_elements());
}



/**
 * Create a view to an entire Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given table, and the number of table entries as the length of the view.
 *
 * This function is used for @p const references to objects of Table type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to the <em>entire</em> table but the order in
 * which the entries are presented in the array is an implementation detail
 * and should not be relied upon.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType> &table)
{
  return ArrayView<const ElementType>(&table[0][0], table.n_elements());
}


/**
 * Create a view to an entire LAPACKFullMatrix object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given object, and the number entries as the length of the view.
 *
 * This function is used for @p non-const references to objects of
 * LAPACKFullMatrix type. Such objects contain elements that can be written to.
 * Consequently, the return type of this function is a view to a set of
 * @p non-const objects.
 *
 * @param[in] matrix The LAPACKFullMatrix for which we want to have an array
 * view object. The array view corresponds to the <em>entire</em> object but
 * the order in which the entries are presented in the array is an
 * implementation detail and should not be relied upon.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(LAPACKFullMatrix<ElementType> &matrix)
{
  return ArrayView<ElementType>(&matrix(0, 0), matrix.n_elements());
}



/**
 * Create a view to an entire LAPACKFullMatrix object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given object, and the number of entries as the length of the view.
 *
 * This function is used for @p const references to objects of LAPACKFullMatrix
 * type because they contain immutable elements. Consequently, the return type
 * of this function is a view to a set of @p const objects.
 *
 * @param[in] matrix The LAPACKFullMatrix for which we want to have an array
 * view object. The array view corresponds to the <em>entire</em> object but
 * the order in which the entries are presented in the array is an
 * implementation detail and should not be relied upon.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const LAPACKFullMatrix<ElementType> &matrix)
{
  return ArrayView<const ElementType>(&matrix(0, 0), matrix.n_elements());
}



/**
 * Create a view to an entire row of a Table<2> object. This is equivalent to
 * initializing an ArrayView object with a pointer to the first element of the
 * given row, and the length of the row as the length of the view.
 *
 * This function is used for @p const references to objects of Table type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType>                    &table,
                const typename Table<2, ElementType>::size_type row)
{
  AssertIndexRange(row, table.size()[0]);
  return ArrayView<const ElementType>(&table[row][0], table.size()[1]);
}



/**
 * Create a view to (a part of) a row of a Table<2> object.
 *
 * This function is used for non-@p const references to objects of Table type.
 * Such objects contain elements that can be written to. Consequently, the
 * return type of this function is a view to a set of writable objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 * @param[in] starting_column The index of the column into the given row of
 * the table that corresponds to the first element of this view.
 * @param[in] size_of_view The number of elements this view should have. This
 * corresponds to the number of columns in the current row to which the view
 * should correspond.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<ElementType>
make_array_view(Table<2, ElementType>                          &table,
                const typename Table<2, ElementType>::size_type row,
                const typename Table<2, ElementType>::size_type starting_column,
                const std::size_t                               size_of_view)
{
  AssertIndexRange(row, table.size()[0]);
  AssertIndexRange(starting_column, table.size()[1]);
  Assert(starting_column + size_of_view <= table.size()[1],
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of a column of the given table."));
  return ArrayView<ElementType>(&table[row][starting_column], size_of_view);
}



/**
 * Create a view to (a part of) a row of a Table<2> object.
 *
 * This function is used for @p const references to objects of Table type
 * because they contain immutable elements. Consequently, the return type of
 * this function is a view to a set of @p const objects.
 *
 * @param[in] table The Table for which we want to have an array view object.
 * The array view corresponds to an <em>entire</em> row.
 * @param[in] row The index of the row into the table to which this view
 * should correspond.
 * @param[in] starting_column The index of the column into the given row of
 * the table that corresponds to the first element of this view.
 * @param[in] size_of_view The number of elements this view should have. This
 * corresponds to the number of columns in the current row to which the view
 * should correspond.
 *
 * @relatesalso ArrayView
 */
template <typename ElementType>
inline ArrayView<const ElementType>
make_array_view(const Table<2, ElementType>                    &table,
                const typename Table<2, ElementType>::size_type row,
                const typename Table<2, ElementType>::size_type starting_column,
                const std::size_t                               size_of_view)
{
  AssertIndexRange(row, table.size()[0]);
  AssertIndexRange(starting_column, table.size()[1]);
  Assert(starting_column + size_of_view <= table.size()[1],
         ExcMessage("The starting index and size of the view you want to "
                    "create would lead to a view that extends beyond the end "
                    "of a column of the given table."));
  return ArrayView<const ElementType>(&table[row][starting_column],
                                      size_of_view);
}



/*
 * Create a view that doesn't allow the container it points to to be modified.
 * This is useful if the object passed in is not `const` already and a function
 * requires a view to constant memory in its signature.
 *
 * This function returns an object of type `ArrayView<const T>` where `T` is the
 * element type of the container.
 *
 * @relatesalso ArrayView
 */
template <typename Container>
inline auto
make_const_array_view(const Container &container)
  -> decltype(make_array_view(container))
{
  return make_array_view(container);
}


DEAL_II_NAMESPACE_CLOSE

#endif
