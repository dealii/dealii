// ---------------------------------------------------------------------
//
// Copyright (C) 2021 by the deal.II authors
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

#ifndef dealii_memory_block_h
#define dealii_memory_block_h

#include <deal.II/base/config.h>

#include <deal.II/base/array_view.h>
#include <deal.II/base/memory_space_utils.h>

#include <type_traits>

DEAL_II_NAMESPACE_OPEN

/**
 * This class allocates a block of memory on the host or the device. The
 * elements of the block must be accessed using ArrayView. Note that when a
 * reinit() function is called, the underlying pointer is changed and thus, one
 * needs to call reinit() on the ArrayView associated with the reinitialized
 * MemoryBlock. The template parameter `ElementType` needs to be trivial.
 */
template <typename ElementType, typename MemorySpaceType>
class MemoryBlock
{
public:
  static_assert(std::is_trivial<ElementType>::value,
                "MemoryBlock only supports trivial ElementType");

  /**
   * An alias that denotes the memory space of this container-like class.
   */
  using memory_space = MemorySpaceType;

  /**
   * Default constructor.
   */
  MemoryBlock() = default;

  /**
   * Constructor. Allocate a block of @p size elements. The data is not
   * initialized.
   */
  MemoryBlock(const unsigned int size);

  /**
   * Copy constructor. A new block of memory is allocated and the data is
   * copied.
   */
  MemoryBlock(const MemoryBlock<ElementType, MemorySpaceType> &other);

  /**
   * Copy the data stored in @p other and move it to the appropriate memory space.
   */
  template <typename MemorySpaceType2>
  MemoryBlock(const MemoryBlock<ElementType, MemorySpaceType2> &other);

  /**
   * Copy the data stored in @p array_view and move it to the appropriate memory space.
   */
  template <typename MemorySpaceType2>
  MemoryBlock(const ArrayView<ElementType, MemorySpaceType2> &array_view);

  /**
   * Release the memory block and allocate a new block. The data is not
   * initialized.
   */
  void
  reinit(const unsigned int size);

  /**
   * Release the memory block, allocate a new block, and copy the data stored in @p other.
   */
  template <typename MemorySpaceType2>
  void
  reinit(const MemoryBlock<ElementType, MemorySpaceType2> &other);

  /**
   * Release the memory block, allocate a new block, and copy the elements in @p
   * array_view.
   */
  template <typename MemorySpaceType2>
  void
  reinit(const ArrayView<ElementType, MemorySpaceType2> &array_view);

  /**
   * Destructor.
   */
  ~MemoryBlock();

  /**
   * Free the allocated memory.
   */
  void
  clear();

  /**
   * Sets all the elements to @p s. This operation is only allowed if @p s is equal to zero.
   */
  MemoryBlock<ElementType, MemorySpaceType> &
  operator=(const ElementType s);

  /**
   * Return the number of elements in the memory block.
   */
  unsigned int
  size() const;

  /**
   * Return a pointer to the underlying array serving as element storage.
   * In case the container is empty a nullptr is returned.
   */
  ElementType *
  data() const;

private:
  template <typename ElementType2, typename MemorySpaceType2>
  friend class MemoryBlock;

  /**
   * Number of elements stored.
   */
  unsigned int n_elements = 0;
  /**
   * Pointer to block of allocated block of memory.
   */
  ElementType *values = nullptr;
};


namespace internal
{
  namespace MemoryBlock
  {
#ifdef DEAL_II_COMPILER_CUDA_AWARE
    template <typename ElementType>
    void
    set_memory_block_to_zero(ElementType *     data_ptr,
                             const std::size_t size,
                             const MemorySpace::CUDA &)
    {
      cudaError_t const error_code =
        cudaMemset(data_ptr, 0, size * sizeof(ElementType));
      AssertCuda(error_code);
    }
#endif

    template <typename ElementType>
    void
    set_memory_block_to_zero(ElementType *     data_ptr,
                             const std::size_t size,
                             const MemorySpace::Host &)
    {
      std::memset(data_ptr, 0, size * sizeof(ElementType));
    }
  } // namespace MemoryBlock
} // namespace internal


template <typename ElementType, typename MemorySpaceType>
MemoryBlock<ElementType, MemorySpaceType>::MemoryBlock(const unsigned int size)
  : n_elements(size)
  , values(Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
      n_elements))
{}



template <typename ElementType, typename MemorySpaceType>
MemoryBlock<ElementType, MemorySpaceType>::MemoryBlock(
  const MemoryBlock<ElementType, MemorySpaceType> &other)
  : n_elements(other.n_elements)
  , values(Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
      n_elements))
{
  Utilities::MemorySpace::copy(
    other.values, MemorySpaceType{}, values, MemorySpaceType{}, n_elements);
}



template <typename ElementType, typename MemorySpaceType>
template <typename MemorySpaceType2>
MemoryBlock<ElementType, MemorySpaceType>::MemoryBlock(
  const MemoryBlock<ElementType, MemorySpaceType2> &other)
  : n_elements(other.n_elements)
  , values(Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
      n_elements))
{
  Utilities::MemorySpace::copy(
    other.values, MemorySpaceType2{}, values, MemorySpaceType{}, n_elements);
}



template <typename ElementType, typename MemorySpaceType>
template <typename MemorySpaceType2>
MemoryBlock<ElementType, MemorySpaceType>::MemoryBlock(
  const ArrayView<ElementType, MemorySpaceType2> &array_view)
  : n_elements(array_view.size())
  , values(Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
      n_elements))
{
  Utilities::MemorySpace::copy(array_view.data(),
                               MemorySpaceType2{},
                               values,
                               MemorySpaceType{},
                               n_elements);
}



template <typename ElementType, typename MemorySpaceType>
void
MemoryBlock<ElementType, MemorySpaceType>::reinit(const unsigned int size)
{
  clear();
  n_elements = size;
  values = Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
    n_elements);
}



template <typename ElementType, typename MemorySpaceType>
template <typename MemorySpaceType2>
void
MemoryBlock<ElementType, MemorySpaceType>::reinit(
  const MemoryBlock<ElementType, MemorySpaceType2> &other)
{
  clear();
  n_elements = other.n_elements;
  values = Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
    n_elements);
  Utilities::MemorySpace::copy(
    other.values, MemorySpaceType2{}, values, MemorySpaceType{}, n_elements);
}



template <typename ElementType, typename MemorySpaceType>
template <typename MemorySpaceType2>
void
MemoryBlock<ElementType, MemorySpaceType>::reinit(
  const ArrayView<ElementType, MemorySpaceType2> &array_view)
{
  clear();
  n_elements = array_view.size();
  values = Utilities::MemorySpace::allocate_data<ElementType, MemorySpaceType>(
    n_elements);
  Utilities::MemorySpace::copy(array_view.data(),
                               MemorySpaceType2{},
                               values,
                               MemorySpaceType{},
                               n_elements);
}



template <typename ElementType, typename MemorySpaceType>
MemoryBlock<ElementType, MemorySpaceType>::~MemoryBlock()
{
  clear();
}



template <typename ElementType, typename MemorySpaceType>
void
MemoryBlock<ElementType, MemorySpaceType>::clear()
{
  if (values != nullptr)
    {
      n_elements = 0;
      Utilities::MemorySpace::delete_data<ElementType, MemorySpaceType>(values);
      values = nullptr;
    }
}



template <typename ElementType, typename MemorySpaceType>
MemoryBlock<ElementType, MemorySpaceType> &
MemoryBlock<ElementType, MemorySpaceType>::operator=(const ElementType s)
{
  Assert(s == ElementType(),
         ExcMessage("Only 0 can be assigned to a MemoryBlock."));
  internal::MemoryBlock::set_memory_block_to_zero(values,
                                                  n_elements,
                                                  MemorySpaceType{});
}



template <typename ElementType, typename MemorySpaceType>
unsigned int
MemoryBlock<ElementType, MemorySpaceType>::size() const
{
  return n_elements;
}



template <typename ElementType, typename MemorySpaceType>
ElementType *
MemoryBlock<ElementType, MemorySpaceType>::data() const
{
  return values;
}

DEAL_II_NAMESPACE_CLOSE

#endif
