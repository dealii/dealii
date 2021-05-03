// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2018 by the deal.II authors
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

#ifndef dealii_tria_iterator_templates_h
#define dealii_tria_iterator_templates_h


#include <deal.II/base/config.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

DEAL_II_NAMESPACE_OPEN


/* Note: This file only contains template definitions and will thus
   not produce an object file. It is rather thought to be included
   into the *_accessor.cc files.
*/


/*------------------------ Functions: TriaRawIterator ------------------*/


template <typename Accessor>
inline TriaRawIterator<Accessor>::TriaRawIterator()
  : accessor(nullptr, -2, -2, nullptr)
{}


template <typename Accessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const TriaRawIterator<Accessor> &i)
  : accessor(i.accessor)
{}



template <typename Accessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
  const int                                                            level,
  const int                                                            index,
  const typename Accessor::AccessorData *local_data)
  : accessor(parent, level, index, local_data)
{}


template <typename Accessor>
inline TriaRawIterator<Accessor>::TriaRawIterator(
  const TriaAccessorBase<Accessor::structure_dimension,
                         Accessor::dimension,
                         Accessor::space_dimension> &tria_accessor,
  const typename Accessor::AccessorData *            local_data)
  : accessor(nullptr, -2, -2, local_data)
{
  accessor.copy_from(tria_accessor);
}


template <typename Accessor>
inline TriaRawIterator<Accessor> &
TriaRawIterator<Accessor>::operator=(const TriaRawIterator<Accessor> &i)
{
  accessor.copy_from(i.accessor);

  return *this;
}



template <typename Accessor>
template <typename OtherAccessor>
inline
  typename std::enable_if<std::is_convertible<OtherAccessor, Accessor>::value,
                          bool>::type
  TriaRawIterator<Accessor>::
  operator==(const TriaRawIterator<OtherAccessor> &other) const
{
  return accessor == other.accessor;
}


template <typename Accessor>
inline bool
TriaRawIterator<Accessor>::
operator!=(const TriaRawIterator<Accessor> &other) const
{
  return !(*this == other);
}


template <typename Accessor>
inline TriaRawIterator<Accessor>
TriaRawIterator<Accessor>::operator++(int)
{
  TriaRawIterator<Accessor> tmp(*this);
                            operator++();

  return tmp;
}


template <typename Accessor>
inline TriaRawIterator<Accessor>
TriaRawIterator<Accessor>::operator--(int)
{
  TriaRawIterator<Accessor> tmp(*this);
                            operator--();

  return tmp;
}


/*-----------------------  functions: TriaIterator ---------------*/


template <typename Accessor>
inline TriaIterator<Accessor>::TriaIterator()
  : TriaRawIterator<Accessor>()
{}


template <typename Accessor>
inline TriaIterator<Accessor>::TriaIterator(const TriaIterator<Accessor> &i)
  : TriaRawIterator<Accessor>(i.accessor)
{}


template <typename Accessor>
inline TriaIterator<Accessor>::TriaIterator(const TriaRawIterator<Accessor> &i)
  : TriaRawIterator<Accessor>(i.accessor)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#endif
}


template <typename Accessor>
inline TriaIterator<Accessor>::TriaIterator(
  const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
  const int                                                            level,
  const int                                                            index,
  const typename Accessor::AccessorData *local_data)
  : TriaRawIterator<Accessor>(parent, level, index, local_data)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#endif
}


template <typename Accessor>
inline TriaIterator<Accessor>::TriaIterator(
  const TriaAccessorBase<Accessor::structure_dimension,
                         Accessor::dimension,
                         Accessor::space_dimension> &tria_accessor,
  const typename Accessor::AccessorData *            local_data)
  : TriaRawIterator<Accessor>(tria_accessor, local_data)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#endif
}


template <typename Accessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator=(const TriaIterator<Accessor> &i)
{
  this->accessor.copy_from(i.accessor);
  return *this;
}


template <typename Accessor>
template <typename OtherAccessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator=(const TriaIterator<OtherAccessor> &i)
{
  this->accessor.copy_from(i.accessor);
  return *this;
}


template <typename Accessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator=(const TriaRawIterator<Accessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#endif
  return *this;
}


template <typename Accessor>
template <typename OtherAccessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator=(const TriaRawIterator<OtherAccessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || used)"
  // used() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used(), ExcAssignmentOfUnusedObject());
#endif
  return *this;
}


template <typename Accessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator++()
{
  while (TriaRawIterator<Accessor>::operator++(),
         (this->state() == IteratorState::valid))
    if (this->accessor.used() == true)
      return *this;
  return *this;
}


template <typename Accessor>
inline TriaIterator<Accessor>
TriaIterator<Accessor>::operator++(int)
{
  TriaIterator<Accessor> tmp(*this);
                         operator++();

  return tmp;
}


template <typename Accessor>
inline TriaIterator<Accessor> &
TriaIterator<Accessor>::operator--()
{
  while (TriaRawIterator<Accessor>::operator--(),
         (this->state() == IteratorState::valid))
    if (this->accessor.used() == true)
      return *this;
  return *this;
}


template <typename Accessor>
inline TriaIterator<Accessor>
TriaIterator<Accessor>::operator--(int)
{
  TriaIterator<Accessor> tmp(*this);
                         operator--();

  return tmp;
}


/*-----------------------  functions: TriaActiveIterator ---------------*/


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator()
  : TriaIterator<Accessor>()
{}


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaActiveIterator<Accessor> &i)
  : TriaIterator<Accessor>(static_cast<TriaIterator<Accessor>>(i))
{}


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaRawIterator<Accessor> &i)
  : TriaIterator<Accessor>(i)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
}


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaIterator<Accessor> &i)
  : TriaIterator<Accessor>(i)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
}


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const Triangulation<Accessor::dimension, Accessor::space_dimension> *parent,
  const int                                                            level,
  const int                                                            index,
  const typename Accessor::AccessorData *local_data)
  : TriaIterator<Accessor>(parent, level, index, local_data)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
}


template <typename Accessor>
inline TriaActiveIterator<Accessor>::TriaActiveIterator(
  const TriaAccessorBase<Accessor::structure_dimension,
                         Accessor::dimension,
                         Accessor::space_dimension> &tria_accessor,
  const typename Accessor::AccessorData *            local_data)
  : TriaIterator<Accessor>(tria_accessor, local_data)
{
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
}


template <typename Accessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator=(const TriaActiveIterator<Accessor> &i)
{
  this->accessor.copy_from(i.accessor);
  return *this;
}


template <typename Accessor>
template <class OtherAccessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::
operator=(const TriaActiveIterator<OtherAccessor> &i)
{
  this->accessor.copy_from(i.accessor);
  return *this;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator=(const TriaRawIterator<Accessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used() && this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
  return *this;
}


template <typename Accessor>
template <class OtherAccessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator=(const TriaRawIterator<OtherAccessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.used() && this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
  return *this;
}


template <typename Accessor>
template <class OtherAccessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator=(const TriaIterator<OtherAccessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
  return *this;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator=(const TriaIterator<Accessor> &i)
{
  this->accessor.copy_from(i.accessor);
#ifdef DEBUG
  // do this like this, because:
  // if we write
  // "Assert (IteratorState::past_the_end || !has_children())"
  // has_children() is called anyway, even if
  // state==IteratorState::past_the_end, and will then
  // throw the exception!
  if (this->state() != IteratorState::past_the_end)
    Assert(this->accessor.has_children() == false,
           ExcAssignmentOfInactiveObject());
#endif
  return *this;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator++()
{
  while (TriaIterator<Accessor>::operator++(),
         (this->state() == IteratorState::valid))
    if (this->accessor.has_children() == false)
      return *this;
  return *this;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor>
TriaActiveIterator<Accessor>::operator++(int)
{
  TriaActiveIterator<Accessor> tmp(*this);
                               operator++();

  return tmp;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor> &
TriaActiveIterator<Accessor>::operator--()
{
  while (TriaIterator<Accessor>::operator--(),
         (this->state() == IteratorState::valid))
    if (this->accessor.has_children() == false)
      return *this;
  return *this;
}


template <typename Accessor>
inline TriaActiveIterator<Accessor>
TriaActiveIterator<Accessor>::operator--(int)
{
  TriaActiveIterator<Accessor> tmp(*this);
                               operator--();

  return tmp;
}

DEAL_II_NAMESPACE_CLOSE

#endif
