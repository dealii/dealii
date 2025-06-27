// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_any_data_h
#define dealii_any_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/types.h>

#include <algorithm>
#include <any>
#include <ostream>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * Store any amount of any type of data accessible by an identifier string.
 *
 * @todo GK: Deprecate access to AnyData by index and change to a map.
 */
class AnyData : public EnableObserverPointer
{
public:
  /// Default constructor for empty object
  AnyData() = default;

  /// Number of stored data objects.
  unsigned int
  size() const;

  /// Add a new data object
  template <typename type>
  void
  add(type entry, const std::string &name);

  /**
   * @brief Merge the data of another AnyData to the end of this object.
   */
  void
  merge(const AnyData &other);

  /**
   * @brief Access to stored data object by name.
   *
   * Find the object with given name, try to convert it to <tt>type</tt> and
   * return it. This function throws an exception if either the name does not
   * exist or if the conversion fails. If such an exception is not desired,
   * use try_read() instead.
   */
  template <typename type>
  type
  entry(const std::string &name);

  /**
   * @brief Read-only access to stored data object by name.
   *
   * Find the object with given name, try to convert it to <tt>type</tt> and
   * return it. This function throws an exception if either the name does not
   * exist or if the conversion fails. If such an exception is not desired,
   * use try_read() instead.
   */
  template <typename type>
  type
  entry(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name.
   *
   * For a constant object, this function equals entry(). For a non-const
   * object, it forces read only access to the data. In particular, it throws
   * an exception if the object is not found or cannot be converted to type.
   * If such an exception is not desired, use try_read() instead.
   *
   * @warning Do not use this function for stored objects which are pointers.
   * Use read_ptr() instead!
   */
  template <typename type>
  type
  read(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name for pointer data.
   *
   * If the stored data object is a pointer to a constant object, the logic of
   * access becomes fairly complicated. Namely, the standard read function may
   * fail, depending on whether it was a const pointer or a regular pointer.
   * This function fixes the logic and ascertains that the object does not
   * become mutable by accident.
   */
  template <typename type>
  const type *
  read_ptr(const std::string &name) const;

  /**
   * Perform the same action as read_ptr(), but do not throw an exception if
   * the pointer does not exist. Return a null pointer instead.
   */
  template <typename type>
  const type *
  try_read_ptr(const std::string &name) const;

  /**
   * @brief Dedicated read only access by name without exceptions.
   *
   * This function tries to find the name in the list and return a pointer to
   * the associated object. If either the name is not found or the object
   * cannot be converted to the return type, a null pointer is returned.
   */
  template <typename type>
  const type *
  try_read(const std::string &name) const;

  /**
   * Access to stored data object by index.
   */
  template <typename type>
  type
  entry(const unsigned int i);

  /// Read-only access to stored data object by index.
  template <typename type>
  type
  entry(const unsigned int i) const;

  /// Dedicated read only access.
  template <typename type>
  type
  read(const unsigned int i) const;

  /// Dedicated read only access to pointer object.
  template <typename type>
  const type *
  read_ptr(const unsigned int i) const;

  /// Dedicated read only access to pointer object without exception.
  template <typename type>
  const type *
  try_read_ptr(const unsigned int i) const;

  /// Dedicated read only access without exception.
  template <typename type>
  const type *
  try_read(const unsigned int i) const;

  /// Name of object at index.
  const std::string &
  name(const unsigned int i) const;

  /**
   * @brief Find index of a named object
   *
   * Try to find the object and return its index in the list. Throw an
   * exception if the object has not been found.
   */
  unsigned int
  find(const std::string &name) const;

  /**
   * @brief Try to find index of a named object
   *
   * Try to find the object and return its index in the list. returns
   * numbers::invalid_unsigned_int if the name was not found.
   */
  unsigned int
  try_find(const std::string &name) const;

  /// Find out if object is of a certain type
  template <typename type>
  bool
  is_type(const unsigned int i) const;

  /// List the contents to a stream
  template <typename StreamType>
  void
  list(StreamType &os) const;

  /// An entry with this name does not exist in the AnyData object.
  DeclException1(ExcNameNotFound,
                 std::string,
                 << "No entry with the name " << arg1 << " exists.");

  /// The requested type and the stored type are different
  DeclException2(ExcTypeMismatch,
                 std::string,
                 std::string,
                 << "The requested type " << arg1 << " and the stored type "
                 << arg2 << " must coincide.");

  /**
   * Exception indicating that a function expected a vector to have a certain
   * name, but we store a different name in that position.
   */
  DeclException2(ExcNameMismatch,
                 int,
                 std::string,
                 << "Name at position " << arg1 << " is not equal to " << arg2
                 << '.');

private:
  /// The stored data
  std::vector<std::any> data;
  /// The names of the stored data
  std::vector<std::string> names;
};


unsigned int inline AnyData::size() const
{
  AssertDimension(data.size(), names.size());
  return data.size();
}


template <typename type>
inline type
AnyData::entry(const unsigned int i)
{
  AssertIndexRange(i, size());
  const type *p = std::any_cast<type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline type
AnyData::entry(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = std::any_cast<type>(&data[i]);
  if (p == nullptr)
    p = std::any_cast<const type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline type
AnyData::read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = std::any_cast<type>(&data[i]);
  if (p == nullptr)
    p = std::any_cast<const type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = std::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = std::any_cast<const type *>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type *).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = std::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = std::any_cast<const type *>(&data[i]);
  if (p == nullptr)
    return nullptr;
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = std::any_cast<type>(&data[i]);
  if (p == 0)
    p = std::any_cast<const type>(&data[i]);
  return p;
}


inline const std::string &
AnyData::name(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return names[i];
}


inline unsigned int
AnyData::try_find(const std::string &n) const
{
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);

  if (it == names.end())
    return numbers::invalid_unsigned_int;

  return it - names.begin();
}


inline unsigned int
AnyData::find(const std::string &n) const
{
  const unsigned int i = try_find(n);
  Assert(i != numbers::invalid_unsigned_int, ExcNameNotFound(n));

  return i;
}


template <typename type>
inline bool
AnyData::is_type(const unsigned int i) const
{
  return data[i].type() == typeid(type);
}


template <typename type>
inline type
AnyData::entry(const std::string &n)
{
  const unsigned int i = find(n);
  const type        *p = std::any_cast<type>(&data[i]);
  Assert(p != 0, ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline type
AnyData::entry(const std::string &n) const
{
  const unsigned int i = find(n);
  const type        *p = std::any_cast<type>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline type
AnyData::read(const std::string &n) const
{
  const unsigned int i = find(n);
  const type        *p = std::any_cast<type>(&data[i]);
  Assert(p != 0, ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::read_ptr(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *const *p = std::any_cast<type *>(&data[i]);
  if (p == nullptr)
    p = std::any_cast<const type *>(&data[i]);
  Assert(p != nullptr,
         ExcTypeMismatch(typeid(type).name(), data[i].type().name()));
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read_ptr(const std::string &n) const
{
  const unsigned int i = try_find(n);
  if (i == numbers::invalid_unsigned_int)
    return 0;

  const type *const *p = std::any_cast<type *>(&data[i]);
  if (p == 0)
    p = std::any_cast<const type *>(&data[i]);
  return *p;
}


template <typename type>
inline const type *
AnyData::try_read(const std::string &n) const
{
  // Try to find name
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);
  // Return null pointer if not found
  if (it == names.end())
    return nullptr;

  // Compute index and return casted pointer
  unsigned int i = it - names.begin();
  const type  *p = std::any_cast<type>(&data[i]);
  return p;
}


template <typename type>
inline void
AnyData::add(type ent, const std::string &n)
{
  std::any e = ent;
  data.push_back(e);
  names.push_back(n);
}


inline void
AnyData::merge(const AnyData &other)
{
  for (unsigned int i = 0; i < other.size(); ++i)
    {
      names.push_back(other.names[i]);
      data.push_back(other.data[i]);
    }
}


template <typename StreamType>
inline void
AnyData::list(StreamType &os) const
{
  for (unsigned int i = 0; i < names.size(); ++i)
    {
      os << i << '\t' << names[i] << '\t' << data[i].type().name() << std::endl;
    }
}


//----------------------------------------------------------------------//



DEAL_II_NAMESPACE_CLOSE

#endif
