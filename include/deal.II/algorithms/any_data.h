// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2014 by the deal.II authors
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

#ifndef __deal2__any_data_h
#define __deal2__any_data_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/named_data.h>

#include <boost/any.hpp>
#include <vector>
#include <algorithm>
#include <typeinfo>

DEAL_II_NAMESPACE_OPEN

/**
 * Store any amount of any type of data accessible by an identifier string.
 *
 * @todo Deprecate access by index after NamedData has been deprecated
 * for long enough, then change to a map.
 */
class AnyData :
  public Subscriptor
{
public:
  /// Default constructor for empty object
  AnyData();

  /// Number of stored data objects.
  unsigned int size() const;

  /// Add a new data object
  template <typename type>
  void add(type entry, const std::string &name);

  /**
   * @brief Merge the data of another NamedData to the end of this object.
   */
  void merge(const AnyData &other);

  /**
   * @brief Access to stored data object by name.
   *
   * Find the object with given name, try to convert it to
   * <tt>type</tt> and return it. This function throws an exception
   * if either the name does not exist or if the conversion
   * fails. If such an exception is not desired, use try_read()
   * instead.
   */
  template <typename type>
  type entry (const std::string &name);

  /**
   * @brief Read-only access to stored data object by name.
   *
   * Find the object with given name, try to convert it to
   * <tt>type</tt> and return it. This function throws an exception
   * if either the name does not exist or if the conversion
   * fails. If such an exception is not desired, use try_read()
   * instead.
   */
  template <typename type>
  const type entry (const std::string &name) const;

  /**
   * @brief Dedicated read only access by name.
   *
   * For a constant object, this function equals entry(). For a
   * non-const object, it forces read only access to the data. In
   * particular, it throws an exception if the object is not found
   * or cannot be converted to type.  If such an exception is not
   * desired, use try_read() instead.
   *
   * @warning Do not use this function for stored objects which are
   * pointers. Use read_ptr() instead!
   */
  template <typename type>
  const type read (const std::string &name) const;

  /**
   * @brief Dedicated read only access by name for pointer data.
   *
   * If the stored data object is a pointer to a constant object, the logic
   * of access becomes fairly complicated. Namely, the standard read
   * function may fail, depending on whether it was a const pointer
   * or a regular pointer. This function fixes the logic and
   * ascertains that the object does not become mutable by accident.
   */
  template <typename type>
  const type *read_ptr (const std::string &name) const;

  /**
   * Perform the same action as read_ptr(), but do not throw an
   * exception if the pointer does not exist. Return a null pointer
   * instead.
   */
  template <typename type>
  const type *try_read_ptr (const std::string &name) const;

  /**
   * @brief Dedicated read only access by name without exceptions.
   *
   * This function tries to find the name in the list and return a
   * pointer to the associated object. If either the name is not
   * found or the object cannot be converted to the return type, a
   * null pointer is returned.
   */
  template <typename type>
  const type *try_read (const std::string &name) const;

  /**
   * Access to stored data object by index.
   */
  template <typename type>
  type entry (const unsigned int i);

  /// Read-only access to stored data object by index.
  template <typename type>
  const type entry (const unsigned int i) const;

  /// Dedicated read only access.
  template <typename type>
  const type read (const unsigned int i) const;

  /// Dedicated read only access to pointer object.
  template <typename type>
  const type *read_ptr (const unsigned int i) const;

  /// Dedicated read only access to pointer object without exception.
  template <typename type>
  const type *try_read_ptr (const unsigned int i) const;

  /// Dedicated read only access without exception.
  template <typename type>
  const type *try_read (const unsigned int i) const;

  /// Name of object at index.
  const std::string &name(const unsigned int i) const;

  /**
   * @brief Find index of a named object
   *
   * Try to find the objecty and return its index in the list. Throw
   * an exception if the object has not been found.
   */
  unsigned int find(const std::string &name) const;

  /**
   * @brief Try to find index of a named object
   *
   * Try to find the objecty and return its index in the
   * list. returns numbers::invalid_unsigned_int if the name was not
   * found.
   */
  unsigned int try_find(const std::string &name) const;

  /// Find out if object is of a certain type
  template <typename type>
  bool is_type(const unsigned int i) const;

  /// List the contents to a stream
  template <class STREAM>
  void list (STREAM &os) const;

  /// Conversion from old NamedData
  template <typename type>
  AnyData(const NamedData<type> &);

  /// An entry with this name does not exist in the AnyData object.
  DeclException1(ExcNameNotFound, std::string &,
                 << "No entry with the name " << arg1 << " exists");

  /// The requested type and the stored type are different
  DeclException2(ExcTypeMismatch,
                 char *, char *,
                 << "The requested type " << arg1
                 << " and the stored type " << arg2
                 << " must coincide");

  /**
   * Exception indicating that a function expected a vector to have a
   * certain name, but NamedData had a different name in that
   * position.
   */
  DeclException2(ExcNameMismatch, int, std::string,
                 << "Name at position " << arg1 << " is not equal to " << arg2);
private:
  /// The stored data
  std::vector<boost::any> data;
  /// The names of the stored data
  std::vector<std::string> names;
};


inline
AnyData::AnyData()
{}


template <typename type>
inline
AnyData::AnyData(const NamedData<type> &other)
{
  for (unsigned int i=0; i<other.size(); ++i)
    add(other(i), other.name(i));
}


unsigned int
inline
AnyData::size () const
{
  AssertDimension(data.size(), names.size());
  return data.size();
}


template <typename type>
inline
type
AnyData::entry (const unsigned int i)
{
  AssertIndexRange(i, size());
  type *p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::entry (const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p==0 )
    p = boost::any_cast<const type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type *
AnyData::read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type *>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type *).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type *
AnyData::try_read_ptr(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type *>(&data[i]);
  return *p;
}


template <typename type>
inline
const type *
AnyData::try_read(const unsigned int i) const
{
  AssertIndexRange(i, size());
  const type *p = boost::any_cast<type>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type>(&data[i]);
  return p;
}


inline
const std::string &
AnyData::name(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return names[i];
}


inline
unsigned int
AnyData::try_find(const std::string &n) const
{
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);

  if (it == names.end())
    return numbers::invalid_unsigned_int;

  return it - names.begin();
}


inline
unsigned int
AnyData::find(const std::string &n) const
{
  const unsigned int i = try_find(n);
  Assert(i != numbers::invalid_unsigned_int, ExcNameNotFound(n));

  return i;
}


template <typename type>
inline
bool
AnyData::is_type(const unsigned int i) const
{
  return data[i].type() == typeid(type);
}


template <typename type>
inline
type
AnyData::entry (const std::string &n)
{
  const unsigned int i = find(n);
  type *p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::entry (const std::string &n) const
{
  const unsigned int i = find(n);
  const type *p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::read(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type *
AnyData::read_ptr(const std::string &n) const
{
  const unsigned int i = find(n);
  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type *>(&data[i]);
  Assert(p != 0,
         ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type *
AnyData::try_read_ptr(const std::string &n) const
{
  const unsigned int i = try_find(n);
  if (i == numbers::invalid_unsigned_int)
    return 0;

  const type *const *p = boost::any_cast<type *>(&data[i]);
  if (p==0)
    p = boost::any_cast<const type *>(&data[i]);
  return *p;
}


template <typename type>
inline
const type *
AnyData::try_read(const std::string &n) const
{
  // Try to find name
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);
  // Return null pointer if not found
  if (it == names.end())
    return 0;

  // Compute index and return casted pointer
  unsigned int i=it-names.begin();
  const type *p = boost::any_cast<type>(&data[i]);
  return p;
}


template <typename type>
inline
void
AnyData::add(type ent, const std::string &n)
{
  boost::any e = ent;
  data.push_back(e);
  names.push_back(n);
}


inline
void
AnyData::merge(const AnyData &other)
{
  for (unsigned int i=0; i<other.size(); ++i)
    {
      names.push_back(other.names[i]);
      data.push_back(other.data[i]);
    }
}


template <class STREAM>
inline
void AnyData::list(STREAM &os) const
{
  for (unsigned int i=0; i<names.size(); ++i)
    {
      os << i
         << '\t' << names[i]
         << '\t' << data[i].type().name()
         << std::endl;
    }
}


//----------------------------------------------------------------------//



DEAL_II_NAMESPACE_CLOSE

#endif










