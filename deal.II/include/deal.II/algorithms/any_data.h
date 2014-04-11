// ---------------------------------------------------------------------
// $Id: named_data.h 30036 2013-07-18 16:55:32Z maier $
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
    /// Number of stored data objects.
    unsigned int size() const;
    
    /// Add a new data object
    template <typename type>
    void add(type entry, const std::string& name);
  
    /**
     * @brief Merge the data of another NamedData to the end of this object.
     */
    void merge(const AnyData& other);

    /**
     * @brief Access to stored data object by name.
     */
    template <typename type>
    type entry (const std::string& name);

    /// Read-only access to stored data object by name.
    template <typename type>
    const type entry (const std::string& name) const;

    /// Dedicated read only access by name.
    template <typename type>
    const type read (const std::string& name) const;
    
  /**
   * @brief Access to stored data object by index.
   */
    template <typename type>
    type entry (const unsigned int i);

    /// Read-only access to stored data object by index.
    template <typename type>
    const type entry (const unsigned int i) const;

    /// Dedicated read only access.
    template <typename type>
    const type read (const unsigned int i) const;
    
    /// Name of object at index.
    const std::string &name(const unsigned int i) const;

    /// Find index of a named object
    unsigned int find(const std::string &name) const;

    /// Find out if object is of a certain type
    template <typename type>
    bool is_type(const unsigned int i) const;
    
    /// The requested type and the stored type are different
  DeclException2(ExcTypeMismatch,
		 char*, char*,
		 << "The requested type " << arg1
		 << " and the stored type " << arg2
		 << " must coincide");
  
  private:
    /// The stored data
    std::vector<boost::any> data;
    /// The names of the stored data
    std::vector<std::string> names;
};


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
  type* p = boost::any_cast<type>(&data[i]);
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
  const type* p = boost::any_cast<type>(&data[i]);
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
  const type* p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
	 ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


inline
const std::string&
AnyData::name(const unsigned int i) const
{
  AssertIndexRange(i, size());
  return names[i];
}


inline
unsigned int
AnyData::find(const std::string& n) const
{
  std::vector<std::string>::const_iterator it =
    std::find(names.begin(), names.end(), n);
  
  Assert(it != names.end(), ExcMessage("An entry with this name does not exist"));
  
  return it - names.begin();
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
AnyData::entry (const std::string& n)
{
  const unsigned int i = find(n);
  type* p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
	 ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::entry (const std::string& n) const
{
  const unsigned int i = find(n);
  const type* p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
	 ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
const type
AnyData::read(const std::string& n) const
{
  const unsigned int i = find(n);
  const type* p = boost::any_cast<type>(&data[i]);
  Assert(p != 0,
	 ExcTypeMismatch(typeid(type).name(),data[i].type().name()));
  return *p;
}


template <typename type>
inline
void
AnyData::add(type ent, const std::string& n)
{
  boost::any e = ent;
  data.push_back(e);
  names.push_back(n);
}


inline
void
AnyData::merge(const AnyData& other)
{
  for (unsigned int i=0; i<other.size(); ++i)
    {
      names.push_back(other.names[i]);
      data.push_back(other.data[i]);
    }
}


//----------------------------------------------------------------------//



DEAL_II_NAMESPACE_CLOSE

#endif
