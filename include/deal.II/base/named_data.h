// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__named_data_h
#define __deal2__named_data_h

#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>

#include <vector>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

class AnyData;

/**
 * @deprecated The use of this class is deprecated and AnyData should
 * be used instead. It is not only more flexible, the way constness
 * was assured in this class never worked the way it was
 * intended. Therefore, the according checks have been disabled.
 *
 * This class is a collection of DATA objects. Each of the pointers
 * has a name associated, enabling identification by this name rather
 * than the index in the vector only.
 *
 * Note that it is the actual data stored in this object. Therefore,
 * for storing vectors and other large objects, it should be
 * considered to use SmartPointer or boost::shared_ptr for DATA.
 *
 * Objects of this kind have a smart way of treating constness: if a
 * const data object is added or a const NamedData is supplied with
 * merge(), the object will henceforth consider itself as const
 * (#is_constant will be true). Thus, any subsequent modification will
 * be illegal and ExcConstantObject will be raised in debug mode.
 *
 * @author Guido Kanschat, 2007, 2008, 2009
 */
template <typename DATA>
class NamedData : public Subscriptor
{
public:
  /** Standard constructor creating
   * an empty object.
   */
  NamedData();

  /**
   * Assignment operator, copying
   * conversible data from another
   * object.
   */
  template <typename DATA2>
  NamedData<DATA> &operator = (const NamedData<DATA2> &other);

  /**
   * \name Adding members
   */
//@{
  /**
   * Add a new data item to the end
   * of the collection.
   */
  void add(DATA &v, const std::string &name);

  /**
   * Add a new constant data item
   * to the end of the collection
   * and make the collection
   * constant.
   */
  void add(const DATA &v, const std::string &name);

  /**
   * Merge the data of another
   * NamedData to the end of this
   * object.
   *
   * If the other object had
   * #is_constant set, so will have
   * this object after merge.
   */
  template <typename DATA2>
  void merge(NamedData<DATA2> &);
  /**
   * Merge the data of another
   * NamedData to the end of this
   * object.
   *
   * After this operation, all
   * data in this object will be
   * treated as const.
   */
  template <typename DATA2>
  void merge(const NamedData<DATA2> &);
//@}
  /**
   * \name Accessing and querying
   * contents
   */
//@{
/// Number of stored data objects.
  unsigned int size() const;

  /**
   * @brief Access to stored data
   * object by index.
   *
   * @note This function throws an
   * exception, if is_const()
   * returns <tt>true</tt>. In such
   * a case, either cast the
   * NamedData object to a const
   * reference, or use the function
   * read() instead of this operator.
   */
  DATA &operator() (unsigned int i);

/// Read-only access to stored data object by index.
  const DATA &operator() (unsigned int i) const;

/// Read only access for a non-const object.
  const DATA &read (unsigned int i) const;

/// Name of object at index.
  const std::string &name(unsigned int i) const;

/// Find index of a named object
  unsigned int find(const std::string &name) const;

/// Returns true if this object contains constant data
  bool is_const () const;

/// List names of stored objects
  template <class OUT>
  void print(OUT &o) const;
//@}
  /**
   * Exception indicating that a
   * function expected a vector
   * to have a certain name, but
   * NamedData had a different
   * name in that position.
   */
  DeclException2(ExcNameMismatch, int, std::string,
                 << "Name at position " << arg1 << " is not equal to " << arg2);

  /**
   * Exception indicating that read
   * access to stored data was
   * attempted although the
   * NamedData object contains
   * const data and #is_constant
   * was true.
   */
  DeclException0(ExcConstantObject);

private:
  /// True if the object is to be treated constant
  bool is_constant;
  /// The actual data stored
  std::vector<DATA> data;

  /// Names for the data
  std::vector<std::string> names;
};


/**
 * Select data from NamedData corresponding to the attached name.
 *
 * Given a list of names to search for (provided by add()), objects of
 * this class provide an index list of the selected data.
 *
 * @author Guido Kanschat, 2009
 */
class NamedSelection
{
public:
  /**
   * Add a new name to be searched
   * for in NamedData.
   *
   * @note Names will be added to
   * the end of the current list.
   */
  void add (const std::string &name);

  /**
   * Create the index vector
   * pointing into the NamedData
   * object.
   */
  template <typename DATA>
  void initialize(const NamedData<DATA> &data);

  /**
   * Create the index vector
   * pointing into the AnyData
   * object.
   */
  void initialize(const AnyData &data);

  /**
   * The number of names in this
   * object. This function may be
   * used whether initialize() was
   * called before or not.
   */
  unsigned int size() const;

  /**
   * Return the corresponding index
   * in the NamedData object
   * supplied to the last
   * initialize(). It is an error
   * if initialize() has not been
   * called before.
   *
   * Indices are in the same order
   * as the calls to add().
   */
  unsigned int operator() (unsigned int i) const;

private:
  /**
   * The selected names.
   */
  std::vector<std::string> names;
  /**
   * The index map generated by
   * initialize() and accessed by
   * operator().
   */
  std::vector<unsigned int> indices;
};


//----------------------------------------------------------------------//

template<typename DATA>
inline
NamedData<DATA>::NamedData()
  :
  is_constant(false)
{}


template<typename DATA>
inline
void
NamedData<DATA>::add(DATA &v, const std::string &n)
{
  // see operator() below

  //  Assert(!is_constant, ExcConstantObject());
  names.push_back(n);
  data.push_back(v);
}


template<typename DATA>
inline
void
NamedData<DATA>::add(const DATA &v, const std::string &n)
{
  // see operator() below

  //  Assert(!is_constant, ExcConstantObject());
  DATA &aux = const_cast<DATA &>(v);
  data.push_back(aux);
  names.push_back(n);
  is_constant = true;
}


template<typename DATA>
template<typename DATA2>
inline
void
NamedData<DATA>::merge(NamedData<DATA2> &other)
{
  // see operator() below

  // Assert(!is_constant, ExcConstantObject());

  for (unsigned int i=0; i<other.size(); ++i)
    {
      names.push_back(other.name(i));
      data.push_back(other.read(i));
    }
  is_constant = other.is_const();
}


template<typename DATA>
template<typename DATA2>
inline
void
NamedData<DATA>::merge(const NamedData<DATA2> &other)
{
  // see operator() below

  //  Assert(!is_constant, ExcConstantObject());
  for (unsigned int i=0; i<other.size(); ++i)
    {
      names.push_back(other.name(i));
      data.push_back(other(i));
    }
  is_constant = true;
}



template<typename DATA>
template<typename DATA2>
inline
NamedData<DATA> &
NamedData<DATA>::operator= (const NamedData<DATA2> &other)
{
  is_constant = false;
  merge(other);
  is_constant = other.is_const();
  return *this;
}



template<typename DATA>
inline
unsigned int
NamedData<DATA>::size() const
{
  return data.size();
}


template<typename DATA>
inline
bool
NamedData<DATA>::is_const() const
{
  return is_constant;
}


template<typename DATA>
inline
DATA &
NamedData<DATA>::operator() (unsigned int i)
{
  // Remove this assertion, since is_const() does not really assert
  // what we would like to. This can only be fixed by using constness
  // in AnyData.

  //  Assert(!is_constant, ExcConstantObject());
  AssertIndexRange(i, size());
  return data[i];
}


template<typename DATA>
inline
const DATA &
NamedData<DATA>::operator() (unsigned int i) const
{
  AssertIndexRange(i, size());
  return data[i];
}


template<typename DATA>
inline
const DATA &
NamedData<DATA>::read (unsigned int i) const
{
  AssertIndexRange(i, size());
  return data[i];
}


template<typename DATA>
inline
const std::string &
NamedData<DATA>::name(unsigned int i) const
{
  AssertIndexRange(i, size());
  return names[i];
}


template<typename DATA>
inline
unsigned int
NamedData<DATA>::find (const std::string &name) const
{
  const std::vector<std::string>::const_iterator
  i = std::find(names.begin(), names.end(), name);
  if (i == names.end())
    return deal_II_numbers::invalid_unsigned_int;
  return i - names.begin();
}


template<typename DATA>
template<class OUT>
inline
void
NamedData<DATA>::print(OUT &o) const
{
  o << "NamedData:";
  for (unsigned int i=0; i<size(); ++i)
    o << ' ' << '\"' << names[i] << '\"';
  o << std::endl;
}


inline
void
NamedSelection::add(const std::string &s)
{
  names.push_back(s);
}


template <typename DATA>
inline void
NamedSelection::initialize(const NamedData<DATA> &data)
{
  indices.resize(names.size());
  for (unsigned int i=0; i<names.size(); ++i)
    indices[i] = data.find(names[i]);
}


inline
unsigned int
NamedSelection::size() const
{
  return names.size();
}


inline
unsigned int
NamedSelection::operator() (unsigned int i) const
{
  Assert (indices.size() == names.size(), ExcNotInitialized());
  AssertIndexRange(i, size());
  return indices[i];
}


DEAL_II_NAMESPACE_CLOSE


#endif
