// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2013 by the deal.II authors
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

#ifndef __deal2__mg_level_object_h
#define __deal2__mg_level_object_h

#include <deal.II/base/subscriptor.h>
#include <vector>

#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN


/**
 * An array with an object for each level.  The purpose of this class
 * is mostly to store objects and allow access by level number, even
 * if the lower levels are not used and therefore have no object at
 * all; this is done by simply shifting the given index by the minimum
 * level we have stored.
 *
 * In most cases, the objects which are stored on each levels, are
 * either matrices or vectors.
 *
 * @ingroup mg
 * @ingroup data
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2005, 2010
 */
template<class Object>
class MGLevelObject : public Subscriptor
{
public:
  /**
   * Constructor allowing to
   * initialize the number of
   * levels. By default, the object
   * is created empty.
   */
  MGLevelObject (const unsigned int minlevel = 0,
                 const unsigned int maxlevel = 0);

  /**
   * Access object on level @p level.
   */
  Object &operator[] (const unsigned int level);

  /**
   * Access object on level
   * @p level. Constant version.
   */
  const Object &operator[] (const unsigned int level) const;

  /**
   * Delete all previous contents
   * of this object and reset its
   * size according to the values
   * of @p new_minlevel and
   * @p new_maxlevel.
   */
  void resize (const unsigned int new_minlevel,
               const unsigned int new_maxlevel);

  /**
   * Call <tt>operator = (s)</tt>
   * on all objects stored by this
   * object.  This is particularly
   * useful for
   * e.g. <tt>Object==Vector@<T@></tt>
   */
  MGLevelObject<Object> &operator = (const double d);

  /**
   * Call @p clear on all objects
   * stored by this object. This
   * function is only implemented
   * for some @p Object classes,
   * e.g. the PreconditionBlockSOR
   * and similar classes.
   */
  void clear();

  /**
   * Coarsest level for multigrid.
   */
  unsigned int min_level () const;

  /**
   * Finest level for multigrid.
   */
  unsigned int max_level () const;

  /**
   * @deprecated Replaced by min_level()
   */
  unsigned int get_minlevel () const DEAL_II_DEPRECATED;

  /**
   * @deprecated Replaced by max_level()
   */
  unsigned int get_maxlevel () const DEAL_II_DEPRECATED;

  /**
   * Memory used by this object.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Level of first component.
   */
  unsigned int minlevel;

  /**
   * Array of the objects to be held.
   */
  std::vector<std_cxx11::shared_ptr<Object> > objects;
};


/* ------------------------------------------------------------------- */


template<class Object>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
                                     const unsigned int max)
  :
  minlevel(0)
{
  resize (min, max);
}


template<class Object>
Object &
MGLevelObject<Object>::operator[] (const unsigned int i)
{
  Assert((i>=minlevel) && (i<minlevel+objects.size()),
         ExcIndexRange (i, minlevel, minlevel+objects.size()));
  return *objects[i-minlevel];
}


template<class Object>
const Object &
MGLevelObject<Object>::operator[] (const unsigned int i) const
{
  Assert((i>=minlevel) && (i<minlevel+objects.size()),
         ExcIndexRange (i, minlevel, minlevel+objects.size()));
  return *objects[i-minlevel];
}


template<class Object>
void
MGLevelObject<Object>::resize (const unsigned int new_minlevel,
                               const unsigned int new_maxlevel)
{
  Assert (new_minlevel <= new_maxlevel, ExcInternalError());
  // note that on clear(), the
  // shared_ptr class takes care of
  // deleting the object it points to
  // by itself
  objects.clear ();

  minlevel = new_minlevel;
  for (unsigned int i=0; i<new_maxlevel-new_minlevel+1; ++i)
    objects.push_back(std_cxx11::shared_ptr<Object> (new Object));
}


template<class Object>
MGLevelObject<Object> &
MGLevelObject<Object>::operator = (const double d)
{
  typename std::vector<std_cxx11::shared_ptr<Object> >::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    **v=d;
  return *this;
}


template<class Object>
void
MGLevelObject<Object>::clear ()
{
  typename std::vector<std_cxx11::shared_ptr<Object> >::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    (*v)->clear();
}


template<class Object>
unsigned int
MGLevelObject<Object>::get_minlevel () const
{
  return minlevel;
}


template<class Object>
unsigned int
MGLevelObject<Object>::get_maxlevel () const
{
  return minlevel + objects.size() - 1;
}


template<class Object>
unsigned int
MGLevelObject<Object>::min_level () const
{
  return minlevel;
}


template<class Object>
unsigned int
MGLevelObject<Object>::max_level () const
{
  return minlevel + objects.size() - 1;
}


template<class Object>
std::size_t
MGLevelObject<Object>::memory_consumption () const
{
  std::size_t result = sizeof(*this);
  typedef typename std::vector<std_cxx11::shared_ptr<Object> >::const_iterator Iter;
  const Iter end = objects.end();
  for (Iter o=objects.begin(); o!=end; ++o)
    result += (*o)->memory_consumption();

  return result;
}

DEAL_II_NAMESPACE_CLOSE

#endif
