// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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

#ifndef dealii__mg_level_object_h
#define dealii__mg_level_object_h

#include <deal.II/base/subscriptor.h>
#include <vector>

#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN


/**
 * This class represents an array with one object for each used level of a
 * multilevel hierarchy, for example for use in the multigrid algorithms.
 * In contrast to just a generic <code>std::vector</code>, this class allows
 * to store objects only between some minimal and maximal index (=level),
 * as one often wants to run a multilevel algorithm only on a subset of
 * the levels of a mesh (e.g., because the second or third coarsest level is
 * already small enough that it is cheaper to run a direct solver there,
 * rather than recurse to even coarser levels). Despite storing objects only
 * for these "interesting" levels, the class allows indexing simply by
 * level. Internally, this is of course done by
 * simply shifting the given index by the minimum level we have stored.
 *
 * In a typical use case for this class, the objects stored on each level
 * are either matrices or vectors.
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
   * Constructor. Create a multilevel object with given minimal and
   * maximal level, and allocate storage for objects on
   * <code>maxlevel-minlevel+1</code> levels.
   *
   * @note Unlike in many other places of the library, the two arguments
   * here do not denote the first level and last-plus-one level, but indeed
   * an <i>inclusive</i> range of levels for which to allocate storage
   * for level objects. Consequently, the defaults for the two arguments
   * will create an array with one level object, rather than an empty
   * array.
   *
   * @param[in] minlevel The lowest level for which to provision memory
   *   for level objects.
   * @param[in] maxlevel The highest level for which to provision memory
   *   for level objects.
   *
   * @pre minlevel <= maxlevel
   */
  MGLevelObject (const unsigned int minlevel = 0,
                 const unsigned int maxlevel = 0);

  /**
   * Access object on level @p level.
   */
  Object &operator[] (const unsigned int level);

  /**
   * Access object on level @p level.
   *
   * This function can be called on a @p const object, and
   * consequently returns a @p const reference.
   */
  const Object &operator[] (const unsigned int level) const;

  /**
   * Delete all previous contents of this object and reset its size according
   * to the values of @p new_minlevel and @p new_maxlevel.
   *
   * @param[in] new_minlevel The lowest level for which to provision memory
   *   for level objects.
   * @param[in] new_maxlevel The highest level for which to provision memory
   *   for level objects.
   *
   * @pre minlevel <= maxlevel
   */
  void resize (const unsigned int new_minlevel,
               const unsigned int new_maxlevel);

  /**
   * Call <tt>operator = (s)</tt> on all objects stored by this object.
   * This clearly requires that the objects stored on each level allow for
   * this operation. This is, in particular, true for vectors and matrices
   * if @p d is zero, thereby zeroing out all vector or matrix entries.
   */
  MGLevelObject<Object> &operator = (const double d);

  /**
   * Call @p clear on all objects stored by this object. This function
   * is only implemented for some @p Object classes, e.g., matrix
   * types or the PreconditionBlockSOR and similar classes. Using this
   * function will fail with a compiler error if the @p Object
   * template type to this class does not provide a
   * <code>clear()</code> member function.
   */
  void clear();

  /**
   * The coarsest level for which this class stores a level object.
   */
  unsigned int min_level () const;

  /**
   * The highest level for which this class stores a level object.
   */
  unsigned int max_level () const;

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
