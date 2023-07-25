// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2022 by the deal.II authors
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

#ifndef dealii_mg_level_object_h
#define dealii_mg_level_object_h

#include <deal.II/base/config.h>

#include <deal.II/base/subscriptor.h>

#include <memory>
#include <vector>

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
 */
template <class Object>
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
   * @param[in] args Optional arguments passed to the constructor of the
   *   underlying object.
   *
   * @pre minlevel <= maxlevel
   */
  template <class... Args>
  MGLevelObject(const unsigned int minlevel,
                const unsigned int maxlevel,
                Args &&...args);

  /**
   * Constructor. Same as above but without arguments to be forwarded to the
   * constructor of the underlying object.
   */
  MGLevelObject(const unsigned int minlevel = 0,
                const unsigned int maxlevel = 0);

  /**
   * Access object on level @p level.
   */
  Object &
  operator[](const unsigned int level);

  /**
   * Access object on level @p level.
   *
   * This function can be called on a @p const object, and
   * consequently returns a @p const reference.
   */
  const Object &
  operator[](const unsigned int level) const;

  /**
   * Return object on level max.
   */
  const Object &
  back() const;

  /**
   * Delete all previous contents of this object and reset its size according
   * to the values of @p new_minlevel and @p new_maxlevel.
   *
   * @param[in] new_minlevel The lowest level for which to provision memory
   *   for level objects.
   * @param[in] new_maxlevel The highest level for which to provision memory
   *   for level objects.
   * @param[in] args Optional arguments passed to the constructor of the
   *   underlying object.
   *
   * @pre minlevel <= maxlevel
   */
  template <class... Args>
  void
  resize(const unsigned int new_minlevel,
         const unsigned int new_maxlevel,
         Args &&...args);

  /**
   * Call <tt>operator = (s)</tt> on all objects stored by this object.
   * This clearly requires that the objects stored on each level allow for
   * this operation. This is, in particular, true for vectors and matrices
   * if @p d is zero, thereby zeroing out all vector or matrix entries.
   */
  MGLevelObject<Object> &
  operator=(const double d);

  /**
   * Clear all data fields and brings the class into a condition similar
   * to after having called the default constructor.
   */
  void
  clear();

  /**
   * Call @p clear on all objects stored by this object. This function
   * is only implemented for some @p Object classes, e.g., matrix
   * types or the PreconditionBlockSOR and similar classes. Using this
   * function will fail with a compiler error if the @p Object
   * template type to this class does not provide a
   * <code>clear()</code> member function.
   */
  void
  clear_elements();

  /**
   * The coarsest level for which this class stores a level object.
   */
  unsigned int
  min_level() const;

  /**
   * The highest level for which this class stores a level object.
   */
  unsigned int
  max_level() const;

  /**
   * Number of levels, i.e., `max_level()-min_level()+1`.
   */
  unsigned int
  n_levels() const;

  /**
   * Apply the action @p action to every object stored in here. The
   * parameter @p action is expected to be a function object that accepts
   * the syntax
   * <code>
   *   action(const unsigned int level, Object &object);
   * </code>
   * This means this function can accept a lambda, a std::function, or a plain
   * function pointer.
   */
  template <typename ActionFunctionObjectType>
  void
  apply(ActionFunctionObjectType action);

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * Level of first component.
   */
  unsigned int minlevel;

  /**
   * Array of the objects to be held.
   */
  std::vector<std::shared_ptr<Object>> objects;
};


/* ------------------------------------------------------------------- */


template <class Object>
template <class... Args>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
                                     const unsigned int max,
                                     Args &&...args)
  : minlevel(0)
{
  resize(min, max, std::forward<Args>(args)...);
}


template <class Object>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
                                     const unsigned int max)
  : minlevel(0)
{
  resize(min, max);
}


template <class Object>
Object &
MGLevelObject<Object>::operator[](const unsigned int i)
{
  Assert((i >= minlevel) && (i < minlevel + objects.size()),
         ExcIndexRange(i, minlevel, minlevel + objects.size()));
  return *objects[i - minlevel];
}


template <class Object>
const Object &
MGLevelObject<Object>::operator[](const unsigned int i) const
{
  Assert((i >= minlevel) && (i < minlevel + objects.size()),
         ExcIndexRange(i, minlevel, minlevel + objects.size()));
  return *objects[i - minlevel];
}


template <class Object>
const Object &
MGLevelObject<Object>::back() const
{
  return this->operator[](this->max_level());
}


template <class Object>
template <class... Args>
void
MGLevelObject<Object>::resize(const unsigned int new_minlevel,
                              const unsigned int new_maxlevel,
                              Args &&...args)
{
  Assert(new_minlevel <= new_maxlevel, ExcInternalError());
  // note that on clear(), the
  // shared_ptr class takes care of
  // deleting the object it points to
  // by itself
  objects.clear();

  minlevel = new_minlevel;
  for (unsigned int i = 0; i < new_maxlevel - new_minlevel + 1; ++i)
    objects.push_back(std::make_shared<Object>(std::forward<Args>(args)...));
}


template <class Object>
MGLevelObject<Object> &
MGLevelObject<Object>::operator=(const double d)
{
  typename std::vector<std::shared_ptr<Object>>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    **v = d;
  return *this;
}


template <class Object>
void
MGLevelObject<Object>::clear()
{
  minlevel = 0;
  objects.clear();
}


template <class Object>
void
MGLevelObject<Object>::clear_elements()
{
  typename std::vector<std::shared_ptr<Object>>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    (*v)->clear();
}


template <class Object>
unsigned int
MGLevelObject<Object>::min_level() const
{
  return minlevel;
}


template <class Object>
unsigned int
MGLevelObject<Object>::max_level() const
{
  return minlevel + objects.size() - 1;
}


template <class Object>
unsigned int
MGLevelObject<Object>::n_levels() const
{
  return objects.size();
}

template <class Object>
template <typename ActionFunctionObjectType>
void
MGLevelObject<Object>::apply(ActionFunctionObjectType action)
{
  for (unsigned int lvl = min_level(); lvl <= max_level(); ++lvl)
    {
      action(lvl, (*this)[lvl]);
    }
}


template <class Object>
std::size_t
MGLevelObject<Object>::memory_consumption() const
{
  std::size_t result = sizeof(*this);
  using Iter = typename std::vector<std::shared_ptr<Object>>::const_iterator;
  const Iter end = objects.end();
  for (Iter o = objects.begin(); o != end; ++o)
    result += (*o)->memory_consumption();

  return result;
}

DEAL_II_NAMESPACE_CLOSE

#endif
