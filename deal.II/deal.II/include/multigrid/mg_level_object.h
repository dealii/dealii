//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------
#ifndef __deal2__mg_level_object_h
#define __deal2__mg_level_object_h

#include <base/subscriptor.h>
#include <vector>


#include <boost_local/shared_ptr.hpp>


/**
 * An array with an object for each level.  The purpose of this class
 * is mostly to allow access by level number, even if the lower levels
 * are not used and therefore have no object at all; this is done by
 * simply shifting the given index by the minimum level we have
 * stored.
 *
 * In most cases, the objects which are stored on each levels, are
 * either matrices or vectors.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
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
				      * Access object on level @p{level}.
				      */
    Object & operator[] (const unsigned int level);
    
				     /**
				      * Access object on level
				      * @p{level}. Constant version.
				      */
    const Object & operator[] (const unsigned int level) const;

				     /**
				      * Delete all previous contents
				      * of this object and reset its
				      * size according to the values
				      * of @p{new_minlevel} and
				      * @p{new_maxlevel}.
				      */
    void resize (const unsigned int new_minlevel,
		 const unsigned int new_maxlevel);
    
				     /**
				      * Call @p{clear} on all objects
				      * stored by this object. This
				      * function is only implemented
				      * for some @p{Object} classes,
				      * most notably for vectors and
				      * matrices. Note that if
				      * @p{Object==Vector<T>}, @p{clear}
				      * will set all entries to zero,
				      * while if
				      * @p{Object==vector<T>},
				      * @p{clear} deletes the elements
				      * of the vectors. This class
				      * might therefore not be useful
				      * for STL vectors.
				      */
    void clear();

				     /**
				      * @brief Coarsest level for multigrid.
				      */
    unsigned int get_minlevel () const;
    
				     /**
				      * Ignored
				      */
    unsigned int get_maxlevel () const;
    
				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    
  private:
				     /**
				      * Level of first component.
				      */
    unsigned int minlevel;

				     /**
				      * Array of the objects to be held.
				      */
    std::vector<boost::shared_ptr<Object> > objects;
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
    objects.push_back(boost::shared_ptr<Object> (new Object)); 
}


template<class Object>
void
MGLevelObject<Object>::clear ()
{
  typename std::vector<boost::shared_ptr<Object> >::iterator v;
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
MGLevelObject<Object>::memory_consumption () const
{
  unsigned int result = sizeof(*this);
  typedef typename std::vector<boost::shared_ptr<Object> >::const_iterator Iter;
  const Iter end = objects.end();
  for (Iter o=objects.begin(); o!=end; ++o)
    result += *o->memory_consumption();
  
  return result;
}


/*-----------------------   mg_level_object.h     -----------------------*/

#endif
/*-----------------------   mg_level_object.h     -----------------------*/
