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
#ifndef __deal2__mg_base_h
#define __deal2__mg_base_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/vector.h>
#include <vector>


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

				     //////
    unsigned int get_minlevel () const;
    unsigned int get_maxlevel () const;
    
  private:
				     /**
				      * Level of first component.
				      */
    unsigned int minlevel;

				     /**
				      * Array of the objects to be held.
				      */
    std::vector<Object> objects;
};


/**
 * Multilevel matrix base. This class sets up the interface needed by
 * multilevel algorithms. It has no relation to the actual matrix type
 * and takes the vector class as only template argument.
 *
 * Usually, the derived class @ref{MGMatrix}, operating on an
 * @ref{MGLevelObject} of matrices will be sufficient for applications.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGMatrixBase : public Subscriptor
{
  public:
				   /*
				    * Virtual destructor.
				    */
  virtual ~MGMatrixBase();

				   /**
				    * Matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Adding matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;

				   /**
				    * Adding transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const = 0;
};


/**
 * Base class for coarse grid solvers.  This defines the virtual
 * parenthesis operator, being the interface used by multigrid
 * methods. Any implementation will be done by derived classes.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGCoarseGrid : public Subscriptor
{
  public:
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGCoarseGrid ();

				     /**
				      * Solver method implemented by
				      * derived classes.
				      */
    virtual void operator() (const unsigned int   level,
			     VECTOR       &dst,
			     const VECTOR &src) const = 0;
};


/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is an abstract one and has no implementations of
 * possible algorithms for these operations.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2002
 */
 template <class VECTOR>
 class MGTransfer : public Subscriptor
 {
   public:
 				     /**
 				      * Destructor. Does nothing here, but
 				      * needs to be declared virtual anyway.
 				      */
     virtual ~MGTransfer();

 				     /**
 				      * Prolongate a vector from level
 				      * @p{to_level-1} to level
 				      * @p{to_level}. The previous
 				      * content of @p{dst} is
 				      * overwritten.
 				      *
 				      * @p{src} is assumed to be a vector with
 				      * as many elements as there are degrees
 				      * of freedom on the coarser level of
 				      * the two involved levels, while @p{src}
 				      * shall have as many elements as there
 				      * are degrees of freedom on the finer
 				      * level.
 				      */
     virtual void prolongate (const unsigned int to_level,
 			     VECTOR&            dst,
 			     const VECTOR&      src) const = 0;

 				     /**
 				      * Restrict a vector from level
 				      * @p{from_level} to level
 				      * @p{from_level-1} and add this
 				      * restriction to
 				      * @p{dst}. Obviously, if the
 				      * refined region on level
 				      * @p{from_level} is smaller than
 				      * that on level @p{from_level-1},
 				      * some degrees of freedom in
 				      * @p{dst} are not covered and will
 				      * not be altered. For the other
 				      * degress of freedom, the result
 				      * of the restriction is added.
 				      *
 				      * @p{src} is assumed to be a vector with
 				      * as many elements as there are degrees
 				      * of freedom on the finer level of
 				      * the two involved levels, while @p{src}
 				      * shall have as many elements as there
 				      * are degrees of freedom on the coarser
 				      * level.
 				      */
     virtual void restrict_and_add (const unsigned int from_level,
 				   VECTOR&            dst,
 				   const VECTOR&      src) const = 0;
 };



/**
 * Base class for multigrid smoothers. Does nothing but defining the
 * interface used by multigrid methods.
 *
 * @author Guido Kanschat, 2002
 */
template <class VECTOR>
class MGSmoother : public Subscriptor
{
  public:
				   /**
				    * Virtual destructor.
				    */
  virtual ~MGSmoother();

				   /**
				    * Smoothing function. This is the
				    * function used in multigrid
				    * methods.
				    */
  virtual void smooth (const unsigned int level,
		       VECTOR&            u,
		       const VECTOR&      rhs) const = 0;  
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
  return objects[i-minlevel];
}


template<class Object>
const Object &
MGLevelObject<Object>::operator[] (const unsigned int i) const
{
  Assert((i>=minlevel) && (i<minlevel+objects.size()),
	 ExcIndexRange (i, minlevel, minlevel+objects.size()));
  return objects[i-minlevel];
}


template<class Object>
void
MGLevelObject<Object>::resize (const unsigned int new_minlevel,
			       const unsigned int new_maxlevel)
{
  Assert (new_minlevel <= new_maxlevel, ExcInternalError());
  objects.clear ();

  minlevel = new_minlevel;
  objects.resize (new_maxlevel - new_minlevel + 1);
}


template<class Object>
void
MGLevelObject<Object>::clear ()
{
  typename std::vector<Object>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    v->clear();  
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

/*----------------------------   mgbase.h     ---------------------------*/

#endif
/*----------------------------   mgbase.h     ---------------------------*/
