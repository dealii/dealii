//----------------------------  mg_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  mg_base.h  ---------------------------
#ifndef __deal2__mg_base_h
#define __deal2__mg_base_h


/*----------------------------   mgbase.h     ---------------------------*/


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
    typename std::vector<Object> objects;
};



/**
 * Coarse grid solver using LAC iterative methods.
 * This is a little wrapper, transforming a triplet of iterative
 * solver, matrix and preconditioner into a coarse grid solver.
 *
 * The type of the matrix (i.e. the template parameter @p{MATRIX})
 * should be derived from @p{Subscriptor} to allow for the use of a
 * smart pointer to it.
 *
 * @author Guido Kanschat, 1999
 */
template<class SOLVER, class MATRIX, class PRECOND, class VECTOR = Vector<double> >
class MGCoarseGridLACIteration :  public Subscriptor
{
  public:
				     /**
				      * Constructor.
				      * Store solver, matrix and
				      * preconditioning method for later
				      * use.
				      */
    MGCoarseGridLACIteration (SOLVER        &,
			      const MATRIX  &,
			      const PRECOND &);
    
				     /**
				      * Implementation of the abstract
				      * function.
				      * Calls the solver method with
				      * matrix, vectors and
				      * preconditioner.
				      */
    void operator() (const unsigned int   level,
		     VECTOR       &dst,
		     const VECTOR &src) const;
  private:
				     /**
				      * Reference to the solver.
				      */
    SOLVER& solver;
    
				     /**
				      * Reference to the matrix.
				      */
    const SmartPointer<const MATRIX> matrix;
    
				     /**
				      * Reference to the preconditioner.
				      */
    const PRECOND& precondition;
};


/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is an abstract one and has no implementations of
 * possible algorithms for these operations.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
//  template <class VECTOR>
//  class MGTransferBase : public Subscriptor
//  {
//    public:
//  				     /**
//  				      * Destructor. Does nothing here, but
//  				      * needs to be declared virtual anyway.
//  				      */
//      virtual ~MGTransferBase();

//  				     /**
//  				      * Prolongate a vector from level
//  				      * @p{to_level-1} to level
//  				      * @p{to_level}. The previous
//  				      * content of @p{dst} is
//  				      * overwritten.
//  				      *
//  				      * @p{src} is assumed to be a vector with
//  				      * as many elements as there are degrees
//  				      * of freedom on the coarser level of
//  				      * the two involved levels, while @p{src}
//  				      * shall have as many elements as there
//  				      * are degrees of freedom on the finer
//  				      * level.
//  				      */
//      virtual void prolongate (const unsigned int to_level,
//  			     VECTOR&            dst,
//  			     const VECTOR&      src) const = 0;

//  				     /**
//  				      * Restrict a vector from level
//  				      * @p{from_level} to level
//  				      * @p{from_level-1} and add this
//  				      * restriction to
//  				      * @p{dst}. Obviously, if the
//  				      * refined region on level
//  				      * @p{from_level} is smaller than
//  				      * that on level @p{from_level-1},
//  				      * some degrees of freedom in
//  				      * @p{dst} are not covered and will
//  				      * not be altered. For the other
//  				      * degress of freedom, the result
//  				      * of the restriction is added.
//  				      *
//  				      * @p{src} is assumed to be a vector with
//  				      * as many elements as there are degrees
//  				      * of freedom on the finer level of
//  				      * the two involved levels, while @p{src}
//  				      * shall have as many elements as there
//  				      * are degrees of freedom on the coarser
//  				      * level.
//  				      */
//      virtual void restrict_and_add (const unsigned int from_level,
//  				   VECTOR&            dst,
//  				   const VECTOR&      src) const = 0;
//  };



/* ------------------------------------------------------------------- */


template<class Object>
MGLevelObject<Object>::MGLevelObject(const unsigned int min,
				     const unsigned int max)
		:
		minlevel(0)
{
  resize (min, max);
};


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
};


template<class Object>
void
MGLevelObject<Object>::clear ()
{
  typename std::vector<Object>::iterator v;
  for (v = objects.begin(); v != objects.end(); ++v)
    v->clear();  
};


template<class Object>
unsigned int
MGLevelObject<Object>::get_minlevel () const
{
  return minlevel;
};


template<class Object>
unsigned int
MGLevelObject<Object>::get_maxlevel () const
{
  return minlevel + objects.size() - 1;
};


/* ------------------ Functions for MGCoarseGridLACIteration ------------ */


template<class SOLVER, class MATRIX, class PRECOND, class VECTOR>
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND, VECTOR>
::MGCoarseGridLACIteration(SOLVER& s,
			   const MATRIX  &m,
			   const PRECOND &p)
		:
		solver(s),
		matrix(&m),
		precondition(p)
{};


template<class SOLVER, class MATRIX, class PRECOND, class VECTOR>
void
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND, VECTOR>
::operator() (const unsigned int    /* level */,
	      VECTOR       &dst,
	      const VECTOR &src) const
{
  solver.solve(*matrix, dst, src, precondition);
}


/*----------------------------   mgbase.h     ---------------------------*/

#endif
/*----------------------------   mgbase.h     ---------------------------*/
