//----------------------------  mg_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
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


#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/vector.h>

#include <vector>

class MGSmootherBase;


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
				      * @p{Object==std::vector<T>},
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
    vector<Object> objects;
};


/**
 * Abstract base class for coarse grid solvers.  The interface of a
 * function call operator is defined to execute coarse grid solution
 * in a derived class.
 *
 * @author Guido Kanschat, 1999
 */
class MGCoarseGridSolver : public Subscriptor
{
  public:
				     /**
				      * Virtual destructor. Does
				      * nothing in particular, but
				      * needs to be declared anyhow.
				      */
    virtual ~MGCoarseGridSolver();
    
				     /**
				      * Coarse grid solving method.
				      * This is only the interface for
				      * a function defined in a derived
				      * class like
				      * @p{MGCoarseGridLACIteration}.
				      *
				      * Note that the information
				      * about the matrix is removed to
				      * that class.
				      */
    virtual void operator() (const unsigned int    level,
			     Vector<double>       &dst,
			     const Vector<double> &src) const = 0;
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
template<class SOLVER, class MATRIX, class PRECOND>
class MGCoarseGridLACIteration :  public MGCoarseGridSolver
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
    virtual void operator() (const unsigned int   level,
			     Vector<double>       &dst,
			     const Vector<double> &src) const;
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
class MGTransferBase : public Subscriptor
{
  public:
				     /**
				      * Destructor. Does nothing here, but
				      * needs to be declared virtual anyway.
				      */
    virtual ~MGTransferBase();

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
    virtual void prolongate (const unsigned int    to_level,
			     Vector<double>       &dst,
			     const Vector<double> &src) const = 0;

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
    virtual void restrict_and_add (const unsigned int    from_level,
				   Vector<double>       &dst,
				   const Vector<double> &src) const = 0;
};


/**
 * Basic class for preconditioning by multigrid.
 *
 * The functionality of the multigrid method is restricted to defect
 * correction. It is @em{not} iterative and the start solution is
 * always zero. Since by this $u^E_l$ and $u^A_l$ (see report on
 * multigrid) are always zero, restriction is simplified a lot and
 * maybe even the seam condition on grids is oblivious. Still, I am
 * not sure that these restrictions on the class might cause numerical
 * inefficiencies.
 *
 * The function @p{precondition} is the actual multigrid method and
 * makes use of several operations to be implemented in derived
 * classes. It takes a defect in @p{src} and the result is the multigrid
 * preconditioned defect in @p{dst}.
 *
 * @author Guido Kanschat, 1999
 */
class MGBase : public Subscriptor
{
  private:
				     /**
				      * Copy constructor. Made private
				      * to prevent use.
				      */
    MGBase(const MGBase&);

				     /**
				      * Copy operator. Made private
				      * to prevent use.
				      */
    const MGBase& operator=(const MGBase&);

  public:
				     /**
				      * Constructor, subject to change.
				      */
    MGBase (const MGTransferBase &transfer,
	    const unsigned int    minlevel,
	    const unsigned int    maxlevel);
    
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGBase();
    
				     /**
				      * Execute one step of the
				      * v-cycle algorithm.  This
				      * function assumes, that the
				      * vector @p{d} is properly filled
				      * with the residual in the outer
				      * defect correction
				      * scheme. After execution of
				      * @p{vcycle()}, the result is in
				      * the vector @p{s}. We propose to
				      * write functions @p{copy_from_mg}
				      * and @p{copy_to_mg} in derived
				      * classes of @p{MGBase} to access
				      * @p{d} and @p{s}.
				      *
				      * The actual work for this
				      * function is done in
				      * @p{level_mgstep}.
				      */
    void vcycle(const MGSmootherBase     &pre_smooth,
		const MGSmootherBase     &post_smooth,
		const MGCoarseGridSolver &cgs);

				     /**
				      * Exception.
				      */
    DeclException2(ExcSwitchedLevels, int, int,
		   << "minlevel and maxlevel switched, should be: "
		   << arg1 << "<=" << arg2);
    
  protected:
				     /**
				      * Highest level of cells.
				      */
    unsigned int maxlevel;

				     /**
				      * Level for coarse grid solution.
				      */
    unsigned int minlevel;
    
				     /**
				      * Auxiliary vector.
				      */
    MGLevelObject<Vector<double> > defect;

				     /**
				      * Auxiliary vector.
				      */
    MGLevelObject<Vector<double> > solution;
    
				     /**
				      * Prolongation and restriction object.
				      */
    SmartPointer<const MGTransferBase> transfer;

				     /**
				      * The actual v-cycle multigrid method.
				      * @p{level} is the level the
				      * function should work on. It
				      * will usually be called for the
				      * highest level from outside,
				      * but will then call itself
				      * recursively for @p{level-1},
				      * unless we are on @p{minlevel}
				      * where instead of recursing
				      * deeper, the coarse grid solver
				      * is used to solve the matrix of
				      * this level.
				      */
    void level_mgstep (const unsigned int        level,
		       const MGSmootherBase     &pre_smooth,
		       const MGSmootherBase     &post_smooth,
		       const MGCoarseGridSolver &cgs);  

				     /**
				      * Apply negative @p{vmult} on all
				      * cells of a level.
				      * This is implemented in a
				      * derived class.
				      */
    virtual void level_vmult (const unsigned int    level,
			      Vector<double>       &dst,
			      const Vector<double> &src,
			      const Vector<double> &rhs) = 0;

				     /**
				      * Print a level vector using
				      * @ref{DoFHandler}.
				      */
    virtual void print_vector (const unsigned int level,
			       const Vector<double>& v,
			       const char* name) const = 0;
    
  private:
				     /**
				      * Auxiliary vector.
				      */
    Vector<double> t;    
};


/**
 * Multi-level preconditioner.
 * Here, we collect all information needed for multi-level preconditioning
 * and provide the standard interface for LAC iterative methods.
 *
 * The template parameter class @p{MG} is required to inherit @p{MGBase}.
 * Furthermore, it needs functions @p{void copy_to_mg(const VECTOR&)}
 * to store @p{src} in the right hand side of the multi-level method and
 * @p{void copy_from_mg(VECTOR&)} to store the result of the v-cycle in @p{dst}.
 *
 * @author Guido Kanschat, 1999
 */
template<class MG, class VECTOR = Vector<double> >
class PreconditionMG
{
  public:
				     /**
				      * Constructor.
				      * Arguments are the multigrid object,
				      * pre-smoother, post-smoother and
				      * coarse grid solver.
				      */
    PreconditionMG(MG                       &mg,
		   const MGSmootherBase     &pre,
		   const MGSmootherBase     &post,
		   const MGCoarseGridSolver &coarse);
    
				     /**
				      * Preconditioning operator,
				      * calling the @p{vcycle} function
				      * of the @p{MG} object passed to
				      * the constructor.
				      *
				      * This is the operator used by
				      * LAC iterative solvers.
				      */
    void vmult (VECTOR       &dst,
		const VECTOR &src) const;
    
  private:
				     /**
				      * The multigrid object.
				      */
    SmartPointer<MG> multigrid;
    
				     /**
				      * The pre-smoothing object.
				      */
    SmartPointer<const MGSmootherBase> pre;

				     /**
				      * The post-smoothing object.
				      */
    SmartPointer<const MGSmootherBase> post;
    
				     /**
				      * The coarse grid solver.
				      */
    SmartPointer<const MGCoarseGridSolver> coarse;
};


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
  typename vector<Object>::iterator v;
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


template<class SOLVER, class MATRIX, class PRECOND>
MGCoarseGridLACIteration<SOLVER,MATRIX,PRECOND>::MGCoarseGridLACIteration(SOLVER        &s,
			 const MATRIX  &m,
			 const PRECOND &p)
		:
		solver(s),
		matrix(&m),
		precondition(p)
{};


template<class SOLVER, class MATRIX, class PRECOND>
void
MGCoarseGridLACIteration<SOLVER,MATRIX,PRECOND>::operator() (const unsigned int    /* level */,
	    Vector<double>       &dst,
	    const Vector<double> &src) const
{
  solver.solve(*matrix, dst, src, precondition);
}


/* ------------------------------------------------------------------- */


template<class MG, class VECTOR>
PreconditionMG<MG, VECTOR>::PreconditionMG(MG                       &mg,
					   const MGSmootherBase     &pre,
					   const MGSmootherBase     &post,
					   const MGCoarseGridSolver &coarse)
		:
		multigrid(&mg),
		pre(&pre),
		post(&post),
		coarse(&coarse)
{}


template<class MG, class VECTOR>
void
PreconditionMG<MG,VECTOR>::vmult (VECTOR       &dst,
				  const VECTOR &src) const
{
  multigrid->copy_to_mg(src);
  multigrid->vcycle(*pre, *post, *coarse);
  multigrid->copy_from_mg(dst);
}


/*----------------------------   mgbase.h     ---------------------------*/

#endif
/*----------------------------   mgbase.h     ---------------------------*/
