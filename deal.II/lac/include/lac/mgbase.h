/*----------------------------   mgbase.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mgbase_H
#define __mgbase_H
/*----------------------------   mgbase.h     ---------------------------*/



#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/forward-declarations.h>
#include <lac/vector.h>

#include <vector>
#include <base/smartpointer.h>

/**
 * Abstract base class for coarse grid solvers.
 * The interface of a function call operator is defined to execute coarse
 * grid solution in a derived class.
 */
class MGCoarseGridSolver
  :
  public Subscriptor
{
  public:
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGCoarseGridSolver();
    
				     /**
				      * Coarse grid solving method.
				      * This is only the interface for
				      * a function defined in a derived
				      * class like
				      * #MGCoarseGridLACIteration#.
				      *
				      * Remark that the information
				      * about the matrix is removed to
				      * that class.
				      */
    virtual void operator() (unsigned int level, Vector<double>& dst,
			     const Vector<double>& src) const = 0;
};



/**
 * Abstract base class for multigrid smoothers.
 * In fact, this class only provides the interface of the smoothing function.
 * Using deal.II grid handling, #MGSmoother# is a good point to start.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGSmootherBase
  :
  public Subscriptor 
{  
  public:
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGSmootherBase();
    
				     /**
				      * Smoothen the residual of #u# on the given
				      * level. This function should keep the interior
				      * level boundary values, so you may want
				      * to call #set_zero_interior_boundary#
				      * in the derived class #MGSmoother#
				      * somewhere in your derived function,
				      * or another function doing similar
				      * things.
				      */
    virtual void smooth (const unsigned int   level,
			 Vector<double>       &u,
			 const Vector<double> &rhs) const = 0;

};

/**
 * Smoother doing nothing.
 * @author Guido Kanschat, 1999
 */
class MGSmootherIdentity
  :
  public MGSmootherBase
{
  public:
				     /**
				      * Implementation of the interface in #MGSmootherBase#.
				      * This function does nothing.
				      */
    virtual void smooth (const unsigned int   level,
			 Vector<double>       &u,
			 const Vector<double> &rhs) const;
};

/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is an abstract one and has no implementations of
 * possible algorithms for these operations.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGTransferBase
  :
  public Subscriptor
{
  public:
				     /**
				      * Destructor. Does nothing here, but
				      * needs to be declared virtual anyway.
				      */
    virtual ~MGTransferBase();

				     /**
				      * Prolongate a vector from level
				      * #to_level-1# to level #to_level#.
				      *
				      * #src# is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the coarser level of
				      * the two involved levels, while #src#
				      * shall have as many elements as there
				      * are degrees of freedom on the finer
				      * level.
				      */
    virtual void prolongate (const unsigned int   to_level,
			     Vector<double>       &dst,
			     const Vector<double> &src) const = 0;

				     /**
				      * Restrict a vector from level
				      * #from_level# to level
				      * #from_level-1#.
				      *
				      * #src# is assumed to be a vector with
				      * as many elements as there are degrees
				      * of freedom on the finer level of
				      * the two involved levels, while #src#
				      * shall have as many elements as there
				      * are degrees of freedom on the coarser
				      * level.
				      */
    virtual void restrict (const unsigned int   from_level,
			   Vector<double>       &dst,
			   const Vector<double> &src) const = 0;
};

/**
 * An array with a vector for each level.
 * The purpose of this class is mostly to provide range checking that is missing in
 * those ultra-fashionable STL vectors.
 * @author Guido Kanschat, 1999
 */
template<class VECTOR>
class MGVector
  :
  public Subscriptor,
  public vector<VECTOR>
{
  public:
				     /**
				      * Constructor allowing to initialize the number of levels.
				      */
    MGVector(unsigned int minlevel, unsigned int maxlevel);
				     /**
				      * Safe access operator.
				      */
    VECTOR& operator[](unsigned int);
				     /**
				      * Safe access operator.
				      */
    const VECTOR& operator[](unsigned int) const;
				     /**
				      * Delete all entries.
				      */
    void clear();
  private:
				     /**
				      * Level of first component.
				      */
    unsigned int minlevel;
};



/**
 * An array of matrices for each level.
 * This class provides the functionality of #vector<MATRIX># combined
 * with a subscriptor for smart pointers.
 * @author Guido Kanschat, 1999
 */
template<class MATRIX>
class MGMatrix
  :
  public Subscriptor,
  public vector<MATRIX>
{
  public:
				     /**
				      * Constructor allowing to initialize the number of levels.
				      */
    MGMatrix(unsigned int minlevel, unsigned int maxlevel);
				     /**
				      * Safe access operator.
				      */
    MATRIX& operator[](unsigned int);
				     /**
				      * Safe access operator.
				      */
    const MATRIX& operator[](unsigned int) const;
    
  private:
				     /**
				      * Level of first component.
				      */
    unsigned int minlevel;
};


/**
 * Coarse grid solver using LAC iterative methods.
 * This is a little wrapper, transforming a triplet of iterative
 * solver, matrix and preconditioner into a coarse grid solver.
 * @author Guido Kanschat, 1999
 */
template<class SOLVER, class MATRIX, class PRECOND>
class MGCoarseGridLACIteration
  :
  public MGCoarseGridSolver
{
  public:
				     /**
				      * Constructor.
				      * Store solver, matrix and
				      * preconditioning method for later
				      * use.
				      */
    MGCoarseGridLACIteration( SOLVER&, const MATRIX&, const PRECOND&);
    
				     /**
				      * Implementation of the abstract
				      * function.
				      * Calls the solver method with
				      * matrix, vectors and
				      * preconditioner.
				      */
    virtual void operator() (unsigned int level, Vector<double>& dst,
			     const Vector<double>& src) const;
  private:
				     /**
				      * Reference to the solver.
				      */
    SOLVER& solver;
				     /**
				      * Reference to the matrix.
				      */
    const MATRIX& matrix;
    
				     /**
				      * Reference to the preconditioner.
				      */
    const PRECOND& precondition;
};

/**
 * Basic class for multigrid preconditioning.
 *
 * The functionality of the multigrid method is restricted to defect
 * correction. It is <B>not</B> iterative and the start solution is
 * always zero. Since by this <I>u<SUP>E</SUP><SUB>l</SUB></I> and
 * <I>u<SUP>A</SUP><SUB>l</SUB></I> (see report on multigrid) are
 * always zero, restriction is simplified a lot and maybe even the
 * seam condition on grids is oblivious. Still, I am not sure that
 * these restrictions on the class might cause numerical
 * inefficiencies.
 *
 * The function #precondition# is the actual multigrid method and
 * makes use of several operations to be implemented in derived
 * classes. It takes a defect in #src# and the result is the multigrid
 * preconditioned defect in #dst#.
 *
 * @author Guido Kanschat, 1999
 * */
class MGBase
  :
  public Subscriptor
{
  private:
    MGBase(const MGBase&);
    const MGBase& operator=(const MGBase&);

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
				      * Auxiliary vector, defect.
				      */
    MGVector<Vector<double> > d;

				     /**
				      * Auxiliary vector, solution.
				      */
    MGVector<Vector<double> > s;
				     /**
				      * Prolongation and restriction object.
				      */
    SmartPointer<const MGTransferBase> transfer;

  private:
				     /**
				      * Auxiliary vector.
				      */
    Vector<double> t;
    
				     /**
				      * Exception.
				      */
    DeclException2(ExcSwitchedLevels, int, int,
		   << "minlevel and maxlevel switched, should be: "
		   << arg1 << "<=" << arg2);
    
  protected:

				     /**
				      * The actual v-cycle multigrid method.
				      * This function is called on the
				      * highest level and recursively
				      * invokes itself down to the
				      * coarsest. There, it calls
				      * #coarse_grid_solver# and
				      * proceeds back up.
				      */
    void level_mgstep(unsigned int level,
		      const MGSmootherBase& pre_smooth,
		      const MGSmootherBase& post_smooth,
		      const MGCoarseGridSolver& cgs);  

				     /**
				      * Apply negative #vmult# on all
				      * cells of a level.
				      * This is implemented in a
				      * derived class.
				      */
    virtual void level_vmult(unsigned int level,
			     Vector<double>& dst,
			     const Vector<double>& src,
			     const Vector<double>& rhs) = 0;  

  
  public:
				     /**
				      * Constructor, subject to change.
				      */
    MGBase(const MGTransferBase& transfer,
		  unsigned int minlevel, unsigned int maxlevel);
				     /**
				      * Virtual destructor.
				      */
    virtual ~MGBase();
				     /**
				      * Execute v-cycle algorithm.
				      * This function assumes, that
				      * the vector #d# is properly
				      * filled with the residual in
				      * the outer defect correction
				      * scheme. After execution of
				      * #vcycle()#, the result is in
				      * the vector #s#. We propose to
				      * write functions #copy_from_mg#
				      * and #copy_to_mg# in derived
				      * classes of #MGBase# to access
				      * #d# and #s#.
				      */
    void vcycle(const MGSmootherBase& pre_smooth,
		const MGSmootherBase& post_smooth,
		const MGCoarseGridSolver& cgs);
};



/**
 * Multi-level preconditioner.
 * Here, we collect all information needed for multi-level preconditioning
 * and provide the standard interface for LAC iterative methods.
 *
 * The template parameter class #MG# is required to inherit #MGBase#.
 * Furthermore, it needs functions #void copy_to_mg(const VECTOR&)#
 * to store #src# in the right hand side of the multi-level method and
 * #void copy_from_mg(VECTOR&)# to store the result of the v-cycle in #dst#.
 * @author Guido Kanschat, 1999
 */
template<class MG, class VECTOR>
class PreconditionMG
{
  public:
				     /**
				      * Constructor.
				      * Arguments are the multigrid object,
				      * pre-smoother, post-smoother and
				      * coarse grid solver.
				      */
    PreconditionMG(MG& mg,
		   const MGSmootherBase& pre,
		   const MGSmootherBase& post,
		   const MGCoarseGridSolver& coarse);
				     /**
				      * Preconditioning operator.
				      * This is the operator used by LAC
				      * iterative solvers.
				      */
    void operator() (VECTOR& dst, const VECTOR& src) const;
  private:
				     /**
				      * The mulrigrid object.
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





template<class SOLVER, class MATRIX, class PRECOND>
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND>::
MGCoarseGridLACIteration(SOLVER& s, const MATRIX& m, const PRECOND& p)
		:
		solver(s), matrix(m), precondition(p)
{}

template<class SOLVER, class MATRIX, class PRECOND>
void
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND>::operator()
(unsigned int, Vector<double>& dst, const Vector<double>& src) const
{
  solver.solve(matrix, dst, src, precondition);
}




/* ------------------------------------------------------------------- */

template<class VECTOR>
MGVector<VECTOR>::MGVector(unsigned int min, unsigned int max)
		:
		vector<VECTOR>(max-min+1),
		minlevel(min)
{}


template<class VECTOR>
VECTOR&
MGVector<VECTOR>::operator[](unsigned int i)
{
  Assert((i>=minlevel)&&(i<minlevel+size()),ExcIndexRange(i,minlevel,minlevel+size()));
  return vector<VECTOR>::operator[](i-minlevel);
}


template<class VECTOR>
const VECTOR&
MGVector<VECTOR>::operator[](unsigned int i) const
{
  Assert((i>=minlevel)&&(i<minlevel+size()),ExcIndexRange(i,minlevel,minlevel+size()));
  return vector<VECTOR>::operator[](i-minlevel);
}

template<class VECTOR>
void
MGVector<VECTOR>::clear()
{
  typename vector<VECTOR>::iterator v;
  for (v = begin(); v != end(); ++v)
    v->clear();
  
}


/* ------------------------------------------------------------------- */


template<class MATRIX>
MGMatrix<MATRIX>::MGMatrix(unsigned int min, unsigned int max)
		:
		vector<MATRIX>(max-min+1),
		minlevel(min)
{}


template<class MATRIX>
MATRIX&
MGMatrix<MATRIX>::operator[](unsigned int i)
{
  Assert((i>=minlevel)&&(i<minlevel+size()),ExcIndexRange(i,minlevel,minlevel+size()));
  return vector<MATRIX>::operator[](i-minlevel);
}


template<class MATRIX>
const MATRIX&
MGMatrix<MATRIX>::operator[](unsigned int i) const
{
  Assert((i>=minlevel)&&(i<minlevel+size()),ExcIndexRange(i,minlevel,minlevel+size()));
  return vector<MATRIX>::operator[](i-minlevel);
}



/* ------------------------------------------------------------------- */


template<class MG, class VECTOR>
PreconditionMG<MG, VECTOR>::PreconditionMG(MG& mg,
					   const MGSmootherBase& pre,
					   const MGSmootherBase& post,
					   const MGCoarseGridSolver& coarse)
		:
		multigrid(&mg),
		pre(&pre),
		post(&post),
		coarse(&coarse)
{}



template<class MG, class VECTOR>
void
PreconditionMG<MG, VECTOR>::operator() (VECTOR& dst,
					const VECTOR& src) const
{
  multigrid->copy_to_mg(src);
  multigrid->vcycle(*pre, *post, *coarse);
  multigrid->copy_from_mg(dst);
}



/*----------------------------   mgbase.h     ---------------------------*/
/* end of #ifndef __mgbase_H */
#endif
/*----------------------------   mgbase.h     ---------------------------*/
