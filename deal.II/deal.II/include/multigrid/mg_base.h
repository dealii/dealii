// $Id$


#ifndef __lac_mg_base_h
#define __lac_mg_base_h

#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/forward-declarations.h>
#include <lac/vector.h>

#include <vector>


/**
 * Abstract base class for coarse grid solvers.
 * The interface of a function call operator is defined to execute coarse
 * grid solution in a derived class.
 */
class MGCoarseGridSolver
{
  public:
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
    virtual void operator() (unsigned int level, Vector<float>& dst,
			     const Vector<float>& src) const = 0;
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
			 Vector<float>       &u,
			 const Vector<float> &rhs) const = 0;

};

/**
 * Base class used to declare the operations needed by a concrete class
 * implementing prolongation and restriction of vectors in the multigrid
 * context. This class is an abstract one and has no implementations of
 * possible algorithms for these operations.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGTransferBase  :  public Subscriptor
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
			     Vector<float>       &dst,
			     const Vector<float> &src) const = 0;

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
			   Vector<float>       &dst,
			   const Vector<float> &src) const = 0;
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
  public Vector<MATRIX>
{
  public:
				     /**
				      * Constructor allowing to initialize the number of levels.
				      */
    MGMatrix(unsigned int n_levels);
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
    MGCoarseGridLACIteration(const SOLVER&, const MATRIX&, const PRECOND&);
    
				     /**
				      * Implementation of the abstract
				      * function.
				      * Calls the solver method with
				      * matrix, vectors and
				      * preconditioner.
				      */
    virtual void operator() (unsigned int level, Vector<float>& dst,
			     const Vector<float>& src) const;
  private:
				     /**
				      * Reference to the solver.
				      */
    const SOLVER& solver;
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
{
    MGBase(const MGBase&);
    const MGBase& operator=(const MGBase&);
    
				     /**
				      * Auxiliary vector, defect.
				      */
    vector<Vector<float> > d;

				     /**
				      * Auxiliary vector, solution.
				      */
    vector<Vector<float> > s;

				     /**
				      * Auxiliary vector.
				      */
    Vector<float> t;
    
				     /**
				      * Prolongation and restriction object.
				      */
    SmartPointer<const MGTransferBase> transfer;
    
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
				      * Apply residual operator on all
				      * cells of a level.
				      * This is implemented in a
				      * derived class.
				      */
    virtual void level_residual(unsigned int level,
			     Vector<float>& dst,
			     const Vector<float>& src,
			     const Vector<float>& rhs) = 0;  

  
  public:
				     /**
				      * Constructor, subject to change.
				      */
    MGBase(const MGTransferBase& transfer,
		  unsigned int maxlevel, unsigned int minlevel);
    virtual ~MGBase();
    
};


template<class SOLVER, class MATRIX, class PRECOND>
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND>
::MGCoarseGridLACIteration(const SOLVER& s, const MATRIX& m, const PRECOND& p)
		:
		solver(s), matrix(m), precondition(p)
{}

template<class SOLVER, class MATRIX, class PRECOND>
void
MGCoarseGridLACIteration<SOLVER, MATRIX, PRECOND>::operator()
(unsigned int level, Vector<float>& dst, const Vector<float>& src) const
{
  solver.solve(matrix, dst, src, precondition);
}

template<class MATRIX>
MGMatrix<MATRIX>::MGMatrix(unsigned int nlevels)
		:
		vector<MATRIX>(nlevels)
{}


#endif
