/*----------------------------   mg_smoother.h     ---------------------------*/
/*      $Id$                 */
#ifndef __mg_smoother_H
#define __mg_smoother_H
/*----------------------------   mg_smoother.h     ---------------------------*/


#include <lac/forward-declarations.h>
#include <lac/mgbase.h>
#include <basic/forward-declarations.h>
#include <base/smartpointer.h>
#include <vector>


/**
 * Base class for multigrid smoothers. It declares the interface
 * to smoothers by inheriting #MGSmootherBase# and implements some functionality for setting
 * the values
 * of vectors at interior boundaries (i.e. boundaries between differing
 * levels of the triangulation) to zero, by building a list of these degrees
 * of freedom's indices at construction time.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGSmoother
  :  public MGSmootherBase
{
  private:
				     /**
				      * Default constructor. Made private
				      * to prevent it being called, which
				      * is necessary since this could
				      * circumpass the set-up of the list
				      * if interior boundary dofs.
				      */
    MGSmoother ();
    
  public:

				     /**
				      * Constructor. This one collects
				      * the indices of the degrees of freedom
				      * on the interior boundaries between
				      * the different levels, which are
				      * needed by the function
				      * #set_zero_interior_boundaries#.
				      *
				      * Since this function is implemented
				      * a bit different in 1d (there are no
				      * faces of cells, just vertices),
				      * there are actually two sets of
				      * constructors, namely this one for 1d
				      * and the following one for all other
				      * dimensions.
				      */
    MGSmoother (const MGDoFHandler<1> &mg_dof);

				     /**
				      * Constructor. This one collects
				      * the indices of the degrees of freedom
				      * on the interior boundaries between
				      * the different levels, which are
				      * needed by the function
				      * #set_zero_interior_boundaries#.
				      *
				      * The parameter steps indicates the number of smoothing
				      * steps to be executed by #smooth#.
				      */
    template <int dim>
    MGSmoother (const MGDoFHandler<dim> &mg_dof, unsigned int steps);    

				     /**
				      * Reset the values of the degrees of
				      * freedom on interior boundaries between
				      * different levels to zero in the given
				      * data vector #u#.
				      *
				      * Since the coarsest level (#level==0#)
				      * has no interior boundaries, this
				      * function does nothing in this case.
				      */
    void set_zero_interior_boundary (const unsigned int  level,
				     Vector<float>      &u) const;
				     /**
				      * Modify the number of smoothing steps.
				      */
    void set_steps(unsigned int steps);
				     /**
				      * How many steps should be used?
				      */
    unsigned int get_steps() const;

  private:
				     /**
				      * Number of smoothing steps.
				      */
    unsigned int steps;
				     /**
				      * For each level, we store a list of
				      * degree of freedom indices which are
				      * located on interior boundaries between
				      * differing levels of the triangulation.
				      * Since the coarsest level has no
				      * interior boundary dofs, the first
				      * entry refers to the second level.
				      *
				      * These arrays are set by the constructor.
				      * The entries for each level are sorted ascendingly.
				      */
    vector<vector<int> > interior_boundary_dofs;
};

/**
 * Smoother doing nothing.
 * @author Guido Kanschat, 1999
 */
class MGSmootherIdentity
  :
  public MGSmoother
{
  public:
				     /**
				      * Constructor.
				      */
    template<int dim>
    MGSmootherIdentity(const MGDoFHandler<dim> &mg_dof);
				     /**
				      * Implementation of the interface in #MGSmootherBase#.
				      * This function does nothing.
				      */
    virtual void smooth (const unsigned int   level,
			 Vector<float>       &u,
			 const Vector<float> &rhs) const;
};

/**
 * Implementation of an SOR smoother.
 *
 * This class is an implementation of a smootin algorithm for multigrid. It
 * knows the actual instance of the multigrid matrix and uses the SOR method of the level
 * matrices.
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGSmootherSOR
  :
  public MGSmoother
{
  public:
				     /**
				      * Constructor.
				      * The constructor uses an #MGDoFHandler# to initialize
				      * data structures for interior level boundary handling
				      * in #MGSmoother#. Furthermore, it takes a pointer to the
				      * matrices and the number of SOR steps required in each
				      * smoothing operation.
				      */

    template<int dim>
    MGSmootherSOR(const MGDoFHandler<dim> &mg_dof,
		  const MGMatrix<SparseMatrix<float> >& matrix,
		  unsigned int steps);
    
				     /**
				      * Implementation of the interface in #MGSmootherBase#.
				      * We use the SOR method in #SparseMatrix# for the real work
				      * and find a way to keep the boundary values.
				      */
    virtual void smooth (const unsigned int   level,
			 Vector<float>       &u,
			 const Vector<float> &rhs) const;
  private:
				     /**
				      * Pointer to the matrices.
				      */
    SmartPointer< const MGMatrix< SparseMatrix<float> > > matrix;
};

inline
void
MGSmoother::set_steps(unsigned int i)
{
  steps = i;
}

inline
unsigned int
MGSmoother::get_steps() const
{
  return steps;
}

/*----------------------------   mg_smoother.h     ---------------------------*/
/* end of #ifndef __mg_smoother_H */
#endif
/*----------------------------   mg_smoother.h     ---------------------------*/
			 
