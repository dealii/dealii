/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

#include <vector>
#include <base/subscriptor.h>
#include <lac/forward-declarations.h>

/**
 * Basic matrix class for multigrid preconditioning.
 *
 * This matrix may be used in the iterative methods of LAC, where the
 * functions #vmult# and #precondition# and possibly their transposed
 * versions are needed.
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
class MultiGridBase
{
    MultiGridBase(const MultiGridBase&);
    const MultiGridBase& operator=(const MultiGridBase&);
    
				     /** Auxiliary vectors.
				      */
    vector<Vector<float> > d;
    vector<Vector<float> > s;
    
				     /**
				      * Highest level of cells.
				      */
    unsigned maxlevel;

				     /**
				      * Level for coarse grid solution.
				      */
    unsigned minlevel;
  
				     /** Tranfer from dVector to
				      * MGVector.
				      *
				      * This function copies data from a
				      * dVector, that is, data on the
				      * locally finest level, into the
				      * corresponding levels of an
				      * MGVector.
				      */
    void copy_to_mg(vector<Vector<float> >& dst, const Vector<double>& src);

				     /**
				      * Transfer from MGVector to
				      * Vector<double>.
				      *
				      * Copies data from active portions
				      * of an MGVector into the
				      * respective positions of a
				      * Vector<double>. All other entries of
				      * #src# are zero.
				      */
    void copy_from_mg(Vector<double>& dst, const vector<Vector<float> >& src);

				     /**
				      * The actual v-cycle multigrid method.
				      * This function is called on the
				      * highest level and recursively
				      * invokes itself down to the
				      * coarsest. There, it calls
				      * #coarse_grid_solution# and
				      * proceeds back up.
				      */
    void level_mgstep(unsigned level);
  
  protected:
				     /**
				      * Number of pre-smoothing steps.
				      */
    unsigned n_pre_smooth;
  
				     /**
				      * Number of post-smoothing steps
				      */
    unsigned n_post_smooth;
				     /**
				      * The (pre-)smoothing algorithm.
				      * This function is required to
				      * perform #steps# iterations to
				      * smoothen the residual $Ax-b$.
				      *
				      * When overloading this function
				      * in derived classes, make sure
				      * that smoothing only applies to
				      * interior degrees of freedom of
				      * the actual level.
				      *
				      */
    virtual void smooth(unsigned level,
			Vector<float>& x,
			const Vector<float>& b,
			unsigned steps);

				     /**
				      * The post-smoothing algorithm.
				      * Defaults to #smooth#.
				      */
    virtual void post_smooth(unsigned level,
			     Vector<float>& dst,
			     const Vector<float>& src,
			     unsigned steps);

				     /**
				      * Apply operator on all
				      * cells of a level.
				      *
				      */
    virtual void level_vmult(unsigned level,
			     Vector<float>& dst,
			     const Vector<float>& src);

				     /**
				      * Solve exactly on coarsest grid.
				      */
    virtual void coarse_grid_solution(unsigned l,
				      Vector<float>& dst,
				      const Vector<float>& src);
  

  

  
  public:
				     /**
				      * Constructor, subject to change.
				      */
    MultiGridBase();
    virtual ~MultiGridBase();
    
};


class MGTransferBase
  :
  public Subscriptor
{
  public:
    virtual ~MGTransferBase();

    virtual void prolongate(unsigned l,
			    Vector<float>& dst,
			    const Vector<float>& src) const = 0;
    virtual void restrict(unsigned l,
			  Vector<float>& dst,
			  const Vector<float>& src) const = 0;
};


/*----------------------------   multigrid.h     ---------------------------*/
/* end of #ifndef __multigrid_H */
#endif
/*----------------------------   multigrid.h     ---------------------------*/
