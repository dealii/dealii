/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

#include <base/smartpointer.h>
#include <lac/forward-declarations.h>
#include <basic/forward-declarations.h>
#include <lac/sparsematrix.h>
#include <lac/vector.h>

#include <vector>
class MGTransferBase;
class MGSmoother;


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
    
				     /**
				      * Auxiliary vector.
				      */
    vector<Vector<float> > d;

				     /**
				      * Auxiliary vector.
				      */
    vector<Vector<float> > s;

				     /**
				      * Auxiliary vector.
				      */
    Vector<float> t;
    
				     /**
				      * Highest level of cells.
				      */
    unsigned maxlevel;

				     /**
				      * Level for coarse grid solution.
				      */
    unsigned minlevel;

				     /**
				      * Prolongation and restriction object.
				      */
    SmartPointer<MGTransferBase> transfer;
    
  protected:
				     /**
				      * Number of pre-smoothing steps. Is used
				      * as a parameter to #pre_smooth#.
				      */
    unsigned n_pre_smooth;
  
				     /**
				      * Number of post-smoothing steps Is used
				      * as a parameter to #post_smooth#.
				      */
    unsigned n_post_smooth;
				     /**
				      * The actual v-cycle multigrid method.
				      * This function is called on the
				      * highest level and recursively
				      * invokes itself down to the
				      * coarsest. There, it calls
				      * #coarse_grid_solution# and
				      * proceeds back up.
				      */
    void level_mgstep(unsigned level, MGSmoother& smoother);  

				     /**
				      * Apply residual operator on all
				      * cells of a level.
				      * This is implemented in a
				      * derived class.
				      */
    virtual void level_residual(unsigned level,
			     Vector<float>& dst,
			     const Vector<float>& src,
			     const Vector<float>& rhs) = 0;

				     /**
				      * Solve exactly on coarsest grid.
				      * Usually, this function should
				      * be overloaded by a more
				      * sophisticated derived
				      * class. Still, there is a
				      * standard implementation doing
				      * #10 * (n_pre_smooth +
				      * n_post_smooth)#
				      * smoothing steps of #pre_smooth#.
				      */
    virtual void coarse_grid_solution(unsigned l,
				      Vector<float>& dst,
				      const Vector<float>& src);
  

  

  
  public:
				     /**
				      * Constructor, subject to change.
				      */
    MultiGridBase(MGTransferBase& transfer,
		  unsigned maxlevel, unsigned minlevel,
		  unsigned pre_smooth, unsigned post_smooth);
    virtual ~MultiGridBase();
    
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
 * Implementation of multigrid with matrices.
 * While #MultiGridBase# was only an abstract framework for the v-cycle,
 * here we have the implementation of the pure virtual functions defined there.
 * Furthermore, level information is obtained from a triangulation object.
 * 
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
template<int dim, class Matrix>
class MultiGrid
{
  public:
    MultiGrid(const MGDoFHandler<dim>&);
    
				     /** Transfer from dVector to
				      * MGVector.
				      *
				      * This function copies data from a
				      * dVector, that is, data on the
				      * locally finest level, into the
				      * corresponding levels of an
				      * MGVector.
				      */
    template<typename number>
    void copy_to_mg(vector<Vector<float> >& dst,
		    const Vector<number>& src) const;

				     /**
				      * Transfer from multi-grid vector to
				      * normal vector.
				      *
				      * Copies data from active portions
				      * of an MGVector into the
				      * respective positions of a
				      * Vector<double>. All other entries of
				      * #src# are zero.
				      */
    template<typename number>
    void copy_from_mg(Vector<number>& dst,
		      const vector<Vector<float> >& src) const;
  private:
				     /**
				      * Associated #MGDoFHandler#.
				      */
    SmartPointer<MGDoFHandler<dim> > dofs;

				     /**
				      * Matrix structures for each level.
				      * These are generated by the constructor
				      * of #MultiGrid# and can be accessed
				      * later.
				      */
    vector<SparseMatrixStruct> matrix_structures;

				     /**
				      * Matrices for each level.
				      * The matrices are prepared by
				      * the constructor of #MultiGrid# and can
				      * be accessed for assembling.
				      */
    vector<SparseMatrix<float> > matrices;
};


/**
 * Implementation of the #MGTransferBase# interface for which the transfer
 * operations are prebuilt upon construction of the object of this class as
 * matrices. This is the fast way, since it only needs to build the operation
 * once by looping over all cells and storing the result in a matrix for
 * each level, but requires additional memory.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGTransferPrebuilt : public MGTransferBase 
{
  public:

				     /**
				      * Destructor.
				      */
    virtual ~MGTransferPrebuilt ();
    
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof);

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
			     const Vector<float> &src) const;

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
			   const Vector<float> &src) const;

  private:

    vector<SparseMatrixStruct>   prolongation_sparsities;

					 /**
					  * The actual prolongation matrix.
					  * column indices belong to the
					  * dof indices of the mother cell,
					  * i.e. the coarse level.
					  * while row indices belong to the
					  * child cell, i.e. the fine level.
					  */
    vector<SparseMatrix<float> > prolongation_matrices;
};



/*----------------------------   multigrid.h     ---------------------------*/
/* end of #ifndef __multigrid_H */
#endif
/*----------------------------   multigrid.h     ---------------------------*/
