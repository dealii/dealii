/*----------------------------   multigrid.h     ---------------------------*/
/*      $Id$                 */
#ifndef __multigrid_H
#define __multigrid_H
/*----------------------------   multigrid.h     ---------------------------*/

#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/forward-declarations.h>
#include <basic/forward-declarations.h>
#include <lac/sparsematrix.h>
#include <lac/vector.h>
#include <lac/mgbase.h>

#include <vector>

/**
 * Implementation of multigrid with matrices.
 * While #MGBase# was only an abstract framework for the v-cycle,
 * here we have the implementation of the pure virtual functions defined there.
 * Furthermore, level information is obtained from a triangulation object.
 * 
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
template<int dim>
class MG
  :
  public MGBase
{
  public:
    MG(const MGDoFHandler<dim>&,
       ConstraintMatrix& hanging_nodes,
       const MGMatrix<SparseMatrix<double> >&,
       const MGTransferBase& transfer,
       unsigned int minlevel = 0, unsigned int maxlevel = 10000);
    
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
    void copy_to_mg(const Vector<number>& src);

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
    void copy_from_mg(Vector<number>& dst) const;

				     /**
				      * Negative #vmult# on a level.
				      * @see MGBase.
				      */
    virtual void level_vmult(unsigned int,
				Vector<double> &,
				const Vector<double> &,
				const Vector<double> &);

				     /**
				      * Read-only access to level matrices.
				      */
    const SparseMatrix<double>& get_matrix(unsigned int level) const;

  private:
				     /**
				      * Associated #MGDoFHandler#.
				      */
    SmartPointer<const MGDoFHandler<dim> > dofs;

				     /**
				      * Matrices for each level.
				      * The matrices are prepared by
				      * the constructor of #MG# and can
				      * be accessed for assembling.
				      */
    SmartPointer<const MGMatrix<SparseMatrix<double> > > matrices;
    ConstraintMatrix& hanging_nodes;
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
				      * Preliminary constructor
				      */
    MGTransferPrebuilt()
      {}
    
    
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
			     Vector<double>       &dst,
			     const Vector<double> &src) const;

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
			   const Vector<double> &src) const;

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


template<int dim>
inline
const SparseMatrix<double>&
MG<dim>::get_matrix(unsigned int level) const
{
  Assert((level>=minlevel) && (level<maxlevel), ExcIndexRange(level, minlevel, maxlevel));
  
  return (*matrices)[level];
}


/*----------------------------   multigrid.h     ---------------------------*/
/* end of #ifndef __multigrid_H */
#endif
/*----------------------------   multigrid.h     ---------------------------*/
