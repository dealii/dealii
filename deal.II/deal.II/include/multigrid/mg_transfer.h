//----------------------------------------------------------------------
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
//----------------------------------------------------------------------
#ifndef __deal2__mg_transfer_h
#define __deal2__mg_transfer_h


#include <lac/sparsity_pattern.h>
#include <lac/block_vector.h>
#include <lac/block_sparsity_pattern.h>
#include <multigrid/mg_base.h>

template <int dim> class MGDoFHandler;


/**
 * Implementation of the @p{MGTransferBase} interface for which the transfer
 * operations are prebuilt upon construction of the object of this class as
 * matrices. This is the fast way, since it only needs to build the operation
 * once by looping over all cells and storing the result in a matrix for
 * each level, but requires additional memory.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999
 */
class MGTransferPrebuilt : public Subscriptor // MGTransferBase<Vector<double> > 
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
			     const Vector<double> &src) const;

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
				   const Vector<double> &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template<int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>& mg_dof,
		MGLevelObject<Vector<double> > &dst,
		const InVector &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<double>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template<int dim, class OutVector>
    void
    copy_from_mg (const MGDoFHandler<dim>& mg_dof,
		  OutVector &dst,
		  const MGLevelObject<Vector<double> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template<int dim, class OutVector>
    void
    copy_from_mg_add (const MGDoFHandler<dim>& mg_dof,
		      OutVector &dst,
		      const MGLevelObject<Vector<double> > &src) const;

  private:

				   /**
				    * Selected component.
				    */
    unsigned int selected;

				   /**
				    * Sizes of the multi-level vectors.
				    */
    std::vector<unsigned int> sizes;


    std::vector<SparsityPattern>   prolongation_sparsities;

				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<SparseMatrix<double> > prolongation_matrices;
};


/**
 * Implementation of matrix generation for @ref{MGTransferBlock} and @p{MGTransferSelected}.
 *
 * @author Guido Kanschat, 2001
 */
class MGTransferBlockBase
{
  protected:  
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      *
				      * If a field selected is given,
				      * only matrices for these
				      * components are to be built. By
				      * default, all matrices are
				      * built.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof,
			 std::vector<bool> selected);

				   /**
				    * Flag of selected components.
				    */
    std::vector<bool> selected;

				   /**
				    * Sizes of the multi-level vectors.
				    */
    std::vector<std::vector<unsigned int> > sizes;
  
				   /**
				    * Start index of each component.
				    */
    std::vector<unsigned int> component_start;
  
				   /**
				    * Start index of eah component on
				    * all levels.
				    */
    std::vector<std::vector<unsigned int> > mg_component_start;

  private:

    std::vector<BlockSparsityPattern>   prolongation_sparsities;

  protected:
  
				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<BlockSparseMatrix<double> > prolongation_matrices;
};



/**
 * Implementation of the @p{MGTransferBase} interface for block
 * matrices and block vectors.  In addition to the functionality of
 * @ref{MGTransferPrebuilt}, the operation may be restricted to
 * certain blocks of the vector.
 *
 * If the restricted mode is chosen, block vectors used in the
 * transfer routines may only have as many components as there are
 * @p{true}s in the selected-field.
 *
 * @author Guido Kanschat, 2001
 */
class MGTransferBlock : public Subscriptor, // MGTransferBase<BlockVector<double> >,
			private MGTransferBlockBase
{
  public:
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferBlock ();
    
				     /**
				      * Build the prolongation
				      * matrices for each level.
				      *
				      * If a field selected is given,
				      * only matrices for these
				      * components are to be built. By
				      * default, all matrices are
				      * built.
				      *
				      * This function is a front-end
				      * for the same function in
				      * @ref{MGTransferBlockBase}.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof,
			 std::vector<bool> selected = std::vector<bool>());

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
			     BlockVector<double>       &dst,
			     const BlockVector<double> &src) const;

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
				   BlockVector<double>       &dst,
				   const BlockVector<double> &src) const;

  protected: 
};



/**
 * Implementation of the @p{MGTransferBase} interface for block
 * matrices and simple vectors. This class uses @ref{MGTransferBlock}
 * selecting a single component. The transfer operators themselves are
 * implemented for simple vectors again.
 *
 * @author Guido Kanschat, 2001
 */
class MGTransferSelect : public Subscriptor, // MGTransferBase<Vector<double> >,
			 private MGTransferBlockBase
{
  public:
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferSelect ();
    
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      *
				      * Select the component you want
				      * to apply multigrid to.
				      *
				      * This function is a front-end
				      * for the same function in
				      * @ref{MGTransferBlockBase}.
				      */
    template <int dim>
    void build_matrices (const MGDoFHandler<dim> &mg_dof,
			 unsigned int selected);

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
			     const Vector<double> &src) const;

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
				   const Vector<double> &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p{MGVector}.
				      */
    template<int dim, class InVector>
    void
    copy_to_mg (const MGDoFHandler<dim>& mg_dof,
		MGLevelObject<Vector<double> > &dst,
		const InVector &src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * @p{Vector<double>}. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template<int dim, class OutVector>
    void
    copy_from_mg (const MGDoFHandler<dim>& mg_dof,
		  OutVector &dst,
		  const MGLevelObject<Vector<double> > &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template<int dim, class OutVector>
    void
    copy_from_mg_add (const MGDoFHandler<dim>& mg_dof,
		      OutVector &dst,
		      const MGLevelObject<Vector<double> > &src) const;

  private:

				   /**
				    * Selected component.
				    */
  unsigned int selected;
};



#endif
