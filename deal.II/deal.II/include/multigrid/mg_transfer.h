//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_transfer_h
#define __deal2__mg_transfer_h

#include <base/config.h>

#include <lac/block_vector.h>
#include <lac/constraint_matrix.h>
#ifdef DEAL_PREFER_MATRIX_EZ
#  include <lac/sparse_matrix_ez.h>
#  include <lac/block_sparse_matrix_ez.h>
#else
#  include <lac/sparsity_pattern.h>
#  include <lac/block_sparsity_pattern.h>
#endif
#include <lac/vector_memory.h>

#include <multigrid/mg_base.h>
#include <multigrid/mg_level_object.h>



#include <dofs/dof_handler.h>

#include <base/std_cxx1x/shared_ptr.h>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim> class MGDoFHandler;

/*
 * MGTransferBase is defined in mg_base.h
 */

/*!@addtogroup mg */
/*@{*/

/**
 * Implementation of the MGTransferBase interface for which the transfer
 * operations are prebuilt upon construction of the object of this class as
 * matrices. This is the fast way, since it only needs to build the operation
 * once by looping over all cells and storing the result in a matrix for
 * each level, but requires additional memory.
 *
 * See MGTransferBase to find out which of the transfer classes
 * is best for your needs.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999-2004
 */
template <class VECTOR>
class MGTransferPrebuilt : public MGTransferBase<VECTOR> 
{
  public:
				     /**
				      * Constructor without constraint
				      * matrices. Use this constructor
				      * only with discontinuous finite
				      * elements or with no local
				      * refinement.
				      */
    MGTransferPrebuilt ();

				     /**
				      * Constructor with constraint matrices.
				      */
    MGTransferPrebuilt (const ConstraintMatrix& constraints);
				     /**
				      * Destructor.
				      */
    virtual ~MGTransferPrebuilt ();
				     /**
				      * Actually build the prolongation
				      * matrices for each level.
				      */
    template <int dim, int spacedim>
    void build_matrices (const MGDoFHandler<dim,spacedim> &mg_dof);

    virtual void prolongate (const unsigned int    to_level,
			     VECTOR       &dst,
			     const VECTOR &src) const;

    virtual void restrict_and_add (const unsigned int    from_level,
				   VECTOR       &dst,
				   const VECTOR &src) const;

    				     /**
				      * Transfer from a vector on the
				      * global grid to vectors defined
				      * on each of the levels
				      * separately, i.a. an @p MGVector.
				      */
    template <int dim, class InVector, int spacedim>
    void
    copy_to_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
		MGLevelObject<VECTOR>& dst,
		const InVector& src) const;

				     /**
				      * Transfer from multi-level vector to
				      * normal vector.
				      *
				      * Copies data from active
				      * portions of an MGVector into
				      * the respective positions of a
				      * <tt>Vector<number></tt>. In order to
				      * keep the result consistent,
				      * constrained degrees of freedom
				      * are set to zero.
				      */
    template <int dim, class OutVector, int spacedim>
    void
    copy_from_mg (const MGDoFHandler<dim,spacedim>& mg_dof,
		  OutVector& dst,
		  const MGLevelObject<VECTOR> &src) const;

				     /**
				      * Add a multi-level vector to a
				      * normal vector.
				      *
				      * Works as the previous
				      * function, but probably not for
				      * continuous elements.
				      */
    template <int dim, class OutVector, int spacedim>
    void
    copy_from_mg_add (const MGDoFHandler<dim,spacedim>& mg_dof,
		      OutVector& dst,
		      const MGLevelObject<VECTOR>& src) const;

				     /**
				      * Finite element does not
				      * provide prolongation matrices.
				      */
    DeclException0(ExcNoProlongation);
    
				     /**
				      * Call @p build_matrices
				      * function first.
				      */
    DeclException0(ExcMatricesNotBuilt);

    				     /**
				      * Memory used by this object.
				      */
    unsigned int memory_consumption () const;
    

  private:

				   /**
				    * Sizes of the multi-level vectors.
				    */
    std::vector<unsigned int> sizes;

				     /**
				      * Sparsity patterns for transfer
				      * matrices.
				      */
    std::vector<std_cxx1x::shared_ptr<SparsityPattern> >   prolongation_sparsities;

				     /**
				      * The actual prolongation matrix.
				      * column indices belong to the
				      * dof indices of the mother cell,
				      * i.e. the coarse level.
				      * while row indices belong to the
				      * child cell, i.e. the fine level.
				      */
    std::vector<std_cxx1x::shared_ptr<SparseMatrix<double> > > prolongation_matrices;    
    
				     /**
				      * Mapping for the
				      * <tt>copy_to/from_mg</tt>-functions.
				      * The data is first the global
				      * index, then the level index.
				     */
    std::vector<std::map<unsigned int, unsigned int> >
    copy_indices;

				     /**
				      * Fill the two boolean vectors
				      * #dofs_on_refinement_edge and
				      * #dofs_on_refinement_boundary.
				      */
    template <int dim, int spacedim>
    void find_dofs_on_refinement_edges (
        const MGDoFHandler<dim,spacedim>& mg_dof);

				     /**
				      * Degrees of freedom on the
				      * refinement edge excluding
				      * those on the boundary.
				      *
				      * @todo Clean up names
				      */
    std::vector<std::vector<bool> > dofs_on_refinement_edge;
				     /**
				      * The degrees of freedom on the
				      * the refinement edges. For each
				      * level (outer vector) and each
				      * dof index (inner vector), this
				      * bool is true if the level
				      * degree of freedom is on the
				      * refinement edge towards the
				      * lower level.
				      *
				      * @todo Clean up names
				      */
    std::vector<std::vector<bool> > dofs_on_refinement_boundary; 
				     /**
				      * The constraints of the global
				      * system.
				      */
    SmartPointer<const ConstraintMatrix> constraints;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
