//----------------------------  mg_dof_tools.h  ---------------------------
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
//----------------------------  mg_dof_tools.h  ---------------------------
#ifndef __deal2__mg_dof_tools_h
#define __deal2__mg_dof_tools_h

template <class Object> class MGLevelObject;
template <int dim> class MGDoFHandler;
template <typename number> class Vector;
template <typename number> class BlockVector;
template <class number> class FullMatrix;

#include <base/config.h>


/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom in a multilevel triangulation. It
 * is similar in purpose and function to the @p{DoFTools} class, but
 * operates on @p{MGDoFHandler} objects instead of @ref{DoFHandler}
 * objects. See there and the documentation of the member functions
 * for more information.
 *
 * All member functions are static, so there is no need to create an
 * object of class @p{MGTools}.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999, 2000
 */
class MGTools
{
  public:
				     /**
				      * Write the sparsity structure
				      * of the matrix belonging to the
				      * specified @p{level}. The sparsity pattern
				      * is not compressed, so before 
				      * creating the actual matrix
				      * you have to compress the
				      * matrix yourself, using
				      * @p{SparseMatrixStruct::compress()}.
				      *
				      * There is no need to consider
				      * hanging nodes here, since only
				      * one level is considered.
				      */
    template <int dim, class SparsityPattern>
    static void
    make_sparsity_pattern (const MGDoFHandler<dim> &dof_handler,
			   SparsityPattern         &sparsity,
			   const unsigned int       level);

				     /**
				      * Make a sparsity pattern including fluxes
				      * of discontinuous Galerkin methods.
				      * @see make_sparsity_pattern
				      * $see DoFTools
				      */
    template <int dim, class SparsityPattern>
    static void
    make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof_handler,
				SparsityPattern         &sparsity,
				const unsigned int       level);

				     /**
				      * Create sparsity pattern for
				      * the fluxes at refinement
				      * edges. 
				      * @see{make_flux_sparsity_pattern}
				      * @see{DoFTools}
				      */
    template <int dim, class SparsityPattern>
    static void
    make_flux_sparsity_pattern_edge (const MGDoFHandler<dim> &dof_handler,
				     SparsityPattern         &sparsity,
				     const unsigned int       level);

				     /**
				      * This function does the same as
				      * the other with the same name,
				      * but it gets two additional
				      * coefficient matrices. A matrix
				      * entry will only be generated
				      * for two basis functions, if
				      * there is a non-zero entry
				      * linking their associated
				      * components in the coefficient
				      * matrix.
				      *
				      * There is one matrix for
				      * couplings in a cell and one
				      * for the couplings occuring in
				      * fluxes.
				      */
    template <int dim, class SparsityPattern>
    static void
    make_flux_sparsity_pattern (const MGDoFHandler<dim> &dof,
				SparsityPattern       &sparsity,
				const unsigned int level,
				const FullMatrix<double>& int_mask,
				const FullMatrix<double>& flux_mask);

				     /**
				      * Count the dofs component-wise
				      * on each level.
				      *
				      * Result is a vector containing
				      * for each level a vector
				      * containing the number of dofs
				      * for each component (access is
				      * @p{result[level][component]}).
				      */
    template <int dim>
      static void count_dofs_per_component (const MGDoFHandler<dim>& mg_dof,
					    std::vector<std::vector<unsigned int> >& result);
    
    
				     /**
				      * Ajust vectors on all levels to
				      * correct size.  Here, we just
				      * count the numbers of degrees
				      * of freedom on each level and
				      * @p{reinit} each level vector
				      * to this length.
				      */
    template <int dim, typename number>
      static void
      reinit_vector (const MGDoFHandler<dim>& mg_dof,
		     MGLevelObject<Vector<number> >& vector);

				     /**
				      * Adjust block-vectors on all
				      * levels to correct size.  Count
				      * the numbers of degrees of
				      * freedom on each level
				      * component-wise. Then, assign
				      * each block of @p{vector} the
				      * corresponding size.
				      *
				      * The boolean field @p{selected}
				      * allows restricting this
				      * operation to certain
				      * components. In this case,
				      * @p{vector} will only have as
				      * many blocks as there are true
				      * values in @p{selected} (no
				      * blocks of length zero are
				      * padded in).
				      *
				      * Degrees of freedom must be
				      * sorted by component in order
				      * to obtain reasonable results
				      * from this function.
				      */
    template <int dim, typename number>
      static void
      reinit_vector (const MGDoFHandler<dim>& mg_dof,
		     MGLevelObject<BlockVector<number> >& vector,
		     const std::vector<bool>& selected);
};


#endif
