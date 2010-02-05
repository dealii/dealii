//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2007, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__mg_tools_h
#define __deal2__mg_tools_h

// This file moved here from mg_dof_tools.h Revision 1.36

#include <base/config.h>
#include <dofs/dof_tools.h>
#include <multigrid/mg_dof_handler.h>

#include <vector>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <class Object> class MGLevelObject;
template <int dim, int spacedim> class MGDoFHandler;
template <typename number> class Vector;
template <typename number> class SparseMatrix;
template <typename number> class BlockVector;
template <typename number> class BlockSparseMatrix;
template <typename number> class FullMatrix;
template <typename number> class BlockSparseMatrix;

/*!@addtogroup mg */
/*@{*/

/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom in a multilevel triangulation. It
 * is similar in purpose and function to the @p DoFTools class, but
 * operates on @p MGDoFHandler objects instead of DoFHandler
 * objects. See there and the documentation of the member functions
 * for more information.
 *
 * All member functions are static, so there is no need to create an
 * object of class MGTools.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999 - 2005
 */
class MGTools
{
  public:
				     /**
				      * Compute row length vector for
				      * multilevel methods.
				      */
    template <int dim, int spacedim>
    static
    void compute_row_length_vector(
      const MGDoFHandler<dim,spacedim>& dofs,
      const unsigned int level,
      std::vector<unsigned int>& row_lengths,
      const DoFTools::Coupling flux_couplings = DoFTools::none);

				     /**
				      * Compute row length vector for
				      * multilevel methods with
				      * optimization for block
				      * couplings.
				      */
    template <int dim, int spacedim>
    static
    void compute_row_length_vector(
      const MGDoFHandler<dim,spacedim>& dofs,
      const unsigned int level,
      std::vector<unsigned int>& row_lengths,
      const Table<2,DoFTools::Coupling>& couplings,
      const Table<2,DoFTools::Coupling>& flux_couplings);

				     /**
				      * Write the sparsity structure
				      * of the matrix belonging to the
				      * specified @p level. The sparsity pattern
				      * is not compressed, so before
				      * creating the actual matrix
				      * you have to compress the
				      * matrix yourself, using
				      * <tt>SparseMatrixStruct::compress()</tt>.
				      *
				      * There is no need to consider
				      * hanging nodes here, since only
				      * one level is considered.
				      */
    template <int dim, class SparsityPattern, int spacedim>
    static void
    make_sparsity_pattern (const MGDoFHandler<dim,spacedim> &dof_handler,
			   SparsityPattern         &sparsity,
			   const unsigned int       level);

				     /**
				      * Make a sparsity pattern including fluxes
				      * of discontinuous Galerkin methods.
				      * @ref make_sparsity_pattern
				      * @ref DoFTools
				      */
    template <int dim, class SparsityPattern, int spacedim>
    static void
    make_flux_sparsity_pattern (const MGDoFHandler<dim,spacedim> &dof_handler,
				SparsityPattern         &sparsity,
				const unsigned int       level);

				     /**
				      * Create sparsity pattern for
				      * the fluxes at refinement
				      * edges. The matrix maps a
				      * function of the fine level
				      * space @p level to the coarser
				      * space.
				      *
				      * make_flux_sparsity_pattern()
				      */
    template <int dim, class SparsityPattern, int spacedim>
    static void
    make_flux_sparsity_pattern_edge (const MGDoFHandler<dim,spacedim> &dof_handler,
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
    template <int dim, class SparsityPattern, int spacedim>
    static void
    make_flux_sparsity_pattern (const MGDoFHandler<dim,spacedim> &dof,
				SparsityPattern       &sparsity,
				const unsigned int level,
				const Table<2,DoFTools::Coupling> &int_mask,
				const Table<2,DoFTools::Coupling> &flux_mask);

				     /**
				      * Create sparsity pattern for
				      * the fluxes at refinement
				      * edges. The matrix maps a
				      * function of the fine level
				      * space @p level to the coarser
				      * space. This is the version
				      * restricting the pattern to the
				      * elements actually needed.
				      *
				      * make_flux_sparsity_pattern()
				      */
    template <int dim, class SparsityPattern, int spacedim>
    static void
    make_flux_sparsity_pattern_edge (const MGDoFHandler<dim,spacedim> &dof_handler,
				     SparsityPattern         &sparsity,
				     const unsigned int       level,
				     const Table<2,DoFTools::Coupling> &flux_mask);

				     /**
				      * Count the dofs block-wise
				      * on each level.
				      *
				      * Result is a vector containing
				      * for each level a vector
				      * containing the number of dofs
				      * for each block (access is
				      * <tt>result[level][block]</tt>).
				      */
    template <int dim, int spacedim>
      static void count_dofs_per_block (
	const MGDoFHandler<dim,spacedim> &mg_dof,
	std::vector<std::vector<unsigned int> > &result,
	std::vector<unsigned int> target_block = std::vector<unsigned int>());

				     /**
				      * Count the dofs component-wise
				      * on each level.
				      *
				      * Result is a vector containing
				      * for each level a vector
				      * containing the number of dofs
				      * for each component (access is
				      * <tt>result[level][component]</tt>).
				      */
    template <int dim, int spacedim>
      static void count_dofs_per_component (
	const MGDoFHandler<dim,spacedim> &mg_dof,
	std::vector<std::vector<unsigned int> > &result,
	const bool only_once = false,
	std::vector<unsigned int> target_component = std::vector<unsigned int>());

				     /**
				      * @deprecated Wrapper for the
				      * other function with same name,
				      * introduced for compatibility.
				      */
    template <int dim, int spacedim>
      static void count_dofs_per_component (
	const MGDoFHandler<dim,spacedim> &mg_dof,
	std::vector<std::vector<unsigned int> > &result,
	std::vector<unsigned int> target_component);

				     /**
				      * Generate a list of those
				      * degrees of freedom at the
				      * boundary which should be
				      * eliminated from the matrix.
				      *
				      * This is the multilevel
				      * equivalent of
				      * VectorTools::interpolate_boundary_values,
				      * but since the multilevel
				      * method does not have its own
				      * right hand side, the function
				      * values are ignored.
				      *
				      * @arg <tt>boundary_indices</tt>
				      * is a vector which on return
				      * contains all indices of
				      * boundary constraint degrees of
				      * freedom for each level. Its
				      * length has to match the number
				      * of levels.
				      */
    template <int dim, int spacedim>
    static void make_boundary_list(
      const MGDoFHandler<dim,spacedim>& mg_dof,
      const typename FunctionMap<dim>::type& function_map,
      std::vector<std::set<unsigned int> >& boundary_indices,
      const std::vector<bool>& component_mask = std::vector<bool>());
                                      /**
                                       * Maybe no longer needed.
                                       */

    template <typename number>
    static void apply_boundary_values (
      const std::set<unsigned int> &boundary_dofs,
      SparseMatrix<number>& matrix,
      const bool preserve_symmetry,
      const bool ignore_zeros = false);

    template <typename number>
    static void apply_boundary_values (
      const std::set<unsigned int>& boundary_dofs,
      BlockSparseMatrix<number>& matrix,
      const bool preserve_symmetry);


				     /**
				      * For each level in a multigrid
				      * hierarchy, produce a boolean
				      * mask that indicates which of
				      * the degrees of freedom are
				      * along interfaces of this level
				      * to cells that only exist on
				      * coarser levels. The function
				      * returns the subset of these
				      * indices in the last argument
				      * that are not only on interior
				      * interfaces (i.e. between cells
				      * of a given level and adjacent
				      * coarser levels) but also on
				      * the external boundary of the
				      * domain.
				      */
    template <int dim, int spacedim>
    static
    void
    extract_inner_interface_dofs (const MGDoFHandler<dim,spacedim> &mg_dof_handler,
				  std::vector<std::vector<bool> >  &interface_dofs,
				  std::vector<std::vector<bool> >  &boundary_interface_dofs);

				     /**
                                      * Does the same as the function above, 
                                      * but fills only the interface_dofs.
				      */
    template <int dim, int spacedim>
    static
    void
    extract_inner_interface_dofs (const MGDoFHandler<dim,spacedim> &mg_dof_handler,
				  std::vector<std::vector<bool> >  &interface_dofs);
};

/*@}*/

DEAL_II_NAMESPACE_CLOSE

#endif
