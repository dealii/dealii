// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2014 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__mg_tools_h
#define __deal2__mg_tools_h

#include <deal.II/base/config.h>
#include <deal.II/base/index_set.h>
#include <deal.II/dofs/dof_tools.h>
#include <deal.II/dofs/dof_handler.h>

#include <vector>
#include <set>


DEAL_II_NAMESPACE_OPEN

template <class Object> class MGLevelObject;
template <int dim, int spacedim> class DoFHandler;
template <typename number> class Vector;
template <typename number> class SparseMatrix;
template <typename number> class BlockVector;
template <typename number> class BlockSparseMatrix;
template <typename number> class FullMatrix;
template <typename number> class BlockSparseMatrix;

/* !@addtogroup mg */
/* @{ */

/**
 * This is a collection of functions operating on, and manipulating
 * the numbers of degrees of freedom in a multilevel triangulation. It
 * is similar in purpose and function to the @p DoFTools class, but
 * operates on levels of DoFHandler
 * objects. See there and the documentation of the member functions
 * for more information.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999 - 2005, 2012
 */
namespace MGTools
{
  /**
   * Compute row length vector for
   * multilevel methods.
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim,spacedim> &dofs,
                            const unsigned int level,
                            std::vector<unsigned int> &row_lengths,
                            const DoFTools::Coupling flux_couplings = DoFTools::none);

  /**
   * Compute row length vector for
   * multilevel methods with
   * optimization for block
   * couplings.
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim,spacedim> &dofs,
                            const unsigned int level,
                            std::vector<unsigned int> &row_lengths,
                            const Table<2,DoFTools::Coupling> &couplings,
                            const Table<2,DoFTools::Coupling> &flux_couplings);

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
  template <class DH, class SparsityPattern>
  void
  make_sparsity_pattern (const DH &dof_handler,
                         SparsityPattern         &sparsity,
                         const unsigned int       level);

  /**
   * Make a sparsity pattern including fluxes
   * of discontinuous Galerkin methods.
   * @ref make_sparsity_pattern
   * @ref DoFTools
   */
  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern (const DoFHandler<dim,spacedim> &dof_handler,
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
  void
  make_flux_sparsity_pattern_edge (const DoFHandler<dim,spacedim> &dof_handler,
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
   * for the couplings occurring in
   * fluxes.
   */
  template <int dim, class SparsityPattern, int spacedim>
  void
  make_flux_sparsity_pattern (const DoFHandler<dim,spacedim> &dof,
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
  void
  make_flux_sparsity_pattern_edge (const DoFHandler<dim,spacedim> &dof_handler,
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
  template <class DH>
  void
  count_dofs_per_block (const DH     &dof_handler,
                        std::vector<std::vector<types::global_dof_index> > &dofs_per_block,
                        std::vector<unsigned int>  target_block = std::vector<unsigned int>());

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
  void
  count_dofs_per_component (const DoFHandler<dim,spacedim> &mg_dof,
                            std::vector<std::vector<types::global_dof_index> > &result,
                            const bool only_once = false,
                            std::vector<unsigned int> target_component = std::vector<unsigned int>());

  /**
   * @deprecated Wrapper for the
   * other function with same name,
   * introduced for compatibility.
   */
  template <int dim, int spacedim>
  void
  count_dofs_per_component (const DoFHandler<dim,spacedim> &mg_dof,
                            std::vector<std::vector<types::global_dof_index> > &result,
                            std::vector<unsigned int> target_component) DEAL_II_DEPRECATED;

  /**
   * Generate a list of those degrees of freedom at the boundary of the domain
   * that should be eliminated from the matrix because they will be
   * constrained by Dirichlet boundary conditions.
   *
   * This is the multilevel equivalent of
   * VectorTools::interpolate_boundary_values, but since the multilevel method
   * does not have its own right hand side, the function values returned by
   * the function object that is part of the function_map argument are
   * ignored.
   *
   * @arg <tt>boundary_indices</tt> is a vector which on return contains all
   * indices of degrees of freedom for each level that are at the part of the
   * boundary identified by the function_map argument. Its length has to match
   * the number of levels in the dof handler object.
   */
  template <int dim, int spacedim>
  void
  make_boundary_list (const DoFHandler<dim,spacedim>      &mg_dof,
                      const typename FunctionMap<dim>::type &function_map,
                      std::vector<std::set<types::global_dof_index> > &boundary_indices,
                      const ComponentMask                   &component_mask = ComponentMask());

  /**
   * The same function as above, but return
   * an IndexSet rather than a
   * std::set<unsigned int> on each level.
   */
  template <int dim, int spacedim>
  void
  make_boundary_list (const DoFHandler<dim,spacedim>      &mg_dof,
                      const typename FunctionMap<dim>::type &function_map,
                      std::vector<IndexSet>                 &boundary_indices,
                      const ComponentMask               &component_mask = ComponentMask());

  /**
   * @deprecated
   */
  template <typename number>
  void
  apply_boundary_values (const std::set<types::global_dof_index> &boundary_dofs,
                         SparseMatrix<number> &matrix,
                         const bool preserve_symmetry,
                         const bool ignore_zeros = false) DEAL_II_DEPRECATED;

  /**
   * @deprecated
   */
  template <typename number>
  void
  apply_boundary_values (const std::set<types::global_dof_index> &boundary_dofs,
                         BlockSparseMatrix<number> &matrix,
                         const bool preserve_symmetry) DEAL_II_DEPRECATED;

  /**
   * For each level in a multigrid hierarchy, produce an IndexSet
   * that indicates which of the degrees of freedom are along
   * interfaces of this level to cells that only exist on coarser
   * levels.
   */
  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                                std::vector<IndexSet>  &interface_dofs);

  /**
   * As above but with a deprecated data structure. This makes one additional copy.
   */
  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                                std::vector<std::vector<bool> >  &interface_dofs) DEAL_II_DEPRECATED;


  template <int dim, int spacedim>
  void
  extract_non_interface_dofs (const DoFHandler<dim,spacedim> &mg_dof_handler,
                              std::vector<std::set<types::global_dof_index> > &non_interface_dofs);
}

/* @} */

DEAL_II_NAMESPACE_CLOSE

#endif
