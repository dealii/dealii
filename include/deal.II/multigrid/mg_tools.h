// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_mg_tools_h
#define dealii_mg_tools_h

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <set>
#include <vector>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class DoFHandler;
class MGConstrainedDoFs;
#endif

/* !@addtogroup mg */
/* @{ */

/**
 * This is a collection of functions operating on, and manipulating the
 * numbers of degrees of freedom in a multilevel triangulation. It is similar
 * in purpose and function to the @p DoFTools namespace, but operates on
 * levels of DoFHandler objects. See there and the documentation of the member
 * functions for more information.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, 1999 - 2005, 2012
 */
namespace MGTools
{
  /**
   * Compute row length vector for multilevel methods.
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(
    const DoFHandler<dim, spacedim> &dofs,
    const unsigned int               level,
    std::vector<unsigned int> &      row_lengths,
    const DoFTools::Coupling         flux_couplings = DoFTools::none);

  /**
   * Compute row length vector for multilevel methods with optimization for
   * block couplings.
   */
  template <int dim, int spacedim>
  void
  compute_row_length_vector(const DoFHandler<dim, spacedim> &   dofs,
                            const unsigned int                  level,
                            std::vector<unsigned int> &         row_lengths,
                            const Table<2, DoFTools::Coupling> &couplings,
                            const Table<2, DoFTools::Coupling> &flux_couplings);

  /**
   * Write the sparsity structure of the matrix belonging to the specified @p
   * level. The sparsity pattern is not compressed, so before creating the
   * actual matrix you have to compress the matrix yourself, using
   * <tt>SparseMatrixStruct::compress()</tt>.
   *
   * There is no need to consider hanging nodes here, since only one level is
   * considered.
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_sparsity_pattern(const DoFHandlerType &dof_handler,
                        SparsityPatternType & sparsity,
                        const unsigned int    level);

  /**
   * Make a sparsity pattern including fluxes of discontinuous Galerkin
   * methods.
   * @see
   * @ref make_sparsity_pattern
   * and
   * @ref DoFTools
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &dof_handler,
                             SparsityPatternType &            sparsity,
                             const unsigned int               level);

  /**
   * Create sparsity pattern for the fluxes at refinement edges. The matrix
   * maps a function of the fine level space @p level to the coarser space.
   *
   * make_flux_sparsity_pattern()
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern_edge(const DoFHandler<dim, spacedim> &dof_handler,
                                  SparsityPatternType &            sparsity,
                                  const unsigned int               level);
  /**
   * This function does the same as the other with the same name, but it gets
   * two additional coefficient matrices. A matrix entry will only be
   * generated for two basis functions, if there is a non-zero entry linking
   * their associated components in the coefficient matrix.
   *
   * There is one matrix for couplings in a cell and one for the couplings
   * occurring in fluxes.
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern(const DoFHandler<dim, spacedim> &   dof,
                             SparsityPatternType &               sparsity,
                             const unsigned int                  level,
                             const Table<2, DoFTools::Coupling> &int_mask,
                             const Table<2, DoFTools::Coupling> &flux_mask);

  /**
   * Create sparsity pattern for the fluxes at refinement edges. The matrix
   * maps a function of the fine level space @p level to the coarser space.
   * This is the version restricting the pattern to the elements actually
   * needed.
   *
   * make_flux_sparsity_pattern()
   */
  template <int dim, typename SparsityPatternType, int spacedim>
  void
  make_flux_sparsity_pattern_edge(
    const DoFHandler<dim, spacedim> &   dof_handler,
    SparsityPatternType &               sparsity,
    const unsigned int                  level,
    const Table<2, DoFTools::Coupling> &flux_mask);


  /**
   * Create sparsity pattern for interface_in/out matrices used in a multigrid
   * computation. These matrices contain an entry representing the coupling of
   * degrees of freedom on a refinement edge to those not on the refinement edge
   * of a certain level.
   */
  template <typename DoFHandlerType, typename SparsityPatternType>
  void
  make_interface_sparsity_pattern(const DoFHandlerType &   dof_handler,
                                  const MGConstrainedDoFs &mg_constrained_dofs,
                                  SparsityPatternType &    sparsity,
                                  const unsigned int       level);


  /**
   * Count the dofs block-wise on each level.
   *
   * Result is a vector containing for each level a vector containing the
   * number of dofs for each block (access is <tt>result[level][block]</tt>).
   */
  template <typename DoFHandlerType>
  void
  count_dofs_per_block(
    const DoFHandlerType &                             dof_handler,
    std::vector<std::vector<types::global_dof_index>> &dofs_per_block,
    std::vector<unsigned int>                          target_block = {});

  /**
   * Count the dofs component-wise on each level.
   *
   * Result is a vector containing for each level a vector containing the
   * number of dofs for each component (access is
   * <tt>result[level][component]</tt>).
   */
  template <int dim, int spacedim>
  void
  count_dofs_per_component(
    const DoFHandler<dim, spacedim> &                  mg_dof,
    std::vector<std::vector<types::global_dof_index>> &result,
    const bool                                         only_once        = false,
    std::vector<unsigned int>                          target_component = {});

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
   *
   * Previous content in @p boundary_indices is not overwritten,
   * but added to.
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(
    const DoFHandler<dim, spacedim> &mg_dof,
    const std::map<types::boundary_id, const Function<spacedim> *>
      &                                             function_map,
    std::vector<std::set<types::global_dof_index>> &boundary_indices,
    const ComponentMask &component_mask = ComponentMask());

  /**
   * The same function as above, but return an IndexSet rather than a
   * std::set<unsigned int> on each level.
   *
   * Previous content in @p boundary_indices is not overwritten,
   * but added to.
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim> &           mg_dof,
                     const std::map<types::boundary_id,
                                    const Function<spacedim> *> &function_map,
                     std::vector<IndexSet> &boundary_indices,
                     const ComponentMask &  component_mask = ComponentMask());

  /**
   * The same function as above, but return an IndexSet rather than a
   * std::set<unsigned int> on each level and use a std::set of boundary_ids
   * as input.
   *
   * Previous content in @p boundary_indices is not overwritten, but added to.
   */
  template <int dim, int spacedim>
  void
  make_boundary_list(const DoFHandler<dim, spacedim> &   mg_dof,
                     const std::set<types::boundary_id> &boundary_ids,
                     std::vector<IndexSet> &             boundary_indices,
                     const ComponentMask &component_mask = ComponentMask());

  /**
   * For each level in a multigrid hierarchy, produce an IndexSet that
   * indicates which of the degrees of freedom are along interfaces of this
   * level to cells that only exist on coarser levels.
   */
  template <int dim, int spacedim>
  void
  extract_inner_interface_dofs(const DoFHandler<dim, spacedim> &mg_dof_handler,
                               std::vector<IndexSet> &          interface_dofs);

  /**
   * For each level in a multigrid hierarchy, produce a std::set of degrees of
   * freedoms that are not located along interfaces of this level to cells that
   * only exist on coarser levels.
   *
   * @deprecated Use extract_inner_interface_dofs() for computing the complement
   * of degrees of freedoms instead.
   */
  template <int dim, int spacedim>
  DEAL_II_DEPRECATED void
  extract_non_interface_dofs(
    const DoFHandler<dim, spacedim> &               mg_dof_handler,
    std::vector<std::set<types::global_dof_index>> &non_interface_dofs);

  /**
   * Return the highest possible level that can be used as the coarsest level in
   * a Multigrid computation, that is, the highest level in the hierarchy whose
   * mesh covers the entire domain. This corresponds to the minimum level of a
   * cell on the active mesh. Since each processor only has a local view of the
   * mesh, each processor must call this function. Note that this is a global
   * minimum over the entire mesh and therefore each processor will return the
   * same value.
   */
  template <int dim, int spacedim>
  unsigned int
  max_level_for_coarse_mesh(const Triangulation<dim, spacedim> &tria);

  /**
   * Return the imbalance of the parallel distribution of the multigrid
   * mesh hierarchy. Ideally this value is equal to 1 (every processor owns
   * the same number of cells on each level, approximately true for most
   * globally refined meshes). Values greater than 1 estimate the slowdown
   * one should see in a geometric multigrid v-cycle as compared with the same
   * computation on a perfectly distributed mesh hierarchy.
   *
   * This function is a collective MPI call between all ranks of the
   * Triangulation and therefore needs to be called from all ranks.
   *
   * @note This function requires that
   * parallel::TriangulationBase::is_multilevel_hierarchy_constructed()
   * is true, which can be controlled by setting the
   * construct_multigrid_hierarchy flag when constructing the
   * Triangulation.
   */
  template <int dim, int spacedim>
  double
  workload_imbalance(const Triangulation<dim, spacedim> &tria);

} // namespace MGTools

/* @} */

DEAL_II_NAMESPACE_CLOSE

#endif
