// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2020 by the deal.II authors
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

#ifndef dealii_mg_transfer_matrix_free_h
#define dealii_mg_transfer_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/multigrid/mg_transfer_internal.h>


DEAL_II_NAMESPACE_OPEN


/*!@addtogroup mg */
/*@{*/

/**
 * Implementation of the MGTransferBase interface for which the transfer
 * operations is implemented in a matrix-free way based on the interpolation
 * matrices of the underlying finite element. This requires considerably less
 * memory than MGTransferPrebuilt and can also be considerably faster than
 * that variant.
 *
 * This class currently only works for tensor-product finite elements based on
 * FE_Q and FE_DGQ elements, including systems involving multiple components
 * of one of these elements. Systems with different elements or other elements
 * are currently not implemented.
 *
 * @author Martin Kronbichler
 * @date 2016
 */
template <int dim, typename Number>
class MGTransferMatrixFree
  : public MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferMatrixFree();

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   */
  MGTransferMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Destructor.
   */
  virtual ~MGTransferMatrixFree() override = default;

  /**
   * Initialize the constraints to be used in build().
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void
  clear();

  /**
   * Actually build the information for the prolongation for each level.
   *
   * The optional second argument of external partitioners allows the user to
   * suggest vector partitioning on the levels. In case the partitioners
   * are found to contain all ghost unknowns that are visited through the
   * transfer, the given partitioners are chosen. This ensures compatibility
   * of vectors during prolongate and restrict with external partitioners as
   * given by the user, which in turn saves some copy operations. However, in
   * case there are unknowns missing -- and this is typically the case at some
   * point during h-coarsening since processors will need to drop out and
   * thus children's unknowns on some processor will be needed as ghosts to a
   * parent cell on another processor -- the provided external partitioners are
   * ignored and internal variants are used instead.
   */
  void
  build(const DoFHandler<dim, dim> &dof_handler,
        const std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>
          &external_partitioners =
            std::vector<std::shared_ptr<const Utilities::MPI::Partitioner>>());

  /**
   * Prolongate a vector from level <tt>to_level-1</tt> to level
   * <tt>to_level</tt> using the embedding matrices of the underlying finite
   * element. The previous content of <tt>dst</tt> is overwritten.
   *
   * @param to_level The index of the level to prolongate to, which is the
   * level of @p dst.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the coarser level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * finer level.
   */
  virtual void
  prolongate(
    const unsigned int                                to_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Restrict a vector from level <tt>from_level</tt> to level
   * <tt>from_level-1</tt> using the transpose operation of the prolongate()
   * method. If the region covered by cells on level <tt>from_level</tt> is
   * smaller than that of level <tt>from_level-1</tt> (local refinement), then
   * some degrees of freedom in <tt>dst</tt> are active and will not be
   * altered. For the other degrees of freedom, the result of the restriction
   * is added.
   *
   * @param from_level The index of the level to restrict from, which is the
   * level of @p src.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the finer level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * coarser level.
   */
  virtual void
  restrict_and_add(
    const unsigned int                                from_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const override;

  /**
   * Restrict fine-mesh field @p src to each multigrid level in @p dof_handler and
   * store the result in @p dst.
   *
   * The argument @p dst has to be initialized with the correct size according
   * to the number of levels of the triangulation.
   *
   * If an inner vector of @p dst is empty or has incorrect locally owned size,
   * it will be resized to locally relevant degrees of freedom on each level.
   */
  template <typename Number2, int spacedim>
  void
  interpolate_to_mg(
    const DoFHandler<dim, spacedim> &                          dof_handler,
    MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
    const LinearAlgebra::distributed::Vector<Number2> &        src) const;

  /**
   * Finite element does not provide prolongation matrices.
   */
  DeclException0(ExcNoProlongation);

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * A variable storing the degree of the finite element contained in the
   * DoFHandler passed to build(). The selection of the computational kernel is
   * based on this number.
   */
  unsigned int fe_degree;

  /**
   * A variable storing whether the element is continuous and there is a joint
   * degree of freedom in the center of the 1D line.
   */
  bool element_is_continuous;

  /**
   * A variable storing the number of components in the finite element contained
   * in the DoFHandler passed to build().
   */
  unsigned int n_components;

  /**
   * A variable storing the number of degrees of freedom on all child cells. It
   * is <tt>2<sup>dim</sup>*fe.dofs_per_cell</tt> for DG elements and somewhat
   * less for continuous elements.
   */
  unsigned int n_child_cell_dofs;

  /**
   * This variable holds the indices for cells on a given level, extracted from
   * DoFHandler for fast access. All DoF indices on a given level are stored as
   * a plain array (since this class assumes constant DoFs per cell). To index
   * into this array, use the cell number times dofs_per_cell.
   *
   * This array first is arranged such that all locally owned level cells come
   * first (found in the variable n_owned_level_cells) and then other cells
   * necessary for the transfer to the next level.
   */
  std::vector<std::vector<unsigned int>> level_dof_indices;

  /**
   * A variable storing the connectivity from parent to child cell numbers for
   * each level.
   */
  std::vector<std::vector<std::pair<unsigned int, unsigned int>>>
    parent_child_connect;

  /**
   * A variable storing the number of cells owned on a given process (sets the
   * bounds for the worker loops) for each level.
   */
  std::vector<unsigned int> n_owned_level_cells;

  /**
   * This variable holds the one-dimensional embedding (prolongation) matrix
   * from mother element to all the children.
   */
  AlignedVector<VectorizedArray<Number>> prolongation_matrix_1d;

  /**
   * This variable holds the temporary values for the tensor evaluation
   */
  mutable AlignedVector<VectorizedArray<Number>> evaluation_data;

  /**
   * For continuous elements, restriction is not additive and we need to
   * weight the result at the end of prolongation (and at the start of
   * restriction) by the valence of the degrees of freedom, i.e., on how many
   * elements they appear. We store the data in vectorized form to allow for
   * cheap access. Moreover, we utilize the fact that we only need to store
   * <tt>3<sup>dim</sup></tt> indices.
   *
   * Data is organized in terms of each level (outer vector) and the cells on
   * each level (inner vector).
   */
  std::vector<AlignedVector<VectorizedArray<Number>>> weights_on_refined;

  /**
   * A variable storing the local indices of Dirichlet boundary conditions on
   * cells for all levels (outer index), the cells within the levels (second
   * index), and the indices on the cell (inner index).
   */
  std::vector<std::vector<std::vector<unsigned short>>> dirichlet_indices;

  /**
   * A vector that holds shared pointers to the partitioners of the
   * transfer. These partitioners might be shared with what was passed in from
   * the outside through build() or be shared with the level vectors inherited
   * from MGLevelGlobalTransfer.
   */
  MGLevelObject<std::shared_ptr<const Utilities::MPI::Partitioner>>
    vector_partitioners;

  /**
   * Perform the prolongation operation.
   */
  template <int degree>
  void
  do_prolongate_add(
    const unsigned int                                to_level,
    LinearAlgebra::distributed::Vector<Number> &      dst,
    const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Performs the restriction operation.
   */
  template <int degree>
  void
  do_restrict_add(const unsigned int                                from_level,
                  LinearAlgebra::distributed::Vector<Number> &      dst,
                  const LinearAlgebra::distributed::Vector<Number> &src) const;
};


/**
 * Implementation of the MGTransferBase interface for which the transfer
 * operations is implemented in a matrix-free way based on the interpolation
 * matrices of the underlying finite element. This requires considerably less
 * memory than MGTransferPrebuilt and can also be considerably faster than
 * that variant.
 *
 * This class works with LinearAlgebra::distributed::BlockVector and
 * performs exactly the same transfer operations for each block as
 * MGTransferMatrixFree.
 * Both the cases that the same DoFHandler is used for all the blocks
 * and that each block uses its own DoFHandler are supported.
 *
 * @author Denis Davydov, Daniel Arndt
 * @date 2017
 */
template <int dim, typename Number>
class MGTransferBlockMatrixFree
  : public MGTransferBase<LinearAlgebra::distributed::BlockVector<Number>>
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferBlockMatrixFree() = default;

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   */
  MGTransferBlockMatrixFree(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  MGTransferBlockMatrixFree(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Destructor.
   */
  virtual ~MGTransferBlockMatrixFree() override = default;

  /**
   * Initialize the constraints to be used in build().
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  void
  initialize_constraints(
    const std::vector<MGConstrainedDoFs> &mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void
  clear();

  /**
   * Actually build the information for the prolongation for each level.
   */
  void
  build(const DoFHandler<dim, dim> &dof_handler);

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  void
  build(const std::vector<const DoFHandler<dim, dim> *> &dof_handler);

  /**
   * Prolongate a vector from level <tt>to_level-1</tt> to level
   * <tt>to_level</tt> using the embedding matrices of the underlying finite
   * element. The previous content of <tt>dst</tt> is overwritten.
   *
   * @param to_level The index of the level to prolongate to, which is the
   * level of @p dst.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the coarser level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * finer level.
   */
  virtual void
  prolongate(
    const unsigned int                                     to_level,
    LinearAlgebra::distributed::BlockVector<Number> &      dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * Restrict a vector from level <tt>from_level</tt> to level
   * <tt>from_level-1</tt> using the transpose operation of the prolongate()
   * method. If the region covered by cells on level <tt>from_level</tt> is
   * smaller than that of level <tt>from_level-1</tt> (local refinement), then
   * some degrees of freedom in <tt>dst</tt> are active and will not be
   * altered. For the other degrees of freedom, the result of the restriction
   * is added.
   *
   * @param from_level The index of the level to restrict from, which is the
   * level of @p src.
   *
   * @param src is a vector with as many elements as there are degrees of
   * freedom on the finer level involved.
   *
   * @param dst has as many elements as there are degrees of freedom on the
   * coarser level.
   */
  virtual void
  restrict_and_add(
    const unsigned int                                     from_level,
    LinearAlgebra::distributed::BlockVector<Number> &      dst,
    const LinearAlgebra::distributed::BlockVector<Number> &src) const override;

  /**
   * Transfer from a block-vector on the global grid to block-vectors defined
   * on each of the levels separately for active degrees of freedom.
   * In particular, for a globally refined mesh only the finest level in @p dst
   * is filled as a plain copy of @p src. All the other level objects are left
   * untouched.
   *
   * This function will initialize @p dst accordingly if needed as required by
   * the Multigrid class.
   */
  template <typename Number2, int spacedim>
  void
  copy_to_mg(
    const DoFHandler<dim, spacedim> &                               dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const LinearAlgebra::distributed::BlockVector<Number2> &        src) const;

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  template <typename Number2, int spacedim>
  void
  copy_to_mg(
    const std::vector<const DoFHandler<dim, spacedim> *> &          dof_handler,
    MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
    const LinearAlgebra::distributed::BlockVector<Number2> &        src) const;

  /**
   * Transfer from multi-level block-vector to normal vector.
   */
  template <typename Number2, int spacedim>
  void
  copy_from_mg(
    const DoFHandler<dim, spacedim> &                 dof_handler,
    LinearAlgebra::distributed::BlockVector<Number2> &dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * Same as above for the case that each block has its own DoFHandler.
   */
  template <typename Number2, int spacedim>
  void
  copy_from_mg(
    const std::vector<const DoFHandler<dim, spacedim> *> &dof_handler,
    LinearAlgebra::distributed::BlockVector<Number2> &    dst,
    const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
    const;

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * This class can both be used with a single DoFHandler
   * or a separate DoFHandler for each block.
   */
  static const bool supports_dof_handler_vector = true;

private:
  /**
   * Non-block matrix-free versions of transfer operation.
   */
  std::vector<MGTransferMatrixFree<dim, Number>> matrix_free_transfer_vector;

  /**
   * A flag to indicate whether the same DoFHandler is used for all
   * the components or if each block has its own DoFHandler.
   */
  const bool same_for_all;
};


/*@}*/


//------------------------ templated functions -------------------------
#ifndef DOXYGEN


template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferMatrixFree<dim, Number>::interpolate_to_mg(
  const DoFHandler<dim, spacedim> &                          dof_handler,
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
  const LinearAlgebra::distributed::Vector<Number2> &        src) const
{
  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  Assert(max_level == dof_handler.get_triangulation().n_global_levels() - 1,
         ExcDimensionMismatch(
           max_level, dof_handler.get_triangulation().n_global_levels() - 1));

  const parallel::TriangulationBase<dim, spacedim> *p_tria =
    (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
      &dof_handler.get_triangulation()));
  MPI_Comm mpi_communicator =
    p_tria != nullptr ? p_tria->get_communicator() : MPI_COMM_SELF;

  // resize the dst vector if it's empty or has incorrect size
  MGLevelObject<IndexSet> relevant_dofs(min_level, max_level);
  for (unsigned int level = min_level; level <= max_level; ++level)
    {
      DoFTools::extract_locally_relevant_level_dofs(dof_handler,
                                                    level,
                                                    relevant_dofs[level]);
      if (dst[level].size() !=
            dof_handler.locally_owned_mg_dofs(level).size() ||
          dst[level].local_size() !=
            dof_handler.locally_owned_mg_dofs(level).n_elements())
        dst[level].reinit(dof_handler.locally_owned_mg_dofs(level),
                          relevant_dofs[level],
                          mpi_communicator);
    }

  const FiniteElement<dim, spacedim> &fe = dof_handler.get_fe();

  // copy fine level vector to active cells in MG hierarchy
  this->copy_to_mg(dof_handler, dst, src, true);

  // FIXME: maybe need to store hanging nodes constraints per level?
  // MGConstrainedDoFs does NOT keep this info right now, only periodicity
  // constraints...
  dst[max_level].update_ghost_values();
  // do the transfer from level to level-1:
  for (unsigned int level = max_level; level > min_level; --level)
    {
      // auxiliary vector which always has ghost elements
      LinearAlgebra::distributed::Vector<Number> ghosted_vector(
        dof_handler.locally_owned_mg_dofs(level),
        relevant_dofs[level],
        mpi_communicator);
      ghosted_vector = dst[level];
      ghosted_vector.update_ghost_values();

      std::vector<Number>                  dof_values_coarse(fe.dofs_per_cell);
      Vector<Number>                       dof_values_fine(fe.dofs_per_cell);
      Vector<Number>                       tmp(fe.dofs_per_cell);
      std::vector<types::global_dof_index> dof_indices(fe.dofs_per_cell);
      typename DoFHandler<dim>::cell_iterator cell =
        dof_handler.begin(level - 1);
      typename DoFHandler<dim>::cell_iterator endc = dof_handler.end(level - 1);
      for (; cell != endc; ++cell)
        if (cell->is_locally_owned_on_level())
          {
            // if we get to a cell without children (== active), we can
            // skip it as there values should be already set by the
            // equivalent of copy_to_mg()
            if (cell->is_active())
              continue;

            std::fill(dof_values_coarse.begin(), dof_values_coarse.end(), 0.);
            for (unsigned int child = 0; child < cell->n_children(); ++child)
              {
                cell->child(child)->get_mg_dof_indices(dof_indices);
                for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                  dof_values_fine(i) = ghosted_vector(dof_indices[i]);
                fe.get_restriction_matrix(child, cell->refinement_case())
                  .vmult(tmp, dof_values_fine);
                for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
                  if (fe.restriction_is_additive(i))
                    dof_values_coarse[i] += tmp[i];
                  else if (tmp(i) != 0.)
                    dof_values_coarse[i] = tmp[i];
              }
            cell->get_mg_dof_indices(dof_indices);
            for (unsigned int i = 0; i < fe.dofs_per_cell; ++i)
              dst[level - 1](dof_indices[i]) = dof_values_coarse[i];
          }

      dst[level - 1].compress(VectorOperation::insert);
      dst[level - 1].update_ghost_values();
    }
}



template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_to_mg(
  const DoFHandler<dim, spacedim> &                               dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const LinearAlgebra::distributed::BlockVector<Number2> &        src) const
{
  AssertDimension(matrix_free_transfer_vector.size(), 1);
  Assert(same_for_all,
         ExcMessage(
           "This object was initialized with support for usage with one "
           "DoFHandler for each block, but this method assumes that "
           "the same DoFHandler is used for all the blocks!"));
  const std::vector<const DoFHandler<dim, spacedim> *> mg_dofs(src.n_blocks(),
                                                               &dof_handler);

  copy_to_mg(mg_dofs, dst, src);
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_to_mg(
  const std::vector<const DoFHandler<dim, spacedim> *> &          dof_handler,
  MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
  const LinearAlgebra::distributed::BlockVector<Number2> &        src) const
{
  const unsigned int n_blocks = src.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  // this function is normally called within the Multigrid class with
  // dst == defect level block vector. At first run this vector is not
  // initialized. Do this below:
  {
    const parallel::TriangulationBase<dim, spacedim> *tria =
      (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
        &(dof_handler[0]->get_triangulation())));
    for (unsigned int i = 1; i < n_blocks; ++i)
      AssertThrow(
        (dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(
           &(dof_handler[0]->get_triangulation())) == tria),
        ExcMessage("The DoFHandler use different Triangulations!"));

    MGLevelObject<bool> do_reinit;
    do_reinit.resize(min_level, max_level);
    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        do_reinit[level] = false;
        if (dst[level].n_blocks() != n_blocks)
          {
            do_reinit[level] = true;
            continue; // level
          }
        for (unsigned int b = 0; b < n_blocks; ++b)
          {
            LinearAlgebra::distributed::Vector<Number> &v = dst[level].block(b);
            if (v.size() !=
                  dof_handler[b]->locally_owned_mg_dofs(level).size() ||
                v.local_size() !=
                  dof_handler[b]->locally_owned_mg_dofs(level).n_elements())
              {
                do_reinit[level] = true;
                break; // b
              }
          }
      }

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        if (do_reinit[level])
          {
            dst[level].reinit(n_blocks);
            for (unsigned int b = 0; b < n_blocks; ++b)
              {
                LinearAlgebra::distributed::Vector<Number> &v =
                  dst[level].block(b);
                v.reinit(dof_handler[b]->locally_owned_mg_dofs(level),
                         tria != nullptr ? tria->get_communicator() :
                                           MPI_COMM_SELF);
              }
            dst[level].collect_sizes();
          }
        else
          dst[level] = 0;
      }
  }

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> dst_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        dst_non_block[l].reinit(dst[l].block(b));
      const unsigned int data_block = same_for_all ? 0 : b;
      matrix_free_transfer_vector[data_block].copy_to_mg(*dof_handler[b],
                                                         dst_non_block,
                                                         src.block(b));

      for (unsigned int l = min_level; l <= max_level; ++l)
        dst[l].block(b) = dst_non_block[l];
    }
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_from_mg(
  const DoFHandler<dim, spacedim> &                 dof_handler,
  LinearAlgebra::distributed::BlockVector<Number2> &dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  AssertDimension(matrix_free_transfer_vector.size(), 1);
  const std::vector<const DoFHandler<dim, spacedim> *> mg_dofs(dst.n_blocks(),
                                                               &dof_handler);

  copy_from_mg(mg_dofs, dst, src);
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim, Number>::copy_from_mg(
  const std::vector<const DoFHandler<dim, spacedim> *> &dof_handler,
  LinearAlgebra::distributed::BlockVector<Number2> &    dst,
  const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src)
  const
{
  const unsigned int n_blocks = dst.n_blocks();
  AssertDimension(dof_handler.size(), n_blocks);

  if (n_blocks == 0)
    return;

  const unsigned int min_level = src.min_level();
  const unsigned int max_level = src.max_level();

  for (unsigned int l = min_level; l <= max_level; ++l)
    AssertDimension(src[l].n_blocks(), dst.n_blocks());

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> src_non_block(
    min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          src_non_block[l].reinit(src[l].block(b));
          src_non_block[l] = src[l].block(b);
        }
      const unsigned int data_block = same_for_all ? 0 : b;
      matrix_free_transfer_vector[data_block].copy_from_mg(*dof_handler[b],
                                                           dst.block(b),
                                                           src_non_block);
    }
}



#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
