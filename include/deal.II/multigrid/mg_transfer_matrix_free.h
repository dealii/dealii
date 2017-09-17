// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2017 by the deal.II authors
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

#ifndef dealii_mg_transfer_matrix_free_h
#define dealii_mg_transfer_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/la_parallel_block_vector.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/dofs/dof_handler.h>


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
class MGTransferMatrixFree : public MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number> >
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferMatrixFree ();

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   */
  MGTransferMatrixFree (const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Destructor.
   */
  virtual ~MGTransferMatrixFree () = default;

  /**
   * Initialize the constraints to be used in build().
   */
  void initialize_constraints (const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void clear ();

  /**
   * Actually build the information for the prolongation for each level.
   */
  void build (const DoFHandler<dim,dim> &mg_dof);

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
  virtual void prolongate (const unsigned int                           to_level,
                           LinearAlgebra::distributed::Vector<Number>       &dst,
                           const LinearAlgebra::distributed::Vector<Number> &src) const;

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
  virtual void restrict_and_add (const unsigned int from_level,
                                 LinearAlgebra::distributed::Vector<Number>       &dst,
                                 const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Finite element does not provide prolongation matrices.
   */
  DeclException0(ExcNoProlongation);

  /**
   * Memory used by this object.
   */
  std::size_t memory_consumption () const;

private:

  /**
   * Stores the degree of the finite element contained in the DoFHandler
   * passed to build(). The selection of the computational kernel is based on
   * this number.
   */
  unsigned int fe_degree;

  /**
   * Stores whether the element is continuous and there is a joint degree of
   * freedom in the center of the 1D line.
   */
  bool element_is_continuous;

  /**
   * Stores the number of components in the finite element contained in the
   * DoFHandler passed to build().
   */
  unsigned int n_components;

  /**
   * Stores the number of degrees of freedom on all child cells. It is
   * <tt>2<sup>dim</sup>*fe.dofs_per_cell</tt> for DG elements and somewhat
   * less for continuous elements.
   */
  unsigned int n_child_cell_dofs;

  /**
   * Holds the indices for cells on a given level, extracted from DoFHandler
   * for fast access. All DoF indices on a given level are stored as a plain
   * array (since this class assumes constant DoFs per cell). To index into
   * this array, use the cell number times dofs_per_cell.
   *
   * This array first is arranged such that all locally owned level cells come
   * first (found in the variable n_owned_level_cells) and then other cells
   * necessary for the transfer to the next level.
   */
  std::vector<std::vector<unsigned int> > level_dof_indices;

  /**
   * Stores the connectivity from parent to child cell numbers for each level.
   */
  std::vector<std::vector<std::pair<unsigned int,unsigned int> > > parent_child_connect;

  /**
   * Stores the number of cells owned on a given process (sets the bounds for
   * the worker loops) for each level.
   */
  std::vector<unsigned int> n_owned_level_cells;

  /**
   * Holds the one-dimensional embedding (prolongation) matrix from mother
   * element to all the children.
   */
  AlignedVector<VectorizedArray<Number> > prolongation_matrix_1d;

  /**
   * Holds the temporary values for the tensor evaluation
   */
  mutable AlignedVector<VectorizedArray<Number> > evaluation_data;

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
  std::vector<AlignedVector<VectorizedArray<Number> > > weights_on_refined;

  /**
   * Stores the local indices of Dirichlet boundary conditions on cells for
   * all levels (outer index), the cells within the levels (second index), and
   * the indices on the cell (inner index).
   */
  std::vector<std::vector<std::vector<unsigned short> > > dirichlet_indices;

  /**
   * Performs templated prolongation operation
   */
  template <int degree>
  void do_prolongate_add(const unsigned int                           to_level,
                         LinearAlgebra::distributed::Vector<Number>       &dst,
                         const LinearAlgebra::distributed::Vector<Number> &src) const;

  /**
   * Performs templated restriction operation
   */
  template <int degree>
  void do_restrict_add(const unsigned int                           from_level,
                       LinearAlgebra::distributed::Vector<Number>       &dst,
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
 * MGTransferMatrixFree. This implies that each block should cover the
 * same index space of DoFs and have the same partitioning.
 *
 * @author Denis Davydov
 * @date 2017
 */
template <int dim, typename Number>
class MGTransferBlockMatrixFree : public MGTransferBase<LinearAlgebra::distributed::BlockVector<Number>>
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferBlockMatrixFree () = default;

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   */
  MGTransferBlockMatrixFree (const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Destructor.
   */
  virtual ~MGTransferBlockMatrixFree () = default;

  /**
   * Initialize the constraints to be used in build().
   */
  void initialize_constraints (const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void clear ();

  /**
   * Actually build the information for the prolongation for each level.
   */
  void build (const DoFHandler<dim,dim> &mg_dof);

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
  virtual void prolongate (const unsigned int                                    to_level,
                           LinearAlgebra::distributed::BlockVector<Number>       &dst,
                           const LinearAlgebra::distributed::BlockVector<Number> &src) const;

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
  virtual void restrict_and_add (const unsigned int from_level,
                                 LinearAlgebra::distributed::BlockVector<Number>       &dst,
                                 const LinearAlgebra::distributed::BlockVector<Number> &src) const;

  /**
   * Transfer from a block-vector on the global grid to block-vectors defined
   * on each of the levels separately.
   *
   * This function will initialize @p dst accordingly if needed as required by
   * the Multigrid class.
   */
  template <typename Number2, int spacedim>
  void
  copy_to_mg (const DoFHandler<dim,spacedim> &mg_dof,
              MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
              const LinearAlgebra::distributed::BlockVector<Number2>         &src) const;

  /**
   * Transfer from multi-level block-vector to normal vector.
   */
  template <typename Number2, int spacedim>
  void
  copy_from_mg (const DoFHandler<dim,spacedim>                                       &mg_dof,
                LinearAlgebra::distributed::BlockVector<Number2>                     &dst,
                const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src) const;

  /**
   * Memory used by this object.
   */
  std::size_t memory_consumption () const;

private:

  /**
   * A non-block matrix-free version of transfer operation.
   */
  MGTransferMatrixFree<dim,Number> matrix_free_transfer;
};


/*@}*/


//------------------------ templated functions -------------------------
#ifndef DOXYGEN

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim,Number>::
copy_to_mg (const DoFHandler<dim,spacedim> &mg_dof,
            MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &dst,
            const LinearAlgebra::distributed::BlockVector<Number2>         &src) const
{
  const unsigned int n_blocks  = src.n_blocks();
  const unsigned int min_level = dst.min_level();
  const unsigned int max_level = dst.max_level();

  // this function is normally called within the Multigrid class with
  // dst == defect level block vector. At first run this vector is not
  // initialized. Do this below:
  {
    const parallel::Triangulation<dim,spacedim> *tria =
      (dynamic_cast<const parallel::Triangulation<dim,spacedim>*>
       (&mg_dof.get_triangulation()));

    for (unsigned int level = min_level; level <= max_level; ++level)
      {
        dst[level].reinit(n_blocks);
        bool collect_size = false;
        for (unsigned int b = 0; b < n_blocks; ++b)
          {
            LinearAlgebra::distributed::Vector<Number> &v = dst[level].block(b);
            if (v.size() != mg_dof.locally_owned_mg_dofs(level).size() ||
                v.local_size() != mg_dof.locally_owned_mg_dofs(level).n_elements())
              {
                v.reinit(mg_dof.locally_owned_mg_dofs(level), tria != nullptr ?
                         tria->get_communicator() : MPI_COMM_SELF);
                collect_size = true;
              }
            else
              v = 0.;
          }
        if (collect_size)
          dst[level].collect_sizes ();
      }
  }

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> dst_non_block(min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        dst_non_block[l].reinit(dst[l].block(b));

      matrix_free_transfer.copy_to_mg(mg_dof, dst_non_block, src.block(b));

      for (unsigned int l = min_level; l <= max_level; ++l)
        dst[l].block(b) = dst_non_block[l];
    }
}

template <int dim, typename Number>
template <typename Number2, int spacedim>
void
MGTransferBlockMatrixFree<dim,Number>::
copy_from_mg (const DoFHandler<dim,spacedim>                        &mg_dof,
              LinearAlgebra::distributed::BlockVector<Number2>      &dst,
              const MGLevelObject<LinearAlgebra::distributed::BlockVector<Number>> &src) const
{
  const unsigned int n_blocks  = dst.n_blocks();
  const unsigned int min_level = src.min_level();
  const unsigned int max_level = src.max_level();

  for (unsigned int l = min_level; l <= max_level; ++l)
    AssertDimension(src[l].n_blocks(), dst.n_blocks());

  // FIXME: this a quite ugly as we need a temporary object:
  MGLevelObject<LinearAlgebra::distributed::Vector<Number>> src_non_block(min_level, max_level);

  for (unsigned int b = 0; b < n_blocks; ++b)
    {
      for (unsigned int l = min_level; l <= max_level; ++l)
        {
          src_non_block[l].reinit(src[l].block(b));
          src_non_block[l] = src[l].block(b);
        }

      matrix_free_transfer.copy_from_mg(mg_dof, dst.block(b), src_non_block);
    }
}



#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
