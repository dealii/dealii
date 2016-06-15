// ---------------------------------------------------------------------
//
// Copyright (C) 2016 - 2016 by the deal.II authors
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

#ifndef dealii__mg_transfer_matrix_free_h
#define dealii__mg_transfer_matrix_free_h

#include <deal.II/base/config.h>

#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/base/mg_level_object.h>
#include <deal.II/multigrid/mg_transfer.h>
#include <deal.II/matrix_free/shape_info.h>

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
  virtual ~MGTransferMatrixFree ();

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
   * element to the children.
   */
  internal::MatrixFreeFunctions::ShapeInfo<Number> shape_info;

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


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
