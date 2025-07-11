// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mg_transfer_component_h
#define dealii_mg_transfer_component_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/component_mask.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>

#include <memory>
#include <set>


DEAL_II_NAMESPACE_OPEN

// Forward declaration
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;
#endif

/*
 * MGTransferBase is defined in mg_base.h
 */

/**
 * @addtogroup mg
 * @{
 */

/**
 * Implementation of matrix generation for component wise multigrid transfer.
 *
 * @note MGTransferBlockBase is probably the more logical class. Still
 * eventually, a class should be developed allowing to select multiple
 * components.
 */
class MGTransferComponentBase
{
public:
  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;


protected:
  /**
   * Actually build the prolongation matrices for each level.
   *
   * This function is only called by derived classes. These can also set the
   * member variables <code>selected_component</code> and
   * <code>mg_selected_component</code> member variables to restrict the
   * transfer matrices to certain components. Furthermore, they use
   * <code>target_component</code> and <code>mg_target_component</code> for
   * re-ordering and grouping of components.
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof_handler);

  /**
   * Flag of selected components.
   *
   * The transfer operators only act on the components having a <tt>true</tt>
   * entry here. If renumbering by #target_component is used, this refers to
   * the <b>renumbered</b> components.
   */
  ComponentMask component_mask;

  /**
   * Flag of selected components.
   *
   * The transfer operators only act on the components having a <tt>true</tt>
   * entry here. If renumbering by #mg_target_component is used, this refers
   * to the <b>renumbered</b> components.
   */
  ComponentMask mg_component_mask;

  /**
   * Target component of the fine-level vector if renumbering is required.
   */
  std::vector<unsigned int> target_component;

  /**
   * Target component if renumbering of level vectors is required.
   */
  std::vector<unsigned int> mg_target_component;

  /**
   * Sizes of the multi-level vectors.
   */
  mutable std::vector<std::vector<types::global_dof_index>> sizes;

  /**
   * Start index of each component.
   */
  std::vector<types::global_dof_index> component_start;

  /**
   * Start index of each component on all levels.
   */
  std::vector<std::vector<types::global_dof_index>> mg_component_start;

  /**
   * Call build() function first.
   */
  DeclException0(ExcMatricesNotBuilt);

private:
  std::vector<std::shared_ptr<BlockSparsityPattern>> prolongation_sparsities;

protected:
  /**
   * The actual prolongation matrix. column indices belong to the dof indices
   * of the parent cell, i.e. the coarse level. while row indices belong to
   * the child cell, i.e. the fine level.
   */
  std::vector<std::shared_ptr<BlockSparseMatrix<double>>> prolongation_matrices;

  /**
   * This variable holds the mapping for the <tt>copy_to/from_mg</tt>-functions.
   * The data is first the global index, then the level index.
   */
  std::vector<std::vector<std::pair<types::global_dof_index, unsigned int>>>
    copy_to_and_from_indices;

  /**
   * Store the boundary_indices. These are needed for the boundary values in
   * the restriction matrix.
   */
  std::vector<std::set<types::global_dof_index>> boundary_indices;
};

// TODO:[GK] Update documentation for copy_* functions

// TODO: Use same kind of template argument as MGTransferSelect

/**
 * Implementation of the MGTransferBase interface for block matrices and
 * simple vectors. This class uses MGTransferComponentBase selecting a single
 * component or grouping several components into a single block. The transfer
 * operators themselves are implemented for Vector and BlockVector objects.
 *
 * See MGTransferBase to find out which of the transfer classes is best for
 * your needs.
 */
template <typename number>
class MGTransferSelect : public MGTransferBase<Vector<number>>,
                         private MGTransferComponentBase
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferSelect();

  /**
   * Constructor with constraint matrices.
   */
  MGTransferSelect(const AffineConstraints<double> &constraints);

  /**
   * Destructor.
   */
  virtual ~MGTransferSelect() override = default;

  // TODO: rewrite docs; make sure defaulted args are actually allowed
  /**
   * Actually build the prolongation matrices for grouped components.
   *
   * This function is a front-end for the same function in
   * MGTransferComponentBase.
   *
   * @arg selected Number of the block of the global vector to be copied from
   * and to the multilevel vector. This number refers to the renumbering by
   * <tt>target_component</tt>.
   *
   * @arg mg_selected Number of the block for which the transfer matrices
   * should be built.
   *
   * If <tt>mg_target_component</tt> is present, this refers to the renumbered
   * components.
   *
   * @arg target_component this argument allows grouping and renumbering of
   * components in the fine-level vector (see DoFRenumbering::component_wise).
   *
   * @arg mg_target_component this argument allows grouping and renumbering
   * of components in the level vectors (see DoFRenumbering::component_wise).
   * It also affects the behavior of the <tt>selected</tt> argument
   *
   * @arg boundary_indices holds the boundary indices on each level.
   */
  template <int dim, int spacedim>
  void
  build(const DoFHandler<dim, spacedim> &dof,
        unsigned int                     selected,
        unsigned int                     mg_selected,
        const std::vector<unsigned int> &target_component =
          std::vector<unsigned int>(),
        const std::vector<unsigned int> &mg_target_component =
          std::vector<unsigned int>(),
        const std::vector<std::set<types::global_dof_index>> &boundary_indices =
          std::vector<std::set<types::global_dof_index>>());

  /**
   * Change selected component. Handle with care!
   */
  void
  select(const unsigned int component,
         const unsigned int mg_component = numbers::invalid_unsigned_int);

  virtual void
  prolongate(const unsigned int    to_level,
             Vector<number>       &dst,
             const Vector<number> &src) const override;

  virtual void
  restrict_and_add(const unsigned int    from_level,
                   Vector<number>       &dst,
                   const Vector<number> &src) const override;

  /**
   * Transfer from a vector on the global grid to a multilevel vector for the
   * active degrees of freedom. In particular, for a globally refined mesh only
   * the finest level in @p dst is filled as a plain copy of @p src. All the
   * other level objects are left untouched.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
             MGLevelObject<Vector<number>>   &dst,
             const Vector<number2>           &src) const;

  /**
   * Transfer from multilevel vector to normal vector.
   *
   * Copies data from active portions of an multilevel vector into the
   * respective positions of a Vector.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim>     &mg_dof,
               Vector<number2>                     &dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * Add a multi-level vector to a normal vector.
   *
   * Works as the previous function, but probably not for continuous elements.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim>     &mg_dof,
                   Vector<number2>                     &dst,
                   const MGLevelObject<Vector<number>> &src) const;

  /**
   * Transfer from a vector on the global grid to a multilevel vector for the
   * active degrees of freedom. In particular, for a globally refined mesh only
   * the finest level in @p dst is filled as a plain copy of @p src. All the
   * other level objects are left untouched.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
             MGLevelObject<Vector<number>>   &dst,
             const BlockVector<number2>      &src) const;

  /**
   * Transfer from multilevel vector to normal vector.
   *
   * Copies data from active portions of a multilevel vector into the
   * respective positions of a global BlockVector.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim>     &mg_dof,
               BlockVector<number2>                &dst,
               const MGLevelObject<Vector<number>> &src) const;

  /**
   * Add a multi-level vector to a normal vector.
   *
   * Works as the previous function, but probably not for continuous elements.
   */
  template <int dim, typename number2, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim>     &mg_dof,
                   BlockVector<number2>                &dst,
                   const MGLevelObject<Vector<number>> &src) const;

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

private:
  /**
   * Implementation of the public function.
   */
  template <int dim, class OutVector, int spacedim>
  void
  do_copy_from_mg(const DoFHandler<dim, spacedim>     &mg_dof,
                  OutVector                           &dst,
                  const MGLevelObject<Vector<number>> &src) const;

  /**
   * Implementation of the public function.
   */
  template <int dim, class OutVector, int spacedim>
  void
  do_copy_from_mg_add(const DoFHandler<dim, spacedim>     &mg_dof,
                      OutVector                           &dst,
                      const MGLevelObject<Vector<number>> &src) const;

  /**
   * Actual implementation of copy_to_mg().
   */
  template <int dim, class InVector, int spacedim>
  void
  do_copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
                MGLevelObject<Vector<number>>   &dst,
                const InVector                  &src) const;
  /**
   * Selected component of global vector.
   */
  unsigned int selected_component;
  /**
   * Selected component inside multigrid.
   */
  unsigned int mg_selected_component;

  /**
   * The degrees of freedom on the refinement edges. For each level the index
   * set denotes which level degrees of freedom are on the refinement edge
   * towards the lower level, excluding boundary dofs.
   */
  std::vector<IndexSet> interface_dofs;

  /**
   * The constraints of the global system.
   */
public:
  ObserverPointer<const AffineConstraints<double>> constraints;
};

/** @} */

//---------------------------------------------------------------------------
template <typename number>
inline void
MGTransferSelect<number>::select(const unsigned int component,
                                 const unsigned int mg_component)
{
  selected_component = component;
  mg_selected_component =
    (mg_component == numbers::invalid_unsigned_int) ? component : mg_component;
}

DEAL_II_NAMESPACE_CLOSE

#endif
