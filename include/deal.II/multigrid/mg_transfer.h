// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2019 by the deal.II authors
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

#ifndef dealii_mg_transfer_h
#define dealii_mg_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/la_parallel_vector.h>
#include <deal.II/lac/petsc_sparse_matrix.h>
#include <deal.II/lac/petsc_vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>
#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>

#include <memory>


DEAL_II_NAMESPACE_OPEN


namespace internal
{
  template <typename VectorType>
  struct MatrixSelector
  {
    using Sparsity = ::dealii::SparsityPattern;
    using Matrix   = ::dealii::SparseMatrix<typename VectorType::value_type>;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &                   matrix,
           Sparsity &                 sparsity,
           int                        level,
           const SparsityPatternType &sp,
           const DoFHandlerType &)
    {
      sparsity.copy_from(sp);
      (void)level;
      matrix.reinit(sparsity);
    }
  };

#ifdef DEAL_II_WITH_TRILINOS
  template <typename Number>
  struct MatrixSelector<LinearAlgebra::distributed::Vector<Number>>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                        level,
           const SparsityPatternType &sp,
           DoFHandlerType &           dh)
    {
      const parallel::TriangulationBase<DoFHandlerType::dimension,
                                        DoFHandlerType::space_dimension>
        *dist_tria = dynamic_cast<
          const parallel::TriangulationBase<DoFHandlerType::dimension,
                                            DoFHandlerType::space_dimension> *>(
          &(dh.get_triangulation()));
      MPI_Comm communicator =
        dist_tria != nullptr ? dist_tria->get_communicator() : MPI_COMM_SELF;

      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };

  template <>
  struct MatrixSelector<dealii::TrilinosWrappers::MPI::Vector>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                        level,
           const SparsityPatternType &sp,
           DoFHandlerType &           dh)
    {
      const parallel::TriangulationBase<DoFHandlerType::dimension,
                                        DoFHandlerType::space_dimension>
        *dist_tria = dynamic_cast<
          const parallel::TriangulationBase<DoFHandlerType::dimension,
                                            DoFHandlerType::space_dimension> *>(
          &(dh.get_triangulation()));
      MPI_Comm communicator =
        dist_tria != nullptr ? dist_tria->get_communicator() : MPI_COMM_SELF;
      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };

#  ifdef DEAL_II_WITH_MPI
#    ifdef DEAL_II_TRILINOS_WITH_TPETRA
  template <typename Number>
  struct MatrixSelector<dealii::LinearAlgebra::TpetraWrappers::Vector<Number>>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                        level,
           const SparsityPatternType &sp,
           DoFHandlerType &           dh)
    {
      const parallel::TriangulationBase<DoFHandlerType::dimension,
                                        DoFHandlerType::space_dimension>
        *dist_tria = dynamic_cast<
          const parallel::TriangulationBase<DoFHandlerType::dimension,
                                            DoFHandlerType::space_dimension> *>(
          &(dh.get_triangulation()));
      MPI_Comm communicator =
        dist_tria != nullptr ? dist_tria->get_communicator() : MPI_COMM_SELF;
      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };
#    endif

  template <>
  struct MatrixSelector<dealii::LinearAlgebra::EpetraWrappers::Vector>
  {
    using Sparsity = ::dealii::TrilinosWrappers::SparsityPattern;
    using Matrix   = ::dealii::TrilinosWrappers::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                        level,
           const SparsityPatternType &sp,
           DoFHandlerType &           dh)
    {
      const parallel::TriangulationBase<DoFHandlerType::dimension,
                                        DoFHandlerType::space_dimension>
        *dist_tria = dynamic_cast<
          const parallel::TriangulationBase<DoFHandlerType::dimension,
                                            DoFHandlerType::space_dimension> *>(
          &(dh.get_triangulation()));
      MPI_Comm communicator =
        dist_tria != nullptr ? dist_tria->get_communicator() : MPI_COMM_SELF;
      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator,
                    true);
    }
  };
#  endif

#else
  // ! DEAL_II_WITH_TRILINOS
  template <typename Number>
  struct MatrixSelector<LinearAlgebra::distributed::Vector<Number>>
  {
    using Sparsity = ::dealii::SparsityPattern;
    using Matrix   = ::dealii::SparseMatrix<Number>;

    static const bool requires_distributed_sparsity_pattern = false;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &,
           Sparsity &,
           int,
           const SparsityPatternType &,
           const DoFHandlerType &)
    {
      AssertThrow(
        false,
        ExcNotImplemented(
          "ERROR: MGTransferPrebuilt with LinearAlgebra::distributed::Vector currently "
          "needs deal.II to be configured with Trilinos."));
    }
  };

#endif

#ifdef DEAL_II_WITH_PETSC
  template <>
  struct MatrixSelector<dealii::PETScWrappers::MPI::Vector>
  {
    using Sparsity = ::dealii::DynamicSparsityPattern;
    using Matrix   = ::dealii::PETScWrappers::MPI::SparseMatrix;

    static const bool requires_distributed_sparsity_pattern = true;

    template <typename SparsityPatternType, typename DoFHandlerType>
    static void
    reinit(Matrix &matrix,
           Sparsity &,
           int                        level,
           const SparsityPatternType &sp,
           const DoFHandlerType &     dh)
    {
      const parallel::TriangulationBase<DoFHandlerType::dimension,
                                        DoFHandlerType::space_dimension>
        *dist_tria = dynamic_cast<
          const parallel::TriangulationBase<DoFHandlerType::dimension,
                                            DoFHandlerType::space_dimension> *>(
          &(dh.get_triangulation()));
      MPI_Comm communicator =
        dist_tria != nullptr ? dist_tria->get_communicator() : MPI_COMM_SELF;
      // Reinit PETSc matrix
      matrix.reinit(dh.locally_owned_mg_dofs(level + 1),
                    dh.locally_owned_mg_dofs(level),
                    sp,
                    communicator);
    }
  };
#endif
} // namespace internal

/*
 * MGTransferBase is defined in mg_base.h
 */

/*!@addtogroup mg */
/*@{*/



/**
 * Implementation of transfer between the global vectors and the multigrid
 * levels for use in the derived class MGTransferPrebuilt and other classes.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Timo Heister, Martin Kronbichler
 * @date 1999, 2000, 2001, 2002, 2003, 2004, 2012, 2015
 */
template <typename VectorType>
class MGLevelGlobalTransfer : public MGTransferBase<VectorType>
{
public:
  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void
  clear();

  /**
   * Transfer from a vector on the global grid to vectors defined on each of
   * the levels separately for the active degrees of freedom. In particular,
   * for a globally refined mesh only the finest level in @p dst is filled as a
   * plain copy of @p src. All the other level objects are left untouched.
   */
  template <int dim, class InVector, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &mg_dof,
             MGLevelObject<VectorType> &      dst,
             const InVector &                 src) const;

  /**
   * Transfer from multi-level vector to normal vector.
   *
   * Copies data from active portions of an MGVector into the respective
   * positions of a <tt>Vector<number></tt>. In order to keep the result
   * consistent, constrained degrees of freedom are set to zero.
   */
  template <int dim, class OutVector, int spacedim>
  void
  copy_from_mg(const DoFHandler<dim, spacedim> &mg_dof,
               OutVector &                      dst,
               const MGLevelObject<VectorType> &src) const;

  /**
   * Add a multi-level vector to a normal vector.
   *
   * Works as the previous function, but probably not for continuous elements.
   */
  template <int dim, class OutVector, int spacedim>
  void
  copy_from_mg_add(const DoFHandler<dim, spacedim> &mg_dof,
                   OutVector &                      dst,
                   const MGLevelObject<VectorType> &src) const;

  /**
   * If this object operates on BlockVector objects, we need to describe how
   * the individual vector components are mapped to the blocks of a vector.
   * For example, for a Stokes system, we have dim+1 vector components for
   * velocity and pressure, but we may want to use block vectors with only two
   * blocks for all velocities in one block, and the pressure variables in the
   * other.
   *
   * By default, if this function is not called, block vectors have as many
   * blocks as the finite element has vector components. However, this can be
   * changed by calling this function with an array that describes how vector
   * components are to be grouped into blocks. The meaning of the argument is
   * the same as the one given to the DoFTools::count_dofs_per_component
   * function.
   */
  void
  set_component_to_block_map(const std::vector<unsigned int> &map);

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Print the copy index fields for debugging purposes.
   */
  void
  print_indices(std::ostream &os) const;

protected:
  /**
   * Internal function to @p fill copy_indices*. Called by derived classes.
   */
  template <int dim, int spacedim>
  void
  fill_and_communicate_copy_indices(const DoFHandler<dim, spacedim> &mg_dof);

  /**
   * Sizes of the multi-level vectors.
   */
  std::vector<types::global_dof_index> sizes;

  /**
   * Mapping for the copy_to_mg() and copy_from_mg() functions. Here only
   * index pairs locally owned
   *
   * The data is organized as follows: one vector per level. Each element of
   * these vectors contains first the global index, then the level index.
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices;

  /**
   * Additional degrees of freedom for the copy_to_mg() function. These are
   * the ones where the global degree of freedom is locally owned and the
   * level degree of freedom is not.
   *
   * Organization of the data is like for @p copy_indices_mine.
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices_global_mine;

  /**
   * Additional degrees of freedom for the copy_from_mg() function. These are
   * the ones where the level degree of freedom is locally owned and the
   * global degree of freedom is not.
   *
   * Organization of the data is like for @p copy_indices_mine.
   */
  std::vector<
    std::vector<std::pair<types::global_dof_index, types::global_dof_index>>>
    copy_indices_level_mine;

  /**
   * This variable stores whether the copy operation from the global to the
   * level vector is actually a plain copy to the finest level. This means that
   * the grid has no adaptive refinement and the numbering on the finest
   * multigrid level is the same as in the global case.
   */
  bool perform_plain_copy;

  /**
   * The vector that stores what has been given to the
   * set_component_to_block_map() function.
   */
  std::vector<unsigned int> component_to_block_map;

  /**
   * The mg_constrained_dofs of the level systems.
   */
  SmartPointer<const MGConstrainedDoFs, MGLevelGlobalTransfer<VectorType>>
    mg_constrained_dofs;
};



/**
 * Implementation of transfer between the global vectors and the multigrid
 * levels for use in the derived class MGTransferPrebuilt and other classes.
 * This class is a specialization for the case of
 * LinearAlgebra::distributed::Vector that requires a few different calling
 * routines as compared to the %parallel vectors in the PETScWrappers and
 * TrilinosWrappers namespaces.
 *
 * @author Martin Kronbichler
 * @date 2016
 */
template <typename Number>
class MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>
  : public MGTransferBase<LinearAlgebra::distributed::Vector<Number>>
{
public:
  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void
  clear();

  /**
   * Transfer from a vector on the global grid to vectors defined on each of
   * the levels separately for the active degrees of freedom. In particular, for
   * a globally refined mesh only the finest level in @p dst is filled as a
   * plain copy of @p src. All the other level objects are left untouched.
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &                          mg_dof,
             MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
             const LinearAlgebra::distributed::Vector<Number2> &src) const;

  /**
   * Transfer from multi-level vector to normal vector.
   *
   * Copies data from active portions of an MGVector into the respective
   * positions of a <tt>Vector<number></tt>. In order to keep the result
   * consistent, constrained degrees of freedom are set to zero.
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_from_mg(
    const DoFHandler<dim, spacedim> &                                mg_dof,
    LinearAlgebra::distributed::Vector<Number2> &                    dst,
    const MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &src) const;

  /**
   * Add a multi-level vector to a normal vector.
   *
   * Works as the previous function, but probably not for continuous elements.
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_from_mg_add(
    const DoFHandler<dim, spacedim> &                                mg_dof,
    LinearAlgebra::distributed::Vector<Number2> &                    dst,
    const MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &src) const;

  /**
   * If this object operates on BlockVector objects, we need to describe how
   * the individual vector components are mapped to the blocks of a vector.
   * For example, for a Stokes system, we have dim+1 vector components for
   * velocity and pressure, but we may want to use block vectors with only two
   * blocks for all velocities in one block, and the pressure variables in the
   * other.
   *
   * By default, if this function is not called, block vectors have as many
   * blocks as the finite element has vector components. However, this can be
   * changed by calling this function with an array that describes how vector
   * components are to be grouped into blocks. The meaning of the argument is
   * the same as the one given to the DoFTools::count_dofs_per_component
   * function.
   */
  void
  set_component_to_block_map(const std::vector<unsigned int> &map);

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Print the copy index fields for debugging purposes.
   */
  void
  print_indices(std::ostream &os) const;

protected:
  /**
   * Internal function to perform transfer of residuals or solutions
   * basesd on the flag @p solution_transfer.
   */
  template <int dim, typename Number2, int spacedim>
  void
  copy_to_mg(const DoFHandler<dim, spacedim> &                          mg_dof,
             MGLevelObject<LinearAlgebra::distributed::Vector<Number>> &dst,
             const LinearAlgebra::distributed::Vector<Number2> &        src,
             const bool solution_transfer) const;

  /**
   * Internal function to @p fill copy_indices*. Called by derived classes.
   */
  template <int dim, int spacedim>
  void
  fill_and_communicate_copy_indices(const DoFHandler<dim, spacedim> &mg_dof);

  /**
   * Sizes of the multi-level vectors.
   */
  std::vector<types::global_dof_index> sizes;

  /**
   * Mapping for the copy_to_mg() and copy_from_mg() functions. Here only
   * index pairs locally owned is stored.
   *
   * The data is organized as follows: one table per level. This table has two
   * rows. The first row contains the global index, the second one the level
   * index.
   */
  std::vector<Table<2, unsigned int>> copy_indices;

  /**
   * Same as above, but used to transfer solution vectors.
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices;

  /**
   * Additional degrees of freedom for the copy_to_mg() function. These are
   * the ones where the global degree of freedom is locally owned and the
   * level degree of freedom is not.
   *
   * Organization of the data is like for @p copy_indices.
   */
  std::vector<Table<2, unsigned int>> copy_indices_global_mine;

  /**
   * Same as above, but used to transfer solution vectors.
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices_global_mine;

  /**
   * Additional degrees of freedom for the copy_from_mg() function. These are
   * the ones where the level degree of freedom is locally owned and the
   * global degree of freedom is not.
   *
   * Organization of the data is like for @p copy_indices.
   */
  std::vector<Table<2, unsigned int>> copy_indices_level_mine;

  /**
   * Same as above, but used to transfer solution vectors.
   */
  std::vector<Table<2, unsigned int>> solution_copy_indices_level_mine;

  /**
   * This variable stores whether the copy operation from the global to the
   * level vector is actually a plain copy to the finest level. This means that
   * the grid has no adaptive refinement and the numbering on the finest
   * multigrid level is the same as in the global case.
   */
  bool perform_plain_copy;

  /**
   * This variable stores whether the copy operation from the global to the
   * level vector is actually a plain copy to the finest level except for a
   * renumbering within the finest level of the degrees of freedom. This means
   * that the grid has no adaptive refinement.
   */
  bool perform_renumbered_plain_copy;

  /**
   * The vector that stores what has been given to the
   * set_component_to_block_map() function.
   */
  std::vector<unsigned int> component_to_block_map;

  /**
   * The mg_constrained_dofs of the level systems.
   */
  SmartPointer<
    const MGConstrainedDoFs,
    MGLevelGlobalTransfer<LinearAlgebra::distributed::Vector<Number>>>
    mg_constrained_dofs;

  /**
   * In the function copy_to_mg, we need to access ghosted entries of the
   * global vector for inserting into the level vectors. This vector is
   * populated with those entries.
   */
  mutable LinearAlgebra::distributed::Vector<Number> ghosted_global_vector;

  /**
   * Same as above but used when working with solution vectors.
   */
  mutable LinearAlgebra::distributed::Vector<Number>
    solution_ghosted_global_vector;

  /**
   * In the function copy_from_mg, we access all level vectors with certain
   * ghost entries for inserting the result into a global vector.
   */
  mutable MGLevelObject<LinearAlgebra::distributed::Vector<Number>>
    ghosted_level_vector;

  /**
   * Same as above but used when working with solution vectors.
   */
  mutable MGLevelObject<LinearAlgebra::distributed::Vector<Number>>
    solution_ghosted_level_vector;
};



/**
 * Implementation of the MGTransferBase interface for which the transfer
 * operations are prebuilt upon construction of the object of this class as
 * matrices. This is the fast way, since it only needs to build the operation
 * once by looping over all cells and storing the result in a matrix for each
 * level, but requires additional memory.
 *
 * See MGTransferBase to find out which of the transfer classes is best for
 * your needs.
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Timo Heister, Martin Kronbichler
 * @date 1999, 2000, 2001, 2002, 2003, 2004, 2012, 2015
 */
template <typename VectorType>
class MGTransferPrebuilt : public MGLevelGlobalTransfer<VectorType>
{
public:
  /**
   * Constructor without constraint matrices. Use this constructor only with
   * discontinuous finite elements or with no local refinement.
   */
  MGTransferPrebuilt() = default;

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   */
  MGTransferPrebuilt(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Constructor with constraints. Equivalent to the default constructor
   * followed by initialize_constraints().
   *
   * @deprecated @p constraints is unused.
   */
  DEAL_II_DEPRECATED
  MGTransferPrebuilt(const AffineConstraints<double> &constraints,
                     const MGConstrainedDoFs &        mg_constrained_dofs);

  /**
   * Destructor.
   */
  virtual ~MGTransferPrebuilt() override = default;

  /**
   * Initialize the constraints to be used in build_matrices().
   */
  void
  initialize_constraints(const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Initialize the constraints to be used in build_matrices().
   *
   * @deprecated @p constraints is unused.
   */
  DEAL_II_DEPRECATED
  void
  initialize_constraints(const AffineConstraints<double> &constraints,
                         const MGConstrainedDoFs &        mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void
  clear();

  /**
   * Actually build the prolongation matrices for each level.
   */
  template <int dim, int spacedim>
  void
  build_matrices(const DoFHandler<dim, spacedim> &mg_dof);

  /**
   * Prolongate a vector from level <tt>to_level-1</tt> to level
   * <tt>to_level</tt> using the embedding matrices of the underlying finite
   * element. The previous content of <tt>dst</tt> is overwritten.
   *
   * @arg src is a vector with as many elements as there are degrees of
   * freedom on the coarser level involved.
   *
   * @arg dst has as many elements as there are degrees of freedom on the
   * finer level.
   */
  virtual void
  prolongate(const unsigned int to_level,
             VectorType &       dst,
             const VectorType & src) const override;

  /**
   * Restrict a vector from level <tt>from_level</tt> to level
   * <tt>from_level-1</tt> using the transpose operation of the @p prolongate
   * method. If the region covered by cells on level <tt>from_level</tt> is
   * smaller than that of level <tt>from_level-1</tt> (local refinement), then
   * some degrees of freedom in <tt>dst</tt> are active and will not be
   * altered. For the other degrees of freedom, the result of the restriction
   * is added.
   *
   * @arg src is a vector with as many elements as there are degrees of
   * freedom on the finer level involved.
   *
   * @arg dst has as many elements as there are degrees of freedom on the
   * coarser level.
   */
  virtual void
  restrict_and_add(const unsigned int from_level,
                   VectorType &       dst,
                   const VectorType & src) const override;

  /**
   * Finite element does not provide prolongation matrices.
   */
  DeclException0(ExcNoProlongation);

  /**
   * You have to call build_matrices() before using this object.
   */
  DeclException0(ExcMatricesNotBuilt);

  /**
   * Memory used by this object.
   */
  std::size_t
  memory_consumption() const;

  /**
   * Print all the matrices for debugging purposes.
   */
  void
  print_matrices(std::ostream &os) const;

private:
  /**
   * Sparsity patterns for transfer matrices.
   */
  std::vector<
    std::shared_ptr<typename internal::MatrixSelector<VectorType>::Sparsity>>
    prolongation_sparsities;

  /**
   * The actual prolongation matrix.  column indices belong to the dof indices
   * of the mother cell, i.e. the coarse level.  while row indices belong to
   * the child cell, i.e. the fine level.
   */
  std::vector<
    std::shared_ptr<typename internal::MatrixSelector<VectorType>::Matrix>>
    prolongation_matrices;

  /**
   * Degrees of freedom on the refinement edge excluding those on the
   * boundary.
   */
  std::vector<std::vector<bool>> interface_dofs;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
