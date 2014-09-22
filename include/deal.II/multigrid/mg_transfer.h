// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__mg_transfer_h
#define __deal2__mg_transfer_h

#include <deal.II/base/config.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/constraint_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparse_matrix.h>

#include <deal.II/lac/vector_memory.h>

#include <deal.II/multigrid/mg_base.h>
#include <deal.II/multigrid/mg_constrained_dofs.h>
#include <deal.II/base/mg_level_object.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>


DEAL_II_NAMESPACE_OPEN


template <int dim, int spacedim> class DoFHandler;

namespace internal
{
  template <class VECTOR>
  struct MatrixSelector
  {
    typedef ::dealii::SparsityPattern Sparsity;
    typedef ::dealii::SparseMatrix<typename VECTOR::value_type> Matrix;

    template <class CSP, class DH>
    static void reinit(Matrix &matrix, Sparsity &sparsity, int level, const CSP &csp, const DH &)
    {
      sparsity.copy_from (csp);
      matrix.reinit (sparsity);
    }
  };

#ifdef DEAL_II_WITH_TRILINOS
  template <>
  struct MatrixSelector<dealii::TrilinosWrappers::MPI::Vector>
  {
    typedef ::dealii::TrilinosWrappers::SparsityPattern Sparsity;
    typedef ::dealii::TrilinosWrappers::SparseMatrix Matrix;

    template <class CSP, class DH>
    static void reinit(Matrix &matrix, Sparsity &sparsity, int level, const CSP &csp, DH &dh)
    {
      matrix.reinit(dh.locally_owned_mg_dofs(level+1),
                    dh.locally_owned_mg_dofs(level),
                    csp, MPI_COMM_WORLD, true);
    }

  };

  template <>
  struct MatrixSelector<dealii::TrilinosWrappers::Vector>
  {
    typedef ::dealii::TrilinosWrappers::SparsityPattern Sparsity;
    typedef ::dealii::TrilinosWrappers::SparseMatrix Matrix;

    template <class CSP, class DH>
    static void reinit(Matrix &matrix, Sparsity &sparsity, int level, const CSP &csp, DH &dh)
    {
    }
  };
#endif
}

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
 * @author Wolfgang Bangerth, Guido Kanschat
 * @date 1999, 2000, 2001, 2002, 2003, 2004, 2012
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
   * Constructor with constraints. Equivalent to the default
   * constructor followed by initialize_constraints().
   */
  MGTransferPrebuilt (const ConstraintMatrix &constraints,
                      const MGConstrainedDoFs &mg_constrained_dofs);
  /**
   * Destructor.
   */
  virtual ~MGTransferPrebuilt ();

  /**
   * Initialize the constraints to be used in build_matrices().
   */
  void initialize_constraints (const ConstraintMatrix &constraints,
                               const MGConstrainedDoFs &mg_constrained_dofs);

  /**
   * Reset the object to the state it had right after the default constructor.
   */
  void clear ();

  /**
   * Actually build the prolongation
   * matrices for each level.
   */
  template <int dim, int spacedim>
  void build_matrices (const DoFHandler<dim,spacedim> &mg_dof);

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
  copy_to_mg (const DoFHandler<dim,spacedim> &mg_dof,
              MGLevelObject<VECTOR> &dst,
              const InVector &src) const;

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
  copy_from_mg (const DoFHandler<dim,spacedim> &mg_dof,
                OutVector &dst,
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
  copy_from_mg_add (const DoFHandler<dim,spacedim> &mg_dof,
                    OutVector &dst,
                    const MGLevelObject<VECTOR> &src) const;

  /**
   * If this object operates on
   * BlockVector objects, we need
   * to describe how the individual
   * vector components are mapped
   * to the blocks of a vector. For
   * example, for a Stokes system,
   * we have dim+1 vector
   * components for velocity and
   * pressure, but we may want to
   * use block vectors with only
   * two blocks for all velocities
   * in one block, and the pressure
   * variables in the other.
   *
   * By default, if this function
   * is not called, block vectors
   * have as many blocks as the
   * finite element has vector
   * components. However, this can
   * be changed by calling this
   * function with an array that
   * describes how vector
   * components are to be grouped
   * into blocks. The meaning of
   * the argument is the same as
   * the one given to the
   * DoFTools::count_dofs_per_component
   * function.
   */
  void
  set_component_to_block_map (const std::vector<unsigned int> &map);

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
  std::size_t memory_consumption () const;

  /**
   * Print all the matrices for debugging purposes.
   */
  void print_matrices(std::ostream &os) const;

  /**
   * Print the copy index fields for debugging purposes.
   */
  void print_indices(std::ostream &os) const;

private:

  /**
   * Sizes of the multi-level vectors.
   */
  std::vector<types::global_dof_index> sizes;

  /**
   * Sparsity patterns for transfer matrices.
   */
  std::vector<std_cxx11::shared_ptr<typename internal::MatrixSelector<VECTOR>::Sparsity> >   prolongation_sparsities;

  /**
   * The actual prolongation matrix.  column indices belong to the dof
   * indices of the mother cell, i.e. the coarse level.  while row
   * indices belong to the child cell, i.e. the fine level.
   */
  std::vector<std_cxx11::shared_ptr<typename internal::MatrixSelector<VECTOR>::Matrix> > prolongation_matrices;

  /**
   * Mapping for the copy_to_mg() and copy_from_mg() functions. Here only
   * index pairs locally owned
   *
   * The data is organized as follows: one vector per level. Each
   * element of these vectors contains first the global index, then
   * the level index.
   */
  std::vector<std::vector<std::pair<types::global_dof_index, unsigned int> > >
  copy_indices;

  /**
   * Additional degrees of freedom for the copy_to_mg()
   * function. These are the ones where the global degree of freedom
   * is locally owned and the level degree of freedom is not.
   *
   * Organization of the data is like for #copy_indices_mine.
   */
  std::vector<std::vector<std::pair<types::global_dof_index, unsigned int> > >
  copy_indices_to_me;

  /**
   * Additional degrees of freedom for the copy_from_mg()
   * function. These are the ones where the level degree of freedom
   * is locally owned and the global degree of freedom is not.
   *
   * Organization of the data is like for #copy_indices_mine.
   */
  std::vector<std::vector<std::pair<types::global_dof_index, unsigned int> > >
  copy_indices_from_me;


  /**
   * The vector that stores what
   * has been given to the
   * set_component_to_block_map()
   * function.
   */
  std::vector<unsigned int> component_to_block_map;

  /**
   * Degrees of freedom on the
   * refinement edge excluding
   * those on the boundary.
   */
  std::vector<std::vector<bool> > interface_dofs;
  /**
   * The constraints of the global
   * system.
   */
  SmartPointer<const ConstraintMatrix, MGTransferPrebuilt<VECTOR> > constraints;
  /**
   * The mg_constrained_dofs of the level
   * systems.
   */

  SmartPointer<const MGConstrainedDoFs, MGTransferPrebuilt<VECTOR> > mg_constrained_dofs;
};


/*@}*/


DEAL_II_NAMESPACE_CLOSE

#endif
