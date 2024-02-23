// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_matrix_tools_h
#define dealii_matrix_tools_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>

#include <deal.II/lac/affine_constraints.h>

#include <deal.II/numerics/matrix_creator.h>

#include <map>

#ifdef DEAL_II_WITH_PETSC
#  include <petscsys.h>
#endif

DEAL_II_NAMESPACE_OPEN


// forward declarations
#ifndef DOXYGEN
template <int dim>
class Quadrature;


template <typename number>
class Vector;
template <typename number>
class FullMatrix;
template <typename number>
class SparseMatrix;

template <typename number>
class BlockSparseMatrix;
template <typename Number>
class BlockVector;

template <int dim, int spacedim>
class Mapping;
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class DoFHandler;

namespace hp
{
  template <int>
  class QCollection;
  template <int, int>
  class MappingCollection;
} // namespace hp


#  ifdef DEAL_II_WITH_PETSC
namespace PETScWrappers
{
  class MatrixBase;
  class VectorBase;
  namespace MPI
  {
    class BlockSparseMatrix;
    class BlockVector;
  } // namespace MPI
} // namespace PETScWrappers
#  endif

#  ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class SparseMatrix;
  class BlockSparseMatrix;
  namespace MPI
  {
    class Vector;
    class BlockVector;
  } // namespace MPI
} // namespace TrilinosWrappers
#  endif
#endif



/**
 * Provide a collection of functions operating on matrices. These include the
 * application of boundary conditions to a linear system of equations and
 * others.
 *
 *
 * <h3>Boundary conditions</h3>
 *
 * The apply_boundary_values() functions modifies a linear system to incorporate
 * the constraints that result from Dirichlet-type boundary conditions (or, more
 * specifically: "strong" boundary conditions). To actually do this, the
 * functions of this name in the current namespace require a list of degree of
 * freedom indices along with the values these degrees of freedom should have.
 * To see how to get such a list, see the discussion of the
 * VectorTools::interpolate_boundary_values() function as one example.
 *
 * There are two ways to incorporate fixed degrees of freedom such as boundary
 * nodes into a linear system, as discussed below. Both operate at either the
 * level of local contributions to the global linear system, or the global
 * system itself. A third way, using
 * AffineConstraints::copy_local_to_global(), performs the same process as part
 * of adding the local contributions of one cell into the global linear system
 * (the "assembly" step) and is the method predominantly used in the tutorial
 * programs today.
 *
 * @dealiiVideoLecture{21.6,21.65}
 *
 *
 *
 * <h3>Global elimination</h3>
 *
 * In the first method, we first assemble the global linear system without
 * respect for fixed degrees of freedom, and in a second step eliminate them
 * again from the linear system. The inclusion into the assembly process is as
 * follows: when the matrix and vectors are set up, a list of nodes subject to
 * Dirichlet boundary conditions is made and matrix and vectors are
 * modified accordingly. This
 * is done by deleting all entries in the matrix in the line of this degree of
 * freedom, setting the main diagonal entry to a suitable positive value and
 * the right hand side element to a value so that the solution of the linear
 * system will have the boundary value at this node. To decouple the remaining
 * linear system of equations and to make the system symmetric again (at least
 * if it was before), one Gauss elimination step is performed with this line,
 * by adding this (now almost empty) line to all other lines which couple with
 * the given degree of freedom and thus eliminating all coupling between this
 * degree of freedom and others. Now the respective column also consists only
 * of zeroes, apart from the main diagonal entry. Alternatively, the functions
 * in this namespace take a boolean parameter that allows to omit this last
 * step, if symmetry of the resulting linear system is not required. Note that
 * usually even CG can cope with a non-symmetric linear system with this
 * particular structure.
 *
 * Finding which rows contain an entry in the column for which we are
 * presently performing a Gauss elimination step is either difficult or very
 * simple, depending on the circumstances. If the sparsity pattern is
 * symmetric (whether the matrix is symmetric is irrelevant here), then we can
 * infer the rows which have a nonzero entry in the present column by looking
 * at which columns in the present row are nonempty. In this case, we only
 * need to look into a fixed number of rows and need not search all rows. On
 * the other hand, if the sparsity pattern is nonsymmetric, then we need to
 * use an iterative solver which can handle nonsymmetric matrices in any case,
 * so there may be no need to do the Gauss elimination anyway. In fact, this
 * is the way the function works: it takes a parameter (@p eliminate_columns)
 * that specifies whether the sparsity pattern is symmetric; if so, then the
 * column is eliminated and the right hand side is also modified accordingly.
 * If not, then only the row is deleted and the column is not touched at all,
 * and all right hand side values apart from the one corresponding to the
 * present row remain unchanged.
 *
 * If the sparsity pattern for your matrix is non-symmetric, you must set the
 * value of this parameter to @p false in any case, since then we can't
 * eliminate the column without searching all rows, which would be too
 * expensive (if @p N be the number of rows, and @p m the number of nonzero
 * elements per row, then eliminating one column is an <tt>O(N*log(m))</tt>
 * operation, since searching in each row takes <tt>log(m)</tt> operations).
 * If your sparsity pattern is symmetric, but your matrix is not, then you
 * might specify @p false as well. If your sparsity pattern and matrix are
 * both symmetric, you might want to specify @p true (the complexity of
 * eliminating one row is then <tt>O(m*log(m))</tt>, since we only have to
 * search @p m rows for the respective element of the column). Given the fact
 * that @p m is roughly constant, irrespective of the discretization, and that
 * the number of boundary nodes is <tt>sqrt(N)</tt> in 2d, the algorithm for
 * symmetric sparsity patterns is <tt>O(sqrt(N)*m*log(m))</tt>, while it would
 * be <tt>O(N*sqrt(N)*log(m))</tt> for the general case; the latter is too
 * expensive to be performed.
 *
 * It seems as if we had to make clear not to overwrite the lines of other
 * boundary nodes when doing the Gauss elimination step. However, since we
 * reset the right hand side when passing such a node, it is not a problem to
 * change the right hand side values of other boundary nodes not yet
 * processed. It would be a problem to change those entries of nodes already
 * processed, but since the matrix entry of the present column on the row of
 * an already processed node is zero, the Gauss step does not change the right
 * hand side. We need therefore not take special care of other boundary nodes.
 *
 * To make solving faster, we preset the solution vector with the right
 * boundary values (as to why this is necessary, see the discussion below in
 * the description of local elimination). It it not clear whether the deletion
 * of coupling between the boundary degree of freedom and other dofs really
 * forces the corresponding entry in the solution vector to have the right
 * value when using iterative solvers, since their search directions may
 * contain components in the direction of the boundary node. For this reason,
 * we perform a very simple line balancing by not setting the main diagonal
 * entry to unity, but rather to the value it had before deleting this line,
 * or to the first nonzero main diagonal entry if it is zero for some reason.
 * Of course we have to change the right hand side appropriately. This is not
 * a very good strategy, but it at least should give the main diagonal entry a
 * value in the right order of dimension, which makes the solution process a
 * bit more stable. A refined algorithm would set the entry to the mean of the
 * other diagonal entries, but this seems to be too expensive.
 *
 * In some cases, it might be interesting to solve several times with the same
 * matrix, but for different right hand sides or boundary values. A typical
 * case would be the solution of a time-dependent problem in which the boundary
 * values or right hand side change, but the matrix itself does not. One
 * may then be tempted to just assemble the matrix once and just call the
 * MatrixTools::apply_boundary_values() function repeatedly on the same
 * matrix object, with a right hand side vector newly formed in each time step.
 * However,
 * since the modification for boundary values of the right hand side vector
 * depends on the original matrix, this is not possible without storing the
 * original matrix somewhere, and in every time step initializing the system
 * matrix with the unmodified matrix stored elsewhere. step-26 does a variation
 * of this process by storing building blocks from which the system matrix is
 * composed, but the general principle is the same. Alternatively, one can
 * use the constrained_linear_operator() function. In its documentation you can
 * also find a formal (mathematical) description of the process of modifying the
 * matrix and right hand side vectors for boundary values.
 *
 *
 * <h3>Local elimination</h3>
 *
 * The second way of handling boundary values is to modify the local matrix
 * and vector contributions appropriately before transferring them into the
 * global sparse matrix and vector. This is what local_apply_boundary_values()
 * does. The advantage is that we save the call to the apply_boundary_values
 * function (which is expensive because it has to work on sparse data
 * structures). On the other hand, the local_apply_boundary_values() function
 * is called many times, even if we only have a very small number of fixed
 * boundary nodes, and the main drawback is that this function doesn't work as
 * expected if there are hanging nodes that also need to be treated. The
 * reason that this function doesn't work is that it is meant to be run before
 * distribution into the global matrix, i.e. before hanging nodes are
 * distributed; since hanging nodes can be constrained to a boundary node, the
 * treatment of hanging nodes can add entries again to rows and columns
 * corresponding to boundary values and that we have already vacated in the
 * local elimination step. To make things worse, in 3d constrained nodes can
 * even lie on the boundary. Thus, it is imperative that boundary node
 * elimination happens @em after hanging node elimination, but this can't be
 * achieved with local elimination of boundary nodes unless there are no
 * hanging node constraints at all.
 *
 * Local elimination has one additional drawback: we don't have access to the
 * solution vector, only to the local contributions to the matrix and right
 * hand side. The problem with this is subtle, but can lead to very hard to
 * find difficulties: when we eliminate a degree of freedom, we delete the row
 * and column of this unknown, and set the diagonal entry to some positive
 * value. To make the problem more or less well-conditioned, we set this
 * diagonal entry to the absolute value of its prior value if that was
 * non-zero, or to the average magnitude of all other nonzero diagonal elements.
 * Then we set the right hand side value such that the resulting solution
 * entry has the right value as given by the boundary values. Since we add
 * these contributions up over all local contributions, the diagonal entry and
 * the respective value in the right hand side are added up correspondingly,
 * so that the entry in the solution of the linear system is still valid.
 *
 * A problem arises, however, if the diagonal entries so chosen are not
 * appropriate for the linear system. Consider, for example, a mixed Laplace
 * problem with matrix <tt>[[A B][C^T 0]]</tt>, where we only specify boundary
 * values for the second component of the solution. In the mixed formulation,
 * the stress-strain tensor only appears in either the matrix @p B or @p C, so
 * one of them may be significantly larger or smaller than the other one. Now,
 * if we eliminate boundary values, we delete some rows and columns, but we
 * also introduce a few entries on the diagonal of the lower right block, so
 * that we get the system <tt>[[A' B'][C'^T X]]</tt>. The diagonal entries in
 * the matrix @p X will be of the same order of magnitude as those in @p A.
 * Now, if we solve this system in the Schur complement formulation, we have
 * to invert the matrix <tt>X-C'^TA'^{-1}B'</tt>. Deleting rows and columns
 * above makes sure that boundary nodes indeed have empty rows and columns in
 * the Schur complement as well, except for the entries in @p X. However, the
 * entries in @p X may be of significantly different orders of magnitude than
 * those in <tt>C'^TA'^{-1}B'</tt>! If this is the case, we may run into
 * trouble with iterative solvers. For example, assume that we start with zero
 * entries in the solution vector and that the entries in @p X are several
 * orders of magnitude too small; in this case, iterative solvers will compute
 * the residual vector in each step and form correction vectors, but since the
 * entries in @p X are so small, the residual contributions for boundary nodes
 * are really small, despite the fact that the boundary nodes are still at
 * values close to zero and not in accordance with the prescribed boundary
 * values. Since the residual is so small, the corrections the iterative
 * solver computes are very small, and in the end the solver will indicate
 * convergence to a small total residual with the boundary values still being
 * significantly wrong.
 *
 * We avoid this problem in the global elimination process described above by
 * 'priming' the solution vector with the correct values for boundary nodes.
 * However, we can't do this for the local elimination process. Therefore, if
 * you experience a problem like the one above, you need to either increase
 * the diagonal entries in @p X to a size that matches those in the other part
 * of the Schur complement, or, simpler, prime the solution vector before you
 * start the solver.
 *
 * In conclusion, local elimination of boundary nodes only works if there are
 * no hanging nodes and even then doesn't always work fully satisfactorily.
 *
 * @ingroup numerics
 */
namespace MatrixTools
{
  /**
   * Import namespace MatrixCreator for backward compatibility with older
   * versions of deal.II in which these namespaces were classes and class
   * MatrixTools was publicly derived from class MatrixCreator.
   */
  using namespace MatrixCreator;

  /**
   * Apply Dirichlet boundary conditions to the system matrix and vectors as
   * described in the general documentation of this namespace.
   */
  template <typename number>
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    SparseMatrix<number>                            &matrix,
    Vector<number>                                  &solution,
    Vector<number>                                  &right_hand_side,
    const bool                                       eliminate_columns = true);

  /**
   * Apply Dirichlet boundary conditions to the system matrix and vectors as
   * described in the general documentation of this namespace. This function
   * works for block sparse matrices and block vectors.
   */
  template <typename number>
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    BlockSparseMatrix<number>                       &matrix,
    BlockVector<number>                             &solution,
    BlockVector<number>                             &right_hand_side,
    const bool                                       eliminate_columns = true);

#ifdef DEAL_II_WITH_PETSC
  /**
   * Apply Dirichlet boundary conditions to the system matrix and vectors as
   * described in the general documentation of this namespace. This function
   * works on the classes that are used to wrap PETSc objects.
   *
   * This function is used in step-17 and step-18.
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, PetscScalar> &boundary_values,
    PETScWrappers::MatrixBase                            &matrix,
    PETScWrappers::VectorBase                            &solution,
    PETScWrappers::VectorBase                            &right_hand_side,
    const bool eliminate_columns = true);

  /**
   * Same as above but for the parallel BlockSparseMatrix.
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, PetscScalar> &boundary_values,
    PETScWrappers::MPI::BlockSparseMatrix                &matrix,
    PETScWrappers::MPI::BlockVector                      &solution,
    PETScWrappers::MPI::BlockVector                      &right_hand_side,
    const bool eliminate_columns = true);

#endif

#ifdef DEAL_II_WITH_TRILINOS
  /**
   * Apply Dirichlet boundary conditions to the system matrix and vectors as
   * described in the general documentation of this namespace. This function
   * works on the classes that are used to wrap Trilinos objects.
   *
   * <b>Important:</b> This function is not very efficient: it needs to
   * alternatingly read and write into the matrix, a situation that Trilinos
   * does not handle well. In addition, we only get rid of rows corresponding
   * to boundary nodes, but the corresponding case of deleting the respective
   * columns (i.e. if @p eliminate_columns is @p true) is not presently
   * implemented, and probably will never because it is too expensive without
   * direct access to the Trilinos data structures. (This leads to the
   * situation where the action indicated by the default value of the last
   * argument is actually not implemented; that argument has <code>true</code>
   * as its default value to stay consistent with the other functions of same
   * name in this namespace.)
   *
   * @note If the matrix is stored in parallel across multiple processors
   * using MPI, this function only touches rows that are locally stored and
   * simply ignores all other rows. In other words, each processor is
   * responsible for its own rows, and the @p boundary_values argument needs
   * to contain all locally owned rows of the matrix that you want to have
   * treated. (But it can also contain entries for degrees of freedom not
   * owned locally; these will simply be ignored.) Further, in the context of
   * parallel computations, you will get into trouble if you treat a row while
   * other processors still have pending writes or additions into the same
   * row. In other words, if another processor still wants to add something to
   * an element of a row and you call this function to zero out the row, then
   * the next time you call compress() may add the remote value to the zero
   * you just created. Consequently, you will want to call compress() after
   * you made the last modifications to a matrix and before starting to clear
   * rows.
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, TrilinosScalar> &boundary_values,
    TrilinosWrappers::SparseMatrix                          &matrix,
    TrilinosWrappers::MPI::Vector                           &solution,
    TrilinosWrappers::MPI::Vector                           &right_hand_side,
    const bool eliminate_columns = true);

  /**
   * This function does the same as the one above, except now working on block
   * structures.
   */
  void
  apply_boundary_values(
    const std::map<types::global_dof_index, TrilinosScalar> &boundary_values,
    TrilinosWrappers::BlockSparseMatrix                     &matrix,
    TrilinosWrappers::MPI::BlockVector                      &solution,
    TrilinosWrappers::MPI::BlockVector                      &right_hand_side,
    const bool eliminate_columns = true);
#endif

  /**
   * Rather than applying boundary values to the global matrix and vector
   * after creating the global matrix, this function does so during assembly,
   * by modifying the local matrix and vector contributions. If you call this
   * function on all local contributions, the resulting matrix will have the
   * same entries, and the final call to apply_boundary_values() on the global
   * system will not be necessary.
   *
   * Since this function does not have to work on the complicated data
   * structures of sparse matrices, it is relatively cheap. It may therefore
   * be a win if you have many fixed degrees of freedom (e.g. boundary nodes),
   * or if access to the sparse matrix is expensive (e.g. for block sparse
   * matrices, or for PETSc or Trilinos matrices). However, it doesn't work as
   * expected if there are also hanging nodes to be considered. More caveats
   * are listed in the general documentation of this namespace.
   *
   * @dealiiVideoLecture{21.6,21.65}
   */
  template <typename number>
  void
  local_apply_boundary_values(
    const std::map<types::global_dof_index, number> &boundary_values,
    const std::vector<types::global_dof_index>      &local_dof_indices,
    FullMatrix<number>                              &local_matrix,
    Vector<number>                                  &local_rhs,
    const bool                                       eliminate_columns);

  /**
   * Exception
   */
  DeclExceptionMsg(ExcBlocksDontMatch,
                   "You are providing a matrix whose subdivision into "
                   "blocks in either row or column direction does not use "
                   "the same blocks sizes as the solution vector or "
                   "right hand side vectors, respectively.");
} // namespace MatrixTools



DEAL_II_NAMESPACE_CLOSE

#endif
