// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__sparse_vanka_h
#define __deal2__sparse_vanka_h



#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/multithread_info.h>

#include <vector>
#include <map>

DEAL_II_NAMESPACE_OPEN

template <typename number> class FullMatrix;
template <typename number> class SparseMatrix;
template <typename number> class Vector;

template <typename number> class SparseVanka;
template <typename number> class SparseBlockVanka;

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Point-wise Vanka preconditioning.
 * This class does Vanka preconditioning  on a point-wise base.
 * Vanka preconditioners are used for saddle point problems like Stoke's
 * problem or problems arising in optimization where Lagrange multiplier
 * occur and let Netwon's matrix have a zero block. With these matrices the
 * application of Jacobi or Gauss-Seidel methods is impossible, because
 * some diagonal elements are zero in the rows of the Lagrange multiplier.
 * The approach of Vanka is to solve a small (usually indefinite) system
 * of equations for each Langrange multiplie variable (we will also call
 * the pressure in Stokes' equation a Langrange multiplier since it
 * can be interpreted as such).
 *
 * Objects of this class are constructed by passing a vector of indices
 * of the degrees of freedom of the Lagrange multiplier. In the actual
 * preconditioning method, these rows are traversed in the order in which
 * the appear in the matrix. Since this is a Gau�-Seidel like procedure,
 * remember to have a good ordering in advance (for transport dominated
 * problems, Cuthill-McKee algorithms are a good means for this, if points
 * on the inflow boundary are chosen as starting points for the renumbering).
 *
 * For each selected degree of freedom, a local system of equations is built
 * by the degree of freedom itself and all other values coupling immediately,
 * i.e. the set of degrees of freedom considered for the local system of
 * equations for degree of freedom @p i is @p i itself and all @p j such that
 * the element <tt>(i,j)</tt> is a nonzero entry in the sparse matrix under
 * consideration. The elements <tt>(j,i)</tt> are not considered. We now pick all
 * matrix entries from rows and columns out of the set of degrees of freedom
 * just described out of the global matrix and put it into a local matrix,
 * which is subsequently inverted. This system may be of different size for
 * each degree of freedom, depending for example on the local neighborhood of
 * the respective node on a computational grid.
 *
 * The right hand side is built up in the same way, i.e. by copying
 * all entries that coupled with the one under present consideration,
 * but it is augmented by all degrees of freedom coupling with the
 * degrees from the set described above (i.e. the DoFs coupling second
 * order to the present one). The reason for this is, that the local
 * problems to be solved shall have Dirichlet boundary conditions on
 * the second order coupling DoFs, so we have to take them into
 * account but eliminate them before actually solving; this
 * elimination is done by the modification of the right hand side, and
 * in the end these degrees of freedom do not occur in the matrix and
 * solution vector any more at all.
 *
 * This local system is solved and the values are updated into the
 * destination vector.
 *
 * Remark: the Vanka method is a non-symmetric preconditioning method.
 *
 *
 * <h3>Example of Use</h3>
 * This little example is taken from a program doing parameter optimization.
 * The Lagrange multiplier is the third component of the finite element
 * used. The system is solved by the GMRES method.
 * @code
 *                        // tag the Lagrange multiplier variable
 *    vector<bool> signature(3);
 *    signature[0] = signature[1] = false;
 *    signature[2] = true;
 *
 *                        // tag all dofs belonging to the
 *                        // Lagrange multiplier
 *    vector<bool> selected_dofs (dof.n_dofs(), false);
 *    DoFTools::extract_dofs(dof, signature, p_select);
 *                        // create the Vanka object
 *    SparseVanka<double> vanka (global_matrix, selected_dofs);
 *
 *                        // create the solver
 *    SolverGMRES<PreconditionedSparseMatrix<double>,
 *                Vector<double> >    gmres(control,memory,504);
 *
 *                        // solve
 *    gmres.solve (global_matrix, solution, right_hand_side,
 *                 vanka);
 * @endcode
 *
 *
 * <h4>Implementor's remark</h4>
 * At present, the local matrices are built up such that the degree of freedom
 * associated with the local Lagrange multiplier is the first one. Thus, usually
 * the upper left entry in the local matrix is zero. It is not clear to me (W.B.)
 * whether this might pose some problems in the inversion of the local matrices.
 * Maybe someone would like to check this.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @author Guido Kanschat, Wolfgang Bangerth; 1999, 2000
 */
template<typename number>
class SparseVanka
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor. Gets the matrix
   * for preconditioning and a bit
   * vector with entries @p true for
   * all rows to be updated. A
   * reference to this vector will
   * be stored, so it must persist
   * longer than the Vanka
   * object. The same is true for
   * the matrix.
   *
   * The matrix @p M which is passed
   * here may or may not be the
   * same matrix for which this
   * object shall act as
   * preconditioner. In particular,
   * it is conceivable that the
   * preconditioner is build up for
   * one matrix once, but is used
   * for subsequent steps in a
   * nonlinear process as well,
   * where the matrix changes in
   * each step slightly.
   *
   * If @p conserve_mem is @p false,
   * then the inverses of the local
   * systems are computed now; if
   * the flag is @p true, then they
   * are computed every time the
   * preconditioner is
   * applied. This saves some
   * memory, but makes
   * preconditioning very
   * slow. Note also, that if the
   * flag is @p false, then the
   * contents of the matrix @p M at
   * the time of calling this
   * constructor are used, while if
   * the flag is @p true, then the
   * values in @p M at the time of
   * preconditioning are used. This
   * may lead to different results,
   * obviously, of @p M changes.
   *
   * The parameter @p n_threads
   * determines how many threads
   * shall be used in parallel when
   * building the inverses of the
   * diagonal blocks. This
   * parameter is ignored if not in
   * multithreaded mode.
   */
  SparseVanka(const SparseMatrix<number> &M,
              const std::vector<bool>    &selected,
              const bool                  conserve_memory = false,
              const unsigned int          n_threads       = multithread_info.n_threads());

  /**
   * Destructor.
   * Delete all allocated matrices.
   */
  ~SparseVanka();

  /**
   * Do the preconditioning.
   * This function takes the residual
   * in @p src and returns the resulting
   * update vector in @p dst.
   */
  template<typename number2>
  void vmult (Vector<number2>       &dst,
              const Vector<number2> &src) const;

protected:
  /**
   * Apply the inverses
   * corresponding to those degrees
   * of freedom that have a @p true
   * value in @p dof_mask, to the
   * @p src vector and move the
   * result into @p dst. Actually,
   * only values for allowed
   * indices are written to @p dst,
   * so the application of this
   * function only does what is
   * announced in the general
   * documentation if the given
   * mask sets all values to zero
   *
   * The reason for providing the
   * mask anyway is that in derived
   * classes we may want to apply
   * the preconditioner to parts of
   * the matrix only, in order to
   * parallelize the
   * application. Then, it is
   * important to only write to
   * some slices of @p dst, in order
   * to eliminate the dependencies
   * of threads of each other.
   *
   * If a null pointer is passed
   * instead of a pointer to the
   * @p dof_mask (as is the default
   * value), then it is assumed
   * that we shall work on all
   * degrees of freedom. This is
   * then equivalent to calling the
   * function with a
   * <tt>vector<bool>(n_dofs,true)</tt>.
   *
   * The @p vmult of this class
   * of course calls this function
   * with a null pointer
   */
  template<typename number2>
  void apply_preconditioner (Vector<number2>         &dst,
                             const Vector<number2>   &src,
                             const std::vector<bool> *const dof_mask = 0) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Pointer to the matrix.
   */
  SmartPointer<const SparseMatrix<number>,SparseVanka<number> > matrix;

  /**
   * Conserve memory flag.
   */
  const bool conserve_mem;

  /**
   * Indices of those degrees of
   * freedom that we shall work on.
   */
  const std::vector<bool> &selected;

  /**
   * Number of threads to be used
   * when building the
   * inverses. Only relevant in
   * multithreaded mode.
   */
  const unsigned int n_threads;

  /**
   * Array of inverse matrices,
   * one for each degree of freedom.
   * Only those elements will be used
   * that are tagged in @p selected.
   */
  mutable std::vector<SmartPointer<FullMatrix<float>,SparseVanka<number> > > inverses;

  /**
   * Compute the inverses of all
   * selected diagonal elements.
   */
  void compute_inverses ();

  /**
   * Compute the inverses at
   * positions in the range
   * <tt>[begin,end)</tt>. In
   * non-multithreaded mode,
   * <tt>compute_inverses()</tt> calls
   * this function with the whole
   * range, but in multithreaded
   * mode, several copies of this
   * function are spawned.
   */
  void compute_inverses (const size_type begin,
                         const size_type end);

  /**
   * Compute the inverse of the
   * block located at position
   * @p row. Since the vector is
   * used quite often, it is
   * generated only once in the
   * caller of this function and
   * passed to this function which
   * first clears it. Reusing the
   * vector makes the process
   * significantly faster than in
   * the case where this function
   * re-creates it each time.
   */
  void compute_inverse (const size_type         row,
                        std::vector<size_type> &local_indices);

  /**
   * Make the derived class a
   * friend. This seems silly, but
   * is actually necessary, since
   * derived classes can only
   * access non-public members
   * through their @p this
   * pointer, but not access these
   * members as member functions of
   * other objects of the type of
   * this base class (i.e. like
   * <tt>x.f()</tt>, where @p x is an
   * object of the base class, and
   * @p f one of it's non-public
   * member functions).
   *
   * Now, in the case of the
   * @p SparseBlockVanka class, we
   * would like to take the address
   * of a function of the base
   * class in order to call it
   * through the multithreading
   * framework, so the derived
   * class has to be a friend.
   */
  template <typename T> friend class SparseBlockVanka;
};



/**
 * Block version of the sparse Vanka preconditioner. This class
 * divides the matrix into blocks and works on the diagonal blocks
 * only, which of course reduces the efficiency as preconditioner, but
 * is perfectly parallelizable. The constructor takes a parameter into
 * how many blocks the matrix shall be subdivided and then lets the
 * underlying class do the work. Division of the matrix is done in
 * several ways which are described in detail below.
 *
 * This class is probably useless if you don't have a multiprocessor
 * system, since then the amount of work per preconditioning step is
 * the same as for the @p SparseVanka class, but preconditioning
 * properties are worse. On the other hand, if you have a
 * multiprocessor system, the worse preconditioning quality (leading
 * to more iterations of the linear solver) usually is well balanced
 * by the increased speed of application due to the parallelization,
 * leading to an overall decrease in elapsed wall-time for solving
 * your linear system. It should be noted that the quality as
 * preconditioner reduces with growing number of blocks, so there may
 * be an optimal value (in terms of wall-time per linear solve) for
 * the number of blocks.
 *
 * To facilitate writing portable code, if the number of blocks into
 * which the matrix is to be subdivided, is set to one, then this
 * class acts just like the @p SparseVanka class. You may therefore
 * want to set the number of blocks equal to the number of processors
 * you have.
 *
 * Note that the parallelization is done if <tt>deal.II</tt> was configured
 * for multithread use and that the number of threads which is spawned
 * equals the number of blocks. This is reasonable since you will not
 * want to set the number of blocks unnecessarily large, since, as
 * mentioned, this reduces the preconditioning properties.
 *
 *
 * <h3>Splitting the matrix into blocks</h3>
 *
 * Splitting the matrix into blocks is always done in a way such that
 * the blocks are not necessarily of equal size, but such that the
 * number of selected degrees of freedom for which a local system is
 * to be solved is equal between blocks. The reason for this strategy
 * to subdivision is load-balancing for multithreading. There are
 * several possibilities to actually split the matrix into blocks,
 * which are selected by the flag @p blocking_strategy that is passed
 * to the constructor. By a block, we will in the sequel denote a list
 * of indices of degrees of freedom; the algorithm will work on each
 * block separately, i.e. the solutions of the local systems
 * corresponding to a degree of freedom of one block will only be used
 * to update the degrees of freedom belonging to the same block, but
 * never to update degrees of freedoms of other blocks. A block can be
 * a consecutive list of indices, as in the first alternative below,
 * or a nonconsecutive list of indices. Of course, we assume that the
 * intersection of each two blocks is empty and that the union of all
 * blocks equals the interval <tt>[0,N)</tt>, where @p N is the number of
 * degrees of freedom of the system of equations.
 *
 * <ul>
 * <li> @p index_intervals:
 *    Here, we chose the blocks to be intervals <tt>[a_i,a_{i+1</tt>)},
 *    i.e. consecutive degrees of freedom are usually also within the
 *    same block. This is a reasonable strategy, if the degrees of
 *    freedom have, for example, be re-numbered using the
 *    Cuthill-McKee algorithm, in which spatially neighboring degrees
 *    of freedom have neighboring indices. In that case, coupling in
 *    the matrix is usually restricted to the vicinity of the diagonal
 *    as well, and we can simply cut the matrix into blocks.
 *
 *    The bounds of the intervals, i.e. the @p a_i above, are chosen
 *    such that the number of degrees of freedom on which we shall
 *    work (i.e. usually the degrees of freedom corresponding to
 *    Lagrange multipliers) is about the same in each block; this does
 *    not mean, however, that the sizes of the blocks are equal, since
 *    the blocks also comprise the other degrees of freedom for which
 *    no local system is solved. In the extreme case, consider that
 *    all Lagrange multipliers are sorted to the end of the range of
 *    DoF indices, then the first block would be very large, since it
 *    comprises all other DoFs and some Lagrange multipliers, while
 *    all other blocks are rather small and comprise only Langrange
 *    multipliers. This strategy therefore does not only depend on the
 *    order in which the Lagrange DoFs are sorted, but also on the
 *    order in which the other DoFs are sorted. It is therefore
 *    necessary to note that this almost renders the capability as
 *    preconditioner useless if the degrees of freedom are numbered by
 *    component, i.e. all Lagrange multipliers en bloc.
 *
 * <li> @p adaptive: This strategy is a bit more clever in cases where
 *    the Langrange DoFs are clustered, as in the example above. It
 *    works as follows: it first groups the Lagrange DoFs into blocks,
 *    using the same strategy as above. However, instead of grouping
 *    the other DoFs into the blocks of Lagrange DoFs with nearest DoF
 *    index, it decides for each non-Lagrange DoF to put it into the
 *    block of Lagrange DoFs which write to this non-Lagrange DoF most
 *    often. This makes it possible to even sort the Lagrange DoFs to
 *    the end and still associate spatially neighboring non-Lagrange
 *    DoFs to the same blocks where the respective Lagrange DoFs are,
 *    since they couple to each other while spatially distant DoFs
 *    don't couple.
 *
 *    The additional computational effort to sorting the non-Lagrange
 *    DoFs is not very large compared with the inversion of the local
 *    systems and applying the preconditioner, so this strategy might
 *    be reasonable if you want to sort your degrees of freedom by
 *    component. If the degrees of freedom are not sorted by
 *    component, the results of the both strategies outlined above
 *    does not differ much. However, unlike the first strategy, the
 *    performance of the second strategy does not deteriorate if the
 *    DoFs are renumbered by component.
 * </ul>
 *
 *
 * <h3>Typical results</h3>
 *
 * As a prototypical test case, we use a nonlinear problem from
 * optimization, which leads to a series of saddle point problems,
 * each of which is solved using GMRES with Vanka as
 * preconditioner. The equation had approx. 850 degrees of
 * freedom. With the non-blocked version @p SparseVanka (or
 * @p SparseBlockVanka with <tt>n_blocks==1</tt>), the following numbers of
 * iterations is needed to solver the linear system in each nonlinear
 * step:
 * @verbatim
 *   101 68 64 53 35 21
 * @endverbatim
 *
 * With four blocks, we need the following numbers of iterations
 * @verbatim
 *   124 88 83 66 44 28
 * @endverbatim
 * As can be seen, more iterations are needed. However, in terms of
 * computing time, the first version needs 72 seconds wall time (and
 * 79 seconds CPU time, which is more than wall time since some other
 * parts of the program were parallelized as well), while the second
 * version needed 53 second wall time (and 110 seconds CPU time) on a
 * four processor machine. The total time is in both cases dominated
 * by the linear solvers. In this case, it is therefore worth while
 * using the blocked version of the preconditioner if wall time is
 * more important than CPU time.
 *
 * The results with the block version above were obtained with the
 * first blocking strategy and the degrees of freedom were not
 * numbered by component. Using the second strategy does not much
 * change the numbers of iterations (at most by one in each step) and
 * they also do not change when the degrees of freedom are sorted
 * by component, while the first strategy significantly deteriorated.
 *
 * @author Wolfgang Bangerth, 2000
 */
template<typename number>
class SparseBlockVanka : public SparseVanka<number>
{
public:
  /**
   * Declate type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Enumeration of the different
   * methods by which the DoFs are
   * distributed to the blocks on
   * which we are to work.
   */
  enum BlockingStrategy
  {
    index_intervals, adaptive
  };

  /**
   * Constructor. Pass all
   * arguments except for
   * @p n_blocks to the base class.
   */
  SparseBlockVanka (const SparseMatrix<number> &M,
                    const std::vector<bool>    &selected,
                    const unsigned int          n_blocks,
                    const BlockingStrategy      blocking_strategy,
                    const bool                  conserve_memory = false,
                    const unsigned int          n_threads       = multithread_info.n_threads());

  /**
   * Apply the preconditioner.
   */
  template<typename number2>
  void vmult (Vector<number2>       &dst,
              const Vector<number2> &src) const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

private:
  /**
   * Store the number of blocks.
   */
  const unsigned int n_blocks;

  /**
   * In this field, we precompute
   * for each block which degrees
   * of freedom belong to it. Thus,
   * if <tt>dof_masks[i][j]==true</tt>,
   * then DoF @p j belongs to block
   * @p i. Of course, no other
   * <tt>dof_masks[l][j]</tt> may be
   * @p true for <tt>l!=i</tt>. This
   * computation is done in the
   * constructor, to avoid
   * recomputing each time the
   * preconditioner is called.
   */
  std::vector<std::vector<bool> > dof_masks;

  /**
   * Compute the contents of the
   * field @p dof_masks. This
   * function is called from the
   * constructor.
   */
  void compute_dof_masks (const SparseMatrix<number> &M,
                          const std::vector<bool>    &selected,
                          const BlockingStrategy      blocking_strategy);
};

/*@}*/
DEAL_II_NAMESPACE_CLOSE

#endif
