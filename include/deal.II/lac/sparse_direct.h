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

#ifndef __deal2__sparse_direct_h
#define __deal2__sparse_direct_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/sparse_matrix_ez.h>
#include <deal.II/lac/block_sparse_matrix.h>

#ifdef DEAL_II_WITH_MUMPS
#  include <deal.II/base/utilities.h>
#  include <dmumps_c.h>
#endif

DEAL_II_NAMESPACE_OPEN

/**
 * This class provides an interface to the sparse direct solver
 * UMFPACK (see <a
 * href="http://www.cise.ufl.edu/research/sparse/umfpack">this
 * link</a>). UMFPACK is a set of routines for solving non-symmetric
 * sparse linear systems, Ax=b, using the Unsymmetric-pattern
 * MultiFrontal method and direct sparse LU factorization. Matrices
 * may have symmetric or unsymmetrix sparsity patterns, and may have
 * unsymmetric entries. The use of this class is explained in the @ref
 * step_22 "step-22" and @ref
 * step_29 "step-29" tutorial programs.
 *
 * This matrix class implements the usual interface of
 * preconditioners, that is a function initialize(const
 * SparseMatrix<double>&matrix,const AdditionalData) for initalizing
 * and the whole set of vmult() functions common to all
 * matrices. Implemented here are only vmult() and vmult_add(), which
 * perform multiplication with the inverse matrix. Furthermore, this
 * class provides an older interface, consisting of the functions
 * factorize() and solve(). Both interfaces are interchangeable.
 *
 * @note This class only exists if support for <a
 * href="http://www.cise.ufl.edu/research/sparse/umfpack">UMFPACK</a> was
 * enabled during configure and if the <a
 * href="http://www.cise.ufl.edu/research/sparse/umfpack">UMFPACK</a> library
 * was configured. The steps to do this are explained in the deal.II ReadMe
 * file. If you do nothing at the time you configure deal.II, then this class
 * will simply not work.
 *
 * @note UMFPACK has its own license, independent of that of deal.II. If you
 * want to use the UMFPACK you have to accept that license. It is linked to
 * from the deal.II ReadMe file. UMFPACK is included courtesy of its author,
 * <a href="http://www.cise.ufl.edu/~davis/">Timothy A. Davis</a>.
 *
 *
 * <h4>Instantiations</h4>
 *
 * There are instantiations of this class for SparseMatrix<double>,
 * SparseMatrix<float>, SparseMatrixEZ<float>, SparseMatrixEZ<double>,
 * BlockSparseMatrix<double>, and BlockSparseMatrix<float>.
 *
 * @ingroup Solvers Preconditioners
 *
 * @author Wolfgang Bangerth, 2004
 */
class SparseDirectUMFPACK : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Dummy class needed for the usual initalization interface of
   * preconditioners.
   */
  class AdditionalData
  {};


  /**
   * Constructor. See the documentation of this class for the meaning of
   * the parameters to this function.
   */
  SparseDirectUMFPACK ();

  /**
   * Destructor.
   */
  ~SparseDirectUMFPACK ();

  /**
   * @name Setting up a sparse factorization
   */
  /**
   * @{
   */

  /**
   * This function does nothing. It is only here to provide a interface
   * consistent with other sparse direct solvers.
   */
  void initialize (const SparsityPattern &sparsity_pattern);

  /**
   * Factorize the matrix. This function may be called multiple times for
   * different matrices, after the object of this class has been
   * initialized for a certain sparsity pattern. You may therefore save
   * some computing time if you want to invert several matrices with the
   * same sparsity pattern. However, note that the bulk of the computing
   * time is actually spent in the factorization, so this functionality may
   * not always be of large benefit.
   *
   * In contrast to the other direct solver classes, the initialisation
   * method does nothing. Therefore initialise is not automatically called
   * by this method, when the initialization step has not been performed
   * yet.
   *
   * This function copies the contents of the matrix into its own storage;
   * the matrix can therefore be deleted after this operation, even if
   * subsequent solves are required.
   */
  template <class Matrix>
  void factorize (const Matrix &matrix);

  /**
   * Initialize memory and call SparseDirectUMFPACK::factorize.
   */
  template <class Matrix>
  void initialize(const Matrix &matrix,
                  const AdditionalData additional_data = AdditionalData());

  /**
   * @}
   */

  /**
   * @name Functions that represent the inverse of a matrix
   */
  /**
   * @{
   */

  /**
   * Preconditioner interface function. Usually, given the source vector,
   * this method returns an approximate solution of <i>Ax = b</i>. As this
   * class provides a wrapper to a direct solver, here it is actually the
   * exact solution (exact within the range of numerical accuracy of
   * course).
   *
   * In other words, this function actually multiplies with the exact
   * inverse of the matrix, $A^{-1}$.
   */
  void vmult (Vector<double> &dst,
              const Vector<double> &src) const;

  /**
   * Same as before, but for block vectors.
   */
  void vmult (BlockVector<double> &dst,
              const BlockVector<double> &src) const;

  /**
   * Same as before, but uses the transpose of the matrix, i.e. this
   * function multiplies with $A^{-T}$.
   */
  void Tvmult (Vector<double> &dst,
               const Vector<double> &src) const;

  /**
   * Same as before, but for block vectors
   */
  void Tvmult (BlockVector<double> &dst,
               const BlockVector<double> &src) const;

  /**
   * Same as vmult(), but adding to the previous solution. Not implemented
   * yet but necessary for compiling certain other classes.
   */
  void vmult_add (Vector<double> &dst,
                  const Vector<double> &src) const;

  /**
   * Same as before, but uses the transpose of the matrix, i.e. this
   * function multiplies with $A^{-T}$.
   */
  void Tvmult_add (Vector<double> &dst,
                   const Vector<double> &src) const;

  /**
   * @}
   */

  /**
   * @name Functions that solve linear systems
   */
  /**
   * @{
   */

  /**
   * Solve for a certain right hand side vector. This function may be
   * called multiple times for different right hand side vectors after the
   * matrix has been factorized. This yields a big saving in computing
   * time, since the actual solution is fast, compared to the factorization
   * of the matrix.
   *
   * The solution will be returned in place of the right hand side vector.
   *
   * If the factorization has not happened before, strange things will
   * happen. Note that we can't actually call the factorize() function from
   * here if it has not yet been called, since we have no access to the
   * actual matrix.
   *
   * If @p transpose is set to true this function solves for the transpose
   * of the matrix, i.e. $x=A^{-T}b$.
   */
  void solve (Vector<double> &rhs_and_solution, bool transpose = false) const;

  /**
   * Same as before, but for block vectors.
   */
  void solve (BlockVector<double> &rhs_and_solution, bool transpose = false) const;

  /**
   * Call the two functions factorize() and solve() in that order, i.e. perform
   * the whole solution process for the given right hand side vector.
   *
   * The solution will be returned in place of the right hand side vector.
   */
  template <class Matrix>
  void solve (const Matrix   &matrix,
              Vector<double> &rhs_and_solution,
              bool            transpose = false);

  /**
   * Same as before, but for block vectors.
   */
  template <class Matrix>
  void solve (const Matrix        &matrix,
              BlockVector<double> &rhs_and_solution,
              bool                 transpose = false);

  /**
   * @}
   */

  /**
   * One of the UMFPack routines threw an error. The error code is included
   * in the output and can be looked up in the UMFPack user manual. The
   * name of the routine is included for reference.
   */
  DeclException2 (ExcUMFPACKError, char *, int,
                  << "UMFPACK routine " << arg1
                  << " returned error status " << arg2
                  << ". See the file <bundled/umfpack/UMFPACK/Include/umfpack.h>"
                  << " for a description of 'status codes'.");

private:
  /**
   * The UMFPACK routines allocate objects in which they store information
   * about symbolic and numeric values of the decomposition. The actual
   * data type of these objects is opaque, and only passed around as void
   * pointers.
   */
  void *symbolic_decomposition;
  void *numeric_decomposition;

  /**
   * Free all memory that hasn't been freed yet.
   */
  void clear ();

  /**
   * Make sure that the arrays Ai and Ap are sorted in each row. UMFPACK
   * wants it this way. We need to have three versions of this function,
   * one for the usual SparseMatrix, one for the SparseMatrixEZ, and one
   * for the BlockSparseMatrix classes
   */
  template <typename number>
  void sort_arrays (const SparseMatrixEZ<number> &);

  template <typename number>
  void sort_arrays (const SparseMatrix<number> &);

  template <typename number>
  void sort_arrays (const BlockSparseMatrix<number> &);

  /**
   * The arrays in which we store the data for the solver.
   */
  std::vector<long int> Ap;
  std::vector<long int> Ai;
  std::vector<double> Ax;

  /**
   * Control and info arrays for the solver routines.
   */
  std::vector<double> control;
};


/**
 * This class provides an interface to the parallel sparse direct solver
 * <a href="http://mumps.enseeiht.fr">MUMPS</a>. MUMPS is direct method
 * based on a multifrontal approach, which performs a direct LU
 * factorization. The matrix coming in may have either symmetric or
 * nonsymmetric sparsity pattern.
 *
 * @note This class is useable if and only if a working installation of <a
 * href="http://mumps.enseeiht.fr">MUMPS</a> exists on your system and was
 * detected during configuration of <code>deal.II</code>.
 *
 * <h4>Instantiations</h4>
 *
 * There are instantiations of this class for SparseMatrix<double>,
 * SparseMatrix<float>, BlockSparseMatrix<double>, and
 * BlockSparseMatrix<float>.
 *
 * @author Markus Buerg, 2010
 */
class SparseDirectMUMPS
{
private:

#ifdef DEAL_II_WITH_MUMPS
  DMUMPS_STRUC_C id;
#endif // DEAL_II_WITH_MUMPS

  double   *a;
  std::vector<double> rhs;
  int      *irn;
  int      *jcn;
  types::global_dof_index n;
  types::global_dof_index nz;

  /**
   * This function initializes a MUMPS instance and hands over the system's
   * matrix <tt>matrix</tt>.
   */
  template<class Matrix>
  void initialize_matrix (const Matrix &matrix);

  /**
   * Copy the computed solution into the solution vector.
   */
  void copy_solution (Vector<double> &vector);

  /**
   *
   */
  void copy_rhs_to_mumps(const Vector<double> &rhs);

  /**
   * Flags storing whether the function <tt>initialize ()</tt> has already
   * been called.
   */
  bool initialize_called;

public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Constructor
   */
  SparseDirectMUMPS ();

  /**
   * Destructor
   */
  ~SparseDirectMUMPS ();

  /**
   * Exception
   */
  DeclException0 (ExcInitializeAlreadyCalled);

  /**
   * This function initializes a MUMPS instance and hands over the system's
   * matrix <tt>matrix</tt> and right-hand side <tt>vector</tt> to the
   * solver.
   */
  template <class Matrix>
  void initialize (const Matrix &matrix,
                   const Vector<double>       &vector);

  /**
   * This function initializes a MUMPS instance and computes the
   * factorization of the system's matrix <tt>matrix</tt>.
   */
  template <class Matrix>
  void initialize (const Matrix &matrix);

  /**
   * A function in which the linear system is solved and the solution
   * vector is copied into the given <tt>vector</tt>. The right-hand side
   * need to be supplied in initialize(matrix, vector);
   */
  void solve (Vector<double> &vector);

  /**
   * A function in which the inverse of the matrix is applied to the input
   * vector <tt>src</tt> and the solution is written into the output vector
   * <tt>dst</tt>.
   */
  void vmult (Vector<double> &dst, const Vector<double> &src);
};

DEAL_II_NAMESPACE_CLOSE

#endif // __deal2__sparse_direct_h
