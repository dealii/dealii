// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2020 by the deal.II authors
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

#ifndef dealii_parpack_solver_h
#define dealii_parpack_solver_h

#include <deal.II/base/config.h>

#include <deal.II/base/index_set.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/base/smartpointer.h>

#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_operation.h>

#include <cstring>


#ifdef DEAL_II_ARPACK_WITH_PARPACK

DEAL_II_NAMESPACE_OPEN

extern "C"
{
  // http://www.mathkeisan.com/usersguide/man/pdnaupd.html
  void
  pdnaupd_(MPI_Fint *comm,
           int *     ido,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdsaupd.html
  void
  pdsaupd_(MPI_Fint *comm,
           int *     ido,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdneupd.html
  void
  pdneupd_(MPI_Fint *comm,
           int *     rvec,
           char *    howmany,
           int *     select,
           double *  d,
           double *  di,
           double *  z,
           int *     ldz,
           double *  sigmar,
           double *  sigmai,
           double *  workev,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // http://www.mathkeisan.com/usersguide/man/pdseupd.html
  void
  pdseupd_(MPI_Fint *comm,
           int *     rvec,
           char *    howmany,
           int *     select,
           double *  d,
           double *  z,
           int *     ldz,
           double *  sigmar,
           char *    bmat,
           int *     n,
           char *    which,
           int *     nev,
           double *  tol,
           double *  resid,
           int *     ncv,
           double *  v,
           int *     nloc,
           int *     iparam,
           int *     ipntr,
           double *  workd,
           double *  workl,
           int *     lworkl,
           int *     info);

  // other resources:
  //    http://acts.nersc.gov/superlu/example5/pnslac.c.html
  //    https://github.com/phpisciuneri/tijo/blob/master/dvr_parpack.cpp
}

/**
 * Interface for using PARPACK. PARPACK is a collection of Fortran77
 * subroutines designed to solve large scale eigenvalue problems. Here we
 * interface to the routines <code>pdneupd</code>, <code>pdseupd</code>,
 * <code>pdnaupd</code>, <code>pdsaupd</code> of PARPACK.  The package is
 * designed to compute a few eigenvalues and corresponding eigenvectors of a
 * general n by n matrix A. It is most appropriate for large sparse matrices
 * A.
 *
 * In this class we make use of the method applied to the generalized
 * eigenspectrum problem $(A-\lambda B)x=0$, for $x\neq0$; where $A$ is a
 * system matrix, $B$ is a mass matrix, and $\lambda, x$ are a set of
 * eigenvalues and eigenvectors respectively.
 *
 * The ArpackSolver can be used in application codes in the following way:
 * @code
 *   SolverControl solver_control (1000, 1e-9);
 *   const unsigned int num_arnoldi_vectors = 2*size_of_spectrum + 2;
 *   PArpackSolver<V>::AdditionalData
 *     additional_data(num_arnoldi_vectors,
 *                     dealii::PArpackSolver<V>::largest_magnitude,
 *                     true);
 *
 *    PArpackSolver<V> eigensolver (solver_control,
 *                                  mpi_communicator,
 *                                  additional_data);
 *    eigensolver.set_shift(sigma);
 *    eigensolver.reinit(locally_owned_dofs);
 *    eigensolver.solve (A,
 *                       B,
 *                       OP,
 *                       lambda,
 *                       x,
 *                       size_of_spectrum);
 * @endcode
 * for the generalized eigenvalue problem $Ax=B\lambda x$, where the variable
 * <code>size_of_spectrum</code> tells PARPACK the number of
 * eigenvector/eigenvalue pairs to solve for. Here, <code>lambda</code> is a
 * vector that will contain the eigenvalues computed, <code>x</code> a vector
 * of objects of type <code>V</code> that will contain the eigenvectors
 * computed.
 *
 * Currently, only three modes of (P)Arpack are implemented. In mode 3
 * (default), <code>OP</code> is an inverse operation for the matrix <code>A -
 * sigma * B</code>, where <code> sigma </code> is a shift value, set to zero
 * by default. Whereas in mode 2, <code>OP</code> is an inverse of
 * <code>M</code>. Finally, mode 1 corresponds to standard eigenvalue problem
 * without spectral transformation $Ax=\lambda x$. The mode can be specified via
 * AdditionalData object. Note that for shift-and-invert (mode=3), the sought
 * eigenpairs are those after the spectral transformation is applied.
 *
 * The <code>OP</code> can be specified by using a LinearOperator:
 * @code
 *   const double shift = 5.0;
 *   const auto op_A = linear_operator<vector_t>(A);
 *   const auto op_B = linear_operator<vector_t>(B);
 *   const auto op_shift = op_A - shift * op_B;
 *   SolverControl solver_control_lin (1000, 1e-10,false,false);
 *
 *   SolverCG<vector_t> cg(solver_control_lin);
 *   const auto op_shift_invert =
 *     inverse_operator(op_shift, cg, PreconditionIdentity ());
 * @endcode
 *
 * The class is intended to be used with MPI and can work on arbitrary vector
 * and matrix distributed classes.  Both symmetric and non-symmetric
 * <code>A</code> are supported.
 *
 * For further information on how the PARPACK routines <code>pdneupd</code>,
 * <code>pdseupd</code>, <code>pdnaupd</code>, <code>pdsaupd</code> work and
 * also how to set the parameters appropriately please take a look into the
 * PARPACK manual.
 *
 * @author Denis Davydov, 2015, 2017
 */
template <typename VectorType>
class PArpackSolver : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  using size_type = types::global_dof_index;

  /**
   * An enum that lists the possible choices for which eigenvalues to compute
   * in the solve() function. Note, that this corresponds to the problem after
   * shift-and-invert (the only currently supported spectral transformation)
   * is applied.
   *
   * A particular choice is limited based on symmetric or non-symmetric matrix
   * <code>A</code> considered.
   */
  enum WhichEigenvalues
  {
    /**
     * The algebraically largest eigenvalues.
     */
    algebraically_largest,
    /**
     * The algebraically smallest eigenvalues.
     */
    algebraically_smallest,
    /**
     * The eigenvalue with the largest magnitudes.
     */
    largest_magnitude,
    /**
     * The eigenvalue with the smallest magnitudes.
     */
    smallest_magnitude,
    /**
     * The eigenvalues with the largest real parts.
     */
    largest_real_part,
    /**
     * The eigenvalues with the smallest real parts.
     */
    smallest_real_part,
    /**
     * The eigenvalues with the largest imaginary parts.
     */
    largest_imaginary_part,
    /**
     * The eigenvalues with the smallest imaginary parts.
     */
    smallest_imaginary_part,
    /**
     * Compute half of the eigenvalues from the high end of the spectrum and
     * the other half from the low end. If the number of requested
     * eigenvectors is odd, then the extra eigenvector comes from the high end
     * of the spectrum.
     */
    both_ends
  };

  /**
   * Standardized data struct to pipe additional data to the solver, should it
   * be needed.
   */
  struct AdditionalData
  {
    const unsigned int     number_of_arnoldi_vectors;
    const WhichEigenvalues eigenvalue_of_interest;
    const bool             symmetric;
    const int              mode;
    AdditionalData(
      const unsigned int     number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest    = largest_magnitude,
      const bool             symmetric                 = false,
      const int              mode                      = 3);
  };

  /**
   * Access to the object that controls convergence.
   */
  SolverControl &
  control() const;

  /**
   * Constructor.
   */
  PArpackSolver(SolverControl &       control,
                const MPI_Comm &      mpi_communicator,
                const AdditionalData &data = AdditionalData());

  /**
   * Initialize internal variables.
   */
  void
  reinit(const IndexSet &locally_owned_dofs);

  /**
   * Initialize internal variables when working with BlockVectors.
   * @p locally_owned_dofs is used to set the dimension of the problem,
   * whereas @p partitioning is used for calling the reinit of the deal.II
   * blockvectors used.
   */
  void
  reinit(const IndexSet &             locally_owned_dofs,
         const std::vector<IndexSet> &partitioning);

  /**
   * Initialize internal variables from the input @p distributed_vector.
   */
  void
  reinit(const VectorType &distributed_vector);

  /**
   * Set initial vector for building Krylov space.
   */
  void
  set_initial_vector(const VectorType &vec);

  /**
   * Set shift @p sigma for shift-and-invert spectral transformation.
   *
   * If this function is not called, the shift is assumed to be zero.
   *
   * @note only relevant for <code>mode=3</code> (see the general documentation of this
   * class for a definition of what the different modes are).
   */
  void
  set_shift(const std::complex<double> sigma);

  /**
   * Solve the generalized eigensprectrum problem $A x=\lambda B x$ by calling
   * the <code>pd(n/s)eupd</code> and <code>pd(n/s)aupd</code> functions of
   * PARPACK.
   *
   * In <code>mode=3</code>, @p inverse should correspond to $[A-\sigma B]^{-1}$,
   * whereas in <code>mode=2</code> it should represent $B^{-1}$. For
   * <code>mode=1</code> both @p B and @p inverse are ignored.
   */
  template <typename MatrixType1, typename MatrixType2, typename INVERSE>
  void
  solve(const MatrixType1 &                A,
        const MatrixType2 &                B,
        const INVERSE &                    inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType> &          eigenvectors,
        const unsigned int                 n_eigenvalues);

  /**
   * Same as above but takes eigenvectors as pointers.
   */
  template <typename MatrixType1, typename MatrixType2, typename INVERSE>
  void
  solve(const MatrixType1 &                A,
        const MatrixType2 &                B,
        const INVERSE &                    inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType *> &        eigenvectors,
        const unsigned int                 n_eigenvalues);

  /**
   * Return the memory consumption of this class in bytes.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * Reference to the object that controls convergence of the iterative
   * solver.
   */
  SolverControl &solver_control;

  /**
   * Store a copy of the flags for this particular solver.
   */
  const AdditionalData additional_data;

  // keep MPI communicator non-const as Arpack functions are not const either:

  /**
   * C++ MPI communicator.
   */
  MPI_Comm mpi_communicator;

  /**
   * Fortran MPI communicator.
   */
  MPI_Fint mpi_communicator_fortran;

  // C++98 guarantees that the elements of a vector are stored contiguously

  /**
   * Length of the work array workl.
   */
  int lworkl;

  /**
   * Double precision  work array of length lworkl
   */
  std::vector<double> workl;

  /**
   * Double precision  work array of length 3*N
   */
  std::vector<double> workd;

  /**
   * Number of local degrees of freedom.
   */
  int nloc;

  /**
   * Number of Arnoldi basis vectors specified in additional_data
   */
  int ncv;


  /**
   * The leading dimension of the array v
   */
  int ldv;

  /**
   * Double precision vector of size ldv by NCV.  Will contains the final set
   * of Arnoldi basis vectors.
   */
  std::vector<double> v;

  /**
   * An auxiliary flag which is set to true when initial vector is provided.
   */
  bool initial_vector_provided;

  /**
   * The initial residual vector, possibly from a previous run.  On output, it
   * contains the final residual vector.
   */
  std::vector<double> resid;

  /**
   * The leading dimension of the array Z equal to nloc.
   */
  int ldz;

  /**
   * A vector of minimum size of nloc by NEV+1.  Z contains the B-orthonormal
   * Ritz vectors of the eigensystem A*z = lambda*B*z corresponding to the
   * Ritz value approximations.
   */
  std::vector<double> z;

  /**
   * The size of the workev array.
   */
  int lworkev;

  /**
   * Double precision  work array of dimension 3* NCV.
   */
  std::vector<double> workev;

  /**
   * A vector of dimension NCV.
   */
  std::vector<int> select;

  /**
   * Temporary vectors used between Arpack and deal.II
   */
  VectorType src, dst, tmp;

  /**
   * Indices of local degrees of freedom.
   */
  std::vector<types::global_dof_index> local_indices;

  /**
   * Real part of the shift
   */
  double sigmar;

  /**
   * Imaginary part of the shift
   */
  double sigmai;

private:
  /**
   * Initialize internal variables which depend on
   * @p locally_owned_dofs.
   *
   * This function is called inside the reinit() functions
   */
  void
  internal_reinit(const IndexSet &locally_owned_dofs);

  /**
   * PArpackExcInfoPdnaupds.
   */
  DeclException2(PArpackExcConvergedEigenvectors,
                 int,
                 int,
                 << arg1 << " eigenpairs were requested, but only " << arg2
                 << " converged");

  DeclException2(PArpackExcInvalidNumberofEigenvalues,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(PArpackExcInvalidEigenvectorSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvectors " << arg2);

  DeclException2(
    PArpackExcInvalidEigenvectorSizeNonsymmetric,
    int,
    int,
    << "To store the real and complex parts of " << arg1
    << " eigenvectors in real-valued vectors, their size (currently set to "
    << arg2 << ") should be greater than or equal to " << arg1 + 1);

  DeclException2(PArpackExcInvalidEigenvalueSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvalues " << arg2);

  DeclException2(PArpackExcInvalidNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(PArpackExcSmallNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is too small to obtain " << arg2 << " eigenvalues");

  DeclException1(PArpackExcIdo,
                 int,
                 << "This ido " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(PArpackExcMode,
                 int,
                 << "This mode " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(PArpackExcInfoPdnaupd,
                 int,
                 << "Error with Pdnaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(PArpackExcInfoPdneupd,
                 int,
                 << "Error with Pdneupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(PArpackExcInfoMaxIt,
                 int,
                 << "Maximum number " << arg1 << " of iterations reached.");

  DeclException1(PArpackExcNoShifts,
                 int,
                 << "No shifts could be applied during implicit"
                 << " Arnoldi update, try increasing the number of"
                 << " Arnoldi vectors.");
};



template <typename VectorType>
std::size_t
PArpackSolver<VectorType>::memory_consumption() const
{
  return MemoryConsumption::memory_consumption(double()) *
           (workl.size() + workd.size() + v.size() + resid.size() + z.size() +
            workev.size()) +
         src.memory_consumption() + dst.memory_consumption() +
         tmp.memory_consumption() +
         MemoryConsumption::memory_consumption(types::global_dof_index()) *
           local_indices.size();
}



template <typename VectorType>
PArpackSolver<VectorType>::AdditionalData::AdditionalData(
  const unsigned int     number_of_arnoldi_vectors,
  const WhichEigenvalues eigenvalue_of_interest,
  const bool             symmetric,
  const int              mode)
  : number_of_arnoldi_vectors(number_of_arnoldi_vectors)
  , eigenvalue_of_interest(eigenvalue_of_interest)
  , symmetric(symmetric)
  , mode(mode)
{
  // Check for possible options for symmetric problems
  if (symmetric)
    {
      Assert(
        eigenvalue_of_interest != largest_real_part,
        ExcMessage(
          "'largest real part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != smallest_real_part,
        ExcMessage(
          "'smallest real part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != largest_imaginary_part,
        ExcMessage(
          "'largest imaginary part' can only be used for non-symmetric problems!"));
      Assert(
        eigenvalue_of_interest != smallest_imaginary_part,
        ExcMessage(
          "'smallest imaginary part' can only be used for non-symmetric problems!"));
    }
  Assert(mode >= 1 && mode <= 3,
         ExcMessage("Currently, only modes 1, 2 and 3 are supported."));
}



template <typename VectorType>
PArpackSolver<VectorType>::PArpackSolver(SolverControl &       control,
                                         const MPI_Comm &      mpi_communicator,
                                         const AdditionalData &data)
  : solver_control(control)
  , additional_data(data)
  , mpi_communicator(mpi_communicator)
  , mpi_communicator_fortran(MPI_Comm_c2f(mpi_communicator))
  , lworkl(0)
  , nloc(0)
  , ncv(0)
  , ldv(0)
  , initial_vector_provided(false)
  , ldz(0)
  , lworkev(0)
  , sigmar(0.0)
  , sigmai(0.0)
{}



template <typename VectorType>
void
PArpackSolver<VectorType>::set_shift(const std::complex<double> sigma)
{
  sigmar = sigma.real();
  sigmai = sigma.imag();
}



template <typename VectorType>
void
PArpackSolver<VectorType>::set_initial_vector(const VectorType &vec)
{
  initial_vector_provided = true;
  Assert(resid.size() == local_indices.size(),
         ExcDimensionMismatch(resid.size(), local_indices.size()));
  vec.extract_subvector_to(local_indices.begin(),
                           local_indices.end(),
                           resid.data());
}



template <typename VectorType>
void
PArpackSolver<VectorType>::internal_reinit(const IndexSet &locally_owned_dofs)
{
  // store local indices to write to vectors
  locally_owned_dofs.fill_index_vector(local_indices);

  // scalars
  nloc = locally_owned_dofs.n_elements();
  ncv  = additional_data.number_of_arnoldi_vectors;

  AssertDimension(local_indices.size(), nloc);

  // vectors
  ldv = nloc;
  v.resize(ldv * ncv, 0.0);

  resid.resize(nloc, 1.0);

  // work arrays for ARPACK
  workd.resize(3 * nloc, 0.0);

  lworkl =
    additional_data.symmetric ? ncv * ncv + 8 * ncv : 3 * ncv * ncv + 6 * ncv;
  workl.resize(lworkl, 0.);

  ldz = nloc;
  z.resize(ldz * ncv, 0.); // TODO we actually need only ldz*nev

  // WORKEV  Double precision  work array of dimension 3*NCV.
  lworkev = additional_data.symmetric ? 0 /*not used in symmetric case*/
                                        :
                                        3 * ncv;
  workev.resize(lworkev, 0.);

  select.resize(ncv, 0);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const IndexSet &locally_owned_dofs)
{
  internal_reinit(locally_owned_dofs);

  // deal.II vectors:
  src.reinit(locally_owned_dofs, mpi_communicator);
  dst.reinit(locally_owned_dofs, mpi_communicator);
  tmp.reinit(locally_owned_dofs, mpi_communicator);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const VectorType &distributed_vector)
{
  internal_reinit(distributed_vector.locally_owned_elements());

  // deal.II vectors:
  src.reinit(distributed_vector);
  dst.reinit(distributed_vector);
  tmp.reinit(distributed_vector);
}



template <typename VectorType>
void
PArpackSolver<VectorType>::reinit(const IndexSet &locally_owned_dofs,
                                  const std::vector<IndexSet> &partitioning)
{
  internal_reinit(locally_owned_dofs);

  // deal.II vectors:
  src.reinit(partitioning, mpi_communicator);
  dst.reinit(partitioning, mpi_communicator);
  tmp.reinit(partitioning, mpi_communicator);
}



template <typename VectorType>
template <typename MatrixType1, typename MatrixType2, typename INVERSE>
void
PArpackSolver<VectorType>::solve(const MatrixType1 &                A,
                                 const MatrixType2 &                B,
                                 const INVERSE &                    inverse,
                                 std::vector<std::complex<double>> &eigenvalues,
                                 std::vector<VectorType> &eigenvectors,
                                 const unsigned int       n_eigenvalues)
{
  std::vector<VectorType *> eigenvectors_ptr(eigenvectors.size());
  for (unsigned int i = 0; i < eigenvectors.size(); ++i)
    eigenvectors_ptr[i] = &eigenvectors[i];
  solve(A, B, inverse, eigenvalues, eigenvectors_ptr, n_eigenvalues);
}



template <typename VectorType>
template <typename MatrixType1, typename MatrixType2, typename INVERSE>
void
PArpackSolver<VectorType>::solve(const MatrixType1 &system_matrix,
                                 const MatrixType2 &mass_matrix,
                                 const INVERSE &    inverse,
                                 std::vector<std::complex<double>> &eigenvalues,
                                 std::vector<VectorType *> &eigenvectors,
                                 const unsigned int         n_eigenvalues)
{
  if (additional_data.symmetric)
    {
      Assert(n_eigenvalues <= eigenvectors.size(),
             PArpackExcInvalidEigenvectorSize(n_eigenvalues,
                                              eigenvectors.size()));
    }
  else
    Assert(n_eigenvalues + 1 <= eigenvectors.size(),
           PArpackExcInvalidEigenvectorSizeNonsymmetric(n_eigenvalues,
                                                        eigenvectors.size()));

  Assert(n_eigenvalues <= eigenvalues.size(),
         PArpackExcInvalidEigenvalueSize(n_eigenvalues, eigenvalues.size()));


  // use eigenvectors to get the problem size so that it is possible to
  // employ LinearOperator for mass_matrix.
  Assert(n_eigenvalues < eigenvectors[0]->size(),
         PArpackExcInvalidNumberofEigenvalues(n_eigenvalues,
                                              eigenvectors[0]->size()));

  Assert(additional_data.number_of_arnoldi_vectors < eigenvectors[0]->size(),
         PArpackExcInvalidNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, eigenvectors[0]->size()));

  Assert(additional_data.number_of_arnoldi_vectors > 2 * n_eigenvalues + 1,
         PArpackExcSmallNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, n_eigenvalues));

  int mode = additional_data.mode;

  // reverse communication parameter
  // must be zero on the first call to pdnaupd
  int ido = 0;

  // 'G' generalized eigenvalue problem
  // 'I' standard eigenvalue problem
  char bmat[2];
  bmat[0] = (mode == 1) ? 'I' : 'G';
  bmat[1] = '\0';

  // Specify the eigenvalues of interest, possible parameters:
  // "LA" algebraically largest
  // "SA" algebraically smallest
  // "LM" largest magnitude
  // "SM" smallest magnitude
  // "LR" largest real part
  // "SR" smallest real part
  // "LI" largest imaginary part
  // "SI" smallest imaginary part
  // "BE" both ends of spectrum simultaneous
  char which[3];
  switch (additional_data.eigenvalue_of_interest)
    {
      case algebraically_largest:
        std::strcpy(which, "LA");
        break;
      case algebraically_smallest:
        std::strcpy(which, "SA");
        break;
      case largest_magnitude:
        std::strcpy(which, "LM");
        break;
      case smallest_magnitude:
        std::strcpy(which, "SM");
        break;
      case largest_real_part:
        std::strcpy(which, "LR");
        break;
      case smallest_real_part:
        std::strcpy(which, "SR");
        break;
      case largest_imaginary_part:
        std::strcpy(which, "LI");
        break;
      case smallest_imaginary_part:
        std::strcpy(which, "SI");
        break;
      case both_ends:
        std::strcpy(which, "BE");
        break;
    }

  // tolerance for ARPACK
  double tol = control().tolerance();

  // information to the routines
  std::vector<int> iparam(11, 0);

  iparam[0] = 1;
  // shift strategy: exact shifts with respect to the current Hessenberg matrix
  // H.

  // maximum number of iterations
  iparam[2] = control().max_steps();

  // Parpack currently works only for NB = 1
  iparam[3] = 1;

  // Sets the mode of dsaupd:
  // 1 is A*x=lambda*x, OP = A, B = I
  // 2 is A*x = lambda*M*x, OP = inv[M]*A, B = M
  // 3 is shift-invert mode, OP = inv[A-sigma*M]*M, B = M
  // 4 is buckling mode,
  // 5 is Cayley mode.

  iparam[6] = mode;
  std::vector<int> ipntr(14, 0);

  // information out of the iteration
  //  If INFO .EQ. 0, a random initial residual vector is used.
  //  If INFO .NE. 0, RESID contains the initial residual vector,
  //  possibly from a previous run.
  // Typical choices in this situation might be to use the final value
  // of the starting vector from the previous eigenvalue calculation
  int info = initial_vector_provided ? 1 : 0;

  // Number of eigenvalues of OP to be computed. 0 < NEV < N.
  int nev             = n_eigenvalues;
  int n_inside_arpack = nloc;

  // IDO = 99: done
  while (ido != 99)
    {
      // call of ARPACK pdnaupd routine
      if (additional_data.symmetric)
        pdsaupd_(&mpi_communicator_fortran,
                 &ido,
                 bmat,
                 &n_inside_arpack,
                 which,
                 &nev,
                 &tol,
                 resid.data(),
                 &ncv,
                 v.data(),
                 &ldv,
                 iparam.data(),
                 ipntr.data(),
                 workd.data(),
                 workl.data(),
                 &lworkl,
                 &info);
      else
        pdnaupd_(&mpi_communicator_fortran,
                 &ido,
                 bmat,
                 &n_inside_arpack,
                 which,
                 &nev,
                 &tol,
                 resid.data(),
                 &ncv,
                 v.data(),
                 &ldv,
                 iparam.data(),
                 ipntr.data(),
                 workd.data(),
                 workl.data(),
                 &lworkl,
                 &info);

      AssertThrow(info == 0, PArpackExcInfoPdnaupd(info));

      // if we converge, we shall not modify anything in work arrays!
      if (ido == 99)
        break;

      // IPNTR(1) is the pointer into WORKD for X,
      // IPNTR(2) is the pointer into WORKD for Y.
      const int shift_x = ipntr[0] - 1;
      const int shift_y = ipntr[1] - 1;
      Assert(shift_x >= 0, dealii::ExcInternalError());
      Assert(shift_x + nloc <= static_cast<int>(workd.size()),
             dealii::ExcInternalError());
      Assert(shift_y >= 0, dealii::ExcInternalError());
      Assert(shift_y + nloc <= static_cast<int>(workd.size()),
             dealii::ExcInternalError());

      src = 0.;

      // switch based on both ido and mode
      if ((ido == -1) || (ido == 1 && mode < 3))
        // compute  Y = OP * X
        {
          src.add(nloc, local_indices.data(), workd.data() + shift_x);
          src.compress(VectorOperation::add);

          if (mode == 3)
            // OP = inv[K - sigma*M]*M
            {
              mass_matrix.vmult(tmp, src);
              inverse.vmult(dst, tmp);
            }
          else if (mode == 2)
            // OP = inv[M]*K
            {
              system_matrix.vmult(tmp, src);
              // store M*X in X
              tmp.extract_subvector_to(local_indices.begin(),
                                       local_indices.end(),
                                       workd.data() + shift_x);
              inverse.vmult(dst, tmp);
            }
          else if (mode == 1)
            {
              system_matrix.vmult(dst, src);
            }
          else
            AssertThrow(false, PArpackExcMode(mode));
        }
      else if (ido == 1 && mode >= 3)
        // compute  Y = OP * X for mode 3, 4 and 5, where
        // the vector B * X is already available in WORKD(ipntr(3)).
        {
          const int shift_b_x = ipntr[2] - 1;
          Assert(shift_b_x >= 0, dealii::ExcInternalError());
          Assert(shift_b_x + nloc <= static_cast<int>(workd.size()),
                 dealii::ExcInternalError());

          // B*X
          src.add(nloc, local_indices.data(), workd.data() + shift_b_x);
          src.compress(VectorOperation::add);

          // solving linear system
          Assert(mode == 3, ExcNotImplemented());
          inverse.vmult(dst, src);
        }
      else if (ido == 2)
        // compute  Y = B * X
        {
          src.add(nloc, local_indices.data(), workd.data() + shift_x);
          src.compress(VectorOperation::add);

          // Multiplication with mass matrix M
          if (mode == 1)
            {
              dst = src;
            }
          else
            // mode 2,3 and 5 have B=M
            {
              mass_matrix.vmult(dst, src);
            }
        }
      else
        AssertThrow(false, PArpackExcIdo(ido));
      // Note: IDO = 3 does not appear to be required for currently
      // implemented modes

      // store the result
      dst.extract_subvector_to(local_indices.begin(),
                               local_indices.end(),
                               workd.data() + shift_y);
    } // end of pd*aupd_ loop

  // 1 - compute eigenvectors,
  // 0 - only eigenvalues
  int rvec = 1;

  // which eigenvectors
  char howmany[4] = "All";

  std::vector<double> eigenvalues_real(n_eigenvalues + 1, 0.);
  std::vector<double> eigenvalues_im(n_eigenvalues + 1, 0.);

  // call of ARPACK pdneupd routine
  if (additional_data.symmetric)
    pdseupd_(&mpi_communicator_fortran,
             &rvec,
             howmany,
             select.data(),
             eigenvalues_real.data(),
             z.data(),
             &ldz,
             &sigmar,
             bmat,
             &n_inside_arpack,
             which,
             &nev,
             &tol,
             resid.data(),
             &ncv,
             v.data(),
             &ldv,
             iparam.data(),
             ipntr.data(),
             workd.data(),
             workl.data(),
             &lworkl,
             &info);
  else
    pdneupd_(&mpi_communicator_fortran,
             &rvec,
             howmany,
             select.data(),
             eigenvalues_real.data(),
             eigenvalues_im.data(),
             v.data(),
             &ldz,
             &sigmar,
             &sigmai,
             workev.data(),
             bmat,
             &n_inside_arpack,
             which,
             &nev,
             &tol,
             resid.data(),
             &ncv,
             v.data(),
             &ldv,
             iparam.data(),
             ipntr.data(),
             workd.data(),
             workl.data(),
             &lworkl,
             &info);

  if (info == 1)
    {
      AssertThrow(false, PArpackExcInfoMaxIt(control().max_steps()));
    }
  else if (info == 3)
    {
      AssertThrow(false, PArpackExcNoShifts(1));
    }
  else if (info != 0)
    {
      AssertThrow(false, PArpackExcInfoPdneupd(info));
    }

  for (int i = 0; i < nev; ++i)
    {
      (*eigenvectors[i]) = 0.0;
      AssertIndexRange(i * nloc + nloc, v.size() + 1);

      eigenvectors[i]->add(nloc, local_indices.data(), &v[i * nloc]);
      eigenvectors[i]->compress(VectorOperation::add);
    }

  for (size_type i = 0; i < n_eigenvalues; ++i)
    eigenvalues[i] =
      std::complex<double>(eigenvalues_real[i], eigenvalues_im[i]);

  // Throw an error if the solver did not converge.
  AssertThrow(iparam[4] >= static_cast<int>(n_eigenvalues),
              PArpackExcConvergedEigenvectors(n_eigenvalues, iparam[4]));

  // both PDNAUPD and PDSAUPD compute eigenpairs of inv[A - sigma*M]*M
  // with respect to a semi-inner product defined by M.

  // resid likely contains residual with respect to M-norm.
  {
    tmp = 0.0;
    tmp.add(nloc, local_indices.data(), resid.data());
    solver_control.check(iparam[2], tmp.l2_norm());
  }
}



template <typename VectorType>
SolverControl &
PArpackSolver<VectorType>::control() const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#endif
#endif
