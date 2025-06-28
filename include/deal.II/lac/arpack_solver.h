// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_arpack_solver_h
#define dealii_arpack_solver_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/types.h>

#include <deal.II/lac/solver_control.h>

#include <cstring>


#ifdef DEAL_II_WITH_ARPACK

DEAL_II_NAMESPACE_OPEN


extern "C" void
dnaupd_(int          *ido,
        char         *bmat,
        unsigned int *n,
        char         *which,
        unsigned int *nev,
        const double *tol,
        double       *resid,
        int          *ncv,
        double       *v,
        int          *ldv,
        int          *iparam,
        int          *ipntr,
        double       *workd,
        double       *workl,
        int          *lworkl,
        int          *info);

extern "C" void
dsaupd_(int          *ido,
        char         *bmat,
        unsigned int *n,
        char         *which,
        unsigned int *nev,
        double       *tol,
        double       *resid,
        int          *ncv,
        double       *v,
        int          *ldv,
        int          *iparam,
        int          *ipntr,
        double       *workd,
        double       *workl,
        int          *lworkl,
        int          *info);

extern "C" void
dneupd_(int          *rvec,
        char         *howmany,
        int          *select,
        double       *d,
        double       *di,
        double       *z,
        int          *ldz,
        double       *sigmar,
        double       *sigmai,
        double       *workev,
        char         *bmat,
        unsigned int *n,
        char         *which,
        unsigned int *nev,
        double       *tol,
        double       *resid,
        int          *ncv,
        double       *v,
        int          *ldv,
        int          *iparam,
        int          *ipntr,
        double       *workd,
        double       *workl,
        int          *lworkl,
        int          *info);

extern "C" void
dseupd_(int          *rvec,
        char         *howmany,
        int          *select,
        double       *d,
        double       *z,
        int          *ldz,
        double       *sigmar,
        char         *bmat,
        unsigned int *n,
        char         *which,
        unsigned int *nev,
        double       *tol,
        double       *resid,
        int          *ncv,
        double       *v,
        int          *ldv,
        int          *iparam,
        int          *ipntr,
        double       *workd,
        double       *workl,
        int          *lworkl,
        int          *info);

/**
 * Interface for using ARPACK. ARPACK is a collection of Fortran77 subroutines
 * designed to solve large scale eigenvalue problems.  Here we interface to
 * the routines <code>dnaupd</code> and <code>dneupd</code> of ARPACK.
 * If the operator is specified to be symmetric we use the symmetric interface
 * <code>dsaupd</code> and <code>dseupd</code> of ARPACK instead.  The
 * package is designed to compute a few eigenvalues and corresponding
 * eigenvectors of a general n by n matrix A. It is most appropriate for large
 * sparse matrices A.
 *
 * In this class we make use of the method applied to the generalized
 * eigenspectrum problem $(A-\lambda B)x=0$, for $x\neq0$; where $A$ is a
 * system matrix, $B$ is a @ref GlossMassMatrix "mass matrix", and $\lambda, x$ are a set of
 * eigenvalues and eigenvectors respectively.
 *
 * The ArpackSolver can be used in application codes with serial objects in
 * the following way:
 * @code
 * SolverControl solver_control(1000, 1e-9);
 * ArpackSolver solver(solver_control);
 * solver.solve(A, B, OP, lambda, x, size_of_spectrum);
 * @endcode
 * for the generalized eigenvalue problem $Ax=B\lambda x$, where the variable
 * <code>size_of_spectrum</code> tells ARPACK the number of
 * eigenvector/eigenvalue pairs to solve for. Here, <code>lambda</code> is a
 * vector that will contain the eigenvalues computed, <code>x</code> a vector
 * that will contain the eigenvectors computed, and <code>OP</code> is an
 * inverse operation for the matrix <code>A</code>. Shift and invert
 * transformation around zero is applied.
 *
 * Through the AdditionalData the user can specify some of the parameters to
 * be set.
 *
 * For further information on how the ARPACK routines <code>dsaupd</code>,
 * <code>dseupd</code>, <code>dnaupd</code> and <code>dneupd</code> work
 * and also how to set the parameters appropriately
 * please take a look into the ARPACK manual.
 *
 * @note Whenever you eliminate degrees of freedom using AffineConstraints,
 * you generate spurious eigenvalues and eigenvectors. If you make sure
 * that the diagonals of eliminated matrix rows are all equal to one, you
 * get a single additional eigenvalue. But beware that some functions in
 * deal.II set these diagonals to rather arbitrary (from the point of view
 * of eigenvalue problems) values. See also
 * @ref step_36 "step-36"
 * for an example.
 */
class ArpackSolver : public EnableObserverPointer
{
public:
  /**
   * Declare the type for container size.
   */
  using size_type = types::global_dof_index;


  /**
   * An enum that lists the possible choices for which eigenvalues to compute
   * in the solve() function.
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
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the number of Arnoldi vectors (Lanczos
     * vectors if the problem is symmetric) to 15. Set the solver to find the
     * eigenvalues of largest magnitude for a non-symmetric problem).
     */
    explicit AdditionalData(
      const unsigned int     number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest    = largest_magnitude,
      const bool             symmetric                 = false);

    /**
     * Number of Arnoldi/Lanczos vectors. This number should be less than the
     * size of the problem but greater than 2 times the number of eigenvalues
     * (or n_eigenvalues if it is set) plus one.
     */
    const unsigned int number_of_arnoldi_vectors;

    /**
     * Specify the eigenvalues of interest.
     */
    const WhichEigenvalues eigenvalue_of_interest;

    /**
     * Specify if the problem is symmetric or not.
     */
    const bool symmetric;
  };

  /**
   * Access to the object that controls convergence.
   */
  SolverControl &
  control() const;

  /**
   * Constructor.
   */
  ArpackSolver(SolverControl        &control,
               const AdditionalData &data = AdditionalData());

  /**
   * Set initial vector for building Krylov space.
   */
  template <typename VectorType>
  void
  set_initial_vector(const VectorType &vec);

  /**
   * Set shift @p sigma for shift-and-invert spectral transformation.
   *
   * If this function is not called, the shift is assumed to be zero.
   */
  void
  set_shift(const std::complex<double> sigma);

  /**
   * Solve the generalized eigensprectrum problem $A x=\lambda B x$ by calling
   * the <code>dsaupd</code> and <code>dseupd</code> or
   * <code>dnaupd</code> and <code>dneupd</code> functions of ARPACK.
   *
   * The function returns a vector of eigenvalues of length <i>n</i> and a
   * vector of eigenvectors of length <i>n</i> in the symmetric case
   * and of length <i>n+1</i> in the non-symmetric case. In the symmetric case
   * all eigenvectors are real. In the non-symmetric case complex eigenvalues
   * always occur as complex conjugate pairs. Therefore the eigenvector for an
   * eigenvalue with nonzero complex part is stored by putting the real and
   * the imaginary parts in consecutive real-valued vectors. The eigenvector
   * of the complex conjugate eigenvalue does not need to be stored, since it
   * is just the complex conjugate of the stored eigenvector. Thus, if the last
   * n-th eigenvalue has a nonzero imaginary part, Arpack needs in total n+1
   * real-valued vectors to store real and imaginary parts of the eigenvectors.
   *
   * @param A The operator for which we want to compute eigenvalues. Actually,
   * this parameter is entirely unused.
   *
   * @param B The inner product of the underlying space, typically the mass
   * matrix. For constrained problems, it can be a partial @ref GlossMassMatrix "mass matrix", like
   * for instance the velocity mass matrix of a Stokes problem. Only its
   * function <code>vmult()</code> is used.
   *
   * @param inverse This is the possibly shifted inverse that is actually used
   * instead of <code>A</code>. Only its function <code>vmult()</code> is
   * used.
   *
   * @param eigenvalues is a vector of complex numbers in which the
   * eigenvalues are returned.
   *
   * @param eigenvectors is a <b>real</b> vector of eigenvectors, containing
   * the real parts of all eigenvectors and the imaginary parts of the
   * eigenvectors corresponding to complex conjugate eigenvalue pairs.
   * Therefore, its length should be <i>n</i> in the symmetric case and
   * <i>n+1</i> in the non-symmetric case. In the non-symmetric case the storage
   * scheme leads for example to the following pattern. Suppose that the first
   * two eigenvalues are real and the third and fourth are a complex conjugate
   * pair. Asking for three eigenpairs results in <i>[real(v1),real(v2),
   * real(v3),imag(v3)]</i>. Note that we get the same pattern if we ask for
   * four eigenpairs in this example, since the fourth eigenvector is simply the
   * complex conjugate of the third one.
   *
   * @param n_eigenvalues The purpose of this parameter is not clear, but it
   * is safe to set it to the size of <code>eigenvalues</code> or greater.
   * Leave it at its default zero, which will be reset to the size of
   * <code>eigenvalues</code> internally.
   */
  template <typename VectorType,
            typename MatrixType1,
            typename MatrixType2,
            typename INVERSE>
  void
  solve(const MatrixType1                 &A,
        const MatrixType2                 &B,
        const INVERSE                     &inverse,
        std::vector<std::complex<double>> &eigenvalues,
        std::vector<VectorType>           &eigenvectors,
        const unsigned int                 n_eigenvalues = 0);

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

  /**
   * Store an initial vector
   */
  bool                initial_vector_provided;
  std::vector<double> resid;

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
   * Exceptions.
   */
  DeclException2(ArpackExcInvalidNumberofEigenvalues,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(ArpackExcInvalidEigenvectorSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvectors " << arg2);

  DeclException2(
    ArpackExcInvalidEigenvectorSizeNonsymmetric,
    int,
    int,
    << "To store the real and complex parts of " << arg1
    << " eigenvectors in real-valued vectors, their size (currently set to "
    << arg2 << ") should be greater than or equal to " << arg1 + 1);

  DeclException2(ArpackExcInvalidEigenvalueSize,
                 int,
                 int,
                 << "Number of wanted eigenvalues " << arg1
                 << " is larger that the size of eigenvalues " << arg2);

  DeclException2(ArpackExcInvalidNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is larger that the size of the matrix " << arg2);

  DeclException2(ArpackExcSmallNumberofArnoldiVectors,
                 int,
                 int,
                 << "Number of Arnoldi vectors " << arg1
                 << " is too small to obtain " << arg2 << " eigenvalues");

  DeclException1(ArpackExcArpackIdo,
                 int,
                 << "This ido " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(ArpackExcArpackMode,
                 int,
                 << "This mode " << arg1
                 << " is not supported. Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodsaupd,
                 int,
                 << "Error with dsaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodnaupd,
                 int,
                 << "Error with dnaupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodseupd,
                 int,
                 << "Error with dseupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfodneupd,
                 int,
                 << "Error with dneupd, info " << arg1
                 << ". Check documentation of ARPACK");

  DeclException1(ArpackExcArpackInfoMaxIt,
                 int,
                 << "Maximum number " << arg1 << " of iterations reached.");

  DeclExceptionMsg(ArpackExcArpackNoShifts,
                   "No shifts could be applied during implicit"
                   " Arnoldi update, try increasing the number of"
                   " Arnoldi vectors.");
};


inline ArpackSolver::AdditionalData::AdditionalData(
  const unsigned int     number_of_arnoldi_vectors,
  const WhichEigenvalues eigenvalue_of_interest,
  const bool             symmetric)
  : number_of_arnoldi_vectors(number_of_arnoldi_vectors)
  , eigenvalue_of_interest(eigenvalue_of_interest)
  , symmetric(symmetric)
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
  // Check for possible options for asymmetric problems
  else
    {
      Assert(
        eigenvalue_of_interest != algebraically_largest,
        ExcMessage(
          "'largest algebraic part' can only be used for symmetric problems!"));
      Assert(
        eigenvalue_of_interest != algebraically_smallest,
        ExcMessage(
          "'smallest algebraic part' can only be used for symmetric problems!"));
      Assert(eigenvalue_of_interest != both_ends,
             ExcMessage(
               "'both ends' can only be used for symmetric problems!"));
    }
}


inline ArpackSolver::ArpackSolver(SolverControl        &control,
                                  const AdditionalData &data)
  : solver_control(control)
  , additional_data(data)
  , initial_vector_provided(false)
  , sigmar(0.0)
  , sigmai(0.0)
{}



inline void
ArpackSolver::set_shift(const std::complex<double> sigma)
{
  sigmar = sigma.real();
  sigmai = sigma.imag();
}



template <typename VectorType>
inline void
ArpackSolver::set_initial_vector(const VectorType &vec)
{
  initial_vector_provided = true;
  resid.resize(vec.size());
  for (size_type i = 0; i < vec.size(); ++i)
    resid[i] = vec[i];
}


template <typename VectorType,
          typename MatrixType1,
          typename MatrixType2,
          typename INVERSE>
inline void
ArpackSolver::solve(const MatrixType1 & /*system_matrix*/,
                    const MatrixType2                 &mass_matrix,
                    const INVERSE                     &inverse,
                    std::vector<std::complex<double>> &eigenvalues,
                    std::vector<VectorType>           &eigenvectors,
                    const unsigned int                 n_eigenvalues)
{
  // Problem size
  unsigned int n = eigenvectors[0].size();

  // Number of eigenvalues
  const unsigned int nev_const =
    (n_eigenvalues == 0) ? eigenvalues.size() : n_eigenvalues;
  // nev for arpack, which might change by plus one during dneupd
  unsigned int nev = nev_const;

  // check input sizes
  if (additional_data.symmetric)
    {
      Assert(nev <= eigenvectors.size(),
             ArpackExcInvalidEigenvectorSize(nev, eigenvectors.size()));
    }
  else
    Assert(nev + 1 <= eigenvectors.size(),
           ArpackExcInvalidEigenvectorSizeNonsymmetric(nev,
                                                       eigenvectors.size()));

  Assert(nev <= eigenvalues.size(),
         ArpackExcInvalidEigenvalueSize(nev, eigenvalues.size()));

  // check large enough problem size
  Assert(nev < n, ArpackExcInvalidNumberofEigenvalues(nev, n));

  Assert(additional_data.number_of_arnoldi_vectors < n,
         ArpackExcInvalidNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, n));

  // check whether we have enough Arnoldi vectors
  Assert(additional_data.number_of_arnoldi_vectors > 2 * nev + 1,
         ArpackExcSmallNumberofArnoldiVectors(
           additional_data.number_of_arnoldi_vectors, nev));

  // ARPACK mode for dsaupd/dnaupd, here only mode 3, i.e. shift-invert mode
  int mode = 3;

  // reverse communication parameter
  int ido = 0;

  // 'G' generalized eigenvalue problem 'I' standard eigenvalue problem
  char bmat[2] = "G";

  // Specify the eigenvalues of interest, possible parameters "LA" algebraically
  // largest "SA" algebraically smallest "LM" largest magnitude "SM" smallest
  // magnitude "LR" largest real part "SR" smallest real part "LI" largest
  // imaginary part "SI" smallest imaginary part "BE" both ends of spectrum
  // simultaneous.
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

  // if the starting vector is used it has to be in resid
  if (!initial_vector_provided || resid.size() != n)
    resid.resize(n, 1.);

  // number of Arnoldi basis vectors specified
  // in additional_data
  int ncv = additional_data.number_of_arnoldi_vectors;

  int                 ldv = n;
  std::vector<double> v(ldv * ncv, 0.0);

  // information to the routines
  std::vector<int> iparam(11, 0);

  iparam[0] = 1; // shift strategy

  // maximum number of iterations
  iparam[2] = control().max_steps();

  // Set the mode of dsaupd. 1 is exact shifting, 2 is user-supplied shifts,
  // 3 is shift-invert mode, 4 is buckling mode, 5 is Cayley mode.

  iparam[6] = mode;
  std::vector<int> ipntr(14, 0);

  // work arrays for ARPACK
  std::vector<double> workd(3 * n, 0.);
  int                 lworkl =
    additional_data.symmetric ? ncv * ncv + 8 * ncv : 3 * ncv * ncv + 6 * ncv;
  std::vector<double> workl(lworkl, 0.);

  // information out of the iteration
  int info = 1;

  while (ido != 99)
    {
      // call of ARPACK dsaupd/dnaupd routine
      if (additional_data.symmetric)
        dsaupd_(&ido,
                bmat,
                &n,
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
        dnaupd_(&ido,
                bmat,
                &n,
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

      if (ido == 99)
        break;

      switch (mode)
        {
          case 3:
            {
              switch (ido)
                {
                  case -1:
                    {
                      VectorType src, dst, tmp;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);
                      tmp.reinit(src);


                      for (size_type i = 0; i < src.size(); ++i)
                        src(i) = workd[ipntr[0] - 1 + i];

                      // multiplication with mass matrix M
                      mass_matrix.vmult(tmp, src);
                      // solving linear system
                      inverse.vmult(dst, tmp);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  case 1:
                    {
                      VectorType src, dst, tmp, tmp2;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);
                      tmp.reinit(src);
                      tmp2.reinit(src);

                      for (size_type i = 0; i < src.size(); ++i)
                        {
                          src(i) = workd[ipntr[2] - 1 + i];
                          tmp(i) = workd[ipntr[0] - 1 + i];
                        }
                      // solving linear system
                      inverse.vmult(dst, src);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  case 2:
                    {
                      VectorType src, dst;
                      src.reinit(eigenvectors[0]);
                      dst.reinit(src);

                      for (size_type i = 0; i < src.size(); ++i)
                        src(i) = workd[ipntr[0] - 1 + i];

                      // Multiplication with mass matrix M
                      mass_matrix.vmult(dst, src);

                      for (size_type i = 0; i < dst.size(); ++i)
                        workd[ipntr[1] - 1 + i] = dst(i);
                    }
                    break;

                  default:
                    Assert(false, ArpackExcArpackIdo(ido));
                    break;
                }
            }
            break;
          default:
            Assert(false, ArpackExcArpackMode(mode));
            break;
        }
    }

  // Set number of used iterations in SolverControl
  control().check(iparam[2], 0.);

  if (info < 0)
    {
      if (additional_data.symmetric)
        {
          Assert(false, ArpackExcArpackInfodsaupd(info));
        }
      else
        Assert(false, ArpackExcArpackInfodnaupd(info));
    }
  else
    {
      // 1 - compute eigenvectors, 0 - only eigenvalues
      int rvec = 1;

      // which eigenvectors
      char howmany = 'A';

      std::vector<int> select(ncv, 1);

      int ldz = n;

      std::vector<double> eigenvalues_real(nev + 1, 0.);
      std::vector<double> eigenvalues_im(nev + 1, 0.);

      // call of ARPACK dseupd/dneupd routine
      if (additional_data.symmetric)
        {
          std::vector<double> z(ldz * nev, 0.);
          dseupd_(&rvec,
                  &howmany,
                  select.data(),
                  eigenvalues_real.data(),
                  z.data(),
                  &ldz,
                  &sigmar,
                  bmat,
                  &n,
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
        }
      else
        {
          std::vector<double> workev(3 * ncv, 0.);
          dneupd_(&rvec,
                  &howmany,
                  select.data(),
                  eigenvalues_real.data(),
                  eigenvalues_im.data(),
                  v.data(),
                  &ldz,
                  &sigmar,
                  &sigmai,
                  workev.data(),
                  bmat,
                  &n,
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
        }

      if (info == 1)
        {
          Assert(false, ArpackExcArpackInfoMaxIt(control().max_steps()));
        }
      else if (info == 3)
        {
          Assert(false, ArpackExcArpackNoShifts());
        }
      else if (info != 0)
        {
          if (additional_data.symmetric)
            {
              Assert(false, ArpackExcArpackInfodseupd(info));
            }
          else
            Assert(false, ArpackExcArpackInfodneupd(info));
        }

      for (unsigned int i = 0; i < nev; ++i)
        for (unsigned int j = 0; j < n; ++j)
          eigenvectors[i](j) = v[i * n + j];

      for (unsigned int i = 0; i < nev_const; ++i)
        eigenvalues[i] =
          std::complex<double>(eigenvalues_real[i], eigenvalues_im[i]);
    }
}


inline SolverControl &
ArpackSolver::control() const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#else

// Make sure the scripts that create the C++20 module input files have
// something to latch on if the preprocessor #ifdef above would
// otherwise lead to an empty content of the file.
DEAL_II_NAMESPACE_OPEN
DEAL_II_NAMESPACE_CLOSE

#endif
#endif
