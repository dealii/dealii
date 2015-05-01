// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2015 by the deal.II authors
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

#ifndef __deal2__arpack_solver_h
#define __deal2__arpack_solver_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector_memory.h>

#include <cstring>

#ifdef DEAL_II_WITH_ARPACK

DEAL_II_NAMESPACE_OPEN


extern "C" void dnaupd_(int *ido, const char *bmat, const unsigned int *n,
                        const char *which, const unsigned int *nev, const
                        double *tol, double *resid, const unsigned int
                        *ncv, double *v, const unsigned int *ldv, int
                        *iparam, int *ipntr, double *workd, double *workl,
                        int *lworkl, int *info);

extern "C" void dneupd_(int *rvec, const char *howmany, const int *select,
                        double *d, double *di, double *z, const unsigned
                        int *ldz, double *sigmar, double *sigmai, double
                        *workev, const char *bmat, const unsigned int *n,
                        const char *which, const unsigned int *nev, const
                        double *tol, double *resid, const unsigned int
                        *ncv, double *v, const unsigned int *ldv, int
                        *iparam, int *ipntr, double *workd, double *workl,
                        int *lworkl, int *info);


/**
 * Interface for using ARPACK. ARPACK is a collection of Fortran77
 * subroutines designed to solve large scale eigenvalue problems.  Here we
 * interface to the routines <code>dneupd</code> and <code>dnaupd</code> of
 * ARPACK. The package is designed to compute a few eigenvalues and
 * corresponding eigenvectors of a general $n\times n$ matrix A. It is most
 * appropriate for large sparse matrices A.
 *
 * In this class we make use of the method applied to the generalized
 * eigenspectrum problem $(A-\lambda B)x=0$, for $x\neq0$; where $A$ is a
 * system matrix, $B$ is a mass matrix, and $\lambda, x$ are a set of
 * eigenvalues and eigenvectors respectively.
 *
 * As an example, the ArpackSolver can be used to compute 10 eigenvalues for
 * the generalized eigenvalue problem $Ax=\lambda Bx$ in the following way:
 * @code
 * SolverControl solver_control(1000, 1e-9);
 * ArpackSolver arpack_solver(solver_control);
 *
 * std::vector<std::complex<double> > eigenvalues(10);
 * std::vector<Vector<double> > eigenvectors(20);
 * for (std::vector::size_type i = 0; i < eigenvectors.size; i++)
 *   {
 *     eigenvectors[i].reinit(...); // system size, e.g. n_dofs
 *   }
 * arpack_solver.solve (B, OP, eigenvalues, eigenvectors);
 * @endcode
 * Here, <code>eigenvalues</code> is a vector that will contain the computed
 * eigenvalues and <code>eigenvectors</code> is a vector that will contain
 * the computed eigenvectors; the first half containing the real parts, the
 * second half containing the imaginary parts of the eigenvectors.
 *
 * <code>OP</code> is an inverse operation for the matrix <code>A</code>.
 * Shift and invert transformation around zero is applied.
 *
 * For further information on how the ARPACK routines <code>dneupd</code> and
 * <code>dnaupd</code> work and also how to set the parameters appropriately
 * please take a look into the ARPACK manual.
 *
 * @note Whenever you eliminate degrees of freedom using ConstraintMatrix,
 * you generate spurious eigenvalues and eigenvectors. If you make sure
 * that the diagonals of eliminated matrix rows are all equal to one, you
 * get a single additional eigenvalue. But beware that some functions in
 * deal.II set these diagonals to rather arbitrary (from the point of view
 * of eigenvalue problems) values. See also step-36 for an example.
 *
 * @author Baerbel Janssen, Agnieszka Miedlar, 2010, Guido Kanschat, Matthias
 *Maier 2015
 */
class ArpackSolver : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;


  /**
   * An enum that lists the possible choices for which eigenvalues to compute
   * in the solve() function.
   */
  enum WhichEigenvalues
  {
    algebraically_largest,
    algebraically_smallest,
    largest_magnitude,
    smallest_magnitude,
    largest_real_part,
    smallest_real_part,
    largest_imaginary_part,
    smallest_imaginary_part,
    both_ends
  };


  /**
   * Access to the object that controls convergence.
   */
  SolverControl &control () const;


  /**
   * Constructor.
   */
  ArpackSolver(SolverControl &control);


  /**
   * Solve the generalized eigensprectrum problem $A x=\lambda B x$ by calling
   * the <code>dneupd</code> and <code>dnaupd</code> functions of
   * ARPACK.
   *
   * The function returns a vector of eigenvalues of length <i>n</i>
   * and a vector of eigenvectors, where the latter should be twice
   * the size of the eigenvalue vector. The first <i>n</i> vectors in
   * <code>eigenvectors</code> will be the real parts of the
   * eigenvectors, the second <i>n</i> the imaginary parts.
   *
   * @param B The inner product of the underlying space, typically the
   * mass matrix. Only its function <code>vmult()</code> is used.
   *
   * @param inverse This is the (possibly shifted) inverse that is actually
   * used instead of <code>A</code>. Only its function <code>vmult()</code>
   * is used.
   *
   * @param eigenvalues is a vector of complex numbers in which the
   * eigenvalues are returned. The size of the vector determines the number
   * of eigenvalues to be computed.
   *
   * @param eigenvectors is a <b>real</b> vector of eigenvectors,
   * the first half containing the real parts and the second half
   * containing the corresponding imaginary parts of the eigenvectors.
   * Therefore, its length should be twice the number of eigenvalues. The
   * vectors have to be initialized to match the matrices.
   *
   * @param eigenvalue_of_interest encodes the eigenvalue of interest given
   * by the enum WhichEigenvalues
   *
   * @param number_of_arnoldi_vectors sets the number of internal Storage
   * for the Arnoldi procedure. If it is left at its default value 0, the
   * minimal required number of Arnoldi vectors is reserved.
   */
  template <typename VECTOR, typename MATRIX, typename INVERSE>
  void solve(
    const MATRIX                       &B,
    const INVERSE                      &inverse,
    std::vector<std::complex<double> > &eigenvalues,
    std::vector<VECTOR>                &eigenvectors,
    const WhichEigenvalues              eigenvalue_of_interest = largest_magnitude,
    unsigned int                        number_of_arnoldi_vectors = 0);

protected:

  /**
   * Reference to the object that controls convergence of the iterative
   * solver.
   */
  SolverControl &solver_control;

private:

  /**
   * Exceptions.
   */
  DeclException2 (ExcInvalidNumberofEigenvalues, int, int,
                  << "Number of wanted eigenvalues " << arg1
                  << " is larger than the size of the matrix " << arg2);

  DeclException2 (ExcInvalidNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is larger than the size of the matrix " << arg2);

  DeclException2 (ExcSmallNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is too small to obtain " << arg2
                  << " eigenvalues");

  DeclException1 (ExcArpackIdo, int, << "This ido " << arg1
                  << " is not supported. Check documentation of ARPACK");

  DeclException1 (ExcArpackInfodsaupd, int,
                  << "Error with dsaupd, info " << arg1
                  << ". Check documentation of ARPACK");

  DeclException1 (ExcArpackInfodneupd, int,
                  << "Error with dneupd, info " << arg1
                  << ". Check documentation of ARPACK");

  DeclException1 (ExcArpackInfoMaxIt, int,
                  << "Maximum number " << arg1
                  << " of iterations reached.");

  DeclException1 (ExcArpackNoShifts, int,
                  << "No shifts could be applied during implicit"
                  << " Arnoldi update, try increasing the number of"
                  << " Arnoldi vectors.");
};


inline
ArpackSolver::ArpackSolver (SolverControl &control)
  :
  solver_control (control)
{}


template <typename VECTOR, typename MATRIX, typename INVERSE>
inline void
ArpackSolver::solve(const MATRIX                      &mass_matrix,
                    const INVERSE                     &inverse,
                    std::vector<std::complex<double>> &eigenvalues,
                    std::vector<VECTOR>               &eigenvectors,
                    const WhichEigenvalues             eigenvalue_of_interest,
                    unsigned int                       number_of_arnoldi_vectors)

{
  static GrowingVectorMemory<Vector<double> > vector_memory;

  // vector size
  const unsigned int n = eigenvectors[0].size();

  // number of eigenvalues to solve for
  const unsigned int n_eigenvalues = eigenvalues.size();

  Assert (n_eigenvalues < n,
          ExcInvalidNumberofEigenvalues(n_eigenvalues, n));

  Assert(eigenvectors.size() >= 2 * n_eigenvalues,
         ExcMessage("The vector eigenvectors must be at least twice the size "
                    "of the vector eigenvalues"));

  Assert(number_of_arnoldi_vectors < n,
         ExcInvalidNumberofArnoldiVectors(number_of_arnoldi_vectors, n));

  if (number_of_arnoldi_vectors == 0)
    number_of_arnoldi_vectors = 2 * n_eigenvalues + 1;

  Assert(number_of_arnoldi_vectors > 2 * n_eigenvalues + 1,
         ExcSmallNumberofArnoldiVectors(number_of_arnoldi_vectors, n_eigenvalues));

  //
  // Set up parameters and intermediate storage for the calls to dnaupd, dneupd:
  //

  // 'G' generalized eigenvalue problem, 'I' standard eigenvalue problem
  std::string bmat = "G";

  // specify the eigenvalues of interest
  std::string which;
  switch (eigenvalue_of_interest)
    {
    case algebraically_largest:
      which = "LA";
      break;
    case algebraically_smallest:
      which = "SA";
      break;
    case largest_magnitude:
      which = "LM";
      break;
    case smallest_magnitude:
      which = "SM";
      break;
    case largest_real_part:
      which = "LR";
      break;
    case smallest_real_part:
      which = "SR";
      break;
    case largest_imaginary_part:
      which = "LI";
      break;
    case smallest_imaginary_part:
      which = "SI";
      break;
    case both_ends:
      which = "BE";
      break;
    }

  std::vector<double> resid (n, 1.); // initialize with 1
  std::vector<double> v (n * number_of_arnoldi_vectors, .0);

  // information for arpack routines are passed with an int array:
  std::vector<int> iparam (11, 0);
  iparam[0] = 1;                     // shift strategy
  iparam[2] = control().max_steps(); // maximum number of iterations
  iparam[6] = 3; // sets the mode of dsaupd. 1 is exact shifting, 2 is user-supplied shifts
  // 3 is shift-invert mode, 4 is buckling mode, 5 is Cayley mode.

  std::vector<int> ipntr (14, 0);
  std::vector<double> workd (3 * n, .0);

  double tol = control().tolerance();

  int lworkl = 3 * number_of_arnoldi_vectors * (number_of_arnoldi_vectors + 6);
  std::vector<double> workl (lworkl, 0.);

  int ido = 0; // reverse communication parameter

  int info = 1; // information about the solver state

  // temporary storage
  VECTOR &src = *vector_memory.alloc();
  VECTOR &dst = *vector_memory.alloc();
  VECTOR &tmp = *vector_memory.alloc();
  VECTOR &tmp2 = *vector_memory.alloc();

  src.reinit(eigenvectors[0]);
  dst.reinit(src);
  tmp.reinit(src);
  tmp2.reinit(src);

  while (ido != 99)
    {
      // call to ARPACK dnaupd routine
      dnaupd_(&ido, bmat.c_str(), &n, which.c_str(), &n_eigenvalues, &tol,
              &resid[0], &number_of_arnoldi_vectors, &v[0], &n, &iparam[0], &ipntr[0],
              &workd[0], &workl[0], &lworkl, &info);

      if (ido == 99)
        break;

      if (ido == -1)
        {
          for (size_type i=0; i<src.size(); ++i)
            src(i) = workd[ipntr[0] - 1 + i];

          // multiplication with mass matrix M
          mass_matrix.vmult(tmp, src);
          // solving linear system
          inverse.vmult(dst, tmp);

          for (size_type i=0; i<dst.size(); ++i)
            workd[ipntr[1] - 1 + i] = dst(i);
        }

      else if (ido == 1)
        {
          for (size_type i=0; i<src.size(); ++i)
            {
              src(i) = workd[ipntr[2] - 1 + i];
              tmp(i) = workd[ipntr[0] - 1 + i];
            }

          // solving linear system
          inverse.vmult(dst,src);

          for (size_type i=0; i<dst.size(); ++i)
            workd[ipntr[1] - 1 + i] = dst(i);
        }

      else if (ido == 2)
        {
          for (size_type i=0; i<src.size(); ++i)
            src(i) = workd[ipntr[0] - 1 + i];

          // Multiplication with mass matrix M
          mass_matrix.vmult(dst, src);

          for (size_type i=0; i<dst.size(); ++i)
            workd[ipntr[1] - 1 + i] = dst(i);
        }

      else
        {
          AssertThrow (false, ExcArpackIdo(ido));
        }
    }

  vector_memory.free(&src);
  vector_memory.free(&dst);
  vector_memory.free(&tmp);
  vector_memory.free(&tmp2);

  // Check whether everything worked out as planned:
  AssertThrow (info >= 0, ExcArpackInfodsaupd(info));


  // 1 - compute eigenvectors, 0 - only eigenvalues
  int rvec = 1;

  // which eigenvectors
  char howmany = 'A';

  std::vector<int> select (number_of_arnoldi_vectors, 1);
  std::vector<double> z (n * number_of_arnoldi_vectors, 0.);

  double sigmar = 0.0; // real part of the shift
  double sigmai = 0.0; // imaginary part of the shift

  int lworkev = 3 * number_of_arnoldi_vectors;
  std::vector<double> workev (lworkev, 0.);

  std::vector<double> eigenvalues_real (n_eigenvalues, 0.);
  std::vector<double> eigenvalues_im (n_eigenvalues, 0.);

  // call to ARPACK dneupd routine
  dneupd_(&rvec, &howmany, &select[0], &eigenvalues_real[0],
          &eigenvalues_im[0], &z[0], &n, &sigmar, &sigmai,
          &workev[0], bmat.c_str(), &n, which.c_str(), &n_eigenvalues, &tol,
          &resid[0], &number_of_arnoldi_vectors, &v[0], &n,
          &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);

  AssertThrow (info != 1, ExcArpackInfoMaxIt(control().max_steps()));

  AssertThrow (info != 3, ExcArpackNoShifts(1));

  AssertThrow (info == 0, ExcArpackInfodneupd(info));

  for (size_type i = 0; i < eigenvectors.size(); ++i)
    for (unsigned int j=0; j<n; ++j)
      eigenvectors[i](j) = z[i*n+j];

  for (size_type i = 0; i < eigenvalues.size(); ++i)
    eigenvalues[i] = std::complex<double> (eigenvalues_real[i],
                                           eigenvalues_im[i]);
}


inline
SolverControl &ArpackSolver::control () const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE

#endif /* DEAL_II_WITH_ARPACK */
#endif
