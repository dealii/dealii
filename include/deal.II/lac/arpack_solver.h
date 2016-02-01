// ---------------------------------------------------------------------
//
// Copyright (C) 2010 - 2016 by the deal.II authors
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

#ifndef dealii__arpack_solver_h
#define dealii__arpack_solver_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/solver_control.h>

#include <cstring>


#ifdef DEAL_II_WITH_ARPACK

DEAL_II_NAMESPACE_OPEN


extern "C" void dnaupd_(int *ido, char *bmat, const unsigned int *n, char *which,
                        const unsigned int *nev, const double *tol, double *resid, int *ncv,
                        double *v, int *ldv, int *iparam, int *ipntr,
                        double *workd, double *workl, int *lworkl,
                        int *info);

extern "C" void dneupd_(int *rvec, char *howmany, int *select, double *d,
                        double *di, double *z, int *ldz, double *sigmar,
                        double *sigmai, double *workev, char *bmat,const unsigned int *n, char *which,
                        const unsigned int *nev, const double *tol, double *resid, int *ncv,
                        double *v, int *ldv, int *iparam, int *ipntr,
                        double *workd, double *workl, int *lworkl, int *info);

/**
 * Interface for using ARPACK. ARPACK is a collection of Fortran77 subroutines
 * designed to solve large scale eigenvalue problems.  Here we interface to
 * the routines <code>dneupd</code> and <code>dnaupd</code> of ARPACK.  The
 * package is designed to compute a few eigenvalues and corresponding
 * eigenvectors of a general n by n matrix A. It is most appropriate for large
 * sparse matrices A.
 *
 * In this class we make use of the method applied to the generalized
 * eigenspectrum problem $(A-\lambda B)x=0$, for $x\neq0$; where $A$ is a
 * system matrix, $B$ is a mass matrix, and $\lambda, x$ are a set of
 * eigenvalues and eigenvectors respectively.
 *
 * The ArpackSolver can be used in application codes with serial objects in
 * the following way:
 * @code
 * SolverControl solver_control (1000, 1e-9);
 * ArpackSolver (solver_control);
 * system.solve (A, B, OP, lambda, x, size_of_spectrum);
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
 * For further information on how the ARPACK routines <code>dneupd</code> and
 * <code>dnaupd</code> work and also how to set the parameters appropriately
 * please take a look into the ARPACK manual.
 *
 * @note Whenever you eliminate degrees of freedom using ConstraintMatrix, you
 * generate spurious eigenvalues and eigenvectors. If you make sure that the
 * diagonals of eliminated matrix rows are all equal to one, you get a single
 * additional eigenvalue. But beware that some functions in deal.II set these
 * diagonals to rather arbitrary (from the point of view of eigenvalue
 * problems) values. See also
 * @ref step_36 "step-36"
 * for an example.
 *
 * @author Baerbel Janssen, Agnieszka Miedlar, 2010, Guido Kanschat 2015
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
   * Standardized data struct to pipe additional data to the solver, should it
   * be needed.
   */
  struct AdditionalData
  {
    const unsigned int number_of_arnoldi_vectors;
    const WhichEigenvalues eigenvalue_of_interest;
    AdditionalData(
      const unsigned int number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest = largest_magnitude);
  };

  /**
   * Access to the object that controls convergence.
   */
  SolverControl &control () const;

  /**
   * Constructor.
   */
  ArpackSolver(SolverControl &control,
               const AdditionalData &data = AdditionalData());

  /**
   * Solve the generalized eigensprectrum problem $A x=\lambda B x$ by calling
   * the <code>dneupd</code> and <code>dnaupd</code> functions of ARPACK.
   *
   * The function returns a vector of eigenvalues of length <i>n</i> and a
   * vector of eigenvectors, where the latter should be twice the size of the
   * eigenvalue vector. The first <i>n</i> vectors in
   * <code>eigenvectors</code> will be the real parts of the eigenvectors, the
   * second <i>n</i> the imaginary parts.
   *
   * @param A The operator for which we want to compute eigenvalues. Actually,
   * this parameter is entirely unused.
   *
   * @param B The inner product of the underlying space, typically the mass
   * matrix. For constrained problems, it can be a partial mass matrix, like
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
   * alternatingly the real parts and the imaginary parts of the eigenvectors.
   * Therefore, its length should be twice the number of eigenvalues. The
   * vectors have to be initialized to match the matrices.
   *
   * @param n_eigenvalues The purpose of this parameter is not clear, but it
   * is safe to set it to the size of <code>eigenvalues</code> or greater.
   * Leave it at its default zero, which will be reset to the size of
   * <code>eigenvalues</code> internally.
   */
  template <typename VectorType, typename MatrixType1,
            typename MatrixType2, typename INVERSE>
  void solve (const MatrixType1                  &A,
              const MatrixType2                  &B,
              const INVERSE                      &inverse,
              std::vector<std::complex<double> > &eigenvalues,
              std::vector<VectorType>            &eigenvectors,
              const unsigned int                  n_eigenvalues = 0);

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

private:

  /**
   * Exceptions.
   */
  DeclException2 (ExcInvalidNumberofEigenvalues, int, int,
                  << "Number of wanted eigenvalues " << arg1
                  << " is larger that the size of the matrix " << arg2);

  DeclException2 (ExcInvalidNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is larger that the size of the matrix " << arg2);

  DeclException2 (ExcSmallNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is too small to obtain " << arg2
                  << " eigenvalues");

  DeclException1 (ExcArpackIdo, int, << "This ido " << arg1
                  << " is not supported. Check documentation of ARPACK");

  DeclException1 (ExcArpackMode, int, << "This mode " << arg1
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

  DeclExceptionMsg (ExcArpackNoShifts,
                    "No shifts could be applied during implicit"
                    " Arnoldi update, try increasing the number of"
                    " Arnoldi vectors.");
};


inline
ArpackSolver::AdditionalData::
AdditionalData (const unsigned int number_of_arnoldi_vectors,
                const WhichEigenvalues eigenvalue_of_interest)
  :
  number_of_arnoldi_vectors(number_of_arnoldi_vectors),
  eigenvalue_of_interest(eigenvalue_of_interest)
{}


inline
ArpackSolver::ArpackSolver (SolverControl &control,
                            const AdditionalData &data)
  :
  solver_control (control),
  additional_data (data)

{}


template <typename VectorType, typename MatrixType1,
          typename MatrixType2, typename INVERSE>
inline
void ArpackSolver::solve (const MatrixType1                  &/*system_matrix*/,
                          const MatrixType2                  &mass_matrix,
                          const INVERSE                      &inverse,
                          std::vector<std::complex<double> > &eigenvalues,
                          std::vector<VectorType>            &eigenvectors,
                          const unsigned int                  n_eigenvalues)
{
  //inside the routines of ARPACK the
  //values change magically, so store
  //them here

  const unsigned int n = eigenvectors[0].size();
  const unsigned int n_inside_arpack = eigenvectors[0].size();

  // Number of eigenvalues for arpack
  const unsigned int nev = (n_eigenvalues == 0) ? eigenvalues.size() : n_eigenvalues;
  AssertIndexRange(eigenvalues.size()-1, nev);
  /*
  if(n < 0 || nev <0 || p < 0 || maxit < 0 )
       std:cout << "All input parameters have to be positive.\n";
       */
  Assert (n_eigenvalues < n,
          ExcInvalidNumberofEigenvalues(nev, n));

  Assert (additional_data.number_of_arnoldi_vectors < n,
          ExcInvalidNumberofArnoldiVectors(
            additional_data.number_of_arnoldi_vectors, n));

  Assert (additional_data.number_of_arnoldi_vectors > 2*nev+1,
          ExcSmallNumberofArnoldiVectors(
            additional_data.number_of_arnoldi_vectors, nev));
  // ARPACK mode for dnaupd, here only mode 3
  int mode = 3;

  // reverse communication parameter
  int ido = 0;

  /**
   * 'G' generalized eigenvalue problem 'I' standard eigenvalue problem
   */
  char bmat[2] = "G";

  /**
   * Specify the eigenvalues of interest, possible parameters "LA"
   * algebraically largest "SA" algebraically smallest "LM" largest magnitude
   * "SM" smallest magnitude "LR" largest real part "SR" smallest real part
   * "LI" largest imaginary part "SI" smallest imaginary part "BE" both ends
   * of spectrum simultaneous
   */
  char which[3];
  switch (additional_data.eigenvalue_of_interest)
    {
    case algebraically_largest:
      std::strcpy (which, "LA");
      break;
    case algebraically_smallest:
      std::strcpy (which, "SA");
      break;
    case largest_magnitude:
      std::strcpy (which, "LM");
      break;
    case smallest_magnitude:
      std::strcpy (which, "SM");
      break;
    case largest_real_part:
      std::strcpy (which, "LR");
      break;
    case smallest_real_part:
      std::strcpy (which, "SR");
      break;
    case largest_imaginary_part:
      std::strcpy (which, "LI");
      break;
    case smallest_imaginary_part:
      std::strcpy (which, "SI");
      break;
    case both_ends:
      std::strcpy (which, "BE");
      break;
    }

  // tolerance for ARPACK
  const double tol = control().tolerance();

  // if the starting vector is used it has to be in resid
  std::vector<double> resid(n, 1.);

  // number of Arnoldi basis vectors specified
  // in additional_data
  int ncv = additional_data.number_of_arnoldi_vectors;

  int ldv = n;
  std::vector<double> v (ldv*ncv, 0.0);

  //information to the routines
  std::vector<int> iparam (11, 0);

  iparam[0] = 1;        // shift strategy

  // maximum number of iterations
  iparam[2] = control().max_steps();

  /**
   * Sets the mode of dsaupd. 1 is exact shifting, 2 is user-supplied shifts,
   * 3 is shift-invert mode, 4 is buckling mode, 5 is Cayley mode.
   */

  iparam[6] = mode;
  std::vector<int> ipntr (14, 0);

  // work arrays for ARPACK
  double *workd;
  workd = new double[3*n];

  for (unsigned int i=0; i<3*n; ++i)
    workd[i] = 0.0;

  int lworkl = 3*ncv*(ncv + 6);
  std::vector<double> workl (lworkl, 0.);
  //information out of the iteration
  int info = 1;

  while (ido != 99)
    {
      // call of ARPACK dnaupd routine
      dnaupd_(&ido, bmat, &n_inside_arpack, which, &nev, &tol,
              &resid[0], &ncv, &v[0], &ldv, &iparam[0], &ipntr[0],
              workd, &workl[0], &lworkl, &info);

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

              VectorType src,dst,tmp;
              src.reinit(eigenvectors[0]);
              dst.reinit(src);
              tmp.reinit(src);


              for (size_type i=0; i<src.size(); ++i)
                src(i) = *(workd+ipntr[0]-1+i);

              // multiplication with mass matrix M
              mass_matrix.vmult(tmp, src);
              // solving linear system
              inverse.vmult(dst,tmp);

              for (size_type i=0; i<dst.size(); ++i)
                *(workd+ipntr[1]-1+i) = dst(i);
            }
            break;

            case  1:
            {

              VectorType src,dst,tmp, tmp2;
              src.reinit(eigenvectors[0]);
              dst.reinit(src);
              tmp.reinit(src);
              tmp2.reinit(src);

              for (size_type i=0; i<src.size(); ++i)
                {
                  src(i) = *(workd+ipntr[2]-1+i);
                  tmp(i) = *(workd+ipntr[0]-1+i);
                }
              // solving linear system
              inverse.vmult(dst,src);

              for (size_type i=0; i<dst.size(); ++i)
                *(workd+ipntr[1]-1+i) = dst(i);
            }
            break;

            case  2:
            {

              VectorType src,dst;
              src.reinit(eigenvectors[0]);
              dst.reinit(src);

              for (size_type i=0; i<src.size(); ++i)
                src(i) = *(workd+ipntr[0]-1+i);

              // Multiplication with mass matrix M
              mass_matrix.vmult(dst, src);

              for (size_type i=0; i<dst.size(); ++i)
                *(workd+ipntr[1]-1+i) = dst(i);

            }
            break;

            default:
              Assert (false, ExcArpackIdo(ido));
              break;
            }
        }
        break;
        default:
          Assert (false, ExcArpackMode(mode));
          break;
        }
    }

  if (info<0)
    {
      Assert (false, ExcArpackInfodsaupd(info));
    }
  else
    {
      /**
       * 1 - compute eigenvectors, 0 - only eigenvalues
       */
      int rvec = 1;

      // which eigenvectors
      char howmany = 'A';

      std::vector<int> select (ncv, 1);

      int ldz = n;

      std::vector<double> z (ldz*ncv, 0.);

      double sigmar = 0.0; // real part of the shift
      double sigmai = 0.0; // imaginary part of the shift

      int lworkev = 3*ncv;
      std::vector<double> workev (lworkev, 0.);

      std::vector<double> eigenvalues_real (nev, 0.);
      std::vector<double> eigenvalues_im (nev, 0.);

      // call of ARPACK dneupd routine
      dneupd_(&rvec, &howmany, &select[0], &eigenvalues_real[0],
              &eigenvalues_im[0], &z[0], &ldz, &sigmar, &sigmai,
              &workev[0], bmat, &n_inside_arpack, which, &nev, &tol,
              &resid[0], &ncv, &v[0], &ldv,
              &iparam[0], &ipntr[0], workd, &workl[0], &lworkl, &info);

      if (info == 1)
        {
          Assert (false, ExcArpackInfoMaxIt(control().max_steps()));
        }
      else if (info == 3)
        {
          Assert (false, ExcArpackNoShifts());
        }
      else if (info!=0)
        {
          Assert (false, ExcArpackInfodneupd(info));
        }


      const unsigned int n_eigenvecs = eigenvectors.size();
      for (size_type i=0; i<n_eigenvecs; ++i)
        for (unsigned int j=0; j<n; ++j)
          eigenvectors[i](j) = v[i*n+j];

      delete[] workd;

      AssertDimension (eigenvalues.size(), eigenvalues_real.size());
      AssertDimension (eigenvalues.size(), eigenvalues_im.size());

      for (size_type i=0; i<eigenvalues.size(); ++i)
        eigenvalues[i] = std::complex<double> (eigenvalues_real[i],
                                               eigenvalues_im[i]);
    }
}


inline
SolverControl &ArpackSolver::control () const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#endif
#endif
