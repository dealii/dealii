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

#ifndef dealii__parpack_solver_h
#define dealii__parpack_solver_h

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/vector.h>
#include <deal.II/base/index_set.h>

#include <cstring>


#ifdef DEAL_II_ARPACK_WITH_PARPACK

DEAL_II_NAMESPACE_OPEN

extern "C" {

  // http://www.mathkeisan.com/usersguide/man/pdnaupd.html
  void pdnaupd_(MPI_Fint *comm, int *ido, char *bmat, int *n, char *which,
                int *nev, double *tol, double *resid, int *ncv,
                double *v, int *nloc, int *iparam, int *ipntr,
                double *workd, double *workl, int *lworkl,
                int *info);

  // http://www.mathkeisan.com/usersguide/man/pdsaupd.html
  void pdsaupd_(MPI_Fint *comm, int *ido, char *bmat, int *n, char *which,
                int *nev, double *tol, double *resid, int *ncv,
                double *v, int *nloc, int *iparam, int *ipntr,
                double *workd, double *workl, int *lworkl,
                int *info);

  // http://www.mathkeisan.com/usersguide/man/pdneupd.html
  void pdneupd_(MPI_Fint *comm, int *rvec, char *howmany, int *select, double *d,
                double *di, double *z, int *ldz, double *sigmar,
                double *sigmai, double *workev, char *bmat, int *n, char *which,
                int *nev, double *tol, double *resid, int *ncv,
                double *v, int *nloc, int *iparam, int *ipntr,
                double *workd, double *workl, int *lworkl, int *info);

  // http://www.mathkeisan.com/usersguide/man/pdseupd.html
  void pdseupd_(MPI_Fint *comm, int *rvec, char *howmany, int *select, double *d,
                double *z, int *ldz, double *sigmar,
                char *bmat, int *n, char *which,
                int *nev, double *tol, double *resid, int *ncv,
                double *v, int *nloc, int *iparam, int *ipntr,
                double *workd, double *workl, int *lworkl, int *info);

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
 * computed. <code>OP</code> is an inverse operation for the matrix <code>A -
 * sigma * B</code>, where <code> sigma </code> is a shift value, set to zero
 * by default.
 *
 * Through the AdditionalData the user can specify some of the parameters to
 * be set.
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
 * @author Denis Davydov, 2015.
 */
template <typename VectorType>
class PArpackSolver : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * An enum that lists the possible choices for which eigenvalues to compute
   * in the solve() function.
   *
   * A particular choice is limited based on symmetric or non-symmetric matrix
   * <code>A</code> considered.
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
   * Auxiliary class to represent <code>A-sigma*B</code> operator.
   */
  template <typename MatrixType>
  class Shift : public dealii::Subscriptor
  {
  public:

    /**
     * Constructor.
     */
    Shift (const MatrixType &A,
           const MatrixType &B,
           const double      sigma)
      :
      A(A),
      B(B),
      sigma(sigma)
    {}

    /**
     * Apply <code>A-sigma * B</code>
     */
    void vmult (VectorType &dst, const VectorType &src) const
    {
      B.vmult(dst,src);
      dst *= (-sigma);
      A.vmult_add(dst,src);
    }

    /**
     * Apply <code>A^T-sigma * B^T</code>
     */
    void Tvmult (VectorType &dst, const VectorType &src) const
    {
      B.Tvmult(dst,src);
      dst *= (-sigma);
      A.Tvmult_add(dst,src);
    }

  private:
    const MatrixType &A;
    const MatrixType &B;
    const double sigma;
  };

  /**
   * Standardized data struct to pipe additional data to the solver, should it
   * be needed.
   */
  struct AdditionalData
  {
    const unsigned int number_of_arnoldi_vectors;
    const WhichEigenvalues eigenvalue_of_interest;
    const bool symmetric;
    AdditionalData(
      const unsigned int number_of_arnoldi_vectors = 15,
      const WhichEigenvalues eigenvalue_of_interest = largest_magnitude,
      const bool symmetric = false);
  };

  /**
   * Access to the object that controls convergence.
   */
  SolverControl &control () const;

  /**
   * Constructor.
   */
  PArpackSolver(SolverControl &control,
                const MPI_Comm &mpi_communicator,
                const AdditionalData &data = AdditionalData());

  /**
   * Initialise internal variables.
   */
  void reinit(const dealii::IndexSet &locally_owned_dofs );

  /**
   * Set desired shift value.
   */
  void set_shift(const double s );

  /**
   * Solve the generalized eigensprectrum problem $A x=\lambda B x$ by calling
   * the <code>pd(n/s)eupd</code> and <code>pd(n/s)aupd</code> functions of
   * PARPACK.
   */
  template <typename MatrixType1,
            typename MatrixType2, typename INVERSE>
  void solve
  (const MatrixType1                  &A,
   const MatrixType2                  &B,
   const INVERSE                      &inverse,
   std::vector<std::complex<double> > &eigenvalues,
   std::vector<VectorType>            &eigenvectors,
   const unsigned int                  n_eigenvalues);

  std::size_t memory_consumption() const;

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
  VectorType src,dst,tmp;

  /**
   * Indices of local degrees of freedom.
   */
  std::vector< types::global_dof_index > local_indices;

  /**
   * The shift value to be applied during solution
   */
  double shift_value;

private:

  /**
   * PArpackExcInfoPdnaupds.
   */
  DeclException2 (PArpackExcConvergedEigenvectors, int, int,
                  << arg1 << "eigenpairs were requested, but only"
                  << arg2 << " converged");

  DeclException2 (PArpackExcInvalidNumberofEigenvalues, int, int,
                  << "Number of wanted eigenvalues " << arg1
                  << " is larger that the size of the matrix " << arg2);

  DeclException2 (PArpackExcInvalidEigenvectorSize, int, int,
                  << "Number of wanted eigenvalues " << arg1
                  << " is larger that the size of eigenvectors " << arg2);

  DeclException2 (PArpackExcInvalidEigenvalueSize, int, int,
                  << "Number of wanted eigenvalues " << arg1
                  << " is larger that the size of eigenvalues " << arg2);

  DeclException2 (PArpackExcInvalidNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is larger that the size of the matrix " << arg2);

  DeclException2 (PArpackExcSmallNumberofArnoldiVectors, int, int,
                  << "Number of Arnoldi vectors " << arg1
                  << " is too small to obtain " << arg2
                  << " eigenvalues");

  DeclException1 (PArpackExcIdo, int, << "This ido " << arg1
                  << " is not supported. Check documentation of ARPACK");

  DeclException1 (PArpackExcMode, int, << "This mode " << arg1
                  << " is not supported. Check documentation of ARPACK");

  DeclException1 (PArpackExcInfoPdnaupd, int,
                  << "Error with Pdnaupd, info " << arg1
                  << ". Check documentation of ARPACK");

  DeclException1 (PArpackExcInfoPdneupd, int,
                  << "Error with Pdneupd, info " << arg1
                  << ". Check documentation of ARPACK");

  DeclException1 (PArpackExcInfoMaxIt, int,
                  << "Maximum number " << arg1
                  << " of iterations reached.");

  DeclException1 (PArpackExcNoShifts, int,
                  << "No shifts could be applied during implicit"
                  << " Arnoldi update, try increasing the number of"
                  << " Arnoldi vectors.");
};

template <typename VectorType>
std::size_t
PArpackSolver<VectorType>::memory_consumption() const
{
  return  MemoryConsumption::memory_consumption (double()) *
          (workl.size()  +
           workd.size()  +
           v.size()      +
           resid.size()  +
           z.size()      +
           workev.size()  )       +
          src.memory_consumption() +
          dst.memory_consumption() +
          tmp.memory_consumption() +
          MemoryConsumption::memory_consumption (types::global_dof_index()) * local_indices.size();
}

template <typename VectorType>
PArpackSolver<VectorType>::AdditionalData::
AdditionalData (const unsigned int     number_of_arnoldi_vectors,
                const WhichEigenvalues eigenvalue_of_interest,
                const bool             symmetric)
  :
  number_of_arnoldi_vectors(number_of_arnoldi_vectors),
  eigenvalue_of_interest(eigenvalue_of_interest),
  symmetric(symmetric)
{}

template <typename VectorType>
PArpackSolver<VectorType>::PArpackSolver (SolverControl        &control,
                                          const MPI_Comm       &mpi_communicator,
                                          const AdditionalData &data)
  :
  solver_control (control),
  additional_data (data),
  mpi_communicator( mpi_communicator ),
  mpi_communicator_fortran ( MPI_Comm_c2f( mpi_communicator ) ),
  shift_value(0.0)

{}

template <typename VectorType>
void PArpackSolver<VectorType>::set_shift(const double s )
{
  shift_value = s;
}

template <typename VectorType>
void PArpackSolver<VectorType>::reinit(const dealii::IndexSet &locally_owned_dofs)
{
  // store local indices to write to vectors
  locally_owned_dofs.fill_index_vector(local_indices);

  // scalars
  nloc = locally_owned_dofs.n_elements ();
  ncv  = additional_data.number_of_arnoldi_vectors;

  Assert ((int)local_indices.size() == nloc, ExcInternalError() );

  // vectors
  ldv = nloc;
  v.resize (ldv*ncv, 0.0);

  // TODO: add optional input for resid
  resid.resize(nloc, 1.0);

  // work arrays for ARPACK
  workd.resize(3*nloc,0.0);

  lworkl = additional_data.symmetric ?
           ncv*ncv + 8*ncv
           :
           3*ncv*ncv+6*ncv;
  workl.resize (lworkl, 0.);

  ldz = nloc;
  z.resize (ldz*ncv, 0.); // TODO we actually need only ldz*nev

  // WORKEV  Double precision  work array of dimension 3*NCV.
  lworkev = additional_data.symmetric ?
            0 /*not used in symmetric case*/
            :
            3*ncv;
  workev.resize (lworkev, 0.);

  select.resize (ncv, 0);

  // deal.II vectors:
  src.reinit (locally_owned_dofs,mpi_communicator);
  dst.reinit (locally_owned_dofs,mpi_communicator);
  tmp.reinit (locally_owned_dofs,mpi_communicator);

}

template <typename VectorType>
template <typename MatrixType1,typename MatrixType2, typename INVERSE>
void PArpackSolver<VectorType>::solve
(const MatrixType1                  &/*system_matrix*/,
 const MatrixType2                  &mass_matrix,
 const INVERSE                      &inverse,
 std::vector<std::complex<double> > &eigenvalues,
 std::vector<VectorType>            &eigenvectors,
 const unsigned int                  n_eigenvalues)
{

  Assert (n_eigenvalues <= eigenvectors.size(),
          PArpackExcInvalidEigenvectorSize(n_eigenvalues, eigenvectors.size()));

  Assert (n_eigenvalues <= eigenvalues.size(),
          PArpackExcInvalidEigenvalueSize(n_eigenvalues, eigenvalues.size()));


  Assert (n_eigenvalues < mass_matrix.m(),
          PArpackExcInvalidNumberofEigenvalues(n_eigenvalues, mass_matrix.m()));

  Assert (additional_data.number_of_arnoldi_vectors < mass_matrix.m(),
          PArpackExcInvalidNumberofArnoldiVectors(
            additional_data.number_of_arnoldi_vectors, mass_matrix.m()));

  Assert (additional_data.number_of_arnoldi_vectors > 2*n_eigenvalues+1,
          PArpackExcSmallNumberofArnoldiVectors(
            additional_data.number_of_arnoldi_vectors, n_eigenvalues));
  // ARPACK mode for dnaupd, here only
  //  Mode 3:  K*x = lambda*M*x, K symmetric, M symmetric positive semi-definite
  //c           ===> OP = (inv[K - sigma*M])*M  and  B = M.
  //c           ===> Shift-and-Invert mode
  int mode = 3;

  // reverse communication parameter
  // must be zero on the first call to pdnaupd
  int ido = 0;

  // 'G' generalized eigenvalue problem
  // 'I' standard eigenvalue problem
  char bmat[2] = "G";

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
  double tol = control().tolerance();

  //information to the routines
  std::vector<int> iparam (11, 0);

  iparam[0] = 1;
  // shift strategy: exact shifts with respect to the current Hessenberg matrix H.

  // maximum number of iterations
  iparam[2] = control().max_steps();

  // Parpack currently works only for NB = 1
  iparam[3] = 1;

  // Sets the mode of dsaupd:
  // 1 is exact shifting,
  // 2 is user-supplied shifts,
  // 3 is shift-invert mode,
  // 4 is buckling mode,
  // 5 is Cayley mode.

  iparam[6] = mode;
  std::vector<int> ipntr (14, 0);

  //information out of the iteration
  //  If INFO .EQ. 0, a random initial residual vector is used.
  //  If INFO .NE. 0, RESID contains the initial residual vector,
  //  possibly from a previous run.
  // Typical choices in this situation might be to use the final value
  // of the starting vector from the previous eigenvalue calculation
  int info = 1;

  // Number of eigenvalues of OP to be computed. 0 < NEV < N.
  int nev = n_eigenvalues;
  int n_inside_arpack = nloc;

  while (ido != 99)
    {
      // call of ARPACK pdnaupd routine
      if (additional_data.symmetric)
        pdsaupd_(&mpi_communicator_fortran,&ido, bmat, &n_inside_arpack, which, &nev, &tol,
                 &resid[0], &ncv, &v[0], &ldv, &iparam[0], &ipntr[0],
                 &workd[0], &workl[0], &lworkl, &info);
      else
        pdnaupd_(&mpi_communicator_fortran,&ido, bmat, &n_inside_arpack, which, &nev, &tol,
                 &resid[0], &ncv, &v[0], &ldv, &iparam[0], &ipntr[0],
                 &workd[0], &workl[0], &lworkl, &info);

      if (ido == 99)
        break;

      switch (mode)
        {
//        OP = (inv[K - sigma*M])*M
        case 3:
        {
          switch (ido)
            {
//            compute  Y = OP * X  where
//            IPNTR(1) is the pointer into WORKD for X,
//            IPNTR(2) is the pointer into WORKD for Y.
            case -1:
            {
              const int shift_x = ipntr[0]-1;
              const int shift_y = ipntr[1]-1;
              Assert (shift_x>=0, dealii::ExcInternalError() );
              Assert (shift_x+nloc <= (int)workd.size(), dealii::ExcInternalError() );
              Assert (shift_y>=0, dealii::ExcInternalError() );
              Assert (shift_y+nloc <= (int)workd.size(), dealii::ExcInternalError() );

              src = 0.0;
              src.add (nloc,
                       &local_indices[0],
                       &workd[0]+shift_x );
              src.compress (VectorOperation::add);

              // multiplication with mass matrix M
              mass_matrix.vmult(tmp, src);
              // solving linear system
              inverse.vmult(dst,tmp);

              // store the result
              dst.extract_subvector_to (local_indices.begin(),
                                        local_indices.end(),
                                        &workd[0]+shift_y  );
            }
            break;

//            compute  Y = OP * X where
//            IPNTR(1) is the pointer into WORKD for X,
//            IPNTR(2) is the pointer into WORKD for Y.
//            In mode 3,4 and 5, the vector B * X is already
//            available in WORKD(ipntr(3)).  It does not
//            need to be recomputed in forming OP * X.
            case  1:
            {
              const int shift_x   = ipntr[0]-1;
              const int shift_y   = ipntr[1]-1;
              const int shift_b_x = ipntr[2]-1;

              Assert (shift_x>=0, dealii::ExcInternalError() );
              Assert (shift_x+nloc <= (int)workd.size(), dealii::ExcInternalError() );
              Assert (shift_y>=0, dealii::ExcInternalError() );
              Assert (shift_y+nloc <= (int)workd.size(), dealii::ExcInternalError() );
              Assert (shift_b_x>=0, dealii::ExcInternalError() );
              Assert (shift_b_x+nloc <= (int)workd.size(), dealii::ExcInternalError() );
              Assert (shift_y>=0, dealii::ExcInternalError() );
              Assert (shift_y+nloc <= (int)workd.size(), dealii::ExcInternalError() );

              src = 0.0; // B*X
              src.add (nloc,
                       &local_indices[0],
                       &workd[0]+shift_b_x );

              tmp = 0.0; // X
              tmp.add (nloc,
                       &local_indices[0],
                       &workd[0]+shift_x);

              src.compress (VectorOperation::add);
              tmp.compress (VectorOperation::add);

              // solving linear system
              inverse.vmult(dst,src);

              // store the result
              dst.extract_subvector_to (local_indices.begin(),
                                        local_indices.end(),
                                        &workd[0]+shift_y  );

            }
            break;

//            compute  Y = B * X  where
//            IPNTR(1) is the pointer into WORKD for X,
//            IPNTR(2) is the pointer into WORKD for Y.
            case  2:
            {

              const int shift_x = ipntr[0]-1;
              const int shift_y = ipntr[1]-1;
              Assert (shift_x>=0, dealii::ExcInternalError() );
              Assert (shift_x+nloc <= (int)workd.size(), dealii::ExcInternalError() );
              Assert (shift_y>=0, dealii::ExcInternalError() );
              Assert (shift_y+nloc <= (int)workd.size(), dealii::ExcInternalError() );

              src = 0.0;
              src.add (nloc,
                       &local_indices[0],
                       &workd[0]+shift_x );
              src.compress (VectorOperation::add);

              // Multiplication with mass matrix M
              mass_matrix.vmult(dst, src);

              // store the result
              dst.extract_subvector_to (local_indices.begin(),
                                        local_indices.end(),
                                        &workd[0]+shift_y);

            }
            break;

            default:
              Assert (false, PArpackExcIdo(ido));
              break;
            }
        }
        break;
        default:
          Assert (false, PArpackExcMode(mode));
          break;
        }
    }

  if (info<0)
    {
      Assert (false, PArpackExcInfoPdnaupd(info));
    }
  else
    {
      // 1 - compute eigenvectors,
      // 0 - only eigenvalues
      int rvec = 1;

      // which eigenvectors
      char howmany[4] = "All";

      double sigmar = shift_value; // real part of the shift
      double sigmai = 0.0; // imaginary part of the shift

      std::vector<double> eigenvalues_real (n_eigenvalues, 0.);
      std::vector<double> eigenvalues_im (n_eigenvalues, 0.);

      // call of ARPACK pdneupd routine
      if (additional_data.symmetric)
        pdseupd_(&mpi_communicator_fortran, &rvec, howmany, &select[0], &eigenvalues_real[0],
                 &z[0], &ldz, &sigmar,
                 bmat, &n_inside_arpack, which, &nev, &tol,
                 &resid[0], &ncv, &v[0], &ldv,
                 &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);
      else
        pdneupd_(&mpi_communicator_fortran, &rvec, howmany, &select[0], &eigenvalues_real[0],
                 &eigenvalues_im[0], &z[0], &ldz, &sigmar, &sigmai,
                 &workev[0], bmat, &n_inside_arpack, which, &nev, &tol,
                 &resid[0], &ncv, &v[0], &ldv,
                 &iparam[0], &ipntr[0], &workd[0], &workl[0], &lworkl, &info);

      if (info == 1)
        {
          Assert (false, PArpackExcInfoMaxIt(control().max_steps()));
        }
      else if (info == 3)
        {
          Assert (false, PArpackExcNoShifts(1));
        }
      else if (info!=0)
        {
          Assert (false, PArpackExcInfoPdneupd(info));
        }

      for (size_type i=0; i<n_eigenvalues; ++i)
        {
          eigenvectors[i] = 0.0;
          Assert (i*nloc + nloc <= v.size(), dealii::ExcInternalError() );

          eigenvectors[i].add (nloc,
                               &local_indices[0],
                               &v[i*nloc] );
          eigenvectors[i].compress (VectorOperation::add);
        }

      for (size_type i=0; i<n_eigenvalues; ++i)
        eigenvalues[i] = std::complex<double> (eigenvalues_real[i],
                                               eigenvalues_im[i]);
    }

  Assert (iparam[4] == (int)n_eigenvalues,
          PArpackExcConvergedEigenvectors(iparam[4], n_eigenvalues));

  // both PDNAUPD and PDSAUPD compute eigenpairs of inv[A - sigma*M]*M
  // with respect to a semi-inner product defined by M.

  // resid likely contains residual with respect to M-norm.
  {

    tmp = 0.0;
    tmp.add (nloc,
             &local_indices[0],
             &resid[0]);
    solver_control.check  ( iparam[2], tmp.l2_norm() );
  }


}

template <typename VectorType>
SolverControl &PArpackSolver<VectorType>::control () const
{
  return solver_control;
}

DEAL_II_NAMESPACE_CLOSE


#endif
#endif
