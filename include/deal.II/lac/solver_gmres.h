// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2016 by the deal.II authors
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

#ifndef dealii__solver_gmres_h
#define dealii__solver_gmres_h



#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/logstream.h>
#include <deal.II/lac/householder.h>
#include <deal.II/lac/solver.h>
#include <deal.II/lac/solver_control.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>
#include <cmath>
#include <algorithm>

DEAL_II_NAMESPACE_OPEN

/*!@addtogroup Solvers */
/*@{*/

namespace internal
{
  /**
   * A namespace for a helper class to the GMRES solver.
   */
  namespace SolverGMRES
  {
    /**
     * Class to hold temporary vectors.  This class automatically allocates a
     * new vector, once it is needed.
     *
     * A future version should also be able to shift through vectors
     * automatically, avoiding restart.
     */

    template <typename VectorType>
    class TmpVectors
    {
    public:
      /**
       * Constructor. Prepares an array of @p VectorType of length @p
       * max_size.
       */
      TmpVectors(const unsigned int       max_size,
                 VectorMemory<VectorType> &vmem);

      /**
       * Delete all allocated vectors.
       */
      ~TmpVectors();

      /**
       * Get vector number @p i. If this vector was unused before, an error
       * occurs.
       */
      VectorType &operator[] (const unsigned int i) const;

      /**
       * Get vector number @p i. Allocate it if necessary.
       *
       * If a vector must be allocated, @p temp is used to reinit it to the
       * proper dimensions.
       */
      VectorType &operator() (const unsigned int i,
                              const VectorType   &temp);

      /**
       * Returns size of data vector. It is used in the solver to store
       * the Arnoldi vectors.
       */
      unsigned int size() const;


    private:
      /**
       * Pool were vectors are obtained from.
       */
      VectorMemory<VectorType> &mem;

      /**
       * Field for storing the vectors.
       */
      std::vector<VectorType *> data;

      /**
       * Offset of the first vector. This is for later when vector rotation
       * will be implemented.
       */
      unsigned int offset;
    };
  }
}

/**
 * Implementation of the Restarted Preconditioned Direct Generalized Minimal
 * Residual Method. The stopping criterion is the norm of the residual.
 *
 * The AdditionalData structure contains the number of temporary vectors used.
 * The size of the Arnoldi basis is this number minus three. Additionally, it
 * allows you to choose between right or left preconditioning. The default is
 * left preconditioning. Finally it includes a flag indicating whether or not
 * the default residual is used as stopping criterion.
 *
 *
 * <h3>Left versus right preconditioning</h3>
 *
 * @p AdditionalData allows you to choose between left and right
 * preconditioning. As expected, this switches between solving for the systems
 * <i>P<sup>-1</sup>A</i> and <i>AP<sup>-1</sup></i>, respectively.
 *
 * A second consequence is the type of residual which is used to measure
 * convergence. With left preconditioning, this is the <b>preconditioned</b>
 * residual, while with right preconditioning, it is the residual of the
 * unpreconditioned system.
 *
 * Optionally, this behavior can be overridden by using the flag
 * AdditionalData::use_default_residual. A <tt>true</tt> value refers to the
 * behavior described in the previous paragraph, while <tt>false</tt> reverts
 * it. Be aware though that additional residuals have to be computed in this
 * case, impeding the overall performance of the solver.
 *
 *
 * <h3>The size of the Arnoldi basis</h3>
 *
 * The maximal basis size is controlled by AdditionalData::max_n_tmp_vectors,
 * and it is this number minus 2. If the number of iteration steps exceeds
 * this number, all basis vectors are discarded and the iteration starts anew
 * from the approximation obtained so far.
 *
 * Note that the minimizing property of GMRes only pertains to the Krylov
 * space spanned by the Arnoldi basis. Therefore, restarted GMRes is
 * <b>not</b> minimizing anymore. The choice of the basis length is a trade-
 * off between memory consumption and convergence speed, since a longer basis
 * means minimization over a larger space.
 *
 * For the requirements on matrices and vectors in order to work with this
 * class, see the documentation of the Solver base class.
 *
 *
 * <h3>Observing the progress of linear solver iterations</h3>
 *
 * The solve() function of this class uses the mechanism described in the
 * Solver base class to determine convergence. This mechanism can also be used
 * to observe the progress of the iteration.
 *
 *
 * <h3>Eigenvalue and condition number estimates</h3>
 *
 * This class can estimate eigenvalues and condition number during the
 * solution process. This is done by creating the Hessenberg matrix during the
 * inner iterations. The eigenvalues are estimated as the eigenvalues of the
 * Hessenberg matrix and the condition number is estimated as the ratio of the
 * largest and smallest singular value of the Hessenberg matrix. The estimates
 * can be obtained by connecting a function as a slot using @p
 * connect_condition_number_slot and @p connect_eigenvalues_slot. These slots
 * will then be called from the solver with the estimates as argument.
 *
 *
 * @author Wolfgang Bangerth, Guido Kanschat, Ralf Hartmann.
 */
template <class VectorType = Vector<double> >
class SolverGMRES : public Solver<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the number of temporary vectors to 30,
     * i.e. do a restart every 28 iterations. Also set preconditioning from
     * left, the residual of the stopping criterion to the default residual,
     * and re-orthogonalization only if necessary.
     */
    explicit
    AdditionalData (const unsigned int max_n_tmp_vectors = 30,
                    const bool right_preconditioning = false,
                    const bool use_default_residual = true,
                    const bool force_re_orthogonalization = false);

    /**
     * Constructor.
     * @deprecated To obtain the estimated eigenvalues instead use:
     * connect_eigenvalues_slot
     */
    AdditionalData (const unsigned int max_n_tmp_vectors,
                    const bool right_preconditioning,
                    const bool use_default_residual,
                    const bool force_re_orthogonalization,
                    const bool compute_eigenvalues) DEAL_II_DEPRECATED;

    /**
     * Maximum number of temporary vectors. This parameter controls the size
     * of the Arnoldi basis, which for historical reasons is
     * #max_n_tmp_vectors-2.
     */
    unsigned int    max_n_tmp_vectors;

    /**
     * Flag for right preconditioning.
     *
     * @note Change between left and right preconditioning will also change
     * the way residuals are evaluated. See the corresponding section in the
     * SolverGMRES.
     */
    bool right_preconditioning;

    /**
     * Flag for the default residual that is used to measure convergence.
     */
    bool use_default_residual;

    /**
     * Flag to force re-orthogonalization of orthonormal basis in every step.
     * If set to false, the solver automatically checks for loss of
     * orthogonality every 5 iterations and enables re-orthogonalization only
     * if necessary.
     */
    bool force_re_orthogonalization;

    /**
     * Compute all eigenvalues of the Hessenberg matrix generated while
     * solving, i.e., the projected system matrix. This gives an approximation
     * of the eigenvalues of the (preconditioned) system matrix. Since the
     * Hessenberg matrix is thrown away at restart, the eigenvalues are
     * printed for every 30 iterations.
     *
     * @note Requires LAPACK support.
     */
    bool compute_eigenvalues;
  };

  /**
   * Constructor.
   */
  SolverGMRES (SolverControl            &cn,
               VectorMemory<VectorType> &mem,
               const AdditionalData     &data=AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverGMRES (SolverControl        &cn,
               const AdditionalData &data=AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template<typename MatrixType, typename PreconditionerType>
  void
  solve (const MatrixType         &A,
         VectorType               &x,
         const VectorType         &b,
         const PreconditionerType &precondition);

  /**
   * Connect a slot to retrieve the estimated condition number. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_condition_number_slot(const std_cxx11::function<void (double)> &slot,
                                const bool every_iteration=false);

  /**
   * Connect a slot to retrieve the estimated eigenvalues. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_eigenvalues_slot(
    const std_cxx11::function<void (const std::vector<std::complex<double> > &)> &slot,
    const bool every_iteration=false);

  /**
   * Connect a slot to retrieve the Hessenberg matrix obtained by the
   * projection of the initial matrix on the Krylov basis. Called on each
   * outer iteration if every_iteration=true, otherwise called once when
   * iterations are ended (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_hessenberg_slot(
    const std_cxx11::function<void (const FullMatrix<double> &)> &slot,
    const bool every_iteration=true);

  /**
   * Connect a slot to retrieve the basis vectors of the Krylov space
   * generated by the Arnoldi algorithm. Called at once when iterations
   * are completed (i.e., either because convergence has been achieved,
   * or because divergence has been detected).
   */
  boost::signals2::connection
  connect_krylov_space_slot(
    const std_cxx11::function<void (const internal::SolverGMRES::TmpVectors<VectorType> &)> &slot);


  DeclException1 (ExcTooFewTmpVectors,
                  int,
                  << "The number of temporary vectors you gave ("
                  << arg1 << ") is too small. It should be at least 10 for "
                  << "any results, and much more for reasonable ones.");

protected:

  /**
   * Includes the maximum number of tmp vectors.
   */
  AdditionalData additional_data;

  /**
   * Signal used to retrieve the estimated condition number. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void (double)> condition_number_signal;

  /**
   * Signal used to retrieve the estimated condition numbers. Called on each
   * outer iteration.
   */
  boost::signals2::signal<void (double)> all_condition_numbers_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called once when all
   * iterations are ended.
   */
  boost::signals2::signal<void (const std::vector<std::complex<double> > &)> eigenvalues_signal;

  /**
   * Signal used to retrieve the estimated eigenvalues. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void (const std::vector<std::complex<double> > &)> all_eigenvalues_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called once when
   * all iterations are ended.
   */
  boost::signals2::signal<void (const FullMatrix<double> &)> hessenberg_signal;

  /**
   * Signal used to retrieve the Hessenberg matrix. Called on each outer
   * iteration.
   */
  boost::signals2::signal<void (const FullMatrix<double> &)> all_hessenberg_signal;

  /**
   * Signal used to retrieve the Krylov space basis vectors. Called once
   * when all iterations are ended.
   */
  boost::signals2::signal<void (const internal::SolverGMRES::TmpVectors<VectorType> &)> krylov_space_signal;

  /**
   * Implementation of the computation of the norm of the residual.
   */
  virtual double criterion();

  /**
   * Transformation of an upper Hessenberg matrix into tridiagonal structure
   * by givens rotation of the last column
   */
  void givens_rotation (Vector<double> &h,  Vector<double> &b,
                        Vector<double> &ci, Vector<double> &si,
                        int col) const;

  /**
   * Orthogonalize the vector @p vv against the @p dim (orthogonal) vectors
   * given by the first argument using the modified Gram-Schmidt algorithm.
   * The factors used for orthogonalization are stored in @p h. The boolean @p
   * re_orthogonalize specifies whether the modified Gram-Schmidt algorithm
   * should be applied twice. The algorithm checks loss of orthogonality in
   * the procedure every fifth step and sets the flag to true in that case.
   * All subsequent iterations use re-orthogonalization.
   */
  static double
  modified_gram_schmidt (const internal::SolverGMRES::TmpVectors<VectorType> &orthogonal_vectors,
                         const unsigned int                                  dim,
                         const unsigned int                                  accumulated_iterations,
                         VectorType                                          &vv,
                         Vector<double>                                      &h,
                         bool                                                &re_orthogonalize);

  /**
   * Estimates the eigenvalues from the Hessenberg matrix, H_orig, generated
   * during the inner iterations. Uses these estimate to compute the condition
   * number. Calls the signals eigenvalues_signal and cond_signal with these
   * estimates as arguments. Outputs the eigenvalues to deallog if
   * log_eigenvalues is true.
   */
  static void
  compute_eigs_and_cond(
    const FullMatrix<double> &H_orig ,
    const unsigned int dim,
    const boost::signals2::signal<void (const std::vector<std::complex<double> > &)> &eigenvalues_signal,
    const boost::signals2::signal<void (const FullMatrix<double> &)> &hessenberg_signal,
    const boost::signals2::signal<void(double)> &cond_signal,
    const bool log_eigenvalues);

  /**
   * Projected system matrix
   */
  FullMatrix<double> H;

  /**
   * Auxiliary matrix for inverting @p H
   */
  FullMatrix<double> H1;


private:
  /**
   * No copy constructor.
   */
  SolverGMRES (const SolverGMRES<VectorType> &);
};

/**
 * Implementation of the Generalized minimal residual method with flexible
 * preconditioning method.
 *
 * This version of the GMRES method allows for the use of a different
 * preconditioner in each iteration step. Therefore, it is also more robust
 * with respect to inaccurate evaluation of the preconditioner. An important
 * application is also the use of a Krylov space method inside the
 * preconditioner. As opposed to SolverGMRES which allows one to choose
 * between left and right preconditioning, this solver always applies the
 * preconditioner from the right.
 *
 * FGMRES needs two vectors in each iteration steps yielding a total of
 * <tt>2*SolverFGMRES::AdditionalData::max_basis_size+1</tt> auxiliary
 * vectors.
 *
 * Caveat: Documentation of this class is not up to date. There are also a few
 * parameters of GMRES we would like to introduce here.
 *
 * @author Guido Kanschat, 2003
 */
template <class VectorType = Vector<double> >
class SolverFGMRES : public Solver<VectorType>
{
public:
  /**
   * Standardized data struct to pipe additional data to the solver.
   */
  struct AdditionalData
  {
    /**
     * Constructor. By default, set the maximum basis size to 30.
     */
    explicit
    AdditionalData(const unsigned int max_basis_size = 30,
                   const bool /*use_default_residual*/ = true)
      :
      max_basis_size(max_basis_size)
    {}

    /**
     * Maximum number of tmp vectors.
     */
    unsigned int    max_basis_size;
  };

  /**
   * Constructor.
   */
  SolverFGMRES (SolverControl            &cn,
                VectorMemory<VectorType> &mem,
                const AdditionalData     &data=AdditionalData());

  /**
   * Constructor. Use an object of type GrowingVectorMemory as a default to
   * allocate memory.
   */
  SolverFGMRES (SolverControl        &cn,
                const AdditionalData &data=AdditionalData());

  /**
   * Solve the linear system $Ax=b$ for x.
   */
  template<typename MatrixType, typename PreconditionerType>
  void
  solve (const MatrixType         &A,
         VectorType               &x,
         const VectorType         &b,
         const PreconditionerType &precondition);

private:

  /**
   * Additional flags.
   */
  AdditionalData additional_data;

  /**
   * Projected system matrix
   */
  FullMatrix<double> H;

  /**
   * Auxiliary matrix for inverting @p H
   */
  FullMatrix<double> H1;
};

/*@}*/
/* --------------------- Inline and template functions ------------------- */


#ifndef DOXYGEN
namespace internal
{
  namespace SolverGMRES
  {
    template <class VectorType>
    inline
    TmpVectors<VectorType>::
    TmpVectors (const unsigned int       max_size,
                VectorMemory<VectorType> &vmem)
      :
      mem(vmem),
      data (max_size, 0),
      offset(0)
    {}


    template <class VectorType>
    inline
    TmpVectors<VectorType>::~TmpVectors ()
    {
      for (typename std::vector<VectorType *>::iterator v = data.begin();
           v != data.end(); ++v)
        if (*v != 0)
          mem.free(*v);
    }


    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator[] (const unsigned int i) const
    {
      Assert (i+offset<data.size(),
              ExcIndexRange(i, -offset, data.size()-offset));

      Assert (data[i-offset] != 0, ExcNotInitialized());
      return *data[i-offset];
    }


    template <class VectorType>
    inline VectorType &
    TmpVectors<VectorType>::operator() (const unsigned int i,
                                        const VectorType       &temp)
    {
      Assert (i+offset<data.size(),
              ExcIndexRange(i,-offset, data.size()-offset));
      if (data[i-offset] == 0)
        {
          data[i-offset] = mem.alloc();
          data[i-offset]->reinit(temp);
        }
      return *data[i-offset];
    }


    template <class VectorType>
    unsigned int
    TmpVectors<VectorType>::size() const
    {
      return (data.size() > 0 ? data.size()-1 : 0);
    }


    // A comparator for better printing eigenvalues
    inline
    bool complex_less_pred(const std::complex<double> &x,
                           const std::complex<double> &y)
    {
      return x.real() < y.real() || (x.real() == y.real() && x.imag() < y.imag());
    }
  }
}



template <class VectorType>
inline
SolverGMRES<VectorType>::AdditionalData::
AdditionalData (const unsigned int max_n_tmp_vectors,
                const bool         right_preconditioning,
                const bool         use_default_residual,
                const bool         force_re_orthogonalization)
  :
  max_n_tmp_vectors(max_n_tmp_vectors),
  right_preconditioning(right_preconditioning),
  use_default_residual(use_default_residual),
  force_re_orthogonalization(force_re_orthogonalization),
  compute_eigenvalues(false)
{}



template <class VectorType>
inline
SolverGMRES<VectorType>::AdditionalData::
AdditionalData (const unsigned int max_n_tmp_vectors,
                const bool         right_preconditioning,
                const bool         use_default_residual,
                const bool         force_re_orthogonalization,
                const bool         compute_eigenvalues)
  :
  max_n_tmp_vectors(max_n_tmp_vectors),
  right_preconditioning(right_preconditioning),
  use_default_residual(use_default_residual),
  force_re_orthogonalization(force_re_orthogonalization),
  compute_eigenvalues(compute_eigenvalues)
{}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES (SolverControl            &cn,
                                      VectorMemory<VectorType> &mem,
                                      const AdditionalData     &data)
  :
  Solver<VectorType> (cn,mem),
  additional_data(data)
{}



template <class VectorType>
SolverGMRES<VectorType>::SolverGMRES (SolverControl        &cn,
                                      const AdditionalData &data) :
  Solver<VectorType> (cn),
  additional_data(data)
{}



template <class VectorType>
inline
void
SolverGMRES<VectorType>::givens_rotation (Vector<double> &h,
                                          Vector<double> &b,
                                          Vector<double> &ci,
                                          Vector<double> &si,
                                          int            col) const
{
  for (int i=0 ; i<col ; i++)
    {
      const double s = si(i);
      const double c = ci(i);
      const double dummy = h(i);
      h(i)   =  c*dummy + s*h(i+1);
      h(i+1) = -s*dummy + c*h(i+1);
    };

  const double r = 1./std::sqrt(h(col)*h(col) + h(col+1)*h(col+1));
  si(col) = h(col+1) *r;
  ci(col) = h(col)   *r;
  h(col)  =  ci(col)*h(col) + si(col)*h(col+1);
  b(col+1)= -si(col)*b(col);
  b(col) *=  ci(col);
}



template <class VectorType>
inline
double
SolverGMRES<VectorType>::modified_gram_schmidt
(const internal::SolverGMRES::TmpVectors<VectorType> &orthogonal_vectors,
 const unsigned int                                  dim,
 const unsigned int                                  accumulated_iterations,
 VectorType                                          &vv,
 Vector<double>                                      &h,
 bool                                                &re_orthogonalize)
{
  Assert(dim > 0, ExcInternalError());
  const unsigned int inner_iteration = dim - 1;

  // need initial norm for detection of re-orthogonalization, see below
  double norm_vv_start = 0;
  if (re_orthogonalize == false && inner_iteration % 5 == 4)
    norm_vv_start = vv.l2_norm();

  // Orthogonalization
  h(0) = vv * orthogonal_vectors[0];
  for (unsigned int i=1 ; i<dim ; ++i)
    h(i) = vv.add_and_dot(-h(i-1), orthogonal_vectors[i-1], orthogonal_vectors[i]);
  double norm_vv = std::sqrt(vv.add_and_dot(-h(dim-1), orthogonal_vectors[dim-1], vv));

  // Re-orthogonalization if loss of orthogonality detected. For the test, use
  // a strategy discussed in C. T. Kelley, Iterative Methods for Linear and
  // Nonlinear Equations, SIAM, Philadelphia, 1995: Compare the norm of vv
  // after orthogonalization with its norm when starting the
  // orthogonalization. If vv became very small (here: less than the square
  // root of the machine precision times 10), it is almost in the span of the
  // previous vectors, which indicates loss of precision.
  if (re_orthogonalize == false && inner_iteration % 5 == 4)
    {
      if (norm_vv > 10. * norm_vv_start *
          std::sqrt(std::numeric_limits<typename VectorType::value_type>::epsilon()))
        return norm_vv;

      else
        {
          re_orthogonalize = true;
          deallog << "Re-orthogonalization enabled at step "
                  << accumulated_iterations << std::endl;
        }
    }

  if (re_orthogonalize == true)
    {
      double htmp = vv * orthogonal_vectors[0];
      h(0) += htmp;
      for (unsigned int i=1 ; i<dim ; ++i)
        {
          htmp = vv.add_and_dot(-htmp, orthogonal_vectors[i-1], orthogonal_vectors[i]);
          h(i) += htmp;
        }
      norm_vv = std::sqrt(vv.add_and_dot(-htmp, orthogonal_vectors[dim-1], vv));
    }

  return norm_vv;
}



template<class VectorType>
inline void
SolverGMRES<VectorType>::compute_eigs_and_cond
(const FullMatrix<double>                     &H_orig,
 const unsigned int                           dim,
 const boost::signals2::signal<void (const std::vector<std::complex<double> > &)> &eigenvalues_signal,
 const boost::signals2::signal<void (const FullMatrix<double> &)> &hessenberg_signal,
 const boost::signals2::signal<void (double)> &cond_signal,
 const bool                                   log_eigenvalues)
{
  //Avoid copying the Hessenberg matrix if it isn't needed.
  if (!eigenvalues_signal.empty() || !hessenberg_signal.empty()
      || !cond_signal.empty() || log_eigenvalues )
    {
      LAPACKFullMatrix<double> mat(dim,dim);
      for (unsigned int i=0; i<dim; ++i)
        for (unsigned int j=0; j<dim; ++j)
          mat(i,j) = H_orig(i,j);
      hessenberg_signal(H_orig);
      //Avoid computing eigenvalues if they are not needed.
      if (!eigenvalues_signal.empty() || log_eigenvalues )
        {
          //Copy mat so that we can compute svd below. Necessary since
          //compute_eigenvalues will leave mat in state LAPACKSupport::unusable.
          LAPACKFullMatrix<double> mat_eig(mat);
          mat_eig.compute_eigenvalues();
          std::vector<std::complex<double> > eigenvalues(dim);
          for (unsigned int i=0; i<mat_eig.n(); ++i)
            eigenvalues[i] = mat_eig.eigenvalue(i);
          //Sort eigenvalues for nicer output.
          std::sort(eigenvalues.begin(), eigenvalues.end(),
                    internal::SolverGMRES::complex_less_pred);
          eigenvalues_signal(eigenvalues);
          if (log_eigenvalues)
            {
              deallog << "Eigenvalue estimate: ";
              for (unsigned int i=0; i<mat_eig.n(); ++i)
                deallog << ' ' << eigenvalues[i];
              deallog << std::endl;
            }
        }
      //Calculate condition number, avoid calculating the svd if a slot
      //isn't connected. Need at least a 2-by-2 matrix to do the estimate.
      if (!cond_signal.empty() && (mat.n()>1))
        {
          mat.compute_svd();
          double condition_number=mat.singular_value(0)/mat.singular_value(mat.n()-1);
          cond_signal(condition_number);
        }
    }
}



template<class VectorType>
template<typename MatrixType, typename PreconditionerType>
void
SolverGMRES<VectorType>::solve (const MatrixType         &A,
                                VectorType               &x,
                                const VectorType         &b,
                                const PreconditionerType &precondition)
{
  // this code was written a very long time ago by people not associated with
  // deal.II. we don't make any guarantees to its optimality or that it even
  // works as expected...

//TODO:[?] Check, why there are two different start residuals.
//TODO:[GK] Make sure the parameter in the constructor means maximum basis size

  deallog.push("GMRES");
  const unsigned int n_tmp_vectors = additional_data.max_n_tmp_vectors;

  // Generate an object where basis vectors are stored.
  internal::SolverGMRES::TmpVectors<VectorType> tmp_vectors (n_tmp_vectors, this->memory);

  // number of the present iteration; this
  // number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  const bool do_eigenvalues=
    !condition_number_signal.empty()
    |!all_condition_numbers_signal.empty()
    |!eigenvalues_signal.empty()
    |!all_eigenvalues_signal.empty()
    |!hessenberg_signal.empty()
    |!all_hessenberg_signal.empty()
    |additional_data.compute_eigenvalues;
  // for eigenvalue computation, need to collect the Hessenberg matrix (before
  // applying Givens rotations)
  FullMatrix<double> H_orig;
  if (do_eigenvalues)
    H_orig.reinit(n_tmp_vectors, n_tmp_vectors-1);

  // matrix used for the orthogonalization process later
  H.reinit(n_tmp_vectors, n_tmp_vectors-1);

  // some additional vectors, also used in the orthogonalization
  dealii::Vector<double>
  gamma(n_tmp_vectors),
        ci   (n_tmp_vectors-1),
        si   (n_tmp_vectors-1),
        h    (n_tmp_vectors-1);


  unsigned int dim = 0;

  SolverControl::State iteration_state = SolverControl::iterate;
  double last_res = -std::numeric_limits<double>::max();

  // switch to determine whether we want a left or a right preconditioner. at
  // present, left is default, but both ways are implemented
  const bool left_precondition = !additional_data.right_preconditioning;

  // Per default the left preconditioned GMRes uses the preconditioned
  // residual and the right preconditioned GMRes uses the unpreconditioned
  // residual as stopping criterion.
  const bool use_default_residual = additional_data.use_default_residual;

  // define two aliases
  VectorType &v = tmp_vectors(0, x);
  VectorType &p = tmp_vectors(n_tmp_vectors-1, x);

  // Following vectors are needed
  // when not the default residuals
  // are used as stopping criterion
  VectorType *r=0;
  VectorType *x_=0;
  dealii::Vector<double> *gamma_=0;
  if (!use_default_residual)
    {
      r=this->memory.alloc();
      x_=this->memory.alloc();
      r->reinit(x);
      x_->reinit(x);

      gamma_ = new dealii::Vector<double> (gamma.size());
    }

  bool re_orthogonalize = additional_data.force_re_orthogonalization;

  ///////////////////////////////////////////////////////////////////////////
  // outer iteration: loop until we either reach convergence or the maximum
  // number of iterations is exceeded. each cycle of this loop amounts to one
  // restart
  do
    {
      // reset this vector to the right size
      h.reinit (n_tmp_vectors-1);

      if (left_precondition)
        {
          A.vmult(p,x);
          p.sadd(-1.,1.,b);
          precondition.vmult(v,p);
        }
      else
        {
          A.vmult(v,x);
          v.sadd(-1.,1.,b);
        };

      double rho = v.l2_norm();

      // check the residual here as well since it may be that we got the exact
      // (or an almost exact) solution vector at the outset. if we wouldn't
      // check here, the next scaling operation would produce garbage
      if (use_default_residual)
        {
          last_res = rho;
          iteration_state = this->iteration_status (accumulated_iterations, rho, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }
      else
        {
          deallog << "default_res=" << rho << std::endl;

          if (left_precondition)
            {
              A.vmult(*r,x);
              r->sadd(-1.,1.,b);
            }
          else
            precondition.vmult(*r,v);

          double res = r->l2_norm();
          last_res = res;
          iteration_state = this->iteration_status (accumulated_iterations, res, x);

          if (iteration_state != SolverControl::iterate)
            break;
        }

      gamma(0) = rho;

      v *= 1./rho;

      // inner iteration doing at most as many steps as there are temporary
      // vectors. the number of steps actually been done is propagated outside
      // through the @p dim variable
      for (unsigned int inner_iteration=0;
           ((inner_iteration < n_tmp_vectors-2)
            &&
            (iteration_state==SolverControl::iterate));
           ++inner_iteration)
        {
          ++accumulated_iterations;
          // yet another alias
          VectorType &vv = tmp_vectors(inner_iteration+1, x);

          if (left_precondition)
            {
              A.vmult(p, tmp_vectors[inner_iteration]);
              precondition.vmult(vv,p);
            }
          else
            {
              precondition.vmult(p, tmp_vectors[inner_iteration]);
              A.vmult(vv,p);
            }

          dim = inner_iteration+1;

          const double s = modified_gram_schmidt(tmp_vectors, dim,
                                                 accumulated_iterations,
                                                 vv, h, re_orthogonalize);
          h(inner_iteration+1) = s;

          //s=0 is a lucky breakdown, the solver will reach convergence,
          //but we must not divide by zero here.
          if (s != 0)
            vv *= 1./s;

          // for eigenvalues, get the resulting coefficients from the
          // orthogonalization process
          if (do_eigenvalues)
            for (unsigned int i=0; i<dim+1; ++i)
              H_orig(i,inner_iteration) = h(i);

          //  Transformation into tridiagonal structure
          givens_rotation(h,gamma,ci,si,inner_iteration);

          //  append vector on matrix
          for (unsigned int i=0; i<dim; ++i)
            H(i,inner_iteration) = h(i);

          //  default residual
          rho = std::fabs(gamma(dim));

          if (use_default_residual)
            {
              last_res = rho;
              iteration_state = this->iteration_status (accumulated_iterations, rho, x);
            }
          else
            {
              deallog << "default_res=" << rho << std::endl;

              dealii::Vector<double> h_(dim);
              *x_=x;
              *gamma_=gamma;
              H1.reinit(dim+1,dim);

              for (unsigned int i=0; i<dim+1; ++i)
                for (unsigned int j=0; j<dim; ++j)
                  H1(i,j) = H(i,j);

              H1.backward(h_,*gamma_);

              if (left_precondition)
                for (unsigned int i=0 ; i<dim; ++i)
                  x_->add(h_(i), tmp_vectors[i]);
              else
                {
                  p = 0.;
                  for (unsigned int i=0; i<dim; ++i)
                    p.add(h_(i), tmp_vectors[i]);
                  precondition.vmult(*r,p);
                  x_->add(1.,*r);
                };
              A.vmult(*r,*x_);
              r->sadd(-1.,1.,b);
              // Now *r contains the unpreconditioned residual!!
              if (left_precondition)
                {
                  const double res=r->l2_norm();
                  last_res = res;

                  iteration_state = this->iteration_status (accumulated_iterations, res, x);
                }
              else
                {
                  precondition.vmult(*x_, *r);
                  const double preconditioned_res=x_->l2_norm();
                  last_res = preconditioned_res;

                  iteration_state = this->iteration_status (accumulated_iterations,
                                                            preconditioned_res, x);
                }
            }
        };
      // end of inner iteration. now calculate the solution from the temporary
      // vectors
      h.reinit(dim);
      H1.reinit(dim+1,dim);

      for (unsigned int i=0; i<dim+1; ++i)
        for (unsigned int j=0; j<dim; ++j)
          H1(i,j) = H(i,j);

      compute_eigs_and_cond(H_orig,dim,all_eigenvalues_signal,
                            all_hessenberg_signal,
                            all_condition_numbers_signal,
                            additional_data.compute_eigenvalues);

      H1.backward(h,gamma);

      if (left_precondition)
        for (unsigned int i=0 ; i<dim; ++i)
          x.add(h(i), tmp_vectors[i]);
      else
        {
          p = 0.;
          for (unsigned int i=0; i<dim; ++i)
            p.add(h(i), tmp_vectors[i]);
          precondition.vmult(v,p);
          x.add(1.,v);
        };
      // end of outer iteration. restart if no convergence and the number of
      // iterations is not exceeded
    }
  while (iteration_state == SolverControl::iterate);

  compute_eigs_and_cond(H_orig,dim,eigenvalues_signal,hessenberg_signal,
                        condition_number_signal,
                        false);

  if (!krylov_space_signal.empty())
    krylov_space_signal(tmp_vectors);

  if (!use_default_residual)
    {
      this->memory.free(r);
      this->memory.free(x_);

      delete gamma_;
    }

  deallog.pop();

  // in case of failure: throw exception
  AssertThrow(iteration_state == SolverControl::success,
              SolverControl::NoConvergence (accumulated_iterations,
                                            last_res));
}



template<class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_condition_number_slot
(const std_cxx11::function<void(double)> &slot,
 const bool every_iteration)
{
  if (every_iteration)
    {
      return all_condition_numbers_signal.connect(slot);
    }
  else
    {
      return condition_number_signal.connect(slot);
    }
}



template<class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_eigenvalues_slot
(const std_cxx11::function<void (const std::vector<std::complex<double> > &)> &slot,
 const bool every_iteration)
{
  if (every_iteration)
    {
      return all_eigenvalues_signal.connect(slot);
    }
  else
    {
      return eigenvalues_signal.connect(slot);
    }
}



template<class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_hessenberg_slot
(const std_cxx11::function<void (const FullMatrix<double> &)> &slot,
 const bool every_iteration)
{
  if (every_iteration)
    {
      return all_hessenberg_signal.connect(slot);
    }
  else
    {
      return hessenberg_signal.connect(slot);
    }
}



template<class VectorType>
boost::signals2::connection
SolverGMRES<VectorType>::connect_krylov_space_slot
(const std_cxx11::function<void (const internal::SolverGMRES::TmpVectors<VectorType> &)> &slot)
{
  return krylov_space_signal.connect(slot);
}



template<class VectorType>
double
SolverGMRES<VectorType>::criterion ()
{
  // dummy implementation. this function is not needed for the present
  // implementation of gmres
  Assert (false, ExcInternalError());
  return 0;
}


//----------------------------------------------------------------------//

template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES (SolverControl            &cn,
                                        VectorMemory<VectorType> &mem,
                                        const AdditionalData     &data)
  :
  Solver<VectorType> (cn, mem),
  additional_data(data)
{}



template <class VectorType>
SolverFGMRES<VectorType>::SolverFGMRES (SolverControl        &cn,
                                        const AdditionalData &data)
  :
  Solver<VectorType> (cn),
  additional_data(data)
{}



template<class VectorType>
template<typename MatrixType, typename PreconditionerType>
void
SolverFGMRES<VectorType>::solve (const MatrixType         &A,
                                 VectorType               &x,
                                 const VectorType         &b,
                                 const PreconditionerType &precondition)
{
  deallog.push("FGMRES");

  SolverControl::State iteration_state = SolverControl::iterate;

  const unsigned int basis_size = additional_data.max_basis_size;

  // Generate an object where basis vectors are stored.
  typename internal::SolverGMRES::TmpVectors<VectorType> v (basis_size, this->memory);
  typename internal::SolverGMRES::TmpVectors<VectorType> z (basis_size, this->memory);

  // number of the present iteration; this number is not reset to zero upon a
  // restart
  unsigned int accumulated_iterations = 0;

  // matrix used for the orthogonalization process later
  H.reinit(basis_size+1, basis_size);

  // Vectors for projected system
  Vector<double> projected_rhs;
  Vector<double> y;

  // Iteration starts here
  double res = -std::numeric_limits<double>::max();

  VectorType *aux = this->memory.alloc();
  aux->reinit(x);
  do
    {
      A.vmult(*aux, x);
      aux->sadd(-1., 1., b);

      double beta = aux->l2_norm();
      res = beta;
      iteration_state = this->iteration_status(accumulated_iterations, res, x);
      if (iteration_state == SolverControl::success)
        break;

      H.reinit(basis_size+1, basis_size);
      double a = beta;

      for (unsigned int j=0; j<basis_size; ++j)
        {
          if (a != 0) // treat lucky breakdown
            v(j,x).equ(1./a, *aux);
          else
            v(j,x) = 0.;


          precondition.vmult(z(j,x), v[j]);
          A.vmult(*aux, z[j]);

          // Gram-Schmidt
          H(0,j) = *aux * v[0];
          for (unsigned int i=1; i<=j; ++i)
            H(i,j) = aux->add_and_dot(-H(i-1,j), v[i-1], v[i]);
          H(j+1,j) = a = std::sqrt(aux->add_and_dot(-H(j,j), v[j], *aux));

          // Compute projected solution

          if (j>0)
            {
              H1.reinit(j+1,j);
              projected_rhs.reinit(j+1);
              y.reinit(j);
              projected_rhs(0) = beta;
              H1.fill(H);

              // check convergence. note that the vector 'x' we pass to the
              // criterion is not the final solution we compute if we
              // decide to jump out of the iteration (we update 'x' again
              // right after the current loop)
              Householder<double> house(H1);
              res = house.least_squares(y, projected_rhs);
              iteration_state = this->iteration_status(++accumulated_iterations, res, x);
              if (iteration_state != SolverControl::iterate)
                break;
            }
        }

      // Update solution vector
      for (unsigned int j=0; j<y.size(); ++j)
        x.add(y(j), z[j]);
    }
  while (iteration_state == SolverControl::iterate);

  this->memory.free(aux);

  deallog.pop();
  // in case of failure: throw exception
  if (iteration_state != SolverControl::success)
    AssertThrow(false, SolverControl::NoConvergence (accumulated_iterations,
                                                     res));
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
