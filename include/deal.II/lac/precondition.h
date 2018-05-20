// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2018 by the deal.II authors
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

#ifndef dealii_precondition_h
#define dealii_precondition_h

// This file contains simple preconditioners.

#include <deal.II/base/config.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/lac/diagonal_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

// forward declarations

template <typename number>
class Vector;
template <typename number>
class SparseMatrix;
namespace LinearAlgebra
{
  namespace distributed
  {
    template <typename number>
    class Vector;
  }
} // namespace LinearAlgebra

/*! @addtogroup Preconditioners
 *@{
 */

/**
 * No preconditioning.  This class helps you, if you want to use a linear
 * solver without preconditioning. All solvers in LAC require a
 * preconditioner. Therefore, you must use the identity provided here to avoid
 * preconditioning. It can be used in the following way:
 *
 * @code
 * SolverControl           solver_control (1000, 1e-12);
 * SolverCG<>              cg (solver_control);
 * cg.solve (system_matrix, solution, system_rhs,
 *          PreconditionIdentity());
 * @endcode
 *
 * See the step-3 tutorial program for an example and additional explanations.
 *
 * Alternatively, the IdentityMatrix class can be used to precondition in this
 * way.
 *
 * @author Guido Kanschat, 1999; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
class PreconditionIdentity : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData() = default;
  };

  /**
   * Constructor, sets the domain and range sizes to their defaults.
   */
  PreconditionIdentity();

  /**
   * The matrix argument is ignored and here just for compatibility with more
   * complex preconditioners.
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType&     matrix,
             const AdditionalData& additional_data = AdditionalData());

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;

  /**
   * Apply preconditioner, adding to the previous value.
   */
  template <class VectorType>
  void
  vmult_add(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType&, const VectorType&) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void
  clear()
  {}

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  n() const;

private:
  /**
   * The dimension of the range space.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space.
   */
  size_type n_columns;
};

/**
 * Preconditioning with Richardson's method. This preconditioner just scales
 * the vector with a constant relaxation factor provided by the AdditionalData
 * object.
 *
 * In Krylov-space methods, this preconditioner should not have any effect.
 * Using SolverRichardson, the two relaxation parameters will be just
 * multiplied. Still, this class is useful in multigrid smoother objects
 * (MGSmootherRelaxation).
 *
 * @author Guido Kanschat, 2005; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
class PreconditionRichardson : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Parameters for Richardson preconditioner.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. Block size must be given since there is no reasonable
     * default parameter.
     */
    AdditionalData(const double relaxation = 1.);

    /**
     * Relaxation parameter.
     */
    double relaxation;
  };

  /**
   * Constructor, sets the relaxation parameter, domain and range sizes to
   * their default.
   */
  PreconditionRichardson();

  /**
   * Change the relaxation parameter.
   */
  void
  initialize(const AdditionalData& parameters);

  /**
   * Change the relaxation parameter in a way consistent with other
   * preconditioners. The matrix argument is ignored and here just for
   * compatibility with more complex preconditioners.
   */
  template <typename MatrixType>
  void
  initialize(const MatrixType& matrix, const AdditionalData& parameters);

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;
  /**
   * Apply preconditioner, adding to the previous value.
   */
  template <class VectorType>
  void
  vmult_add(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template <class VectorType>
  void
  Tvmult_add(VectorType&, const VectorType&) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void
  clear()
  {}

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   *
   * @note This function should only be called if the preconditioner has been
   * initialized.
   */
  size_type
  n() const;

private:
  /**
   * The relaxation parameter multiplied with the vectors.
   */
  double relaxation;

  /**
   * The dimension of the range space.
   */
  size_type n_rows;

  /**
   * The dimension of the domain space.
   */
  size_type n_columns;
};

/**
 * Preconditioner using a matrix-builtin function.  This class forms a
 * preconditioner suitable for the LAC solver classes. Since many
 * preconditioning methods are based on matrix entries, these have to be
 * implemented as member functions of the underlying matrix implementation.
 * This class now is intended to allow easy access to these member functions
 * from LAC solver classes.
 *
 * It seems that all builtin preconditioners have a relaxation parameter, so
 * please use PreconditionRelaxation for these.
 *
 * You will usually not want to create a named object of this type, although
 * possible. The most common use is like this:
 * @code
 *    SolverGMRES<SparseMatrix<double>,
 *                Vector<double> >      gmres(control,memory,500);
 *
 *    gmres.solve (matrix, solution, right_hand_side,
 *                 PreconditionUseMatrix<SparseMatrix<double>,Vector<double> >
 *                 (matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 * This creates an unnamed object to be passed as the fourth parameter to the
 * solver function of the SolverGMRES class. It assumes that the SparseMatrix
 * class has a function <tt>precondition_Jacobi</tt> taking two vectors
 * (source and destination) as parameters (Actually, there is no function like
 * that, the existing function takes a third parameter, denoting the
 * relaxation parameter; this example is therefore only meant to illustrate
 * the general idea).
 *
 * Note that due to the default template parameters, the above example could
 * be written shorter as follows:
 * @code
 *    ...
 *    gmres.solve (matrix, solution, right_hand_side,
 *                 PreconditionUseMatrix<>
 *                   (matrix,&SparseMatrix<double>::template precondition_Jacobi<double>));
 * @endcode
 *
 * @author Guido Kanschat, Wolfgang Bangerth, 1999
 */
template <typename MatrixType = SparseMatrix<double>,
          class VectorType    = Vector<double>>
class PreconditionUseMatrix : public Subscriptor
{
public:
  /**
   * Type of the preconditioning function of the matrix.
   */
  typedef void (MatrixType::*function_ptr)(VectorType&,
                                           const VectorType&) const;

  /**
   * Constructor.  This constructor stores a reference to the matrix object
   * for later use and selects a preconditioning method, which must be a
   * member function of that matrix.
   */
  PreconditionUseMatrix(const MatrixType& M, const function_ptr method);

  /**
   * Execute preconditioning. Calls the function passed to the constructor of
   * this object with the two arguments given here.
   */
  void
  vmult(VectorType& dst, const VectorType& src) const;

private:
  /**
   * Pointer to the matrix in use.
   */
  const MatrixType& matrix;

  /**
   * Pointer to the preconditioning function.
   */
  const function_ptr precondition;
};

/**
 * Base class for other preconditioners. Here, only some common features
 * Jacobi, SOR and SSOR preconditioners are implemented. For preconditioning,
 * refer to derived classes.
 *
 * @author Guido Kanschat, 2000; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionRelaxation : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef typename MatrixType::size_type size_type;

  /**
   * Class for parameters.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor.
     */
    AdditionalData(const double relaxation = 1.);

    /**
     * Relaxation parameter.
     */
    double relaxation;
  };

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object. The relaxation parameter should be larger than
   * zero and smaller than 2 for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType&     A,
             const AdditionalData& parameters = AdditionalData());

  /**
   * Release the matrix and reset its pointer.
   */
  void
  clear();

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  n() const;

protected:
  /**
   * Pointer to the matrix object.
   */
  SmartPointer<const MatrixType, PreconditionRelaxation<MatrixType>> A;

  /**
   * Relaxation parameter.
   */
  double relaxation;
};

/**
 * Jacobi preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have a function
 * <tt>precondition_Jacobi(VectorType&, const VectorType&, double</tt>). This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner:
 *
 * PreconditionJacobi<SparseMatrix<double> > precondition;
 * precondition.initialize (A, PreconditionJacobi<SparseMatrix<double> >::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionJacobi : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * A typedef to the base class AdditionalData.
   */
  typedef
    typename PreconditionRelaxation<MatrixType>::AdditionalData AdditionalData;

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner. Since this is a symmetric preconditioner,
   * this function is the same as vmult().
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  template <class VectorType>
  void
  step(VectorType& x, const VectorType& rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template <class VectorType>
  void
  Tstep(VectorType& x, const VectorType& rhs) const;
};

/**
 * SOR preconditioner using matrix built-in function.
 *
 * Assuming the matrix <i>A = D + L + U</i> is split into its diagonal
 * <i>D</i> as well as the strict lower and upper triangles <i>L</i> and
 * <i>U</i>, then the SOR preconditioner with relaxation parameter <i>r</i> is
 * @f[
 *  P^{-1} = r (D+rL)^{-1}.
 * @f]
 * It is this operator <i>P<sup>-1</sup></i>, which is implemented by vmult()
 * through forward substitution. Analogously, Tvmult() implements the
 * operation of <i>r(D+rU)<sup>-1</sup></i>.
 *
 * The SOR iteration itself can be directly written as
 * @f[
 *  x^{k+1} = x^k - r D^{-1} \bigl(L x^{k+1} + U x^k - b\bigr).
 * @f]
 * Using the right hand side <i>b</i> and the previous iterate <i>x</i>, this
 * is the operation implemented by step().
 *
 * The MatrixType class used is required to have functions
 * <tt>precondition_SOR(VectorType&, const VectorType&, double)</tt> and
 * <tt>precondition_TSOR(VectorType&, const VectorType&, double)</tt>. This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 * // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, PreconditionSOR<SparseMatrix<double> >::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * A typedef to the base class AdditionalData.
   */
  typedef
    typename PreconditionRelaxation<MatrixType>::AdditionalData AdditionalData;

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner.
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  template <class VectorType>
  void
  step(VectorType& x, const VectorType& rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template <class VectorType>
  void
  Tstep(VectorType& x, const VectorType& rhs) const;
};

/**
 * SSOR preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have a function
 * <tt>precondition_SSOR(VectorType&, const VectorType&, double)</tt>. This
 * class satisfies the
 * @ref ConceptRelaxationType "relaxation concept".
 *
 * @code
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * // Define and initialize preconditioner
 *
 * PreconditionSSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, PreconditionSSOR<SparseMatrix<double> >::AdditionalData(.6));
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionSSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * A typedef to the base class AdditionalData.
   */
  typedef
    typename PreconditionRelaxation<MatrixType>::AdditionalData AdditionalData;

  /**
   * Declare type for container size.
   */
  typedef typename MatrixType::size_type size_type;

  /**
   * A typedef to the base class.
   */
  typedef PreconditionRelaxation<MatrixType> BaseClass;

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object. The relaxation parameter should be larger than
   * zero and smaller than 2 for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType&                         A,
             const typename BaseClass::AdditionalData& parameters
             = typename BaseClass::AdditionalData());

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner. Since this is a symmetric preconditioner,
   * this function is the same as vmult().
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;

  /**
   * Perform one step of the preconditioned Richardson iteration
   */
  template <class VectorType>
  void
  step(VectorType& x, const VectorType& rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template <class VectorType>
  void
  Tstep(VectorType& x, const VectorType& rhs) const;

private:
  /**
   * An array that stores for each matrix row where the first position after
   * the diagonal is located.
   */
  std::vector<std::size_t> pos_right_of_diagonal;
};

/**
 * Permuted SOR preconditioner using matrix built-in function.  The
 * <tt>MatrixType</tt> class used is required to have functions
 * <tt>PSOR(VectorType&, const VectorType&, double)</tt> and
 * <tt>TPSOR(VectorType&, const VectorType&, double)</tt>.
 *
 * @code
 *     // Declare related objects
 *
 * SparseMatrix<double> A;
 * Vector<double> x;
 * Vector<double> b;
 * SolverCG<> solver(...);
 *
 * //...initialize and build A
 *
 * std::vector<unsigned int> permutation(x.size());
 * std::vector<unsigned int> inverse_permutation(x.size());
 *
 * //...fill permutation and its inverse with reasonable values
 *
 *     // Define and initialize preconditioner
 *
 * PreconditionPSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, permutation, inverse_permutation, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2003; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
template <typename MatrixType = SparseMatrix<double>>
class PreconditionPSOR : public PreconditionRelaxation<MatrixType>
{
public:
  /**
   * Declare type for container size.
   */
  typedef typename MatrixType::size_type size_type;

  /**
   * Parameters for PreconditionPSOR.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor. For the parameters' description, see below.
     *
     * The permutation vectors are stored as a reference. Therefore, it has to
     * be assured that the lifetime of the vector exceeds the lifetime of the
     * preconditioner.
     *
     * The relaxation parameter should be larger than zero and smaller than 2
     * for numerical reasons. It defaults to 1.
     */
    AdditionalData(
      const std::vector<size_type>& permutation,
      const std::vector<size_type>& inverse_permutation,
      const typename PreconditionRelaxation<MatrixType>::AdditionalData&
        parameters
      = typename PreconditionRelaxation<MatrixType>::AdditionalData());

    /**
     * Storage for the permutation vector.
     */
    const std::vector<size_type>& permutation;
    /**
     * Storage for the inverse permutation vector.
     */
    const std::vector<size_type>& inverse_permutation;
    /**
     * Relaxation parameters
     */
    typename PreconditionRelaxation<MatrixType>::AdditionalData parameters;
  };

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object.
   *
   * The permutation vector is stored as a pointer. Therefore, it has to be
   * assured that the lifetime of the vector exceeds the lifetime of the
   * preconditioner.
   *
   * The relaxation parameter should be larger than zero and smaller than 2
   * for numerical reasons. It defaults to 1.
   */
  void
  initialize(const MatrixType&             A,
             const std::vector<size_type>& permutation,
             const std::vector<size_type>& inverse_permutation,
             const typename PreconditionRelaxation<MatrixType>::AdditionalData&
               parameters
             = typename PreconditionRelaxation<MatrixType>::AdditionalData());

  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object.
   *
   * For more detail about possible parameters, see the class documentation
   * and the documentation of the PreconditionPSOR::AdditionalData class.
   *
   * After this function is called the preconditioner is ready to be used
   * (using the <code>vmult</code> function of derived classes).
   */
  void
  initialize(const MatrixType& A, const AdditionalData& additional_data);

  /**
   * Apply preconditioner.
   */
  template <class VectorType>
  void
  vmult(VectorType&, const VectorType&) const;

  /**
   * Apply transpose preconditioner.
   */
  template <class VectorType>
  void
  Tvmult(VectorType&, const VectorType&) const;

private:
  /**
   * Storage for the permutation vector.
   */
  const std::vector<size_type>* permutation;
  /**
   * Storage for the inverse permutation vector.
   */
  const std::vector<size_type>* inverse_permutation;
};

/**
 * Preconditioning with a Chebyshev polynomial for symmetric positive definite
 * matrices. This preconditioner is based on an iteration of an inner
 * preconditioner of type @p PreconditionerType with coefficients that are
 * adapted to optimally cover an eigenvalue range between the largest
 * eigenvalue down to a given lower eigenvalue specified by the optional
 * parameter @p smoothing_range. The typical use case for the preconditioner
 * is a Jacobi preconditioner specified through DiagonalMatrix, which is also
 * the default value for the preconditioner. Note that if the degree variable
 * is set to zero, the Chebyshev iteration corresponds to a Jacobi
 * preconditioner (or the underlying preconditioner type) with relaxation
 * parameter according to the specified smoothing range.
 *
 * Besides the default choice of a pointwise Jacobi preconditioner, this class
 * also allows for more advanced types of preconditioners, for example
 * iterating block-Jacobi preconditioners in DG methods.
 *
 * Apart from the inner preconditioner object, this iteration does not need
 * access to matrix entries, which makes it an ideal ingredient for
 * matrix-free computations. In that context, this class can be used as a
 * multigrid smoother that is trivially %parallel (assuming that matrix-vector
 * products are %parallel and the inner preconditioner is %parallel). Its use
 * is demonstrated in the step-37 tutorial program.
 *
 * <h4>Algorithm execution</h4>
 *
 * The Chebyshev method relies on an estimate of the eigenvalues of the matrix
 * which are computed during the first invocation of vmult(). The algorithm
 * invokes a conjugate gradient solver so symmetry and positive definiteness
 * of the (preconditioned) matrix system are strong requirements. The
 * computation of eigenvalues needs to be deferred until the first vmult()
 * invocation because temporary vectors of the same layout as the source and
 * destination vectors are necessary for these computations and this
 * information gets only available through vmult().
 *
 * The estimation of eigenvalues can also be bypassed by setting
 * PreconditionChebyshev::AdditionalData::eig_cg_n_iterations to zero and
 * providing sensible values for the largest eigenvalues in the field
 * PreconditionChebyshev::AdditionalData::max_eigenvalue. If the range
 * <tt>[max_eigenvalue/smoothing_range, max_eigenvalue]</tt> contains all
 * eigenvalues of the preconditioned matrix system and the degree (i.e.,
 * number of iterations) is high enough, this class can also be used as a
 * direct solver. For an error estimation of the Chebyshev iteration that can
 * be used to determine the number of iteration, see Varga (2009).
 *
 * In order to use Chebyshev as a solver, set the degree to
 * numbers::invalid_unsigned_int to force the automatic computation of the
 * number of iterations needed to reach a given target tolerance. In this
 * case, the target tolerance is read from the variable
 * PreconditionChebyshev::AdditionalData::smoothing_range (it needs to be a
 * number less than one to force any iterations obviously).
 *
 * For details on the algorithm, see section 5.1 of
 * @code{.bib}
 * @Book{Varga2009,
 *   Title       = {Matrix iterative analysis},
 *   Author      = {Varga, R. S.},
 *   Publisher   = {Springer},
 *   Address     = {Berlin},
 *   Edition     = {2nd},
 *   Year        = {2009},
 * }
 * @endcode
 *
 * <h4>Requirements on the templated classes</h4>
 *
 * The class MatrixType must be derived from Subscriptor because a
 * SmartPointer to MatrixType is held in the class. In particular, this means
 * that the matrix object needs to persist during the lifetime of
 * PreconditionChebyshev. The preconditioner is held in a shared_ptr that is
 * copied into the AdditionalData member variable of the class, so the
 * variable used for initialization can safely be discarded after calling
 * initialize(). Both the matrix and the preconditioner need to provide @p
 * vmult functions for the matrix-vector product and @p m functions for
 * accessing the number of rows in the (square) matrix. Furthermore, the
 * matrix must provide <tt>el(i,i)</tt> methods for accessing the matrix
 * diagonal in case the preconditioner type is a diagonal matrix. Even though
 * it is highly recommended to pass the inverse diagonal entries inside a
 * separate preconditioner object for implementing the Jacobi method (which is
 * the only possible way to operate this class when computing in %parallel
 * with MPI because there is no knowledge about the locally stored range of
 * entries that would be needed from the matrix alone), there is a backward
 * compatibility function that can extract the diagonal in case of a serial
 * computation.
 *
 * @author Martin Kronbichler, 2009, 2016; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
template <typename MatrixType         = SparseMatrix<double>,
          typename VectorType         = Vector<double>,
          typename PreconditionerType = DiagonalMatrix<VectorType>>
class PreconditionChebyshev : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  // avoid warning about use of deprecated variables

  /**
   * Standardized data struct to pipe additional parameters to the
   * preconditioner.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData(const unsigned int degree              = 0,
                   const double       smoothing_range     = 0.,
                   const bool         nonzero_starting    = false,
                   const unsigned int eig_cg_n_iterations = 8,
                   const double       eig_cg_residual     = 1e-2,
                   const double       max_eigenvalue      = 1);

    /**
     * This determines the degree of the Chebyshev polynomial. The degree of
     * the polynomial gives the number of matrix-vector products to be
     * performed for one application of the vmult() operation. Degree zero
     * corresponds to a damped Jacobi method.
     *
     * If the degree is set to numbers::invalid_unsigned_int, the algorithm
     * will automatically determine the number of necessary iterations based
     * on the usual Chebyshev error formula as mentioned in the discussion of
     * the main class.
     */
    unsigned int degree;

    /**
     * This sets the range between the largest eigenvalue in the matrix and
     * the smallest eigenvalue to be treated. If the parameter is set to a
     * number less than 1, an estimate for the largest and for the smallest
     * eigenvalue will be calculated internally. For a smoothing range larger
     * than one, the Chebyshev polynomial will act in the interval
     * $[\lambda_\mathrm{max}/ \tt{smoothing\_range}, \lambda_\mathrm{max}]$,
     * where $\lambda_\mathrm{max}$ is an estimate of the maximum eigenvalue
     * of the matrix. A choice of <tt>smoothing_range</tt> between 5 and 20 is
     * useful in case the preconditioner is used as a smoother in multigrid.
     */
    double smoothing_range;

    /**
     * When this flag is set to <tt>true</tt>, it enables the method
     * <tt>vmult(dst, src)</tt> to use non-zero data in the vector
     * <tt>dst</tt>, appending to it the Chebyshev corrections. This can be
     * useful in some situations (e.g. when used for high-frequency error
     * smoothing in a multigrid algorithm), but not the way the solver classes
     * expect a preconditioner to work (where one ignores the content in
     * <tt>dst</tt> for the preconditioner application).
     *
     * @deprecated For non-zero starting, use the step() and Tstep()
     * interfaces, whereas vmult() provides the preconditioner interface.
     */
    bool nonzero_starting DEAL_II_DEPRECATED;

    /**
     * Maximum number of CG iterations performed for finding the maximum
     * eigenvalue. If set to zero, no computations are performed and the
     * eigenvalues according to the given input are used instead.
     */
    unsigned int eig_cg_n_iterations;

    /**
     * Tolerance for CG iterations performed for finding the maximum
     * eigenvalue.
     */
    double eig_cg_residual;

    /**
     * Maximum eigenvalue to work with. Only in effect if @p
     * eig_cg_n_iterations is set to zero, otherwise this parameter is
     * ignored.
     */
    double max_eigenvalue;

    /**
     * Stores the inverse of the diagonal of the underlying matrix.
     *
     * @deprecated Set the variable @p preconditioner defined below instead.
     */
    VectorType matrix_diagonal_inverse DEAL_II_DEPRECATED;

    /**
     * Stores the preconditioner object that the Chebyshev is wrapped around.
     */
    std::shared_ptr<PreconditionerType> preconditioner;
  };

  PreconditionChebyshev();

  /**
   * Initialize function. Takes the matrix which is used to form the
   * preconditioner, and additional flags if there are any. This function
   * works only if the input matrix has an operator <tt>el(i,i)</tt> for
   * accessing all the elements in the diagonal. Alternatively, the diagonal
   * can be supplied with the help of the AdditionalData field.
   *
   * This function calculates an estimate of the eigenvalue range of the
   * matrix weighted by its diagonal using a modified CG iteration in case the
   * given number of iterations is positive.
   */
  void
  initialize(const MatrixType&     matrix,
             const AdditionalData& additional_data = AdditionalData());

  /**
   * Compute the action of the preconditioner on <tt>src</tt>, storing the
   * result in <tt>dst</tt>.
   */
  void
  vmult(VectorType& dst, const VectorType& src) const;

  /**
   * Compute the action of the transposed preconditioner on <tt>src</tt>,
   * storing the result in <tt>dst</tt>.
   */
  void
  Tvmult(VectorType& dst, const VectorType& src) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  void
  step(VectorType& dst, const VectorType& src) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  void
  Tstep(VectorType& dst, const VectorType& src) const;

  /**
   * Resets the preconditioner.
   */
  void
  clear();

  /**
   * Return the dimension of the codomain (or range) space. Note that the
   * matrix is of dimension $m \times n$.
   */
  size_type
  m() const;

  /**
   * Return the dimension of the domain space. Note that the matrix is of
   * dimension $m \times n$.
   */
  size_type
  n() const;

private:
  /**
   * A pointer to the underlying matrix.
   */
  SmartPointer<
    const MatrixType,
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>>
    matrix_ptr;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType update1;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType update2;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VectorType update3;

  /**
   * Stores the additional data passed to the initialize function, obtained
   * through a copy operation.
   */
  AdditionalData data;

  /**
   * Average of the largest and smallest eigenvalue under consideration.
   */
  double theta;

  /**
   * Half the interval length between the largest and smallest eigenvalue
   * under consideration.
   */
  double delta;

  /**
   * Stores whether the preconditioner has been set up and eigenvalues have
   * been computed.
   */
  bool eigenvalues_are_initialized;

  /**
   * A mutex to avoid that multiple vmult() invocations by different threads
   * overwrite the temporary vectors.
   */
  mutable Threads::Mutex mutex;

  /**
   * Runs the inner loop of the Chebyshev preconditioner that is the same for
   * vmult() and step() methods.
   */
  void
  do_chebyshev_loop(VectorType& dst, const VectorType& src) const;

  /**
   * Runs the inner loop of the Chebyshev preconditioner that is the same for
   * vmult() and step() methods. Uses a separate function to not force users
   * to provide both vmult() and Tvmult() in case only one variant is
   * requested in subsequent calls.
   */
  void
  do_transpose_chebyshev_loop(VectorType& dst, const VectorType& src) const;

  /**
   * Initializes the factors theta and delta based on an eigenvalue
   * computation. If the user set provided values for the largest eigenvalue
   * in AdditionalData, no computation is performed and the information given
   * by the user is used.
   */
  void
  estimate_eigenvalues(const VectorType& src) const;
};

/*@}*/
/* ---------------------------------- Inline functions ------------------- */

#ifndef DOXYGEN

inline PreconditionIdentity::PreconditionIdentity() : n_rows(0), n_columns(0)
{}

template <typename MatrixType>
inline void
PreconditionIdentity::initialize(const MatrixType& matrix,
                                 const PreconditionIdentity::AdditionalData&)
{
  n_rows    = matrix.m();
  n_columns = matrix.n();
}

template <class VectorType>
inline void
PreconditionIdentity::vmult(VectorType& dst, const VectorType& src) const
{
  dst = src;
}

template <class VectorType>
inline void
PreconditionIdentity::Tvmult(VectorType& dst, const VectorType& src) const
{
  dst = src;
}

template <class VectorType>
inline void
PreconditionIdentity::vmult_add(VectorType& dst, const VectorType& src) const
{
  dst += src;
}

template <class VectorType>
inline void
PreconditionIdentity::Tvmult_add(VectorType& dst, const VectorType& src) const
{
  dst += src;
}

inline PreconditionIdentity::size_type
PreconditionIdentity::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionIdentity::size_type
PreconditionIdentity::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

inline PreconditionRichardson::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}

inline PreconditionRichardson::PreconditionRichardson()
  : relaxation(0), n_rows(0), n_columns(0)
{
  AdditionalData add_data;
  relaxation = add_data.relaxation;
}

inline void
PreconditionRichardson::initialize(
  const PreconditionRichardson::AdditionalData& parameters)
{
  relaxation = parameters.relaxation;
}

template <typename MatrixType>
inline void
PreconditionRichardson::initialize(
  const MatrixType&                             matrix,
  const PreconditionRichardson::AdditionalData& parameters)
{
  relaxation = parameters.relaxation;
  n_rows     = matrix.m();
  n_columns  = matrix.n();
}

template <class VectorType>
inline void
PreconditionRichardson::vmult(VectorType& dst, const VectorType& src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}

template <class VectorType>
inline void
PreconditionRichardson::Tvmult(VectorType& dst, const VectorType& src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.equ(relaxation, src);
}

template <class VectorType>
inline void
PreconditionRichardson::vmult_add(VectorType& dst, const VectorType& src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}

template <class VectorType>
inline void
PreconditionRichardson::Tvmult_add(VectorType& dst, const VectorType& src) const
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionRichardson and VectorType must have the same size_type.");

  dst.add(relaxation, src);
}

inline PreconditionRichardson::size_type
PreconditionRichardson::m() const
{
  Assert(n_rows != 0, ExcNotInitialized());
  return n_rows;
}

inline PreconditionRichardson::size_type
PreconditionRichardson::n() const
{
  Assert(n_columns != 0, ExcNotInitialized());
  return n_columns;
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionRelaxation<MatrixType>::initialize(const MatrixType&     rA,
                                               const AdditionalData& parameters)
{
  A          = &rA;
  relaxation = parameters.relaxation;
}

template <typename MatrixType>
inline void
PreconditionRelaxation<MatrixType>::clear()
{
  A = nullptr;
}

template <typename MatrixType>
inline typename PreconditionRelaxation<MatrixType>::size_type
PreconditionRelaxation<MatrixType>::m() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->m();
}

template <typename MatrixType>
inline typename PreconditionRelaxation<MatrixType>::size_type
PreconditionRelaxation<MatrixType>::n() const
{
  Assert(A != nullptr, ExcNotInitialized());
  return A->n();
}

//---------------------------------------------------------------------------

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::vmult(VectorType&       dst,
                                      const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_Jacobi(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::Tvmult(VectorType&       dst,
                                       const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_Jacobi(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::step(VectorType&       dst,
                                     const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->Jacobi_step(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionJacobi<MatrixType>::Tstep(VectorType&       dst,
                                      const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionJacobi<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionJacobi and VectorType must have the same size_type.");

  step(dst, src);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::vmult(VectorType& dst, const VectorType& src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SOR(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::Tvmult(VectorType&       dst,
                                    const VectorType& src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_TSOR(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::step(VectorType& dst, const VectorType& src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->SOR_step(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSOR<MatrixType>::Tstep(VectorType& dst, const VectorType& src) const
{
  static_assert(std::is_same<typename PreconditionSOR<MatrixType>::size_type,
                             typename VectorType::size_type>::value,
                "PreconditionSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->TSOR_step(dst, src, this->relaxation);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionSSOR<MatrixType>::initialize(
  const MatrixType&                         rA,
  const typename BaseClass::AdditionalData& parameters)
{
  this->PreconditionRelaxation<MatrixType>::initialize(rA, parameters);

  // in case we have a SparseMatrix class, we can extract information about
  // the diagonal.
  const SparseMatrix<typename MatrixType::value_type>* mat
    = dynamic_cast<const SparseMatrix<typename MatrixType::value_type>*>(
      &*this->A);

  // calculate the positions first after the diagonal.
  if(mat != nullptr)
    {
      const size_type n = this->A->n();
      pos_right_of_diagonal.resize(n, static_cast<std::size_t>(-1));
      for(size_type row = 0; row < n; ++row)
        {
          // find the first element in this line which is on the right of the
          // diagonal.  we need to precondition with the elements on the left
          // only. note: the first entry in each line denotes the diagonal
          // element, which we need not check.
          typename SparseMatrix<typename MatrixType::value_type>::const_iterator
            it
            = mat->begin(row) + 1;
          for(; it < mat->end(row); ++it)
            if(it->column() > row)
              break;
          pos_right_of_diagonal[row] = it - mat->begin();
        }
    }
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::vmult(VectorType&       dst,
                                    const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SSOR(dst, src, this->relaxation, pos_right_of_diagonal);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::Tvmult(VectorType&       dst,
                                     const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->precondition_SSOR(dst, src, this->relaxation, pos_right_of_diagonal);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::step(VectorType& dst, const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  this->A->SSOR_step(dst, src, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionSSOR<MatrixType>::Tstep(VectorType&       dst,
                                    const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionSSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionSSOR and VectorType must have the same size_type.");

  step(dst, src);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(
  const MatrixType&                                                  rA,
  const std::vector<size_type>&                                      p,
  const std::vector<size_type>&                                      ip,
  const typename PreconditionRelaxation<MatrixType>::AdditionalData& parameters)
{
  permutation         = &p;
  inverse_permutation = &ip;
  PreconditionRelaxation<MatrixType>::initialize(rA, parameters);
}

template <typename MatrixType>
inline void
PreconditionPSOR<MatrixType>::initialize(const MatrixType&     A,
                                         const AdditionalData& additional_data)
{
  initialize(A,
             additional_data.permutation,
             additional_data.inverse_permutation,
             additional_data.parameters);
}

template <typename MatrixType>
template <typename VectorType>
inline void
PreconditionPSOR<MatrixType>::vmult(VectorType&       dst,
                                    const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionPSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionPSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  dst = src;
  this->A->PSOR(dst, *permutation, *inverse_permutation, this->relaxation);
}

template <typename MatrixType>
template <class VectorType>
inline void
PreconditionPSOR<MatrixType>::Tvmult(VectorType&       dst,
                                     const VectorType& src) const
{
  static_assert(
    std::is_same<typename PreconditionPSOR<MatrixType>::size_type,
                 typename VectorType::size_type>::value,
    "PreconditionPSOR and VectorType must have the same size_type.");

  Assert(this->A != nullptr, ExcNotInitialized());
  dst = src;
  this->A->TPSOR(dst, *permutation, *inverse_permutation, this->relaxation);
}

template <typename MatrixType>
PreconditionPSOR<MatrixType>::AdditionalData::AdditionalData(
  const std::vector<size_type>& permutation,
  const std::vector<size_type>& inverse_permutation,
  const typename PreconditionRelaxation<MatrixType>::AdditionalData& parameters)
  : permutation(permutation),
    inverse_permutation(inverse_permutation),
    parameters(parameters)
{}

//---------------------------------------------------------------------------

template <typename MatrixType, class VectorType>
PreconditionUseMatrix<MatrixType, VectorType>::PreconditionUseMatrix(
  const MatrixType&  M,
  const function_ptr method)
  : matrix(M), precondition(method)
{}

template <typename MatrixType, class VectorType>
void
PreconditionUseMatrix<MatrixType, VectorType>::vmult(
  VectorType&       dst,
  const VectorType& src) const
{
  (matrix.*precondition)(dst, src);
}

//---------------------------------------------------------------------------

template <typename MatrixType>
inline PreconditionRelaxation<MatrixType>::AdditionalData::AdditionalData(
  const double relaxation)
  : relaxation(relaxation)
{}

//---------------------------------------------------------------------------

namespace internal
{
  namespace PreconditionChebyshevImplementation
  {
    // for deal.II vectors, perform updates for Chebyshev preconditioner all
    // at once to reduce memory transfer. Here, we select between general
    // vectors and deal.II vectors where we expand the loop over the (local)
    // size of the vector

    // generic part for non-deal.II vectors
    template <typename VectorType, typename PreconditionerType>
    inline void
    vector_updates(const VectorType&         src,
                   const PreconditionerType& preconditioner,
                   const bool                start_zero,
                   const double              factor1,
                   const double              factor2,
                   VectorType&               update1,
                   VectorType&               update2,
                   VectorType&               update3,
                   VectorType&               dst)
    {
      if(start_zero)
        {
          update1.equ(factor2, src);
          preconditioner.vmult(dst, update1);
          update1.equ(-1., dst);
        }
      else
        {
          update2 -= src;
          preconditioner.vmult(update3, update2);
          update2 = update3;
          if(factor1 == 0.)
            update1.equ(factor2, update2);
          else
            update1.sadd(factor1, factor2, update2);
          dst -= update1;
        }
    }

    // worker routine for deal.II vectors. Because of vectorization, we need
    // to put the loop into an extra structure because the virtual function of
    // VectorUpdatesRange prevents the compiler from applying vectorization.
    template <typename Number>
    struct VectorUpdater
    {
      VectorUpdater(const Number* src,
                    const Number* matrix_diagonal_inverse,
                    const bool    start_zero,
                    const Number  factor1,
                    const Number  factor2,
                    Number*       update1,
                    Number*       update2,
                    Number*       dst)
        : src(src),
          matrix_diagonal_inverse(matrix_diagonal_inverse),
          do_startup(factor1 == Number()),
          start_zero(start_zero),
          factor1(factor1),
          factor2(factor2),
          update1(update1),
          update2(update2),
          dst(dst)
      {}

      void
      apply_to_subrange(const std::size_t begin, const std::size_t end) const
      {
        // To circumvent a bug in gcc
        // (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63945), we create copies
        // of the variables factor1 and factor2 and do not check based on
        // factor1.
        const Number factor1 = this->factor1;
        const Number factor2 = this->factor2;
        if(do_startup)
          {
            if(start_zero)
              DEAL_II_OPENMP_SIMD_PRAGMA
            for(std::size_t i = begin; i < end; ++i)
              {
                dst[i]     = factor2 * src[i] * matrix_diagonal_inverse[i];
                update1[i] = -dst[i];
              }
            else DEAL_II_OPENMP_SIMD_PRAGMA for(std::size_t i = begin; i < end;
                                                ++i)
            {
              update1[i] = ((update2[i] - src[i]) * factor2
                            * matrix_diagonal_inverse[i]);
              dst[i] -= update1[i];
            }
          }
        else
          DEAL_II_OPENMP_SIMD_PRAGMA
        for(std::size_t i = begin; i < end; ++i)
          {
            const Number update
              = factor1 * update1[i]
                + factor2
                    * ((update2[i] - src[i]) * matrix_diagonal_inverse[i]);
            update1[i] = update;
            dst[i] -= update;
          }
      }

      const Number*   src;
      const Number*   matrix_diagonal_inverse;
      const bool      do_startup;
      const bool      start_zero;
      const Number    factor1;
      const Number    factor2;
      mutable Number* update1;
      mutable Number* update2;
      mutable Number* dst;
    };

    template <typename Number>
    struct VectorUpdatesRange : public parallel::ParallelForInteger
    {
      VectorUpdatesRange(const VectorUpdater<Number>& updater,
                         const std::size_t            size)
        : updater(updater)
      {
        if(size < internal::VectorImplementation::minimum_parallel_grain_size)
          apply_to_subrange(0, size);
        else
          apply_parallel(
            0,
            size,
            internal::VectorImplementation::minimum_parallel_grain_size);
      }

      ~VectorUpdatesRange() override = default;

      virtual void
      apply_to_subrange(const std::size_t begin,
                        const std::size_t end) const override
      {
        updater.apply_to_subrange(begin, end);
      }

      const VectorUpdater<Number>& updater;
    };

    // selection for diagonal matrix around deal.II vector
    template <typename Number>
    inline void
    vector_updates(const ::dealii::Vector<Number>&                 src,
                   const DiagonalMatrix<::dealii::Vector<Number>>& jacobi,
                   const bool                                      start_zero,
                   const double                                    factor1,
                   const double                                    factor2,
                   ::dealii::Vector<Number>&                       update1,
                   ::dealii::Vector<Number>&                       update2,
                   ::dealii::Vector<Number>&,
                   ::dealii::Vector<Number>& dst)
    {
      VectorUpdater<Number> upd(src.begin(),
                                jacobi.get_vector().begin(),
                                start_zero,
                                factor1,
                                factor2,
                                update1.begin(),
                                update2.begin(),
                                dst.begin());
      VectorUpdatesRange<Number>(upd, src.size());
    }

    // selection for diagonal matrix around parallel deal.II vector
    template <typename Number>
    inline void
    vector_updates(
      const LinearAlgebra::distributed::Vector<Number>&                 src,
      const DiagonalMatrix<LinearAlgebra::distributed::Vector<Number>>& jacobi,
      const bool                                  start_zero,
      const double                                factor1,
      const double                                factor2,
      LinearAlgebra::distributed::Vector<Number>& update1,
      LinearAlgebra::distributed::Vector<Number>& update2,
      LinearAlgebra::distributed::Vector<Number>&,
      LinearAlgebra::distributed::Vector<Number>& dst)
    {
      VectorUpdater<Number> upd(src.begin(),
                                jacobi.get_vector().begin(),
                                start_zero,
                                factor1,
                                factor2,
                                update1.begin(),
                                update2.begin(),
                                dst.begin());
      VectorUpdatesRange<Number>(upd, src.local_size());
    }

    template <typename MatrixType,
              typename VectorType,
              typename PreconditionerType>
    inline void
    initialize_preconditioner(
      const MatrixType&                    matrix,
      std::shared_ptr<PreconditionerType>& preconditioner,
      VectorType&)
    {
      (void) matrix;
      (void) preconditioner;
      AssertThrow(preconditioner.get() != nullptr, ExcNotInitialized());
    }

    template <typename MatrixType, typename VectorType>
    inline void
    initialize_preconditioner(
      const MatrixType&                            matrix,
      std::shared_ptr<DiagonalMatrix<VectorType>>& preconditioner,
      VectorType&                                  diagonal_inverse)
    {
      if(preconditioner.get() == nullptr || preconditioner->m() != matrix.m())
        {
          if(preconditioner.get() == nullptr)
            preconditioner = std::make_shared<DiagonalMatrix<VectorType>>();

          Assert(
            preconditioner->m() == 0,
            ExcMessage(
              "Preconditioner appears to be initialized but not sized correctly"));

          // Check if we can initialize from vector that then gets set to zero
          // as the matrix will own the memory
          preconditioner->reinit(diagonal_inverse);
          {
            VectorType empty_vector;
            diagonal_inverse.reinit(empty_vector);
          }

          // This part only works in serial
          if(preconditioner->m() != matrix.m())
            {
              preconditioner->get_vector().reinit(matrix.m());
              for(typename VectorType::size_type i = 0; i < matrix.m(); ++i)
                preconditioner->get_vector()(i) = 1. / matrix.el(i, i);
            }
        }
    }

    template <typename VectorType>
    void
    set_initial_guess(VectorType& vector)
    {
      vector = 1. / std::sqrt(static_cast<double>(vector.size()));
      if(vector.locally_owned_elements().is_element(0))
        vector(0) = 0.;
    }

    template <typename Number>
    void
    set_initial_guess(::dealii::Vector<Number>& vector)
    {
      // Choose a high-frequency mode consisting of numbers between 0 and 1
      // that is cheap to compute (cheaper than random numbers) but avoids
      // obviously re-occurring numbers in multi-component systems by choosing
      // a period of 11
      for(unsigned int i = 0; i < vector.size(); ++i)
        vector(i) = i % 11;

      const Number mean_value = vector.mean_value();
      vector.add(-mean_value);
    }

    template <typename Number>
    void
    set_initial_guess(
      ::dealii::LinearAlgebra::distributed::Vector<Number>& vector)
    {
      // Choose a high-frequency mode consisting of numbers between 0 and 1
      // that is cheap to compute (cheaper than random numbers) but avoids
      // obviously re-occurring numbers in multi-component systems by choosing
      // a period of 11.
      // Make initial guess robust with respect to number of processors
      // by operating on the global index.
      types::global_dof_index first_local_range = 0;
      if(!vector.locally_owned_elements().is_empty())
        first_local_range = vector.locally_owned_elements().nth_index_in_set(0);
      for(unsigned int i = 0; i < vector.local_size(); ++i)
        vector.local_element(i) = (i + first_local_range) % 11;

      const Number mean_value = vector.mean_value();
      vector.add(-mean_value);
    }

    struct EigenvalueTracker
    {
    public:
      void
      slot(const std::vector<double>& eigenvalues)
      {
        values = eigenvalues;
      }

      std::vector<double> values;
    };
  } // namespace PreconditionChebyshevImplementation
} // namespace internal

// avoid warning about deprecated variable nonzero_starting

template <typename MatrixType, class VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  AdditionalData::AdditionalData(const unsigned int degree,
                                 const double       smoothing_range,
                                 const bool         nonzero_starting,
                                 const unsigned int eig_cg_n_iterations,
                                 const double       eig_cg_residual,
                                 const double       max_eigenvalue)
  : degree(degree),
    smoothing_range(smoothing_range),
    nonzero_starting(nonzero_starting),
    eig_cg_n_iterations(eig_cg_n_iterations),
    eig_cg_residual(eig_cg_residual),
    max_eigenvalue(max_eigenvalue)
{}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  PreconditionChebyshev()
  : theta(1.), delta(1.), eigenvalues_are_initialized(false)
{
  static_assert(
    std::is_same<size_type, typename VectorType::size_type>::value,
    "PreconditionChebyshev and VectorType must have the same size_type.");
}

// avoid warning about deprecated variable AdditionalData::matrix_diagonal_inverse

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::initialize(
  const MatrixType&     matrix,
  const AdditionalData& additional_data)
{
  matrix_ptr = &matrix;
  data       = additional_data;
  internal::PreconditionChebyshevImplementation::initialize_preconditioner(
    matrix, data.preconditioner, data.matrix_diagonal_inverse);
  eigenvalues_are_initialized = false;
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::clear()
{
  eigenvalues_are_initialized = false;
  theta = delta = 1.0;
  matrix_ptr    = nullptr;
  {
    VectorType empty_vector;
    data.matrix_diagonal_inverse.reinit(empty_vector);
    update1.reinit(empty_vector);
    update2.reinit(empty_vector);
    update3.reinit(empty_vector);
  }
  data.preconditioner.reset();
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  estimate_eigenvalues(const VectorType& src) const
{
  Assert(eigenvalues_are_initialized == false, ExcInternalError());
  Assert(data.preconditioner.get() != nullptr, ExcNotInitialized());

  update1.reinit(src);
  update2.reinit(src, true);

  // calculate largest eigenvalue using a hand-tuned CG iteration on the
  // matrix weighted by its diagonal. we start with a vector that consists of
  // ones only, weighted by the length.
  double max_eigenvalue, min_eigenvalue;
  if(data.eig_cg_n_iterations > 0)
    {
      Assert(
        data.eig_cg_n_iterations > 2,
        ExcMessage("Need to set at least two iterations to find eigenvalues."));

      // set a very strict tolerance to force at least two iterations
      ReductionControl control(
        data.eig_cg_n_iterations,
        std::sqrt(
          std::numeric_limits<typename VectorType::value_type>::epsilon()),
        1e-10,
        false,
        false);

      internal::PreconditionChebyshevImplementation::EigenvalueTracker
                           eigenvalue_tracker;
      SolverCG<VectorType> solver(control);
      solver.connect_eigenvalues_slot(std::bind(
        &internal::PreconditionChebyshevImplementation::EigenvalueTracker::slot,
        &eigenvalue_tracker,
        std::placeholders::_1));

      // set an initial guess which is close to the constant vector but where
      // one entry is different to trigger high frequencies
      internal::PreconditionChebyshevImplementation::set_initial_guess(update2);

      try
        {
          solver.solve(*matrix_ptr, update1, update2, *data.preconditioner);
        }
      catch(SolverControl::NoConvergence&)
        {}

      // read the eigenvalues from the attached eigenvalue tracker
      if(eigenvalue_tracker.values.empty())
        min_eigenvalue = max_eigenvalue = 1;
      else
        {
          min_eigenvalue = eigenvalue_tracker.values.front();

          // include a safety factor since the CG method will in general not
          // be converged
          max_eigenvalue = 1.2 * eigenvalue_tracker.values.back();
        }
    }
  else
    {
      max_eigenvalue = data.max_eigenvalue;
      min_eigenvalue = data.max_eigenvalue / data.smoothing_range;
    }

  const double alpha = (data.smoothing_range > 1. ?
                          max_eigenvalue / data.smoothing_range :
                          std::min(0.9 * max_eigenvalue, min_eigenvalue));

  // in case the user set the degree to invalid unsigned int, we have to
  // determine the number of necessary iterations from the Chebyshev error
  // estimate, given the target tolerance specified by smoothing_range. This
  // estimate is based on the error formula given in section 5.1 of
  // R. S. Varga, Matrix iterative analysis, 2nd ed., Springer, 2009
  if(data.degree == numbers::invalid_unsigned_int)
    {
      const double actual_range = max_eigenvalue / alpha;
      const double sigma        = (1. - std::sqrt(1. / actual_range))
                           / (1. + std::sqrt(1. / actual_range));
      const double eps = data.smoothing_range;
      const_cast<
        PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>*>(
        this)
        ->data.degree
        = 1
          + std::log(1. / eps + std::sqrt(1. / eps / eps - 1))
              / std::log(1. / sigma);
    }

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>*>(this)
    ->delta
    = (max_eigenvalue - alpha) * 0.5;
  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>*>(this)
    ->theta
    = (max_eigenvalue + alpha) * 0.5;

  // We do not need the third auxiliary vector in case we have a
  // DiagonalMatrix as preconditioner and use deal.II's own vectors
  if(std::is_same<PreconditionerType, DiagonalMatrix<VectorType>>::value
       == false
     || (std::is_same<VectorType,
                      dealii::Vector<typename VectorType::value_type>>::value
           == false
         && std::is_same<VectorType,
                         LinearAlgebra::distributed::Vector<
                           typename VectorType::value_type>>::value
              == false))
    update3.reinit(src, true);

  const_cast<
    PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>*>(this)
    ->eigenvalues_are_initialized
    = true;
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  do_chebyshev_loop(VectorType& dst, const VectorType& src) const
{
  // if delta is zero, we do not need to iterate because the updates will be
  // zero
  if(std::abs(delta) < 1e-40)
    return;

  double rhok = delta / theta, sigma = theta / delta;
  for(unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->vmult(update2, dst);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        src,
        *data.preconditioner,
        false,
        factor1,
        factor2,
        update1,
        update2,
        update3,
        dst);
    }
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::
  do_transpose_chebyshev_loop(VectorType& dst, const VectorType& src) const
{
  double rhok = delta / theta, sigma = theta / delta;
  for(unsigned int k = 0; k < data.degree; ++k)
    {
      matrix_ptr->Tvmult(update2, dst);
      const double rhokp   = 1. / (2. * sigma - rhok);
      const double factor1 = rhokp * rhok, factor2 = 2. * rhokp / delta;
      rhok = rhokp;
      internal::PreconditionChebyshevImplementation::vector_updates(
        src,
        *data.preconditioner,
        false,
        factor1,
        factor2,
        update1,
        update2,
        update3,
        dst);
    }
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::vmult(
  VectorType&       dst,
  const VectorType& src) const
{
  Threads::Mutex::ScopedLock lock(mutex);
  if(eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  internal::PreconditionChebyshevImplementation::vector_updates(
    src,
    *data.preconditioner,
    true,
    0.,
    1. / theta,
    update1,
    update2,
    update3,
    dst);

  do_chebyshev_loop(dst, src);
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tvmult(
  VectorType&       dst,
  const VectorType& src) const
{
  Threads::Mutex::ScopedLock lock(mutex);
  if(eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  internal::PreconditionChebyshevImplementation::vector_updates(
    src,
    *data.preconditioner,
    true,
    0.,
    1. / theta,
    update1,
    update2,
    update3,
    dst);

  do_transpose_chebyshev_loop(dst, src);
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::step(
  VectorType&       dst,
  const VectorType& src) const
{
  Threads::Mutex::ScopedLock lock(mutex);
  if(eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  matrix_ptr->vmult(update2, dst);
  internal::PreconditionChebyshevImplementation::vector_updates(
    src,
    *data.preconditioner,
    false,
    0.,
    1. / theta,
    update1,
    update2,
    update3,
    dst);

  do_chebyshev_loop(dst, src);
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline void
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::Tstep(
  VectorType&       dst,
  const VectorType& src) const
{
  Threads::Mutex::ScopedLock lock(mutex);
  if(eigenvalues_are_initialized == false)
    estimate_eigenvalues(src);

  matrix_ptr->Tvmult(update2, dst);
  internal::PreconditionChebyshevImplementation::vector_updates(
    src,
    *data.preconditioner,
    false,
    0.,
    1. / theta,
    update1,
    update2,
    update3,
    dst);

  do_transpose_chebyshev_loop(dst, src);
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::m() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->m();
}

template <typename MatrixType, typename VectorType, typename PreconditionerType>
inline typename PreconditionChebyshev<MatrixType,
                                      VectorType,
                                      PreconditionerType>::size_type
PreconditionChebyshev<MatrixType, VectorType, PreconditionerType>::n() const
{
  Assert(matrix_ptr != nullptr, ExcNotInitialized());
  return matrix_ptr->n();
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
