// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2015 by the deal.II authors
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

#ifndef dealii__precondition_h
#define dealii__precondition_h

// This file contains simple preconditioners.

#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/parallel.h>
#include <deal.II/base/template_constraints.h>
#include <deal.II/lac/tridiagonal_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/vector_memory.h>

DEAL_II_NAMESPACE_OPEN

// forward declarations

template <typename number> class Vector;
template <typename number> class SparseMatrix;
namespace parallel
{
  namespace distributed
  {
    template <typename number> class Vector;
  }
}



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
 * @author Guido Kanschat, 1999
 */
class PreconditionIdentity : public Subscriptor
{
public:
  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData () {}
  };

  /**
   * The matrix argument is ignored and here just for compatibility with more
   * complex preconditioners.
   */
  template <class MATRIX>
  void initialize (const MATRIX         &matrix,
                   const AdditionalData &additional_data = AdditionalData());

  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply preconditioner, adding to the previous value.
   */
  template<class VECTOR>
  void vmult_add (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template<class VECTOR>
  void Tvmult_add (VECTOR &, const VECTOR &) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void clear () {}
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
 * @author Guido Kanschat, 2005
 */
class PreconditionRichardson : public Subscriptor
{
public:
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
    AdditionalData (const double relaxation = 1.);

    /**
     * Relaxation parameter.
     */
    double relaxation;
  };

  /**
   * Constructor, sets the relaxation parameter to its default.
   */
  PreconditionRichardson();

  /**
   * Change the relaxation parameter.
   */
  void initialize (const AdditionalData &parameters);

  /**
   * Change the relaxation parameter in a way consistent with other
   * preconditioners. The matrix argument is ignored and here just for
   * compatibility with more complex preconditioners.
   */
  template <class MATRIX>
  void initialize (const MATRIX &,
                   const AdditionalData &parameters);

  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner. Since this is the identity, this function
   * is the same as vmult().
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;
  /**
   * Apply preconditioner, adding to the previous value.
   */
  template<class VECTOR>
  void vmult_add (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner, adding. Since this is the identity, this
   * function is the same as vmult_add().
   */
  template<class VECTOR>
  void Tvmult_add (VECTOR &, const VECTOR &) const;

  /**
   * This function is only present to provide the interface of a
   * preconditioner to be handed to a smoother.  This does nothing.
   */
  void clear () {}

private:
  /**
   * The relaxation parameter multiplied with the vectors.
   */
  double relaxation;
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
template<class MATRIX = SparseMatrix<double>, class VECTOR = Vector<double> >
class PreconditionUseMatrix : public Subscriptor
{
public:
  /**
   * Type of the preconditioning function of the matrix.
   */
  typedef void ( MATRIX::* function_ptr)(VECTOR &, const VECTOR &) const;

  /**
   * Constructor.  This constructor stores a reference to the matrix object
   * for later use and selects a preconditioning method, which must be a
   * member function of that matrix.
   */
  PreconditionUseMatrix(const MATRIX      &M,
                        const function_ptr method);

  /**
   * Execute preconditioning. Calls the function passed to the constructor of
   * this object with the two arguments given here.
   */
  void vmult (VECTOR       &dst,
              const VECTOR &src) const;

private:
  /**
   * Pointer to the matrix in use.
   */
  const MATRIX &matrix;

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
 * @author Guido Kanschat, 2000
 */
template<class MATRIX = SparseMatrix<double> >
class PreconditionRelaxation : public Subscriptor
{
public:
  /**
   * Class for parameters.
   */
  class AdditionalData
  {
  public:
    /**
     * Constructor.
     */
    AdditionalData (const double relaxation = 1.);

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
  void initialize (const MATRIX &A,
                   const AdditionalData &parameters = AdditionalData());

  /**
   * Release the matrix and reset its pointer.
   */
  void clear();

protected:
  /**
   * Pointer to the matrix object.
   */
  SmartPointer<const MATRIX, PreconditionRelaxation<MATRIX> > A;

  /**
   * Relaxation parameter.
   */
  double relaxation;
};



/**
 * Jacobi preconditioner using matrix built-in function.  The <tt>MATRIX</tt>
 * class used is required to have a function <tt>precondition_Jacobi(VECTOR&,
 * const VECTOR&, double</tt>)
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
 *     // Define and initialize preconditioner
 *
 * PreconditionJacobi<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionJacobi : public PreconditionRelaxation<MATRIX>
{
public:
  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner. Since this is a symmetric preconditioner,
   * this function is the same as vmult().
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  template<class VECTOR>
  void step (VECTOR &x, const VECTOR &rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template<class VECTOR>
  void Tstep (VECTOR &x, const VECTOR &rhs) const;
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
 * The MATRIX class used is required to have functions
 * <tt>precondition_SOR(VECTOR&, const VECTOR&, double)</tt> and
 * <tt>precondition_TSOR(VECTOR&, const VECTOR&, double)</tt>.
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
 *     // Define and initialize preconditioner
 *
 * PreconditionSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionSOR : public PreconditionRelaxation<MATRIX>
{
public:
  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner.
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;

  /**
   * Perform one step of the preconditioned Richardson iteration.
   */
  template<class VECTOR>
  void step (VECTOR &x, const VECTOR &rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template<class VECTOR>
  void Tstep (VECTOR &x, const VECTOR &rhs) const;
};



/**
 * SSOR preconditioner using matrix built-in function.  The <tt>MATRIX</tt>
 * class used is required to have a function <tt>precondition_SSOR(VECTOR&,
 * const VECTOR&, double)</tt>
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
 *     // Define and initialize preconditioner
 *
 * PreconditionSSOR<SparseMatrix<double> > precondition;
 * precondition.initialize (A, .6);
 *
 * solver.solve (A, x, b, precondition);
 * @endcode
 *
 * @author Guido Kanschat, 2000
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionSSOR : public PreconditionRelaxation<MATRIX>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * A typedef to the base class.
   */
  typedef PreconditionRelaxation<MATRIX> BaseClass;


  /**
   * Initialize matrix and relaxation parameter. The matrix is just stored in
   * the preconditioner object. The relaxation parameter should be larger than
   * zero and smaller than 2 for numerical reasons. It defaults to 1.
   */
  void initialize (const MATRIX &A,
                   const typename BaseClass::AdditionalData &parameters = typename BaseClass::AdditionalData());

  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner. Since this is a symmetric preconditioner,
   * this function is the same as vmult().
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;


  /**
   * Perform one step of the preconditioned Richardson iteration
   */
  template<class VECTOR>
  void step (VECTOR &x, const VECTOR &rhs) const;

  /**
   * Perform one transposed step of the preconditioned Richardson iteration.
   */
  template<class VECTOR>
  void Tstep (VECTOR &x, const VECTOR &rhs) const;

private:
  /**
   * An array that stores for each matrix row where the first position after
   * the diagonal is located.
   */
  std::vector<std::size_t> pos_right_of_diagonal;
};


/**
 * Permuted SOR preconditioner using matrix built-in function.  The
 * <tt>MATRIX</tt> class used is required to have functions <tt>PSOR(VECTOR&,
 * const VECTOR&, double)</tt> and <tt>TPSOR(VECTOR&, const VECTOR&,
 * double)</tt>.
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
 * @author Guido Kanschat, 2003
 */
template <class MATRIX = SparseMatrix<double> >
class PreconditionPSOR : public PreconditionRelaxation<MATRIX>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

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
  void initialize (const MATRIX &A,
                   const std::vector<size_type> &permutation,
                   const std::vector<size_type> &inverse_permutation,
                   const typename PreconditionRelaxation<MATRIX>::AdditionalData &
                   parameters = typename PreconditionRelaxation<MATRIX>::AdditionalData());

  /**
   * Apply preconditioner.
   */
  template<class VECTOR>
  void vmult (VECTOR &, const VECTOR &) const;

  /**
   * Apply transpose preconditioner.
   */
  template<class VECTOR>
  void Tvmult (VECTOR &, const VECTOR &) const;
private:
  /**
   * Storage for the permutation vector.
   */
  const std::vector<size_type> *permutation;
  /**
   * Storage for the inverse permutation vector.
   */
  const std::vector<size_type> *inverse_permutation;
};



/**
 * Preconditioning with a Chebyshev polynomial for symmetric positive definite
 * matrices. This preconditioner is similar to a Jacobi preconditioner if the
 * degree variable is set to one, otherwise some higher order polynomial
 * corrections are used. This preconditioner needs access to the diagonal of
 * the matrix it acts on and needs a respective <tt>vmult</tt>
 * implementation. However, it does not need to explicitly know the matrix
 * entries.
 *
 * This class is useful e.g. in multigrid smoother objects, since it is
 * trivially %parallel (assuming that matrix-vector products are %parallel).
 *
 * @author Martin Kronbichler, 2009
 */
template <class MATRIX=SparseMatrix<double>, class VECTOR=Vector<double> >
class PreconditionChebyshev : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Standardized data struct to pipe additional parameters to the
   * preconditioner.
   */
  struct AdditionalData
  {
    /**
     * Constructor.
     */
    AdditionalData (const unsigned int degree              = 0,
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
     */
    unsigned int degree;

    /**
     * This sets the range between the largest eigenvalue in the matrix and
     * the smallest eigenvalue to be treated. If the parameter is zero, an
     * estimate for the largest and for the smallest eigenvalue will be
     * calculated internally. Otherwise, the Chebyshev polynomial will act in
     * the interval $[\lambda_\mathrm{max}/ \tt{smoothing\_range},
     * \lambda_\mathrm{max}]$, where $\lambda_\mathrm{max}$ is an estimate of
     * the maximum eigenvalue of the matrix. A choice of
     * <tt>smoothing_range</tt> between 5 and 20 is useful in case the
     * preconditioner is used as a smoother in multigrid.
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
     */
    bool nonzero_starting;

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
     */
    VECTOR matrix_diagonal_inverse;
  };

  PreconditionChebyshev ();

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
  void initialize (const MATRIX         &matrix,
                   const AdditionalData &additional_data = AdditionalData());

  /**
   * Computes the action of the preconditioner on <tt>src</tt>, storing the
   * result in <tt>dst</tt>.
   */
  void vmult (VECTOR       &dst,
              const VECTOR &src) const;

  /**
   * Computes the action of the transposed preconditioner on <tt>src</tt>,
   * storing the result in <tt>dst</tt>.
   */
  void Tvmult (VECTOR       &dst,
               const VECTOR &src) const;

  /**
   * Resets the preconditioner.
   */
  void clear ();

private:

  /**
   * A pointer to the underlying matrix.
   */
  SmartPointer<const MATRIX,PreconditionChebyshev<MATRIX,VECTOR> > matrix_ptr;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VECTOR update1;

  /**
   * Internal vector used for the <tt>vmult</tt> operation.
   */
  mutable VECTOR update2;

  /**
   * Stores the additional data provided to the initialize function.
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
   * Stores whether the preconditioner has been set up.
   */
  bool is_initialized;
};



/*@}*/
/* ---------------------------------- Inline functions ------------------- */

#ifndef DOXYGEN

template <class MATRIX>
inline void
PreconditionIdentity::initialize (
  const MATRIX &,
  const PreconditionIdentity::AdditionalData &)
{}


template<class VECTOR>
inline void
PreconditionIdentity::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst = src;
}



template<class VECTOR>
inline void
PreconditionIdentity::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst = src;
}

template<class VECTOR>
inline void
PreconditionIdentity::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(src);
}



template<class VECTOR>
inline void
PreconditionIdentity::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(src);
}

//---------------------------------------------------------------------------

inline
PreconditionRichardson::AdditionalData::AdditionalData (
  const double relaxation)
  :
  relaxation(relaxation)
{}


inline
PreconditionRichardson::PreconditionRichardson ()
  :
  relaxation(0)
{
  AdditionalData add_data;
  relaxation=add_data.relaxation;
}



inline void
PreconditionRichardson::initialize (
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template <class MATRIX>
inline void
PreconditionRichardson::initialize (
  const MATRIX &,
  const PreconditionRichardson::AdditionalData &parameters)
{
  relaxation = parameters.relaxation;
}



template<class VECTOR>
inline void
PreconditionRichardson::vmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(relaxation,src);
}



template<class VECTOR>
inline void
PreconditionRichardson::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  dst.equ(relaxation,src);
}

template<class VECTOR>
inline void
PreconditionRichardson::vmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(relaxation,src);
}



template<class VECTOR>
inline void
PreconditionRichardson::Tvmult_add (VECTOR &dst, const VECTOR &src) const
{
  dst.add(relaxation,src);
}

//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionRelaxation<MATRIX>::initialize (const MATRIX &rA,
                                            const AdditionalData &parameters)
{
  A = &rA;
  relaxation = parameters.relaxation;
}


template <class MATRIX>
inline void
PreconditionRelaxation<MATRIX>::clear ()
{
  A = 0;
}


//---------------------------------------------------------------------------

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_Jacobi (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_Jacobi (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->Jacobi_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionJacobi<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  step (dst, src);
}



//---------------------------------------------------------------------------

template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SOR (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_TSOR (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->SOR_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSOR<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->TSOR_step (dst, src, this->relaxation);
}



//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionSSOR<MATRIX>::initialize (const MATRIX &rA,
                                      const typename BaseClass::AdditionalData &parameters)
{
  this->PreconditionRelaxation<MATRIX>::initialize (rA, parameters);

  // in case we have a SparseMatrix class, we can extract information about
  // the diagonal.
  const SparseMatrix<typename MATRIX::value_type> *mat =
    dynamic_cast<const SparseMatrix<typename MATRIX::value_type> *>(&*this->A);

  // calculate the positions first after the diagonal.
  if (mat != 0)
    {
      const size_type n = this->A->n();
      pos_right_of_diagonal.resize(n, static_cast<std::size_t>(-1));
      for (size_type row=0; row<n; ++row)
        {
          // find the first element in this line which is on the right of the
          // diagonal.  we need to precondition with the elements on the left
          // only. note: the first entry in each line denotes the diagonal
          // element, which we need not check.
          typename SparseMatrix<typename MATRIX::value_type>::const_iterator
          it = mat->begin(row)+1;
          for ( ; it < mat->end(row); ++it)
            if (it->column() > row)
              break;
          pos_right_of_diagonal[row] = it - mat->begin();
        }
    }
}


template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SSOR (dst, src, this->relaxation, pos_right_of_diagonal);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->precondition_SSOR (dst, src, this->relaxation, pos_right_of_diagonal);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::step (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  this->A->SSOR_step (dst, src, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionSSOR<MATRIX>::Tstep (VECTOR &dst, const VECTOR &src) const
{
  step (dst, src);
}



//---------------------------------------------------------------------------

template <class MATRIX>
inline void
PreconditionPSOR<MATRIX>::initialize (
  const MATRIX &rA,
  const std::vector<size_type> &p,
  const std::vector<size_type> &ip,
  const typename PreconditionRelaxation<MATRIX>::AdditionalData &parameters)
{
  permutation = &p;
  inverse_permutation = &ip;
  PreconditionRelaxation<MATRIX>::initialize(rA, parameters);
}


template <class MATRIX>
template<class VECTOR>
inline void
PreconditionPSOR<MATRIX>::vmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  dst = src;
  this->A->PSOR (dst, *permutation, *inverse_permutation, this->relaxation);
}



template <class MATRIX>
template<class VECTOR>
inline void
PreconditionPSOR<MATRIX>::Tvmult (VECTOR &dst, const VECTOR &src) const
{
  Assert (this->A!=0, ExcNotInitialized());
  dst = src;
  this->A->TPSOR (dst, *permutation, *inverse_permutation, this->relaxation);
}


//---------------------------------------------------------------------------


template<class MATRIX, class VECTOR>
PreconditionUseMatrix<MATRIX,VECTOR>::PreconditionUseMatrix(const MATRIX       &M,
                                                            const function_ptr  method)
  :
  matrix(M), precondition(method)
{}



template<class MATRIX, class VECTOR>
void
PreconditionUseMatrix<MATRIX,VECTOR>::vmult (VECTOR &dst,
                                             const VECTOR &src) const
{
  (matrix.*precondition)(dst, src);
}

//---------------------------------------------------------------------------

template<class MATRIX>
inline
PreconditionRelaxation<MATRIX>::AdditionalData::
AdditionalData (const double relaxation)
  :
  relaxation (relaxation)
{}



//---------------------------------------------------------------------------

namespace internal
{
  namespace PreconditionChebyshev
  {
    // for deal.II vectors, perform updates for Chebyshev preconditioner all
    // at once to reduce memory transfer. Here, we select between general
    // vectors and deal.II vectors where we expand the loop over the (local)
    // size of the vector

    // generic part for non-deal.II vectors
    template <typename VECTOR>
    inline
    void
    vector_updates (const VECTOR &src,
                    const VECTOR &matrix_diagonal_inverse,
                    const bool    start_zero,
                    const double  factor1,
                    const double  factor2,
                    VECTOR &update1,
                    VECTOR &update2,
                    VECTOR &dst)
    {
      if (start_zero)
        {
          dst.equ (factor2, src);
          dst.scale (matrix_diagonal_inverse);
          update1.equ(-1.,dst);
        }
      else
        {
          update2 -= src;
          update2.scale (matrix_diagonal_inverse);
          if (factor1 == 0.)
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
      VectorUpdater (const Number *src,
                     const Number *matrix_diagonal_inverse,
                     const bool    start_zero,
                     const Number  factor1,
                     const Number  factor2,
                     Number       *update1,
                     Number       *update2,
                     Number       *dst)
        :
        src (src),
        matrix_diagonal_inverse (matrix_diagonal_inverse),
        do_startup (factor1 == Number()),
        start_zero (start_zero),
        factor1 (factor1),
        factor2 (factor2),
        update1 (update1),
        update2 (update2),
        dst (dst)
      {}

      void
      apply_to_subrange (const std::size_t begin,
                         const std::size_t end) const
      {
        // To circumvent a bug in gcc
        // (https://gcc.gnu.org/bugzilla/show_bug.cgi?id=63945), we create copies
        // of the variables factor1 and factor2 and do not check based on
        // factor1.
        const Number factor1 = this->factor1;
        const Number factor2 = this->factor2;
        if (do_startup)
          {
            if (start_zero)
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i=begin; i<end; ++i)
                {
                  dst[i] = factor2 * src[i] * matrix_diagonal_inverse[i];
                  update1[i] = -dst[i];
                }
            else
              DEAL_II_OPENMP_SIMD_PRAGMA
              for (std::size_t i=begin; i<end; ++i)
                {
                  update1[i] = ((update2[i]-src[i]) *
                                factor2*matrix_diagonal_inverse[i]);
                  dst[i] -= update1[i];
                }
          }
        else
          DEAL_II_OPENMP_SIMD_PRAGMA
          for (std::size_t i=begin; i<end; ++i)
            {
              const Number update =
                factor1 * update1[i] + factor2 *
                ((update2[i] - src[i]) * matrix_diagonal_inverse[i]);
              update1[i] = update;
              dst[i] -= update;
            }
      }

      const Number *src;
      const Number *matrix_diagonal_inverse;
      const bool do_startup;
      const bool start_zero;
      const Number factor1;
      const Number factor2;
      mutable Number *update1;
      mutable Number *update2;
      mutable Number *dst;
    };

    template<typename Number>
    struct VectorUpdatesRange : public parallel::ParallelForInteger
    {
      VectorUpdatesRange(const VectorUpdater<Number> &updater,
                         const std::size_t size)
        :
        updater (updater)
      {
        if (size < internal::Vector::minimum_parallel_grain_size)
          apply_to_subrange (0, size);
        else
          apply_parallel (0, size,
                          internal::Vector::minimum_parallel_grain_size);
      }

      ~VectorUpdatesRange() {}

      virtual void
      apply_to_subrange (const std::size_t begin,
                         const std::size_t end) const
      {
        updater.apply_to_subrange(begin, end);
      }

      const VectorUpdater<Number> &updater;
    };

    // selection for deal.II vector
    template <typename Number>
    inline
    void
    vector_updates (const ::dealii::Vector<Number> &src,
                    const ::dealii::Vector<Number> &matrix_diagonal_inverse,
                    const bool    start_zero,
                    const double  factor1,
                    const double  factor2,
                    ::dealii::Vector<Number> &update1,
                    ::dealii::Vector<Number> &update2,
                    ::dealii::Vector<Number> &dst)
    {
      VectorUpdater<Number> upd(src.begin(), matrix_diagonal_inverse.begin(),
                                start_zero, factor1, factor2,
                                update1.begin(), update2.begin(), dst.begin());
      VectorUpdatesRange<Number>(upd, src.size());
    }

    // selection for parallel deal.II vector
    template <typename Number>
    inline
    void
    vector_updates (const parallel::distributed::Vector<Number> &src,
                    const parallel::distributed::Vector<Number> &matrix_diagonal_inverse,
                    const bool    start_zero,
                    const double  factor1,
                    const double  factor2,
                    parallel::distributed::Vector<Number> &update1,
                    parallel::distributed::Vector<Number> &update2,
                    parallel::distributed::Vector<Number> &dst)
    {
      VectorUpdater<Number> upd(src.begin(), matrix_diagonal_inverse.begin(),
                                start_zero, factor1, factor2,
                                update1.begin(), update2.begin(), dst.begin());
      VectorUpdatesRange<Number>(upd, src.local_size());
    }

    template <typename VECTOR>
    struct DiagonalPreconditioner
    {
      DiagonalPreconditioner (const VECTOR &vector)
        :
        diagonal_vector(vector)
      {}

      void vmult (VECTOR &dst,
                  const VECTOR &src) const
      {
        dst = src;
        dst.scale(diagonal_vector);
      }

      const VECTOR &diagonal_vector;
    };

    struct EigenvalueTracker
    {
    public:
      void slot(const std::vector<double> &eigenvalues)
      {
        values = eigenvalues;
      }

      std::vector<double> values;
    };
  }
}



template <class MATRIX, class VECTOR>
inline
PreconditionChebyshev<MATRIX,VECTOR>::AdditionalData::
AdditionalData (const unsigned int degree,
                const double       smoothing_range,
                const bool         nonzero_starting,
                const unsigned int eig_cg_n_iterations,
                const double       eig_cg_residual,
                const double       max_eigenvalue)
  :
  degree  (degree),
  smoothing_range (smoothing_range),
  nonzero_starting (nonzero_starting),
  eig_cg_n_iterations (eig_cg_n_iterations),
  eig_cg_residual (eig_cg_residual),
  max_eigenvalue (max_eigenvalue)
{}



template <class MATRIX, class VECTOR>
inline
PreconditionChebyshev<MATRIX,VECTOR>::PreconditionChebyshev ()
  :
  is_initialized (false)
{}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::initialize (const MATRIX &matrix,
                                                  const AdditionalData &additional_data)
{
  matrix_ptr = &matrix;
  data = additional_data;
  if (data.matrix_diagonal_inverse.size() != matrix.m())
    {
      Assert(data.matrix_diagonal_inverse.size() == 0,
             ExcMessage("Matrix diagonal vector set but not sized correctly"));
      data.matrix_diagonal_inverse.reinit(matrix.m());
      for (unsigned int i=0; i<matrix.m(); ++i)
        data.matrix_diagonal_inverse(i) = 1./matrix.el(i,i);
    }


  // calculate largest eigenvalue using a hand-tuned CG iteration on the
  // matrix weighted by its diagonal. we start with a vector that consists of
  // ones only, weighted by the length.
  double max_eigenvalue, min_eigenvalue;
  if (data.eig_cg_n_iterations > 0)
    {
      Assert (additional_data.eig_cg_n_iterations > 2,
              ExcMessage ("Need to set at least two iterations to find eigenvalues."));

      // set a very strict tolerance to force at least two iterations
      ReductionControl control (data.eig_cg_n_iterations, 1e-20, 1e-20);
      GrowingVectorMemory<VECTOR> memory;
      VECTOR *rhs = memory.alloc();
      VECTOR *dummy = memory.alloc();
      rhs->reinit(data.matrix_diagonal_inverse, true);
      dummy->reinit(data.matrix_diagonal_inverse);
      *rhs = 1./std::sqrt(static_cast<double>(matrix.m()));

      internal::PreconditionChebyshev::EigenvalueTracker eigenvalue_tracker;
      SolverCG<VECTOR> solver (control, memory);
      solver.connect_eigenvalues_slot(std_cxx11::bind(&internal::PreconditionChebyshev::EigenvalueTracker::slot,
                                                      &eigenvalue_tracker,
                                                      std_cxx11::_1));
      internal::PreconditionChebyshev::DiagonalPreconditioner<VECTOR>
      preconditioner(data.matrix_diagonal_inverse);
      try
        {
          solver.solve(matrix, *dummy, *rhs, preconditioner);
        }
      catch (SolverControl::NoConvergence &)
        {
        }
      Assert(control.last_step() >= 2,
             ExcMessage("Could not find eigenvalues"));

      memory.free(dummy);
      memory.free(rhs);

      // read the eigenvalues from the attached eigenvalue tracker
      if (eigenvalue_tracker.values.empty())
        min_eigenvalue = max_eigenvalue = 1;
      else
        {
          min_eigenvalue = eigenvalue_tracker.values.front();
          max_eigenvalue = eigenvalue_tracker.values.back();
        }

      // include a safety factor since the CG method will in general not be
      // converged
      max_eigenvalue *= 1.2;
    }
  else
    {
      max_eigenvalue = data.max_eigenvalue;
      min_eigenvalue = data.max_eigenvalue/data.smoothing_range;
    }

  const double alpha = (data.smoothing_range > 1. ?
                        max_eigenvalue / data.smoothing_range :
                        std::min(0.9*max_eigenvalue, min_eigenvalue));
  delta = (max_eigenvalue-alpha)*0.5;
  theta = (max_eigenvalue+alpha)*0.5;

  update1.reinit (data.matrix_diagonal_inverse, true);
  update2.reinit (data.matrix_diagonal_inverse, true);

  is_initialized = true;
}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::vmult (VECTOR &dst,
                                             const VECTOR &src) const
{
  Assert (is_initialized, ExcMessage("Preconditioner not initialized"));
  double rhok  = delta / theta,  sigma = theta / delta;
  if (data.nonzero_starting && !dst.all_zero())
    {
      matrix_ptr->vmult (update2, dst);
      internal::PreconditionChebyshev::vector_updates
      (src, data.matrix_diagonal_inverse, false, 0., 1./theta, update1,
       update2, dst);
    }
  else
    internal::PreconditionChebyshev::vector_updates
    (src, data.matrix_diagonal_inverse, true, 0., 1./theta, update1,
     update2, dst);

  for (unsigned int k=0; k<data.degree; ++k)
    {
      matrix_ptr->vmult (update2, dst);
      const double rhokp = 1./(2.*sigma-rhok);
      const double factor1 = rhokp * rhok, factor2 = 2.*rhokp/delta;
      rhok = rhokp;
      internal::PreconditionChebyshev::vector_updates
      (src, data.matrix_diagonal_inverse, false, factor1, factor2, update1,
       update2, dst);
    }
}



template <class MATRIX, class VECTOR>
inline
void
PreconditionChebyshev<MATRIX,VECTOR>::Tvmult (VECTOR &dst,
                                              const VECTOR &src) const
{
  Assert (is_initialized, ExcMessage("Preconditioner not initialized"));
  double rhok  = delta / theta,  sigma = theta / delta;
  if (data.nonzero_starting && !dst.all_zero())
    {
      matrix_ptr->Tvmult (update2, dst);
      internal::PreconditionChebyshev::vector_updates
      (src, data.matrix_diagonal_inverse, false, 0., 1./theta, update1,
       update2, dst);
    }
  else
    internal::PreconditionChebyshev::vector_updates
    (src, data.matrix_diagonal_inverse, true, 0., 1./theta, update1,
     update2, dst);

  for (unsigned int k=0; k<data.degree; ++k)
    {
      matrix_ptr->Tvmult (update2, dst);
      const double rhokp = 1./(2.*sigma-rhok);
      const double factor1 = rhokp * rhok, factor2 = 2.*rhokp/delta;
      rhok = rhokp;
      internal::PreconditionChebyshev::vector_updates
      (src, data.matrix_diagonal_inverse, false, factor1, factor2, update1,
       update2, dst);
    }
}



template <class MATRIX, class VECTOR>
inline
void PreconditionChebyshev<MATRIX,VECTOR>::clear ()
{
  is_initialized = false;
  matrix_ptr = 0;
  data.matrix_diagonal_inverse.reinit(0);
  update1.reinit(0);
  update2.reinit(0);
}



#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
