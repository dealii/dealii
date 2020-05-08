// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
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

#ifndef dealii_precondition_selector_h
#define dealii_precondition_selector_h


#include <deal.II/base/config.h>

#include <deal.II/base/smartpointer.h>
#include <deal.II/base/subscriptor.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <class number>
class Vector;
template <class number>
class SparseMatrix;
#endif


/*! @addtogroup Preconditioners
 *@{
 */

/**
 * Selects the preconditioner. The constructor of this class takes the name of
 * the preconditioning and the damping parameter @p omega of the
 * preconditioning and the @p use_matrix function takes the matrix that is
 * used by the matrix-builtin precondition functions. Each time, the
 * <tt>operator()</tt> function is called, this preselected preconditioner,
 * this matrix and this @p omega is used for the preconditioning. This class
 * is designed for being used as argument of the @p solve function of a @p
 * Solver and it covers the selection of all matrix-builtin precondition
 * functions. The selection of other preconditioners, like BlockSOR or ILU
 * should be handled in derived classes by the user.
 *
 * <h3>Usage</h3> The simplest use of this class is the following:
 * @code
 * // generate a @p SolverControl and a @p VectorMemory
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *
 * // generate a solver
 * SolverCG<SparseMatrix<double>, Vector<double> > solver(control, memory);
 *
 * // generate a @p PreconditionSelector
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *
 * // give a matrix whose diagonal entries are to be used for the
 * // preconditioning. Generally the matrix of the linear equation system Ax=b.
 * preconditioning.use_matrix(A);
 *
 * // call the @p solve function with this preconditioning as last argument
 * solver.solve(A,x,b,preconditioning);
 * @endcode
 * The same example where also the @p SolverSelector class is used reads
 * @code
 * // generate a @p SolverControl and a @p VectorMemory
 * SolverControl control;
 * VectorMemory<Vector<double> > memory;
 *
 * // generate a @p SolverSelector that calls the @p SolverCG
 * SolverSelector<SparseMatrix<double>, Vector<double> >
 *   solver_selector("cg", control, memory);
 *
 * // generate a @p PreconditionSelector
 * PreconditionSelector<SparseMatrix<double>, Vector<double> >
 *   preconditioning("jacobi", 1.);
 *
 * preconditioning.use_matrix(A);
 *
 * solver_selector.solve(A,x,b,preconditioning);
 * @endcode
 * Now the use of the @p SolverSelector in combination with the @p
 * PreconditionSelector allows the user to select both, the solver and the
 * preconditioner, at the beginning of their program and each time the solver is
 * started (that is several times e.g. in a nonlinear iteration) this
 * preselected solver and preconditioner is called.
 *
 * @author Ralf Hartmann, 1999; extension for full compatibility with
 * LinearOperator class: Jean-Paul Pelteret, 2015
 */
template <typename MatrixType = SparseMatrix<double>,
          typename VectorType = dealii::Vector<double>>
class PreconditionSelector : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  using size_type = typename MatrixType::size_type;

  /**
   * Constructor. @p omega denotes the damping parameter of the
   * preconditioning.
   */
  PreconditionSelector(const std::string &                    preconditioning,
                       const typename VectorType::value_type &omega = 1.);

  /**
   * Destructor.
   */
  virtual ~PreconditionSelector() override;

  /**
   * Takes the matrix that is needed for preconditionings that involves a
   * matrix. e.g. for @p precondition_jacobi, <tt>~_sor</tt>, <tt>~_ssor</tt>.
   */
  void
  use_matrix(const MatrixType &M);

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

  /**
   * Precondition procedure. Calls the preconditioning that was specified in
   * the constructor.
   */
  virtual void
  vmult(VectorType &dst, const VectorType &src) const;

  /**
   * Transpose precondition procedure. Calls the preconditioning that was
   * specified in the constructor.
   */
  virtual void
  Tvmult(VectorType &dst, const VectorType &src) const;

  /**
   * Get the names of all implemented preconditionings. The list of possible
   * options includes:
   * <ul>
   * <li>  "none" </li>
   * <li>  "jacobi" </li>
   * <li>  "sor" </li>
   * <li>  "ssor" </li>
   * </ul>
   */
  static std::string
  get_precondition_names();

  /**
   * @addtogroup Exceptions
   * @{
   */


  /**
   * Exception.
   */
  DeclException0(ExcNoMatrixGivenToUse);

  //@}
protected:
  /**
   * Stores the name of the preconditioning.
   */
  std::string preconditioning;

private:
  /**
   * Matrix that is used for the matrix-builtin preconditioning function. cf.
   * also @p PreconditionUseMatrix.
   */
  SmartPointer<const MatrixType, PreconditionSelector<MatrixType, VectorType>>
    A;

  /**
   * Stores the damping parameter of the preconditioner.
   */
  const typename VectorType::value_type omega;
};

/*@}*/
/* --------------------- Inline and template functions ------------------- */


template <typename MatrixType, typename VectorType>
PreconditionSelector<MatrixType, VectorType>::PreconditionSelector(
  const std::string &                    preconditioning,
  const typename VectorType::value_type &omega)
  : preconditioning(preconditioning)
  , omega(omega)
{}


template <typename MatrixType, typename VectorType>
PreconditionSelector<MatrixType, VectorType>::~PreconditionSelector()
{
  // release the matrix A
  A = nullptr;
}


template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::use_matrix(const MatrixType &M)
{
  A = &M;
}


template <typename MatrixType, typename VectorType>
inline typename PreconditionSelector<MatrixType, VectorType>::size_type
PreconditionSelector<MatrixType, VectorType>::m() const
{
  Assert(A != nullptr, ExcNoMatrixGivenToUse());
  return A->m();
}


template <typename MatrixType, typename VectorType>
inline typename PreconditionSelector<MatrixType, VectorType>::size_type
PreconditionSelector<MatrixType, VectorType>::n() const
{
  Assert(A != nullptr, ExcNoMatrixGivenToUse());
  return A->n();
}



template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::vmult(VectorType &      dst,
                                                    const VectorType &src) const
{
  if (preconditioning == "none")
    {
      dst = src;
    }
  else
    {
      Assert(A != nullptr, ExcNoMatrixGivenToUse());

      if (preconditioning == "jacobi")
        {
          A->precondition_Jacobi(dst, src, omega);
        }
      else if (preconditioning == "sor")
        {
          A->precondition_SOR(dst, src, omega);
        }
      else if (preconditioning == "ssor")
        {
          A->precondition_SSOR(dst, src, omega);
        }
      else
        Assert(false, ExcNotImplemented());
    }
}


template <typename MatrixType, typename VectorType>
void
PreconditionSelector<MatrixType, VectorType>::Tvmult(
  VectorType &      dst,
  const VectorType &src) const
{
  if (preconditioning == "none")
    {
      dst = src;
    }
  else
    {
      Assert(A != nullptr, ExcNoMatrixGivenToUse());

      if (preconditioning == "jacobi")
        {
          A->precondition_Jacobi(dst, src, omega); // Symmetric operation
        }
      else if (preconditioning == "sor")
        {
          A->precondition_TSOR(dst, src, omega);
        }
      else if (preconditioning == "ssor")
        {
          A->precondition_SSOR(dst, src, omega); // Symmetric operation
        }
      else
        Assert(false, ExcNotImplemented());
    }
}


template <typename MatrixType, typename VectorType>
std::string
PreconditionSelector<MatrixType, VectorType>::get_precondition_names()
{
  return "none|jacobi|sor|ssor";
}


DEAL_II_NAMESPACE_CLOSE

#endif
