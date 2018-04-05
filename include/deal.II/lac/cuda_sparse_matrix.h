// ---------------------------------------------------------------------
//
// Copyright (C) 2018 by the deal.II authors
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

#ifndef dealii_cuda_sparse_matrix_h
#define dealii_cuda_sparse_matrix_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>

#ifdef DEAL_II_WITH_CUDA
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/cuda_vector.h>
#include <cusparse.h>

DEAL_II_NAMESPACE_OPEN

namespace CUDAWrappers
{
  /**
   * This class is a wrapper around cuSPARSE csr sparse matrix. Unlike deal.II's
   * own SparseMatrix all elements within each row are stored in increasing
   * column index order.
   *
   * @note Instantiations for this template are provided for <tt>@<float@> and
   * @<double@></tt>.
   *
   * @ingroup Matrix1
   * @author Bruno Turcksin
   * @date 2018
   */
  template <typename Number>
  class SparseMatrix: public virtual Subscriptor
  {
  public:
    /**
     *  Declare type for container size.
     */
    typedef unsigned int size_type;

    /**
     * Type of the matrix entries.
     */
    typedef Number value_type;

    /**
     * Declare a type that holds real-valued numbers with the same precision
     * as the template argument to this class.
     */
    typedef Number real_type;

    /**
     * @name Constructors and initialization
     */
    //@{
    /**
     * Constructor. Initialize the matrix to be empty, without any structure,
     * i.e., the matrix is not usable at all. This constructor is therefore
     * only useful for matrices which are members of a class.
     *
     * You have to initialize the matrix before usage with reinit.
     */
    SparseMatrix();

    /**
     * Constructor. Takes a cuSPARSE handle and a sparse matrix on the host.
     * The sparse matrix on the host is copied on the device and the elements
     * are reordered according to the format supported by cuSPARSE.
     */
    SparseMatrix(cusparseHandle_t handle,
                 const ::dealii::SparseMatrix<Number> &sparse_matrix_host);

    /**
     * Move constructor. Create a new SparseMatrix by stealing the internal
     * data.
     */
    SparseMatrix(CUDAWrappers::SparseMatrix<Number> &&);

    /**
     * Copy constructor is deleted.
     */
    SparseMatrix(const CUDAWrappers::SparseMatrix<Number> &) = delete;

    /**
     * Destructor. Free all memory.
     */
    ~SparseMatrix();

    /**
     * Reinitialize the sparse matrix. The sparse matrix on the host is copied
     * to the device and the elementes are reordered according to the format
     * supported by cuSPARSE.
     */
    void reinit(cusparseHandle_t handle,
                const ::dealii::SparseMatrix<Number> &sparse_matrix_host);
    //@}

    /**
     * @name Information on the matrix
     */
    //@{
    /**
     * Return the dimension of the codomain (or range) space. Note that the
     * matrix is of dimension $m \times n$.
     */
    size_type m() const;

    /**
     * Return the dimension of the domain space. Note that the matrix is of
     * dimension $m \times n$.
     */
    size_type n() const;

    /**
     * Return the number of nonzero elements of this matrix. Actually, it
     * returns the number of entries in the sparsity pattern; if any of the
     * entries should happen to be zero, it is counted anyway.
     */
    std::size_t n_nonzero_elements() const;
    //@}

    /**
     * @name Modifying entries
     */
    //@{
    /**
     * Multiply the entire matrix by a fixed factor.
     */
    SparseMatrix &operator*= (const Number factor);

    /**
     * Divide the entrie matrix by a fixed factor.
     */
    SparseMatrix &operator/= (const Number factor);
    //@}

    /**
     * @name Multiplications
     */
    //@{
    /**
     * Matrix-vector multiplication: let $dst = M \cdot src$ with $M$
     * being this matrix.
     */
    void vmult(LinearAlgebra::CUDAWrappers::Vector<Number> &dst,
               const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Matrix-vector multiplication: let $dst = M^T \cdot src$ with
     * $M$ being this matrix. This function does the same as vmult() but
     * takes this transposed matrix.
     */
    void Tvmult(LinearAlgebra::CUDAWrappers::Vector<Number> &dst,
                const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Adding matrix-vector multiplication. Add $M \cdot src$ on $dst$
     * with $M$ being this matrix.
     */
    void vmult_add(LinearAlgebra::CUDAWrappers::Vector<Number> &dst,
                   const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Adding matrix-vector multiplication. Add $M^T \cdot src$ to
     * $dst$ with $M$ being this matrix. This function foes the same
     * as vmult_add() but takes the transposed matrix.
     */
    void Tvmult_add(LinearAlgebra::CUDAWrappers::Vector<Number> &dst,
                    const LinearAlgebra::CUDAWrappers::Vector<Number> &src) const;

    /**
     * Return the square of the norm of the vector $v$ with respect to the
     * norm induced by this matrix, i.e., $\left(v,Mv\right)$. This is useful,
     * e.g., in the finite context, where the $L_2$ norm of a function equals
     * the matrix norm with respect to the mass matrix of the vector
     * representing the nodal values of the finite element function.
     *
     * Obviously, the matrix needs to be quadratic for this operation.
     */
    Number matrix_norm_square(const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const;

    /**
     * Compute the matrix scalar product $\left(u,Mv\right)$.
     */
    Number matrix_scalar_product(const LinearAlgebra::CUDAWrappers::Vector<Number> &u,
                                 const LinearAlgebra::CUDAWrappers::Vector<Number> &v) const;

    /**
     * Compute the residual of an equation $M \cdot x=b$, where the residual is
     * defined to be $r=b-M \cdot x$. Write the residual into $dst$. The
     * $l_2$ norm of the residual vector is returned.
     *
     * Source $x$ and destination $dst$ must not be the same vector.
     */
    Number residual(LinearAlgebra::CUDAWrappers::Vector<Number>       &dst,
                    const LinearAlgebra::CUDAWrappers::Vector<Number> &x,
                    const LinearAlgebra::CUDAWrappers::Vector<Number> &b) const;
    //@}

    /**
     * @name Matrix norms
     */
    //@{
    /**
     * Return the $l_1$-norm of the matrix, that is $|M|_1=\max_{\mathrm{all\
     * columns\ }j}\sum_{\mathrm{all\ rows\ }i} |M_{ij}|$, (max. sum of
     * columns). This is the natural matrix norm that is compatible to the
     * $l_1$-norm for vectors, i.e., $|Mv|_1\leq |M|_1 |v|_1$.
     */
    Number l1_norm() const;

    /**
     * Return the $l_\infty$-norm of the matrix, that is
     * $|M|_\infty=\max_{\mathrm{all\ rows\ }i}\sum_{\mathrm{all\ columns\ }j}
     * |M_{ij}|$, (max. sum of rows). This is the natural norm that is
     * compatible to the $l_\infty$-norm of vectors, i.e., $|Mv|_\infty \leq
     * |M|_\infty |v|_\infty$.
     */
    Number linfty_norm() const;

    /**
     * Return the frobenius norm of the matrix, i.e., the square root of the
     * sum of squares of all entries in the matrix.
     */
    Number frobenius_norm() const;
    //@}

    /**
     *@name Access to underlying CUDA data
     */
    //@{
    /**
     * Return a tuple containing the pointer to the values of matrix, the
     * pointer to the columns indices, the pointer to the rows pointer, and
     * the cuSPARSE matrix description.
     */
    std::tuple<Number *, int *, int *, cusparseMatDescr_t>
    get_cusparse_matrix();
    //*}

  private:
    /**
     * cuSPARSE used to call cuSPARSE function. The cuSPARSE handle needs to
     * be mutable to be called in a const function.
     */
    mutable cusparseHandle_t cusparse_handle;

    /**
     * Number of non-zero elements in the sparse matrix.
     */
    int nnz;

    /**
     * Number of rows of the sparse matrix.
     */
    int n_rows;

    /**
     * Number of columns of the sparse matrix.
     */
    int n_cols;

    /**
     * Pointer to the values (on the device) of the sparse matrix.
     */
    Number *val_dev;

    /**
     * Pointer to the column indices (on the device) of the sparse matrix.
     */
    int *column_index_dev;

    /**
     * Pointer to the row pointer (on the device) of the sparse matrix.
     */
    int *row_ptr_dev;

    /**
     * cuSPARSE description of the sparse matrix.
     */
    cusparseMatDescr_t descr;
  };



  template <typename Number>
  inline
  unsigned int SparseMatrix<Number>::m() const
  {
    return n_rows;
  }



  template <typename Number>
  inline
  unsigned int SparseMatrix<Number>::n() const
  {
    return n_cols;
  }



  template <typename Number>
  inline
  std::size_t SparseMatrix<Number>::n_nonzero_elements() const
  {
    return nnz;
  }
}

DEAL_II_NAMESPACE_CLOSE

#endif
#endif
