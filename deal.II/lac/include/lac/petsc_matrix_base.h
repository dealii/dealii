//----------------------------  petsc_matrix_base.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  petsc_matrix_base.h  ---------------------------
#ifndef __deal2__petsc_matrix_base_h
#define __deal2__petsc_matrix_base_h


#include <base/config.h>
#include <base/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <petscmat.h>


namespace PETScWrappers
{
                                   // forward declaration
  class VectorBase;
  
  
/**
 * Base class for all matrix classes that are implemented on top of the PETSc
 * matrix types. Since in PETSc all matrix types (i.e. sequential and
 * parallel, sparse, blocked, etc.)  are built by filling the contents of an
 * abstract object that is only referenced through a pointer of a type that is
 * independent of the actual matrix type, we can implement almost all
 * functionality of matrices in this base class. Derived classes will then only
 * have to provide the functionality to create one or the other kind of
 * matrix.
 *
 * The interface of this class is modeled after the existing
 * @ref{SparseMatrix} class in deal.II. It has almost the same member
 * functions, and is often exchangable. However, since PETSc only supports a
 * single scalar type (either double, float, or a complex data type), it is
 * not templated, and only works with whatever your PETSc installation has
 * defined the data type @p{PetscScalar} to.
 *
 * Note that PETSc only guarantees that operations do what you expect if the
 * functions @p{MatAssemblyBegin} and @p{MatAssemblyEnd} have been called
 * after matrix assembly. Therefore, you need to call
 * @ref{SparseMatrix::compress} before you actually use the matrix. This also
 * calls @p{MatCompress} that compresses the storage format for sparse
 * matrices by discarding unused elements. PETSc allows to continue with
 * assembling the matrix after calls to these functions, but since there are
 * no more free entries available after that any more, it is better to only
 * call @ref{SparseMatrix::compress} once at the end of the assembly stage and
 * before the matrix is actively used.
 * 
 * @author Wolfgang Bangerth, 2004
 */
  class MatrixBase 
  {
    public:
                                       /**
                                        * Default constructor.
                                        */
      MatrixBase ();

                                       /**
                                        * Destructor. Made virtual so that one
                                        * can use pointers to this class.
                                        */
      virtual ~MatrixBase ();
      
                                       /**
                                        * Set all matrix entries to zero, but
                                        * retain the sparsity pattern.
                                        */
      void reinit ();

                                       /**
                                        * Set the element (<i>i,j</i>)
                                        * to <tt>value</tt>. Throws an
                                        * error if the entry does not
                                        * exist. Still, it is allowed to
                                        * store zero values in
                                        * non-existent fields.
                                        */
      void set (const unsigned int i,
                const unsigned int j,
                const PetscScalar value);

                                       /**
                                        * Add <tt>value</tt> to the
                                        * element (<i>i,j</i>).  Throws
                                        * an error if the entry does not
                                        * exist. Still, it is allowed to
                                        * store zero values in
                                        * non-existent fields.
                                        */
      void add (const unsigned int i,
                const unsigned int j,
                const PetscScalar value);

                                       /**
                                        * PETSc matrices store their own
                                        * sparsity patterns. So, in analogy to
                                        * our own @ref{SparsityPattern} class,
                                        * this function compresses the
                                        * sparsity pattern and allows the
                                        * resulting matrix to be used in all
                                        * other operations where before only
                                        * assembly functions were
                                        * allowed. This function must
                                        * therefore be called once you have
                                        * assembled the matrix.
                                        */
      void compress ();
      
                                       /**
                                        * Return the value of the entry
                                        * (<i>i,j</i>).  This may be an
                                        * expensive operation and you should
                                        * always take care where to call this
                                        * function. In contrast to the
                                        * respective function in the
                                        * @p{MatrixBase} class, we don't
                                        * throw an exception if the respective
                                        * entry doesn't exist in the sparsity
                                        * pattern of this class, since PETSc
                                        * does not transmit this information.
                                        *
                                        * This function is therefore exactly
                                        * equivalent to the @p{el()} function.
                                        */
      PetscScalar operator () (const unsigned int i,
                               const unsigned int j) const;

                                       /**
                                        * Return the value of the matrix entry
                                        * (<i>i,j</i>). If this entry does not
                                        * exist in the sparsity pattern, then
                                        * zero is returned. While this may be
                                        * convenient in some cases, note that
                                        * it is simple to write algorithms
                                        * that are slow compared to an optimal
                                        * solution, since the sparsity of the
                                        * matrix is not used.
                                        */
      PetscScalar el (const unsigned int i,
                      const unsigned int j) const;

                                       /**
                                        * Return the number of rows in this
                                        * matrix.
                                        */
      unsigned int m () const;

                                       /**
                                        * Return the number of columns in this
                                        * matrix.
                                        */
      unsigned int n () const;

                                       /**
                                        * Return the number of nonzero
                                        * elements of this
                                        * matrix. Actually, it returns
                                        * the number of entries in the
                                        * sparsity pattern; if any of
                                        * the entries should happen to
                                        * be zero, it is counted anyway.
                                        */
      unsigned int n_nonzero_elements () const;
      
                                       /**
                                        * Return the l1-norm of the matrix, that is
                                        * $|M|_1=max_{all columns j}\sum_{all 
                                        * rows i} |M_ij|$,
                                        * (max. sum of columns).
                                        * This is the
                                        * natural matrix norm that is compatible
                                        * to the l1-norm for vectors, i.e.
                                        * $|Mv|_1\leq |M|_1 |v|_1$.
                                        * (cf. Haemmerlin-Hoffmann:
                                        * Numerische Mathematik)
                                        */
      PetscReal l1_norm () const;

                                       /**
                                        * Return the linfty-norm of the
                                        * matrix, that is
                                        * $|M|_infty=max_{all rows i}\sum_{all 
                                        * columns j} |M_ij|$,
                                        * (max. sum of rows).
                                        * This is the
                                        * natural matrix norm that is compatible
                                        * to the linfty-norm of vectors, i.e.
                                        * $|Mv|_infty \leq |M|_infty |v|_infty$.
                                        * (cf. Haemmerlin-Hoffmann:
                                        * Numerische Mathematik)
                                        */
      PetscReal linfty_norm () const;

                                       /**
                                        * Return the frobenius norm of the
                                        * matrix, i.e. the square root of the
                                        * sum of squares of all entries in the
                                        * matrix.
                                        */
      PetscReal frobenius_norm () const;
      
                                       /**
                                        * Multiply the entire matrix by a
                                        * fixed factor.
                                        */
      MatrixBase & operator *= (const PetscScalar factor);
    
                                       /**
                                        * Divide the entire matrix by a
                                        * fixed factor.
                                        */
      MatrixBase & operator /= (const PetscScalar factor);

                                       /**
                                        * Matrix-vector multiplication:
                                        * let <i>dst = M*src</i> with
                                        * <i>M</i> being this matrix.
                                        *
                                        * Source and destination must
                                        * not be the same vector.
                                        */      
      void vmult (VectorBase       &dst,
                  const VectorBase &src) const;

                                       /**
                                        * Matrix-vector multiplication: let
                                        * <i>dst = M<sup>T</sup>*src</i> with
                                        * <i>M</i> being this matrix. This
                                        * function does the same as vmult()
                                        * but takes the transposed matrix.
                                        *
                                        * Source and destination must
                                        * not be the same vector.
                                        */
      void Tvmult (VectorBase       &dst,
                   const VectorBase &src) const;

                                       /**
                                        * Adding Matrix-vector
                                        * multiplication. Add
                                        * <i>M*src</i> on <i>dst</i>
                                        * with <i>M</i> being this
                                        * matrix.
                                        *
                                        * Source and destination must
                                        * not be the same vector.
                                        */
      void vmult_add (VectorBase       &dst,
                      const VectorBase &src) const;

                                       /**
                                        * Adding Matrix-vector
                                        * multiplication. Add
                                        * <i>M<sup>T</sup>*src</i> to
                                        * <i>dst</i> with <i>M</i> being
                                        * this matrix. This function
                                        * does the same as vmult_add()
                                        * but takes the transposed
                                        * matrix.
                                        *
                                        * Source and destination must
                                        * not be the same vector.
                                        */
      void Tvmult_add (VectorBase       &dst,
                       const VectorBase &src) const;

                                       /**
                                        * Return the square of the norm
                                        * of the vector $v$ with respect
                                        * to the norm induced by this
                                        * matrix,
                                        * i.e. $\left(v,Mv\right)$. This
                                        * is useful, e.g. in the finite
                                        * element context, where the
                                        * $L_2$ norm of a function
                                        * equals the matrix norm with
                                        * respect to the mass matrix of
                                        * the vector representing the
                                        * nodal values of the finite
                                        * element function.
                                        *
                                        * Obviously, the matrix needs to
                                        * be quadratic for this operation.
                                        *
                                        * The implementation of this function
                                        * is not as efficient as the one in
                                        * the @p{MatrixBase} class used in
                                        * deal.II (i.e. the original one, not
                                        * the PETSc wrapper class) since PETSc
                                        * doesn't support this operation and
                                        * needs a temporary vector.
                                        */
      PetscScalar matrix_norm_square (const VectorBase &v) const;

                                       /**
                                        * Compute the matrix scalar
                                        * product $\left(u,Mv\right)$.
                                        *
                                        * The implementation of this function
                                        * is not as efficient as the one in
                                        * the @p{MatrixBase} class used in
                                        * deal.II (i.e. the original one, not
                                        * the PETSc wrapper class) since PETSc
                                        * doesn't support this operation and
                                        * needs a temporary vector.
                                        */
      PetscScalar matrix_scalar_product (const VectorBase &u,
                                         const VectorBase &v) const;

                                       /**
                                        * Compute the residual of an
                                        * equation <i>Mx=b</i>, where
                                        * the residual is defined to be
                                        * <i>r=b-Mx</i>. Write the
                                        * residual into
                                        * <tt>dst</tt>. The
                                        * <i>l<sub>2</sub></i> norm of
                                        * the residual vector is
                                        * returned.
                                        *
                                        * Source <i>x</i> and destination
                                        * <i>dst</i> must not be the same
                                        * vector.
                                        */
      PetscScalar residual (VectorBase       &dst,
                            const VectorBase &x,
                            const VectorBase &b) const;

                                       /**
                                        * Conversion operator to gain access
                                        * to the underlying PETSc type. If you
                                        * do this, you cut this class off some
                                        * information it may need, so this
                                        * conversion operator should only be
                                        * used if you know what you do. In
                                        * particular, it should only be used
                                        * for read-only operations into the
                                        * matrix.
                                        */
      operator const Mat () const;

                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcPETScError,
                      int,
                      << "An error with error number " << arg1
                      << " occured while calling a PETSc function");
                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcSourceEqualsDestination);
      
    protected:
                                       /**
                                        * A generic matrix object in
                                        * PETSc. The actual type, a sparse
                                        * matrix, is set in the constructor.
                                        */
      Mat matrix;

                                       /**
                                        * PETSc doesn't allow to mix additions
                                        * to matrix entries and overwriting
                                        * them (to make synchronisation of
                                        * parallel computations
                                        * simpler). Since the interface of the
                                        * existing classes don't support the
                                        * notion of not interleaving things,
                                        * we have to emulate this
                                        * ourselves. The way we do it is to,
                                        * for each access operation, store
                                        * whether it is an insertion or an
                                        * addition. If the previous one was of
                                        * different type, then we first have
                                        * to flush the PETSc buffers;
                                        * otherwise, we can simply go on.
                                        *
                                        * The following structure and variable
                                        * declare and store the previous
                                        * state.
                                        */
      struct LastAction
      {
          enum Values { none, insert, add };
      };

                                       /**
                                        * Store whether the last action was a
                                        * write or add operation.
                                        */
      LastAction::Values last_action;            
  };




// -------------------------- inline and template functions ----------------------

  inline
  PetscScalar
  MatrixBase::operator() (const unsigned int i,
                          const unsigned int j) const
  {
    return el(i,j);
  }
  
  
}

#endif // DEAL_II_USE_PETSC


/*----------------------------   petsc_matrix_base.h     ---------------------------*/

#endif
/*----------------------------   petsc_matrix_base.h     ---------------------------*/
