//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__petsc_matrix_base_h
#define __deal2__petsc_matrix_base_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>

#ifdef DEAL_II_USE_PETSC

#include <petscmat.h>
#include <boost/shared_ptr.hpp>
#include <vector>

namespace PETScWrappers
{
                                   // forward declarations
  class VectorBase;
  class MatrixBase;

  namespace MatrixIterators
  {
/**
 * STL conforming iterator. This class acts as an iterator walking over the
 * elements of PETSc matrices. Since PETSc offers a uniform interface for all
 * types of matrices, this iterator can be used to access both sparse and full
 * matrices.
 *
 * Note that PETSc does not give any guarantees as to the order of elements
 * within each row. Note also that accessing the elements of a full matrix
 * surprisingly only shows the nonzero elements of the matrix, not all
 * elements.
 *
 * @ingroup PETScWrappers
 * @author Guido Kanschat, Roy Stogner, Wolfgang Bangerth, 2004
 */
    class const_iterator
    {
      private:
                                         /**
                                          * Accessor class for iterators
                                          */
        class Accessor
        {
          public:
                                             /**
                                              * Constructor. Since we use
                                              * accessors only for read
                                              * access, a const matrix
                                              * pointer is sufficient.
                                              */
            Accessor (const MatrixBase    *matrix,
                      const unsigned int   row,
                      const unsigned int   index);

                                             /**
                                              * Row number of the element
                                              * represented by this
                                              * object.
                                              */
            unsigned int row() const;

                                             /**
                                              * Index in row of the element
                                              * represented by this
                                              * object.
                                              */
            unsigned int index() const;

                                             /**
                                              * Column number of the
                                              * element represented by
                                              * this object.
                                              */
            unsigned int column() const;

                                             /**
                                              * Value of this matrix entry.
                                              */
            PetscScalar value() const;

                                             /**
                                              * Exception
                                              */
            DeclException0 (ExcBeyondEndOfMatrix);
                                             /**
                                              * Exception
                                              */
            DeclException3 (ExcAccessToNonlocalRow,
                            int, int, int,
                            << "You tried to access row " << arg1
                            << " of a distributed matrix, but only rows "
                            << arg2 << " through " << arg3
                            << " are stored locally and can be accessed.");
            
          private:
                                             /**
                                              * The matrix accessed.
                                              */
            mutable MatrixBase *matrix;

                                             /**
                                              * Current row number.
                                              */
            unsigned int a_row;

                                             /**
                                              * Current index in row.
                                              */
            unsigned int a_index;

                                             /**
                                              * Cache where we store the
                                              * column indices of the present
                                              * row. This is necessary, since
                                              * PETSc makes access to the
                                              * elements of its matrices
                                              * rather hard, and it is much
                                              * more efficient to copy all
                                              * column entries of a row once
                                              * when we enter it than
                                              * repeatedly asking PETSc for
                                              * individual ones. This also
                                              * makes some sense since it is
                                              * likely that we will access
                                              * them sequentially anyway.
                                              *
                                              * In order to make copying of
                                              * iterators/accessor of
                                              * acceptable performance, we
                                              * keep a shared pointer to these
                                              * entries so that more than one
                                              * accessor can access this data
                                              * if necessary.
                                              */
            boost::shared_ptr<const std::vector<unsigned int> > colnum_cache;

                                             /**
                                              * Similar cache for the values
                                              * of this row.
                                              */
            boost::shared_ptr<const std::vector<PetscScalar> > value_cache;
            
                                             /**
                                              * Discard the old row caches
                                              * (they may still be used by
                                              * other accessors) and generate
                                              * new ones for the row pointed
                                              * to presently by this accessor.
                                              */
            void visit_present_row ();

                                             /**
                                              * Make enclosing class a
                                              * friend.
                                              */
            friend class const_iterator;
        };
        
      public:
          
                                         /**
                                          * Constructor. Create an iterator
                                          * into the matrix @p matrix for the
                                          * given row and the index within it.
                                          */ 
        const_iterator (const MatrixBase   *matrix,
                        const unsigned int  row,
                        const unsigned int  index);
          
                                         /**
                                          * Prefix increment.
                                          */
        const_iterator& operator++ ();

                                         /**
                                          * Postfix increment.
                                          */
        const_iterator operator++ (int);

                                         /**
                                          * Dereferencing operator.
                                          */
        const Accessor& operator* () const;

                                         /**
                                          * Dereferencing operator.
                                          */
        const Accessor* operator-> () const;

                                         /**
                                          * Comparison. True, if
                                          * both iterators point to
                                          * the same matrix
                                          * position.
                                          */
        bool operator == (const const_iterator&) const;
                                         /**
                                          * Inverse of <tt>==</tt>.
                                          */
        bool operator != (const const_iterator&) const;

                                         /**
                                          * Comparison
                                          * operator. Result is true
                                          * if either the first row
                                          * number is smaller or if
                                          * the row numbers are
                                          * equal and the first
                                          * index is smaller.
                                          */
        bool operator < (const const_iterator&) const;

                                         /**
                                          * Exception
                                          */
        DeclException2 (ExcInvalidIndexWithinRow,
                        int, int,
                        << "Attempt to access element " << arg2
                        << " of row " << arg1
                        << " which doesn't have that many elements.");
        
      private:
                                         /**
                                          * Store an object of the
                                          * accessor class.
                                          */
        Accessor accessor;
    };
    
  }
  
  
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
 * SparseMatrix class in deal.II. It has almost the same member
 * functions, and is often exchangable. However, since PETSc only supports a
 * single scalar type (either double, float, or a complex data type), it is
 * not templated, and only works with whatever your PETSc installation has
 * defined the data type @p PetscScalar to.
 *
 * Note that PETSc only guarantees that operations do what you expect if the
 * functions @p MatAssemblyBegin and @p MatAssemblyEnd have been called
 * after matrix assembly. Therefore, you need to call
 * SparseMatrix::compress() before you actually use the matrix. This also
 * calls @p MatCompress that compresses the storage format for sparse
 * matrices by discarding unused elements. PETSc allows to continue with
 * assembling the matrix after calls to these functions, but since there are
 * no more free entries available after that any more, it is better to only
 * call SparseMatrix::compress() once at the end of the assembly stage and
 * before the matrix is actively used.
 * 
 * @ingroup PETScWrappers
 * @author Wolfgang Bangerth, 2004
 */
  class MatrixBase : public Subscriptor
  {
    public:
                                       /**
                                        * Declare a typedef for the iterator
                                        * class.
                                        */
      typedef MatrixIterators::const_iterator const_iterator;

                                       /**
                                        * Declare a typedef in analogy to all
                                        * the other container classes.
                                        */
      typedef PetscScalar value_type;
      
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
                                        * This operator assigns a scalar to a
                                        * matrix. Since this does usually not
                                        * make much sense (should we set all
                                        * matrix entries to this value? Only
                                        * the nonzero entries of the sparsity
                                        * pattern?), this operation is only
                                        * allowed if the actual value to be
                                        * assigned is zero. This operator only
                                        * exists to allow for the obvious
                                        * notation <tt>matrix=0</tt>, which
                                        * sets all elements of the matrix to
                                        * zero, but keeps the sparsity pattern
                                        * previously used.
                                        */
      MatrixBase &
      operator = (const double d);      
                                       /**
                                        * Release all memory and return
                                        * to a state just like after
                                        * having called the default
                                        * constructor.
                                        */
      void clear ();

                                       /**
                                        * Set the element (<i>i,j</i>)
                                        * to @p value.
					*
					* If the present object (from
					* a derived class of this one)
					* happens to be a sparse
					* matrix, then this function
					* adds a new entry to the
					* matrix if it didn't exist
					* before, very much in
					* contrast to the SparseMatrix
					* class which throws an error
					* if the entry does not exist.
                                        */
      void set (const unsigned int i,
                const unsigned int j,
                const PetscScalar value);

                                       /**
                                        * Add @p value to the
                                        * element (<i>i,j</i>).
					*
					* If the present object (from
					* a derived class of this one)
					* happens to be a sparse
					* matrix, then this function
					* adds a new entry to the
					* matrix if it didn't exist
					* before, very much in
					* contrast to the SparseMatrix
					* class which throws an error
					* if the entry does not exist.
                                        */
      void add (const unsigned int i,
                const unsigned int j,
                const PetscScalar value);

                                       /**
                                        * PETSc matrices store their own
                                        * sparsity patterns. So, in analogy to
                                        * our own SparsityPattern class,
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
                                        * @p MatrixBase class, we don't
                                        * throw an exception if the respective
                                        * entry doesn't exist in the sparsity
                                        * pattern of this class, since PETSc
                                        * does not transmit this information.
                                        *
                                        * This function is therefore exactly
                                        * equivalent to the <tt>el()</tt> function.
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
                                        * Return the main diagonal
                                        * element in the <i>i</i>th
                                        * row. This function throws an
                                        * error if the matrix is not
                                        * quadratic.
                                        *
                                        * Since we do not have direct access
                                        * to the underlying data structure,
                                        * this function is no faster than the
                                        * elementwise access using the el()
                                        * function. However, we provide this
                                        * function for compatibility with the
                                        * SparseMatrix class.
                                        */
      PetscScalar diag_element (const unsigned int i) const;
      
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
                                        * Return the local dimension of the
                                        * matrix, i.e. the number of rows
                                        * stored on the present MPI
                                        * process. For sequential matrices,
                                        * this number is the same as m(),
                                        * but for parallel matrices it may be
                                        * smaller.
					*
					* To figure out which elements
					* exactly are stored locally,
					* use local_range().
                                        */
      unsigned int local_size () const;

                                       /**
					* Return a pair of indices
					* indicating which rows of
					* this matrix are stored
					* locally. The first number is
					* the index of the first
					* row stored, the second
					* the index of the one past
					* the last one that is stored
					* locally. If this is a
					* sequential matrix, then the
					* result will be the pair
					* (0,m()), otherwise it will be
					* a pair (i,i+n), where
					* <tt>n=local_size()</tt>.
					*/
      std::pair<unsigned int, unsigned int>
      local_range () const;

				       /**
					* Return whether @p index is
					* in the local range or not,
					* see also local_range().
					*/
      bool in_local_range (const unsigned int index) const;

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
                                        * Number of entries in a specific row.
                                        */
      unsigned int row_length (const unsigned int row) const;
      
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
                                        * the @p MatrixBase class used in
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
                                        * the @p MatrixBase class used in
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
                                        * @p dst. The
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
                                        * STL-like iterator with the
                                        * first entry.
                                        */
      const_iterator begin () const;

                                       /**
                                        * Final iterator.
                                        */
      const_iterator end () const;

                                       /**
                                        * STL-like iterator with the
                                        * first entry of row @p r.
                                        *
                                        * Note that if the given row is empty,
                                        * i.e. does not contain any nonzero
                                        * entries, then the iterator returned by
                                        * this function equals
                                        * <tt>end(r)</tt>. Note also that the
                                        * iterator may not be dereferencable in
                                        * that case.
                                        */
      const_iterator begin (const unsigned int r) const;

                                       /**
                                        * Final iterator of row <tt>r</tt>. It
                                        * points to the first element past the
                                        * end of line @p r, or past the end of
                                        * the entire sparsity pattern.
                                        *
                                        * Note that the end iterator is not
                                        * necessarily dereferencable. This is in
                                        * particular the case if it is the end
                                        * iterator for the last row of a matrix.
                                        */
      const_iterator end (const unsigned int r) const;
      
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



/// @if NoDoc
// -------------------------- inline and template functions ----------------------


  namespace MatrixIterators
  {

    inline
    const_iterator::Accessor::
    Accessor (const MatrixBase   *matrix,
              const unsigned int  row,
              const unsigned int  index)
                    :
                    matrix(const_cast<MatrixBase*>(matrix)),
                    a_row(row),
                    a_index(index)
    {
      visit_present_row ();
    }


    inline
    unsigned int
    const_iterator::Accessor::row() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_row;
    }


    inline
    unsigned int
    const_iterator::Accessor::column() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*colnum_cache)[a_index];
    }


    inline
    unsigned int
    const_iterator::Accessor::index() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return a_index;
    }


    inline
    PetscScalar
    const_iterator::Accessor::value() const
    {
      Assert (a_row < matrix->m(), ExcBeyondEndOfMatrix());
      return (*value_cache)[a_index];
    }


    inline
    const_iterator::
    const_iterator(const MatrixBase   *matrix,
                   const unsigned int  row,
                   const unsigned int  index)
                    :
                    accessor(matrix, row, index)
    {}



    inline
    const_iterator &
    const_iterator::operator++ ()
    {
      Assert (accessor.a_row < accessor.matrix->m(), ExcIteratorPastEnd());

      ++accessor.a_index;

                                       // if at end of line: do one step, then
                                       // cycle until we find a row with a
                                       // nonzero number of entries
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;
      
          while (accessor.a_index >= accessor.matrix->row_length(accessor.a_row))
            {
              ++accessor.a_row;

                                           // if we happened to find the end
                                           // of the matrix, then stop here
              if (accessor.a_row == accessor.matrix->m())
                break;
            }

          accessor.visit_present_row();
        }
      return *this;
    }


    inline
    const_iterator
    const_iterator::operator++ (int)
    {
      const const_iterator old_state = *this;
      ++(*this);
      return old_state;
    }


    inline
    const const_iterator::Accessor &
    const_iterator::operator* () const
    {
      return accessor;
    }


    inline
    const const_iterator::Accessor *
    const_iterator::operator-> () const
    {
      return &accessor;
    }


    inline
    bool
    const_iterator::
    operator == (const const_iterator& other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }


    inline
    bool
    const_iterator::
    operator != (const const_iterator& other) const
    {
      return ! (*this == other);
    }


    inline
    bool
    const_iterator::
    operator < (const const_iterator& other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }
    
  }
  
  
  inline
  PetscScalar
  MatrixBase::operator() (const unsigned int i,
                          const unsigned int j) const
  {
    return el(i,j);
  }

  

  inline
  MatrixBase::const_iterator
  MatrixBase::begin() const
  {
    return const_iterator(this, 0, 0);
  }


  inline
  MatrixBase::const_iterator
  MatrixBase::end() const
  {
    return const_iterator(this, m(), 0);
  }


  inline
  MatrixBase::const_iterator
  MatrixBase::begin(const unsigned int r) const
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));
    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end (r);
  }


  inline
  MatrixBase::const_iterator
  MatrixBase::end(const unsigned int r) const
  {
    Assert (r < m(), ExcIndexRange(r, 0, m()));

                                     // place the iterator on the first entry
                                     // past this line, or at the end of the
                                     // matrix
    for (unsigned int i=r+1; i<m(); ++i)
      if (row_length(i) > 0)
        return const_iterator(this, i, 0);
    
                                     // if there is no such line, then take the
                                     // end iterator of the matrix
    return end();
  }


  
  inline
  bool
  MatrixBase::in_local_range (const unsigned int index) const
  {
    int begin, end;
    const int ierr = MatGetOwnershipRange (static_cast<const Mat &>(matrix),
					   &begin, &end);
    AssertThrow (ierr == 0, ExcPETScError(ierr));
    
    return ((index >= static_cast<unsigned int>(begin)) &&
            (index < static_cast<unsigned int>(end)));
  }

/// @endif      
}


#endif // DEAL_II_USE_PETSC


/*----------------------------   petsc_matrix_base.h     ---------------------------*/

#endif
/*----------------------------   petsc_matrix_base.h     ---------------------------*/
