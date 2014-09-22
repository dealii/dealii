// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2013 by the deal.II authors
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

#ifndef __deal2__block_matrix_base_h
#define __deal2__block_matrix_base_h


#include <deal.II/base/config.h>
#include <deal.II/base/table.h>
#include <deal.II/base/thread_management.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/memory_consumption.h>
#include <deal.II/lac/block_indices.h>
#include <deal.II/lac/exceptions.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/matrix_iterator.h>
#include <deal.II/lac/vector.h>

#include <cmath>

DEAL_II_NAMESPACE_OPEN


template <typename> class MatrixIterator;



/*! @addtogroup Matrix1
 *@{
 */

/**
 * Namespace in which iterators in block matrices are implemented.
 *
 * @author Wolfgang Bangerth, 2004
 */
namespace BlockMatrixIterators
{
  /**
   * Base class for block matrix
   * accessors, implementing the
   * stepping through a matrix.
   */
  template <class BlockMatrix>
  class AccessorBase
  {
  public:
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * Typedef the value type of the
     * matrix we point into.
     */
    typedef typename BlockMatrix::value_type value_type;

    /**
     * Initialize data fields to
     * default values.
     */
    AccessorBase ();

    /**
     * Block row of the
     * element represented by
     * this object.
     */
    unsigned int block_row() const;

    /**
     * Block column of the
     * element represented by
     * this object.
     */
    unsigned int block_column() const;

  protected:
    /**
     * Block row into which we presently
     * point.
     */
    unsigned int row_block;

    /**
     * Block column into which we
     * presently point.
     */
    unsigned int col_block;

    /**
     * Let the iterator class be a
     * friend.
     */
    template <typename>
    friend class MatrixIterator;
  };



  /**
   * Accessor classes in
   * block matrices.
   */
  template <class BlockMatrix, bool ConstNess>
  class Accessor;


  /**
   * Block matrix accessor for non
   * const matrices.
   */
  template <class BlockMatrix>
  class Accessor<BlockMatrix, false>
    :
    public AccessorBase<BlockMatrix>
  {
  public:
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * Type of the matrix used in
     * this accessor.
     */
    typedef BlockMatrix MatrixType;

    /**
     * Typedef the value type of the
     * matrix we point into.
     */
    typedef typename BlockMatrix::value_type value_type;

    /**
     * Constructor. Since we use
     * accessors only for read
     * access, a const matrix
     * pointer is sufficient.
     *
     * Place the iterator at the
     * beginning of the given row of the
     * matrix, or create the end pointer
     * if @p row equals the total number
     * of rows in the matrix.
     */
    Accessor (BlockMatrix *m,
              const size_type row,
              const size_type col);

    /**
     * Row number of the element
     * represented by this
     * object.
     */
    size_type row() const;

    /**
     * Column number of the
     * element represented by
     * this object.
     */
    size_type column() const;

    /**
     * Value of the entry at the
     * current position.
     */
    value_type value() const;

    /**
     * Set new value.
     */
    void set_value(value_type newval) const;

  protected:
    /**
     * The matrix accessed.
     */
    BlockMatrix *matrix;

    /**
     * Iterator of the underlying matrix
     * class.
     */
    typename BlockMatrix::BlockType::iterator base_iterator;

    /**
     * Move ahead one element.
     */
    void advance ();

    /**
     * Compare this accessor with another
     * one for equality.
     */
    bool operator == (const Accessor &a) const;

    template <typename> friend class MatrixIterator;
    friend class Accessor<BlockMatrix, true>;
  };

  /**
   * Block matrix accessor for
   * constant matrices, implementing
   * the stepping through a matrix.
   */
  template <class BlockMatrix>
  class Accessor<BlockMatrix, true>
    :
    public AccessorBase<BlockMatrix>
  {
  public:
    /**
     * Declare type for container size.
     */
    typedef types::global_dof_index size_type;

    /**
     * Type of the matrix used in
     * this accessor.
     */
    typedef const BlockMatrix MatrixType;

    /**
     * Typedef the value type of the
     * matrix we point into.
     */
    typedef typename BlockMatrix::value_type value_type;

    /**
     * Constructor. Since we use
     * accessors only for read
     * access, a const matrix
     * pointer is sufficient.
     *
     * Place the iterator at the
     * beginning of the given row of the
     * matrix, or create the end pointer
     * if @p row equals the total number
     * of rows in the matrix.
     */
    Accessor (const BlockMatrix *m,
              const size_type row,
              const size_type col);

    /**
     * Initalize const accessor
     * from non const accessor.
     */
    Accessor(const Accessor<BlockMatrix, false> &);

    /**
     * Row number of the element
     * represented by this
     * object.
     */
    size_type row() const;

    /**
     * Column number of the
     * element represented by
     * this object.
     */
    size_type column() const;

    /**
     * Value of the entry at the
     * current position.
     */
    value_type value() const;
  protected:
    /**
     * The matrix accessed.
     */
    const BlockMatrix *matrix;

    /**
     * Iterator of the underlying matrix
     * class.
     */
    typename BlockMatrix::BlockType::const_iterator base_iterator;

    /**
     * Move ahead one element.
     */
    void advance ();

    /**
     * Compare this accessor with another
     * one for equality.
     */
    bool operator == (const Accessor &a) const;

    /**
     * Let the iterator class be a
     * friend.
     */
    template <typename>
    friend class dealii::MatrixIterator;
  };
}



/**
 * Blocked matrix class. The behaviour of objects of this type is almost as
 * for the usual matrix objects, with most of the functions being implemented
 * in both classes. The main difference is that the matrix represented by this
 * object is composed of an array of matrices (e.g. of type
 * SparseMatrix<number>) and all accesses to the elements of this object are
 * relayed to accesses of the base matrices. The actual type of the individual
 * blocks of this matrix is the type of the template argument, and can, for
 * example be the usual SparseMatrix or PETScWrappers::SparseMatrix.
 *
 * In addition to the usual matrix access and linear algebra
 * functions, there are functions block() which allow access to the
 * different blocks of the matrix. This may, for example, be of help
 * when you want to implement Schur complement methods, or block
 * preconditioners, where each block belongs to a specific component
 * of the equation you are presently discretizing.
 *
 * Note that the numbers of blocks and rows are implicitly determined
 * by the sparsity pattern objects used.
 *
 * Objects of this type are frequently used when a system of differential
 * equations has solutions with variables that fall into different
 * classes. For example, solutions of the Stokes or Navier-Stokes equations
 * have @p dim velocity components and one pressure component. In this case,
 * it may make sense to consider the linear system of equations as a system of
 * 2x2 blocks, and one can construct preconditioners or solvers based on this
 * 2x2 block structure. This class can help you in these cases, as it allows
 * to view the matrix alternatively as one big matrix, or as a number of
 * individual blocks.
 *
 *
 * <h3>Inheriting from this class</h3>
 *
 * Since this class simply forwards its calls to the subobjects (if necessary
 * after adjusting indices denoting which subobject is meant), this class is
 * completely independent of the actual type of the subobject. The functions
 * that set up block matrices and destroy them, however, have to be
 * implemented in derived classes. These functions also have to fill the data
 * members provided by this base class, as they are only used passively in
 * this class.
 *
 *
 * Most of the functions take a vector or block vector argument. These
 * functions can, in general, only successfully be compiled if the
 * individual blocks of this matrix implement the respective functions
 * operating on the vector type in question. For example, if you have
 * a block sparse matrix over deal.II SparseMatrix objects, then you
 * will likely not be able to form the matrix-vector multiplication
 * with a block vector over PETScWrappers::SparseMatrix objects. If
 * you attempt anyway, you will likely get a number of compiler
 * errors.
 *
 * @note Instantiations for this template are provided for <tt>@<float@> and
 * @<double@></tt>; others can be generated in application programs (see the
 * section on @ref Instantiations in the manual).
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, 2000, 2004
 */
template <typename MatrixType>
class BlockMatrixBase : public Subscriptor
{
public:
  /**
   * Typedef the type of the underlying
   * matrix.
   */
  typedef MatrixType BlockType;

  /**
   * Type of matrix entries. In analogy to
   * the STL container classes.
   */
  typedef typename BlockType::value_type value_type;
  typedef value_type             *pointer;
  typedef const value_type       *const_pointer;
  typedef value_type             &reference;
  typedef const value_type       &const_reference;
  typedef types::global_dof_index size_type;

  typedef
  MatrixIterator<BlockMatrixIterators::Accessor<BlockMatrixBase, false> >
  iterator;

  typedef
  MatrixIterator<BlockMatrixIterators::Accessor<BlockMatrixBase, true> >
  const_iterator;


  /**
   * Default constructor.
   */
  BlockMatrixBase ();

  /**
   * Destructor.
   */
  ~BlockMatrixBase ();

  /**
   * Copy the given matrix to this
   * one.  The operation throws an
   * error if the sparsity patterns
   * of the two involved matrices
   * do not point to the same
   * object, since in this case the
   * copy operation is
   * cheaper. Since this operation
   * is notheless not for free, we
   * do not make it available
   * through operator=(), since
   * this may lead to unwanted
   * usage, e.g. in copy arguments
   * to functions, which should
   * really be arguments by
   * reference.
   *
   * The source matrix may be a
   * matrix of arbitrary type, as
   * long as its data type is
   * convertible to the data type
   * of this matrix.
   *
   * The function returns a
   * reference to <tt>this</tt>.
   */
  template <class BlockMatrixType>
  BlockMatrixBase &
  copy_from (const BlockMatrixType &source);

  /**
   * Access the block with the
   * given coordinates.
   */
  BlockType &
  block (const unsigned int row,
         const unsigned int column);


  /**
   * Access the block with the
   * given coordinates. Version for
   * constant objects.
   */
  const BlockType &
  block (const unsigned int row,
         const unsigned int column) const;

  /**
   * Return the dimension of the
   * image space.  To remember: the
   * matrix is of dimension
   * $m \times n$.
   */
  size_type m () const;

  /**
   * Return the dimension of the
   * range space.  To remember: the
   * matrix is of dimension
   * $m \times n$.
   */
  size_type n () const;


  /**
   * Return the number of blocks in
   * a column. Returns zero if no
   * sparsity pattern is presently
   * associated to this matrix.
   */
  unsigned int n_block_rows () const;

  /**
   * Return the number of blocks in
   * a row. Returns zero if no
   * sparsity pattern is presently
   * associated to this matrix.
   */
  unsigned int n_block_cols () const;

  /**
   * Set the element <tt>(i,j)</tt>
   * to <tt>value</tt>. Throws an
   * error if the entry does not
   * exist or if <tt>value</tt> is
   * not a finite number. Still, it
   * is allowed to store zero
   * values in non-existent fields.
   */
  void set (const size_type i,
            const size_type j,
            const value_type value);

  /**
   * Set all elements given in a
   * FullMatrix into the sparse matrix
   * locations given by
   * <tt>indices</tt>. In other words,
   * this function writes the elements
   * in <tt>full_matrix</tt> into the
   * calling matrix, using the
   * local-to-global indexing specified
   * by <tt>indices</tt> for both the
   * rows and the columns of the
   * matrix. This function assumes a
   * quadratic sparse matrix and a
   * quadratic full_matrix, the usual
   * situation in FE calculations.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be set anyway or
   * they should be filtered away (and
   * not change the previous content in
   * the respective element if it
   * exists). The default value is
   * <tt>false</tt>, i.e., even zero
   * values are treated.
   */
  template <typename number>
  void set (const std::vector<size_type> &indices,
            const FullMatrix<number>     &full_matrix,
            const bool                    elide_zero_values = false);

  /**
   * Same function as before, but now
   * including the possibility to use
   * rectangular full_matrices and
   * different local-to-global indexing
   * on rows and columns, respectively.
   */
  template <typename number>
  void set (const std::vector<size_type> &row_indices,
            const std::vector<size_type> &col_indices,
            const FullMatrix<number>     &full_matrix,
            const bool                    elide_zero_values = false);

  /**
   * Set several elements in the
   * specified row of the matrix with
   * column indices as given by
   * <tt>col_indices</tt> to the
   * respective value.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be set anyway or
   * they should be filtered away (and
   * not change the previous content in
   * the respective element if it
   * exists). The default value is
   * <tt>false</tt>, i.e., even zero
   * values are treated.
   */
  template <typename number>
  void set (const size_type row,
            const std::vector<size_type> &col_indices,
            const std::vector<number>    &values,
            const bool                    elide_zero_values = false);

  /**
   * Set several elements to values
   * given by <tt>values</tt> in a
   * given row in columns given by
   * col_indices into the sparse
   * matrix.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be inserted anyway
   * or they should be filtered
   * away. The default value is
   * <tt>false</tt>, i.e., even zero
   * values are inserted/replaced.
   */
  template <typename number>
  void set (const size_type  row,
            const size_type  n_cols,
            const size_type *col_indices,
            const number    *values,
            const bool       elide_zero_values = false);

  /**
   * Add <tt>value</tt> to the
   * element (<i>i,j</i>).  Throws
   * an error if the entry does not
   * exist or if <tt>value</tt> is
   * not a finite number. Still, it
   * is allowed to store zero
   * values in non-existent fields.
   */
  void add (const size_type i,
            const size_type j,
            const value_type value);

  /**
   * Add all elements given in a
   * FullMatrix<double> into sparse
   * matrix locations given by
   * <tt>indices</tt>. In other words,
   * this function adds the elements in
   * <tt>full_matrix</tt> to the
   * respective entries in calling
   * matrix, using the local-to-global
   * indexing specified by
   * <tt>indices</tt> for both the rows
   * and the columns of the
   * matrix. This function assumes a
   * quadratic sparse matrix and a
   * quadratic full_matrix, the usual
   * situation in FE calculations.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const std::vector<size_type> &indices,
            const FullMatrix<number>     &full_matrix,
            const bool                    elide_zero_values = true);

  /**
   * Same function as before, but now
   * including the possibility to use
   * rectangular full_matrices and
   * different local-to-global indexing
   * on rows and columns, respectively.
   */
  template <typename number>
  void add (const std::vector<size_type> &row_indices,
            const std::vector<size_type> &col_indices,
            const FullMatrix<number>     &full_matrix,
            const bool                    elide_zero_values = true);

  /**
   * Set several elements in the
   * specified row of the matrix with
   * column indices as given by
   * <tt>col_indices</tt> to the
   * respective value.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const size_type row,
            const std::vector<size_type> &col_indices,
            const std::vector<number>    &values,
            const bool                    elide_zero_values = true);

  /**
   * Add an array of values given by
   * <tt>values</tt> in the given
   * global matrix row at columns
   * specified by col_indices in the
   * sparse matrix.
   *
   * The optional parameter
   * <tt>elide_zero_values</tt> can be
   * used to specify whether zero
   * values should be added anyway or
   * these should be filtered away and
   * only non-zero data is added. The
   * default value is <tt>true</tt>,
   * i.e., zero values won't be added
   * into the matrix.
   */
  template <typename number>
  void add (const size_type  row,
            const size_type  n_cols,
            const size_type *col_indices,
            const number    *values,
            const bool       elide_zero_values = true,
            const bool       col_indices_are_sorted = false);

  /**
   * Add <tt>matrix</tt> scaled by
   * <tt>factor</tt> to this matrix,
   * i.e. the matrix
   * <tt>factor*matrix</tt> is added to
   * <tt>this</tt>. If the sparsity
   * pattern of the calling matrix does
   * not contain all the elements in
   * the sparsity pattern of the input
   * matrix, this function will throw
   * an exception.
   */
  void add (const value_type                   factor,
            const BlockMatrixBase<MatrixType> &matrix);

  /**
   * Return the value of the entry
   * (i,j).  This may be an
   * expensive operation and you
   * should always take care where
   * to call this function.  In
   * order to avoid abuse, this
   * function throws an exception
   * if the wanted element does not
   * exist in the matrix.
   */
  value_type operator () (const size_type i,
                          const size_type j) const;

  /**
   * This function is mostly like
   * operator()() in that it
   * returns the value of the
   * matrix entry <tt>(i,j)</tt>. The only
   * difference is that if this
   * entry does not exist in the
   * sparsity pattern, then instead
   * of raising an exception, zero
   * is returned. While this may be
   * convenient in some cases, note
   * that it is simple to write
   * algorithms that are slow
   * compared to an optimal
   * solution, since the sparsity
   * of the matrix is not used.
   */
  value_type el (const size_type i,
                 const size_type j) const;

  /**
   * Return the main diagonal element in
   * the <i>i</i>th row. This function
   * throws an error if the matrix is not
   * quadratic and also if the diagonal
   * blocks of the matrix are not
   * quadratic.
   *
   * This function is considerably
   * faster than the operator()(),
   * since for quadratic matrices, the
   * diagonal entry may be the
   * first to be stored in each row
   * and access therefore does not
   * involve searching for the
   * right column number.
   */
  value_type diag_element (const size_type i) const;

  /**
   * Call the compress() function on all
   * the subblocks of the matrix.
  *
  *
  * See @ref GlossCompress "Compressing
  * distributed objects" for more
  * information.
   */
  void compress (::dealii::VectorOperation::values operation);

  /**
   * @deprecated: use compress() with VectorOperation instead.
   */
  void compress () DEAL_II_DEPRECATED;

  /**
   * Multiply the entire matrix by a
   * fixed factor.
   */
  BlockMatrixBase &operator *= (const value_type factor);

  /**
   * Divide the entire matrix by a
   * fixed factor.
   */
  BlockMatrixBase &operator /= (const value_type factor);

  /**
   * Adding Matrix-vector
   * multiplication. Add $M*src$ on
   * $dst$ with $M$ being this
   * matrix.
   */
  template <class BlockVectorType>
  void vmult_add (BlockVectorType       &dst,
                  const BlockVectorType &src) const;

  /**
   * Adding Matrix-vector
   * multiplication. Add
   * <i>M<sup>T</sup>src</i> to
   * <i>dst</i> with <i>M</i> being
   * this matrix. This function
   * does the same as vmult_add()
   * but takes the transposed
   * matrix.
   */
  template <class BlockVectorType>
  void Tvmult_add (BlockVectorType       &dst,
                   const BlockVectorType &src) const;

  /**
   * Return the norm of the vector
   * <i>v</i> with respect to the
   * norm induced by this matrix,
   * i.e. <i>v<sup>T</sup>Mv)</i>. This
   * is useful, e.g. in the finite
   * element context, where the
   * <i>L<sup>T</sup></i>-norm of a
   * function equals the matrix
   * norm with respect to the mass
   * matrix of the vector
   * representing the nodal values
   * of the finite element
   * function. Note that even
   * though the function's name
   * might suggest something
   * different, for historic
   * reasons not the norm but its
   * square is returned, as defined
   * above by the scalar product.
   *
   * Obviously, the matrix needs to
   * be square for this operation.
   */
  template <class BlockVectorType>
  value_type
  matrix_norm_square (const BlockVectorType &v) const;

  /**
   * Compute the matrix scalar
   * product $\left(u,Mv\right)$.
   */
  template <class BlockVectorType>
  value_type
  matrix_scalar_product (const BlockVectorType &u,
                         const BlockVectorType &v) const;

  /**
   * Compute the residual
   * <i>r=b-Ax</i>. Write the
   * residual into <tt>dst</tt>.
   */
  template <class BlockVectorType>
  value_type residual (BlockVectorType       &dst,
                       const BlockVectorType &x,
                       const BlockVectorType &b) const;

  /**
   * Print the matrix to the given
   * stream, using the format
   * <tt>(line,col) value</tt>, i.e. one
   * nonzero entry of the matrix per
   * line. The optional flag outputs the
   * sparsity pattern in a different style
   * according to the underlying
   * sparsematrix type.
   */
  void print (std::ostream &out,
              const bool    alternative_output = false) const;

  /**
   * STL-like iterator with the
   * first entry.
   */
  iterator begin ();

  /**
   * Final iterator.
   */
  iterator end ();

  /**
   * STL-like iterator with the
   * first entry of row <tt>r</tt>.
   */
  iterator begin (const size_type r);

  /**
   * Final iterator of row <tt>r</tt>.
   */
  iterator end (const size_type r);
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
   * first entry of row <tt>r</tt>.
   */
  const_iterator begin (const size_type r) const;

  /**
   * Final iterator of row <tt>r</tt>.
   */
  const_iterator end (const size_type r) const;

  /**
   * Return a reference to the underlying
   * BlockIndices data of the rows.
   */
  const BlockIndices &get_row_indices () const;

  /**
   * Return a reference to the underlying
   * BlockIndices data of the columns.
   */
  const BlockIndices &get_column_indices () const;

  /**
   * Determine an estimate for the memory
   * consumption (in bytes) of this
   * object. Note that only the memory
   * reserved on the current processor is
   * returned in case this is called in
   * an MPI-based program.
   */
  std::size_t memory_consumption () const;

  /** @addtogroup Exceptions
   * @{ */

  /**
   * Exception
   */
  DeclException4 (ExcIncompatibleRowNumbers,
                  int, int, int, int,
                  << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                  << arg3 << ',' << arg4 << "] have differing row numbers.");
  /**
   * Exception
   */
  DeclException4 (ExcIncompatibleColNumbers,
                  int, int, int, int,
                  << "The blocks [" << arg1 << ',' << arg2 << "] and ["
                  << arg3 << ',' << arg4 << "] have differing column numbers.");
  //@}
protected:
  /**
   * Release all memory and return
   * to a state just like after
   * having called the default
   * constructor. It also forgets
   * the sparsity pattern it was
   * previously tied to.
   *
   * This calls clear for all
   * sub-matrices and then resets this
   * object to have no blocks at all.
   *
   * This function is protected
   * since it may be necessary to
   * release additional structures.
   * A derived class can make it
   * public again, if it is
   * sufficient.
   */
  void clear ();

  /**
   * Index arrays for rows and columns.
   */
  BlockIndices row_block_indices;
  BlockIndices column_block_indices;

  /**
   * Array of sub-matrices.
   */
  Table<2,SmartPointer<BlockType, BlockMatrixBase<MatrixType> > > sub_objects;

  /**
   * This function collects the
   * sizes of the sub-objects and
   * stores them in internal
   * arrays, in order to be able to
   * relay global indices into the
   * matrix to indices into the
   * subobjects. You *must* call
   * this function each time after
   * you have changed the size of
   * the sub-objects.
   *
   * Derived classes should call this
   * function whenever the size of the
   * sub-objects has changed and the @p
   * X_block_indices arrays need to be
   * updated.
   *
   * Note that this function is not public
   * since not all derived classes need to
   * export its interface. For example, for
   * the usual deal.II SparseMatrix class,
   * the sizes are implicitly determined
   * whenever reinit() is called, and
   * individual blocks cannot be
   * resized. For that class, this function
   * therefore does not have to be
   * public. On the other hand, for the
   * PETSc classes, there is no associated
   * sparsity pattern object that
   * determines the block sizes, and for
   * these the function needs to be
   * publicly available. These classes
   * therefore export this function.
   */
  void collect_sizes ();

  /**
   * Matrix-vector multiplication:
   * let $dst = M*src$ with $M$
   * being this matrix.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType>
  void vmult_block_block (BlockVectorType       &dst,
                          const BlockVectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block column.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType,
            class VectorType>
  void vmult_block_nonblock (BlockVectorType          &dst,
                             const VectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block row.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType,
            class VectorType>
  void vmult_nonblock_block (VectorType    &dst,
                             const BlockVectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class VectorType>
  void vmult_nonblock_nonblock (VectorType       &dst,
                                const VectorType &src) const;

  /**
   * Matrix-vector multiplication:
   * let $dst = M^T*src$ with $M$
   * being this matrix. This
   * function does the same as
   * vmult() but takes the
   * transposed matrix.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType>
  void Tvmult_block_block (BlockVectorType       &dst,
                           const BlockVectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block row.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType,
            class VectorType>
  void Tvmult_block_nonblock (BlockVectorType  &dst,
                              const VectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block column.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class BlockVectorType,
            class VectorType>
  void Tvmult_nonblock_block (VectorType    &dst,
                              const BlockVectorType &src) const;

  /**
   * Matrix-vector
   * multiplication. Just like the
   * previous function, but only
   * applicable if the matrix has
   * only one block.
   *
   * Due to problems with deriving template
   * arguments between the block and
   * non-block versions of the vmult/Tvmult
   * functions, the actual functions are
   * implemented in derived classes, with
   * implementations forwarding the calls
   * to the implementations provided here
   * under a unique name for which template
   * arguments can be derived by the
   * compiler.
   */
  template <class VectorType>
  void Tvmult_nonblock_nonblock (VectorType       &dst,
                                 const VectorType &src) const;


protected:

  /**
   * Some matrix types, in particular PETSc,
   * need to synchronize set and add
   * operations. This has to be done for all
   * matrices in the BlockMatrix.
   * This routine prepares adding of elements
   * by notifying all blocks. Called by all
   * internal routines before adding
   * elements.
   */
  void prepare_add_operation();

  /**
   * Notifies all blocks to let them prepare
   * for setting elements, see
   * prepare_add_operation().
   */
  void prepare_set_operation();


private:

  /**
   * A structure containing some fields used by the
   * set() and add() functions that is used to pre-sort
   * the input fields. Since one can reasonably expect
   * to call set() and add() from multiple threads at once
   * as long as the matrix indices that are touched are
   * disjoint, these temporary data fields need to be
   * guarded by a mutex; the structure therefore contains such
   * a mutex as a member variable.
   */
  struct TemporaryData
  {
    /**
     * Temporary vector for counting the
     * elements written into the
     * individual blocks when doing a
     * collective add or set.
     */
    std::vector<size_type> counter_within_block;

    /**
     * Temporary vector for column
     * indices on each block when writing
     * local to global data on each
     * sparse matrix.
     */
    std::vector<std::vector<size_type> > column_indices;

    /**
     * Temporary vector for storing the
     * local values (they need to be
     * reordered when writing local to
     * global).
     */
    std::vector<std::vector<double> > column_values;

    /**
     * A mutex variable used to guard access to the member
     * variables of this structure;
     */
    Threads::Mutex mutex;

    /**
     * Copy operator. This is needed because the default copy
     * operator of this class is deleted (since Threads::Mutex is
     * not copyable) and hence the default copy operator of the
     * enclosing class is also deleted.
     *
     * The implementation here simply does nothing -- TemporaryData
     * objects are just scratch objects that are resized at the
     * beginning of their use, so there is no point actually copying
     * anything.
     */
    TemporaryData &operator = (const TemporaryData &)
    {
      return *this;
    }
  };

  /**
   * A set of scratch arrays that can be used by the add()
   * and set() functions that take pointers to data to
   * pre-sort indices before use. Access from multiple threads
   * is synchronized via the mutex variable that is part of the
   * structure.
   */
  TemporaryData temporary_data;

  /**
   * Make the iterator class a
   * friend. We have to work around
   * a compiler bug here again.
   */
  template <typename, bool>
  friend class BlockMatrixIterators::Accessor;

  template <typename>
  friend class MatrixIterator;
};


/*@}*/

#ifndef DOXYGEN
/* ------------------------- Template functions ---------------------- */


namespace BlockMatrixIterators
{
  template <class BlockMatrix>
  inline
  AccessorBase<BlockMatrix>::AccessorBase()
    :
    row_block(0),
    col_block(0)
  {}


  template <class BlockMatrix>
  inline
  unsigned int
  AccessorBase<BlockMatrix>::block_row() const
  {
    Assert (row_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    return row_block;
  }


  template <class BlockMatrix>
  inline
  unsigned int
  AccessorBase<BlockMatrix>::block_column() const
  {
    Assert (col_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    return col_block;
  }


  template <class BlockMatrix>
  inline
  Accessor<BlockMatrix, true>::Accessor (
    const BlockMatrix  *matrix,
    const size_type     row,
    const size_type     col)
    :
    matrix(matrix),
    base_iterator(matrix->block(0,0).begin())
  {
    Assert(col==0, ExcNotImplemented());

    // check if this is a regular row or
    // the end of the matrix
    if (row < matrix->m())
      {
        const std::pair<unsigned int,size_type> indices
          = matrix->row_block_indices.global_to_local(row);

        // find the first block that does
        // have an entry in this row
        for (unsigned int bc=0; bc<matrix->n_block_cols(); ++bc)
          {
            base_iterator
              = matrix->block(indices.first, bc).begin(indices.second);
            if (base_iterator !=
                matrix->block(indices.first, bc).end(indices.second))
              {
                this->row_block = indices.first;
                this->col_block = bc;
                return;
              }
          }

        // hm, there is no block that has
        // an entry in this column. we need
        // to take the next entry then,
        // which may be the first entry of
        // the next row, or recursively the
        // next row, or so on
        *this = Accessor (matrix, row+1, 0);
      }
    else
      {
        // we were asked to create the end
        // iterator for this matrix
        this->row_block = numbers::invalid_unsigned_int;
        this->col_block = numbers::invalid_unsigned_int;
      }
  }


//   template <class BlockMatrix>
//   inline
//   Accessor<BlockMatrix, true>::Accessor (const Accessor<BlockMatrix, true>& other)
//                :
//                matrix(other.matrix),
//                base_iterator(other.base_iterator)
//   {
//     this->row_block = other.row_block;
//     this->col_block = other.col_block;
//   }


  template <class BlockMatrix>
  inline
  Accessor<BlockMatrix, true>::Accessor (const Accessor<BlockMatrix, false> &other)
    :
    matrix(other.matrix),
    base_iterator(other.base_iterator)
  {
    this->row_block = other.row_block;
    this->col_block = other.col_block;
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, true>::size_type
  Accessor<BlockMatrix, true>::row() const
  {
    Assert (this->row_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    return (matrix->row_block_indices.local_to_global(this->row_block, 0) +
            base_iterator->row());
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, true>::size_type
  Accessor<BlockMatrix, true>::column() const
  {
    Assert (this->col_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    return (matrix->column_block_indices.local_to_global(this->col_block,0) +
            base_iterator->column());
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, true>::value_type
  Accessor<BlockMatrix, true>::value () const
  {
    Assert (this->row_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());
    Assert (this->col_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    return base_iterator->value();
  }



  template <class BlockMatrix>
  inline
  void
  Accessor<BlockMatrix, true>::advance ()
  {
    Assert (this->row_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());
    Assert (this->col_block != numbers::invalid_unsigned_int,
            ExcIteratorPastEnd());

    // Remember current row inside block
    size_type local_row = base_iterator->row();

    // Advance one element inside the
    // current block
    ++base_iterator;

    // while we hit the end of the row of a
    // block (which may happen multiple
    // times if rows inside a block are
    // empty), we have to jump to the next
    // block and take the
    while (base_iterator ==
           matrix->block(this->row_block, this->col_block).end(local_row))
      {
        // jump to next block in this block
        // row, if possible, otherwise go
        // to next row
        if (this->col_block < matrix->n_block_cols()-1)
          {
            ++this->col_block;
            base_iterator
              = matrix->block(this->row_block, this->col_block).begin(local_row);
          }
        else
          {
            // jump back to next row in
            // first block column
            this->col_block = 0;
            ++local_row;

            // see if this has brought us
            // past the number of rows in
            // this block. if so see
            // whether we've just fallen
            // off the end of the whole
            // matrix
            if (local_row == matrix->block(this->row_block, this->col_block).m())
              {
                local_row = 0;
                ++this->row_block;
                if (this->row_block == matrix->n_block_rows())
                  {
                    this->row_block = numbers::invalid_unsigned_int;
                    this->col_block = numbers::invalid_unsigned_int;
                    return;
                  }
              }

            base_iterator
              = matrix->block(this->row_block, this->col_block).begin(local_row);
          }
      }
  }


  template <class BlockMatrix>
  inline
  bool
  Accessor<BlockMatrix, true>::operator == (const Accessor &a) const
  {
    if (matrix != a.matrix)
      return false;

    if (this->row_block == a.row_block
        && this->col_block == a.col_block)
      // end iterators do not necessarily
      // have to have the same
      // base_iterator representation, but
      // valid iterators have to
      return (((this->row_block == numbers::invalid_unsigned_int)
               &&
               (this->col_block == numbers::invalid_unsigned_int))
              ||
              (base_iterator == a.base_iterator));

    return false;
  }

//----------------------------------------------------------------------//


  template <class BlockMatrix>
  inline
  Accessor<BlockMatrix, false>::Accessor (
    BlockMatrix  *matrix,
    const size_type row,
    const size_type col)
    :
    matrix(matrix),
    base_iterator(matrix->block(0,0).begin())
  {
    Assert(col==0, ExcNotImplemented());
    // check if this is a regular row or
    // the end of the matrix
    if (row < matrix->m())
      {
        const std::pair<unsigned int,size_type> indices
          = matrix->row_block_indices.global_to_local(row);

        // find the first block that does
        // have an entry in this row
        for (size_type bc=0; bc<matrix->n_block_cols(); ++bc)
          {
            base_iterator
              = matrix->block(indices.first, bc).begin(indices.second);
            if (base_iterator !=
                matrix->block(indices.first, bc).end(indices.second))
              {
                this->row_block = indices.first;
                this->col_block = bc;
                return;
              }
          }

        // hm, there is no block that has
        // an entry in this column. we need
        // to take the next entry then,
        // which may be the first entry of
        // the next row, or recursively the
        // next row, or so on
        *this = Accessor (matrix, row+1, 0);
      }
    else
      {
        // we were asked to create the end
        // iterator for this matrix
        this->row_block = numbers::invalid_size_type;
        this->col_block = numbers::invalid_size_type;
      }
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, false>::size_type
  Accessor<BlockMatrix, false>::row() const
  {
    Assert (this->row_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());

    return (matrix->row_block_indices.local_to_global(this->row_block, 0) +
            base_iterator->row());
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, false>::size_type
  Accessor<BlockMatrix, false>::column() const
  {
    Assert (this->col_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());

    return (matrix->column_block_indices.local_to_global(this->col_block,0) +
            base_iterator->column());
  }


  template <class BlockMatrix>
  inline
  typename Accessor<BlockMatrix, false>::value_type
  Accessor<BlockMatrix, false>::value () const
  {
    Assert (this->row_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());
    Assert (this->col_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());

    return base_iterator->value();
  }



  template <class BlockMatrix>
  inline
  void
  Accessor<BlockMatrix, false>::set_value (typename Accessor<BlockMatrix, false>::value_type newval) const
  {
    Assert (this->row_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());
    Assert (this->col_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());

    base_iterator->value() = newval;
  }



  template <class BlockMatrix>
  inline
  void
  Accessor<BlockMatrix, false>::advance ()
  {
    Assert (this->row_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());
    Assert (this->col_block != numbers::invalid_size_type,
            ExcIteratorPastEnd());

    // Remember current row inside block
    size_type local_row = base_iterator->row();

    // Advance one element inside the
    // current block
    ++base_iterator;

    // while we hit the end of the row of a
    // block (which may happen multiple
    // times if rows inside a block are
    // empty), we have to jump to the next
    // block and take the
    while (base_iterator ==
           matrix->block(this->row_block, this->col_block).end(local_row))
      {
        // jump to next block in this block
        // row, if possible, otherwise go
        // to next row
        if (this->col_block < matrix->n_block_cols()-1)
          {
            ++this->col_block;
            base_iterator
              = matrix->block(this->row_block, this->col_block).begin(local_row);
          }
        else
          {
            // jump back to next row in
            // first block column
            this->col_block = 0;
            ++local_row;

            // see if this has brought us
            // past the number of rows in
            // this block. if so see
            // whether we've just fallen
            // off the end of the whole
            // matrix
            if (local_row == matrix->block(this->row_block, this->col_block).m())
              {
                local_row = 0;
                ++this->row_block;
                if (this->row_block == matrix->n_block_rows())
                  {
                    this->row_block = numbers::invalid_size_type;
                    this->col_block = numbers::invalid_size_type;
                    return;
                  }
              }

            base_iterator
              = matrix->block(this->row_block, this->col_block).begin(local_row);
          }
      }
  }



  template <class BlockMatrix>
  inline
  bool
  Accessor<BlockMatrix, false>::operator == (const Accessor &a) const
  {
    if (matrix != a.matrix)
      return false;

    if (this->row_block == a.row_block
        && this->col_block == a.col_block)
      // end iterators do not necessarily
      // have to have the same
      // base_iterator representation, but
      // valid iterators have to
      return (((this->row_block == numbers::invalid_size_type)
               &&
               (this->col_block == numbers::invalid_size_type))
              ||
              (base_iterator == a.base_iterator));

    return false;
  }
}


//---------------------------------------------------------------------------


template <typename MatrixType>
inline
BlockMatrixBase<MatrixType>::BlockMatrixBase ()
{}

template <typename MatrixType>
inline
BlockMatrixBase<MatrixType>::~BlockMatrixBase ()
{
  clear ();
}


template <class MatrixType>
template <class BlockMatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::
copy_from (const BlockMatrixType &source)
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).copy_from (source.block(r,c));

  return *this;
}


template <class MatrixType>
std::size_t
BlockMatrixBase<MatrixType>::memory_consumption () const
{
  std::size_t mem =
    MemoryConsumption::memory_consumption(row_block_indices)+
    MemoryConsumption::memory_consumption(column_block_indices)+
    MemoryConsumption::memory_consumption(sub_objects)+
    MemoryConsumption::memory_consumption(temporary_data.counter_within_block)+
    MemoryConsumption::memory_consumption(temporary_data.column_indices)+
    MemoryConsumption::memory_consumption(temporary_data.column_values)+
    sizeof(temporary_data.mutex);

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      {
        MatrixType *p = this->sub_objects[r][c];
        mem += MemoryConsumption::memory_consumption(*p);
      }

  return mem;
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::clear ()
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      {
        MatrixType *p = this->sub_objects[r][c];
        this->sub_objects[r][c] = 0;
        delete p;
      }
  sub_objects.reinit (0,0);

  // reset block indices to empty
  row_block_indices = column_block_indices = BlockIndices ();
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block (const unsigned int row,
                                    const unsigned int column)
{
  Assert (row<n_block_rows(),
          ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(),
          ExcIndexRange (column, 0, n_block_cols()));

  return *sub_objects[row][column];
}



template <class MatrixType>
inline
const typename BlockMatrixBase<MatrixType>::BlockType &
BlockMatrixBase<MatrixType>::block (const unsigned int row,
                                    const unsigned int column) const
{
  Assert (row<n_block_rows(),
          ExcIndexRange (row, 0, n_block_rows()));
  Assert (column<n_block_cols(),
          ExcIndexRange (column, 0, n_block_cols()));

  return *sub_objects[row][column];
}


template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::size_type
BlockMatrixBase<MatrixType>::m () const
{
  return row_block_indices.total_size();
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::size_type
BlockMatrixBase<MatrixType>::n () const
{
  return column_block_indices.total_size();
}



template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::n_block_cols () const
{
  return column_block_indices.size();
}



template <class MatrixType>
inline
unsigned int
BlockMatrixBase<MatrixType>::n_block_rows () const
{
  return row_block_indices.size();
}



// Write the single set manually,
// since the other function has a lot
// of overhead in that case.
template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::set (const size_type i,
                                  const size_type j,
                                  const value_type value)
{
  prepare_set_operation();

  Assert (numbers::is_finite(value), ExcNumberNotFinite());

  const std::pair<unsigned int,size_type>
  row_index = row_block_indices.global_to_local (i),
  col_index = column_block_indices.global_to_local (j);
  block(row_index.first,col_index.first).set (row_index.second,
                                              col_index.second,
                                              value);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::set (const std::vector<size_type> &row_indices,
                                  const std::vector<size_type> &col_indices,
                                  const FullMatrix<number>        &values,
                                  const bool                       elide_zero_values)
{
  Assert (row_indices.size() == values.m(),
          ExcDimensionMismatch(row_indices.size(), values.m()));
  Assert (col_indices.size() == values.n(),
          ExcDimensionMismatch(col_indices.size(), values.n()));

  for (size_type i=0; i<row_indices.size(); ++i)
    set (row_indices[i], col_indices.size(), &col_indices[0], &values(i,0),
         elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::set (const std::vector<size_type> &indices,
                                  const FullMatrix<number>        &values,
                                  const bool                       elide_zero_values)
{
  Assert (indices.size() == values.m(),
          ExcDimensionMismatch(indices.size(), values.m()));
  Assert (values.n() == values.m(), ExcNotQuadratic());

  for (size_type i=0; i<indices.size(); ++i)
    set (indices[i], indices.size(), &indices[0], &values(i,0),
         elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::set (const size_type               row,
                                  const std::vector<size_type> &col_indices,
                                  const std::vector<number>    &values,
                                  const bool                    elide_zero_values)
{
  Assert (col_indices.size() == values.size(),
          ExcDimensionMismatch(col_indices.size(), values.size()));

  set (row, col_indices.size(), &col_indices[0], &values[0],
       elide_zero_values);
}



// This is a very messy function, since
// we need to calculate to each position
// the location in the global array.
template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::set (const size_type  row,
                                  const size_type  n_cols,
                                  const size_type *col_indices,
                                  const number    *values,
                                  const bool       elide_zero_values)
{
  prepare_set_operation();

  // lock access to the temporary data structure to
  // allow multiple threads to call this function concurrently
  Threads::Mutex::ScopedLock lock (temporary_data.mutex);

  // Resize scratch arrays
  if (temporary_data.column_indices.size() < this->n_block_cols())
    {
      temporary_data.column_indices.resize (this->n_block_cols());
      temporary_data.column_values.resize (this->n_block_cols());
      temporary_data.counter_within_block.resize (this->n_block_cols());
    }

  // Resize sub-arrays to n_cols. This
  // is a bit wasteful, but we resize
  // only a few times (then the maximum
  // row length won't increase that
  // much any more). At least we know
  // that all arrays are going to be of
  // the same size, so we can check
  // whether the size of one is large
  // enough before actually going
  // through all of them.
  if (temporary_data.column_indices[0].size() < n_cols)
    {
      for (unsigned int i=0; i<this->n_block_cols(); ++i)
        {
          temporary_data.column_indices[i].resize(n_cols);
          temporary_data.column_values[i].resize(n_cols);
        }
    }

  // Reset the number of added elements
  // in each block to zero.
  for (unsigned int i=0; i<this->n_block_cols(); ++i)
    temporary_data.counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be set in which
  // block of the matrix. We need to
  // touch all the data, since we can't
  // be sure that the data of one block
  // is stored contiguously (in fact,
  // indices will be intermixed when it
  // comes from an element matrix).
  for (size_type j=0; j<n_cols; ++j)
    {
      double value = values[j];

      if (value == 0 && elide_zero_values == true)
        continue;

      const std::pair<unsigned int, size_type>
      col_index = this->column_block_indices.global_to_local(col_indices[j]);

      const size_type local_index = temporary_data.counter_within_block[col_index.first]++;

      temporary_data.column_indices[col_index.first][local_index] = col_index.second;
      temporary_data.column_values[col_index.first][local_index] = value;
    }

#ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (unsigned int i=0; i<this->n_block_cols(); ++i)
    length += temporary_data.counter_within_block[i];
  Assert (length <= n_cols, ExcInternalError());
#endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<unsigned int,size_type>
  row_index = this->row_block_indices.global_to_local (row);
  for (unsigned int block_col=0; block_col<n_block_cols(); ++block_col)
    {
      if (temporary_data.counter_within_block[block_col] == 0)
        continue;

      block(row_index.first, block_col).set
      (row_index.second,
       temporary_data.counter_within_block[block_col],
       &temporary_data.column_indices[block_col][0],
       &temporary_data.column_values[block_col][0],
       false);
    }
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::add (const size_type  i,
                                  const size_type  j,
                                  const value_type value)
{

  Assert (numbers::is_finite(value), ExcNumberNotFinite());

  prepare_add_operation();

  // save some cycles for zero additions, but
  // only if it is safe for the matrix we are
  // working with
  typedef typename MatrixType::Traits MatrixTraits;
  if ((MatrixTraits::zero_addition_can_be_elided == true)
      &&
      (value == value_type()))
    return;

  const std::pair<unsigned int,size_type>
  row_index = row_block_indices.global_to_local (i),
  col_index = column_block_indices.global_to_local (j);
  block(row_index.first,col_index.first).add (row_index.second,
                                              col_index.second,
                                              value);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::add (const std::vector<size_type> &row_indices,
                                  const std::vector<size_type> &col_indices,
                                  const FullMatrix<number>     &values,
                                  const bool                    elide_zero_values)
{
  Assert (row_indices.size() == values.m(),
          ExcDimensionMismatch(row_indices.size(), values.m()));
  Assert (col_indices.size() == values.n(),
          ExcDimensionMismatch(col_indices.size(), values.n()));

  for (size_type i=0; i<row_indices.size(); ++i)
    add (row_indices[i], col_indices.size(), &col_indices[0], &values(i,0),
         elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::add (const std::vector<size_type> &indices,
                                  const FullMatrix<number>     &values,
                                  const bool                    elide_zero_values)
{
  Assert (indices.size() == values.m(),
          ExcDimensionMismatch(indices.size(), values.m()));
  Assert (values.n() == values.m(), ExcNotQuadratic());

  for (size_type i=0; i<indices.size(); ++i)
    add (indices[i], indices.size(), &indices[0], &values(i,0),
         elide_zero_values);
}



template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::add (const size_type               row,
                                  const std::vector<size_type> &col_indices,
                                  const std::vector<number>    &values,
                                  const bool                    elide_zero_values)
{
  Assert (col_indices.size() == values.size(),
          ExcDimensionMismatch(col_indices.size(), values.size()));

  add (row, col_indices.size(), &col_indices[0], &values[0],
       elide_zero_values);
}



// This is a very messy function, since
// we need to calculate to each position
// the location in the global array.
template <class MatrixType>
template <typename number>
inline
void
BlockMatrixBase<MatrixType>::add (const size_type  row,
                                  const size_type  n_cols,
                                  const size_type *col_indices,
                                  const number    *values,
                                  const bool       elide_zero_values,
                                  const bool       col_indices_are_sorted)
{
  prepare_add_operation();

  // TODO: Look over this to find out
  // whether we can do that more
  // efficiently.
  if (col_indices_are_sorted == true)
    {
#ifdef DEBUG
      // check whether indices really are
      // sorted.
      size_type before = col_indices[0];
      for (size_type i=1; i<n_cols; ++i)
        if (col_indices[i] <= before)
          Assert (false, ExcMessage ("Flag col_indices_are_sorted is set, but "
                                     "indices appear to not be sorted."))
          else
            before = col_indices[i];
#endif
      const std::pair<unsigned int,size_type>
      row_index = this->row_block_indices.global_to_local (row);

      if (this->n_block_cols() > 1)
        {
          const size_type *first_block = Utilities::lower_bound (col_indices,
                                                                 col_indices+n_cols,
                                                                 this->column_block_indices.block_start(1));

          const size_type n_zero_block_indices = first_block - col_indices;
          block(row_index.first, 0).add (row_index.second,
                                         n_zero_block_indices,
                                         col_indices,
                                         values,
                                         elide_zero_values,
                                         col_indices_are_sorted);

          if (n_zero_block_indices < n_cols)
            this->add(row, n_cols - n_zero_block_indices, first_block,
                      values + n_zero_block_indices, elide_zero_values,
                      false);
        }
      else
        {
          block(row_index.first, 0). add (row_index.second,
                                          n_cols,
                                          col_indices,
                                          values,
                                          elide_zero_values,
                                          col_indices_are_sorted);
        }

      return;
    }

  // Lock scratch arrays, then resize them
  Threads::Mutex::ScopedLock lock (temporary_data.mutex);

  if (temporary_data.column_indices.size() < this->n_block_cols())
    {
      temporary_data.column_indices.resize (this->n_block_cols());
      temporary_data.column_values.resize (this->n_block_cols());
      temporary_data.counter_within_block.resize (this->n_block_cols());
    }

  // Resize sub-arrays to n_cols. This
  // is a bit wasteful, but we resize
  // only a few times (then the maximum
  // row length won't increase that
  // much any more). At least we know
  // that all arrays are going to be of
  // the same size, so we can check
  // whether the size of one is large
  // enough before actually going
  // through all of them.
  if (temporary_data.column_indices[0].size() < n_cols)
    {
      for (unsigned int i=0; i<this->n_block_cols(); ++i)
        {
          temporary_data.column_indices[i].resize(n_cols);
          temporary_data.column_values[i].resize(n_cols);
        }
    }

  // Reset the number of added elements
  // in each block to zero.
  for (unsigned int i=0; i<this->n_block_cols(); ++i)
    temporary_data.counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be written into
  // which block of the matrix. We need
  // to touch all the data, since we
  // can't be sure that the data of one
  // block is stored contiguously (in
  // fact, data will be intermixed when
  // it comes from an element matrix).
  for (size_type j=0; j<n_cols; ++j)
    {
      double value = values[j];

      if (value == 0 && elide_zero_values == true)
        continue;

      const std::pair<unsigned int, size_type>
      col_index = this->column_block_indices.global_to_local(col_indices[j]);

      const size_type local_index = temporary_data.counter_within_block[col_index.first]++;

      temporary_data.column_indices[col_index.first][local_index] = col_index.second;
      temporary_data.column_values[col_index.first][local_index] = value;
    }

#ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (unsigned int i=0; i<this->n_block_cols(); ++i)
    length += temporary_data.counter_within_block[i];
  Assert (length <= n_cols, ExcInternalError());
#endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<unsigned int,size_type>
  row_index = this->row_block_indices.global_to_local (row);
  for (unsigned int block_col=0; block_col<n_block_cols(); ++block_col)
    {
      if (temporary_data.counter_within_block[block_col] == 0)
        continue;

      block(row_index.first, block_col).add
      (row_index.second,
       temporary_data.counter_within_block[block_col],
       &temporary_data.column_indices[block_col][0],
       &temporary_data.column_values[block_col][0],
       false,
       col_indices_are_sorted);
    }
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::add (const value_type                   factor,
                                  const BlockMatrixBase<MatrixType> &matrix)
{
  Assert (numbers::is_finite(factor), ExcNumberNotFinite());

  prepare_add_operation();

  // save some cycles for zero additions, but
  // only if it is safe for the matrix we are
  // working with
  typedef typename MatrixType::Traits MatrixTraits;
  if ((MatrixTraits::zero_addition_can_be_elided == true)
      &&
      (factor == 0))
    return;

  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      // This function should throw if the sparsity
      // patterns of the two blocks differ
      block(row, col).add(factor, matrix.block(row,col));
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::operator () (const size_type i,
                                          const size_type j) const
{
  const std::pair<unsigned int,size_type>
  row_index = row_block_indices.global_to_local (i),
  col_index = column_block_indices.global_to_local (j);
  return block(row_index.first,col_index.first) (row_index.second,
                                                 col_index.second);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::el (const size_type i,
                                 const size_type j) const
{
  const std::pair<unsigned int,size_type>
  row_index = row_block_indices.global_to_local (i),
  col_index = column_block_indices.global_to_local (j);
  return block(row_index.first,col_index.first).el (row_index.second,
                                                    col_index.second);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::diag_element (const size_type i) const
{
  Assert (n_block_rows() == n_block_cols(),
          ExcNotQuadratic());

  const std::pair<unsigned int,size_type>
  index = row_block_indices.global_to_local (i);
  return block(index.first,index.first).diag_element(index.second);
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::compress (::dealii::VectorOperation::values operation)
{
  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c).compress (operation);
}

template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::compress ()
{
  compress(::dealii::VectorOperation::unknown);
}



template <class MatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator *= (const value_type factor)
{
  Assert (n_block_cols() != 0, ExcNotInitialized());
  Assert (n_block_rows() != 0, ExcNotInitialized());

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) *= factor;

  return *this;
}



template <class MatrixType>
inline
BlockMatrixBase<MatrixType> &
BlockMatrixBase<MatrixType>::operator /= (const value_type factor)
{
  Assert (n_block_cols() != 0, ExcNotInitialized());
  Assert (n_block_rows() != 0, ExcNotInitialized());
  Assert (factor !=0, ExcDivideByZero());

  const value_type factor_inv = 1. / factor;

  for (unsigned int r=0; r<n_block_rows(); ++r)
    for (unsigned int c=0; c<n_block_cols(); ++c)
      block(r,c) *= factor_inv;

  return *this;
}



template <class MatrixType>
const BlockIndices &
BlockMatrixBase<MatrixType>::get_row_indices () const
{
  return this->row_block_indices;
}



template <class MatrixType>
const BlockIndices &
BlockMatrixBase<MatrixType>::get_column_indices () const
{
  return this->column_block_indices;
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::
vmult_block_block (BlockVectorType       &dst,
                   const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (size_type row=0; row<n_block_rows(); ++row)
    {
      block(row,0).vmult (dst.block(row),
                          src.block(0));
      for (size_type col=1; col<n_block_cols(); ++col)
        block(row,col).vmult_add (dst.block(row),
                                  src.block(col));
    };
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::
vmult_nonblock_block (VectorType    &dst,
                      const BlockVectorType &src) const
{
  Assert (n_block_rows() == 1,
          ExcDimensionMismatch(1, n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  block(0,0).vmult (dst, src.block(0));
  for (size_type col=1; col<n_block_cols(); ++col)
    block(0,col).vmult_add (dst, src.block(col));
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::
vmult_block_nonblock (BlockVectorType  &dst,
                      const VectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (1 == n_block_cols(),
          ExcDimensionMismatch(1, n_block_cols()));

  for (size_type row=0; row<n_block_rows(); ++row)
    block(row,0).vmult (dst.block(row),
                        src);
}



template <class MatrixType>
template <class VectorType>
void
BlockMatrixBase<MatrixType>::
vmult_nonblock_nonblock (VectorType       &dst,
                         const VectorType &src) const
{
  Assert (1 == n_block_rows(),
          ExcDimensionMismatch(1, n_block_rows()));
  Assert (1 == n_block_cols(),
          ExcDimensionMismatch(1, n_block_cols()));

  block(0,0).vmult (dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::vmult_add (BlockVectorType       &dst,
                                        const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (src.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(src.n_blocks(), n_block_cols()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      block(row,col).vmult_add (dst.block(row),
                                src.block(col));
}




template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::
Tvmult_block_block (BlockVectorType       &dst,
                    const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  dst = 0.;

  for (unsigned int  row=0; row<n_block_rows(); ++row)
    {
      for (unsigned int col=0; col<n_block_cols(); ++col)
        block(row,col).Tvmult_add (dst.block(col),
                                   src.block(row));
    };
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::
Tvmult_block_nonblock (BlockVectorType  &dst,
                       const VectorType &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (1 == n_block_rows(),
          ExcDimensionMismatch(1, n_block_rows()));

  dst = 0.;

  for (unsigned int col=0; col<n_block_cols(); ++col)
    block(0,col).Tvmult_add (dst.block(col), src);
}



template <class MatrixType>
template <class BlockVectorType,
          class VectorType>
void
BlockMatrixBase<MatrixType>::
Tvmult_nonblock_block (VectorType    &dst,
                       const BlockVectorType &src) const
{
  Assert (1 == n_block_cols(),
          ExcDimensionMismatch(1, n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  block(0,0).Tvmult (dst, src.block(0));

  for (size_type row=1; row<n_block_rows(); ++row)
    block(row,0).Tvmult_add (dst, src.block(row));
}



template <class MatrixType>
template <class VectorType>
void
BlockMatrixBase<MatrixType>::
Tvmult_nonblock_nonblock (VectorType       &dst,
                          const VectorType &src) const
{
  Assert (1 == n_block_cols(),
          ExcDimensionMismatch(1, n_block_cols()));
  Assert (1 == n_block_rows(),
          ExcDimensionMismatch(1, n_block_rows()));

  block(0,0).Tvmult (dst, src);
}



template <class MatrixType>
template <class BlockVectorType>
void
BlockMatrixBase<MatrixType>::Tvmult_add (BlockVectorType &dst,
                                         const BlockVectorType &src) const
{
  Assert (dst.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_cols()));
  Assert (src.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(src.n_blocks(), n_block_rows()));

  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      block(row,col).Tvmult_add (dst.block(col),
                                 src.block(row));
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::matrix_norm_square (const BlockVectorType &v) const
{
  Assert (n_block_rows() == n_block_cols(), ExcNotQuadratic());
  Assert (v.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(v.n_blocks(), n_block_rows()));

  value_type norm_sqr = 0;
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      if (row==col)
        norm_sqr += block(row,col).matrix_norm_square (v.block(row));
      else
        norm_sqr += block(row,col).matrix_scalar_product (v.block(row),
                                                          v.block(col));
  return norm_sqr;
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::
matrix_scalar_product (const BlockVectorType    &u,
                       const BlockVectorType &v) const
{
  Assert (u.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(u.n_blocks(), n_block_rows()));
  Assert (v.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(v.n_blocks(), n_block_cols()));

  value_type result = 0;
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      result += block(row,col).matrix_scalar_product (u.block(row),
                                                      v.block(col));
  return result;
}



template <class MatrixType>
template <class BlockVectorType>
typename BlockMatrixBase<MatrixType>::value_type
BlockMatrixBase<MatrixType>::
residual (BlockVectorType          &dst,
          const BlockVectorType &x,
          const BlockVectorType    &b) const
{
  Assert (dst.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(dst.n_blocks(), n_block_rows()));
  Assert (b.n_blocks() == n_block_rows(),
          ExcDimensionMismatch(b.n_blocks(), n_block_rows()));
  Assert (x.n_blocks() == n_block_cols(),
          ExcDimensionMismatch(x.n_blocks(), n_block_cols()));
  // in block notation, the residual is
  // r_i = b_i - \sum_j A_ij x_j.
  // this can be written as
  // r_i = b_i - A_i0 x_0 - \sum_{j>0} A_ij x_j.
  //
  // for the first two terms, we can
  // call the residual function of
  // A_i0. for the other terms, we
  // use vmult_add. however, we want
  // to subtract, so in order to
  // avoid a temporary vector, we
  // perform a sign change of the
  // first two term before, and after
  // adding up
  for (unsigned int row=0; row<n_block_rows(); ++row)
    {
      block(row,0).residual (dst.block(row),
                             x.block(0),
                             b.block(row));

      for (size_type i=0; i<dst.block(row).size(); ++i)
        dst.block(row)(i) = -dst.block(row)(i);

      for (unsigned int col=1; col<n_block_cols(); ++col)
        block(row,col).vmult_add (dst.block(row),
                                  x.block(col));

      for (size_type i=0; i<dst.block(row).size(); ++i)
        dst.block(row)(i) = -dst.block(row)(i);
    };

  value_type res = 0;
  for (size_type row=0; row<n_block_rows(); ++row)
    res += dst.block(row).norm_sqr ();
  return std::sqrt(res);
}



template <class MatrixType>
inline
void
BlockMatrixBase<MatrixType>::print (std::ostream &out,
                                    const bool    alternative_output) const
{
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      {
        if (!alternative_output)
          out << "Block (" << row << ", " << col << ")" << std::endl;

        block(row, col).print(out, alternative_output);
      }
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin () const
{
  return const_iterator(this, 0);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end () const
{
  return const_iterator(this, m());
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::begin (const size_type r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::const_iterator
BlockMatrixBase<MatrixType>::end (const size_type r) const
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return const_iterator(this, r+1);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::begin ()
{
  return iterator(this, 0);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::end ()
{
  return iterator(this, m());
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::begin (const size_type r)
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return iterator(this, r);
}



template <class MatrixType>
inline
typename BlockMatrixBase<MatrixType>::iterator
BlockMatrixBase<MatrixType>::end (const size_type r)
{
  Assert (r<m(), ExcIndexRange(r,0,m()));
  return iterator(this, r+1);
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::collect_sizes ()
{
  std::vector<size_type> row_sizes (this->n_block_rows());
  std::vector<size_type> col_sizes (this->n_block_cols());

  // first find out the row sizes
  // from the first block column
  for (unsigned int r=0; r<this->n_block_rows(); ++r)
    row_sizes[r] = sub_objects[r][0]->m();
  // then check that the following
  // block columns have the same
  // sizes
  for (unsigned int c=1; c<this->n_block_cols(); ++c)
    for (unsigned int r=0; r<this->n_block_rows(); ++r)
      Assert (row_sizes[r] == sub_objects[r][c]->m(),
              ExcIncompatibleRowNumbers (r,0,r,c));

  // finally initialize the row
  // indices with this array
  this->row_block_indices.reinit (row_sizes);


  // then do the same with the columns
  for (unsigned int c=0; c<this->n_block_cols(); ++c)
    col_sizes[c] = sub_objects[0][c]->n();
  for (unsigned int r=1; r<this->n_block_rows(); ++r)
    for (unsigned int c=0; c<this->n_block_cols(); ++c)
      Assert (col_sizes[c] == sub_objects[r][c]->n(),
              ExcIncompatibleRowNumbers (0,c,r,c));

  // finally initialize the row
  // indices with this array
  this->column_block_indices.reinit (col_sizes);
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::prepare_add_operation ()
{
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      block(row, col).prepare_add();
}



template <class MatrixType>
void
BlockMatrixBase<MatrixType>::prepare_set_operation ()
{
  for (unsigned int row=0; row<n_block_rows(); ++row)
    for (unsigned int col=0; col<n_block_cols(); ++col)
      block(row, col).prepare_set();
}

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif    // __deal2__block_matrix_base_h
