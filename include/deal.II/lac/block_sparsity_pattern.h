// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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

#ifndef __deal2__block_sparsity_pattern_h
#define __deal2__block_sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/table.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/trilinos_sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/block_indices.h>

DEAL_II_NAMESPACE_OPEN


template <typename number> class BlockSparseMatrix;
class BlockSparsityPattern;
class BlockCompressedSparsityPattern;
class BlockCompressedSimpleSparsityPattern;
class BlockCompressedSetSparsityPattern;
#ifdef DEAL_II_WITH_TRILINOS
namespace TrilinosWrappers
{
  class BlockSparsityPattern;
}
#endif

/*! @addtogroup Sparsity
 *@{
 */


/**
 * This is the base class for block versions of the sparsity pattern and
 * compressed sparsity pattern classes. It has not much functionality, but
 * only administrates an array of sparsity pattern objects and delegates work
 * to them. It has mostly the same interface as has the SparsityPattern,
 * CompressedSparsityPattern, and CompressedSetSparsityPattern classes, and
 * simply transforms calls to its member functions to calls to the respective
 * member functions of the member sparsity patterns.
 *
 * The largest difference between the SparsityPattern and
 * CompressedSparsityPattern classes and this class is that
 * mostly, the matrices have different properties and you will want to
 * work on the blocks making up the matrix rather than the whole
 * matrix. You can access the different blocks using the
 * <tt>block(row,col)</tt> function.
 *
 * Attention: this object is not automatically notified if the size of
 * one of its subobjects' size is changed. After you initialize the
 * sizes of the subobjects, you will therefore have to call the
 * <tt>collect_sizes()</tt> function of this class! Note that, of course, all
 * sub-matrices in a (block-)row have to have the same number of rows, and
 * that all sub-matrices in a (block-)column have to have the same number of
 * columns.
 *
 * You will in general not want to use this class, but one of the
 * derived classes.
 *
 * @todo Handle optimization of diagonal elements of the underlying
 * SparsityPattern correctly.
 *
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Wolfgang Bangerth, 2000, 2001
 */
template <class SparsityPatternBase>
class BlockSparsityPatternBase : public Subscriptor
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Define a value which is used
   * to indicate that a certain
   * value in the @p colnums array
   * is unused, i.e. does not
   * represent a certain column
   * number index.
   *
   * This value is only an alias to
   * the respective value of the
   * SparsityPattern class.
   */
  static const size_type invalid_entry = SparsityPattern::invalid_entry;

  /**
   * Initialize the matrix empty,
   * that is with no memory
   * allocated. This is useful if
   * you want such objects as
   * member variables in other
   * classes. You can make the
   * structure usable by calling
   * the reinit() function.
   */
  BlockSparsityPatternBase ();

  /**
   * Initialize the matrix with the
   * given number of block rows and
   * columns. The blocks themselves
   * are still empty, and you have
   * to call collect_sizes() after
   * you assign them sizes.
   */
  BlockSparsityPatternBase (const size_type n_block_rows,
                            const size_type n_block_columns);

  /**
   * Copy constructor. This
   * constructor is only allowed to
   * be called if the sparsity pattern to be
   * copied is empty, i.e. there
   * are no block allocated at
   * present. This is for the same
   * reason as for the
   * SparsityPattern, see there
   * for the details.
   */
  BlockSparsityPatternBase (const BlockSparsityPatternBase &bsp);

  /**
   * Destructor.
   */
  ~BlockSparsityPatternBase ();

  /**
   * Resize the matrix, by setting
   * the number of block rows and
   * columns. This deletes all
   * blocks and replaces them by
   * unitialized ones, i.e. ones
   * for which also the sizes are
   * not yet set. You have to do
   * that by calling the reinit()
   * functions of the blocks
   * themselves. Do not forget to
   * call collect_sizes() after
   * that on this object.
   *
   * The reason that you have to
   * set sizes of the blocks
   * yourself is that the sizes may
   * be varying, the maximum number
   * of elements per row may be
   * varying, etc. It is simpler
   * not to reproduce the interface
   * of the SparsityPattern
   * class here but rather let the
   * user call whatever function
   * she desires.
   */
  void reinit (const size_type n_block_rows,
               const size_type n_block_columns);

  /**
   * Copy operator. For this the
   * same holds as for the copy
   * constructor: it is declared,
   * defined and fine to be called,
   * but the latter only for empty
   * objects.
   */
  BlockSparsityPatternBase &operator = (const BlockSparsityPatternBase &);

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
   */
  void collect_sizes ();

  /**
   * Access the block with the
   * given coordinates.
   */
  SparsityPatternBase &
  block (const size_type row,
         const size_type column);


  /**
   * Access the block with the
   * given coordinates. Version for
   * constant objects.
   */
  const SparsityPatternBase &
  block (const size_type row,
         const size_type column) const;

  /**
   * Grant access to the object
   * describing the distribution of
   * row indices to the individual
   * blocks.
   */
  const BlockIndices &
  get_row_indices () const;

  /**
   * Grant access to the object
   * describing the distribution of
   * column indices to the individual
   * blocks.
   */
  const BlockIndices &
  get_column_indices () const;

  /**
   * This function compresses the
   * sparsity structures that this
   * object represents. It simply
   * calls @p compress for all
   * sub-objects.
   */
  void compress ();

  /**
   * Return the number of blocks in a
   * column.
   */
  size_type n_block_rows () const;

  /**
   * Return the number of blocks in a
   * row.
   */
  size_type n_block_cols () const;

  /**
   * Return whether the object is
   * empty. It is empty if no
   * memory is allocated, which is
   * the same as that both
   * dimensions are zero. This
   * function is just the
   * concatenation of the
   * respective call to all
   * sub-matrices.
   */
  bool empty () const;

  /**
   * Return the maximum number of
   * entries per row. It returns
   * the maximal number of entries
   * per row accumulated over all
   * blocks in a row, and the
   * maximum over all rows.
   */
  size_type max_entries_per_row () const;

  /**
   * Add a nonzero entry to the matrix.
   * This function may only be called
   * for non-compressed sparsity patterns.
   *
   * If the entry already exists, nothing
   * bad happens.
   *
   * This function simply finds out
   * to which block <tt>(i,j)</tt> belongs
   * and then relays to that block.
   */
  void add (const size_type i, const size_type j);

  /**
   * Add several nonzero entries to the
   * specified matrix row.  This function
   * may only be called for
   * non-compressed sparsity patterns.
   *
   * If some of the entries already
   * exist, nothing bad happens.
   *
   * This function simply finds out to
   * which blocks <tt>(row,col)</tt> for
   * <tt>col</tt> in the iterator range
   * belong and then relays to those
   * blocks.
   */
  template <typename ForwardIterator>
  void add_entries (const size_type  row,
                    ForwardIterator  begin,
                    ForwardIterator  end,
                    const bool       indices_are_sorted = false);

  /**
   * Return number of rows of this
   * matrix, which equals the
   * dimension of the image
   * space. It is the sum of rows
   * of the (block-)rows of
   * sub-matrices.
   */
  size_type n_rows () const;

  /**
   * Return number of columns of
   * this matrix, which equals the
   * dimension of the range
   * space. It is the sum of
   * columns of the (block-)columns
   * of sub-matrices.
   */
  size_type n_cols () const;

  /**
   * Check if a value at a certain
   * position may be non-zero.
   */
  bool exists (const size_type i, const size_type j) const;

  /**
   * Number of entries in a
   * specific row, added up over
   * all the blocks that form this
   * row.
   */
  unsigned int row_length (const size_type row) const;

  /**
   * Return the number of nonzero
   * elements of this
   * matrix. Actually, it returns
   * the number of entries in the
   * sparsity pattern; if any of
   * the entries should happen to
   * be zero, it is counted anyway.
   *
   * This function may only be
   * called if the matrix struct is
   * compressed. It does not make
   * too much sense otherwise
   * anyway.
   *
   * In the present context, it is
   * the sum of the values as
   * returned by the sub-objects.
   */
  size_type n_nonzero_elements () const;

  /**
   * Print the sparsity of the
   * matrix. The output consists of
   * one line per row of the format
   * <tt>[i,j1,j2,j3,...]</tt>. <i>i</i>
   * is the row number and
   * <i>jn</i> are the allocated
   * columns in this row.
   */
  void print (std::ostream &out) const;

  /**
   * Print the sparsity of the matrix
   * in a format that <tt>gnuplot</tt> understands
   * and which can be used to plot the
   * sparsity pattern in a graphical
   * way. This is the same functionality
   * implemented for usual sparsity
   * patterns, see @ref SparsityPattern.
   */
  void print_gnuplot (std::ostream &out) const;

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
  /**
   * Exception
   */
  DeclException0 (ExcInvalidConstructorCall);
  //@}

protected:

  /**
   * Number of block rows.
   */
  size_type rows;

  /**
   * Number of block columns.
   */
  size_type columns;

  /**
   * Array of sparsity patterns.
   */
  Table<2,SmartPointer<SparsityPatternBase, BlockSparsityPatternBase<SparsityPatternBase> > > sub_objects;

  /**
   * Object storing and managing
   * the transformation of row
   * indices to indices of the
   * sub-objects.
   */
  BlockIndices    row_indices;

  /**
   * Object storing and managing
   * the transformation of column
   * indices to indices of the
   * sub-objects.
   */
  BlockIndices    column_indices;

private:
  /**
   * Temporary vector for counting the
   * elements written into the
   * individual blocks when doing a
   * collective add or set.
   */
  std::vector<size_type > counter_within_block;

  /**
   * Temporary vector for column
   * indices on each block when writing
   * local to global data on each
   * sparse matrix.
   */
  std::vector<std::vector<size_type > > block_column_indices;

  /**
   * Make the block sparse matrix a
   * friend, so that it can use our
   * #row_indices and
   * #column_indices objects.
   */
  template <typename number>
  friend class BlockSparseMatrix;
};



/**
 * This class extends the base class to implement an array of sparsity
 * patterns that can be used by block sparse matrix objects. It only
 * adds a few additional member functions, but the main interface
 * stems from the base class, see there for more information.
 *
 * This class is an example of the "static" type of @ref Sparsity.
 *
 * @author Wolfgang Bangerth, 2000, 2001
 */
class BlockSparsityPattern : public BlockSparsityPatternBase<SparsityPattern>
{
public:

  /**
   * Initialize the matrix empty,
   * that is with no memory
   * allocated. This is useful if
   * you want such objects as
   * member variables in other
   * classes. You can make the
   * structure usable by calling
   * the reinit() function.
   */
  BlockSparsityPattern ();

  /**
   * Initialize the matrix with the
   * given number of block rows and
   * columns. The blocks themselves
   * are still empty, and you have
   * to call collect_sizes() after
   * you assign them sizes.
   */
  BlockSparsityPattern (const size_type n_rows,
                        const size_type n_columns);

  /**
   * Forwarding to
   * BlockSparsityPatternBase::reinit().
   */
  void reinit (const size_type n_block_rows,
               const size_type n_block_columns);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns as well as a row
   * length vector.
   *
   * The row length vector should
   * be in the format produced by
   * DoFTools. Alternatively, there
   * is a simplified version,
   * where each of the inner
   * vectors has length one. Then,
   * the corresponding entry is
   * used as the maximal row length.
   *
   * For the diagonal blocks, the
   * inner SparsityPattern is
   * initialized with optimized
   * diagonals, while this is not
   * done for the off-diagonal blocks.
   */
  void reinit (const BlockIndices &row_indices,
               const BlockIndices &col_indices,
               const std::vector<std::vector<unsigned int> > &row_lengths);


  /**
   * Return whether the structure
   * is compressed or not,
   * i.e. whether all sub-matrices
   * are compressed.
   */
  bool is_compressed () const;

  /**
   * Determine an estimate for the
   * memory consumption (in bytes)
   * of this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Copy data from an object of
   * type
   * BlockCompressedSparsityPattern,
   * i.e. resize this object to the
   * size of the given argument,
   * and copy over the contents of
   * each of the
   * subobjects. Previous content
   * of this object is lost.
   */
  void copy_from (const BlockCompressedSparsityPattern &csp);

  /**
   * Copy data from an object of
   * type
   * BlockCompressedSetSparsityPattern,
   * i.e. resize this object to the
   * size of the given argument,
   * and copy over the contents of
   * each of the
   * subobjects. Previous content
   * of this object is lost.
   */
  void copy_from (const BlockCompressedSetSparsityPattern &csp);

  /**
   * Copy data from an object of
   * type
   * BlockCompressedSimpleSparsityPattern,
   * i.e. resize this object to the
   * size of the given argument,
   * and copy over the contents of
   * each of the
   * subobjects. Previous content
   * of this object is lost.
   */
  void copy_from (const BlockCompressedSimpleSparsityPattern &csp);



};



/**
 * This class extends the base class to implement an array of compressed
 * sparsity patterns that can be used to initialize objects of type
 * BlockSparsityPattern. It does not add additional member functions, but
 * rather acts as a @p typedef to introduce the name of this class, without
 * requiring the user to specify the templated name of the base class. For
 * information on the interface of this class refer to the base class. The
 * individual blocks are based on the CompressedSparsityPattern class.
 *
 * This class is an example of the "dynamic" type of @ref Sparsity.
 *
 * <b>Note:</b>
 * There are several, exchangeable variations of this class, see @ref Sparsity,
 * section 'Dynamic block sparsity patterns' for more information.
 *
 * <b>Note:</b> This class used to be called
 * CompressedBlockSparsityPattern. However, since it's a block wrapper around
 * the CompressedSparsityPattern class, this is a misnomer and the class has
 * been renamed.
 *
 * <h3>Example</h3>
 *
 * Usage of this class is very similar to CompressedSparsityPattern,
 * but since the use of block indices causes some additional
 * complications, we give a short example.
 *
 * @dontinclude compressed_block_sparsity_pattern.cc
 *
 * After the the DoFHandler <tt>dof</tt> and the ConstraintMatrix
 * <tt>constraints</tt> have been set up with a system element, we
 * must count the degrees of freedom in each matrix block:
 *
 * @skipline dofs_per_block
 * @until count
 *
 * Now, we are ready to set up the BlockCompressedSparsityPattern.
 *
 * @until collect
 *
 * It is filled as if it were a normal pattern
 *
 * @until condense
 *
 * In the end, it is copied to a normal BlockSparsityPattern for later
 * use.
 *
 * @until copy
 *
 * @author Wolfgang Bangerth, 2000, 2001, Guido Kanschat, 2006, 2007
 */
class BlockCompressedSparsityPattern : public BlockSparsityPatternBase<CompressedSparsityPattern>
{
public:

  /**
   * Initialize the matrix empty,
   * that is with no memory
   * allocated. This is useful if
   * you want such objects as
   * member variables in other
   * classes. You can make the
   * structure usable by calling
   * the reinit() function.
   */
  BlockCompressedSparsityPattern ();

  /**
   * Initialize the matrix with the
   * given number of block rows and
   * columns. The blocks themselves
   * are still empty, and you have
   * to call collect_sizes() after
   * you assign them sizes.
   */
  BlockCompressedSparsityPattern (const size_type n_rows,
                                  const size_type n_columns);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns. This function is
   * equivalent to calling the
   * previous constructor with the
   * length of the two index vector
   * and then entering the index
   * values.
   */
  BlockCompressedSparsityPattern (const std::vector<size_type> &row_block_sizes,
                                  const std::vector<size_type> &col_block_sizes);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns.
   */
  BlockCompressedSparsityPattern (const BlockIndices &row_indices,
                                  const BlockIndices &col_indices);

  /**
   * Resize the matrix to a tensor
   * product of matrices with
   * dimensions defined by the
   * arguments.
   *
   * The matrix will have as many
   * block rows and columns as
   * there are entries in the two
   * arguments. The block at
   * position (<i>i,j</i>) will
   * have the dimensions
   * <tt>row_block_sizes[i]</tt>
   * times <tt>col_block_sizes[j]</tt>.
   */
  void reinit (const std::vector<size_type> &row_block_sizes,
               const std::vector<size_type> &col_block_sizes);

  /**
   * Resize the matrix to a tensor
   * product of matrices with
   * dimensions defined by the
   * arguments. The two
   * BlockIndices objects must be
   * initialized and the sparsity
   * pattern will have the
   * same block structure afterwards.
   */
  void reinit (const BlockIndices &row_indices, const BlockIndices &col_indices);

  /**
   * Allow the use of the reinit
   * functions of the base class as
   * well.
   */
  using BlockSparsityPatternBase<CompressedSparsityPattern>::reinit;
};


/**
 * A typdef to preserve the old name CompressedBlockSparsityPattern even after
 * we changed the naming of the class to BlockCompressedSparsityPattern. The
 * old name is now deprecated and user codes should use the name
 * BlockCompressedSparsityPattern instead.
 *
 * @deprecated
 */
typedef BlockCompressedSparsityPattern CompressedBlockSparsityPattern DEAL_II_DEPRECATED;



/**
 * This class extends the base class to implement an array of compressed
 * sparsity patterns that can be used to initialize objects of type
 * BlockSparsityPattern. It is used in the same way as the
 * BlockCompressedSparsityPattern except that it builds upon the
 * CompressedSetSparsityPattern instead of the CompressedSparsityPattern. See
 * the documentation of the BlockCompressedSparsityPattern for examples.
 *
 * This class is an example of the "dynamic" type of @ref Sparsity.
 *
 * <b>Note:</b>
 * There are several, exchangeable variations of this class, see @ref Sparsity,
 * section 'Dynamic block sparsity patterns' for more information.
 *
 * @author Wolfgang Bangerth, 2007
 */
class BlockCompressedSetSparsityPattern : public BlockSparsityPatternBase<CompressedSetSparsityPattern>
{
public:

  /**
   * Initialize the matrix empty,
   * that is with no memory
   * allocated. This is useful if
   * you want such objects as
   * member variables in other
   * classes. You can make the
   * structure usable by calling
   * the reinit() function.
   */
  BlockCompressedSetSparsityPattern ();

  /**
   * Initialize the matrix with the
   * given number of block rows and
   * columns. The blocks themselves
   * are still empty, and you have
   * to call collect_sizes() after
   * you assign them sizes.
   */
  BlockCompressedSetSparsityPattern (const size_type n_rows,
                                     const size_type n_columns);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns. This function is
   * equivalent to calling the
   * previous constructor with the
   * length of the two index vector
   * and then entering the index
   * values.
   */
  BlockCompressedSetSparsityPattern (const std::vector<size_type> &row_block_sizes,
                                     const std::vector<size_type> &col_block_sizes);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns.
   */
  BlockCompressedSetSparsityPattern (const BlockIndices &row_indices,
                                     const BlockIndices &col_indices);

  /**
   * Resize the matrix to a tensor
   * product of matrices with
   * dimensions defined by the
   * arguments.
   *
   * The matrix will have as many
   * block rows and columns as
   * there are entries in the two
   * arguments. The block at
   * position (<i>i,j</i>) will
   * have the dimensions
   * <tt>row_block_sizes[i]</tt>
   * times <tt>col_block_sizes[j]</tt>.
   */
  void reinit (const std::vector<size_type> &row_block_sizes,
               const std::vector<size_type> &col_block_sizes);

  /**
   * Resize the matrix to a tensor
   * product of matrices with
   * dimensions defined by the
   * arguments. The two
   * BlockIndices objects must be
   * initialized and the sparsity
   * pattern will have the
   * same block structure afterwards.
   */
  void reinit (const BlockIndices &row_indices, const BlockIndices &col_indices);

  /**
   * Allow the use of the reinit
   * functions of the base class as
   * well.
   */
  using BlockSparsityPatternBase<CompressedSetSparsityPattern>::reinit;
};





/**
 * This class extends the base class to implement an array of compressed
 * sparsity patterns that can be used to initialize objects of type
 * BlockSparsityPattern. It is used in the same way as the
 * BlockCompressedSparsityPattern except that it builds upon the
 * CompressedSimpleSparsityPattern instead of the CompressedSparsityPattern. See
 * the documentation of the BlockCompressedSparsityPattern for examples.
 *
 * This class is an example of the "dynamic" type of @ref Sparsity.
 *
 * <b>Note:</b>
 * There are several, exchangeable variations of this class, see @ref Sparsity,
 * section 'Dynamic block sparsity patterns' for more information.
 *
 * This class is used in step-22 and step-31.
 *
 * @author Timo Heister, 2008
 */
class BlockCompressedSimpleSparsityPattern : public BlockSparsityPatternBase<CompressedSimpleSparsityPattern>
{
public:

  /**
   * Initialize the matrix empty,
   * that is with no memory
   * allocated. This is useful if
   * you want such objects as
   * member variables in other
   * classes. You can make the
   * structure usable by calling
   * the reinit() function.
   */
  BlockCompressedSimpleSparsityPattern ();

  /**
   * Initialize the matrix with the
   * given number of block rows and
   * columns. The blocks themselves
   * are still empty, and you have
   * to call collect_sizes() after
   * you assign them sizes.
   */
  BlockCompressedSimpleSparsityPattern (const size_type n_rows,
                                        const size_type n_columns);

  /**
   * Initialize the pattern with
   * two BlockIndices for the block
   * structures of matrix rows and
   * columns. This function is
   * equivalent to calling the
   * previous constructor with the
   * length of the two index vector
   * and then entering the index
   * values.
   */
  BlockCompressedSimpleSparsityPattern (const std::vector<size_type> &row_block_sizes,
                                        const std::vector<size_type> &col_block_sizes);

  /**
   * Initialize the pattern with symmetric
   * blocks. The number of IndexSets in the
   * vector determine the number of rows
   * and columns of blocks. The size of
   * each block is determined by the size()
   * of the respective IndexSet. Each block
   * only stores the rows given by the
   * values in the IndexSet, which is
   * useful for distributed memory parallel
   * computations and usually corresponds
   * to the locally owned DoFs.
   */
  BlockCompressedSimpleSparsityPattern (const std::vector<IndexSet> &partitioning);

  /**
   * Resize the pattern to a tensor product
   * of matrices with dimensions defined by
   * the arguments.
   *
   * The matrix will have as many
   * block rows and columns as
   * there are entries in the two
   * arguments. The block at
   * position (<i>i,j</i>) will
   * have the dimensions
   * <tt>row_block_sizes[i]</tt>
   * times <tt>col_block_sizes[j]</tt>.
   */
  void reinit (const std::vector<size_type> &row_block_sizes,
               const std::vector<size_type> &col_block_sizes);

  /**
   * Resize the pattern with symmetric
   * blocks determined by the size() of
   * each IndexSet. See the constructor
   * taking a vector of IndexSets for
   * details.
   */
  void reinit(const std::vector<IndexSet> &partitioning);

  /**
   * Access to column number field.
   * Return the column number of
   * the @p index th entry in row @p row.
   */
  size_type column_number (const size_type row,
                           const unsigned int index) const;

  /**
   * Allow the use of the reinit
   * functions of the base class as
   * well.
   */
  using BlockSparsityPatternBase<CompressedSimpleSparsityPattern>::reinit;
};




#ifdef DEAL_II_WITH_TRILINOS


/**
 * This class extends the base class to implement an array of Trilinos
 * sparsity patterns that can be used to initialize Trilinos block sparse
 * matrices that can be distributed among different processors. It is used
 * in the same way as the BlockSparsityPattern except that it builds upon
 * the TrilinosWrappers::SparsityPattern instead of the
 * dealii::SparsityPattern. See the documentation of the
 * BlockSparsityPattern for examples.
 *
 * This class is has properties of the "dynamic" type of @ref Sparsity (in
 * the sense that it can extend the memory if too little elements were
 * allocated), but otherwise is more like the basic deal.II SparsityPattern
 * (in the sense that the method compress() needs to be called before the
 * pattern can be used).
 *
 * This class is used in step-32.
 *
 * @author Martin Kronbichler, 2008, 2009
 */
namespace TrilinosWrappers
{
  class BlockSparsityPattern :
    public dealii::BlockSparsityPatternBase<SparsityPattern>
  {
  public:

    /**
     * Initialize the matrix empty, that is with no memory allocated. This is
     * useful if you want such objects as member variables in other
     * classes. You can make the structure usable by calling the reinit()
     * function.
     */
    BlockSparsityPattern ();

    /**
     * Initialize the matrix with the given number of block rows and
     * columns. The blocks themselves are still empty, and you have to call
     * collect_sizes() after you assign them sizes.
     */
    BlockSparsityPattern (const size_type n_rows,
                          const size_type n_columns);

    /**
     * Initialize the pattern with two BlockIndices for the block structures
     * of matrix rows and columns. This function is equivalent to calling the
     * previous constructor with the length of the two index vector and then
     * entering the index values.
     */
    BlockSparsityPattern (const std::vector<size_type> &row_block_sizes,
                          const std::vector<size_type> &col_block_sizes);

    /**
     * Initialize the pattern with an array Epetra_Map that specifies both
     * rows and columns of the matrix (so the final matrix will be a square
     * matrix), where the Epetra_Map specifies the parallel distribution of
     * the degrees of freedom on the individual block.  This function is
     * equivalent to calling the second constructor with the length of the
     * mapping vector and then entering the index values.
     */
    BlockSparsityPattern (const std::vector<Epetra_Map> &parallel_partitioning);

    /**
     * Initialize the pattern with an array of index sets that specifies both
     * rows and columns of the matrix (so the final matrix will be a square
     * matrix), where the size() of the IndexSets specifies the size of the
     * blocks and the values in each IndexSet denotes the rows that are going
     * to be saved in each block.
     */
    BlockSparsityPattern (const std::vector<IndexSet> &parallel_partitioning,
                          const MPI_Comm &communicator = MPI_COMM_WORLD);

    /**
     * Initialize the pattern with two arrays of index sets that specify rows
     * and columns of the matrix, where the size() of the IndexSets specifies
     * the size of the blocks and the values in each IndexSet denotes the rows
     * that are going to be saved in each block. The additional index set
     * writable_rows is used to set all rows that we allow to write
     * locally. This constructor is used to create matrices that allow several
     * threads to write simultaneously into the matrix (to different rows, of
     * course), see the method TrilinosWrappers::SparsityPattern::reinit
     * method with three index set arguments for more details.
     */
    BlockSparsityPattern (const std::vector<IndexSet> &row_parallel_partitioning,
                          const std::vector<IndexSet> &column_parallel_partitioning,
                          const std::vector<IndexSet> &writeable_rows,
                          const MPI_Comm              &communicator = MPI_COMM_WORLD);

    /**
     * Resize the matrix to a tensor product of matrices with dimensions
     * defined by the arguments.
     *
     * The matrix will have as many block rows and columns as there are
     * entries in the two arguments. The block at position (<i>i,j</i>) will
     * have the dimensions <tt>row_block_sizes[i]</tt> times
     * <tt>col_block_sizes[j]</tt>.
     */
    void reinit (const std::vector<size_type> &row_block_sizes,
                 const std::vector<size_type> &col_block_sizes);

    /**
     * Resize the matrix to a square tensor product of matrices with parallel
     * distribution according to the specifications in the array of
     * Epetra_Maps.
     */
    void reinit (const std::vector<Epetra_Map> &parallel_partitioning);

    /**
     * Resize the matrix to a square tensor product of matrices. See the
     * constructor that takes a vector of IndexSets for details.
     */
    void reinit (const std::vector<IndexSet> &parallel_partitioning,
                 const MPI_Comm              &communicator = MPI_COMM_WORLD);

    /**
     * Resize the matrix to a rectangular block matrices. This method allows
     * rows and columns to be different, both in the outer block structure and
     * within the blocks.
     */
    void reinit (const std::vector<IndexSet> &row_parallel_partitioning,
                 const std::vector<IndexSet> &column_parallel_partitioning,
                 const MPI_Comm              &communicator = MPI_COMM_WORLD);

    /**
     * Resize the matrix to a rectangular block matrices that furthermore
     * explicitly specify the writable rows in each of the blocks. This method
     * is used to create matrices that allow several threads to write
     * simultaneously into the matrix (to different rows, of course), see the
     * method TrilinosWrappers::SparsityPattern::reinit method with three
     * index set arguments for more details.
     */
    void reinit (const std::vector<IndexSet> &row_parallel_partitioning,
                 const std::vector<IndexSet> &column_parallel_partitioning,
                 const std::vector<IndexSet> &writeable_rows,
                 const MPI_Comm              &communicator = MPI_COMM_WORLD);

    /**
     * Allow the use of the reinit functions of the base class as well.
     */
    using BlockSparsityPatternBase<SparsityPattern>::reinit;
  };
}

#endif


/*@}*/
/*---------------------- Template functions -----------------------------------*/



template <class SparsityPatternBase>
inline
SparsityPatternBase &
BlockSparsityPatternBase<SparsityPatternBase>::block (const size_type row,
                                                      const size_type column)
{
  Assert (row<rows, ExcIndexRange(row,0,rows));
  Assert (column<columns, ExcIndexRange(column,0,columns));
  return *sub_objects[row][column];
}



template <class SparsityPatternBase>
inline
const SparsityPatternBase &
BlockSparsityPatternBase<SparsityPatternBase>::block (const size_type row,
                                                      const size_type column) const
{
  Assert (row<rows, ExcIndexRange(row,0,rows));
  Assert (column<columns, ExcIndexRange(column,0,columns));
  return *sub_objects[row][column];
}



template <class SparsityPatternBase>
inline
const BlockIndices &
BlockSparsityPatternBase<SparsityPatternBase>::get_row_indices () const
{
  return row_indices;
}



template <class SparsityPatternBase>
inline
const BlockIndices &
BlockSparsityPatternBase<SparsityPatternBase>::get_column_indices () const
{
  return column_indices;
}



template <class SparsityPatternBase>
inline
void
BlockSparsityPatternBase<SparsityPatternBase>::add (const size_type i,
                                                    const size_type j)
{
  // if you get an error here, are
  // you sure you called
  // <tt>collect_sizes()</tt> before?
  const std::pair<size_type,size_type>
  row_index = row_indices.global_to_local (i),
  col_index = column_indices.global_to_local (j);
  sub_objects[row_index.first][col_index.first]->add (row_index.second,
                                                      col_index.second);
}



template <class SparsityPatternBase>
template <typename ForwardIterator>
void
BlockSparsityPatternBase<SparsityPatternBase>::add_entries (const size_type row,
                                                            ForwardIterator begin,
                                                            ForwardIterator end,
                                                            const bool      indices_are_sorted)
{
  // Resize scratch arrays
  if (block_column_indices.size() < this->n_block_cols())
    {
      block_column_indices.resize (this->n_block_cols());
      counter_within_block.resize (this->n_block_cols());
    }

  const size_type n_cols = static_cast<size_type>(end - begin);

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
  if (block_column_indices[0].size() < n_cols)
    for (size_type i=0; i<this->n_block_cols(); ++i)
      block_column_indices[i].resize(n_cols);

  // Reset the number of added elements
  // in each block to zero.
  for (size_type i=0; i<this->n_block_cols(); ++i)
    counter_within_block[i] = 0;

  // Go through the column indices to
  // find out which portions of the
  // values should be set in which
  // block of the matrix. We need to
  // touch all the data, since we can't
  // be sure that the data of one block
  // is stored contiguously (in fact,
  // indices will be intermixed when it
  // comes from an element matrix).
  for (ForwardIterator it = begin; it != end; ++it)
    {
      const size_type col = *it;

      const std::pair<size_type , size_type>
      col_index = this->column_indices.global_to_local(col);

      const size_type local_index = counter_within_block[col_index.first]++;

      block_column_indices[col_index.first][local_index] = col_index.second;
    }

#ifdef DEBUG
  // If in debug mode, do a check whether
  // the right length has been obtained.
  size_type length = 0;
  for (size_type i=0; i<this->n_block_cols(); ++i)
    length += counter_within_block[i];
  Assert (length == n_cols, ExcInternalError());
#endif

  // Now we found out about where the
  // individual columns should start and
  // where we should start reading out
  // data. Now let's write the data into
  // the individual blocks!
  const std::pair<size_type , size_type>
  row_index = this->row_indices.global_to_local (row);
  for (size_type block_col=0; block_col<n_block_cols(); ++block_col)
    {
      if (counter_within_block[block_col] == 0)
        continue;
      sub_objects[row_index.first][block_col]->
      add_entries (row_index.second,
                   block_column_indices[block_col].begin(),
                   block_column_indices[block_col].begin()+counter_within_block[block_col],
                   indices_are_sorted);
    }
}



template <class SparsityPatternBase>
inline
bool
BlockSparsityPatternBase<SparsityPatternBase>::exists (const size_type i,
                                                       const size_type j) const
{
  // if you get an error here, are
  // you sure you called
  // <tt>collect_sizes()</tt> before?
  const std::pair<size_type , size_type>
  row_index = row_indices.global_to_local (i),
  col_index = column_indices.global_to_local (j);
  return sub_objects[row_index.first][col_index.first]->exists (row_index.second,
         col_index.second);
}



template <class SparsityPatternBase>
inline
unsigned int
BlockSparsityPatternBase<SparsityPatternBase>::
row_length (const size_type row) const
{
  const std::pair<size_type , size_type>
  row_index = row_indices.global_to_local (row);

  unsigned int c = 0;

  for (size_type b=0; b<rows; ++b)
    c += sub_objects[row_index.first][b]->row_length (row_index.second);

  return c;
}



template <class SparsityPatternBase>
inline
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::n_block_cols () const
{
  return columns;
}



template <class SparsityPatternBase>
inline
typename BlockSparsityPatternBase<SparsityPatternBase>::size_type
BlockSparsityPatternBase<SparsityPatternBase>::n_block_rows () const
{
  return rows;
}


inline
BlockCompressedSimpleSparsityPattern::size_type
BlockCompressedSimpleSparsityPattern::column_number (const size_type row,
                                                     const unsigned int index) const
{
  // .first= ith block, .second = jth row in that block
  const std::pair<size_type ,size_type >
  row_index = row_indices.global_to_local (row);

  Assert(index<row_length(row), ExcIndexRange(index, 0, row_length(row)));

  size_type c = 0;
  size_type block_columns = 0; //sum of n_cols for all blocks to the left
  for (unsigned int b=0; b<columns; ++b)
    {
      unsigned int rowlen = sub_objects[row_index.first][b]->row_length (row_index.second);
      if (index<c+rowlen)
        return block_columns+sub_objects[row_index.first][b]->column_number(row_index.second, index-c);
      c += rowlen;
      block_columns += sub_objects[row_index.first][b]->n_cols();
    }

  Assert(false, ExcInternalError());
  return 0;
}


inline
void
BlockSparsityPattern::reinit (
  const size_type n_block_rows,
  const size_type n_block_columns)
{
  BlockSparsityPatternBase<SparsityPattern>::reinit (
    n_block_rows, n_block_columns);
}


DEAL_II_NAMESPACE_CLOSE

#endif
