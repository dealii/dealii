//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_sparsity_pattern_h
#define __deal2__block_sparsity_pattern_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/table.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>
#include <lac/sparsity_pattern.h>
#include <lac/trilinos_sparsity_pattern.h>
#include <lac/compressed_sparsity_pattern.h>
#include <lac/compressed_set_sparsity_pattern.h>
#include <lac/compressed_simple_sparsity_pattern.h>
#include <lac/block_indices.h>

DEAL_II_NAMESPACE_OPEN


template <typename number> class BlockSparseMatrix;
class BlockSparsityPattern;
class BlockCompressedSparsityPattern;
class BlockCompressedSimpleSparsityPattern;
class BlockCompressedSetSparsityPattern;
#ifdef DEAL_II_USE_TRILINOS
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
 * @author Wolfgang Bangerth, 2000, 2001
 */
template <class SparsityPatternBase>
class BlockSparsityPatternBase : public Subscriptor
{
  public:
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
    static const unsigned int invalid_entry = SparsityPattern::invalid_entry;
    
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
    BlockSparsityPatternBase (const unsigned int n_block_rows,
			      const unsigned int n_block_columns);
    
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
    void reinit (const unsigned int n_block_rows,
		 const unsigned int n_block_columns);

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    BlockSparsityPatternBase & operator = (const BlockSparsityPatternBase &);

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
    block (const unsigned int row,
	   const unsigned int column);
    
    
				     /**
				      * Access the block with the
				      * given coordinates. Version for
				      * constant objects.
				      */
    const SparsityPatternBase &
    block (const unsigned int row,
	   const unsigned int column) const;    

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
    unsigned int n_block_rows () const;
    
				     /**
				      * Return the number of blocks in a
				      * row.
				      */
    unsigned int n_block_cols () const;
  
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
    unsigned int max_entries_per_row () const;

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
    void add (const unsigned int i, const unsigned int j);
    
				     /**
				      * Return number of rows of this
				      * matrix, which equals the
				      * dimension of the image
				      * space. It is the sum of rows
				      * of the (block-)rows of
				      * sub-matrices.
				      */
    unsigned int n_rows () const;

				     /**
				      * Return number of columns of
				      * this matrix, which equals the
				      * dimension of the range
				      * space. It is the sum of
				      * columns of the (block-)columns
				      * of sub-matrices.
				      */
    unsigned int n_cols () const;

				     /**
				      * Check if a value at a certain
				      * position may be non-zero.
				      */
    bool exists (const unsigned int i, const unsigned int j) const;

				     /**
				      * Number of entries in a
				      * specific row, added up over
				      * all the blocks that form this
				      * row.
				      */
    unsigned int row_length (const unsigned int row) const;
    
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
    unsigned int n_nonzero_elements () const;

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
    DeclException2 (ExcIncompatibleSizes,
		    int, int,
		    << "The number of blocks " << arg1 << " and " << arg2
		    << " are different.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCall);
				     //@}
    
  protected:

				     /**
				      * Number of block rows.
				      */
    unsigned int rows;

				     /**
				      * Number of block columns.
				      */
    unsigned int columns;
    
				     /**
				      * Array of sparsity patterns.
				      */
    Table<2,SmartPointer<SparsityPatternBase> > sub_objects;

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
    BlockSparsityPattern (const unsigned int n_rows,
			  const unsigned int n_columns);

				     /**
				      * Forwarding to
				      * BlockSparsityPatternBase::reinit().
				      */
    void reinit (const unsigned int n_block_rows,
		 const unsigned int n_block_columns);
    
    				     /**
				      * Initialize the pattern with
				      * two BlockIndices for the block
				      * structures of matrix rows and
				      * columns as well as a row
				      * length vector in the format
				      * produced by DoFTools.
				      */
    void reinit (const BlockIndices& row_indices,
		 const BlockIndices& col_indices,
		 const std::vector<std::vector<unsigned int> >& row_lengths);
    

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
    unsigned int memory_consumption () const;

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
    BlockCompressedSparsityPattern (const unsigned int n_rows,
				    const unsigned int n_columns);

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
    BlockCompressedSparsityPattern (const std::vector<unsigned int>& row_block_sizes,
				    const std::vector<unsigned int>& col_block_sizes);
    
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
    void reinit (const std::vector< unsigned int > &row_block_sizes,
		 const std::vector< unsigned int > &col_block_sizes);

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
typedef BlockCompressedSparsityPattern CompressedBlockSparsityPattern;



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
    BlockCompressedSetSparsityPattern (const unsigned int n_rows,
				    const unsigned int n_columns);

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
    BlockCompressedSetSparsityPattern (const std::vector<unsigned int>& row_block_sizes,
				    const std::vector<unsigned int>& col_block_sizes);
    
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
    void reinit (const std::vector< unsigned int > &row_block_sizes,
		 const std::vector< unsigned int > &col_block_sizes);

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
 * This class is used in @ref step_22 "step-22" and @ref step_31 "step-31".
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
    BlockCompressedSimpleSparsityPattern (const unsigned int n_rows,
					  const unsigned int n_columns);

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
    BlockCompressedSimpleSparsityPattern (const std::vector<unsigned int>& row_block_sizes,
				    const std::vector<unsigned int>& col_block_sizes);

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
    void reinit (const std::vector< unsigned int > &row_block_sizes,
		 const std::vector< unsigned int > &col_block_sizes);

				     /**
				      * Allow the use of the reinit
				      * functions of the base class as
				      * well.
				      */
    using BlockSparsityPatternBase<CompressedSimpleSparsityPattern>::reinit;
};




#ifdef DEAL_II_USE_TRILINOS


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
 * This class is used in @ref step_32 "step-32".
 *
 * @author Martin Kronbichler, 2008
 */
namespace TrilinosWrappers
{
  class BlockSparsityPattern : 
    public dealii::BlockSparsityPatternBase<SparsityPattern>
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
      BlockSparsityPattern (const unsigned int n_rows,
			    const unsigned int n_columns);

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
      BlockSparsityPattern (const std::vector<unsigned int>& row_block_sizes,
			    const std::vector<unsigned int>& col_block_sizes);

    				     /**
				      * Initialize the pattern with an array
				      * Epetra_Map that specifies both rows
				      * and columns of the matrix (so the
				      * final matrix will be a square
				      * matrix), where the Epetra_Map
				      * specifies the parallel distribution
				      * of the degrees of freedom on the
				      * individual block.  This function is
				      * equivalent to calling the second
				      * constructor with the length of the
				      * mapping vector and then entering the
				      * index values.
				      */
      BlockSparsityPattern (const std::vector<Epetra_Map>& input_maps);

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
      void reinit (const std::vector< unsigned int > &row_block_sizes,
		   const std::vector< unsigned int > &col_block_sizes);

				     /**
				      * Resize the matrix to a square tensor
				      * product of matrices with parallel
				      * distribution according to the
				      * specifications in the array of
				      * Epetra_Maps.
				      */
      void reinit (const std::vector<Epetra_Map>& input_maps);


				     /**
				      * Allow the use of the reinit
				      * functions of the base class as
				      * well.
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
BlockSparsityPatternBase<SparsityPatternBase>::block (const unsigned int row,
                                                      const unsigned int column)
{
  Assert (row<rows, ExcIndexRange(row,0,rows));
  Assert (column<columns, ExcIndexRange(column,0,columns));
  return *sub_objects[row][column];
}



template <class SparsityPatternBase>
inline
const SparsityPatternBase &
BlockSparsityPatternBase<SparsityPatternBase>::block (const unsigned int row,
                                                      const unsigned int column) const
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
BlockSparsityPatternBase<SparsityPatternBase>::add (const unsigned int i,
						    const unsigned int j)
{
				   // if you get an error here, are
				   // you sure you called
				   // <tt>collect_sizes()</tt> before?
  const std::pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (i),
    col_index = column_indices.global_to_local (j);
  sub_objects[row_index.first][col_index.first]->add (row_index.second,
						      col_index.second);
}



template <class SparsityPatternBase>
inline
bool
BlockSparsityPatternBase<SparsityPatternBase>::exists (const unsigned int i,
						       const unsigned int j) const
{
				   // if you get an error here, are
				   // you sure you called
				   // <tt>collect_sizes()</tt> before?
  const std::pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (i),
    col_index = column_indices.global_to_local (j);
  return sub_objects[row_index.first][col_index.first]->exists (row_index.second,
								col_index.second);
}



template <class SparsityPatternBase>
inline
unsigned int
BlockSparsityPatternBase<SparsityPatternBase>::
row_length (const unsigned int row) const
{
  const std::pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (row);

  unsigned int c = 0;
  
  for (unsigned int b=0; b<rows; ++b)
    c += sub_objects[row_index.first][b]->row_length (row_index.second);

  return c;
}



template <class SparsityPatternBase>
inline
unsigned int
BlockSparsityPatternBase<SparsityPatternBase>::n_block_cols () const
{
  return columns;
}



template <class SparsityPatternBase>
inline
unsigned int
BlockSparsityPatternBase<SparsityPatternBase>::n_block_rows () const
{
  return rows;
}


inline
void
BlockSparsityPattern::reinit (
  const unsigned int n_block_rows,
  const unsigned int n_block_columns)
{
  BlockSparsityPatternBase<SparsityPattern>::reinit (
    n_block_rows, n_block_columns);
}


DEAL_II_NAMESPACE_CLOSE

#endif
