//----------------------------  block_sparsity_pattern.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  block_sparsity_pattern.h  ---------------------------
#ifndef __deal2__block_sparsity_pattern_h
#define __deal2__block_sparsity_pattern_h


#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <lac/forward_declarations.h>
#include <lac/sparsity_pattern.h>
#include <lac/block_indices.h>


/**
 * This is a block version of the sparsity pattern class. It has not
 * much functionality, but only administrated an array of sparsity
 * pattern objects and delegates work to them. It has mostly the same
 * interface as has the #SparsityPattern# class, and simply transforms
 * calls to its member functions to calls to the respective member
 * functions of the member sparsity patterns.
 *
 * The largest difference between the #SparsityPattern# class and this
 * class is probably the absence of several of the constructors as
 * well as the absenceof the #reinit# functions. The reason is that
 * mostly, the matrices have different properties and you will want to
 * initialize the blocks making up the matrix separately. You can
 * access the different blocks using the #block(row,col)# function.
 *
 * Attention: this object is not automatically notified if the size of
 * one of its subobjects' size is changed. After you initialize the
 * sizes of the subobjects, you will therefore have to call the
 * #collect_sizes()# function of this class! Note that, of course, all
 * sub-matrices in a row have to have the same number of rows, and
 * that all sub-matrices in a column have to have the same number of
 * columns.
 *
 *
 * \section{On template instantiations}
 *
 * Member functions of this class are either implemented in this file
 * or in a file of the same name with suffix ``.templates.h''. For the
 * most common combinations of the template parameters, instantiations
 * of this class are provided in a file with suffix ``.cc'' in the
 * ``source'' directory. If you need an instantiation that is not
 * listed there, you have to include this file along with the
 * corresponding ``.templates.h'' file and instantiate the respective
 * class yourself.
 *
 * @author Wolfgang Bangerth, 2000
 */
template <int rows, int columns=rows>
class BlockSparsityPattern : public Subscriptor
{
  public:
				     /**
				      * Define a value which is used
				      * to indicate that a certain
				      * value in the #colnums# array
				      * is unused, i.e. does not
				      * represent a certain column
				      * number index.
				      *
				      * This value is only an alias to
				      * the respective value of the
				      * #SparsityPattern# class.
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
				      * the #reinit# function.
				      */
    BlockSparsityPattern ();

				     /**
				      * Copy operator. For this the same holds
				      * as for the copy constructor: it is
				      * declared, defined and fine to be called,
				      * but the latter only for empty objects.
				      */
    BlockSparsityPattern & operator = (const BlockSparsityPattern &);

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
    SparsityPattern &
    block (const unsigned int row,
	   const unsigned int column);
    
    
				     /**
				      * Access the block with the
				      * given coordinates. Version for
				      * constant objects.
				      */
    const SparsityPattern &
    block (const unsigned int row,
	   const unsigned int column) const;    
    
				     /**
				      * This function compresses the
				      * sparsity structures that this
				      * object represents. It simply
				      * calls #compress# for all
				      * sub-objects.
				      */
    void compress ();

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
				      * to which block #(i,j)# belongs
				      * and then relays to that block.
				      */
    void add (const unsigned int i, const unsigned int j);
	
				     /**
				      * Return number of rows of this
				      * matrix, which equals the
				      * dimension of the image
				      * space. It is the sum of rows
				      * of the rows of sub-matrices.
				      */
    unsigned int n_rows () const;

				     /**
				      * Return number of columns of
				      * this matrix, which equals the
				      * dimension of the range
				      * space. It is the sum of
				      * columns of the columns of
				      * sub-matrices.
				      */
    unsigned int n_cols () const;

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
				      * Return whether the structure
				      * is compressed or not,
				      * i.e. whether all sub-matrices
				      * are compressed.
				      */
    bool is_compressed () const;

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidSize,
		    int,
		    << "The number of blocks must be greater or equal to one, "
		    << "but is " << arg1);
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
    
    
  private:
				     /**
				      * Array of sparsity patterns.
				      */
    SparsityPattern sub_objects[rows][columns];

				     /**
				      * Object storing and managing
				      * the transformation of row
				      * indices to indices of the
				      * sub-objects.
				      */
    BlockIndices<rows>    row_indices;

				     /**
				      * Object storing and managing
				      * the transformation of column
				      * indices to indices of the
				      * sub-objects.
				      */
    BlockIndices<columns> column_indices;

				     /**
				      * Make the block sparse matrix a
				      * friend, so that it can use our
				      * #row_indices# and
				      * #column_indices# objects.
				      */
    template <typename number, int r, int c>
    friend class BlockSparseMatrix;
};



/*---------------------- Template functions -----------------------------------*/



template <int rows, int columns>
inline
SparsityPattern &
BlockSparsityPattern<rows,columns>::block (const unsigned int row,
					   const unsigned int column)
{
  return sub_objects[row][column];
};



template <int rows, int columns>
inline
const SparsityPattern &
BlockSparsityPattern<rows,columns>::block (const unsigned int row,
					   const unsigned int column) const
{
  return sub_objects[row][column];
};



template <int rows, int columns>
inline
void
BlockSparsityPattern<rows,columns>::add (const unsigned int i,
					 const unsigned int j)
{
				   // if you get an error here, are
				   // you sure you called
				   // #collect_sizes()# before?
  const pair<unsigned int,unsigned int>
    row_index = row_indices.global_to_local (i),
    col_index = column_indices.global_to_local (j);
  sub_objects[row_index.first][col_index.first].add (row_index.second,
						     col_index.second);
};



#endif
