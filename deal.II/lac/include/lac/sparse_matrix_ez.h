//----------------------------  sparse_matrix.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix.h  ---------------------------
#ifndef __deal2__sparse_matrix_ez_h
#define __deal2__sparse_matrix_ez_h


#include <base/config.h>
#include <base/exceptions.h>
#include <base/subscriptor.h>
#include <base/smartpointer.h>

#include <vector>

template<typename number> class Vector;
template<typename number> class FullMatrix;

/**
 * Sparse matrix without sparsity pattern.
 *
 * Instead of using a pre-assembled sparsity pattern, this matrix
 * builds the pattern on the fly. Filling the matrix may consume more
 * time as with @p{SparseMatrix}, since large memory movements may be
 * involved. To help optimizing things, an expected row-length may be
 * provided to the constructor, as well as a mininmum size increment
 * for rows.
 *
 * The name of this matrix is in reverence to a publication of the
 * Internal Revenue Service of the United States of America. I hope
 * some other aliens will appreciate it. By the way, the suffix makes
 * sense by pronouncing it the American way.
 *
 * @author Guido Kanschat, 2002
 */
template <typename number>
class SparseMatrixEZ : public Subscriptor
{
  private:
				     /**
				      * The class for storing the
				      * column number of an entry
				      * together with its value.
				      */
    struct Entry
    {
					 /**
					  * Standard constructor. Sets
					  * @p{column} to
					  * @p{invalid}.
					  */
	Entry();

					 /**
					  * Constructor. Fills column
					  * and value.
					  */
	Entry(unsigned int column,
	      const number& value);
	
					 /**
					  * The column number.
					  */
	unsigned int column;
					 /**
					  * The value there.
					  */
	number value;
					 /**
					  * Comparison operator for finding.
					  */
//	bool operator==(const Entry&) const;

					 /**
					  * Less than operator for sorting.
					  */
//	bool operator < (const Entry&) const;
					 /**
					  * Non-existent column number.
					  */
	static const unsigned int invalid = static_cast<unsigned int>(-1);
    };

				     /**
				      * Structure for storing
				      * information on a matrix
				      * row. One object for each row
				      * will be stored in the matrix.
				      */
    struct RowInfo
    {
					 /**
					  * Constructor.
					  */
	RowInfo (unsigned int start = Entry::invalid);
	
					 /**
					  * Index of first entry of
					  * the row in the data field.
					  */
	unsigned int start;
					 /**
					  * Number of entries in this
					  * row.
					  */
	unsigned short length;
					 /**
					  * Position of the diagonal
					  * element relative tor the
					  * start index.
					  */
	unsigned short diagonal;
					 /**
					  * Value for non-existing diagonal.
					  */
	static const unsigned short
	invalid_diagonal = static_cast<unsigned short>(-1);
    };
    
    
  public:
				     /**
				      * Type of matrix entries. In analogy to
				      * the STL container classes.
				      */
    typedef number value_type;
    
				     /**
				      * Constructor. Initializes an
				      * empty matrix of dimension zero
				      * times zero.
				      */
    SparseMatrixEZ ();

				     /**
				      * Copy constructor. This constructor is
				      * only allowed to be called if the matrix
				      * to be copied is empty. This is for the
				      * same reason as for the
				      * @p{SparsityPattern}, see there for the
				      * details.
				      *
				      * If you really want to copy a whole
				      * matrix, you can do so by using the
				      * @p{copy_from} function.
				      */
    SparseMatrixEZ (const SparseMatrixEZ &);

				     /**
				      * Constructor. Generates a
				      * matrix of the given size,
				      * ready to be filled. The
				      * optional parameters
				      * @p{default_row_length} and
				      * @p{default_increment} allow
				      * for preallocating
				      * memory. Providing these
				      * properly is essential for an
				      * efficient assembling of the
				      * matrix.
				      */
    explicit SparseMatrixEZ (unsigned int n_rows,
			     unsigned int n_columns = n_rows,
			     unsigned int default_row_length = Entry::invalid,
			     unsigned int default_increment = Entry::invalid);
    
				     /**
				      * Destructor. Free all memory, but do not
				      * release the memory of the sparsity
				      * structure.
				      */
    virtual ~SparseMatrixEZ ();

				     /** 
				      * Pseudo operator only copying
				      * empty objects.
				      */
    SparseMatrixEZ<number>& operator = (const SparseMatrixEZ<number> &);

				     /**
				      * Reinitialize the sparse matrix
				      * to the dimensions provided.
				      * The matrix will have no
				      * entries at this point. The
				      * optional parameters
				      * @p{default_row_length} and
				      * @p{default_increment} allow
				      * for preallocating
				      * memory. Providing these
				      * properly is essential for an
				      * efficient assembling of the
				      * matrix.
				      */
    virtual void reinit (unsigned int n_rows,
			 unsigned int n_columns = n_rows,
			 unsigned int default_row_length = Entry::invalid,
			 unsigned int default_increment = Entry::invalid);

				     /**
				      * Release all memory and return
				      * to a state just like after
				      * having called the default
				      * constructor. It also forgets
				      * the sparsity pattern it was
				      * previously tied to.
				      */
    virtual void clear ();
    
				     /**
				      * Return whether the object is
				      * empty. It is empty if
				      * both dimensions are zero.
				      */
    bool empty () const;

				     /**
				      * Return the dimension of the
				      * image space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int m () const;
    
				     /**
				      * Return the dimension of the
				      * range space.  To remember: the
				      * matrix is of dimension
				      * $m \times n$.
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
//    unsigned int n_nonzero_elements () const;

				     /**
				      * Return the number of actually
				      * nonzero elements of this
				      * matrix.
				      *
				      * Note, that this function does
				      * (in contrary to the
				      * @p{n_nonzero_elements}) NOT
				      * count all entries of the
				      * sparsity pattern but only the
				      * ones that are nonzero.
				      */
//    unsigned int n_actually_nonzero_elements () const;
    
				     /**
				      * Set the element @p{(i,j)} to
				      * @p{value}. Allocates the entry
				      * if it does not exist. Filters
				      * out zero automatically.
				      */
    void set (const unsigned int i, const unsigned int j,
	      const number value);
    
				     /**
				      * Add @p{value} to the element
				      * @p{(i,j)}. Allocates the entry
				      * if it does not exist. Filters
				      * out zero automatically.
				      */
    void add (const unsigned int i, const unsigned int j,
	      const number value);

				     /**
				      * Symmetrize the matrix by
				      * forming the mean value between
				      * the existing matrix and its
				      * transpose, $A = \frac 12(A+A^T)$.
				      *
				      * This operation assumes that
				      * the underlying sparsity
				      * pattern represents a symmetric
				      * object. If this is not the
				      * case, then the result of this
				      * operation will not be a
				      * symmetric matrix, since it
				      * only explicitly symmetrizes
				      * by looping over the lower left
				      * triangular part for efficiency
				      * reasons; if there are entries
				      * in the upper right triangle,
				      * then these elements are missed
				      * in the
				      * symmetrization. Symmetrization
				      * of the sparsity pattern can be
				      * obtain by the
				      * @ref{SparsityPattern}@p{::symmetrize}
				      * function.
				      */
//    void symmetrize ();
    
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
				      * through @p{operator =}, since
				      * this may lead to unwanted
				      * usage, e.g. in copy arguments
				      * to functions, which should
				      * really be arguments by
				      * reference.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      *
				      * The function returns a reference to
				      * @p{this}.
				      */
//    template <typename somenumber>
//    SparseMatrix<number> &
//    copy_from (const SparseMatrix<somenumber> &source);

				     /**
				      * This function is complete
				      * analogous to the
				      * @ref{SparsityPattern}@p{::copy_from}
				      * function in that it allows to
				      * initialize a whole matrix in
				      * one step. See there for more
				      * information on argument types
				      * and their meaning. You can
				      * also find a small example on
				      * how to use this function
				      * there.
				      *
				      * The only difference to the
				      * cited function is that the
				      * objects which the inner
				      * iterator points to need to be
				      * of type @p{std::pair<unsigned int, value},
				      * where @p{value}
				      * needs to be convertible to the
				      * element type of this class, as
				      * specified by the @p{number}
				      * template argument.
				      *
				      * Previous content of the matrix
				      * is overwritten. Note that the
				      * entries specified by the input
				      * parameters need not
				      * necessarily cover all elements
				      * of the matrix. Elements not
				      * covered remain untouched.
				      */
//    template <typename ForwardIterator>
//    void copy_from (const ForwardIterator begin,
//		    const ForwardIterator end);    

				     /**
				      * Copy the nonzero entries of a
				      * full matrix into this
				      * object. Previous content is
				      * deleted. Note that the
				      * underlying sparsity pattern
				      * must be appropriate to hold
				      * the nonzero entries of the
				      * full matrix.
				      */
//    template <typename somenumber>
//    void copy_from (const FullMatrix<somenumber> &matrix);
    
				     /**
				      * Add @p{matrix} scaled by
				      * @p{factor} to this matrix. The
				      * function throws an error if
				      * the sparsity patterns of the
				      * two involved matrices do not
				      * point to the same object,
				      * since in this case the
				      * operation is cheaper.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      */
//    template <typename somenumber>
//    void add_scaled (const number factor,
//		     const SparseMatrix<somenumber> &matrix);
    
				     /**
				      * Return the value of the entry
				      * (i,j).  This may be an
				      * expensive operation and you
				      * should always take care where
				      * to call this function.  In
				      * order to avoid abuse, this
				      * function throws an exception
				      * if the required element does
				      * not exist in the matrix.
				      *
				      * In case you want a function
				      * that returns zero instead (for
				      * entries that are not in the
				      * sparsity pattern of the
				      * matrix), use the @p{el}
				      * function.
				      */
    number operator () (const unsigned int i,
			const unsigned int j) const;

				     /**
				      * Return the main diagonal element in
				      * the @p{i}th row. This function throws an
				      * error if the matrix is not square.
				      *
				      * This function is considerably
				      * faster than the @p{operator()},
				      * since for square matrices, the
				      * diagonal entry is always the
				      * first to be stored in each row
				      * and access therefore does not
				      * involve searching for the
				      * right column number.
				      */
    number diag_element (const unsigned int i) const;

				     /**
				      * Same as above, but return a
				      * writeable reference. You're
				      * sure you know what you do?
				      */
    number & diag_element (const unsigned int i);
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M*src$ with $M$
				      * being this matrix.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>       &dst,
		const Vector<somenumber> &src) const;
    
				     /**
				      * Matrix-vector multiplication:
				      * let $dst = M^T*src$ with $M$
				      * being this matrix. This
				      * function does the same as
				      * @p{vmult} but takes the
				      * transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult (Vector<somenumber>       &dst,
		 const Vector<somenumber> &src) const;
  
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M*src$ on
				      * $dst$ with $M$ being this
				      * matrix.
				      */
    template <typename somenumber>
    void vmult_add (Vector<somenumber>       &dst,
		    const Vector<somenumber> &src) const;
    
				     /**
				      * Adding Matrix-vector
				      * multiplication. Add $M^T*src$
				      * to $dst$ with $M$ being this
				      * matrix. This function does the
				      * same as @p{vmult_add} but takes
				      * the transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult_add (Vector<somenumber>       &dst,
		     const Vector<somenumber> &src) const;
  
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
				      * be square for this operation.
				      */
    template <typename somenumber>
    somenumber matrix_norm_square (const Vector<somenumber> &v) const;

				     /**
				      * Compute the matrix scalar
				      * product $\left(u,Mv\right)$.
				      */
    template <typename somenumber>
    somenumber matrix_scalar_product (const Vector<somenumber> &u,
				      const Vector<somenumber> &v) const;
    
    				     /**
				      * Return the l1-norm of the matrix, that is
				      * $|M|_1=max_{all columns j}\sum_{all 
				      * rows i} |M_ij|$,
				      * (max. sum of columns).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the l1-norm for vectors, i.e.
				      * $|Mv|_1\leq |M|_1 |v|_1$.
				      * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
				      */
    number l1_norm () const;

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
				      * (cf. Haemmerlin-Hoffmann : Numerische Mathematik)
				      */
    number linfty_norm () const;

				     /**
				      * Apply the Jacobi
				      * preconditioner, which
				      * multiplies every element of
				      * the @p{src} vector by the
				      * inverse of the respective
				      * diagonal element and
				      * multiplies the result with the
				      * damping factor @p{omega}.
				      */
    template <typename somenumber>
    void precondition_Jacobi (Vector<somenumber>       &dst,
			      const Vector<somenumber> &src,
			      const number              omega = 1.) const;

				     /**
				      * Apply SSOR preconditioning to
				      * @p{src}.
				      */
    template <typename somenumber>
    void precondition_SSOR (Vector<somenumber>       &dst,
			    const Vector<somenumber> &src,
			    const number              om = 1.) const;

				     /**
				      * Apply SOR preconditioning matrix to @p{src}.
				      * The result of this method is
				      * $dst = (om D - L)^{-1} src$.
				      */
    template <typename somenumber>
    void precondition_SOR (Vector<somenumber>       &dst,
			   const Vector<somenumber> &src,
 			   const number              om = 1.) const;
    
				     /**
				      * Apply transpose SOR preconditioning matrix to @p{src}.
				      * The result of this method is
				      * $dst = (om D - U)^{-1} src$.
				      */
    template <typename somenumber>
    void precondition_TSOR (Vector<somenumber>       &dst,
			    const Vector<somenumber> &src,
			    const number              om = 1.) const;
    
				     /**
				      * Print the matrix to the given
				      * stream, using the format
				      * @p{(line,col) value}, i.e. one
				      * nonzero entry of the matrix
				      * per line.
				      */
    void print (std::ostream &out) const;

				     /**
				      * Print the matrix in the usual
				      * format, i.e. as a matrix and
				      * not as a list of nonzero
				      * elements. For better
				      * readability, elements not in
				      * the matrix are displayed as
				      * empty space, while matrix
				      * elements which are explicitly
				      * set to zero are displayed as
				      * such.
				      *
				      * The parameters allow for a
				      * flexible setting of the output
				      * format: @p{precision} and
				      * @p{scientific} are used to
				      * determine the number format,
				      * where @p{scientific} = @p{false}
				      * means fixed point notation.  A
				      * zero entry for @p{width} makes
				      * the function compute a width,
				      * but it may be changed to a
				      * positive value, if output is
				      * crude.
				      *
				      * Additionally, a character for
				      * an empty value may be
				      * specified.
				      *
				      * Finally, the whole matrix can
				      * be multiplied with a common
				      * denominator to produce more
				      * readable output, even
				      * integers.
				      *
				      * This function
				      * may produce @em{large} amounts of
				      * output if applied to a large matrix!
				      */
    void print_formatted (std::ostream       &out,
			  const unsigned int  precision   = 3,
			  const bool          scientific  = true,
			  const unsigned int  width       = 0,
			  const char         *zero_string = " ",
			  const double        denominator = 1.) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception for applying
				      * inverse-type operators to
				      * rectangular matrices.
				      */
    DeclException0(ExcNoSquareMatrix);
    
				     /**
				      * Exception for missing diagonal entry.
				      */
    DeclException0(ExcNoDiagonal);
    
  private:

				     /**
				      * Find an entry. Return a
				      * zero-pointer if the entry does
				      * not exist.
				      */
    const Entry* locate (unsigned int row,
			 unsigned int col) const;

				     /**
				      * Find an entry or generate it.
				      */
    Entry* allocate (unsigned int row,
		     unsigned int col);
    
				     /**
				      * Version of @p{vmult} which only
				      * performs its actions on the
				      * region defined by
				      * @p{[begin_row,end_row)}. This
				      * function is called by @p{vmult}
				      * in the case of enabled
				      * multithreading.
				      */
    template <typename somenumber>
    void threaded_vmult (Vector<somenumber>       &dst,
			 const Vector<somenumber> &src,
			 const unsigned int        begin_row,
			 const unsigned int        end_row) const;

				     /**
				      * Version of
				      * @p{matrix_norm_square} which
				      * only performs its actions on
				      * the region defined by
				      * @p{[begin_row,end_row)}. This
				      * function is called by
				      * @p{matrix_norm_square} in the
				      * case of enabled
				      * multithreading.
				      */
    template <typename somenumber>
    void threaded_matrix_norm_square (const Vector<somenumber> &v,
				      const unsigned int        begin_row,
				      const unsigned int        end_row,
				      somenumber               *partial_sum) const;

    				     /**
				      * Version of
				      * @p{matrix_scalar_product} which
				      * only performs its actions on
				      * the region defined by
				      * @p{[begin_row,end_row)}. This
				      * function is called by
				      * @p{matrix_scalar_product} in the
				      * case of enabled
				      * multithreading.
				      */
    template <typename somenumber>
    void threaded_matrix_scalar_product (const Vector<somenumber> &u,
					 const Vector<somenumber> &v,
					 const unsigned int        begin_row,
					 const unsigned int        end_row,
					 somenumber               *partial_sum) const;

				     /**
				      * Version of @p{residual} which
				      * only performs its actions on
				      * the region defined by
				      * @p{[begin_row,end_row)} (these
				      * numbers are the components of
				      * @p{interval}). This function is
				      * called by @p{residual} in the
				      * case of enabled
				      * multithreading.
				      */
    template <typename somenumber>
    void threaded_residual (Vector<somenumber>       &dst,
			    const Vector<somenumber> &u,
			    const Vector<somenumber> &b,
			    const std::pair<unsigned int,unsigned int> interval,
			    somenumber               *partial_norm) const;

				     /**
				      * Number of columns. This is
				      * used to check vector
				      * dimensions only.
				      */
    unsigned int n_columns;

				     /**
				      * Info structure for each row.
				      */
    std::vector<RowInfo> row_info;
    
				     /**
				      * Data storage.
				      */
    std::vector<Entry> data;

				     /**
				      * Increment when a row grows.
				      */
    unsigned int increment;
    
				     // make all other sparse matrices
				     // friends
    template <typename somenumber> friend class SparseMatrix;
};


/*---------------------- Inline functions -----------------------------------*/

template <typename number>
inline
SparseMatrixEZ<number>::Entry::Entry(unsigned int column,
				     const number& value)
		:
		column(column),
  value(value)
{}



template <typename number>
inline
SparseMatrixEZ<number>::Entry::Entry()
		:
		column(invalid),
  value(0)
{}


// template <typename number>
// inline
// bool
// SparseMatrixEZ<number>::Entry::operator==(const Entry& e) const
// {
//   return column == e.column;
// }


// template <typename number>
// inline
// bool
// SparseMatrixEZ<number>::Entry::operator<(const Entry& e) const
// {
//   return column < e.column;
// }

template <typename number>
inline
SparseMatrixEZ<number>::RowInfo::RowInfo(unsigned int start)
		:
		start(start), length(0), diagonal(invalid_diagonal)
{}


//----------------------------------------------------------------------//
template <typename number>
inline
unsigned int SparseMatrixEZ<number>::m () const
{
  return row_info.size();
};


template <typename number>
inline
unsigned int SparseMatrixEZ<number>::n () const
{
  return n_columns;
};


template <typename number>
inline
const SparseMatrixEZ<number>::Entry* SparseMatrixEZ<number>::locate (
  const unsigned int row,
  const unsigned int col) const
{
  Assert (row<m(), ExcIndexRange(row,0,m()));
  Assert (col<n(), ExcIndexRange(col,0,n()));

  const RowInfo& r = row_info[row];
  const unsigned int end = r.start + r.length;
  for (unsigned int i=r.start;i<end;++i)
    {
      const Entry * const entry = &data[i];
      if (entry->column == col)
	return entry;
      if (entry->column == Entry::invalid)
	return 0;
    }
  return 0;
}



template <typename number>
inline
SparseMatrixEZ<number>::Entry* SparseMatrixEZ<number>::allocate (
  const unsigned int row,
  const unsigned int col)
{
  Assert (row<m(), ExcIndexRange(row,0,m()));
  Assert (col<n(), ExcIndexRange(col,0,n()));

  RowInfo& r = row_info[row];
  const unsigned int end = r.start + r.length;

  unsigned int i = r.start;
  if (r.diagonal != RowInfo::invalid_diagonal && col >= row)
    i += r.diagonal;
				   // Find position of entry
  while (i<end && data[i].column < col) ++i;
  
  Entry* entry = &data[i];
				   // entry found
  if (entry->column == col)
    return entry;

				   // Now, we must insert the new
				   // entry and move all successive
				   // entries back.
  
				   // If no more space is available
				   // for this row, insert new
				   // elements into the vector.
  if (row != row_info.size()-1)
    {
      if (end >= row_info[row+1].start)
	{
					   // Insert new entries
	  data.insert(data.begin()+end, increment, Entry());
	  entry = &data[i];
					   // Update starts of
					   // following rows
	  for (unsigned int rn=row+1;rn<row_info.size();++rn)
	    row_info[rn].start += increment;
	}
    } else {
      if (end >= data.size())
	{
					   // Here, appending a block
					   // does not increase
					   // performance.
	  data.push_back(Entry());
	  entry = &data[i];
	}
    }
				   // Save original entry
  Entry temp = *entry;
				   // Insert new entry here to
				   // make sure all entries
				   // are ordered by column
				   // index
  entry->column = col;
  entry->value = 0;
				   // Update row_info
  ++r.length;
  if (col == row)
    r.diagonal = i - r.start;
  else if (col<row && r.diagonal!= RowInfo::invalid_diagonal)
    ++r.diagonal;

  if (i == end)
      return entry;
  
				   // Move all entries in this
				   // row up by one
  for (unsigned int j = i+1; j < end; ++j)
    {
				       // There should be no invalid
				       // entry below end
      Assert (data[j].column != Entry::invalid, ExcInternalError());
      Entry temp2 = data[j];
      data[j] = temp;
      temp = temp2;
    }
  Assert (data[end].column == Entry::invalid, ExcInternalError());
  data[end] = temp;

  return entry;
}



template <typename number>
inline
void SparseMatrixEZ<number>::set (const unsigned int i,
				  const unsigned int j,
				  const number value)
{
  Assert (i<m(), ExcIndexRange(i,0,m()));
  Assert (j<n(), ExcIndexRange(j,0,n()));
  Entry* entry = allocate(i,j);
  entry->value = value;
};



template <typename number>
inline
void SparseMatrixEZ<number>::add (const unsigned int i,
				  const unsigned int j,
				  const number value)
{
  Assert (i<m(), ExcIndexRange(i,0,m()));
  Assert (j<n(), ExcIndexRange(j,0,n()));
  Entry* entry = allocate(i,j);
  entry->value += value;
};


#endif
/*----------------------------   sparse_matrix.h     ---------------------------*/
