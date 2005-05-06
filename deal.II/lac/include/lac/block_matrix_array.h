//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__block_matrix_array_h
#define __deal2__block_matrix_array_h

#include <base/config.h>
#include <base/subscriptor.h>
#include <base/table.h>

#include <lac/pointer_matrix.h>
#include <lac/vector_memory.h>

#include <vector>
#include <map>
#include <string>
#include <memory>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


template <typename> class BlockVector;
template <typename> class Vector;

/*! @addtogroup Matrix2
 *@{
 */


/**
 * Block matrix composed of different single matrices; these matrices
 * may even be of different types.
 *
 * Given a set of arbitrary matrices <i>A<sub>i</sub></i>, this class
 * implements a block matrix with block entries of the form
 * <i>M<sub>jk</sub> = s<sub>jk</sub>A<sub>i</sub></i>.  Each
 * <i>A<sub>i</sub></i> may be used several times with different
 * prefix. The matrices are not copied into the BlockMatrixArray
 * object, but rather a PointerMatrix referencing each of them will be
 * stored along with factors and transposition flags.
 *
 * Non-zero entries are registered by the function enter(), zero
 * entries are not stored at all. Using enter() with the same location
 * <tt>(i,j)</tt> several times will add the corresponding matrices in
 * matrix-vector multiplications. These matrices will not be actually
 * added, but the multiplications with them will be summed up.
 *
 * @note This mechanism makes it impossible to access single entries
 * of BlockMatrixArray. In particular, (block) relaxation
 * preconditioners based on PreconditionRelaxation or
 * PreconditionBlock <b>cannot</b> be used with this class. If you
 * need a preconditioner for a BlockMatrixArray object, use
 * BlockTrianglePrecondition.
 *
 * <h3>Requirements on MATRIX</h3>
 *
 * The template argument <tt>MATRIX</tt> is a class providing the
 * matrix-vector multiplication functions vmult(), Tvmult(),
 * vmult_add() and Tvmult_add() used in this class, but with arguments
 * of type Vector&lt;number&gt; instead of
 * BlockVector&lt;number&gt;. Every matrix which can be used by
 * PointerMatrix is allowed, in particular SparseMatrix is a possible
 * entry type.
 *
 * <h3>Example program</h3>
 * We document the relevant parts of <tt>examples/doxygen/block_matrix_array.cc</tt>.
 *
 * @dontinclude block_matrix_array.cc
 *
 * Obviously, we have to include the header file containing the definition
 * of BlockMatrixArray:
 * @skipline block_matrix_array.h
 *
 * First, we set up some matrices to be entered into the blocks.
 * @skip main
 * @until C.fill
 *
 * The BlockMatrixArray needs a VectorMemory&lt;Vector&lt;number&gt;
 * &gt; object to allocate auxiliary vectors. <tt>number</tt> must
 * equal the second template argument of BlockMatrixArray and also the
 * number type of the BlockVector used later. We use the
 * GrowingVectorMemory type, since it remembers the vector and avoids
 * reallocating.
 *
 * @ line Growing
 *
 * Now, we are ready to build a <i>2x2</i> BlockMatrixArray.
 * @line Block
 * First, we enter the matrix <tt>A</tt> multiplied by 2 in the upper left block
 * @line enter
 * Now -1 times <tt>B1</tt> in the upper right block.
 * @line enter
 * We add the transpose of <tt>B2</tt> to the upper right block and
 * continue in a similar fashion. In the end, the block matrix
 * structure is printed into an LaTeX table.
 * @until latex
 *
 * Now, we set up vectors to be multiplied with this matrix and do a
 * multiplication.
 * @until vmult
 *
 * Finally, we solve a linear system with BlockMatrixArray, using no
 * preconditioning and the conjugate gradient method.
 * @until Error
 *
 * The remaining code of this sample program concerns preconditioning
 * and is described in the documentation of
 * BlockTrianglePrecondition.
 *
 * @author Guido Kanschat, 2000 - 2005
 */
template <typename number = double>
class BlockMatrixArray : public Subscriptor
{
  public:
				     /**
				      * Default constructor creating a
				      * useless object. initialize()
				      * must be called before using
				      * it.
				      */
    BlockMatrixArray ();
    
				     /**
				      * Constructor fixing the
				      * dimensions.
				      */
    BlockMatrixArray (const unsigned int n_block_rows,
		      const unsigned int n_block_cols,
		      VectorMemory<Vector<number> >& mem);

				     /**
				      * Initialize object
				      * completely. This is the
				      * function to call for an object
				      * created by the default
				      * constructor.
				      */
    void initialize (const unsigned int n_block_rows,
		     const unsigned int n_block_cols,
		     VectorMemory<Vector<number> >& mem);
    
				     /**
				      * Adjust the matrix to a new
				      * size and delete all blocks.
				      */
    void reinit (const unsigned int n_block_rows,
		 const unsigned int n_block_cols);
    
				     /**
				      * Add a block matrix entry. The
				      * <tt>matrix</tt> is entered
				      * into a list of blocks for
				      * multiplication, together with
				      * its coordinates <tt>row</tt>
				      * and <tt>col</tt> as well as
				      * optional multiplication factor
				      * <tt>prefix</tt> and transpose
				      * flag <tt>transpose</tt>.
				      *
				      * @note No check for consistency
				      * of block sizes is
				      * made. Therefore, entering a
				      * block of wrong dimension here
				      * will only lead to a
				      * ExcDimensionMismatch in one of
				      * the multiplication functions.
				      */
    template <class MATRIX>
    void enter (const MATRIX      &matrix,
		const unsigned int row,
		const unsigned int col,
		const double       prefix = 1.,
		const bool         transpose = false);
    
				     /**
				      * Delete all entries, i.e. reset
				      * the matrix to an empty state.
				      */
    void clear();
    
				     /**
				      * Number of block-entries per
				      * column.
				      */
    unsigned int n_block_rows () const;

				     /**
				      * Number of block-entries per
				      * row.
				      */
    unsigned int n_block_cols () const;

				     /**
				      * Matrix-vector multiplication.
				      */
    void vmult (BlockVector<number>& dst,
		const BlockVector<number>& src) const;
  
				     /**
				      * Matrix-vector multiplication
				      * adding to <tt>dst</tt>.
				      */
    void vmult_add (BlockVector<number>& dst,
		    const BlockVector<number>& src) const;
  
				     /**
				      * Transposed matrix-vector
				      * multiplication.
				      */
    void Tvmult (BlockVector<number>& dst,
		 const BlockVector<number>& src) const;
  
				     /**
				      * Transposed matrix-vector
				      * multiplication adding to
				      * <tt>dst</tt>.
				      */
    void Tvmult_add (BlockVector<number>& dst,
		     const BlockVector<number>& src) const;

				     /**
				      * Matrix scalar product between
				      * two vectors (at least for a
				      * symmetric matrix).
				      */
    number matrix_scalar_product (const BlockVector<number>& u,
				  const BlockVector<number>& v) const;
    
				     /**
				      * Square of the matrix norm of a
				      * vector (at least for a
				      * symmetric matrix).
				      */
    number matrix_norm_square (const BlockVector<number>& u) const;
    
				     /**
				      * Print the block structure as a
				      * LaTeX-array. This output will
				      * not be very intuitive, since
				      * the matrix object lacks
				      * important information. What
				      * you see is an entry for each
				      * block showing all the matrices
				      * with their multiplicaton
				      * factors and possibly transpose
				      * mark. The matrices itself are
				      * numbered successively upon
				      * being entred. If the same
				      * matrix is entered several
				      * times, it will be listed with
				      * the same number everytime.
				      */
    template <class STREAM>
    void print_latex (STREAM& out) const;
    
  protected:
				     /**
				      * Internal data structure.
				      *
				      * For each entry of a
				      * BlockMatrixArray, its
				      * position, matrix, prefix and
				      * optional transposition must be
				      * stored. This structure
				      * encapsulates all of them.
				      *
				      * @author Guido Kanschat, 2000, 2001
				      */
    class Entry
    {
      public:
					 /**
					  * Constructor initializing all
					  * data fields.
					  */
	template<class MATRIX> 
	Entry (const MATRIX& matrix,
	       unsigned row, unsigned int col,
	       double prefix, bool transpose);
					 /**
					  * Copy constructor
					  * invalidating the old
					  * object. Since it is only
					  * used for entering
					  * temporary objects into a
					  * vector, this is ok.
					  *
					  * For a deep copy, we would
					  * need a reproduction
					  * operator in
					  * PointerMatixBase.
					  */
	Entry(const Entry&);
	
					 /**
					  * Destructor, where we
					  * delete the PointerMatrix
					  * created by the
					  * constructor.
					  */
	~Entry();
	
					 /**
					  * Row number in the block
					  * matrix.
					  */
	unsigned int row;

					 /**
					  * Column number in the block
					  * matrix.
					  */
	unsigned int col;

					 /**
					  * Factor in front of the matrix
					  * block.
					  */
	double prefix;

					 /**
					  * Indicates that matrix block
					  * must be transposed for
					  * multiplication.
					  */
	bool transpose;

					 /**
					  * The matrix block itself.
					  */
	PointerMatrixBase<Vector<number> >* matrix;
    };
  
				     /**
				      * Array of block entries in the
				      * matrix.
				      */
    std::vector<Entry> entries;

  private:
				     /**
				      * Number of blocks per column.
				      */
    unsigned int block_rows;
				     /**
				      * number of blocks per row.
				      */
    unsigned int block_cols;
  protected:
				     /**
				      * The memory used for auxiliary
				      * vectors.
				      */
    mutable SmartPointer<VectorMemory<Vector<number> > > mem;
};

/*@}*/

/*! @addtogroup Preconditioners
 *@{
 */


/**
 * Inversion of a block-triangular matrix.
 *
 * In this block matrix, the inverses of the diagonal blocks are
 * stored together with the off-diagonal blocks of a block
 * matrix. Then, forward or backward insertion is performed
 * block-wise. The diagonal blocks are NOT inverted for this purpose!
 *
 * Like for all preconditioners, the preconditioning operation is
 * performed by the vmult() member function.
 *
 * @note While block indices may be duplicated (see BlockMatrixArray)
 * to add blocks, this is not allowed for diagonal blocks, since
 * summing up the inverse of two blocks does not yield the inverse of
 * the sum.
 *
 * The implementation may be a little clumsy, but it should be
 * sufficient as long as the block sizes are much larger than the
 * number of blocks.
 *
 * <h3>Example</h3>
 * Here, we document the second part of
 * <tt>examples/doxygen/block_matrix_array.cc</tt>. For the beginning
 * of this file, see BlockMatrixArray.
 *
 * In order to set up the preconditioner, we have to compute the
 * inverses of the diagonal blocks ourselves. Since we used FullMatrix
 * objects, this is fairly easy.
 * @dontinclude block_matrix_array.cc
 * @skip Error
 * @until Cinv.invert
 *
 * After creating a <i>2x2</i> BlockTrianglePrecondition object, we
 * only fill its diagonals. The scaling factor <i>1/2</i> used for
 * <tt>A</tt> is the reciprocal of the scaling factor used for the
 * <tt>matrix</tt> itself. Remember, this preconditioner actually
 * <b>multiplies</b> with the diagonal blocks.
 * @until Cinv
 *
 * Now, we have a block Jacobi preconditioner, which is still
 * symmetric, since the blocks are symmetric. Therefore, we can still
 * use the preconditioned conjugate gradient method.
 * @until Error
 *
 * Now, we enter the subdiagonal block. This is the same as in
 * <tt>matrix</tt>.
 * @until B2
 *
 * Since the preconditioner is not symmetric anymore, we use the GMRES
 * method for solving.
 * @until Error
 *
 * @author Guido Kanschat, 2001, 2005
 */
template <typename number = double>
class BlockTrianglePrecondition
  : private BlockMatrixArray<number>
{
  public:
				     /**
				      * Default constructor creating a
				      * useless object. initialize()
				      * must be called before using
				      * it.
				      */
    BlockTrianglePrecondition ();
    
				     /**
				      * Constructor. This matrix must be
				      * block-quadratic. The additional
				      * parameter allows for backward
				      * insertion instead of forward.
				      */
    BlockTrianglePrecondition (unsigned int n_block_rows,
			       VectorMemory<Vector<number> >& mem,
			       bool backward = false);
    
				     /**
				      * Initialize object
				      * completely. This is the
				      * function to call for an object
				      * created by the default
				      * constructor.
				      */
    void initialize (const unsigned int n_block_rows,
		     VectorMemory<Vector<number> >& mem,
		     bool backward = false);
    
				     /**
				      * Resize preconditioner to a new
				      * size and clear all blocks.
				      */
    void reinit(const unsigned int n_block_rows);
    
    
				     /**
				      * Enter a block. This calls
				      * BlockMatrixArray::enter(). Remember
				      * that the diagonal blocks
				      * should actually be inverse
				      * matrices or preconditioners.
				      */
    template <class MATRIX>
    void enter (const MATRIX      &matrix,
		const unsigned int row,
		const unsigned int col,
		const double       prefix = 1.,
		const bool         transpose = false);

                                     /**
				      * Preconditioning.
				      */
    void vmult (BlockVector<number>& dst,
		const BlockVector<number>& src) const;
  
				     /**
				      * Preconditioning
				      * adding to <tt>dst</tt>.
				      */
    void vmult_add (BlockVector<number>& dst,
		    const BlockVector<number>& src) const;
  
				     /**
				      * Transposed preconditioning
				      */
    void Tvmult (BlockVector<number>& dst,
		 const BlockVector<number>& src) const;
  
				     /**
				      * Transposed preconditioning
				      * adding to <tt>dst</tt>.
				      */
    void Tvmult_add (BlockVector<number>& dst,
		     const BlockVector<number>& src) const;

				     /**
				      * Make function of base class available.
				      */
    BlockMatrixArray<number>::print_latex;

				     /**
				      * Make function of base class available.
				      */
    BlockMatrixArray<number>::n_block_rows;

				     /**
				      * Make function of base class available.
				      */
    BlockMatrixArray<number>::n_block_cols;
    BlockMatrixArray<number>::clear;
    BlockMatrixArray<number>::Subscriptor::subscribe;
    BlockMatrixArray<number>::Subscriptor::unsubscribe;

      				     /** @addtogroup Exceptions
				      * @{ */

				     /**
				      * Multiple diagonal element.
				      */
    DeclException1(ExcMultipleDiagonal,
		   unsigned int,
		   << "Inverse diagonal entries may not be added in block "
		   << arg1);
				     //@}    
  private:
				     /**
				      * Add all off-diagonal
				      * contributions and return the
				      * entry of the diagonal element
				      * for one row.
				      */
    void do_row (BlockVector<number>& dst,
		 unsigned int row_num) const;
    
				     /**
				      * Flag for backward insertion.
				      */
    bool backward;
};

/*@}*/

///@if NoDoc
//---------------------------------------------------------------------------

template <typename number>
template <class MATRIX>
inline
BlockMatrixArray<number>::Entry::Entry (const MATRIX& matrix,
					unsigned row, unsigned int col,
					double prefix, bool transpose)
		:
		row (row),
		col (col),
		prefix (prefix),
		transpose (transpose),
		matrix (new PointerMatrix<MATRIX, Vector<number> >(&matrix, typeid(*this).name()))
{}



template <typename number>
template <class MATRIX>
inline
void
BlockMatrixArray<number>::enter (const MATRIX& matrix,
				 unsigned row, unsigned int col,
				 double prefix, bool transpose)
{
  Assert (mem != 0, ExcNotInitialized());
  Assert(row<n_block_rows(), ExcIndexRange(row, 0, n_block_rows()));
  Assert(col<n_block_cols(), ExcIndexRange(col, 0, n_block_cols()));
  entries.push_back(Entry(matrix, row, col, prefix, transpose));
}


template<typename number>
struct BlockMatrixArrayPointerMatrixLess
{
    bool operator () (const PointerMatrixBase<Vector<number> >*a,
		      const PointerMatrixBase<Vector<number> >* b) const
      {
	return *a < *b;
      }
};
  

template <typename number>
template <class STREAM>
inline
void
BlockMatrixArray<number>::print_latex (STREAM& out) const
{
  Assert (mem != 0, ExcNotInitialized());
  out << "\\begin{array}{"
      << std::string(n_block_cols(), 'c')
      << "}" << std::endl;

  Table<2,std::string> array(n_block_rows(), n_block_cols());
  
  typedef std::map<const PointerMatrixBase<Vector<number> >*,
    std::string, BlockMatrixArrayPointerMatrixLess<number> > NameMap;
  NameMap matrix_names;
  
  typename std::vector<Entry>::const_iterator m = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();

  unsigned int matrix_number = 0;
  for (; m != end ; ++m)
    {
      if (matrix_names.find(m->matrix) == matrix_names.end())
	{
	  std::pair<typename NameMap::iterator, bool> x =
	    matrix_names.insert(
	      std::pair<const PointerMatrixBase<Vector<number> >*, std::string> (m->matrix,
                                                     std::string("M")));
#ifdef HAVE_STD_STRINGSTREAM
          std::ostringstream stream;
#else
          std::ostrstream stream;
#endif
	  
	  stream << matrix_number++;
	  
#ifndef HAVE_STD_STRINGSTREAM
          stream << std::ends;
#endif
	  x.first->second += stream.str();
	}
      
#ifdef HAVE_STD_STRINGSTREAM
      std::ostringstream stream;
#else
      std::ostrstream stream;
#endif
      if (array(m->row, m->col) != "" && m->prefix >= 0)
	stream << "+";
      if (m->prefix != 1.)
	stream << m->prefix << 'x';
      stream << matrix_names.find(m->matrix)->second;
//      stream << '(' << m->matrix << ')';
      if (m->transpose)
	stream << "^T";

#ifndef HAVE_STD_STRINGSTREAM
      stream << std::ends;
#endif
      array(m->row, m->col) += stream.str();
    }
  for (unsigned int i=0;i<n_block_rows();++i)
    for (unsigned int j=0;j<n_block_cols();++j)
      {
	out << '\t' << array(i,j);
	if (j==n_block_cols()-1)
	  out << "\\\\" << std::endl;
	else
	  out << " &";
      }
  out << "\\end{array}" << std::endl;
}

template <typename number>
template <class MATRIX>
inline
void
BlockTrianglePrecondition<number>::enter (const MATRIX& matrix,
					  unsigned row, unsigned int col,
					  double prefix, bool transpose)
{
  BlockMatrixArray<number>::enter(matrix, row, col, prefix, transpose);
}



///@endif

#endif
