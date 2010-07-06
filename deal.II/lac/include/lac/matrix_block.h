//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2007, 2008, 2009, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__matrix_block_h
#define __deal2__matrix_block_h

#include <base/config.h>
#include <base/named_data.h>
#include <base/smartpointer.h>
#include <base/std_cxx1x/shared_ptr.h>
#include <lac/block_indices.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/sparse_matrix.h>
#include <lac/full_matrix.h>


DEAL_II_NAMESPACE_OPEN

/**
 * A wrapper around a matrix object, storing the coordinates in a
 * block matrix as well.
 *
 * This class is an alternative to BlockMatrixBase, if you only want
 * to generate a single block of the system, not the whole
 * system. Using the add() functions of this class, it is possible to
 * use the standard assembling functions used for block matrices, but
 * only enter in one of the blocks and still avoiding the index
 * computations innvolved.

 * The reason for this class is, that we may need a different number
 * of matrices for different blocks in a block system. For example, a
 * preconditioner for the Oseen system can be built as a block system,
 * where the pressure block is of the form
 * <b>M</b><sup>-1</sup><b>FA</b><sup>-1</sup> with <b>M</b> the
 * pressure mass matrix, <b>A</b> the pressure Laplacian and <b>F</b>
 * the advection diffusion operator applied to the pressure
 * space. Since only a single matrix is needed for the other blocks,
 * using BlockSparseMatrix or similar would be a waste of memory.
 *
 * While the add() functions make a MatrixBlock appear like a block
 * matrix for assembling, the functions vmult(), Tvmult(),
 * vmult_add(), and Tvmult_add() make it behave like a MATRIX, when it
 * comes to applying it to a vector. This behavior allows us to store
 * MatrixBlock objects in vectors, for instance in MGLevelObject
 * without extracting the #matrix first.
 *
 * MatrixBlock comes handy when using BlockMatrixArray. Once the
 * MatrixBlock has been properly initalized and filled, it can be used
 * in the simplest case as:
 * @code
 * MatrixBlockVector<SparseMatrix<double> > > blocks;
 *
 * ...
 *
 * BlockMatrixArray matrix (n_blocks, n_blocks);
 *
 * for (unsigned int i=0;i<blocks.size;++i)
 *   matrix.enter(blocks.block(i).row, blocks.block(i).column, blocks.matrix(i));
 * @endcode
 *
 * Here, we have not gained very much, except that we do not need to
 * set up empty blocks in the block system.
 *
 * @todo Example for the product preconditioner of the pressure Schur
 * complement.
 *
 * @ingroup Matrix2
 * @ingroup vector_valued
 * @see @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat, 2006
 */
template <class MATRIX>
class MatrixBlock
{
  public:
				     /**
				      * Constructor rendering an
				      * uninitialized object.
				      */
    MatrixBlock();

				     /**
				      * Copy constructor.
				      */
    MatrixBlock(const MatrixBlock<MATRIX>& M);

				     /**
				      * Constructor setting block
				      * coordinates, but not
				      * initializing the matrix.
				      */

    MatrixBlock(unsigned int i, unsigned int j,
		const BlockIndices* block_indices = 0);

				     /**
				      * Add <tt>value</tt> to the
				      * element (<i>i,j</i>). Throws
				      * an error if the entry does not
				      * exist or if it is in a
				      * different block.
				      */
    void add (const unsigned int i,
              const unsigned int j,
	      const typename MATRIX::value_type value);

                                       /**
                                        * Add all elements in a
                                        * FullMatrix into sparse
                                        * matrix locations given by
                                        * <tt>indices</tt>. This function
                                        * assumes a quadratic sparse
                                        * matrix and a quadratic
                                        * full_matrix.  The global
                                        * locations are translated
                                        * into locations in this block
                                        * and ExcBlockIndexMismatch is
                                        * thrown, if the global index
                                        * does not point into the
                                        * block refered to by #row and
                                        * #column.
					*
					* @todo
					* <tt>elide_zero_values</tt> is
					* currently ignored.
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
    void add (const std::vector<unsigned int>& indices,
	      const FullMatrix<number>&        full_matrix,
	      const bool                       elide_zero_values = true);

                                       /**
					* Add all elements in a
                                        * FullMatrix into global
                                        * locations given by
                                        * <tt>row_indices</tt> and
                                        * <tt>col_indices</tt>,
                                        * respectively. The global
                                        * locations are translated
                                        * into locations in this block
                                        * and ExcBlockIndexMismatch is
                                        * thrown, if the global index
                                        * does not point into the
                                        * block refered to by #row and
                                        * #column.
					*
					* @todo
					* <tt>elide_zero_values</tt> is
					* currently ignored.
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
    void add (const std::vector<unsigned int>& row_indices,
	      const std::vector<unsigned int>& col_indices,
	      const FullMatrix<number>&        full_matrix,
	      const bool                       elide_zero_values = true);

                                       /**
                                        * Set several elements in the
                                        * specified row of the matrix
                                        * with column indices as given
                                        * by <tt>col_indices</tt> to
                                        * the respective value. This
                                        * is the function doing thye
                                        * actual work for the ones
                                        * adding full matrices. The
                                        * global locations
                                        * <tt>row_index</tt> and
                                        * <tt>col_indices</tt> are
                                        * translated into locations in
                                        * this block and
                                        * ExcBlockIndexMismatch is
                                        * thrown, if the global index
                                        * does not point into the
                                        * block refered to by #row and
                                        * #column.
					*
					* @todo
					* <tt>elide_zero_values</tt> is
					* currently ignored.
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
    void add (const unsigned int               row_index,
	      const std::vector<unsigned int>& col_indices,
	      const std::vector<number>&       values,
	      const bool                       elide_zero_values = true);

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
    void add (const unsigned int  row,
	      const unsigned int  n_cols,
	      const unsigned int *col_indices,
	      const number       *values,
	      const bool          elide_zero_values = true,
	      const bool          col_indices_are_sorted = false);

				     /**
				      * Matrix-vector-multiplication,
				      * forwarding to the same
				      * function in MATRIX. No index
				      * computations are done, thus,
				      * the vectors need to have sizes
				      * matching #matrix.
				      */
    template<class VECTOR>
    void vmult (VECTOR& w, const VECTOR& v) const;

				     /**
				      * Matrix-vector-multiplication,
				      * forwarding to the same
				      * function in MATRIX. No index
				      * computations are done, thus,
				      * the vectors need to have sizes
				      * matching #matrix.
				      */
    template<class VECTOR>
    void vmult_add (VECTOR& w, const VECTOR& v) const;

				     /**
				      * Matrix-vector-multiplication,
				      * forwarding to the same
				      * function in MATRIX. No index
				      * computations are done, thus,
				      * the vectors need to have sizes
				      * matching #matrix.
				      */
    template<class VECTOR>
    void Tvmult (VECTOR& w, const VECTOR& v) const;

				     /**
				      * Matrix-vector-multiplication,
				      * forwarding to the same
				      * function in MATRIX. No index
				      * computations are done, thus,
				      * the vectors need to have sizes
				      * matching #matrix.
				      */
    template<class VECTOR>
    void Tvmult_add (VECTOR& w, const VECTOR& v) const;

				     /**
				      * The block number computed from
				      * an index by using
				      * #block_indices does not match
				      * the block coordinates stored
				      * in this object.
				      */
    DeclException2(ExcBlockIndexMismatch, unsigned int, unsigned int,
		   << "Block index " << arg1 << " does not match " << arg2);

				     /**
				      * Row coordinate.  This is the
				      * position of the data member
				      * matrix on the global matrix.
				      */
    unsigned int row;
				     /**
				      * Column coordinate.  This is
				      * the position of the data
				      * member matrix on the global
				      * matrix.
				      */
    unsigned int column;

				     /**
				      * The matrix itself
				      */
    MATRIX matrix;

				     /**
				      * The BlockIndices of the whole
				      * system. Using row() and
				      * column(), this allows us to
				      * find the index of the first
				      * row and column degrees of
				      * freedom for this block.
				      */
    SmartPointer<const BlockIndices, MatrixBlock<MATRIX> > block_indices;
};


/**
 * A vector of MatrixBlock, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each
 * matrix block can be identified by name.
 *
 * @relates MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <class MATRIX>
class MatrixBlockVector : public NamedData<std_cxx1x::shared_ptr<MatrixBlock<MATRIX> > >
{
  public:
				     /**
				      * Add a new matrix block at the
				      * position <tt>(row,column)</tt>
				      * in the block system.
				      */
    void add(unsigned int row, unsigned int column,
	     const std::string& name,
	     const BlockIndices* block_indices);

				     /**
				      * For matrices using a
				      * SparsityPattern, this function
				      * reinitializes each matrix in
				      * the vector with the correct
				      * pattern from the block system.
				      */
    void reinit(const BlockSparsityPattern& sparsity);
    
				     /**
				      * Access a constant reference to
				      * the block at position <i>i</i>.
				      */
    const MatrixBlock<MATRIX>& block(unsigned int i) const;
    
				     /**
				      * Access a reference to
				      * the block at position <i>i</i>.
				      */
    MatrixBlock<MATRIX>& block(unsigned int i);
    
				     /**
				      * Access the matrix at position
				      * <i>i</i> for read and write
				      * access.
				      */
    MATRIX& matrix(unsigned int i);


};


/**
 * A vector of MGLevelObject<MatrixBlock>, which is implemented using shared
 * pointers, in order to allow for copying and rearranging. Each
 * matrix block can be identified by name.
 *
 * @relates MatrixBlock
 * @ingroup vector_valued
 * @author Baerbel Janssen, Guido Kanschat, 2010
 */
template <class MATRIX>
class MGMatrixBlockVector : public NamedData<std_cxx1x::shared_ptr<MatrixBlock<MATRIX> > >
{
  public:
				     /**
				      * Add a new matrix block at the
				      * position <tt>(row,column)</tt>
				      * in the block system. The third
				      * argument allows to give the
				      * matrix a name for later
				      * identification.
				      *
				      * @deprecated The final
				      * argument is ignored and will
				      * be removed in a future release.
				      */
    void add(unsigned int row, unsigned int column,
	     const std::string& name,
	     const BlockIndices* block_indices = 0);
    
				     /**
				      * Access a constant reference to
				      * the block at position <i>i</i>.
				      */
    const MatrixBlock<MATRIX>& block(unsigned int i) const;
    
				     /**
				      * Access a reference to
				      * the block at position <i>i</i>.
				      */
    MatrixBlock<MATRIX>& block(unsigned int i);
    
				     /**
				      * Access the matrix at position
				      * <i>i</i> for read and write
				      * access.
				      */
    MATRIX& matrix(unsigned int i);


};


//----------------------------------------------------------------------//

namespace internal
{
  template <class MATRIX>
  void
  reinit(MatrixBlockVector<MATRIX>& v, const BlockSparsityPattern& p)
  {
    Assert(false, ExcNotImplemented());
  }
  
  template <typename number>
  void
  reinit(MatrixBlockVector<dealii::SparseMatrix<number> >& v, const BlockSparsityPattern& p)
  {
    for (unsigned int i=0;i<v.size();++i)
      v(i)->matrix.reinit(p.block(v(i)->row, v(i)->column));
  }
  
}

template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock()
		:
		row(deal_II_numbers::invalid_unsigned_int),
		column(deal_II_numbers::invalid_unsigned_int)
{}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock(const MatrixBlock<MATRIX>& M)
		:
		row(M.row),
		column(M.column),
		matrix(M.matrix),
		block_indices(M.block_indices)
{}


template <class MATRIX>
inline
MatrixBlock<MATRIX>::MatrixBlock(unsigned int i, unsigned int j,
				 const BlockIndices* block_indices)
		:
		row(i), column(j), block_indices(block_indices)
{}


template <class MATRIX>
inline void
MatrixBlock<MATRIX>::add (
  const unsigned int gi,
  const unsigned int gj,
  const typename MATRIX::value_type value)
{
  Assert(block_indices != 0, ExcNotInitialized());

  const std::pair<unsigned int, unsigned int> bi
    = block_indices->global_to_local(gi);
  const std::pair<unsigned int, unsigned int> bj
    = block_indices->global_to_local(gj);

  Assert (bi.first == row, ExcBlockIndexMismatch(bi.first, row));
  Assert (bj.first == column, ExcBlockIndexMismatch(bj.first, column));

  matrix.add(bi.second, bj.second, value);
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const std::vector<unsigned int>&         row_indices,
				  const std::vector<unsigned int>& col_indices,
				  const FullMatrix<number>&        values,
				  const bool                       elide_zero_values)
{
  Assert(block_indices != 0, ExcNotInitialized());

  AssertDimension (row_indices.size(), values.m());
  AssertDimension (col_indices.size(), values.n());

  for (unsigned int i=0; i<row_indices.size(); ++i)
    add (row_indices[i], col_indices.size(), &col_indices[0], &values(i,0),
	 elide_zero_values);
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const unsigned int   b_row,
			  const unsigned int   n_cols,
			  const unsigned int  *col_indices,
			  const number        *values,
			  const bool,
			  const bool)
{
  Assert(block_indices != 0, ExcNotInitialized());
  const std::pair<unsigned int, unsigned int> bi
    = block_indices->global_to_local(b_row);

				   // In debug mode, we check whether
				   // all indices are in the correct
				   // block.

				   // Actually, for the time being, we
				   // leave it at this. While it may
				   // not be the most efficient way,
				   // it is at least thread safe.
//#ifdef DEBUG
  Assert(bi.first == row, ExcBlockIndexMismatch(bi.first, row));

  for (unsigned int j=0;j<n_cols;++j)
    {
      const std::pair<unsigned int, unsigned int> bj
	= block_indices->global_to_local(col_indices[j]);
      Assert(bj.first == column, ExcBlockIndexMismatch(bj.first, column));

      matrix.add(bi.second, bj.second, values[j]);
    }
//#endif
}


template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const std::vector<unsigned int> &indices,
			  const FullMatrix<number>        &values,
			  const bool                       elide_zero_values)
{
  Assert(block_indices != 0, ExcNotInitialized());

  AssertDimension (indices.size(), values.m());
  Assert (values.n() == values.m(), ExcNotQuadratic());

  for (unsigned int i=0; i<indices.size(); ++i)
    add (indices[i], indices.size(), &indices[0], &values(i,0),
	 elide_zero_values);
}



template <class MATRIX>
template <typename number>
inline
void
MatrixBlock<MATRIX>::add (const unsigned int               row,
			  const std::vector<unsigned int> &col_indices,
			  const std::vector<number>       &values,
			  const bool                       elide_zero_values)
{
  Assert(block_indices != 0, ExcNotInitialized());
  AssertDimension (col_indices.size(), values.size());
  add (row, col_indices.size(), &col_indices[0], &values[0],
       elide_zero_values);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::vmult (VECTOR& w, const VECTOR& v) const
{
  matrix.vmult(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::vmult_add (VECTOR& w, const VECTOR& v) const
{
  matrix.vmult_add(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::Tvmult (VECTOR& w, const VECTOR& v) const
{
  matrix.Tvmult(w,v);
}


template <class MATRIX>
template <class VECTOR>
inline
void
MatrixBlock<MATRIX>::Tvmult_add (VECTOR& w, const VECTOR& v) const
{
  matrix.Tvmult_add(w,v);
}



//----------------------------------------------------------------------//

template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::add(
  unsigned int row, unsigned int column,
  const std::string& name,
  const BlockIndices* block_indices)
{
  std_cxx1x::shared_ptr<MatrixBlock<MATRIX> > p(new MatrixBlock<MATRIX>(row, column, block_indices));
  NamedData<std_cxx1x::shared_ptr<MatrixBlock<MATRIX> > >::add(p, name);
}


template <class MATRIX>
inline void
MatrixBlockVector<MATRIX>::reinit(const BlockSparsityPattern& sparsity)
{
  internal::reinit(*this, sparsity);
}



template <class MATRIX>
inline const MatrixBlock<MATRIX>&
MatrixBlockVector<MATRIX>::block(unsigned int i) const
{
  return *this->read(i);
}


template <class MATRIX>
inline MatrixBlock<MATRIX>&
MatrixBlockVector<MATRIX>::block(unsigned int i)
{
  return *(*this)(i);
}


template <class MATRIX>
inline MATRIX&
MatrixBlockVector<MATRIX>::matrix(unsigned int i)
{
  return (*this)(i)->matrix;
}



DEAL_II_NAMESPACE_CLOSE

#endif
