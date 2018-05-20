// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2018 by the deal.II authors
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

#ifndef dealii_block_matrix_array_h
#define dealii_block_matrix_array_h

#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/table.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/pointer_matrix.h>
#include <deal.II/lac/vector_memory.h>

#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/*! @addtogroup Matrix2
 *@{
 */

/**
 * Block matrix composed of different single matrices; these matrices may even
 * be of different types.
 *
 * Given a set of arbitrary matrices <i>A<sub>i</sub></i>, this class
 * implements a block matrix with block entries of the form <i>M<sub>jk</sub>
 * = s<sub>jk</sub>A<sub>i</sub></i>.  Each <i>A<sub>i</sub></i> may be used
 * several times with different prefix. The matrices are not copied into the
 * BlockMatrixArray object, but rather a PointerMatrix referencing each of
 * them will be stored along with factors and transposition flags.
 *
 * Non-zero entries are registered by the function enter(), zero entries are
 * not stored at all. Using enter() with the same location <tt>(i,j)</tt>
 * several times will add the corresponding matrices in matrix-vector
 * multiplications. These matrices will not be actually added, but the
 * multiplications with them will be summed up.
 *
 * @note This mechanism makes it impossible to access single entries of
 * BlockMatrixArray. In particular, (block) relaxation preconditioners based
 * on PreconditionRelaxation or PreconditionBlock <b>cannot</b> be used with
 * this class. If you need a preconditioner for a BlockMatrixArray object, use
 * BlockTrianglePrecondition.
 *
 * <h3>Requirements on MatrixType</h3>
 *
 * The template argument <tt>MatrixType</tt> is a class providing the matrix-
 * vector multiplication functions vmult(), Tvmult(), vmult_add() and
 * Tvmult_add() used in this class, but with arguments of type
 * Vector&lt;number&gt; instead of BlockVector&lt;number&gt;. Every matrix
 * which can be used by PointerMatrix is allowed, in particular SparseMatrix
 * is a possible entry type.
 *
 * <h3>Example program</h3> We document the relevant parts of
 * <tt>examples/doxygen/block_matrix_array.cc</tt>.
 *
 * @dontinclude block_matrix_array.cc
 *
 * Obviously, we have to include the header file containing the definition of
 * BlockMatrixArray:
 * @skipline block_matrix_array.h
 *
 * First, we set up some matrices to be entered into the blocks.
 * @skip main
 * @until C.fill
 *
 * Now, we are ready to build a <i>2x2</i> BlockMatrixArray.
 * @line Block
 * First, we enter the matrix <tt>A</tt> multiplied by 2 in the upper left
 * block
 * @line enter
 * Now -1 times <tt>B1</tt> in the upper right block.
 * @line enter
 * We add the transpose of <tt>B2</tt> to the upper right block and continue
 * in a similar fashion. In the end, the block matrix structure is printed
 * into an LaTeX table.
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
 * The remaining code of this sample program concerns preconditioning and is
 * described in the documentation of BlockTrianglePrecondition.
 *
 * @deprecated This class has been superseded by BlockLinearOperator.
 *
 * @see
 * @ref GlossBlockLA "Block (linear algebra)"
 * @author Guido Kanschat
 * @date 2000-2005, 2010
 */
template <typename number          = double,
          typename BlockVectorType = BlockVector<number>>
class DEAL_II_DEPRECATED BlockMatrixArray : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Default constructor creating a useless object. initialize() must be
   * called before using it.
   */
  BlockMatrixArray();

  /**
   * Constructor fixing the dimensions.
   */
  BlockMatrixArray(const unsigned int n_block_rows,
                   const unsigned int n_block_cols);

  /**
   * Initialize object completely. This is the function to call for an object
   * created by the default constructor.
   */
  void
  initialize(const unsigned int n_block_rows, const unsigned int n_block_cols);

  /**
   * Adjust the matrix to a new size and delete all blocks.
   */
  void
  reinit(const unsigned int n_block_rows, const unsigned int n_block_cols);

  /**
   * Add a block matrix entry. The <tt>matrix</tt> is entered into a list of
   * blocks for multiplication, together with its coordinates <tt>row</tt> and
   * <tt>col</tt> as well as optional multiplication factor <tt>prefix</tt>
   * and transpose flag <tt>transpose</tt>.
   *
   * @note No check for consistency of block sizes is made. Therefore,
   * entering a block of wrong dimension here will only lead to a
   * ExcDimensionMismatch in one of the multiplication functions.
   */
  template <typename MatrixType>
  void
  enter(const MatrixType&  matrix,
        const unsigned int row,
        const unsigned int col,
        const number       prefix    = 1.,
        const bool         transpose = false);

  /**
   * Delete all entries, i.e. reset the matrix to an empty state.
   */
  void
  clear();

  /**
   * Number of block-entries per column.
   */
  unsigned int
  n_block_rows() const;

  /**
   * Number of block-entries per row.
   */
  unsigned int
  n_block_cols() const;

  /**
   * Matrix-vector multiplication.
   */
  void
  vmult(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Matrix-vector multiplication adding to <tt>dst</tt>.
   */
  void
  vmult_add(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Transposed matrix-vector multiplication.
   */
  void
  Tvmult(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Transposed matrix-vector multiplication adding to <tt>dst</tt>.
   */
  void
  Tvmult_add(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Matrix scalar product between two vectors (at least for a symmetric
   * matrix).
   */
  number
  matrix_scalar_product(const BlockVectorType& u,
                        const BlockVectorType& v) const;

  /**
   * Compute $u^T M u$. This is the square of the norm induced by the matrix
   * assuming the matrix is symmetric positive definitive.
   */
  number
  matrix_norm_square(const BlockVectorType& u) const;

  /**
   * Print the block structure as a LaTeX-array. This output will not be very
   * intuitive, since the current object lacks knowledge about what the
   * individual blocks represent or how they should be named. Instead, what
   * you will see is an entry for each block showing all the matrices with
   * their multiplication factors and possibly transpose marks. The matrices
   * itself are named successively as they are encountered. If the same matrix
   * is entered several times, it will be listed with different names.
   *
   * As an example, consider the following code:
   * @code
   *   FullMatrix<double> A1(4,4);
   *  FullMatrix<double> A2(4,4);
   *  FullMatrix<double> B(4,3);
   *  FullMatrix<double> C(3,3);
   *
   *  BlockMatrixArray<double> block(2,2);
   *
   *  block.enter(A1,0,0);
   *  block.enter(A2,0,0,2,true);
   *  block.enter(B,0,1,-3.);
   *  block.enter(B,0,1,-3.,true);
   *  block.enter(C,1,1,1.,true);
   *
   *  block.print_latex(std::cout);
   * @endcode
   * The current function will then produce output of the following kind:
   * @code
   * \begin{array}{cc}
   *    M0+2xM1^T &     -3xM2-3xM3^T\\
   *    &      M4^T
   * \end{array}
   * @endcode
   * Note how the individual blocks here are just numbered successively as
   * <code>M0</code> to <code>M4</code> and that the output misses the fact
   * that <code>M2</code> and <code>M3</code> are, in fact, the same matrix.
   * Nevertheless, the output at least gives some kind of idea of the block
   * structure of this matrix.
   */
  template <class StreamType>
  void
  print_latex(StreamType& out) const;

protected:
  /**
   * Internal data structure.
   *
   * For each entry of a BlockMatrixArray, its position, matrix, prefix and
   * optional transposition must be stored. This structure encapsulates all of
   * them.
   *
   * @author Guido Kanschat, 2000, 2001
   */
  class Entry
  {
  public:
    /**
     * Constructor initializing all data fields. A PointerMatrix object is
     * generated for <tt>matrix</tt>.
     */
    template <typename MatrixType>
    Entry(const MatrixType& matrix,
          size_type         row,
          size_type         col,
          number            prefix,
          bool              transpose);

    /**
     * Copy constructor invalidating the old object. Since it is only used for
     * entering temporary objects into a vector, this is ok.
     *
     * For a deep copy, we would need a reproduction operator in
     * PointerMatixBase.
     */
    Entry(const Entry&);

    /**
     * Destructor, where we delete the PointerMatrix created by the
     * constructor.
     */
    ~Entry();

    /**
     * Row number in the block matrix.
     */
    size_type row;

    /**
     * Column number in the block matrix.
     */
    size_type col;

    /**
     * Factor in front of the matrix block.
     */
    number prefix;

    /**
     * Indicates that matrix block must be transposed for multiplication.
     */
    bool transpose;

    /**
     * The matrix block itself.
     */
    PointerMatrixBase<typename BlockVectorType::BlockType>* matrix;

    /**
     * Assignment operator.
     *
     * @note Since the copy constructor is destructive (see its documentation)
     * and only exists for convenience there is no reasonable way to implement
     * this, so it is explicitly deleted.
     */
    Entry&
    operator=(const Entry&)
      = delete;
  };

  /**
   * Array of block entries in the matrix.
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
};

/*@}*/

/**
 * Inversion of a block-triangular matrix.
 *
 * In this block matrix, the inverses of the diagonal blocks are stored
 * together with the off-diagonal blocks of a block matrix. Then, forward or
 * backward insertion is performed block-wise. The diagonal blocks are NOT
 * inverted for this purpose!
 *
 * Like for all preconditioners, the preconditioning operation is performed by
 * the vmult() member function.
 *
 * @note While block indices may be duplicated (see BlockMatrixArray) to add
 * blocks, this has to be used with caution, since summing up the inverse of
 * two blocks does not yield the inverse of the sum. While the latter would be
 * desirable, we can only perform the first.
 *
 * The implementation may be a little clumsy, but it should be sufficient as
 * long as the block sizes are much larger than the number of blocks.
 *
 * <h3>Example</h3> Here, we document the second part of
 * <tt>examples/doxygen/block_matrix_array.cc</tt>. For the beginning of this
 * file, see BlockMatrixArray.
 *
 * In order to set up the preconditioner, we have to compute the inverses of
 * the diagonal blocks ourselves. Since we used FullMatrix objects, this is
 * fairly easy.
 * @dontinclude block_matrix_array.cc
 * @skip Error
 * @until Cinv.invert
 *
 * After creating a <i>2x2</i> BlockTrianglePrecondition object, we only fill
 * its diagonals. The scaling factor <i>1/2</i> used for <tt>A</tt> is the
 * reciprocal of the scaling factor used for the <tt>matrix</tt> itself.
 * Remember, this preconditioner actually <b>multiplies</b> with the diagonal
 * blocks.
 * @until Cinv
 *
 * Now, we have a block Jacobi preconditioner, which is still symmetric, since
 * the blocks are symmetric. Therefore, we can still use the preconditioned
 * conjugate gradient method.
 * @until Error
 *
 * Now, we enter the subdiagonal block. This is the same as in
 * <tt>matrix</tt>.
 * @until B2
 *
 * Since the preconditioner is not symmetric anymore, we use the GMRES method
 * for solving.
 * @until Error
 *
 * @deprecated This class has been superseded by block_back_substitution and
 * block_forward_substitution, which build on BlockLinearOperator.
 *
 * @ingroup Preconditioners
 * @author Guido Kanschat, 2001, 2005
 */
template <typename number          = double,
          typename BlockVectorType = BlockVector<number>>
class DEAL_II_DEPRECATED BlockTrianglePrecondition
  : private BlockMatrixArray<number, BlockVectorType>
{
public:
  /**
   * Declare type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * Default constructor creating a useless object. initialize() must be
   * called before using it.
   */
  BlockTrianglePrecondition();

  /**
   * Constructor. This matrix must be block-quadratic, and <tt>n_blocks</tt>
   * is the number of blocks in each direction.
   */
  BlockTrianglePrecondition(const unsigned int n_blocks);

  /**
   * Resize preconditioner to a new size and clear all blocks.
   */
  void
  reinit(const unsigned int n_block_rows);

  /**
   * Enter a block. This calls BlockMatrixArray::enter(). Remember that the
   * diagonal blocks should actually be inverse matrices or preconditioners.
   */
  template <typename MatrixType>
  void
  enter(const MatrixType& matrix,
        const size_type   row,
        const size_type   col,
        const number      prefix    = 1.,
        const bool        transpose = false);

  /**
   * Preconditioning.
   */
  void
  vmult(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Preconditioning adding to <tt>dst</tt>.
   */
  void
  vmult_add(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Transposed preconditioning
   */
  void
  Tvmult(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Transposed preconditioning adding to <tt>dst</tt>.
   */
  void
  Tvmult_add(BlockVectorType& dst, const BlockVectorType& src) const;

  /**
   * Make function of base class available.
   */
  using BlockMatrixArray<number, BlockVectorType>::print_latex;

  /**
   * Make function of base class available.
   */
  using BlockMatrixArray<number, BlockVectorType>::n_block_rows;

  /**
   * Make function of base class available.
   */
  using BlockMatrixArray<number, BlockVectorType>::n_block_cols;
  using BlockMatrixArray<number, BlockVectorType>::clear;
  using BlockMatrixArray<number, BlockVectorType>::Subscriptor::subscribe;
  using BlockMatrixArray<number, BlockVectorType>::Subscriptor::unsubscribe;

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Each diagonal block must contain one and only one matrix. If this
   * exception is thrown, you did not enter a matrix here.
   */
  DeclException1(ExcNoDiagonal,
                 size_type,
                 << "No diagonal entry was added for block " << arg1);

  /**
   * Each diagonal block must contain one and only one matrix. If this
   * exception is thrown, you entered a second matrix here.
   */
  DeclException1(ExcMultipleDiagonal,
                 size_type,
                 << "Inverse diagonal entries may not be added in block "
                 << arg1);
  //@}
private:
  /**
   * Add all off-diagonal contributions and return the entry of the diagonal
   * element for one row.
   */
  void
  do_row(BlockVectorType& dst, size_type row_num) const;

  /**
   * Flag for backward insertion.
   */
  bool backward;
};

#ifndef DOXYGEN
//---------------------------------------------------------------------------

template <typename number, typename BlockVectorType>
template <typename MatrixType>
inline BlockMatrixArray<number, BlockVectorType>::Entry::Entry(
  const MatrixType& m,
  size_type         row,
  size_type         col,
  number            prefix,
  bool              transpose)
  : row(row),
    col(col),
    prefix(prefix),
    transpose(transpose),
    matrix(new_pointer_matrix_base(m,
                                   typename BlockVectorType::BlockType(),
                                   typeid(*this).name()))
{}

template <typename number, typename BlockVectorType>
template <typename MatrixType>
inline void
BlockMatrixArray<number, BlockVectorType>::enter(const MatrixType& matrix,
                                                 unsigned int      row,
                                                 unsigned int      col,
                                                 number            prefix,
                                                 bool              transpose)
{
  Assert(row < n_block_rows(), ExcIndexRange(row, 0, n_block_rows()));
  Assert(col < n_block_cols(), ExcIndexRange(col, 0, n_block_cols()));
  entries.push_back(Entry(matrix, row, col, prefix, transpose));
}

template <typename number, typename BlockVectorType>
template <class StreamType>
inline void
BlockMatrixArray<number, BlockVectorType>::print_latex(StreamType& out) const
{
  out << "\\begin{array}{" << std::string(n_block_cols(), 'c') << "}"
      << std::endl;

  Table<2, std::string> array(n_block_rows(), n_block_cols());

  typedef std::map<
    const PointerMatrixBase<typename BlockVectorType::BlockType>*,
    std::string>
          NameMap;
  NameMap matrix_names;

  typename std::vector<Entry>::const_iterator m   = entries.begin();
  typename std::vector<Entry>::const_iterator end = entries.end();

  size_type matrix_number = 0;
  for(; m != end; ++m)
    {
      if(matrix_names.find(m->matrix) == matrix_names.end())
        {
          std::pair<typename NameMap::iterator, bool> x = matrix_names.insert(
            std::pair<
              const PointerMatrixBase<typename BlockVectorType::BlockType>*,
              std::string>(m->matrix, std::string("M")));
          std::ostringstream stream;
          stream << matrix_number++;

          x.first->second += stream.str();
        }

      std::ostringstream stream;

      if(array(m->row, m->col) != "" && m->prefix >= 0)
        stream << "+";
      if(m->prefix != 1.)
        stream << m->prefix << 'x';
      stream << matrix_names.find(m->matrix)->second;
      //      stream << '(' << m->matrix << ')';
      if(m->transpose)
        stream << "^T";

      array(m->row, m->col) += stream.str();
    }
  for(unsigned int i = 0; i < n_block_rows(); ++i)
    for(unsigned int j = 0; j < n_block_cols(); ++j)
      {
        out << '\t' << array(i, j);
        if(j == n_block_cols() - 1)
          {
            if(i != n_block_rows() - 1)
              out << "\\\\" << std::endl;
            else
              out << std::endl;
          }
        else
          out << " &";
      }
  out << "\\end{array}" << std::endl;
}

template <typename number, typename BlockVectorType>
template <typename MatrixType>
inline void
BlockTrianglePrecondition<number, BlockVectorType>::enter(
  const MatrixType& matrix,
  size_type         row,
  size_type         col,
  number            prefix,
  bool              transpose)
{
  BlockMatrixArray<number, BlockVectorType>::enter(
    matrix, row, col, prefix, transpose);
}

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
