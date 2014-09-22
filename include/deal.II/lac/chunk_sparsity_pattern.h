// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
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

#ifndef __deal2__chunk_sparsity_pattern_h
#define __deal2__chunk_sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/base/vector_slice.h>

#include <deal.II/lac/sparsity_pattern.h>

#include <vector>
#include <iostream>

DEAL_II_NAMESPACE_OPEN


template <typename> class ChunkSparseMatrix;


/*! @addtogroup Sparsity
 *@{
 */



/**
 * Iterators on sparsity patterns
 */
namespace ChunkSparsityPatternIterators
{
  // forward declaration
  class Iterator;

  /**
   * Accessor class for iterators into sparsity patterns. This class is
   * also the base class for both const and non-const accessor classes
   * into sparse matrices.
   *
   * Note that this class only allows read access to elements, providing
   * their row and column number. It does not allow modifying the
   * sparsity pattern itself.
   *
   * @author Martin Kronbichler
   * @date 2013
   */
  class Accessor
  {
  public:
    /**
     * Constructor.
     */
    Accessor (const ChunkSparsityPattern *matrix,
              const unsigned int          row);

    /**
     * Constructor. Construct the end accessor for the given sparsity pattern.
     */
    Accessor (const ChunkSparsityPattern *matrix);

    /**
     * Row number of the element represented by this object. This function can
     * only be called for entries for which is_valid_entry() is true.
     */
    unsigned int row () const;

    /**
     * Returns the global index from the reduced sparsity pattern.
     */
    std::size_t reduced_index() const;

    /**
     * Column number of the element represented by this object. This function
     * can only be called for entries for which is_valid_entry() is true.
     */
    unsigned int column () const;

    /**
     * Return whether the sparsity pattern entry pointed to by this iterator
     * is valid or not. Note that after compressing the sparsity pattern, all
     * entries are valid. However, before compression, the sparsity pattern
     * allocated some memory to be used while still adding new nonzero
     * entries; if you create iterators in this phase of the sparsity
     * pattern's lifetime, you will iterate over elements that are not
     * valid. If this is so, then this function will return false.
     */
    bool is_valid_entry () const;


    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool operator == (const Accessor &) const;


    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * sparsity pattern.
     */
    bool operator < (const Accessor &) const;

  protected:
    /**
     * The sparsity pattern we operate on accessed.
     */
    const ChunkSparsityPattern *sparsity_pattern;

    /**
     * The accessor of the (reduced) sparsity pattern.
     */
    SparsityPatternIterators::Accessor reduced_accessor;

    /**
     * Current chunk row number.
     */
    unsigned int chunk_row;

    /**
     * Current chunk col number.
     */
    unsigned int chunk_col;

    /**
     * Move the accessor to the next nonzero entry in the matrix.
     */
    void advance ();

    /**
     * Grant access to iterator class.
     */
    friend class Iterator;
  };



  /**
   * STL conforming iterator walking over the elements of a sparsity pattern.
   */
  class Iterator
  {
  public:
    /**
     * Constructor. Create an iterator into the sparsity pattern @p sp for the
     * given row and the index within it.
     */
    Iterator (const ChunkSparsityPattern *sp,
              const unsigned int          row);

    /**
     * Prefix increment.
     */
    Iterator &operator++ ();

    /**
     * Postfix increment.
     */
    Iterator operator++ (int);

    /**
     * Dereferencing operator.
     */
    const Accessor &operator* () const;

    /**
     * Dereferencing operator.
     */
    const Accessor *operator-> () const;

    /**
     * Comparison. True, if both iterators point to the same matrix position.
     */
    bool operator == (const Iterator &) const;

    /**
     * Inverse of <tt>==</tt>.
     */
    bool operator != (const Iterator &) const;

    /**
     * Comparison operator. Result is true if either the first row number is
     * smaller or if the row numbers are equal and the first index is smaller.
     *
     * This function is only valid if both iterators point into the same
     * matrix.
     */
    bool operator < (const Iterator &) const;

  private:
    /**
     * Store an object of the accessor class.
     */
    Accessor accessor;
  };
}



/**
 * Structure representing the sparsity pattern of a sparse matrix.
 * This class is an example of the "static" type of @ref Sparsity.
 * It uses the compressed row storage (CSR) format to store data.
 *
 * The use of this class is demonstrated in step-51.
 *
 * @author Wolfgang Bangerth, 2008
 */
class ChunkSparsityPattern : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;
  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   */
  typedef ChunkSparsityPatternIterators::Iterator const_iterator;

  /**
   * Typedef an iterator class that allows to walk over all nonzero elements
   * of a sparsity pattern.
   *
   * Since the iterator does not allow to modify the sparsity pattern, this
   * type is the same as that for @p const_iterator.
   */
  typedef ChunkSparsityPatternIterators::Iterator iterator;

  /**
   * Define a value which is used to indicate that a certain value in the
   * colnums array is unused, i.e. does not represent a certain column number
   * index.
   *
   * Indices with this invalid value are used to insert new entries to the
   * sparsity pattern using the add() member function, and are removed when
   * calling compress().
   *
   * You should not assume that the variable declared here has a certain
   * value. The initialization is given here only to enable the compiler to
   * perform some optimizations, but the actual value of the variable may
   * change over time.
   */
  static const size_type invalid_entry = SparsityPattern::invalid_entry;

  /**
   * Initialize the matrix empty, that is with no memory allocated. This is
   * useful if you want such objects as member variables in other classes. You
   * can make the structure usable by calling the reinit() function.
   */
  ChunkSparsityPattern ();

  /**
   * Copy constructor. This constructor is only allowed to be called if the
   * matrix structure to be copied is empty. This is so in order to prevent
   * involuntary copies of objects for temporaries, which can use large
   * amounts of computing time.  However, copy constructors are needed if yo
   * want to use the STL data types on classes like this, e.g. to write such
   * statements like <tt>v.push_back (ChunkSparsityPattern());</tt>, with
   * <tt>v</tt> a vector of ChunkSparsityPattern objects.
   *
   * Usually, it is sufficient to use the explicit keyword to disallow
   * unwanted temporaries, but for the STL vectors, this does not work. Since
   * copying a structure like this is not useful anyway because multiple
   * matrices can use the same sparsity structure, copies are only allowed for
   * empty objects, as described above.
   */
  ChunkSparsityPattern (const ChunkSparsityPattern &);

  /**
   * Initialize a rectangular matrix.
   *
   * @arg m number of rows
   * @arg n number of columns
   * @arg max_per_row maximum number of nonzero entries per row
   */
  ChunkSparsityPattern (const size_type m,
                        const size_type n,
                        const size_type max_chunks_per_row,
                        const size_type chunk_size);

  /**
   * @deprecated This constructor is deprecated. Use the version
   * without the last argument
   */
  ChunkSparsityPattern (const size_type m,
                        const size_type n,
                        const size_type max_chunks_per_row,
                        const size_type chunk_size,
                        const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Initialize a rectangular matrix.
   *
   * @arg m number of rows
   * @arg n number of columns
   * @arg row_lengths possible number of nonzero entries for each row.  This
   * vector must have one entry for each row.
   */
  ChunkSparsityPattern (const size_type m,
                        const size_type n,
                        const std::vector<size_type> &row_lengths,
                        const size_type chunk_size);

  /**
   * @deprecated This constructor is deprecated. Use the version
   * without the last argument
   */
  ChunkSparsityPattern (const size_type               m,
                        const size_type               n,
                        const std::vector<size_type> &row_lengths,
                        const size_type chunk_size,
                        const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Initialize a quadratic matrix of dimension <tt>n</tt> with at most
   * <tt>max_per_row</tt> nonzero entries per row.
   *
   * This constructor automatically enables optimized storage of diagonal
   * elements. To avoid this, use the constructor taking row and column
   * numbers separately.
   */
  ChunkSparsityPattern (const size_type n,
                        const size_type max_per_row,
                        const size_type chunk_size);

  /**
   * Initialize a quadratic matrix.
   *
   * @arg m number of rows and columns
   * @arg row_lengths possible number of nonzero entries for each row.  This
   * vector must have one entry for each row.
   */
  ChunkSparsityPattern (const size_type                m,
                        const std::vector<size_type>  &row_lengths,
                        const size_type                chunk_size);

  /**
   * @deprecated This constructor is deprecated. Use the version
   * without the last argument
   */
  ChunkSparsityPattern (const size_type               m,
                        const std::vector<size_type> &row_lengths,
                        const size_type               chunk_size,
                        const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Destructor.
   */
  ~ChunkSparsityPattern ();

  /**
   * Copy operator. For this the same holds as for the copy constructor: it is
   * declared, defined and fine to be called, but the latter only for empty
   * objects.
   */
  ChunkSparsityPattern &operator = (const ChunkSparsityPattern &);

  /**
   * Reallocate memory and set up data structures for a new matrix with <tt>m
   * </tt>rows and <tt>n</tt> columns, with at most <tt>max_per_row</tt>
   * nonzero entries per row.
   *
   * This function simply maps its operations to the other <tt>reinit</tt>
   * function.
   */
  void reinit (const size_type m,
               const size_type n,
               const size_type max_per_row,
               const size_type chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  void reinit (const size_type m,
               const size_type n,
               const size_type max_per_row,
               const size_type chunk_size,
               const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Reallocate memory for a matrix of size <tt>m x n</tt>. The number of
   * entries for each row is taken from the array <tt>row_lengths</tt> which
   * has to give this number of each row <tt>i=1...m</tt>.
   *
   * If <tt>m*n==0</tt> all memory is freed, resulting in a total
   * reinitialization of the object. If it is nonzero, new memory is only
   * allocated if the new size extends the old one. This is done to save time
   * and to avoid fragmentation of the heap.
   *
   * If the number of rows equals
   * the number of columns then
   * diagonal elements are stored
   * first in each row to allow
   * optimized access in relaxation
   * methods of SparseMatrix.
   */
  void reinit (const size_type m,
               const size_type n,
               const std::vector<size_type> &row_lengths,
               const size_type chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  void reinit (const size_type                m,
               const size_type                n,
               const std::vector<size_type > &row_lengths,
               const size_type                chunk_size,
               const bool                     optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Same as above, but with a VectorSlice argument instead.
   */
  void reinit (const size_type m,
               const size_type n,
               const VectorSlice<const std::vector<size_type> > &row_lengths,
               const size_type chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  void reinit (const size_type m,
               const size_type n,
               const VectorSlice<const std::vector<size_type> > &row_lengths,
               const size_type chunk_size,
               const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * This function compresses the sparsity structure that this object
   * represents.  It does so by eliminating unused entries and sorting the
   * remaining ones to allow faster access by usage of binary search
   * algorithms. A special sorting scheme is used for the diagonal entry of
   * quadratic matrices, which is always the first entry of each row.
   *
   * The memory which is no more needed is released.
   *
   * SparseMatrix objects require the ChunkSparsityPattern objects they are
   * initialized with to be compressed, to reduce memory requirements.
   */
  void compress ();

  /**
   * This function can be used as a replacement for reinit(), subsequent calls
   * to add() and a final call to close() if you know exactly in advance the
   * entries that will form the matrix sparsity pattern.
   *
   * The first two parameters determine the size of the matrix. For the two
   * last ones, note that a sparse matrix can be described by a sequence of
   * rows, each of which is represented by a sequence of pairs of column
   * indices and values. In the present context, the begin() and end()
   * parameters designate iterators (of forward iterator type) into a
   * container, one representing one row. The distance between begin() and
   * end() should therefore be equal to n_rows(). These iterators may be
   * iterators of <tt>std::vector</tt>, <tt>std::list</tt>, pointers into a
   * C-style array, or any other iterator satisfying the requirements of a
   * forward iterator. The objects pointed to by these iterators (i.e. what we
   * get after applying <tt>operator*</tt> or <tt>operator-></tt> to one of
   * these iterators) must be a container itself that provides functions
   * <tt>begin</tt> and <tt>end</tt> designating a range of iterators that
   * describe the contents of one line. Dereferencing these inner iterators
   * must either yield a pair of an unsigned integer as column index and a
   * value of arbitrary type (such a type would be used if we wanted to
   * describe a sparse matrix with one such object), or simply an unsigned
   * integer (of we only wanted to describe a sparsity pattern). The function
   * is able to determine itself whether an unsigned integer or a pair is what
   * we get after dereferencing the inner iterators, through some template
   * magic.
   *
   * While the order of the outer iterators denotes the different rows of the
   * matrix, the order of the inner iterator denoting the columns does not
   * matter, as they are sorted internal to this function anyway.
   *
   * Since that all sounds very complicated, consider the following example
   * code, which may be used to fill a sparsity pattern:
   * @code
   * std::vector<std::vector<size_type> > column_indices (n_rows);
   * for (size_type row=0; row<n_rows; ++row)
   *         // generate necessary columns in this row
   *   fill_row (column_indices[row]);
   *
   * sparsity.copy_from (n_rows, n_cols,
   *                     column_indices.begin(),
   *                     column_indices.end());
   * @endcode
   *
   * Note that this example works since the iterators dereferenced yield
   * containers with functions <tt>begin</tt> and <tt>end</tt> (namely
   * <tt>std::vector</tt>s), and the inner iterators dereferenced yield
   * unsigned integers as column indices. Note that we could have replaced
   * each of the two <tt>std::vector</tt> occurrences by <tt>std::list</tt>,
   * and the inner one by <tt>std::set</tt> as well.
   *
   * Another example would be as follows, where we initialize a whole matrix,
   * not only a sparsity pattern:
   * @code
   * std::vector<std::map<size_type,double> > entries (n_rows);
   * for (size_type row=0; row<n_rows; ++row)
   *         // generate necessary pairs of columns
   *         // and corresponding values in this row
   *   fill_row (entries[row]);
   *
   * sparsity.copy_from (n_rows, n_cols,
   *                     column_indices.begin(),
   *                     column_indices.end());
   * matrix.reinit (sparsity);
   * matrix.copy_from (column_indices.begin(),
   *                   column_indices.end());
   * @endcode
   *
   * This example works because dereferencing iterators of the inner type
   * yields a pair of unsigned integers and a value, the first of which we
   * take as column index. As previously, the outer <tt>std::vector</tt> could
   * be replaced by <tt>std::list</tt>, and the inner <tt>std::map<unsigned
   * int,double></tt> could be replaced by <tt>std::vector<std::pair<unsigned
   * int,double> ></tt>, or a list or set of such pairs, as they all return
   * iterators that point to such pairs.
   */
  template <typename ForwardIterator>
  void copy_from (const size_type n_rows,
                  const size_type n_cols,
                  const ForwardIterator begin,
                  const ForwardIterator end,
                  const size_type chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  template <typename ForwardIterator>
  void copy_from (const size_type       n_rows,
                  const size_type       n_cols,
                  const ForwardIterator begin,
                  const ForwardIterator end,
                  const size_type       chunk_size,
                  const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Copy data from an object of type CompressedSparsityPattern,
   * CompressedSetSparsityPattern or CompressedSimpleSparsityPattern.
   * Previous content of this object is lost, and the sparsity pattern is in
   * compressed mode afterwards.
   */
  template <typename SparsityType>
  void copy_from (const SparsityType &csp,
                  const size_type     chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  template <typename SparsityType>
  void copy_from (const SparsityType &csp,
                  const size_type     chunk_size,
                  const bool          optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Take a full matrix and use its nonzero entries to generate a sparse
   * matrix entry pattern for this object.
   *
   * Previous content of this object is lost, and the sparsity pattern is in
   * compressed mode afterwards.
   */
  template <typename number>
  void copy_from (const FullMatrix<number> &matrix,
                  const size_type chunk_size);

  /**
   * @deprecated This function is deprecated. Use the function
   * without the last argument
   */
  template <typename number>
  void copy_from (const FullMatrix<number> &matrix,
                  const size_type chunk_size,
                  const bool optimize_diagonal) DEAL_II_DEPRECATED;

  /**
   * Set the sparsity pattern of the chunk sparsity pattern to be given by
   * <tt>chunk_size*chunksize</tt> blocks of the sparsity pattern for chunks
   * specified. Note that the final number of rows <tt>m</tt> of the sparsity
   * pattern will be approximately <tt>sparsity_pattern_for_chunks.n_rows() *
   * chunk_size</tt> (modulo padding elements in the last chunk) and similarly
   * for the number of columns <tt>n</tt>.
   *
   * This is a special initialization option in case you can tell the position
   * of the chunk already from the beginning without generating the sparsity
   * pattern using <tt>make_sparsity_pattern</tt> calls. This bypasses the
   * search for chunks but of course needs to be handled with care in order to
   * give a correct sparsity pattern.
   *
   * Previous content of this object is lost, and the sparsity pattern is in
   * compressed mode afterwards.
   */
  template <typename Sparsity>
  void create_from (const unsigned int  m,
                    const unsigned int  n,
                    const Sparsity     &sparsity_pattern_for_chunks,
                    const unsigned int  chunk_size,
                    const bool          optimize_diagonal = true);

  /**
   * Return whether the object is empty. It is empty if no memory is
   * allocated, which is the same as that both dimensions are zero.
   */
  bool empty () const;

  /**
   * Return the chunk size given as argument when constructing this object.
   */
  size_type get_chunk_size () const;

  /**
   * Return the maximum number of entries per row. Before compression, this
   * equals the number given to the constructor, while after compression, it
   * equals the maximum number of entries actually allocated by the user.
   */
  size_type max_entries_per_row () const;

  /**
   * Add a nonzero entry to the matrix. This function may only be called for
   * non-compressed sparsity patterns.
   *
   * If the entry already exists, nothing bad happens.
   */
  void add (const size_type i,
            const size_type j);

  /**
   * Make the sparsity pattern symmetric by adding the sparsity pattern of the
   * transpose object.
   *
   * This function throws an exception if the sparsity pattern does not
   * represent a quadratic matrix.
   */
  void symmetrize ();

  /**
   * Return number of rows of this matrix, which equals the dimension of the
   * image space.
   */
  inline size_type n_rows () const;

  /**
   * Return number of columns of this matrix, which equals the dimension of
   * the range space.
   */
  inline size_type n_cols () const;

  /**
   * Check if a value at a certain position may be non-zero.
   */
  bool exists (const size_type i,
               const size_type j) const;

  /**
   * Number of entries in a specific row.
   */
  size_type row_length (const size_type row) const;

  /**
   * Compute the bandwidth of the matrix represented by this structure. The
   * bandwidth is the maximum of $|i-j|$ for which the index pair $(i,j)$
   * represents a nonzero entry of the matrix. Consequently, the maximum
   * bandwidth a $n\times m$ matrix can have is $\max\{n-1,m-1\}$.
   */
  size_type bandwidth () const;

  /**
   * Return the number of nonzero elements of this matrix. Actually, it
   * returns the number of entries in the sparsity pattern; if any of the
   * entries should happen to be zero, it is counted anyway.
   *
   * This function may only be called if the matrix struct is compressed. It
   * does not make too much sense otherwise anyway.
   */
  size_type n_nonzero_elements () const;

  /**
   * Return whether the structure is compressed or not.
   */
  bool is_compressed () const;

  /**
   * Determine whether the matrix uses the special convention for quadratic
   * matrices that the diagonal entries are stored first in each row.
   *
   * @deprecated The function always returns true if the matrix is
   * square and false if it is not.
   */
  bool optimize_diagonal () const DEAL_II_DEPRECATED;

  /**
   * Return whether this object stores only those entries that have been added
   * explicitly, or if the sparsity pattern contains elements that have been
   * added through other means (implicitly) while building it. For the current
   * class, the result is true if and only if it is square because it then
   * unconditionally stores the diagonal entries whether they have been added explicitly or not.
   *
   * This function mainly serves the purpose of describing the current class
   * in cases where several kinds of sparsity patterns can be passed as
   * template arguments.
   */
  bool stores_only_added_elements () const;

  /**
   * STL-like iterator with the first entry of the matrix. The resulting
   * iterator can be used to walk over all nonzero entries of the sparsity
   * pattern.
   */
  iterator begin () const;

  /**
   * Final iterator.
   */
  iterator end () const;

  /**
   * STL-like iterator with the first entry of row <tt>r</tt>.
   *
   * Note that if the given row is empty, i.e. does not contain any nonzero
   * entries, then the iterator returned by this function equals
   * <tt>end(r)</tt>. Note also that the iterator may not be dereferencable in
   * that case.
   */
  iterator begin (const unsigned int r) const;

  /**
   * Final iterator of row <tt>r</tt>. It points to the first element past the
   * end of line @p r, or past the end of the entire sparsity pattern.
   *
   * Note that the end iterator is not necessarily dereferencable. This is in
   * particular the case if it is the end iterator for the last row of a
   * matrix.
   */
  iterator end (const unsigned int r) const;

  /**
   * Write the data of this object en bloc to a file. This is done in a binary
   * mode, so the output is neither readable by humans nor (probably) by other
   * computers using a different operating system of number format.
   *
   * The purpose of this function is that you can swap out matrices and
   * sparsity pattern if you are short of memory, want to communicate between
   * different programs, or allow objects to be persistent across different
   * runs of the program.
   */
  void block_write (std::ostream &out) const;

  /**
   * Read data that has previously been written by block_write() from a
   * file. This is done using the inverse operations to the above function, so
   * it is reasonably fast because the bitstream is not interpreted except for
   * a few numbers up front.
   *
   * The object is resized on this operation, and all previous contents are
   * lost.
   *
   * A primitive form of error checking is performed which will recognize the
   * bluntest attempts to interpret some data as a vector stored bitwise to a
   * file, but not more.
   */
  void block_read (std::istream &in);

  /**
   * Print the sparsity of the matrix. The output consists of one line per row
   * of the format <tt>[i,j1,j2,j3,...]</tt>. <i>i</i> is the row number and
   * <i>jn</i> are the allocated columns in this row.
   */
  void print (std::ostream &out) const;

  /**
   * Print the sparsity of the matrix in a format that <tt>gnuplot</tt>
   * understands and which can be used to plot the sparsity pattern in a
   * graphical way. The format consists of pairs <tt>i j</tt> of nonzero
   * elements, each representing one entry of this matrix, one per line of the
   * output file. Indices are counted from zero on, as usual. Since sparsity
   * patterns are printed in the same way as matrices are displayed, we print
   * the negative of the column index, which means that the <tt>(0,0)</tt>
   * element is in the top left rather than in the bottom left corner.
   *
   * Print the sparsity pattern in gnuplot by setting the data style to dots
   * or points and use the <tt>plot</tt> command.
   */
  void print_gnuplot (std::ostream &out) const;

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object. See MemoryConsumption.
   */
  std::size_t memory_consumption () const;

  /** @addtogroup Exceptions
   * @{ */
  /**
   * Exception
   */
  DeclException1 (ExcInvalidNumber,
                  int,
                  << "The provided number is invalid here: " << arg1);
  /**
   * Exception
   */
  DeclException2 (ExcInvalidIndex,
                  int, int,
                  << "The given index " << arg1
                  << " should be less than " << arg2 << ".");
  /**
   * Exception
   */
  DeclException2 (ExcNotEnoughSpace,
                  int, int,
                  << "Upon entering a new entry to row " << arg1
                  << ": there was no free entry any more. " << std::endl
                  << "(Maximum number of entries for this row: "
                  << arg2 << "; maybe the matrix is already compressed?)");
  /**
   * Exception
   */
  DeclException0 (ExcNotCompressed);
  /**
   * Exception
   */
  DeclException0 (ExcMatrixIsCompressed);
  /**
   * Exception
   */
  DeclException0 (ExcEmptyObject);
  /**
   * Exception
   */
  DeclException0 (ExcInvalidConstructorCall);
  /**
   * Exception
   */
  DeclException2 (ExcIteratorRange,
                  int, int,
                  << "The iterators denote a range of " << arg1
                  << " elements, but the given number of rows was " << arg2);
  /**
   * Exception
   */
  DeclException0 (ExcMETISNotInstalled);
  /**
   * Exception
   */
  DeclException1 (ExcInvalidNumberOfPartitions,
                  int,
                  << "The number of partitions you gave is " << arg1
                  << ", but must be greater than zero.");
  /**
   * Exception
   */
  DeclException2 (ExcInvalidArraySize,
                  int, int,
                  << "The array has size " << arg1 << " but should have size "
                  << arg2);
  //@}
private:
  /**
   * Number of rows that this sparsity structure shall represent.
   */
  size_type rows;

  /**
   * Number of columns that this sparsity structure shall represent.
   */
  size_type cols;

  /**
   * The size of chunks.
   */
  size_type chunk_size;

  /**
   * The reduced sparsity pattern. We store only which chunks exist, with each
   * chunk a block in the matrix of size chunk_size by chunk_size.
   */
  SparsityPattern sparsity_pattern;

  /**
   * Make all the chunk sparse matrix kinds friends.
   */
  template <typename> friend class ChunkSparseMatrix;

  /**
   * Make the accessor class a friend.
   */
  friend class ChunkSparsityPatternIterators::Accessor;
};


/*@}*/
/*---------------------- Inline functions -----------------------------------*/

#ifndef DOXYGEN

namespace ChunkSparsityPatternIterators
{
  inline
  Accessor::
  Accessor (const ChunkSparsityPattern *sparsity_pattern,
            const unsigned int          row)
    :
    sparsity_pattern(sparsity_pattern),
    reduced_accessor(row==sparsity_pattern->n_rows() ?
                     *sparsity_pattern->sparsity_pattern.end() :
                     *sparsity_pattern->sparsity_pattern.
                     begin(row/sparsity_pattern->get_chunk_size())),
    chunk_row (row==sparsity_pattern->n_rows() ? 0 :
               row%sparsity_pattern->get_chunk_size()),
    chunk_col (0)
  {}



  inline
  Accessor::
  Accessor (const ChunkSparsityPattern *sparsity_pattern)
    :
    sparsity_pattern(sparsity_pattern),
    reduced_accessor(*sparsity_pattern->sparsity_pattern.end()),
    chunk_row (0),
    chunk_col (0)
  {}



  inline
  bool
  Accessor::is_valid_entry () const
  {
    return reduced_accessor.is_valid_entry()
           &&
           sparsity_pattern->get_chunk_size()*reduced_accessor.row()+chunk_row <
           sparsity_pattern->n_rows()
           &&
           sparsity_pattern->get_chunk_size()*reduced_accessor.column()+chunk_col <
           sparsity_pattern->n_cols();
  }



  inline
  unsigned int
  Accessor::row() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    return sparsity_pattern->get_chunk_size()*reduced_accessor.row() +
           chunk_row;
  }



  inline
  unsigned int
  Accessor::column() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    return sparsity_pattern->get_chunk_size()*reduced_accessor.column() +
           chunk_col;
  }



  inline
  std::size_t
  Accessor::reduced_index() const
  {
    Assert (is_valid_entry() == true, ExcInvalidIterator());

    return reduced_accessor.index_within_sparsity;
  }




  inline
  bool
  Accessor::operator == (const Accessor &other) const
  {
    // no need to check for equality of sparsity patterns as this is done in
    // the reduced case already and every ChunkSparsityPattern has its own
    // reduced sparsity pattern
    return (reduced_accessor == other.reduced_accessor &&
            chunk_row == other.chunk_row &&
            chunk_col == other.chunk_col);
  }



  inline
  bool
  Accessor::operator < (const Accessor &other) const
  {
    Assert (sparsity_pattern == other.sparsity_pattern,
            ExcInternalError());

    if (chunk_row != other.chunk_row)
      {
        if (reduced_accessor.index_within_sparsity ==
            reduced_accessor.sparsity_pattern->n_nonzero_elements())
          return false;
        if (other.reduced_accessor.index_within_sparsity ==
            reduced_accessor.sparsity_pattern->n_nonzero_elements())
          return true;

        const unsigned int
        global_row = sparsity_pattern->get_chunk_size()*reduced_accessor.row()
                     +chunk_row,
                     other_global_row = sparsity_pattern->get_chunk_size()*
                                        other.reduced_accessor.row()+other.chunk_row;
        if (global_row < other_global_row)
          return true;
        else if (global_row > other_global_row)
          return false;
      }

    return (reduced_accessor.index_within_sparsity <
            other.reduced_accessor.index_within_sparsity ||
            (reduced_accessor.index_within_sparsity ==
             other.reduced_accessor.index_within_sparsity &&
             chunk_col < other.chunk_col));
  }


  inline
  void
  Accessor::advance ()
  {
    const unsigned int chunk_size = sparsity_pattern->get_chunk_size();
    Assert (chunk_row < chunk_size && chunk_col < chunk_size,
            ExcIteratorPastEnd());
    Assert (reduced_accessor.row() * chunk_size + chunk_row <
            sparsity_pattern->n_rows()
            &&
            reduced_accessor.column() * chunk_size + chunk_col <
            sparsity_pattern->n_cols(),
            ExcIteratorPastEnd());
    if (chunk_size == 1)
      {
        reduced_accessor.advance();
        return;
      }

    ++chunk_col;

    // end of chunk
    if (chunk_col == chunk_size
        ||
        reduced_accessor.column() * chunk_size + chunk_col ==
        sparsity_pattern->n_cols())
      {
        const unsigned int reduced_row = reduced_accessor.row();
        // end of row
        if (reduced_accessor.index_within_sparsity + 1 ==
            reduced_accessor.sparsity_pattern->rowstart[reduced_row+1])
          {
            ++chunk_row;

            chunk_col = 0;

            // end of chunk rows or end of matrix
            if (chunk_row == chunk_size ||
                (reduced_row * chunk_size + chunk_row ==
                 sparsity_pattern->n_rows()))
              {
                chunk_row = 0;
                reduced_accessor.advance();
              }
            // go back to the beginning of the same reduced row but with
            // chunk_row increased by one
            else
              reduced_accessor.index_within_sparsity =
                reduced_accessor.sparsity_pattern->rowstart[reduced_row];
          }
        // advance within chunk
        else
          {
            reduced_accessor.advance();
            chunk_col = 0;
          }
      }
  }



  inline
  Iterator::Iterator (const ChunkSparsityPattern *sparsity_pattern,
                      const unsigned int          row)
    :
    accessor(sparsity_pattern, row)
  {}



  inline
  Iterator &
  Iterator::operator++ ()
  {
    accessor.advance ();
    return *this;
  }



  inline
  Iterator
  Iterator::operator++ (int)
  {
    const Iterator iter = *this;
    accessor.advance ();
    return iter;
  }



  inline
  const Accessor &
  Iterator::operator* () const
  {
    return accessor;
  }



  inline
  const Accessor *
  Iterator::operator-> () const
  {
    return &accessor;
  }


  inline
  bool
  Iterator::operator == (const Iterator &other) const
  {
    return (accessor == other.accessor);
  }



  inline
  bool
  Iterator::operator != (const Iterator &other) const
  {
    return ! (accessor == other.accessor);
  }


  inline
  bool
  Iterator::operator < (const Iterator &other) const
  {
    return accessor < other.accessor;
  }

}



inline
ChunkSparsityPattern::iterator
ChunkSparsityPattern::begin () const
{
  return iterator(this, 0);
}


inline
ChunkSparsityPattern::iterator
ChunkSparsityPattern::end () const
{
  return iterator(this, n_rows());
}



inline
ChunkSparsityPattern::iterator
ChunkSparsityPattern::begin (const unsigned int r) const
{
  Assert (r<n_rows(), ExcIndexRange(r,0,n_rows()));
  return iterator(this, r);
}



inline
ChunkSparsityPattern::iterator
ChunkSparsityPattern::end (const unsigned int r) const
{
  Assert (r<n_rows(), ExcIndexRange(r,0,n_rows()))
  return iterator(this, r+1);
}



inline
ChunkSparsityPattern::size_type
ChunkSparsityPattern::n_rows () const
{
  return rows;
}


inline
ChunkSparsityPattern::size_type
ChunkSparsityPattern::n_cols () const
{
  return cols;
}



inline
ChunkSparsityPattern::size_type
ChunkSparsityPattern::get_chunk_size () const
{
  return chunk_size;
}



inline
bool
ChunkSparsityPattern::is_compressed () const
{
  return sparsity_pattern.compressed;
}



inline
bool
ChunkSparsityPattern::optimize_diagonal () const
{
  return sparsity_pattern.store_diagonal_first_in_row;
}



template <typename ForwardIterator>
inline
void
ChunkSparsityPattern::copy_from (const size_type n_rows,
                                 const size_type n_cols,
                                 const ForwardIterator begin,
                                 const ForwardIterator end,
                                 const size_type chunk_size,
                                 const bool)
{
  copy_from (n_rows, n_cols, begin, end, chunk_size);
}



template <typename ForwardIterator>
void
ChunkSparsityPattern::copy_from (const size_type       n_rows,
                                 const size_type       n_cols,
                                 const ForwardIterator begin,
                                 const ForwardIterator end,
                                 const size_type       chunk_size)
{
  Assert (static_cast<size_type>(std::distance (begin, end)) == n_rows,
          ExcIteratorRange (std::distance (begin, end), n_rows));

  // first determine row lengths for each row. if the matrix is quadratic,
  // then we might have to add an additional entry for the diagonal, if that
  // is not yet present. as we have to call compress anyway later on, don't
  // bother to check whether that diagonal entry is in a certain row or not
  const bool is_square = (n_rows == n_cols);
  std::vector<size_type> row_lengths;
  row_lengths.reserve(n_rows);
  for (ForwardIterator i=begin; i!=end; ++i)
    row_lengths.push_back (std::distance (i->begin(), i->end())
                           +
                           (is_square ? 1 : 0));
  reinit (n_rows, n_cols, row_lengths, chunk_size);

  // now enter all the elements into the matrix
  size_type row = 0;
  typedef typename std::iterator_traits<ForwardIterator>::value_type::const_iterator inner_iterator;
  for (ForwardIterator i=begin; i!=end; ++i, ++row)
    {
      const inner_iterator end_of_row = i->end();
      for (inner_iterator j=i->begin(); j!=end_of_row; ++j)
        {
          const size_type col
            = internal::SparsityPatternTools::get_column_index_from_iterator(*j);
          Assert (col < n_cols, ExcInvalidIndex(col,n_cols));

          add (row, col);
        }
    }

  // finally compress everything. this also sorts the entries within each row
  compress ();
}


#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
