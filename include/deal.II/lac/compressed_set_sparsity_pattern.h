// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2013 by the deal.II authors
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

#ifndef __deal2__compressed_set_sparsity_pattern_h
#define __deal2__compressed_set_sparsity_pattern_h


#include <deal.II/base/config.h>
#include <deal.II/base/subscriptor.h>
#include <deal.II/lac/exceptions.h>

#include <vector>
#include <algorithm>
#include <set>

DEAL_II_NAMESPACE_OPEN

template <typename number> class SparseMatrix;

/*! @addtogroup Sparsity
 *@{
 */


/**
 * This class acts as an intermediate form of the
 * SparsityPattern class. From the interface it mostly
 * represents a SparsityPattern object that is kept compressed
 * at all times. However, since the final sparsity pattern is not
 * known while constructing it, keeping the pattern compressed at all
 * times can only be achieved at the expense of either increased
 * memory or run time consumption upon use. The main purpose of this
 * class is to avoid some memory bottlenecks, so we chose to implement
 * it memory conservative, but the chosen data format is too unsuited
 * to be used for actual matrices. It is therefore necessary to first
 * copy the data of this object over to an object of type
 * SparsityPattern before using it in actual matrices.
 *
 * Another viewpoint is that this class does not need up front allocation of a
 * certain amount of memory, but grows as necessary.  An extensive description
 * of sparsity patterns can be found in the documentation of the @ref Sparsity
 * module.
 *
 * This class is an example of the "dynamic" type of @ref Sparsity. It
 * is discussed in the step-27 and @ref step_22
 * "step-22" tutorial programs.
 *
 * <h3>Interface</h3>
 *
 * Since this class is intended as an intermediate replacement of the
 * SparsityPattern class, it has mostly the same interface, with
 * small changes where necessary. In particular, the add()
 * function, and the functions inquiring properties of the sparsity
 * pattern are the same.
 *
 *
 * <h3>Usage</h3>
 *
 * Use this class as follows:
 * @code
 * CompressedSetSparsityPattern compressed_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  compressed_pattern);
 * constraints.condense (compressed_pattern);
 *
 * SparsityPattern sp;
 * sp.copy_from (compressed_pattern);
 * @endcode
 *
 * See also step-11 and step-18 for usage
 * patterns of the related CompressedSparsityPattern class, and
 * step-27 of the current class.
 *
 * <h3>Notes</h3>
 *
 * There are several, exchangeable variations of this class, see @ref Sparsity,
 * section '"Dynamic" or "compressed" sparsity patterns' for more information.
 *
 * This class is a variation of the CompressedSparsityPattern class.
 * Instead of using sorted vectors together with a caching algorithm
 * for storing the column indices of nonzero entries, the std::set
 * container is used. This solution might not be the fastest in all
 * situations, but seems to work much better than the
 * CompressedSparsityPattern in the context of hp-adaptivity (see for
 * example step-27), or generally when there are many
 * nonzero entries in each row of a matrix (see @ref step_22
 * "step-22").  On the other hand, a benchmark where nonzero entries
 * were randomly inserted into the sparsity pattern revealed that this
 * class is slower by a factor 4-6 in this situation. Hence, currently
 * the suggestion is to carefully analyze which of the
 * CompressedSparsityPattern classes works best in a certain
 * setting. An algorithm which performs equally well in all situations
 * still has to be found.
 *
 *
 * @author Oliver Kayser-Herold, 2007
 */
class CompressedSetSparsityPattern : public Subscriptor
{
public:
  /**
   * Declare the type for container size.
   */
  typedef types::global_dof_index size_type;

  /**
   * An iterator that can be used to
   * iterate over the elements of a single
   * row. The result of dereferencing such
   * an iterator is a column index.
   */
  typedef std::set<size_type>::const_iterator row_iterator;


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
  CompressedSetSparsityPattern ();

  /**
   * Copy constructor. This constructor is
   * only allowed to be called if the
   * matrix structure to be copied is
   * empty. This is so in order to prevent
   * involuntary copies of objects for
   * temporaries, which can use large
   * amounts of computing time.  However,
   * copy constructors are needed if yo
   * want to use the STL data types on
   * classes like this, e.g. to write such
   * statements like <tt>v.push_back
   * (CompressedSetSparsityPattern());</tt>,
   * with @p v a vector of @p
   * CompressedSetSparsityPattern objects.
   */
  CompressedSetSparsityPattern (const CompressedSetSparsityPattern &);

  /**
   * Initialize a rectangular
   * matrix with @p m rows and
   * @p n columns.
   */
  CompressedSetSparsityPattern (const size_type m,
                                const size_type n);

  /**
   * Initialize a square matrix of
   * dimension @p n.
   */
  CompressedSetSparsityPattern (const size_type n);

  /**
   * Copy operator. For this the
   * same holds as for the copy
   * constructor: it is declared,
   * defined and fine to be called,
   * but the latter only for empty
   * objects.
   */
  CompressedSetSparsityPattern &operator = (const CompressedSetSparsityPattern &);

  /**
   * Reallocate memory and set up
   * data structures for a new
   * matrix with @p m rows and
   * @p n columns, with at most
   * max_entries_per_row() nonzero
   * entries per row.
   */
  void reinit (const size_type m,
               const size_type n);

  /**
   * Since this object is kept
   * compressed at all times anway,
   * this function does nothing,
   * but is declared to make the
   * interface of this class as
   * much alike as that of the
   * SparsityPattern class.
   */
  void compress ();

  /**
   * Return whether the object is
   * empty. It is empty if no
   * memory is allocated, which is
   * the same as that both
   * dimensions are zero.
   */
  bool empty () const;

  /**
   * Return the maximum number of
   * entries per row. Note that
   * this number may change as
   * entries are added.
   */
  size_type max_entries_per_row () const;

  /**
   * Add a nonzero entry to the
   * matrix. If the entry already
   * exists, nothing bad happens.
   */
  void add (const size_type i,
            const size_type j);

  /**
   * Add several nonzero entries to the
   * specified row of the matrix. If the
   * entries already exist, nothing bad
   * happens.
   */
  template <typename ForwardIterator>
  void add_entries (const size_type row,
                    ForwardIterator begin,
                    ForwardIterator end,
                    const bool      indices_are_sorted = false);

  /**
   * Check if a value at a certain
   * position may be non-zero.
   */
  bool exists (const size_type i,
               const size_type j) const;

  /**
   * Make the sparsity pattern
   * symmetric by adding the
   * sparsity pattern of the
   * transpose object.
   *
   * This function throws an
   * exception if the sparsity
   * pattern does not represent a
   * square matrix.
   */
  void symmetrize ();

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
   * Print the sparsity of the matrix in a
   * format that @p gnuplot understands and
   * which can be used to plot the sparsity
   * pattern in a graphical way. The format
   * consists of pairs <tt>i j</tt> of
   * nonzero elements, each representing
   * one entry of this matrix, one per line
   * of the output file. Indices are
   * counted from zero on, as usual. Since
   * sparsity patterns are printed in the
   * same way as matrices are displayed, we
   * print the negative of the column
   * index, which means that the
   * <tt>(0,0)</tt> element is in the top
   * left rather than in the bottom left
   * corner.
   *
   * Print the sparsity pattern in
   * gnuplot by setting the data style
   * to dots or points and use the
   * @p plot command.
   */
  void print_gnuplot (std::ostream &out) const;

  /**
   * Return number of rows of this
   * matrix, which equals the dimension
   * of the image space.
   */
  size_type n_rows () const;

  /**
   * Return number of columns of this
   * matrix, which equals the dimension
   * of the range space.
   */
  size_type n_cols () const;

  /**
   * Number of entries in a specific row.
   */
  size_type row_length (const size_type row) const;

  /**
   * Return an iterator that can loop over
   * all entries in the given
   * row. Dereferencing the iterator yields
   * a column index.
   */
  row_iterator row_begin (const size_type row) const;

  /**
   * End iterator for the given row.
   */
  row_iterator row_end (const size_type row) const;


  /**
   * Compute the bandwidth of the matrix
   * represented by this structure. The
   * bandwidth is the maximum of
   * $|i-j|$ for which the index pair
   * $(i,j)$ represents a nonzero entry
   * of the matrix.
   */
  size_type bandwidth () const;

  /**
   * Return the number of nonzero elements
   * allocated through this sparsity
   * pattern.
   */
  size_type n_nonzero_elements () const;

  /**
   * Return whether this object stores only
   * those entries that have been added
   * explicitly, or if the sparsity pattern
   * contains elements that have been added
   * through other means (implicitly) while
   * building it. For the current class,
   * the result is always true.
   *
   * This function mainly serves the
   * purpose of describing the current
   * class in cases where several kinds of
   * sparsity patterns can be passed as
   * template arguments.
   */
  static
  bool stores_only_added_elements ();

private:
  /**
   * Number of rows that this sparsity
   * structure shall represent.
   */
  size_type rows;

  /**
   * Number of columns that this sparsity
   * structure shall represent.
   */
  size_type cols;

  /**
   * For each row of the matrix, store the
   * allocated non-zero entries as a
   * std::set of column indices. For a
   * discussion of storage schemes see the
   * CompressedSparsityPattern::Line class.
   */
  struct Line
  {
    std::set<size_type> entries;

    /**
     * Constructor.
     */
    Line ();

    /**
     * Add the given column number to
     * this line.
     */
    void add (const size_type col_num);

    /**
     * Add the columns specified by the
     * iterator range to this line.
     */
    template <typename ForwardIterator>
    void add_entries (ForwardIterator begin,
                      ForwardIterator end);
  };


  /**
   * Actual data: store for each
   * row the set of nonzero
   * entries.
   */
  std::vector<Line> lines;
};

/*@}*/
/*---------------------- Inline functions -----------------------------------*/


inline
CompressedSetSparsityPattern::Line::Line ()
{}



inline
void
CompressedSetSparsityPattern::Line::add (const size_type j)
{
  entries.insert (j);
}



template <typename ForwardIterator>
inline
void
CompressedSetSparsityPattern::Line::add_entries (ForwardIterator begin,
                                                 ForwardIterator end)
{
  entries.insert (begin, end);
}



inline
CompressedSetSparsityPattern::size_type
CompressedSetSparsityPattern::n_rows () const
{
  return rows;
}



inline
CompressedSetSparsityPattern::size_type
CompressedSetSparsityPattern::n_cols () const
{
  return cols;
}



inline
void
CompressedSetSparsityPattern::add (const size_type i,
                                   const size_type j)
{
  Assert (i<rows, ExcIndexRangeType<size_type>(i, 0, rows));
  Assert (j<cols, ExcIndexRangeType<size_type>(j, 0, cols));

  lines[i].add (j);
}



template <typename ForwardIterator>
inline
void
CompressedSetSparsityPattern::add_entries (const size_type row,
                                           ForwardIterator begin,
                                           ForwardIterator end,
                                           const bool         /*indices_are_sorted*/)
{
  Assert (row < rows, ExcIndexRangeType<size_type> (row, 0, rows));

  lines[row].add_entries (begin, end);
}



inline
CompressedSetSparsityPattern::size_type
CompressedSetSparsityPattern::row_length (const size_type row) const
{
  Assert (row < n_rows(), ExcIndexRangeType<size_type> (row, 0, n_rows()));

  return lines[row].entries.size();
}



inline
CompressedSetSparsityPattern::row_iterator
CompressedSetSparsityPattern::row_begin (const size_type row) const
{
  return (lines[row].entries.begin ());
}



inline
CompressedSetSparsityPattern::row_iterator
CompressedSetSparsityPattern::row_end (const size_type row) const
{
  return (lines[row].entries.end ());
}



inline
bool
CompressedSetSparsityPattern::stores_only_added_elements ()
{
  return true;
}


DEAL_II_NAMESPACE_CLOSE

#endif
