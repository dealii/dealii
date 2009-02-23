//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__compressed_simple_sparsity_pattern_h
#define __deal2__compressed_simple_sparsity_pattern_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>

#include <vector>
#include <algorithm>

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
 * This class is an example of the "dynamic" type of @ref Sparsity.
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
 * @verbatim
 * CompressedSimpleSparsityPattern compressed_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  compressed_pattern);
 * constraints.condense (compressed_pattern);
 *
 * SparsityPattern sp;
 * sp.copy_from (compressed_pattern);
 * @endverbatim
 *
 *
 * <h3>Notes</h3>
 *
 * There are several, exchangeable variations of this class, see @ref Sparsity,
 * section '"Dynamic" or "compressed" sparsity patterns' for more information.
 *
 * @author Timo Heister, 2008
 */
class CompressedSimpleSparsityPattern : public Subscriptor
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
    CompressedSimpleSparsityPattern ();

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
				      * (CompressedSparsityPattern());</tt>,
				      * with @p v a vector of @p
				      * CompressedSparsityPattern objects.
				      */
    CompressedSimpleSparsityPattern (const CompressedSimpleSparsityPattern &);

				     /**
				      * Initialize a rectangular
				      * matrix with @p m rows and
				      * @p n columns.
				      */
    CompressedSimpleSparsityPattern (const unsigned int m,
			       const unsigned int n);

				     /**
				      * Initialize a square matrix of
				      * dimension @p n.
				      */
    CompressedSimpleSparsityPattern (const unsigned int n);

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    CompressedSimpleSparsityPattern & operator = (const CompressedSimpleSparsityPattern &);

				     /**
				      * Reallocate memory and set up
				      * data structures for a new
				      * matrix with @p m rows and
				      * @p n columns, with at most
				      * max_entries_per_row() nonzero
				      * entries per row.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n);

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
    unsigned int max_entries_per_row () const;

				     /**
				      * Add a nonzero entry to the
				      * matrix. If the entry already
				      * exists, nothing bad happens.
				      */
    void add (const unsigned int i,
	      const unsigned int j);

				     /**
				      * Check if a value at a certain
				      * position may be non-zero.
				      */
    bool exists (const unsigned int i,
                 const unsigned int j) const;

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
    unsigned int n_rows () const;

				     /**
				      * Return number of columns of this
				      * matrix, which equals the dimension
				      * of the range space.
				      */
    unsigned int n_cols () const;

				     /**
				      * Number of entries in a specific row.
				      */
    unsigned int row_length (const unsigned int row) const;

				     /**
				      * Access to column number field.
				      * Return the column number of
				      * the @p indexth entry in @p row.
				      */
    unsigned int column_number (const unsigned int row,
				const unsigned int index) const;

				     /**
				      * Compute the bandwidth of the matrix
				      * represented by this structure. The
				      * bandwidth is the maximum of
				      * $|i-j|$ for which the index pair
				      * $(i,j)$ represents a nonzero entry
				      * of the matrix.
				      */
    unsigned int bandwidth () const;

				     /**
				      * Return the number of nonzero elements
				      * allocated through this sparsity pattern.
				      */
    unsigned int n_nonzero_elements () const;

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
    unsigned int rows;

				     /**
				      * Number of columns that this sparsity
				      * structure shall represent.
				      */
    unsigned int cols;

                                     /**
                                      * Store some data for each row
                                      * describing which entries of this row
                                      * are nonzero. Data is stored sorted in
                                      * the @p entries std::vector.
                                      * The vector per row is dynamically
                                      * growing upon insertion doubling its
                                      * memory each time.
                                      */
    struct Line
    {
      public:
                                         /**
                                          * Storage for the column indices of
                                          * this row. This array is always
                                          * kept sorted.
                                          */
        std::vector<unsigned int> entries;

                                         /**
                                          * Constructor.
                                          */
        Line ();

                                         /**
                                          * Add the given column number to
                                          * this line.
                                          */
        void add (const unsigned int col_num);
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
void
CompressedSimpleSparsityPattern::Line::add (const unsigned int j)
{
				   // first check the last element (or if line
				   // is still empty)
  if ( (entries.size()==0) || ( entries.back() < j) )
    {
      entries.push_back(j);
      return;
    }

				   // do a binary search to find the place
				   // where to insert:
  std::vector<unsigned int>::iterator it = std::lower_bound(entries.begin(),
							    entries.end(),
							    j);

				   // If this entry is a duplicate, exit
				   // immediately
  if (*it == j)
    return;

				   // Insert at the right place in the
				   // vector. Vector grows automatically to
				   // fit elements. Always doubles its size.
  entries.insert(it, j);
}



inline
unsigned int
CompressedSimpleSparsityPattern::n_rows () const
{
  return rows;
}



inline
unsigned int
CompressedSimpleSparsityPattern::n_cols () const
{
  return cols;
}



inline
void
CompressedSimpleSparsityPattern::add (const unsigned int i,
				      const unsigned int j)
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));

  lines[i].add (j);
}



inline
CompressedSimpleSparsityPattern::Line::Line ()
{}



inline
unsigned int
CompressedSimpleSparsityPattern::row_length (const unsigned int row) const
{
  Assert (row < n_rows(), ExcIndexRange (row, 0, n_rows()));

  return lines[row].entries.size();
}



inline
unsigned int
CompressedSimpleSparsityPattern::column_number (const unsigned int row,
						const unsigned int index) const
{
  Assert (row < n_rows(), ExcIndexRange (row, 0, n_rows()));
  Assert (index < lines[row].entries.size(),
	  ExcIndexRange (index, 0, lines[row].entries.size()));

  return lines[row].entries[index];
}


inline
bool
CompressedSimpleSparsityPattern::stores_only_added_elements ()
{
  return true;
}


DEAL_II_NAMESPACE_CLOSE

#endif
