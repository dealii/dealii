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
#ifndef __deal2__compressed_sparsity_pattern_h
#define __deal2__compressed_sparsity_pattern_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>

template <typename number> class SparseMatrix;

#include <vector>
#include <algorithm>

/*! @addtogroup Matrix1
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
 * Another viewpoint is that this class does not need up front
 * allocation of a certain amount of memory, but grows as necessary.
 *
 *
 * <h3>Rationale</h3>
 *
 * When constructing the sparsity pattern of a matrix, you usually
 * first have to provide an empty sparsity pattern object with a fixed
 * maximal number of entries per row. To find out about this maximal
 * row length, one usually calls the function
 * DoFHandler@p ::max_couplings_per_dof which returns an
 * estimate for that quantity. While this estimate is usually quite
 * good in 2d and exact in 1d, it is often significantly too large in
 * 3d and especially for higher order elements. Furthermore, normally
 * only a small fraction of the rows of a matrix will end up having
 * the maximal number of nonzero entries per row (usually those nodes
 * adjacent to hanging nodes), most have much less. In effect, the
 * empty SparsityPattern object has allocated much too much
 * memory. Although this unnecessarily allocated memory is later freed
 * when SparsityPattern@p ::compress is called, this
 * overallocation has, with higher order elements and in 3d, sometimes
 * been so large that the program aborted due to lack of memory.
 *
 * This class therefore provides an alternative representation of a
 * sparsity pattern: we don't specify a maximal row length initially,
 * but store a set of column indices indicating possible nonzero
 * entries in the sparsity pattern for each row. This is very much
 * like the final "compressed" format used in the
 * SparsityPattern object after compression, but uses a less
 * compact memory storage format, since the exact number of entries
 * per row is only known a posteriori and since it may change (for the
 * SparsityPattern class, no more changes are allowed after
 * compressing it). We can therefore not store all the column indices
 * in a big array, but have to use a vector of sets. This can later be
 * used to actually initialize a SparsityPattern object with the
 * then final set of necessary indices.
 *
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
 * CompressedSparsityPattern compressed_pattern (dof_handler.n_dofs());
 * DoFTools::make_sparsity_pattern (dof_handler,
 *                                  compressed_pattern);
 * constraints.condense (compressed_pattern);
 *
 * SparsityPattern sp;
 * sp.copy_from (compressed_pattern);
 * @endverbatim
 *
 *
 * @author Wolfgang Bangerth, 2001
 */
class CompressedSparsityPattern : public Subscriptor
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
				      * the @p reinit function.
				      */
    CompressedSparsityPattern ();
    
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
				      *
				      * Usually, it is sufficient to
				      * use the explicit keyword to
				      * disallow unwanted temporaries,
				      * but for the STL vectors, this
				      * does not work. Since copying a
				      * structure like this is not
				      * useful anyway because multiple
				      * matrices can use the same
				      * sparsity structure, copies are
				      * only allowed for empty
				      * objects, as described above.
				      */
    CompressedSparsityPattern (const CompressedSparsityPattern &);

				     /**
				      * Initialize a rectangular
				      * matrix with @p m rows and
				      * @p n columns.
				      */
    CompressedSparsityPattern (const unsigned int m,
			       const unsigned int n);
    
				     /**
				      * Initialize a square matrix of
				      * dimension @p n.
				      */
    CompressedSparsityPattern (const unsigned int n);

				     /**
				      * Copy operator. For this the
				      * same holds as for the copy
				      * constructor: it is declared,
				      * defined and fine to be called,
				      * but the latter only for empty
				      * objects.
				      */
    CompressedSparsityPattern & operator = (const CompressedSparsityPattern &);
    
				     /**
				      * Reallocate memory and set up
				      * data structures for a new
				      * matrix with @p m rows and
				      * @p n columns, with at most
				      * @p max_per_row nonzero
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
				      * of this matrix. Actually, it returns
				      * the number of entries in the sparsity
				      * pattern; if any of the entries should
				      * happen to be zero, it is counted
				      * anyway.
				      *
				      * This function may only be called if
				      * the matrix struct is compressed. It
				      * does not make too much sense otherwise
				      * anyway.
				      */
    unsigned int n_nonzero_elements () const;

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
                                      * are nonzero. Data is organized as
                                      * follows: if an entry is added to a
                                      * row, it is first added to the @p cache
                                      * variable, irrespective of whether an
                                      * entry with same column number has
                                      * already been added. Only if the cache
                                      * is full do we flush it by removing
                                      * duplicates, removing entries that are
                                      * already stored in the @p entries
                                      * array, sorting everything, and merging
                                      * the two arrays.
                                      *
                                      * The reasoning behind this scheme is
                                      * that memory allocation is expensive,
                                      * and we only want to do it when really
                                      * necessary. Previously (in deal.II
                                      * versions up to 5.0), we used to store
                                      * the column indices inside a std::set,
                                      * but this would allocate 20 bytes each
                                      * time we added an entry. Using the
                                      * present scheme, we only need to
                                      * allocate memory once for every 8 added
                                      * entries, and we waste a lot less
                                      * memory by not using a balanced tree
                                      * for storing column indices.
                                      *
                                      * Since some functions that are @p const
                                      * need to access the data of this
                                      * object, but need to flush caches
                                      * before, the @p flush_cache function is
                                      * marked const, and the data members are
                                      * marked @p mutable.
                                      *
                                      * A small testseries about the size of
                                      * the cache showed that the run time of
                                      * a small program just testing the
                                      * compressed sparsity pattern element
                                      * insertion routine ran for 3.6 seconds
                                      * with a cache size of 8, and 4.2
                                      * seconds with a cache size of 16. We
                                      * deem even smaller cache sizes
                                      * undesirable, since they lead to more
                                      * memory allocations, while larger cache
                                      * sizes lead to waste of memory. The
                                      * original version of this class, with
                                      * one std::set per row took 8.2 seconds
                                      * on the same program.
                                      */
    struct Line
    {
                                         /**
                                          * Storage for the column indices of
                                          * this row, unless they are still in
                                          * the cache. This array is always
                                          * kept sorted.
                                          */
        mutable std::vector<unsigned int> entries;

                                         /**
                                          * Size of the cache.
                                          */
        static const unsigned int cache_size = 8;
        
                                         /**
                                          * Cache of entries that have not yet
                                          * been written to @p entries;
                                          */
        mutable unsigned int cache[cache_size];

                                         /**
                                          * Number of entries in the cache.
                                          */
        mutable unsigned int cache_entries;

                                         /**
                                          * Constructor.
                                          */
        Line ();

                                         /**
                                          * Add the given column number to
                                          * this line.
                                          */
        void add (const unsigned int col_num);

                                         /**
                                          * Flush the cache my merging it with
                                          * the @p entries array.
                                          */
        void flush_cache () const;
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
CompressedSparsityPattern::Line::flush_cache () const
{
                                   // do nothing if the cache is empty
  if (cache_entries == 0)
    return;
  
                                   // first sort the entries in the cache, so
                                   // that it is easier to merge it with the
                                   // main array. note that due to the way
                                   // add() inserts elements, there can be no
                                   // duplicates in the cache
                                   //
                                   // do the sorting in a way that is fast for
                                   // the small cache sizes we have
                                   // here. basically, use bubble sort
  switch (cache_entries)
    {
      case 1:
      {
        break;
      }

      case 2:
      {
        if (cache[1] < cache[0])
          std::swap (cache[0], cache[1]);
        break;
      }

      case 3:
      {
        if (cache[1] < cache[0])
          std::swap (cache[0], cache[1]);
        if (cache[2] < cache[1])
          std::swap (cache[1], cache[2]);
        if (cache[1] < cache[0])
          std::swap (cache[0], cache[1]);
        break;
      }

      case 4:
      case 5:
      case 6:
      case 7:
      {
        for (unsigned int i=0; i<cache_entries; ++i)
          for (unsigned int j=i+1; j<cache_entries; ++j)
            if (cache[j] < cache[i])
              std::swap (cache[i], cache[j]);
        break;
      }

      default:
      {
        std::sort (&cache[0], &cache[cache_entries]);
        break;
      }
    }

                                   // next job is to merge the two
                                   // arrays. special case the case that the
                                   // original array is empty.
  if (entries.size() == 0)
    {
      entries.resize (cache_entries);
      for (unsigned int i=0; i<cache_entries; ++i)
        entries[i] = cache[i];
    }
  else
    {
                                       // first count how many of the cache
                                       // entries are already in the main
                                       // array, so that we can efficiently
                                       // allocate memory
      unsigned int n_new_entries = 0;
      {
        unsigned int cache_position = 0;
        unsigned int entry_position = 0;
        while ((entry_position<entries.size()) &&
               (cache_position<cache_entries))
          {
            ++n_new_entries;
            if (entries[entry_position] < cache[cache_position])
              ++entry_position;
            else if (entries[entry_position] == cache[cache_position])
              {
                ++entry_position;
                ++cache_position;
              }
            else
              ++cache_position;
          }

                                         // scoop up leftovers in arrays
        n_new_entries += (entries.size() - entry_position) +
                         (cache_entries - cache_position);
      }

                                       // then allocate new memory and merge
                                       // arrays
      std::vector<unsigned int> new_entries;
      new_entries.reserve (n_new_entries);
      unsigned int cache_position = 0;
      unsigned int entry_position = 0;
      while ((entry_position<entries.size()) && (cache_position<cache_entries))
        if (entries[entry_position] < cache[cache_position])
          {
            new_entries.push_back (entries[entry_position]);
            ++entry_position;
          }
        else if (entries[entry_position] == cache[cache_position])
          {
            new_entries.push_back (entries[entry_position]);
            ++entry_position;
            ++cache_position;
          }
        else
          {
            new_entries.push_back (cache[cache_position]);
            ++cache_position;
          }

                                       // copy remaining elements from the
                                       // array that we haven't finished. note
                                       // that at most one of the following
                                       // loops will run at all
      for (; entry_position < entries.size(); ++entry_position)
        new_entries.push_back (entries[entry_position]);
      for (; cache_position < cache_entries; ++cache_position)
        new_entries.push_back (cache[cache_position]);
      
      Assert (new_entries.size() == n_new_entries,
              ExcInternalError());

                                       // finally swap old and new array, and
                                       // set cache size to zero
      new_entries.swap (entries);
    }
  
  cache_entries = 0;
}



inline
void
CompressedSparsityPattern::Line::add (const unsigned int j)
{
                                   // first check whether this entry is
                                   // already in the cache. if so, we can
                                   // safely return
  for (unsigned int i=0; i<cache_entries; ++i)
    if (cache[i] == j)
      return;

                                   // if not, see whether there is still some
                                   // space in the cache. if not, then flush
                                   // the cache first
  if (cache_entries == cache_size)
    flush_cache ();
  
  cache[cache_entries] = j;
  ++cache_entries;
}



inline
unsigned int
CompressedSparsityPattern::n_rows () const
{
  return rows;
}



inline
unsigned int
CompressedSparsityPattern::n_cols () const
{
  return cols;
}



inline
void
CompressedSparsityPattern::add (const unsigned int i,
				const unsigned int j)
{
  Assert (i<rows, ExcIndexRange(i, 0, rows));
  Assert (j<cols, ExcIndexRange(j, 0, cols));

  lines[i].add (j);
}



inline
CompressedSparsityPattern::Line::Line ()
                :
                cache_entries (0)
{}



inline
unsigned int
CompressedSparsityPattern::row_length (const unsigned int row) const
{
  Assert (row < n_rows(), ExcIndexRange (row, 0, n_rows()));
  
  lines[row].flush_cache ();
  return lines[row].entries.size();
}



inline
unsigned int
CompressedSparsityPattern::column_number (const unsigned int row,
					  const unsigned int index) const
{
  Assert (row < n_rows(), ExcIndexRange (row, 0, n_rows()));
  Assert (index < lines[row].entries.size(),
	  ExcIndexRange (index, 0, lines[row].entries.size()));

  lines[row].flush_cache ();
  return lines[row].entries[index];
}




#endif
