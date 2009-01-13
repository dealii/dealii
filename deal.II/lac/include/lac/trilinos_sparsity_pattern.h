//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__trilinos_sparsity_pattern_h
#define __deal2__trilinos_sparsity_pattern_h


#include <base/config.h>
#include <base/subscriptor.h>
#include <lac/exceptions.h>

#include <base/std_cxx0x/shared_ptr.h>
#include <vector>
#include <cmath>
#include <memory>

#ifdef DEAL_II_USE_TRILINOS

#  include <Epetra_FECrsGraph.h>
#  include <Epetra_Map.h>
#  ifdef DEAL_II_COMPILER_SUPPORTS_MPI
#    include <Epetra_MpiComm.h>
#    include "mpi.h"
#  else
#    include "Epetra_SerialComm.h"
#  endif

DEAL_II_NAMESPACE_OPEN

                                   // forward declarations
class SparsityPattern;
class CompressedSparsityPattern;
class CompressedSetSparsityPattern;
class CompressedSimpleSparsityPattern;

namespace TrilinosWrappers
{
                                   // forward declarations
  class SparsityPattern;

  namespace SparsityPatternIterators
  {
/**
 * STL conforming iterator. This class acts as an iterator walking
 * over the elements of Trilinos sparsity pattern.
 *
 * @ingroup TrilinosWrappers
 * @author Martin Kronbichler, Wolfgang Bangerth, 2008
 */
    class const_iterator
    {
      private:
                                       /**
					* Accessor class for iterators
					*/
        class Accessor
        {
          public:
                                       /**
					* Constructor. Since we use
					* accessors only for read
					* access, a const matrix
					* pointer is sufficient.
					*/
            Accessor (const SparsityPattern *sparsity_pattern,
                      const unsigned int     row,
                      const unsigned int     index);

                                       /**
					* Row number of the element
					* represented by this object.
					*/
            unsigned int row() const;

                                       /**
					* Index in row of the element
					* represented by this object.
					*/
            unsigned int index() const;

                                       /**
					* Column number of the element
					* represented by this object.
					*/
            unsigned int column() const;

                                       /**
					* Exception
					*/
            DeclException0 (ExcBeyondEndOfSparsityPattern);

	                               /**
					* Exception
					*/
            DeclException3 (ExcAccessToNonlocalRow,
                            int, int, int,
                            << "You tried to access row " << arg1
                            << " of a distributed sparsity pattern, "
			    << " but only rows " << arg2 << " through " << arg3
                            << " are stored locally and can be accessed.");
            
          private:
                                       /**
					* The matrix accessed.
					*/
            mutable SparsityPattern *sparsity_pattern;

                                       /**
					* Current row number.
					*/
            unsigned int a_row;

                                       /**
					* Current index in row.
					*/
            unsigned int a_index;

                                       /**
					* Cache where we store the
					* column indices of the
					* present row. This is
					* necessary, since Trilinos
					* makes access to the elements
					* of its matrices rather hard,
					* and it is much more
					* efficient to copy all column
					* entries of a row once when
					* we enter it than repeatedly
					* asking Trilinos for
					* individual ones. This also
					* makes some sense since it is
					* likely that we will access
					* them sequentially anyway.
					*
					* In order to make copying of
					* iterators/accessor of
					* acceptable performance, we
					* keep a shared pointer to
					* these entries so that more
					* than one accessor can access
					* this data if necessary.
					*/
            std_cxx0x::shared_ptr<const std::vector<unsigned int> > colnum_cache;
            
	                               /**
					* Discard the old row caches
					* (they may still be used by
					* other accessors) and
					* generate new ones for the
					* row pointed to presently by
					* this accessor.
					*/
            void visit_present_row ();

                                       /**
					* Make enclosing class a
					* friend.
					*/
            friend class const_iterator;
        };
        
      public:
          
                                       /**
					* Constructor. Create an
					* iterator into the matrix @p
					* matrix for the given row and
					* the index within it.
					*/ 
        const_iterator (const SparsityPattern *sparsity_pattern,
                        const unsigned int     row,
                        const unsigned int     index);
          
                                       /**
					* Prefix increment.
					*/
        const_iterator& operator++ ();

                                       /**
					* Postfix increment.
					*/
        const_iterator operator++ (int);

                                       /**
					* Dereferencing operator.
					*/
        const Accessor& operator* () const;

                                       /**
					* Dereferencing operator.
					*/
        const Accessor* operator-> () const;

                                       /**
					* Comparison. True, if both
					* iterators point to the same
					* matrix position.
					*/
        bool operator == (const const_iterator&) const;

                                       /**
					* Inverse of <tt>==</tt>.
					*/
        bool operator != (const const_iterator&) const;

                                       /**
					* Comparison operator. Result
					* is true if either the first
					* row number is smaller or if
					* the row numbers are equal
					* and the first index is
					* smaller.
					*/
        bool operator < (const const_iterator&) const;

	                               /**
                                        * Exception
					*/
        DeclException2 (ExcInvalidIndexWithinRow,
                        int, int,
                        << "Attempt to access element " << arg2
                        << " of row " << arg1
                        << " which doesn't have that many elements.");
        
      private:
                                       /**
                                        * Store an object of the
                                        * accessor class.
					*/
        Accessor accessor;
    };
    
  }
  
  
/**
 * This class implements a wrapper class to use the Trilinos distributed
 * sparsity pattern class Epetra_FECrsGraph. This class is designed to be
 * used for construction of parallel Trilinos matrices. The functionality of
 * this class is modeled after the existing sparsity pattern classes, with
 * the difference that this class can work fully in parallel according to a
 * partitioning of the sparsity pattern rows.
 *
 * This class has many similarities to the compressed sparsity pattern
 * classes of deal.II (i.e., the classes CompressedSparsityPattern,
 * CompressedSetSparsityPattern, and CompressedSimpleSparsityPattern), since
 * it can dynamically add elements to the pattern without any memory being
 * previously reserved for it. However, it also has a method
 * SparsityPattern::compress(), that finalizes the pattern and enables its
 * use with Trilinos sparse matrices.
 * 
 * @ingroup TrilinosWrappers
 * @ingroup Sparsity
 * @author Martin Kronbichler, 2008
 */
  class SparsityPattern : public Subscriptor
  {
    public:

                                       /**
                                        * Declare a typedef for the
                                        * iterator class.
                                        */
      typedef SparsityPatternIterators::const_iterator const_iterator;
      
/**
 * @name Constructors and initalization.
 */
//@{
                                       /**
                                        * Default constructor. Generates an
                                        * empty (zero-size) sparsity
                                        * pattern.
                                        */
      SparsityPattern ();

                                       /**
                                        * Constructor for a square sparsity
				        * pattern using an Epetra_Map and
				        * the number of nonzero entries in
				        * the rows of the sparsity
				        * pattern. Note that this number
				        * does not need to be exact, and it
				        * is even allowed that the actual
				        * sparsity structure has more
				        * nonzero entries than specified in
				        * the constructor. However it is
				        * still advantageous to provide good
				        * estimates here since this will
				        * considerably increase the
				        * performance when creating the
				        * sparsity pattern.
                                        */
      SparsityPattern (const Epetra_Map   &InputMap,
		       const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * Same as before, but now use the
				        * exact number of nonzeros in each m
				        * row. Since we know the number of
				        * elements in the sparsity pattern
				        * exactly in this case, we can
				        * already allocate the right amount
				        * of memory, which makes the
				        * creation process by the respective
				        * SparsityPattern::reinit call
				        * considerably faster. However, this
				        * is a rather unusual situation,
				        * since knowing the number of
				        * entries in each row is usually
				        * connected to knowing the indices
				        * of nonzero entries, which the
				        * sparsity pattern is designed to
				        * describe.
                                        */
      SparsityPattern (const Epetra_Map                &InputMap,
		       const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * This constructor is similar to the
                                        * one above, but it now takes two
                                        * different Epetra maps for rows and
                                        * columns. This interface is meant
                                        * to be used for generating
                                        * rectangular sparsity pattern,
                                        * where one map describes the
                                        * parallel partitioning of the dofs
                                        * associated with the sparsity
                                        * pattern rows and the other one of
                                        * the sparsity pattern columns. Note
                                        * that there is no real parallelism
                                        * along the columns &ndash; the
                                        * processor that owns a certain row
                                        * always owns all the column
                                        * elements, no matter how far they
                                        * might be spread out. The second
                                        * Epetra_Map is only used to specify
                                        * the number of columns and for
                                        * internal arragements when doing
                                        * matrix-vector products with
                                        * vectors based on that column map.
					*
					* The number of columns entries
					* per row is specified as the
					* maximum number of entries
					* argument.
                                        */
      SparsityPattern (const Epetra_Map   &InputRowMap,
		       const Epetra_Map   &InputColMap,
		       const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * This constructor is similar to the
				        * one above, but it now takes two
				        * different Epetra maps for rows and
				        * columns. This interface is meant
				        * to be used for generating
				        * rectangular matrices, where one
				        * map specifies the parallel
				        * distribution of rows and the
				        * second one specifies the
				        * distribution of degrees of freedom
				        * associated with matrix
				        * columns. This second map is
				        * however not used for the
				        * distribution of the columns
				        * themselves &ndash; rather, all
				        * column elements of a row are
				        * stored on the same processor. The
				        * vector <tt>n_entries_per_row</tt>
				        * specifies the number of entries in
				        * each row of the newly generated
				        * matrix.
                                        */
      SparsityPattern (const Epetra_Map                &InputRowMap,
		       const Epetra_Map                &InputColMap,
		       const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * Generate a sparsity pattern that
                                        * is completely stored locally,
                                        * having $m$ rows and $n$ columns. The
                                        * resulting matrix will be
                                        * completely stored locally.
					*
					* The number of columns entries per
					* row is specified as the maximum
					* number of entries argument. As
					* above, this does not need to be an
					* accurate number since the entries
					* are allocated dynamically in a
					* similar manner as for the deal.II
					* CompressedSparsityPattern classes,
					* but a good estimate will reduce
					* the setup time of the sparsity
					* pattern.
                                        */
      SparsityPattern (const unsigned int  m,
		       const unsigned int  n,
		       const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * Generate a sparsity pattern that
                                        * is completely stored locally,
                                        * having $m$ rows and $n$ columns. The
                                        * resulting matrix will be
                                        * completely stored locally.
					*
					* The vector
					* <tt>n_entries_per_row</tt>
					* specifies the number of entries in
					* each row.
                                        */
      SparsityPattern (const unsigned int               m,
		       const unsigned int               n,
		       const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * Copy constructor. Sets the calling
                                        * sparsity pattern to be the same as
                                        * the input sparsity pattern.
                                        */
      SparsityPattern (const SparsityPattern &SP);

                                       /**
                                        * Destructor. Made virtual so that
                                        * one can use pointers to this
                                        * class.
                                        */
      virtual ~SparsityPattern ();

                                       /**
                                        * Reinitialization function for
				        * generating a square sparsity
				        * pattern using an Epetra_Map and
				        * the number of nonzero entries in
				        * the rows of the sparsity
				        * pattern. Note that this number
				        * does not need to be exact, and it
				        * is even allowed that the actual
				        * sparsity structure has more
				        * nonzero entries than specified in
				        * the constructor. However it is
				        * still advantageous to provide good
				        * estimates here since this will
				        * considerably increase the
				        * performance when creating the
				        * sparsity pattern.
					*
					* This function does not create any
					* entries by itself, but provides
					* the correct data structures that
					* can be used by the respective
					* add() function.
                                        */
      void  
      reinit (const Epetra_Map   &InputMap,
	      const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * Same as before, but now use the
				        * exact number of nonzeros in each m
				        * row. Since we know the number of
				        * elements in the sparsity pattern
				        * exactly in this case, we can
				        * already allocate the right amount
				        * of memory, which makes process of
				        * adding entries to the sparsity
				        * pattern considerably
				        * faster. However, this is a rather
				        * unusual situation, since knowing
				        * the number of entries in each row
				        * is usually connected to knowing
				        * the indices of nonzero entries,
				        * which the sparsity pattern is
				        * designed to describe.
                                        */
      void  
      reinit (const Epetra_Map                &InputMap,
	      const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * This reinit function is similar to
                                        * the one above, but it now takes
                                        * two different Epetra maps for rows
                                        * and columns. This interface is
                                        * meant to be used for generating
                                        * rectangular sparsity pattern,
                                        * where one map describes the
                                        * parallel partitioning of the dofs
                                        * associated with the sparsity
                                        * pattern rows and the other one of
                                        * the sparsity pattern columns. Note
                                        * that there is no real parallelism
                                        * along the columns &ndash; the
                                        * processor that owns a certain row
                                        * always owns all the column
                                        * elements, no matter how far they
                                        * might be spread out. The second
                                        * Epetra_Map is only used to specify
                                        * the number of columns and for
                                        * internal arragements when doing
                                        * matrix-vector products with
                                        * vectors based on that column map.
					*
					* The number of columns entries per
					* row is specified by the argument
					* <tt>n_entries_per_row</tt>.
                                        */
      void  
      reinit (const Epetra_Map   &InputRowMap,
	      const Epetra_Map   &InputColMap,
	      const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * This reinit function is similar to
				        * the one above, but it now takes
				        * two different Epetra maps for rows
				        * and columns. This interface is
				        * meant to be used for generating
				        * rectangular matrices, where one
				        * map specifies the parallel
				        * distribution of rows and the
				        * second one specifies the
				        * distribution of degrees of freedom
				        * associated with matrix
				        * columns. This second map is
				        * however not used for the
				        * distribution of the columns
				        * themselves &ndash; rather, all
				        * column elements of a row are
				        * stored on the same processor. The
				        * vector <tt>n_entries_per_row</tt>
				        * specifies the number of entries in
				        * each row of the newly generated
				        * matrix.
                                        */
      void  
      reinit (const Epetra_Map                &InputRowMap,
	      const Epetra_Map                &InputColMap,
	      const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * Initialize a sparsity pattern that
                                        * is completely stored locally,
                                        * having $m$ rows and $n$ columns. The
                                        * resulting matrix will be
                                        * completely stored locally.
					*
					* The number of columns entries per
					* row is specified as the maximum
					* number of entries argument. As
					* above, this does not need to be an
					* accurate number since the entries
					* are allocated dynamically in a
					* similar manner as for the deal.II
					* CompressedSparsityPattern classes,
					* but a good estimate will reduce
					* the setup time of the sparsity
					* pattern.
                                        */
      void  
      reinit (const unsigned int  m,
	      const unsigned int  n,
	      const unsigned int  n_entries_per_row = 1);

                                       /**
                                        * Initialize a sparsity pattern that
                                        * is completely stored locally,
                                        * having $m$ rows and $n$ columns. The
                                        * resulting matrix will be
                                        * completely stored locally.
					*
					* The vector
					* <tt>n_entries_per_row</tt>
					* specifies the number of entries in
					* each row.
                                        */
      void  
      reinit (const unsigned int               m,
	      const unsigned int               n,
	      const std::vector<unsigned int> &n_entries_per_row);

                                       /**
                                        * Reinit function. Takes one of the
                                        * deal.II sparsity patterns and a
                                        * parallel partitioning of the rows
                                        * and columns for initializing the
                                        * current Trilinos sparsity pattern.
                                        */
      template<typename SparsityType>
      void  
      reinit (const Epetra_Map   &InputRowMap,
	      const Epetra_Map   &InputColMap,
	      const SparsityType &SP);

                                       /**
                                        * Reinit function. Takes one of the
                                        * deal.II sparsity patterns and a
                                        * parallel partitioning of the rows
                                        * and columns for initializing the
                                        * current Trilinos sparsity pattern.
                                        */
      template<typename SparsityType>
      void  
      reinit (const Epetra_Map   &InputMap,
	      const SparsityType &SP);

                                       /**
                                        * Copy function. Sets the calling
                                        * sparsity pattern to be the same as
                                        * the input sparsity pattern.
                                        */
      void  
      copy_from (const SparsityPattern &SP);

                                       /**
                                        * Copy function from one of the
                                        * deal.II sparsity patterns. If used
                                        * in parallel, this function uses an
                                        * ad-hoc partitioning of the rows
                                        * and columns.
                                        */
      template<typename SparsityType>
      void  
      copy_from (const SparsityType &SP);

                                       /**
                                        * Release all memory and
                                        * return to a state just like
                                        * after having called the
                                        * default constructor.
					*
					* This is a
				        * collective operation that
				        * needs to be called on all
				        * processors in order to avoid a
				        * dead lock.
                                        */
      void clear ();

                                       /**
                                        * In analogy to our own
                                        * SparsityPattern class, this
                                        * function compresses the sparsity
                                        * pattern and allows the resulting
                                        * pattern to be used for actually
                                        * generating a matrix. This function
                                        * also exchanges non-local data that
                                        * might have accumulated during the
                                        * addition of new elements. This
                                        * function must therefore be called
                                        * once the structure is fixed. This
                                        * is a collective operation, i.e.,
                                        * it needs to be run on all
                                        * processors when used in parallel.
                                        */
      void compress ();
//@}
/**
 * @name Information on the sparsity pattern
 */
//@{

				       /**
					* Returns the state of the sparsity
					* pattern, i.e., whether compress()
					* needs to be called after an
					* operation requiring data
					* exchange.
					*/
      bool is_compressed () const;

                                       /**
                                        * Gives the maximum number of
                                        * entries per row on the current
                                        * processor.
                                        */
      unsigned int max_entries_per_row () const;
      
                                       /**
                                        * Return the number of rows in this
                                        * sparsity pattern.
                                        */
      unsigned int n_rows () const;

                                       /**
                                        * Return the number of columns in
                                        * this sparsity pattern.
                                        */
      unsigned int n_cols () const;

                                       /**
                                        * Return the local dimension of the
                                        * sparsity pattern, i.e. the number
                                        * of rows stored on the present MPI
                                        * process. In the sequential case,
                                        * this number is the same as
                                        * n_rows(), but for parallel
                                        * matrices it may be smaller.
					*
					* To figure out which elements
					* exactly are stored locally,
					* use local_range().
                                        */
      unsigned int local_size () const;

                                       /**
					* Return a pair of indices
					* indicating which rows of this
					* sparsity pattern are stored
					* locally. The first number is the
					* index of the first row stored, the
					* second the index of the one past
					* the last one that is stored
					* locally. If this is a sequential
					* matrix, then the result will be
					* the pair (0,n_rows()), otherwise
					* it will be a pair (i,i+n), where
					* <tt>n=local_size()</tt>.
					*/
      std::pair<unsigned int, unsigned int>
	local_range () const;

				       /**
					* Return whether @p index is
					* in the local range or not,
					* see also local_range().
					*/
      bool in_local_range (const unsigned int index) const;

                                       /**
                                        * Return the number of nonzero
                                        * elements of this sparsity pattern.
                                        */
      unsigned int n_nonzero_elements () const;

                                       /**
                                        * Number of entries in a
                                        * specific row.
                                        */
      unsigned int row_length (const unsigned int row) const;

                                       /**
					* Compute the bandwidth of the
					* matrix represented by this
					* structure. The bandwidth is the
					* maximum of $|i-j|$ for which the
					* index pair $(i,j)$ represents a
					* nonzero entry of the
					* matrix. Consequently, the maximum
					* bandwidth a $n\times m$ matrix can
					* have is $\max\{n-1,m-1\}$.
                                        */
      unsigned int bandwidth () const;

                                       /**
                                        * Return whether the object is
                                        * empty. It is empty if no memory is
                                        * allocated, which is the same as
                                        * when both dimensions are zero.
                                        */
      bool empty () const;

                                       /**
                                        * Return whether the index
                                        * (<i>i,j</i>) exists in the
                                        * sparsity pattern (i.e., it may be
                                        * non-zero) or not.
                                        */
      bool exists (const unsigned int i,
		   const unsigned int j) const;

				       /**
					* Determine an estimate for the
					* memory consumption (in bytes)
					* of this object. Currently not
					* implemented for this class.
					*/
      unsigned int memory_consumption () const;

//@}
/**
 * @name Adding entries
 */
//@{
                                       /**
                                        * Add the element (<i>i,j</i>) to
                                        * the sparsity pattern.
                                        */
      void add (const unsigned int i,
                const unsigned int j);


                                       /**
                                        * Add several elements in one row to
                                        * the sparsity pattern.
                                        */
      void add (const unsigned int    row,
		const unsigned int    n_cols,
		const unsigned int   *col_indices);
//@}
/**
 * @name Iterators
 */
//@{

                                       /**
                                        * STL-like iterator with the
                                        * first entry.
                                        */
      const_iterator begin () const;

                                       /**
                                        * Final iterator.
                                        */
      const_iterator end () const;

                                       /**
                                        * STL-like iterator with the
                                        * first entry of row @p r.
                                        *
                                        * Note that if the given row
                                        * is empty, i.e. does not
                                        * contain any nonzero entries,
                                        * then the iterator returned
                                        * by this function equals
                                        * <tt>end(r)</tt>. Note also
                                        * that the iterator may not be
                                        * dereferencable in that case.
                                        */
      const_iterator begin (const unsigned int r) const;

                                       /**
                                        * Final iterator of row
                                        * <tt>r</tt>. It points to the
                                        * first element past the end
                                        * of line @p r, or past the
                                        * end of the entire sparsity
                                        * pattern.
                                        *
                                        * Note that the end iterator
                                        * is not necessarily
                                        * dereferencable. This is in
                                        * particular the case if it is
                                        * the end iterator for the
                                        * last row of a matrix.
                                        */
      const_iterator end (const unsigned int r) const;

//@}
/**
 * @name Input/Output
 */
//@{

                                       /**
					* Abstract Trilinos object
					* that helps view in ASCII
					* other Trilinos
					* objects. Currently this
					* function is not
					* implemented.  TODO: Not
					* implemented.
					*/
      void write_ascii ();

				       /**
					* Print the sparsity pattern to the
					* given stream, using the format
					* <tt>(line,col)</tt>.
					*/
      void print (std::ostream &out) const;

				       /**
					* Print the sparsity of the matrix
					* in a format that <tt>gnuplot</tt>
					* understands and which can be used
					* to plot the sparsity pattern in a
					* graphical way. The format consists
					* of pairs <tt>i j</tt> of nonzero
					* elements, each representing one
					* entry of this matrix, one per line
					* of the output file. Indices are
					* counted from zero on, as
					* usual. Since sparsity patterns are
					* printed in the same way as
					* matrices are displayed, we print
					* the negative of the column index,
					* which means that the
					* <tt>(0,0)</tt> element is in the
					* top left rather than in the bottom
					* left corner.
					*
					* Print the sparsity pattern in
					* gnuplot by setting the data style
					* to dots or points and use the
					* <tt>plot</tt> command.
					*/
      void print_gnuplot (std::ostream &out) const;
    
                                        // TODO: Write an overloading
                                        // of the operator << for output.
                                        // Since the underlying Trilinos 
                                        // object supports it, this should 
                                        // be very easy.

//@}
    				     /** @addtogroup Exceptions
				      * @{ */
                                       /**
                                        * Exception
                                        */
      DeclException1 (ExcTrilinosError,
                      int,
                      << "An error with error number " << arg1
                      << " occured while calling a Trilinos function");

				       /**
					* Exception
					*/
      DeclException2 (ExcInvalidIndex,
		      int, int,
		      << "The entry with index <" << arg1 << ',' << arg2
		      << "> does not exist.");

                                       /**
                                        * Exception
                                        */
      DeclException0 (ExcSourceEqualsDestination);

                                       /**
                                        * Exception
                                        */
      DeclException4 (ExcAccessToNonLocalElement,
		      int, int, int, int,
		      << "You tried to access element (" << arg1
		      << "/" << arg2 << ")"
		      << " of a distributed matrix, but only rows "
		      << arg3 << " through " << arg4
		      << " are stored locally and can be accessed.");

                                       /**
                                        * Exception
                                        */
      DeclException2 (ExcAccessToNonPresentElement,
		      int, int,
		      << "You tried to access element (" << arg1
		      << "/" << arg2 << ")"
		      << " of a sparse matrix, but it appears to not"
		      << " exist in the Trilinos sparsity pattern.");
				     //@}    
    private:
                                       /**
				        * Epetra Trilinos
				        * mapping of the matrix rows
				        * that assigns parts of the
				        * matrix to the individual
				        * processes. This map is
				        * provided either via the
				        * constructor or in a reinit
				        * function.
				        */
      Epetra_Map row_map;

                                       /**
                                        * Pointer to the user-supplied
				        * Epetra Trilinos mapping of
				        * the matrix columns that
				        * assigns parts of the matrix
				        * to the individual processes.
				        */
      Epetra_Map col_map;

				       /**
					* A boolean variable to hold
					* information on whether the
					* vector is compressed or not.
					*/
      bool compressed;

                                       /**
                                        * A sparsity pattern object in
                                        * Trilinos to be used for finite
                                        * element based problems which
                                        * allows for adding non-local
                                        * elements to the pattern.
                                        */
      std::auto_ptr<Epetra_FECrsGraph> graph;

                                       /**
					* Scratch array that holds several
					* indices that should be written
					* into the same row of the sparsity
					* pattern. This is to increase the
					* speed of this function.
					*/
      std::vector<unsigned int> cached_row_indices;

				       /**
					* A number that tells how many
					* indices currently are active in
					* the cache.
					*/
      unsigned int n_cached_elements;

				       /**
					* The row that is currently in the
					* cache.
					*/
      unsigned int row_in_cache;

      friend class SparseMatrix;
      friend class SparsityPatternIterators::const_iterator;
  };



// -------------------------- inline and template functions ----------------------


#ifndef DOXYGEN

  namespace SparsityPatternIterators
  {

    inline
    const_iterator::Accessor::
    Accessor (const SparsityPattern *sp,
              const unsigned int     row,
              const unsigned int     index)
                    :
                    sparsity_pattern(const_cast<SparsityPattern*>(sp)),
                    a_row(row),
                    a_index(index)
    {
      visit_present_row ();
    }


    inline
    unsigned int
    const_iterator::Accessor::row() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return a_row;
    }



    inline
    unsigned int
    const_iterator::Accessor::column() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return (*colnum_cache)[a_index];
    }



    inline
    unsigned int
    const_iterator::Accessor::index() const
    {
      Assert (a_row < sparsity_pattern->n_rows(), ExcBeyondEndOfSparsityPattern());
      return a_index;
    }



    inline
    const_iterator::
    const_iterator(const SparsityPattern *sp,
                   const unsigned int     row,
                   const unsigned int     index)
                    :
                    accessor(sp, row, index)
    {}



    inline
    const_iterator &
    const_iterator::operator++ ()
    {
      Assert (accessor.a_row < accessor.sparsity_pattern->n_rows(), 
	      ExcIteratorPastEnd());

      ++accessor.a_index;

				        // If at end of line: do one
				        // step, then cycle until we 
				        // find a row with a nonzero
				        // number of entries.
      if (accessor.a_index >= accessor.colnum_cache->size())
        {
          accessor.a_index = 0;
          ++accessor.a_row;
      
          while ((accessor.a_row < accessor.sparsity_pattern->n_rows())
                 &&
                 (accessor.sparsity_pattern->row_length(accessor.a_row) == 0))
            ++accessor.a_row;

          accessor.visit_present_row();
        }
      return *this;
    }



    inline
    const_iterator
    const_iterator::operator++ (int)
    {
      const const_iterator old_state = *this;
      ++(*this);
      return old_state;
    }



    inline
    const const_iterator::Accessor &
    const_iterator::operator* () const
    {
      return accessor;
    }



    inline
    const const_iterator::Accessor *
    const_iterator::operator-> () const
    {
      return &accessor;
    }



    inline
    bool
    const_iterator::
    operator == (const const_iterator& other) const
    {
      return (accessor.a_row == other.accessor.a_row &&
              accessor.a_index == other.accessor.a_index);
    }



    inline
    bool
    const_iterator::
    operator != (const const_iterator& other) const
    {
      return ! (*this == other);
    }



    inline
    bool
    const_iterator::
    operator < (const const_iterator& other) const
    {
      return (accessor.row() < other.accessor.row() ||
              (accessor.row() == other.accessor.row() &&
               accessor.index() < other.accessor.index()));
    }
    
  }
 


  inline
  SparsityPattern::const_iterator
  SparsityPattern::begin() const
  {
    return const_iterator(this, 0, 0);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::end() const
  {
    return const_iterator(this, n_rows(), 0);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::begin(const unsigned int r) const
  {
    Assert (r < n_rows(), ExcIndexRange(r, 0, n_rows()));
    if (row_length(r) > 0)
      return const_iterator(this, r, 0);
    else
      return end (r);
  }



  inline
  SparsityPattern::const_iterator
  SparsityPattern::end(const unsigned int r) const
  {
    Assert (r < n_rows(), ExcIndexRange(r, 0, n_rows()));

                                     // place the iterator on the first entry
                                     // past this line, or at the end of the
                                     // matrix
    for (unsigned int i=r+1; i<n_rows(); ++i)
      if (row_length(i) > 0)
        return const_iterator(this, i, 0);
    
                                     // if there is no such line, then take the
                                     // end iterator of the matrix
    return end();
  }



  inline
  bool
  SparsityPattern::in_local_range (const unsigned int index) const
  {
    int begin, end;
    begin = graph->RowMap().MinMyGID();
    end = graph->RowMap().MaxMyGID()+1;
    
    return ((index >= static_cast<unsigned int>(begin)) &&
            (index < static_cast<unsigned int>(end)));
  }



  inline
  bool
  SparsityPattern::is_compressed () const
  {
    return compressed;
  }



  inline
  bool
  SparsityPattern::empty () const
  {
    return ((n_rows() == 0) && (n_cols() == 0));
  }



  inline
  void
  SparsityPattern::add (const unsigned int i,
			const unsigned int j)
  {
				   // if we want to write an element to the
				   // row the cache is currently pointed to,
				   // we just append the data to the cache
    if (i == row_in_cache)
      {
				   // if the size is too small, extend the
				   // cache by 10 elements
	if (n_cached_elements >= cached_row_indices.size())
	  cached_row_indices.resize(cached_row_indices.size() + 10);

	cached_row_indices[n_cached_elements] = j;
	++n_cached_elements;
	return;
      }

				   // if we are to write another row data,
				   // we write the cache data into the
				   // sparsity pattern, and then call this
				   // function again
    add (row_in_cache, n_cached_elements, &cached_row_indices[0]);
    row_in_cache = i;
    n_cached_elements = 0;
    add (i,j);
  }



  inline
  void
  SparsityPattern::add (const unsigned int    row,
			const unsigned int    n_cols,
			const unsigned int   *col_indices)
  {
    if (n_cols == 0)
      return;

    int * col_index_ptr = (int*)col_indices;
    compressed = false;

    int ierr;

				   // If the calling sparsity pattern owns
				   // the row to which we want to add
				   // values, we can directly call the
				   // Epetra_CrsGraph input function, which
				   // is much faster than the
				   // Epetra_FECrsGraph function.
    if (row_map.MyGID(row) == true)
      ierr = graph->Epetra_CrsGraph::InsertGlobalIndices(row, 
							 n_cols,
							 col_index_ptr);
    else
      {
				   // When we're at off-processor data, we
				   // have to stick with the standard
				   // SumIntoGlobalValues
				   // function. Nevertheless, the way we
				   // call it is the fastest one (any other
				   // will lead to repeated allocation and
				   // deallocation of memory in order to
				   // call the function we already use,
				   // which is very unefficient if writing
				   // one element at a time).

	ierr = graph->InsertGlobalIndices (1, (int*)&row, n_cols, 
					   col_index_ptr);
      }

    //Assert (ierr <= 0, ExcAccessToNonPresentElement(row, col_index_ptr[0]));
    AssertThrow (ierr >= 0, ExcTrilinosError(ierr));
  }


#endif // DOXYGEN      
}


DEAL_II_NAMESPACE_CLOSE


#endif // DEAL_II_USE_TRILINOS


/*--------------------   trilinos_sparsity_pattern.h     --------------------*/

#endif
/*--------------------   trilinos_sparsity_pattern.h     --------------------*/
