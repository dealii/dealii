/*----------------------------   sparsematrix.h     ---------------------------*/
/*      $Id$                 */
#ifndef __sparsematrix_H
#define __sparsematrix_H
/*----------------------------   sparsematrix.h     ---------------------------*/


// This file was once part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier
// Revised, modified and extended by Wolfgang Bangerth


#include <base/exceptions.h>


//forward declarations
template <typename number> class Vector;
template <typename number> class SparseMatrix;

class iVector;
class ostream;



/*

   @author Original version by Roland Becker, Guido Kanschat, Franz-Theo Suttmeier; lots of enhancements, reorganisation and documentation by Wolfgang Bangerth
   */
class SparseMatrixStruct
{
  public:
				     /**
				      * Initialize the matrix empty, i.e. with
				      * no memory allocated. This is useful if
				      * you want such objects as member
				      * variables in other classes. You can make
				      * the structure usable by calling the
				      * #reinit# function.
				      */
    SparseMatrixStruct ();
    
				     /**
				      * Copy constructor. This constructor is
				      * only allowed to be called if the matrix
				      * structure to be copied is empty. This is
				      * so in order to prevent involuntary
				      * copies of objects for temporaries, which
				      * can use large amounts of computing time.
				      * However, copy constructors are needed
				      * if yo want to use the STL data types
				      * on classes like this, e.g. to write
				      * such statements like
				      * #v.push_back (SparseMatrixStruct());#,
				      * with #v# a vector of #SparseMatrixStruct#
				      * objects.
				      *
				      * Usually, it is sufficient to use the
				      * explicit keyword to disallow unwanted
				      * temporaries, but for the STL vectors,
				      * this does not work. Since copying a
				      * structure like this is not useful
				      * anyway because multiple matrices can
				      * use the same sparsity structure, copies
				      * are only allowed for empty objects, as
				      * described above.
				      */
    SparseMatrixStruct (const SparseMatrixStruct &);

				     /**
				      * Initialize a rectangular matrix with
				      * #m# rows and #n# columns.
				      * The matrix may contain at most #max_per_row#
				      * nonzero entries per row.
				      */
    SparseMatrixStruct (const unsigned int m,
			const unsigned int n,
			const unsigned int max_per_row);
    
				     /**
				      * Initialize a square matrix of dimension
				      * #n# with at most #max_per_row#
				      * nonzero entries per row.
				      */
    SparseMatrixStruct (const unsigned int n,
			const unsigned int max_per_row);

				     /**
				      * Make a copy with extra off-diagonals.
				      *
				      * This constructs objects intended for
				      * the application of the ILU(n)-method
				      * or other incomplete decompositions.
				      * Therefore, additional to the original
				      * entry structure, space for
				      * #extra_off_diagonals#
				      * side-diagonals is provided on both
				      * sides of the main diagonal.
				      *
				      * #max_per_row# is the maximum number of
				      * nonzero elements per row which this
				      * structure is to hold. It is assumed
				      * that this number is sufficiently large
				      * to accomodate both the elements in
				      * #original# as well as the new
				      * off-diagonal elements created by this
				      * constructor. You will usually want to
				      * give the same number as you gave for
				      * #original# plus the number of side
				      * diagonals times two. You may however
				      * give a larger value if you wish to add
				      * further nonzero entries for the
				      * decomposition based on other criteria
				      * than their being on side-diagonals.
				      *
				      * This function requires that #original#
				      * refer to a square matrix structure.
				      * It shall be compressed. The matrix 
				      * structure is not compressed
				      * after this function finishes.
				      */
    SparseMatrixStruct(const SparseMatrixStruct& original,
		       const unsigned int        max_per_row,
		       const unsigned int        extra_off_diagonals);
    
				     /**
				      * Destructor.
				      */
    ~SparseMatrixStruct ();
    
				     /**
				      * Reallocate memory and set up data
				      * structures for a new matrix with
				      * #m# rows and #n# columns,
				      * with at most #max_per_row#
				      * nonzero entries per row.
				      *
				      * If #m*n==0# all memory is freed,
				      * resulting in a total reinitialization
				      * of the object. If it is nonzero, new
				      * memory is only allocated if the new
				      * size extends the old one. This is done
				      * to save time and to avoid fragmentation
				      * of the heap.
				      */
    void reinit (const unsigned int m,
		 const unsigned int n,
		 const unsigned int max_per_row);

				     /**
				      * This function compresses the sparsity
				      * structure that this object represents.
				      * It does so by eliminating unused
				      * entries and sorting the remaining
				      * ones to allow faster access by usage
				      * of binary search algorithms. A special
				      * sorting scheme is used for the diagonal
				      * entry of square matrices, which is
				      * always the first entry of each row.
				      *
				      * #SparseMatrix# objects require the
				      * #SparseMatrixStruct# objects they are
				      * initialized with to be compressed, to
				      * reduce memory requirements.
				      */
    void compress ();

				     /**
				      * Return whether the object is empty. It
				      * is empty if no memory is allocated,
				      * which is the same as that both
				      * dimensions are zero.
				      */
    bool empty () const;
    

				     /**
				      * Return the index of the matrix
				      * element with row number #i# and
				      * column number #j#. If the matrix
				      * element is not a nonzero one,
				      * return -1.
				      *
				      * This function is usually called
				      * by the #operator()# of the
				      * #SparseMatrix#. It shall only be
				      * called for compressed sparsity
				      * patterns, since in this case
				      * searching whether the entry
				      * exists can be done quite fast
				      * with a binary sort algorithm
				      * because the column numbers are
				      * sorted.
				      */
    int operator() (const unsigned int i, const unsigned int j) const;

				     /**
				      * Add a nonzero entry to the matrix.
				      * This function may only be called
				      * for non-compressed sparsity patterns.
				      *
				      * If the entry already exists, nothing
				      * bad happens.
				      */
    void add (const unsigned int i, const unsigned int j);
    
				     /**
				      * This matrix adds a whole connectivity
				      * list to the sparsity structure
				      * respresented by this object. It assumes
				      * the #rowcols# array to be a list of
				      * indices which are all linked together,
				      * i.e. all entries
				      * #(rowcols[i], rowcols[j])# for all
				      * #i,j=0...n# for this sparsity pattern
				      * are created. #n# is assumed to be the
				      * number of elements pointed to by
				      * #rowcols#.
				      */
    void add_matrix (const unsigned int n, const int* rowcols);

				     //////////
    void add_matrix (const unsigned int m, const unsigned int n,
		     const int* rows, const int* cols);
				     //////////
    void add_matrix (const iVector& rowcols);
				     //////////
    void add_matrix (const iVector& rows, const iVector& cols);

				     /**
				      * Print the sparsity of the matrix
				      * in a format that #gnuplot# understands
				      * and which can be used to plot the
				      * sparsity pattern in a graphical
				      * way. The format consists of pairs
				      * #i j# of nonzero elements, each
				      * representing one entry of this
				      * matrix, one per line of the output
				      * file. Indices are counted from
				      * zero on, as usual. Since sparsity
				      * patterns are printed in the same
				      * way as matrices are displayed, we
				      * print the negative of the column
				      * index, which means that the
				      * #(0,0)# element is in the top left
				      * rather than in the bottom left
				      * corner.
				      *
				      * Print the sparsity pattern in
				      * gnuplot by setting the data style
				      * to dots or points and use the
				      * #plot# command.
				      */
    void print_gnuplot (ostream &out) const;

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
				      * Compute the bandwidth of the matrix
				      * represented by this structure. The
				      * bandwidth is the maximum of
				      * $|i-j|$ for which the index pair
				      * $(i,j)$ represents a nonzero entry
				      * of the matrix.
				      */
    unsigned int bandwidth () const;

				     /**
				      * Return the number of nonzero elements of
				      * this matrix. Actually, it returns the
				      * number of entries in the sparsity
				      * pattern; if any of the entries should
				      * happen to be zero, it is counted
				      * anyway.
				      *
				      * This function may only be called if the
				      * matrix struct is compressed. It does not
				      * make too much sense otherwise anyway.
				      */
    unsigned int n_nonzero_elements () const;

				     /**
				      * Return whether the structure is
				      * compressed or not.
				      */
    bool is_compressed () const;
    
				     /**
				      * This is kind of an expert mode: get
				      * access to the rowstart array, but
				      * readonly.
				      *
				      * Though the return value is declared
				      * #const#, you should be aware that it
				      * may change if you call any nonconstant
				      * function of objects which operate on
				      * it.
				      *
				      * You should use this interface very
				      * carefully and only if you are absolutely
				      * sure to know what you do. You should
				      * also note that the structure of these
				      * arrays may change over time.
				      * If you change the layout yourself, you
				      * should also rename this function to
				      * avoid programs relying on outdated
				      * information!
				      */
    const unsigned int * get_rowstart_indices () const;

				     /**
				      * This is kind of an expert mode: get
				      * access to the colnums array, but
				      * readonly.
				      *
				      * Though the return value is declared
				      * #const#, you shoudl be aware that it
				      * may change if you call any nonconstant
				      * function of objects which operate on
				      * it.
				      *
				      * You should use this interface very
				      * carefully and only if you are absolutely
				      * sure to know what you do. You should
				      * also note that the structure of these
				      * arrays may change over time.
				      * If you change the layout yourself, you
				      * should also rename this function to
				      * avoid programs relying on outdated
				      * information!
				      */
    const int * get_column_numbers () const;
    
    
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
		    << ": there was no free entry any more. " << endl
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
    DeclException0 (ExcInternalError);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCall);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotSquare);
    
  private:
    unsigned int max_dim;
    unsigned int rows, cols;
    unsigned int vec_len, max_vec_len;
    unsigned int max_row_len;
    unsigned int* rowstart;
    int* colnums;

				     /**
				      * Store whether the #compress# function
				      * was called for this object.
				      */
    bool compressed;
    
    template <typename number> friend class SparseMatrix;
};




/*
CLASS
   SparseMatrix

   @author Original version by Roland Becker, Guido Kanschat, Franz-Theo Suttmeier; lots of enhancements, reorganisation and documentation by Wolfgang Bangerth 1998
   */
template <typename number>
class SparseMatrix
{
  public:
				     /**
				      * Type of matrix entries.
				      */
    typedef number entry_type;
    
				     /**
				      * Constructor; initializes the matrix to
				      * be empty, without any structure, i.e.
				      * the matrix is not usable at all. This
				      * constructor is therefore only useful
				      * for matrices which are members of a
				      * class. All other matrices should be
				      * created at a point in the data flow
				      * where all necessary information is
				      * available.
				      *
				      * You have to initialize
				      * the matrix before usage with
				      * #reinit(SparseMatrixStruct)#.
				      */
    SparseMatrix ();

				     /**
				      * Copy constructor. This constructor is
				      * only allowed to be called if the matrix
				      * to be copied is empty. This is for the
				      * same reason as for the
				      * #SparseMatrixStruct#, see there for the
				      * details.
				      *
				      * If you really want to copy a whole
				      * matrix, you can do so by using the
				      * #copy_from# function.
				      */
    SparseMatrix (const SparseMatrix &);
    
    
				     /**
				      * Constructor. Takes the given matrix
				      * sparsity structure to represent the
				      * sparsity pattern of this matrix. You
				      * can change the sparsity pattern later
				      * on by calling the #reinit# function.
				      *
				      * You have to make sure that the lifetime
				      * of the sparsity structure is at least
				      * as long as that of this matrix or as
				      * long as #reinit# is not called with a
				      * new sparsity structure.
				      */
    SparseMatrix (const SparseMatrixStruct &sparsity);
    
				     /**
				      * Destructor. Free all memory, but do not
				      * release the memory of the sparsity
				      * structure.
				      */
    virtual ~SparseMatrix ();
    

				     /**
				      * Reinitialize the object but keep to
				      * the sparsity pattern previously used.
				      * This may be necessary if you #reinit#'d
				      * the sparsity structure and want to
				      * update the size of the matrix.
				      *
				      * Note that memory is only reallocated if
				      * the new size exceeds the old size. If
				      * that is not the case, the allocated
				      * memory is not reduced. However, if the
				      * sparsity structure is empty (i.e. the
				      * dimensions are zero), then all memory
				      * is freed.
				      */
    virtual void reinit ();

				     /**
				      * Reinitialize the sparse matrix with the
				      * given sparsity pattern. The latter tells
				      * the matrix how many nonzero elements
				      * there need to be reserved.
				      *
				      * Regarding memory allocation, the same
				      * applies as said above.
				      *
				      * You have to make sure that the lifetime
				      * of the sparsity structure is at least
				      * as long as that of this matrix or as
				      * long as #reinit# is not called with a
				      * new sparsity structure.
				      */
    virtual void reinit (const SparseMatrixStruct &sparsity);

				     /**
				      * Release all memory and return to a state
				      * just like after having called the
				      * default constructor. It also forgets the
				      * sparsity pattern it was previously tied
				      * to.
				      */
    virtual void clear ();
    
				     /**
				      * Return whether the object is empty. It
				      * is empty if either both dimensions
				      * are zero or no #SparseMatrixStruct#
				      * is associated.
				      */
    bool empty () const;

				     /**
				      * Return the dimension of the image space.
				      * To remember: the matrix is of dimension
				      * $m \times n$.
				      */
    unsigned int m () const;
    
				     /**
				      * Return the dimension of the range space.
				      * To remember: the matrix is of dimension
				      * $m \times n$.
				      */

    unsigned int n () const;

				     /**
				      * Return the number of nonzero elements of
				      * this matrix. Actually, it returns the
				      * number of entries in the sparsity
				      * pattern; if any of the entries should
				      * happen to be zero, it is counted
				      * anyway.
				      */
    unsigned int n_nonzero_elements () const;
    
				     /**
				      * Set the element #(i,j)# to #value#.
				      * Throws an error if the entry does
				      * not exist. Still, it is allowed to store
				      * zero values in non-existent fields.
				      */
    void set (const unsigned int i, const unsigned int j,
	      const number value);
    
				     /**
				      * Add #value# to the element #(i,j)#.
				      * Throws an error if the entry does
				      * not exist. Still, it is allowed to store
				      * zero values in non-existent fields.
				      */
    void add (const unsigned int i, const unsigned int j,
	      const number value);

				     /**
				      * Copy the given matrix to this one.
				      * The operation throws an error if the
				      * sparsity patterns of the two involved
				      * matrices do not point to the same
				      * object, since in this case the copy
				      * operation is cheaper. Since this
				      * operation is notheless not for free,
				      * we do not make it available through
				      * #operator =#, since this may lead
				      * to unwanted usage, e.g. in copy
				      * arguments to functions, which should
				      * really be arguments by reference.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      *
				      * The function returns a reference to
				      * #this#.
				      */
    template <typename somenumber>
    SparseMatrix<number> & copy_from (const SparseMatrix<somenumber> &source);

				     /**
				      * Generate ILU.
				      * The matrix entries will contain the
				      * incomplete LU-factorization of
				      * the matrix #source#. Having a matrix
				      * #source# and a matrix structure object
				      * #struct#, the code for generating
				      * an ILU factorization reads
				      * \begin{verbatim}
				      * SparseMatrix<float> ilu(struct);
				      * ilu.ILU(source);
				      * \end{verbatim}
				      *
				      * If additional side diagonals are
				      * needed (ILU(n)-algorithm), you have to
				      * construct a second matrix structure:
				      * \begin{verbatim}
				      * SparseMatrixStruct ilustruct(struct,n);
				      * SparseMatrix<float> ilu(ilustruct);
				      * ilu.ILU(source);
				      * \end{verbatim}
				      *
				      * After generating the ILU-decomposition,
				      * it can be applied to a vector by
				      * #backward_forward#.
				      */
    template <typename somenumber>
    SparseMatrix<number> & ILU (const SparseMatrix<somenumber> &source);

				     /**
				      * Add #matrix# scaled by #factor# to this
				      * matrix. The function throws an error
				      * if the sparsity patterns of the two
				      * involved matrices do not point to the
				      * same object, since in this case the
				      * operation is cheaper.
				      *
				      * The source matrix may be a matrix
				      * of arbitrary type, as long as its
				      * data type is convertible to the
				      * data type of this matrix.
				      */
    template <typename somenumber>
    void add_scaled (const number factor,
		     const SparseMatrix<somenumber> &matrix);
    
				     /**
				      * Return the value of the entry (i,j).
				      * This may be an expensive operation
				      * and you should always take care
				      * where to call this function.
				      * In order to avoid abuse, this function
				      * throws an exception if the wanted
				      * element does not exist in the matrix.
				      */
    number operator () (const unsigned int i, const unsigned int j) const;

				     /**
				      * Return the main diagonal element in
				      * the #i#th row. This function throws an
				      * error if the matrix is not square.
				      *
				      * This function is considerably faster
				      * than the #operator()#, since for
				      * square matrices, the diagonal entry is
				      * always the first to be stored in each
				      * row and access therefore does not
				      * involve searching for the right column
				      * number.
				      */
    number diag_element (const unsigned int i) const;

				     /**
				      * Same as above, but return a writeable
				      * reference. You're sure what you do?
				      */
    number & diag_element (const unsigned int i);
    
    				     /**
				      * This is kind of an expert mode: get
				      * access to the #i#th element of this
				      * matrix. The elements are stored in
				      * a consecutive way, refer to the
				      * #SparseMatrixStruct# class for more details.
				      *
				      * You should use this interface very
				      * carefully and only if you are absolutely
				      * sure to know what you do. You should
				      * also note that the structure of these
				      * arrays may change over time.
				      * If you change the layout yourself, you
				      * should also rename this function to
				      * avoid programs relying on outdated
				      * information!
				      */
    number global_entry (const unsigned int i) const;

				     /**
				      * Same as above, but with write access.
				      * You certainly know what you do?
				      */
    number & global_entry (const unsigned int i);

				     /**
				      * Matrix-vector multiplication: let
				      * $dst = M*src$ with $M$ being this matrix.
				      */
    template <typename somenumber>
    void vmult (Vector<somenumber>& dst, const Vector<somenumber>& src) const;
    
				     /**
				      * Matrix-vector multiplication: let
				      * $dst = M^T*src$ with $M$ being this
				      * matrix. This function does the same as
				      * #vmult# but takes the transposed matrix.
				      */
    template <typename somenumber>
    void Tvmult (Vector<somenumber>& dst, const Vector<somenumber>& src) const;
  
				     /**
				      * Do backward-forward solution of a
				      * previously generated ILU.
				      */
    template <typename somenumber>
    void backward_forward (Vector<somenumber>& v);
    
				     /**
				      * Return the norm of the vector $v$ with
				      * respect to the norm induced by this
				      * matrix, i.e. $\left(v,Mv\right)$. This
				      * is useful, e.g. in the finite element
				      * context, where the $L_2$ norm of a
				      * function equals the matrix norm with
				      * respect to the mass matrix of the vector
				      * representing the nodal values of the
				      * finite element function.
				      *
				      * Note the order in which the matrix
				      * appears. For non-symmetric matrices
				      * there is a difference whether the
				      * matrix operates on the first
				      * or on the second operand of the
				      * scalar product.
				      *
				      * Obviously, the matrix needs to be square
				      * for this operation.
				      */
    template <typename somenumber>
    double matrix_norm (const Vector<somenumber> &v) const;

    				     /**
				      * Return the l1-norm of the matrix, i.e.
				      * $|M|_1=max_{all columns j}\sum_{all 
				      * rows i} |M_ij|$,
				      * (max. sum of columns).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the l1-norm for vectors, i.e.
				      * $|Mv|_1\leq |M|_1 |v|_1$.
				      * (cf. Rannacher Numerik0)
				      */
    number l1_norm () const;

    				     /**
				      * Return the linfty-norm of the
				      * matrix, i.e.
				      * $|M|_infty=max_{all rows i}\sum_{all 
				      * columns j} |M_ij|$,
				      * (max. sum of rows).
				      * This is the
				      * natural matrix norm that is compatible
				      * to the linfty-norm of vectors, i.e.
				      * $|Mv|_infty \leq |M|_infty |v|_infty$.
				      * (cf. Rannacher Numerik0)
				      */
    number linfty_norm () const;

				     /**
				      * Compute the residual of an equation
				      * #Ax=b#, where the residual is defined
				      * to be #r=b-Ax# with #x# typically being
				      * an approximate of the true solution of
				      * the equation. Write the residual into
				      * #dst#.
				      */
    template <typename somenumber>
    double residual (Vector<somenumber>& dst, const Vector<somenumber>& x,
		     const Vector<somenumber>& b) const;
				     //
    template <typename somenumber>
    void precondition_Jacobi (Vector<somenumber>& dst, const Vector<somenumber>& src,
			      const number om = 1.) const;
				     //
    template <typename somenumber>
    void precondition_SSOR (Vector<somenumber>& dst, const Vector<somenumber>& src,
			    const number om = 1.) const;
				     //
    template <typename somenumber>
    void precondition_SOR (Vector<somenumber>& dst, const Vector<somenumber>& src,
			   const number om = 1.) const;
				     //
    template <typename somenumber>
    void SSOR (Vector<somenumber>& dst, const number om = 1.) const;
				     //
    template <typename somenumber>
    void SOR (Vector<somenumber>& dst, const number om = 1.) const;

				     /**
				      * Return a (constant) reference to the
				      * underlying sparsity pattern of this
				      * matrix.
				      *
				      * Though the return value is declared
				      * #const#, you shoudl be aware that it
				      * may change if you call any nonconstant
				      * function of objects which operate on
				      * it.
				      */
    const SparseMatrixStruct & get_sparsity_pattern () const;

				     /**
				      * Print the matrix to the given stream,
				      * using the format
				      * #(line,col) value#, i.e. one
				      * nonzero entry of the matrix per line.
				      */
    void print (ostream &out) const;

				     /**
				      * Print the matrix in the usual format,
				      * i.e. as a matrix and not as a list of
				      * nonzero elements. For better
				      * readability, elements not in the matrix
				      * are displayed as empty space, while
				      * matrix elements which are explicitely
				      * set to zero are displayed as such.
				      *
				      * Each entry is printed in scientific
				      * format, with one pre-comma digit and
				      * the number of digits given by
				      * #precision# after the comma, with one
				      * space following.
				      * The precision defaults to four, which
				      * suffices for most cases. The precision
				      * and output format are {\it not}
				      * properly reset to the old values
				      * when the function exits.
				      *
				      * You should be aware that this function
				      * may produce {\bf large} amounts of
				      * output if applied to a large matrix!
				      * Be careful with it.
				      */
    void print_formatted (ostream &out,
			  const unsigned int presicion=3) const;
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotCompressed);
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotInitialized);
				     /**
				      * Exception
				      */
    DeclException2 (ExcDimensionsDontMatch,
		    int, int,
		    << "The dimensions " << arg1 << " and " << arg2
		    << " do not match properly.");
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
    DeclException1 (ExcInvalidIndex1,
		    int,
		    << "The index " << arg1 << " is not in the allowed range.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcMatrixNotSquare);
				     /**
				      * Exception
				      */
    DeclException0 (ExcDifferentSparsityPatterns);
				     /**
				      * Exception
				      */
    DeclException0 (ExcIO);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoILU);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidConstructorCall);
    
  private:
    const SparseMatrixStruct * cols;
    number* val;
    unsigned int max_len;
    bool is_ilu;
    

				     // make all other sparse matrices
				     // friends
    template <typename somenumber> friend class SparseMatrix<somenumber>;
};





/*---------------------- Inline functions -----------------------------------*/


inline
unsigned int SparseMatrixStruct::n_rows () const {
  return rows;
};



inline
unsigned int SparseMatrixStruct::n_cols () const {
  return cols;
};



inline
bool SparseMatrixStruct::is_compressed () const {
  return compressed;
};



inline
const unsigned int * SparseMatrixStruct::get_rowstart_indices () const {
  return rowstart;
};



inline
const int * SparseMatrixStruct::get_column_numbers () const {
  return colnums;
};



template <typename number>
inline
unsigned int SparseMatrix<number>::m () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->rows;
};



template <typename number>
inline
unsigned int SparseMatrix<number>::n () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->cols;
};



template <typename number>
inline
void SparseMatrix<number>::set (const unsigned int i, const unsigned int j,
				const number value) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert ((cols->operator()(i,j) != -1) || (value == 0.),
	  ExcInvalidIndex(i,j));

  const int index = cols->operator()(i,j);

  if (index >= 0) val[index] = value;
};



template <typename number>
inline
void SparseMatrix<number>::add (const unsigned int i, const unsigned int j,
				const number value) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert ((cols->operator()(i,j) != -1) || (value == 0.),
	  ExcInvalidIndex(i,j));

  const int index = cols->operator()(i,j);
  
  if (index >= 0) val[index] += value;
};





template <typename number>
inline
number SparseMatrix<number>::operator () (const unsigned int i, const unsigned int j) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  return val[cols->operator()(i,j)];
};



template <typename number>
inline
number SparseMatrix<number>::diag_element (const unsigned int i) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (i<max_len, ExcInvalidIndex1(i));
  
				   // Use that the first element in each
				   // row of a square matrix is the main
				   // diagonal
  return val[cols->rowstart[i]];
};



template <typename number>
inline
number & SparseMatrix<number>::diag_element (const unsigned int i) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (i<max_len, ExcInvalidIndex1(i));
  
				   // Use that the first element in each
				   // row of a square matrix is the main
				   // diagonal
  return val[cols->rowstart[i]];
};



template <typename number>
inline
number SparseMatrix<number>::global_entry (const unsigned int j) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return val[j];
};



template <typename number>
inline
number & SparseMatrix<number>::global_entry (const unsigned int j) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return val[j];
};



/*----------------------------   sparsematrix.h     ---------------------------*/
/* end of #ifndef __sparsematrix_H */
#endif
/*----------------------------   sparsematrix.h     ---------------------------*/


