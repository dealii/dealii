/*----------------------------   dsmatrix.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dsmatrix_H
#define __dsmatrix_H
/*----------------------------   dsmatrix.h     ---------------------------*/


// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier
// Revised by Wolfgang Bangerth


#include <base/exceptions.h>


//forward declarations
class dVector;
class iVector;
class ostream;



/*
CLASS
   dSMatrixStruct

   @author Original version by Roland Becker, Guido Kanschat, Franz-Theo Suttmeier; lots of enhancements, reorganisation and documentation by Wolfgang Bangerth
   */
class dSMatrixStruct
{
  private:
				     /**
				      * Copy constructor, made private in order to
				      * prevent copying such an object which does
				      * not make much sense because you can use
				      * a structure like this for more than one
				      * matrix.
				      *
				      * Because it is not needed, this function
				      * is not implemented.
				      */
    dSMatrixStruct (const dSMatrixStruct &);
    
  public:
				     /**
				      * Initialize the matrix empty, i.e. with
				      * no memory allocated. This is useful if
				      * you want such objects as member
				      * variables in other classes. You can make
				      * the structure usable by calling the
				      * #reinit# function.
				      */
    dSMatrixStruct ();
    
				     /**
				      * Initialize a rectangular matrix with
				      * #m# rows and #n# columns,
				      * with at most #max_per_row#
				      * nonzero entries per row.
				      */
    dSMatrixStruct (const unsigned int m,
		    const unsigned int n,
		    const unsigned int max_per_row);
    
				     /**
				      * Initialize a square matrix of dimension
				      * #n# with at most #max_per_row#
				      * nonzero entries per row.
				      */
    dSMatrixStruct (const unsigned int n,
		    const unsigned int max_per_row);

				     /**
				      * Destructor.
				      */
    ~dSMatrixStruct ();
    
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
				      * #dSMatrix# objects require the
				      * #dSMatrixStruct# objects they are
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
				      * #dSMatrix#. It shall only be
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
    
    friend class dSMatrix;
};




/*
CLASS
   dSMatrix

   @author Original version by Roland Becker, Guido Kanschat, Franz-Theo Suttmeier; lots of enhancements, reorganisation and documentation by Wolfgang Bangerth 1998
   */
class dSMatrix
{
  public:

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
				      * #reinit(dSMatrixStruct)#.
				      */
    dSMatrix ();
    
				     /**
				      * Constructor. Takes the given matrix
				      * sparisty structure to represent the
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
    dSMatrix (const dSMatrixStruct &sparsity);
    
				     /**
				      * Destructor. Free all memory, but do not
				      * release the memory of the sparsity
				      * structure.
				      */
    virtual ~dSMatrix ();
    

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
    virtual void reinit (const dSMatrixStruct &sparsity);

				     /**
				      * Release all memory and return to a state
				      * just like after having called the
				      * default constructor. It also forgets the
				      * sparsity pattern it was previously tied
				      * to.
				      */
    virtual void clear ();
    
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
				      * not exist.
				      */
    void set (const unsigned int i, const unsigned int j,
	      const double value);
    
				     /**
				      * Add #value# to the element #(i,j)#.
				      * Throws an error if the entry does
				      * not exist.
				      */
    void add (const unsigned int i, const unsigned int j,
	      const double value);

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
				      * The function returns a reference to
				      * #this#.
				      */
    dSMatrix & copy_from (const dSMatrix &);

				     /**
				      * Add #matrix# scaled by #factor# to this
				      * matrix. The function throws an error
				      * if the sparsity patterns of the two
				      * involved matrices do not point to the
				      * same object, since in this case the
				      * operation is cheaper.
				      */
    void add_scaled (const double factor, const dSMatrix &matrix);
    
				     /**
				      * Return the value of the entry (i,j).
				      * This may be an expensive operation
				      * and you should always take care
				      * where to call this function.
				      * In order to avoid abuse, this function
				      * throws an exception if the wanted
				      * element does not exist in the matrix.
				      */
    double operator () (const unsigned int i, const unsigned int j) const;

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
    double diag_element (const unsigned int i) const;

    				     /**
				      * This is kind of an expert mode: get
				      * access to the #i#th element of this
				      * matrix. The elements are stored in
				      * a consecutive way, refer to the
				      * #dSMatrixStruct# class for more details.
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
    double global_entry (const unsigned int i) const;

				     /**
				      * Same as above, but with write access.
				      * You certainly know what you do?
				      */
    double & global_entry (const unsigned int i);

				     /**
				      * Matrix-vector multiplication: let
				      * #dst = M*src# with #M# being this matrix.
				      */
    void vmult (dVector& dst, const dVector& src) const;
    
				     /**
				      * Matrix-vector multiplication: let
				      * #dst = M^T*src# with #M# being this
				      * matrix. This function does the same as
				      * #vmult# but takes the transposed matrix.
				      */
    void Tvmult (dVector& dst, const dVector& src) const;
  

				     /**
				      * Return the norm of the vector #v# with
				      * respect to the norm induced by this
				      * matrix, i.e. $\left<v,Mv\right>$. This
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
    double matrix_norm (const dVector &v) const;
    
				     //
    double residual (dVector& dst, const dVector& x,
		     const dVector& b) const;
				     //
    void precondition_Jacobi (dVector& dst, const dVector& src,
			      const double om = 1.) const;
				     //
    void precondition_SSOR (dVector& dst, const dVector& src,
			    const double om = 1.) const;
				     //
    void precondition_SOR (dVector& dst, const dVector& src,
			   const double om = 1.) const;
				     //
    void SSOR (dVector& dst, const double om = 1.) const;
				     //
    void SOR (dVector& dst, const double om = 1.) const;
				     //
    void precondition (dVector& dst, const dVector& src) const;

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
    const dSMatrixStruct & get_sparsity_pattern () const;

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
    
  private:
    const dSMatrixStruct * cols;
    double* val;
    unsigned int max_len;
};





/*---------------------- Inline functions -----------------------------------*/

inline
unsigned int dSMatrixStruct::n_rows () const {
  return rows;
};



inline
unsigned int dSMatrixStruct::n_cols () const {
  return cols;
};



inline
bool dSMatrixStruct::is_compressed () const {
  return compressed;
};



inline
const unsigned int * dSMatrixStruct::get_rowstart_indices () const {
  return rowstart;
};



inline
const int * dSMatrixStruct::get_column_numbers () const {
  return colnums;
};



inline
unsigned int dSMatrix::m () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->rows;
};



inline
unsigned int dSMatrix::n () const
{
  Assert (cols != 0, ExcMatrixNotInitialized());
  return cols->cols;
};



inline
void dSMatrix::set (const unsigned int i, const unsigned int j,
		    const double value) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)] = value;
};



inline
void dSMatrix::add (const unsigned int i, const unsigned int j,
		    const double value) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)] += value;
};





inline
double dSMatrix::operator () (const unsigned int i, const unsigned int j) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  return val[cols->operator()(i,j)];
};



inline
double dSMatrix::diag_element (const unsigned int i) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (i<max_len, ExcInvalidIndex1(i));
  
				   // Use that the first element in each
				   // row of a square matrix is the main
				   // diagonal
  return val[cols->rowstart[i]];
};



inline
double dSMatrix::global_entry (const unsigned int j) const {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return val[j];
};



inline
double & dSMatrix::global_entry (const unsigned int j) {
  Assert (cols != 0, ExcMatrixNotInitialized());
  return val[j];
};



/*----------------------------   dsmatrix.h     ---------------------------*/
/* end of #ifndef __dsmatrix_H */
#endif
/*----------------------------   dsmatrix.h     ---------------------------*/


