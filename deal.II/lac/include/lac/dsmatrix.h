/*----------------------------   dsmatrix.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dsmatrix_H
#define __dsmatrix_H
/*----------------------------   dsmatrix.h     ---------------------------*/


// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_dvector_h
#include <lac/dvector.h>
#endif
#ifndef __lac_ivector_h
#include <lac/ivector.h>
#endif

#include <base/exceptions.h>



class ostream;
template <int dim> class DoFHandler;

/*
CLASS
   dSMatrixStruct
   */
class dSMatrixStruct
{
  friend class dSMatrix;
  unsigned int max_dim;
  unsigned int rows, cols;
  unsigned int vec_len, max_vec_len;
  unsigned int max_row_len;
  unsigned int* rowstart;
  int* colnums;
  bool compressed;

public:
  //////////
  void reinit (const unsigned int m, const unsigned int n,
	       const unsigned int max_per_row);
				     //////////
  dSMatrixStruct (const unsigned int m, const unsigned int n,
		  const unsigned int max_per_row);
    
				     //////////
  dSMatrixStruct (const unsigned int n, const unsigned int max_per_row);
				     //////////
  ~dSMatrixStruct ();
  //////////
  void compress ();
  //////////
  int operator() (const unsigned int i, const unsigned int j);
  //////////
  void add (const unsigned int i, const unsigned int j);
  //////////
  void add_matrix (const unsigned int n, const int* rowcols);
  //////////
  void add_matrix (const unsigned int m, const unsigned int n,
		   const int* rows, const int* cols);
  //////////
  void add_matrix (const iVector& rowcols);
  //////////
  void add_matrix (const iVector& rows, const iVector& cols);

    void print_gnuplot (ostream &) const;
    unsigned int n_rows () const;
    unsigned int n_cols () const;
    unsigned int bandwidth () const;

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
};




/*
CLASS
   dSMatrix
   */
class dSMatrix
{
    dSMatrixStruct * cols;
    double* val;
    unsigned int max_len;
  public:

				     /**
				      * Constructor; initializes the matrix to be
				      * empty, without any structure, i.e. the
				      * matrix is not usable at all. This constructor
				      * is therefore only useful for matrices which
				      * are members of a class. You have to initialize
				      * the matrix before usage with
				      * #reinit(dSMatrixStruct)#.
				      */
    dSMatrix ();
    
				     //
    dSMatrix (dSMatrixStruct& c);
    
				     //
    ~dSMatrix ();
    

				     //
    void reinit ();
				     //
    void reinit (dSMatrixStruct &);

				     /**
				      * Release all memory and return to a state
				      * just like after having called the
				      * default constructor.
				      */
    void clear ();
    
				     //
    unsigned int m () const;
				     //
    unsigned int n () const;

				     //
    void set (const unsigned int i, const unsigned int j,
	      const double value);
				     //
    void add (const unsigned int i, const unsigned int j,
	      const double value);

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

				     //
    void vmult (dVector& dst,const dVector& src) const;
				     //
    void Tvmult (dVector& dst,const dVector& src) const;
  
				     //
    double residual (dVector& dst, const dVector& x, const dVector& b);
				     //
    void Jacobi_precond (dVector& dst, const dVector& src,
			 const double om = 1.);
				     //
    void SSOR_precond (dVector& dst, const dVector& src,
		       const double om = 1.);
				     //
    void SOR_precond (dVector& dst, const dVector& src,
		      const double om = 1.);
				     //
    void SSOR (dVector& dst, const double om = 1.);
				     //
    void SOR (dVector& dst, const double om = 1.);
				     //
    void precondition (dVector& dst, const dVector& src) { dst=src; }

    const dSMatrixStruct & get_sparsity_pattern () const;
    
    void print (ostream &) const;

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
  return cols->rows;
};



inline
unsigned int dSMatrix::n () const
{
  return cols->cols;
};



inline
void dSMatrix::set (const unsigned int i, const unsigned int j,
		    const double value) {
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)] = value;
};



inline
void dSMatrix::add (const unsigned int i, const unsigned int j,
		    const double value) {
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)] += value;
};





inline
double dSMatrix::operator () (const unsigned int i, const unsigned int j) const {
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  return val[cols->operator()(i,j)];
};



inline
double dSMatrix::diag_element (const unsigned int i) const {
  Assert (m() == n(), ExcMatrixNotSquare());
  Assert (i<max_len, ExcInvalidIndex1(i));
  
				   // Use that the first element in each
				   // row of a square matrix is the main
				   // diagonal
  return val[cols->rowstart[i]];
};



inline
double dSMatrix::global_entry (const unsigned int j) const {
  return val[j];
};



inline
double & dSMatrix::global_entry (const unsigned int j) {
  return val[j];
};



/*----------------------------   dsmatrix.h     ---------------------------*/
/* end of #ifndef __dsmatrix_H */
#endif
/*----------------------------   dsmatrix.h     ---------------------------*/


