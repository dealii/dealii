/*----------------------------   dsmatrix.h     ---------------------------*/
/*      $Id$                 */
#ifndef __dsmatrix_H
#define __dsmatrix_H
/*----------------------------   dsmatrix.h     ---------------------------*/


// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#include <math.h>
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
  int max_dim;
  int rows, cols;
  int vec_len, max_vec_len;
  int max_row_len;
  int* rowstart;
  int* colnums;
  int compressed;

public:
  //////////
  void reinit(int m, int n, int max_per_row);
				     //////////
  dSMatrixStruct(int m, int n, int max_per_row);
    
				     //////////
  dSMatrixStruct(int n, int max_per_row);
				     //////////
  ~dSMatrixStruct();
  //////////
  void compress();
  //////////
  int operator () (int i, int j);
  //////////
  void add(int i, int j);
  //////////
  void add_matrix(int n, int* rowcols);
  //////////
  void add_matrix(int m, int n, int* rows, int* cols);
  //////////
  void add_matrix(const iVector& rowcols);
  //////////
  void add_matrix(const iVector& rows, const iVector& cols);

    void print_gnuplot (ostream &) const;
    int n_rows () const;
    int n_cols () const;
    int bandwidth () const;

    friend class ConstraintMatrix;
    friend class DoFHandler<1>;
    friend class DoFHandler<2>;

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
    DeclException1 (ExcNotEnoughSpace,
		    int,
		    << "Upon entering a new entry to row " << arg1
		    << ": there was no free entry any more.");
};




/*
CLASS
   dSMatrix
   */
class dSMatrix
{
    dSMatrixStruct * cols;
    double* val;
    int max_len;
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
    dSMatrix(dSMatrixStruct& c);
    
				     //
    ~dSMatrix();
    

				     //
    void reinit();
				     //
    void reinit (dSMatrixStruct &);

				     //
    int m() const;
				     //
    int n() const;

				     //
    void set (int i, int j, double value);
				     //
    void add (int i, int j, double value);
    
				     //
    void vmult (dVector& dst,const dVector& src) const;
				     //
    void Tvmult(dVector& dst,const dVector& src) const;
  
				     //
    double residual (dVector& dst,const dVector& x,const dVector& b);
				     //
    void Jacobi_precond(dVector& dst,const dVector& src,double om = 1.);
				     //
    void SSOR_precond(dVector& dst,const dVector& src,double om = 1.);
				     //
    void SOR_precond(dVector& dst,const dVector& src,double om = 1.);
				     //
    void SSOR(dVector& dst,double om = 1.);
				     //
    void SOR(dVector& dst,double om = 1.);
				     //
    void precondition(dVector& dst,const dVector& src) { dst=src; }

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

  friend class ConstraintMatrix;
};





/*---------------------- Inline functions -----------------------------------*/

inline
int dSMatrixStruct::n_rows () const {
  return rows;
};



inline
int dSMatrixStruct::n_cols () const {
  return cols;
};





inline
int dSMatrix::m () const
{
  return cols->rows;
};



inline
int dSMatrix::n () const
{
  return cols->cols;
};



inline
void dSMatrix::set (int i, int j, double value) {
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)] = value;
};



inline
void dSMatrix::add (int i, int j, double value) {
  Assert (cols->operator()(i,j) != -1,
	  ExcInvalidIndex(i,j));
  val[cols->operator()(i,j)]+= value;
};












/*----------------------------   dsmatrix.h     ---------------------------*/
/* end of #ifndef __dsmatrix_H */
#endif
/*----------------------------   dsmatrix.h     ---------------------------*/
