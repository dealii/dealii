// $Id$

// This file is part of the DEAL Library
// DEAL is Copyright(1995) by
// Roland Becker, Guido Kanschat, Franz-Theo Suttmeier


#ifndef __lac_dsmatrix_h
#define __lac_dsmatrix_h

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
  dSMatrixStruct(int m, int n, int max_per_row) 
    : max_dim(0), max_vec_len(0), rowstart(0), colnums(0)
      {
	reinit(m,n,max_per_row);
      }
  //////////
  dSMatrixStruct(int n, int max_per_row)
    : max_dim(0), max_vec_len(0), rowstart(0), colnums(0)
      {
	reinit(n,n,max_per_row);
      }
  //////////
  ~dSMatrixStruct()
    {
      delete[] rowstart;
      delete[] colnums;
    }
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
  int n_rows () const {  return rows;   };
  int n_cols () const {  return cols;   };
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
    dSMatrix(dSMatrixStruct& c)
		    : cols(&c), val(0), max_len(0)
      {
	reinit();
      }
				     //
    ~dSMatrix()
      {
	delete[] val;
      }


				     //
    void reinit();
				     //
    void reinit (dSMatrixStruct &);

				     //
    int m() const { return cols->rows; }
				     //
    int n() const { return cols->cols; }

				     //
    void set(int i,int j,double value) { val[cols->operator()(i,j)] = value; }
				     //
    void add(int i,int j,double value) { val[cols->operator()(i,j)]+= value; }
  
				     //
    void vmult (dVector& dst,const dVector& src);
				     //
    void Tvmult(dVector& dst,const dVector& src);
  
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
};
#endif
