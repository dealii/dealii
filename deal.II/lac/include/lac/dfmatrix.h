// $Id$

#ifndef __lac_dfmatrix_h
#define __lac_dfmatrix_h

#ifndef __stdio_h
#include <stdio.h>
#endif
#ifndef __lac_dvector_h
#include <lac/dvector.h>
#endif
#ifndef __lac_ivector_h
#include <lac/ivector.h>
#endif

/// Double precision full matrix
class dFMatrix
{
  double* val;
  int dim_range, dim_image, val_size;
  void init(int n, int m);
//  dFMatrix(const dFMatrix&);

  double& el(int i, int j)  { return val[i*dim_range+j]; }
  double el(int i, int j) const { return val[i*dim_range+j]; }
  double el(int i) const { return val[i]; }
 public:
  /// Number of rows
  int m() const { return dim_image; }
  /// Number of columns
  int n() const { return dim_range; }


  /// copy constructor. Be very careful with this constructor, since
  // it may take a hige amount of computing time for large matrices!!
  dFMatrix(const dFMatrix&);
  /// Constructor for quadratic n x n matrices
  dFMatrix(int n = 1) { init(n,n); }
  /// Constructor for rectangular matrices with m rows and n columns.
  dFMatrix(int m,int n) { init(m,n); }
  /// Destructor.
  ~dFMatrix();

  /// Reinitialization of rectangular matrix.
  void reinit(int m, int n);
  /// Reinitialization of quadratic matrix
  void reinit(int n) { reinit(n,n); }
  /// Reinitialization to the same dimensions of another matrix.
  void reinit(const dFMatrix& A) { reinit(A.m(), A.n()); }

  /// Read access to matrix elements.
  double operator() (int i, int j) const
  {
    //THROW2((i<0) || (i>=dim_image), IntError(IntError::Range,i));
    //THROW2((j<0) || (j>=dim_range), IntError(IntError::Range,j));
    return el(i,j);
  }

  /// Read-Write access to matrix elements
  double& operator() (int i, int j)
  {
    //THROW2((i<0) || (i>=dim_image), IntError(IntError::Range,i));
    //THROW2((j<0) || (j>=dim_range), IntError(IntError::Range,j));
    return el(i,j);
  }

  /// Assignment operator
  dFMatrix& operator = (const dFMatrix&);
  /** Copying matrix entries.
   * Copy the elements of matrix src into this beginning at
   * element (i,j)
   */

  void fill(const dFMatrix& src, int i = 0, int j = 0);
  void add(double s, const dFMatrix& src);
  void add_diag(double s, const dFMatrix& src);
  void Tadd(double s, const dFMatrix& src);

  ///
    /* Add different rows of a matrix
     * a(i,.) += s * a(j,.)
     */
  void add_row(int i, double s, int j);
  ///
    /**
     * Add different rows of a matrix
     * a(i,.) += s * a(j,.) + t * a(k,.)
     */
  void add_row(int i, double s, int j, double t, int k);
  ///
    /**
     * Add different columns of a matrix
     * a(.,i) += s * a(.,j)
     */
  void add_col(int i, double s, int j);
  ///
    /**
     * Add different columns of a matrix
     * a(.,i) += s * a(.,j) + t * a(.,k)
     */
  void add_col(int i, double s, int j, double t, int k);
  /// Exchange contents of rows i and j
  void swap_row(int i, int j);
  /// Exchange contents of columns i and j
  void swap_col(int i, int j);
  /// Adding a scalar value on the diagonal
  void diagadd(const double& src);

  /// Matrix-matrix-multiplication dst = this * src
  void mmult(dFMatrix& dst, const dFMatrix& src) const;
  ///
  void Tmmult(dFMatrix& dst, const dFMatrix& src) const;
  ///
    /*
     * Application of a matrix to a vector.
     *
     * [flag] adding determines if the result is copied to dst or added.
     *
     * dst (+)= this * src
     */
  void vmult(dVector& dst, const dVector& src,const int adding = 0) const;
  ///
  void gsmult(dVector& dst, const dVector& src,const iVector& gl) const;
  /// Application of the transpose matrix to a vector.
  void Tvmult(dVector& dst,const dVector& src,const int adding=0) const;
  /// 
  double residual(dVector& dst, const dVector& src, const dVector& right) const;
  /// Inversion of lower triangle
  void forward(dVector& dst, const dVector& src) const;
  /// Inversion of upper triangle
  void backward(dVector& dst, const dVector& src) const;
  /// Replace this by its inverse matrix calculated with Gau&szlig;-Jordan algorithm.
  void gauss_jordan();

  ///
  /*
   * QR - factorization of a matrix.
   * The orthogonal transformation Q is applied to the vector y and the
   * matrix. After execution of householder, the upper triangle contains
   * the resulting matrix R, the lower the incomplete factorization matrices.
   */
  void householder(dVector& y);

  /// Least - Squares - Approximation by QR-factorization.
  double least_squares(dVector& dst, dVector& src);

  /// Output of the matrix in user-defined format.
  void print(FILE* fp, const char* format = 0) const;


				     /**
				      * Comparison operator. Be careful with
				      * this thing, it may eat up huge amounts
				      * of computing time!
				      */
    bool operator == (const dFMatrix &) const;

				     /**
                                      * Computes the determinant of a matrix.
                                      * This is only implemented for one two and
                                      * three dimensions, since for higher
                                      * dimensions the numerical work explodes.
                                      * Obviously, the matrix needs to be square
                                      * for this function.
                                      */
    double determinant () const;
};
#endif
