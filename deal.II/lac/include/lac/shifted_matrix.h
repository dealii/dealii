//----------------------------  eigen.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  eigen.h  ---------------------------
#ifndef __deal2__shifted_matrix_h
#define __deal2__shifted_matrix_h


#include <base/smartpointer.h>

/**
 * Matrix with shifted diagonal values.
 *
 * Given a matrix p{A}, this class implements a matrix-vector product with
 * @p{A+\sigma I}, where sigma is a provided shift parameter.
 *
 * @author Guido Kanschat, 2000, 2001
 */
template<class MATRIX>
class ShiftedMatrix
{
  public:
				     /**
				      * Constructor.
				      * Provide the base matrix and a shift parameter.
				      */
    ShiftedMatrix (const MATRIX& A, const double sigma);

				     /**
				      * Set the shift parameter.
				      */
    void shift (const double sigma);

				     /**
				      * Access to the shift parameter.
				      */
    double shift () const;

				     /**
				      * Matrix-vector-product.
				      */
    template <class VECTOR>
    void vmult (VECTOR& dst, const VECTOR& src) const;

				     /**
				      * Residual.
				      */
    template <class VECTOR>
    double residual (VECTOR& dst, const VECTOR& src, const VECTOR& rhs) const;
    
  private:
				     /**
				      * Storage for base matrix.
				      */
    SmartPointer<const MATRIX> A;

				     /**
				      * Shift parameter.
				      */
    double sigma;
};
    


/**
 * Matrix with shifted diagonal values with respect to a certain scalar product.
 *
 * Given a matrix @p{A}, this class implements a matrix-vector product
 * with @p{A+\sigma M}, where sigma is a provided shift parameter and
 * @p{M} is the matrix representing the identity
 *
 * @author Guido Kanschat, 2001
 */
template<class MATRIX, class MASSMATRIX>
class ShiftedMatrixGeneralized
{
  public:
				     /**
				      * Constructor.
				      * Provide the base matrix and a shift parameter.
				      */
    ShiftedMatrixGeneralized (const MATRIX& A,
			      const MASSMATRIX& M,
			      const double sigma);

				     /**
				      * Set the shift parameter.
				      */
    void shift (const double sigma);

				     /**
				      * Access to the shift parameter.
				      */
    double shift () const;

				     /**
				      * Matrix-vector-product.
				      */
    template <class VECTOR>
    void vmult (VECTOR& dst, const VECTOR& src) const;

				     /**
				      * Residual.
				      */
    template <class VECTOR>
    double residual (VECTOR& dst, const VECTOR& src, const VECTOR& rhs) const;
    
  private:
				     /**
				      * Storage for base matrix.
				      */
    SmartPointer<const MATRIX> A;
				     /**
				      * Storage for mass matrix.
				      */
    SmartPointer<const MASSMATRIX> M;

				     /**
				      * Shift parameter.
				      */
    double sigma;
};
    


//----------------------------------------------------------------------//

template <class MATRIX>
inline
ShiftedMatrix<MATRIX>::ShiftedMatrix (const MATRIX& A, const double sigma)
		:
		A(&A), sigma(sigma)
{}



template <class MATRIX>
inline void
ShiftedMatrix<MATRIX>::shift (const double s)
{
  sigma = s;
}


template <class MATRIX>
inline double
ShiftedMatrix<MATRIX>::shift () const
{
  return sigma;
}



template <class MATRIX>
template <class VECTOR>
inline void
ShiftedMatrix<MATRIX>::vmult (VECTOR& dst, const VECTOR& src) const
{
  A.vmult(dst, src);
  dst.add(sigma, src);
}


template <class MATRIX>
template <class VECTOR>
inline double
ShiftedMatrix<MATRIX>::residual (VECTOR& dst,
				 const VECTOR& src,
				 const VECTOR& rhs) const
{
  A.vmult(dst, src);
  dst.add(sigma, src);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


//----------------------------------------------------------------------//

template <class MATRIX, class MASSMATRIX>
inline
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX>
::ShiftedMatrixGeneralized (const MATRIX& A,
			    const MASSMATRIX& M,
			    const double sigma)
		:
		A(&A), M(&M), sigma(sigma)
{}



template <class MATRIX, class MASSMATRIX>
inline void
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX>::shift (const double s)
{
  sigma = s;
}


template <class MATRIX, class MASSMATRIX>
inline double
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX>::shift () const
{
  return sigma;
}



template <class MATRIX, class MASSMATRIX>
template <class VECTOR>
inline void
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX>::vmult (VECTOR& dst,
						     const VECTOR& src) const
{
  A.vmult(dst, src);
  dst.add(sigma, src);
}


template <class MATRIX, class MASSMATRIX>
template <class VECTOR>
inline double
ShiftedMatrixGeneralized<MATRIX, MASSMATRIX>::residual (VECTOR& dst,
							const VECTOR& src,
							const VECTOR& rhs) const
{
  A.vmult(dst, src);
  dst.add(sigma, src);
  dst.sadd(-1.,1.,rhs);
  return dst.l2_norm ();
}


#endif
