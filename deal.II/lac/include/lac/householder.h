//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__householder_h
#define __deal2__householder_h


#include <base/config.h>
#include <lac/full_matrix.h>

#include <vector>

// forward declarations
template<typename number> class Vector;


/*! @addtogroup Matrix2
 *@{
 */


/**
 * QR-decomposition of a full matrix.
 *
 * @ref Instantiations: some (<tt>@<float@> @<double@></tt>)
 *
 * @author Guido Kanschat, 2005
 */
template<typename number>
class Householder : private FullMatrix<number>
{
  public:
				     /**
				      * Create an object holding the
				      * QR-decomposition of a matrix.
				      */
    template<typename number2>
    Householder (const FullMatrix<number2>&);

				     /**
				      * Solve the least-squares
				      * problem for the right hand
				      * side <tt>src</tt>.  The return
				      * value is the Euclidean norm of
				      * the approximation error.
				      */
    template<typename number2>
    double least_squares (Vector<number2> &dst,
			  Vector<number2> &src);

  private:
				     /**
				      * Storage for the diagonal
				      * elements of the orthogonal
				      * transformation.
				      */
    std::vector<number> diagonal;
};

/*@}*/

/// @if NoDoc
/*-------------------------Inline functions -------------------------------*/

// QR-transformation cf. Stoer 1 4.8.2 (p. 191)

template <typename number>
template <typename number2>
Householder<number>::Householder(const FullMatrix<number2>& M)
		:
		FullMatrix<number>(M),
		diagonal(M.n_rows())
{
//  Assert (!this->empty(), ExcEmptyMatrix());
  
				   // m > n, src.n() = m
  Assert (this->n_cols() <= this->n_rows(),
	  ExcDimensionMismatch(this->n_cols(), this->n_rows()));

  for (unsigned int j=0 ; j<n() ; ++j)
  {
    number2 sigma = 0;
    unsigned int i;
    for (i=j ; i<m() ; ++i) sigma += this->el(i,j)*this->el(i,j);
    if (std::fabs(sigma) < 1.e-15) return;
    number2 s = this->el(j,j);
    s = (s<0) ? std::sqrt(sigma) : -std::sqrt(sigma);
    number2 dj = s;

    number2 beta = 1./(s*this->el(j,j)-sigma);
    this->el(j,j) -= s;
    
    for (unsigned int k=j+1 ; k<n() ; ++k)
    {
      number2 sum = 0.;
      for (i=j ; i<m() ; ++i) sum += this->el(i,j)*this->el(i,k);
      sum *= beta;

      for (i=j ; i<m() ; ++i) this->el(i,k) += sum*this->el(i,j);
    }

    diagonal[j] = this->el(j,j);
    this->el(j,j) = dj;    
  }
}


template <typename number>
template <typename number2>
double
Householder<number>::least_squares (Vector<number2>& dst,
				    Vector<number2>& src)
{
//  Assert (!this->empty(), ExcEmptyMatrix());
  
				   // m > n, m = src.n, n = dst.n
  
  for (unsigned int j=0;j<n();++j)
    {
      number2 sum = diagonal[j]*src(j);
      for (unsigned int i=j+1 ; i<m() ; ++i)
	sum += this->el(i,j)*src(i);
// F*** what is beta???
//      sum *= beta;

      src(j) += sum*diagonal[j];
      for (unsigned int i=j+1 ; i<m() ; ++i)
	src(i) += sum*this->el(i,j);
    }
  
  backward(dst, src);
  
  number2 sum = 0.;
  for (unsigned int i=n() ; i<m() ; ++i) sum += src(i) * src(i);
  return std::sqrt(sum);
}



///@endif

#endif

