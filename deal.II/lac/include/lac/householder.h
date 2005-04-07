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
 * Whenever an object of this class is created, it copies the matrix
 * given and computes its QR-decomposition by Householder
 * algorithm. Then, the function least_squares() can be used to
 * compute the vector <i>x</i> minimizing <i>||Ax-b||</i> for a given
 * vector <i>b</i>.
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
				      * side <tt>src</tt>. The return
				      * value is the Euclidean norm of
				      * the approximation error.
				      *
				      * @arg @c dst contains the
				      * solution of the least squares
				      * problem on return.
				      *
				      * @arg @c src contains the right
				      * hand side <i>b</i> of the
				      * least squares problem. It will
				      * be changed uring the algorithm
				      * and is unusable on return.
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

  for (unsigned int j=0 ; j<this->n() ; ++j)
  {
    number2 sigma = 0;
    unsigned int i;
				     // sigma = ||v||^2
    for (i=j ; i<this->m() ; ++i)
      sigma += this->el(i,j)*this->el(i,j);
				     // We are ready if the column is
				     // empty. Are we?
    if (std::fabs(sigma) < 1.e-15) return;
    
    number2 s = (this->el(j,j) < 0) ? std::sqrt(sigma) : -std::sqrt(sigma);
				     // 
    number2 beta = std::sqrt(1./(sigma-s*this->el(j,j)));

				     // Make column j the Householder
				     // vector, store first entry in
				     // diagonal
    diagonal[j] = beta*(this->el(j,j) - s);
    this->el(j,j) = s;

    for (i=j+1 ; i<this->m() ; ++i)
      this->el(i,j) *= beta;


				     // For all subsequent columns do
				     // the Householder reflexion
    for (unsigned int k=j+1 ; k<this->n() ; ++k)
    {
      number2 sum = diagonal[j]*this->el(j,k);
      for (i=j+1 ; i<this->m() ; ++i)
	sum += this->el(i,j)*this->el(i,k);

      this->el(j,k) -= sum*this->diagonal[j];
      for (i=j+1 ; i<this->m() ; ++i)
	this->el(i,k) -= sum*this->el(i,j);
    }    
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

				   // Multiply Q_n ... Q_2 Q_1 src
				   // Where Q_i = I-v_i v_i^T
  for (unsigned int j=0;j<this->n();++j)
    {
				       // sum = v_i^T src
      number2 sum = diagonal[j]*src(j);
      for (unsigned int i=j+1 ; i<this->m() ; ++i)
	sum += this->el(i,j)*src(i);
				       // src -= v * sum
      src(j) -= sum*diagonal[j];
      for (unsigned int i=j+1 ; i<this->m() ; ++i)
	src(i) -= sum*this->el(i,j);
    }
  
  backward(dst, src);
  
  number2 sum = 0.;
  for (unsigned int i=this->n() ; i<this->m() ; ++i) sum += src(i) * src(i);
  return std::sqrt(sum);
}



///@endif

#endif

