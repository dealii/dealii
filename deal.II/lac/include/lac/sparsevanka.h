// $Id$
// Copyright Guido Kanschat, 1999

#ifndef __lac_sparsematrix_H
#define __lac_sparsematrix_H

#include <base/smartpointer.h>
#include <lac/forward-declarations.h>

#include <vector>

/**
 * Point-wise Vanka preconditioning.
 * This class does Vanka preconditioning  on a point-wise base.
 * Vanka preconditioners are used for saddle point problems. There the
 * application of Jacobi or Gauﬂ-Seidel methods is impossible, because
 * the diagonal elements are zero in the rows of the Lagrange multiplier.
 *
 * It is constructed initializing a vector of indices to the degrees of
 * freedom of the Lagrange multiplier.
 *
 * In the actual preconditioning method, these rows are traversed in
 * original order. Since this is a Gauﬂ-Seidel like procedure,
 * remember to have a good ordering in advance.
 *
 * For each row, a local system of equations is built by the degree of
 * freedom itself and all other values coupling immediately. The right
 * hand side is augmented by all further couplings.
 *
 * This local system is solved and the values are updated into the
 * destination vector.
 * @author Guido Kanschat
 */
template<typename number>
class SparseVanka
{
  public:
				     /**
				      * Constructor.
				      * Take a vector of the indices
				      * of the Lagrange multiplier as
				      * argument. A reference to this
				      * vector will be stored, so it
				      * must persist longer than the
				      * Vanka object. The same is true
				      * for the matrix.
				      */
    SparseVanka(const SparseMatrix<number>& M,
		const vector<int>& indices);
				     /**
				      * Do the preconditioning.
				      */
    template<typename number2>
    void operator() (Vector<number2>& dst,
		     const Vector<number2>& src) const;

				     /**
				      * In-place application of the
				      * method.
				      */
    template<typename number2>
    void apply(Vector<number2>& dst) const;
    
  private:
				     /**
				      * Pointer to the matrix.
				      */
    SmartPointer<const SparseMatrix<number> > matrix;
    
				     /**
				      * Indices of Lagrange
				      * multipliers.
				      */
    const vector<int>& indices;
};

template<typename number>
template<typename number2>
inline
void
SparseVanka<number>::operator() (Vector<number2>& dst,
				 const Vector<number2>& src) const
{
  dst = src;
  apply(dst);
}


#endif
