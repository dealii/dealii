// $Id$
// Copyright Guido Kanschat, 1999

#ifndef __lac_sparse_vanka_H
#define __lac_sparse_vanka_H

#include <base/smartpointer.h>
#include <lac/forward-declarations.h>

#include <bvector.h>
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
 *
 * Remark: the Vanka method is a non-symmetric preconditioning method.
 * @author Guido Kanschat */
template<typename number>
class SparseVanka
{
  public:
				     /**
				      * Constructor. Gets the matrix
				      * for preconditioning and a bit
				      * vector with entries #true# for
				      * all rows to be updated. A
				      * reference to this vector will
				      * be stored, so it must persist
				      * longer than the Vanka
				      * object. The same is true for
				      * the matrix.
				      */
    SparseVanka(const SparseMatrix<number>& M,
		const bit_vector& selected);
				     /**
				      * Destructor.
				      * Delete all allocated matrices.
				      */
    ~SparseVanka();
    
				     /**
				      * Do the preconditioning. This
				      * function contains a dispatch
				      * mechanism to use the
				      * multiplicative version by
				      * default and the additive version
				      * if requested by #set_additive#.
				      */
    template<typename number2>
    void operator() (Vector<number2>& dst,
		     const Vector<number2>& src) const;

				     /**
				      * Application of the Vanka operator.
				      * This function takes two vector
				      * arguments, the residual in #src#
				      * and the resulting update vector
				      * in #dst#.
				      */
    template<typename number2>
    void forward(Vector<number2>& dst, const Vector<number2>& src) const;
				     /**
				      * Application of the transpose
				      * Vanka operator.
				      * This function takes two vector
				      * arguments, the residual in #src#
				      * and the resulting update vector
				      * in #dst#.
				      */
    template<typename number2>
    void backward(Vector<number2>& dst, const Vector<number2>& src) const;
				     /**
				      * Minimize memory consumption.
				      * Activating this option reduces
				      * memory needs of the Vanka object
				      * to nealy zero. You pay for this
				      * by a high increase of computing
				      * time, since all local matrices
				      * are built up and inverted every time.
				      */
    void conserve_memory();
    
  private:
				     /**
				      * Pointer to the matrix.
				      */
    SmartPointer<const SparseMatrix<number> > matrix;
    
				     /**
				      * Indices of Lagrange
				      * multipliers.
				      */
    const bit_vector& selected;
				     /**
				      * Conserve memory flag.
				      */
    bool conserve_mem;
				     /**
				      * Array of inverse matrices.
				      */
    mutable vector<SmartPointer<FullMatrix<float> > > inverses;
};

template<typename number>
template<typename number2>
inline
void
SparseVanka<number>::operator() (Vector<number2>& dst,
				 const Vector<number2>& src) const
{
  dst = 0.;
  forward(dst, src);
}


#endif
