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
#ifndef __deal2__lapack_full_matrix_h
#define __deal2__lapack_full_matrix_h

#include <base/table.h>
#include <lac/lapack_support.h>

// forward declarations
template<typename number> class Vector;
template<typename number> class FullMatrix;

/*! @addtogroup Matrix1
 *@{
 */


/**
 * A variant of FullMatrix using LAPACK functions whereever possible.
 *
 * @author Guido Kanschat, 2005
 */
template <typename number>
class LAPACKFullMatrix : public TransposeTable<number>
{
  public:
				     /**
				      * Constructor. Initialize the
				      * matrix as a square matrix with
				      * dimension <tt>n</tt>.
				      *
				      * In order to avoid the implicit
				      * conversion of integers and
				      * other types to a matrix, this
				      * constructor is declared
				      * <tt>explicit</tt>.
				      *
				      * By default, no memory is
				      * allocated.
				      */
    explicit LAPACKFullMatrix (const unsigned int n = 0);
    
				     /**
				      * Constructor. Initialize the
				      * matrix as a rectangular
				      * matrix.
				      */
    LAPACKFullMatrix (const unsigned int rows,
		      const unsigned int cols);
    
				     /** 
				      * Copy constructor. This
				      * constructor does a deep copy
				      * of the matrix. Therefore, it
				      * poses a possible efficiency
				      * problem, if for example,
				      * function arguments are passed
				      * by value rather than by
				      * reference. Unfortunately, we
				      * can't mark this copy
				      * constructor <tt>explicit</tt>,
				      * since that prevents the use of
				      * this class in containers, such
				      * as <tt>std::vector</tt>. The
				      * responsibility to check
				      * performance of programs must
				      * therefore remain with the
				      * user of this class.
				      */
    LAPACKFullMatrix (const LAPACKFullMatrix&);

				     /**
				      * Assignment operator.
				      */
    LAPACKFullMatrix<number> &
    operator = (const LAPACKFullMatrix<number>&);
    
				     /**
				      * Assignment operator for a
				      * regular FullMatrix.
				      */
    template <typename number2>
    LAPACKFullMatrix<number> &
    operator = (const FullMatrix<number2>&);
    
				     /**
				      * This operator assigns a scalar
				      * to a matrix. To avoid
				      * confusion with constructors,
				      * zero is the only value allowed
				      * for <tt>d</tt>
				      */
    LAPACKFullMatrix<number> &
    operator = (const double d);

                                     /**
				      * Assignment from different
				      * matrix classes. This
				      * assignment operator uses
				      * iterators of the class
				      * MATRIX. Therefore, sparse
				      * matrices are possible sources.
				      */
    template <class MATRIX>
    void copy_from (const MATRIX&);
    
				     /**
				      * Fill rectangular block.
				      *
				      * A rectangular block of the
				      * matrix <tt>src</tt> is copied into
				      * <tt>this</tt>. The upper left
				      * corner of the block being
				      * copied is
				      * <tt>(src_offset_i,src_offset_j)</tt>.
				      * The upper left corner of the
				      * copied block is
				      * <tt>(dst_offset_i,dst_offset_j)</tt>.
				      * The size of the rectangular
				      * block being copied is the
				      * maximum size possible,
				      * determined either by the size
				      * of <tt>this</tt> or <tt>src</tt>.
				      */
    template<class MATRIX>
    void fill (const MATRIX &src,
	       const unsigned int dst_offset_i = 0,
	       const unsigned int dst_offset_j = 0,
	       const unsigned int src_offset_i = 0,
	       const unsigned int src_offset_j = 0);
    
				     /**
				      * Matrix-vector-multiplication.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>w</tt> or added
				      * to <tt>w</tt>.
				      *
				      * if (adding)
				      *  <i>w += A*v</i>
				      *
				      * if (!adding)
				      *  <i>w = A*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    void vmult (Vector<number>   &w,
		const Vector<number> &v,
		const bool            adding=false) const;
				     /**
				      * Adding Matrix-vector-multiplication.
				      *  <i>w += A*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    void vmult_add (Vector<number>       &w,
		    const Vector<number> &v) const;
    
				     /**
				      * Transpose
				      * matrix-vector-multiplication.
				      *
				      * The optional parameter
				      * <tt>adding</tt> determines, whether the
				      * result is stored in <tt>w</tt> or added
				      * to <tt>w</tt>.
				      *
				      * if (adding)
				      *  <i>w += A<sup>T</sup>*v</i>
				      *
				      * if (!adding)
				      *  <i>w = A<sup>T</sup>*v</i>
                                      *
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    void Tvmult (Vector<number>       &w,
		 const Vector<number> &v,
		 const bool            adding=false) const;

				     /**
				      * Adding transpose
				      * matrix-vector-multiplication.
				      *  <i>w += A<sup>T</sup>*v</i>
                                      *
                                      * Source and destination must
                                      * not be the same vector.
				      */
    void Tvmult_add (Vector<number>       &w,
		     const Vector<number> &v) const;

};


template <typename number>
template <class MATRIX>
void
LAPACKFullMatrix<number>::copy_from (const MATRIX& M)
{
  this->reinit (M.m(), M.n());
  const typename MATRIX::const_iterator end = M.end();
  for (typename MATRIX::const_iterator entry = M.begin();
       entry != end; ++entry)
    this->el(entry->row(), entry->column()) = entry->value();
}



template <typename number>
template <class MATRIX>
void
LAPACKFullMatrix<number>::fill (
  const MATRIX& M,
  const unsigned int dst_offset_i,
  const unsigned int dst_offset_j,
  const unsigned int src_offset_i,
  const unsigned int src_offset_j)
{
  const unsigned int endrow = src_offset_i + this->m();
  const unsigned int endcol = src_offset_j + this->n();

  for (unsigned int i=0;i<this->m();++i)
    {
      const typename MATRIX::const_iterator end = M.end(src_offset_i+i);
      for (typename MATRIX::const_iterator entry = M.begin(src_offset_i+i);
	   entry != end; ++entry)
	if (entry->column() < endcol)
	  this->el(dst_offset_i+i,
		   dst_offset_j+entry->column()-src_offset_j)
	    = entry->value();
    }
}




#endif
