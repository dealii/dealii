//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------
#ifndef __deal2__mg_matrix_h
#define __deal2__mg_matrix_h

#include <lac/vector.h>

/**
 * Multilevel matrix. This class implements the interface defined by
 * @ref{MGMatrixBase}, using @ref{MGLevelObject} of an arbitrary
 * matrix class.
 *
 * @author Guido Kanschat, 2002
 */
template <class MATRIX, class VECTOR>
class MGMatrix : public MGMatrixBase<VECTOR>,
  public SmartPointer<MGLevelObject<MATRIX> >
{
  public:
				     /**
				      * Constructor. The argument is
				      * handed over to the
				      * @p{SmartPointer} constructor.
				      */
    MGMatrix (MGLevelObject<MATRIX>* = 0);
    
				     /**
				      * Matrix-vector-multiplication on
				      * a certain level.
				      */
    virtual void vmult(unsigned int level, VECTOR& dst,
		       const VECTOR& src) const;
    
				   /**
				    * Adding matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const;

				   /**
				    * Transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const;

				   /**
				    * Adding transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult_add(unsigned int level, VECTOR& dst,
		     const VECTOR& src) const;    
};


/**
 * Multilevel matrix selecting from block matrices. This class
 * implements the interface defined by @ref{MGMatrixBase}.  The
 * template parameter @p{MATRIX} should be a block matrix class like
 * @ref{BlockSparseMatrix} or @p{BlockSparseMatrixEZ}. Then, this
 * class stores a pointer to a @ref{MGLevelObject} of this matrix
 * class. In each @p{vmult}, the block selected on initialization will
 * be multiplied with the vector provided.
 *
 * @author Guido Kanschat, 2002
 */
template <class MATRIX, typename number>
class MGMatrixSelect : public MGMatrixBase<Vector<number> >,
  public SmartPointer<MGLevelObject<MATRIX> >
{
  public:
				     /**
				      * Constructor. @p{row} and
				      * @p{col} are the coordinate of
				      * the selected block. The other
				      * argument is handed over to the
				      * @p{SmartPointer} constructor.
				      */
    MGMatrixSelect (const unsigned int row = 0,
		    const unsigned int col = 0,
		    MGLevelObject<MATRIX>* = 0);

				     /**
				      * Select the block for
				      * multiplication.
				      */
    void select_block (const unsigned int row,
		       const unsigned int col);
    
				     /**
				      * Matrix-vector-multiplication on
				      * a certain level.
				      */
    virtual void vmult(unsigned int level, Vector<number>& dst,
		       const Vector<number>& src) const;
    
				   /**
				    * Adding matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void vmult_add(unsigned int level, Vector<number>& dst,
			 const Vector<number>& src) const;

				   /**
				    * Transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult(unsigned int level, Vector<number>& dst,
		      const Vector<number>& src) const;

				   /**
				    * Adding transpose
				    * matrix-vector-multiplication on
				    * a certain level.
				    */
  virtual void Tvmult_add(unsigned int level, Vector<number>& dst,
			  const Vector<number>& src) const;    

  private:
				     /**
				      * Row coordinate of selected block.
				      */
    unsigned int row;
				     /**
				      * Column coordinate of selected block.
				      */
    unsigned int col;
    
};


/*----------------------------------------------------------------------*/

template <class MATRIX, class VECTOR>
MGMatrix<MATRIX, VECTOR>::MGMatrix (MGLevelObject<MATRIX>* p)
		:
		SmartPointer<MGLevelObject<MATRIX> > (p)
{}



template <class MATRIX, class VECTOR>
void
MGMatrix<MATRIX, VECTOR>::vmult (unsigned int level,
				 VECTOR& dst,
				 const VECTOR& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].vmult(dst, src);
}


template <class MATRIX, class VECTOR>
void
MGMatrix<MATRIX, VECTOR>::vmult_add (unsigned int level,
				 VECTOR& dst,
				 const VECTOR& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].vmult_add(dst, src);
}


template <class MATRIX, class VECTOR>
void
MGMatrix<MATRIX, VECTOR>::Tvmult (unsigned int level,
				 VECTOR& dst,
				 const VECTOR& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].Tvmult(dst, src);
}


template <class MATRIX, class VECTOR>
void
MGMatrix<MATRIX, VECTOR>::Tvmult_add (unsigned int level,
				 VECTOR& dst,
				 const VECTOR& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].Tvmult_add(dst, src);
}


/*----------------------------------------------------------------------*/

template <class MATRIX, typename number>
MGMatrixSelect<MATRIX, number>::
MGMatrixSelect (const unsigned int row,
		const unsigned int col,
		MGLevelObject<MATRIX>* p)
		:
		SmartPointer<MGLevelObject<MATRIX> > (p),
  row(row),
  col(col)
{}



template <class MATRIX, typename number>
void
MGMatrixSelect<MATRIX, number>::
select_block (const unsigned int brow,
	      const unsigned int bcol)
{
  row = brow;
  col = bcol;
}


template <class MATRIX, typename number>
void
MGMatrixSelect<MATRIX, number>::
vmult (unsigned int level,
       Vector<number>& dst,
       const Vector<number>& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].block(row, col).vmult(dst, src);
}


template <class MATRIX, typename number>
void
MGMatrixSelect<MATRIX, number>::
vmult_add (unsigned int level,
	   Vector<number>& dst,
	   const Vector<number>& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].block(row, col).vmult_add(dst, src);
}


template <class MATRIX, typename number>
void
MGMatrixSelect<MATRIX, number>::
Tvmult (unsigned int level,
	Vector<number>& dst,
	const Vector<number>& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].block(row, col).Tvmult(dst, src);
}


template <class MATRIX, typename number>
void
MGMatrixSelect<MATRIX, number>::
Tvmult_add (unsigned int level,
	    Vector<number>& dst,
	    const Vector<number>& src) const
{
  Assert((SmartPointer<MGLevelObject<MATRIX> >) *this != 0,
	 ExcNotInitialized());
  
  const MGLevelObject<MATRIX>& m = **this;
  m[level].block(row, col).Tvmult_add(dst, src);
}

#endif
