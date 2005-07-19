//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__matrix_lib_templates_h
#define __deal2__matrix_lib_templates_h

#include <lac/matrix_lib.h>
#include <lac/vector.h>
#include <lac/block_vector.h>


template <typename number>
void
MeanValueFilter::filter(Vector<number>& v) const
{
  number mean = v.mean_value();

  for (unsigned int i=0;i<v.size();++i)
    v(i) -= mean;
}



template <typename number>
void
MeanValueFilter::vmult(Vector<number>& dst,
		       const Vector<number>& src) const
{
  Assert (dst.size() == src.size(),
	  ExcDimensionMismatch(dst.size(), src.size()));
  
  number mean = src.mean_value();

  for (unsigned int i=0;i<dst.size();++i)
    dst(i) = src(i) - mean;
}



template <typename number>
void
MeanValueFilter::vmult_add(Vector<number>& dst,
			   const Vector<number>& src) const
{
  Assert (dst.size() == src.size(),
	  ExcDimensionMismatch(dst.size(), src.size()));
  
  number mean = src.mean_value();
  
  for (unsigned int i=0;i<dst.size();++i)
    dst(i) += src(i) - mean;
}



template <typename number>
void
MeanValueFilter::filter(BlockVector<number>& v) const
{
  Assert (component != deal_II_numbers::invalid_unsigned_int,
	  ExcNotInitialized());
  
  for (unsigned int i=0;i<v.n_blocks();++i)
    if (i == component)
      vmult(v.block(i), v.block(i));
    else
      v.block(i) = v.block(i);
}



template <typename number>
void
MeanValueFilter::vmult(BlockVector<number>& dst,
			   const BlockVector<number>& src) const
{
  Assert (component != deal_II_numbers::invalid_unsigned_int,
	  ExcNotInitialized());
  
  Assert (dst.n_blocks() == src.n_blocks(),
	  ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));
  
  for (unsigned int i=0;i<dst.n_blocks();++i)
    if (i == component)
      vmult(dst.block(i), src.block(i));
    else
      dst.block(i) = src.block(i);
}



template <typename number>
void
MeanValueFilter::vmult_add(BlockVector<number>& dst,
			   const BlockVector<number>& src) const
{
  Assert (component != deal_II_numbers::invalid_unsigned_int,
	  ExcNotInitialized());
  
  Assert (dst.n_blocks() == src.n_blocks(),
	  ExcDimensionMismatch(dst.n_blocks(), src.n_blocks()));
  
  for (unsigned int i=0;i<dst.n_blocks();++i)
    if (i == component)
      vmult_add(dst.block(i), src.block(i));
    else
      dst.block(i).add(src.block(i));
}


//----------------------------------------------------------------------//


template <class VECTOR>
InverseMatrixRichardson<VECTOR>::InverseMatrixRichardson(
  SolverControl& c,
  VectorMemory<VECTOR>& m)
		:
		mem(m),
		solver(c,m),
		matrix(0),
		precondition(0)
{}




template <class VECTOR>
void
InverseMatrixRichardson<VECTOR>::vmult(VECTOR& dst, const VECTOR& src) const
{
  Assert (matrix != 0, ExcNotInitialized());
  Assert (precondition != 0, ExcNotInitialized());
  solver.solve(*matrix, dst, src, *precondition);
}



template <class VECTOR>
void
InverseMatrixRichardson<VECTOR>::vmult_add(VECTOR& dst, const VECTOR& src) const
{
  Assert (matrix != 0, ExcNotInitialized());
  Assert (precondition != 0, ExcNotInitialized());
  VECTOR* aux = mem.alloc();
  aux->reinit(dst);
  solver.solve(*matrix, *aux, src, *precondition);
  dst += *aux;
  mem.free(aux);
}



template <class VECTOR>
void
InverseMatrixRichardson<VECTOR>::Tvmult(VECTOR& dst, const VECTOR& src) const
{
  Assert (matrix != 0, ExcNotInitialized());
  Assert (precondition != 0, ExcNotInitialized());
  solver.Tsolve(*matrix, dst, src, *precondition);
}



template <class VECTOR>
void
InverseMatrixRichardson<VECTOR>::Tvmult_add(VECTOR& dst, const VECTOR& src) const
{
  Assert (matrix != 0, ExcNotInitialized());
  Assert (precondition != 0, ExcNotInitialized());
  VECTOR* aux = mem.alloc();
  aux->reinit(dst);
  solver.Tsolve(*matrix, *aux, src, *precondition);
  dst += *aux;
  mem.free(aux);
}




#endif
