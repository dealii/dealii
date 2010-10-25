//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2005, 2006 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <lac/matrix_lib.templates.h>
#include <lac/sparse_matrix.h>

DEAL_II_NAMESPACE_OPEN

MeanValueFilter::MeanValueFilter(unsigned int component)
		:
		component(component)
{}


template<typename number, typename vnumber>
ProductSparseMatrix<number, vnumber>::ProductSparseMatrix(
  const MatrixType& mat1,
  const MatrixType& mat2,
  VectorMemory<VectorType>& mem)
		:
		m1(&mat1, typeid(*this).name()),
		m2(&mat2, typeid(*this).name()),
		mem(&mem, typeid(*this).name())
{
  Assert(mat1.n() == mat2.m(), ExcDimensionMismatch(mat1.n(),mat2.m()));
}


template<typename number, typename vnumber>
ProductSparseMatrix<number, vnumber>::ProductSparseMatrix()
		:
		m1(0, typeid(*this).name()),
		m2(0, typeid(*this).name()),
		mem(0, typeid(*this).name())
{}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::initialize(
  const MatrixType& mat1,
  const MatrixType& mat2,
  VectorMemory<VectorType>& memory)
{
  Assert(mat1.n() == mat2.m(), ExcDimensionMismatch(mat1.n(),mat2.m()));
  mem = &memory;
  m1 = &mat1;
  m2 = &mat2;
}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::clear()
{
  m1 = 0;
  m2 = 0;
}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::vmult (VectorType& dst, const VectorType& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Assert(m1 != 0, ExcNotInitialized());
  Assert(m2 != 0, ExcNotInitialized());
  
  VectorType* v = mem->alloc();
  v->reinit(m1->n());
  m2->vmult (*v, src);
  m1->vmult (dst, *v);
  mem->free(v);
}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::vmult_add (VectorType& dst, const VectorType& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Assert(m1 != 0, ExcNotInitialized());
  Assert(m2 != 0, ExcNotInitialized());
  
  VectorType* v = mem->alloc();
  v->reinit(m1->n());
  m2->vmult (*v, src);
  m1->vmult_add (dst, *v);
  mem->free(v);
}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::Tvmult (VectorType& dst, const VectorType& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Assert(m1 != 0, ExcNotInitialized());
  Assert(m2 != 0, ExcNotInitialized());
  
  VectorType* v = mem->alloc();
  v->reinit(m1->n());
  m1->Tvmult (*v, src);
  m2->Tvmult (dst, *v);
  mem->free(v);
}


template<typename number, typename vnumber>
void
ProductSparseMatrix<number, vnumber>::Tvmult_add (VectorType& dst, const VectorType& src) const
{
  Assert(mem != 0, ExcNotInitialized());
  Assert(m1 != 0, ExcNotInitialized());
  Assert(m2 != 0, ExcNotInitialized());
  
  VectorType* v = mem->alloc();
  v->reinit(m1->n());
  m1->Tvmult (*v, src);
  m2->Tvmult_add (dst, *v);
  mem->free(v);
}


template<typename number, typename vnumber>
const void*
ProductSparseMatrix<number, vnumber>::get () const
{
  return &*m1;
}


template class ProductSparseMatrix<double, double>;
template class ProductSparseMatrix<double, float>;
template class ProductSparseMatrix<float, double>;
template class ProductSparseMatrix<float, float>;

template void MeanValueFilter::filter(Vector<float>&) const;
template void MeanValueFilter::filter(Vector<double>&) const;
template void MeanValueFilter::filter(BlockVector<float>&) const;
template void MeanValueFilter::filter(BlockVector<double>&) const;
template void MeanValueFilter::vmult(Vector<float>&,
				     const Vector<float>&) const;
template void MeanValueFilter::vmult(Vector<double>&,
				     const Vector<double>&) const;
template void MeanValueFilter::vmult(BlockVector<float>&,
				     const BlockVector<float>&) const;
template void MeanValueFilter::vmult(BlockVector<double>&,
				     const BlockVector<double>&) const;

template void MeanValueFilter::vmult_add(Vector<float>&,
					 const Vector<float>&) const;
template void MeanValueFilter::vmult_add(Vector<double>&,
					 const Vector<double>&) const;
template void MeanValueFilter::vmult_add(BlockVector<float>&,
					 const BlockVector<float>&) const;
template void MeanValueFilter::vmult_add(BlockVector<double>&,
					 const BlockVector<double>&) const;

template class InverseMatrixRichardson<Vector<float> >;
template class InverseMatrixRichardson<Vector<double> >;
template class InverseMatrixRichardson<BlockVector<float> >;
template class InverseMatrixRichardson<BlockVector<double> >;

DEAL_II_NAMESPACE_CLOSE
