#include <lac/sparse_matrix_ez.h>

template <typename number>
void
SparseMatrixEZ<number>::Row::set(unsigned int column,
				 const number& value)
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
void
SparseMatrixEZ<number>::Row::add(unsigned int column,
				 const number& value)
{
  Assert(false, ExcNotImplemented());
}


//----------------------------------------------------------------------//

template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ()
{
  n_columns = 0;
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(const SparseMatrixEZ&)
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(unsigned int n_rows,
				       unsigned int n_cols)
{
  reinit(n_rows, n_cols);
}


template <typename number>
SparseMatrixEZ<number>::~SparseMatrixEZ()
{}


template <typename number>
SparseMatrixEZ<number>&
SparseMatrixEZ<number>::operator= (const SparseMatrixEZ<number>&)
{
  Assert (false, ExcNotImplemented());
  return *this;
}


template <typename number>
void
SparseMatrixEZ<number>::reinit(unsigned int n_rows,
			       unsigned int n_cols)
{
  n_columns = n_cols;
  rows.resize(n_rows);
}


template <typename number>
void
SparseMatrixEZ<number>::clear()
{
  n_columns = 0;
  rows.resize(0);
}


template <typename number>
bool
SparseMatrixEZ<number>::empty()
{
  return ((n_columns == 0) && (rows.size()==0));
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::vmult (Vector<somenumber>& dst,
			       const Vector<somenumber>& src) const
{
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));
}

