#include <lac/sparse_matrix_ez.h>

#include <algorithm>

//----------------------------------------------------------------------//

template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ()
{
  n_columns = 0;
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(const SparseMatrixEZ<number>&)
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(unsigned int n_rows,
				       unsigned int n_cols,
				       unsigned int default_row_length,
				       unsigned int default_increment)
{
  reinit(n_rows, n_cols, default_row_length, default_increment);
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
			       unsigned int n_cols,
			       unsigned int default_row_length,
			       unsigned int default_increment)
{
  if (default_row_length == Entry::invalid_entry)
    default_row_length = 5;
  if (default_increment == Entry::invalid_entry)
    default_increment = 4;
  if (default_increment == 0)
    default_increment = 4;
  increment = default_increment;
  
  n_columns = n_cols;
  row_start.resize(n_rows+1);
  data.reserve(default_row_length * n_rows + n_rows * increment);
  data.resize(default_row_length * n_rows);

  for (unsigned int i=0;i<=n_rows;++i)
    row_start[i] = i * default_row_length;
}


template <typename number>
void
SparseMatrixEZ<number>::clear()
{
  n_columns = 0;
  row_start.resize(0);
  data.resize(0);
}


template <typename number>
bool
SparseMatrixEZ<number>::empty() const
{
  return ((n_columns == 0) && (row_start.size()==0));
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::vmult (Vector<somenumber>& dst,
			       const Vector<somenumber>& src) const
{
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  const unsigned int end_row = row_start.size() - 1;
  for (unsigned int i = 0; i < end_row;++i)
    {
      unsigned int index = row_start[i];
      unsigned int end = row_start[i+1];
      double s = 0.;
      for (;index != end && data[index].column != Entry::invalid_entry;++index)
	{
	  const Entry& entry = data[index];
	  s += entry.value * src(entry.column);
	}
      dst(i) = s;
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::Tvmult (Vector<somenumber>& dst,
				const Vector<somenumber>& src) const
{
  dst = 0.;
  Tvmult_add(dst, src);
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::vmult_add (Vector<somenumber>& dst,
				   const Vector<somenumber>& src) const
{
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  const unsigned int end_row = row_start.size() - 1;
  for (unsigned int i = 0; i < end_row;++i)
    {
      unsigned int index = row_start[i];
      unsigned int end = row_start[i+1];
      double s = 0.;
      for (;index != end && data[index].column != Entry::invalid_entry;++index)
	{
	  const Entry& entry = data[index];
	  s += entry.value * src(entry.column);
	}
      dst(i) += s;
    }
}

template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::Tvmult_add (Vector<somenumber>& dst,
				    const Vector<somenumber>& src) const
{
  Assert(n() == dst.size(), ExcDimensionMismatch(n(),dst.size()));
  Assert(m() == src.size(), ExcDimensionMismatch(m(),src.size()));

  const unsigned int end_row = row_start.size() - 1;
  for (unsigned int i = 0; i < end_row;++i)
    {
      unsigned int index = row_start[i];
      unsigned int end = row_start[i+1];
      for (;index != end && data[index].column != Entry::invalid_entry;++index)
	{
	  const Entry& entry = data[index];
	  dst(entry.column) += entry.value * src(i);
	}
    }
}

template <typename number>
unsigned int
SparseMatrixEZ<number>::memory_consumption() const
{
  unsigned int result =
    sizeof (*this)
    + sizeof(unsigned int) * row_start.capacity();
    + sizeof(SparseMatrixEZ<number>::Entry) * data.capacity();
  return result;
}
