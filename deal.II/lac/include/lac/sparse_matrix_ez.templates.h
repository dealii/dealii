#include <lac/sparse_matrix_ez.h>

#include <algorithm>

template <typename number>
void
SparseMatrixEZ<number>::Row::set(unsigned int column,
				 const number& value)
{
				   // Store end of vector
  const typename std::vector<Entry>::iterator end_col = end();

				   // Create Entry for inserting
  const Entry e(column, value);

				   // Find position for inserting
				   // should return first Entry with
				   // higher column number.
  typename std::vector<Entry>::iterator col = lower_bound(begin(), end_col, e);

				   // Insert Entry if column did not exist.
				   // Edit existing entry otherwise.
  if (col==end_col)
    values.push_back(e);
  else
    if (col->column == column)
      col->value = value;
    else
      values.insert(col, e);
}


template <typename number>
void
SparseMatrixEZ<number>::Row::add(unsigned int column,
				 const number& value)
{
  const typename std::vector<Entry>::const_iterator end_col = end();

				   // Create Entry for inserting
  const typename SparseMatrixEZ<number>::Entry e(column, value);
  
				   // Find position for inserting
				   // should return first Entry with
				   // higher column number.
  typename std::vector<Entry>::iterator col = lower_bound(begin(), end, e);

				   // Insert Entry if column did not exist.
				   // Edit existing entry otherwise.
  if (col==end_col)
    values.push_back(Entry(column, value));
  else
    if (col->column == column)
      col->values += value;
    else
      values.insert(col, Entry(column, value));
}


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
SparseMatrixEZ<number>::empty() const
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

  typename std::vector<Row>::const_iterator row = rows.begin();
  const typename std::vector<Row>::const_iterator end_row = rows.end();
  for (unsigned int i=0; row != end_row; ++i, ++row)
    {
      typename std::vector<Entry>::const_iterator col = row->begin();
      const typename std::vector<Entry>::const_iterator end_col = row->end();

      double s = 0.;
      for (;col != end_col; ++col)
	s += col->value * src(col->column);
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

  typename std::vector<Row>::const_iterator row = rows.begin();
  const typename std::vector<Row>::const_iterator end_row = rows.end();
  for (unsigned int i=0; row != end_row; ++i, ++row)
    {
      typename std::vector<Entry>::const_iterator col = row->begin();
      const typename std::vector<Entry>::const_iterator end_col = row->end();

      double s = 0.;
      for (;col != end_col; ++col)
	s += col->value * src(col->column);
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

  typename std::vector<Row>::const_iterator row = rows.begin();
  const typename std::vector<Row>::const_iterator end_row = rows.end();
  for (unsigned int i=0; row != end_row; ++i, ++row)
    {
      typename std::vector<Entry>::const_iterator col = row->begin();
      const typename std::vector<Entry>::const_iterator end_col = row->end();

      for (;col != end_col; ++col)
	dst(col->column) += col->value * src(i);
    }
}

template <typename number>
unsigned int
SparseMatrixEZ<number>::memory_consumption() const
{
  unsigned int result = sizeof (*this)
			+ sizeof(Row) * sizeof (rows);

  for (typename std::vector<Row>::const_iterator r = rows.begin();
       r != rows.end(); ++r)
    result += r->size() * sizeof(Entry);

  return result;
}
