//----------------------------  sparse_matrix_ez.templates.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  sparse_matrix_ez.templates.h  ---------------------------
#ifndef __deal2__sparse_matrix_ez_templates_h
#define __deal2__sparse_matrix_ez_templates_h


#include <lac/sparse_matrix_ez.h>
#include <lac/vector.h>

#include <iostream>
#include <algorithm>
#include <cmath>

//----------------------------------------------------------------------//

template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ()
{
  n_columns = 0;
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(const SparseMatrixEZ<number> &m)
                :
                Subscriptor (m)
{
  Assert(false, ExcNotImplemented());
}


template <typename number>
SparseMatrixEZ<number>::SparseMatrixEZ(const unsigned int n_rows,
				       const unsigned int n_cols,
				       const unsigned int default_row_length,
				       const unsigned int default_increment)
{
  reinit(n_rows, n_cols, default_row_length, default_increment);
}


template <typename number>
SparseMatrixEZ<number>::~SparseMatrixEZ()
{}


template <typename number>
SparseMatrixEZ<number>&
SparseMatrixEZ<number>::operator= (const SparseMatrixEZ<number>& m)
{
  Assert (m.empty(), ExcInvalidConstructorCall());
  return *this;
}


template <typename number>
void
SparseMatrixEZ<number>::reinit(const unsigned int n_rows,
			       const unsigned int n_cols,
			       unsigned int default_row_length,
			       unsigned int default_increment)
{
  clear();

  increment = default_increment;
  
  n_columns = n_cols;
  row_info.resize(n_rows);
//TODO:[GK] allow for flexible memory reservation in later version
  data.reserve(default_row_length * n_rows + n_rows * increment);
  data.resize(default_row_length * n_rows);

  for (unsigned int i=0;i<n_rows;++i)
    row_info[i].start = i * default_row_length;
}


template <typename number>
void
SparseMatrixEZ<number>::clear()
{
  n_columns = 0;
  row_info.resize(0);
  data.resize(0);
}


template <typename number>
bool
SparseMatrixEZ<number>::empty() const
{
  return ((n_columns == 0) && (row_info.size()==0));
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::vmult (Vector<somenumber>& dst,
			       const Vector<somenumber>& src) const
{
  Assert(m() == dst.size(), ExcDimensionMismatch(m(),dst.size()));
  Assert(n() == src.size(), ExcDimensionMismatch(n(),src.size()));

  const unsigned int end_row = row_info.size();
  for (unsigned int row = 0; row < end_row; ++row)
    {
      const RowInfo& ri = row_info[row];
      typename std::vector<Entry>::const_iterator
	entry = data.begin() + ri.start;
      double s = 0.;
      for (unsigned short i=0;i<ri.length;++i,++entry)
	{
	  Assert (entry->column != Entry::invalid,
		  ExcInternalError());
	  s += entry->value * src(entry->column);
	}
      dst(row) = s;
    }
}


template <typename number>
number
SparseMatrixEZ<number>::l2_norm () const
{
  number sum = 0.;
  const_iterator start = begin();
  const_iterator final = end();

  while (start != final)
    {
      const double value = start->value();
      sum += value*value;
      ++start;
    }
  return std::sqrt(sum);
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

  const unsigned int end_row = row_info.size();
  for (unsigned int row = 0; row < end_row; ++row)
    {
      const RowInfo& ri = row_info[row];
      typename std::vector<Entry>::const_iterator
	entry = data.begin() + ri.start;
      double s = 0.;
      for (unsigned short i=0;i<ri.length;++i,++entry)
	{
	  Assert (entry->column != Entry::invalid,
		  ExcInternalError());
	  s += entry->value * src(entry->column);
	}
      dst(row) += s;
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

  const unsigned int end_row = row_info.size();
  for (unsigned int row = 0; row < end_row; ++row)
    {
      const RowInfo& ri = row_info[row];
      typename std::vector<Entry>::const_iterator
	entry = data.begin() + ri.start;
      for (unsigned short i=0;i<ri.length;++i,++entry)
	{
	  Assert (entry->column != Entry::invalid,
		  ExcInternalError());
	  dst(entry->column) += entry->value * src(row);
	}
    }
}


template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::precondition_Jacobi (Vector<somenumber>       &dst,
					     const Vector<somenumber> &src,
					     const number              om) const
{
  Assert (m() == n(), ExcNoSquareMatrix());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  somenumber              *dst_ptr = dst.begin();
  const somenumber        *src_ptr = src.begin();
  typename std::vector<RowInfo>::const_iterator ri = row_info.begin();
  const typename std::vector<RowInfo>::const_iterator end = row_info.end();
  
  for (; ri != end; ++dst_ptr, ++src_ptr, ++ri)
    {
      Assert (ri->diagonal != RowInfo::invalid_diagonal, ExcNoDiagonal());
      *dst_ptr = om * *src_ptr / data[ri->start + ri->diagonal].value;
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::precondition_SOR (Vector<somenumber>       &dst,
					  const Vector<somenumber> &src,
					  const number              om) const
{
  Assert (m() == n(), ExcNoSquareMatrix());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  somenumber       *dst_ptr = dst.begin();
  const somenumber *src_ptr = src.begin();
  typename std::vector<RowInfo>::const_iterator ri = row_info.begin();
  const typename std::vector<RowInfo>::const_iterator end = row_info.end();
  
  for (; ri != end; ++dst_ptr, ++src_ptr, ++ri)
    {
      Assert (ri->diagonal != RowInfo::invalid_diagonal, ExcNoDiagonal());
      number s = *src_ptr;
      const unsigned int end_row = ri->start + ri->diagonal;
      for (unsigned int i=ri->start;i<end_row;++i)
	s -= data[i].value * dst(data[i].column);
      
      *dst_ptr = om * s / data[ri->start + ri->diagonal].value;
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::precondition_TSOR (Vector<somenumber>       &dst,
					  const Vector<somenumber> &src,
					  const number              om) const
{
  Assert (m() == n(), ExcNoSquareMatrix());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  somenumber       *dst_ptr = dst.begin()+dst.size()-1;
  const somenumber *src_ptr = src.begin()+src.size()-1;
  typename std::vector<RowInfo>::const_reverse_iterator
    ri = row_info.rbegin();
  const typename std::vector<RowInfo>::const_reverse_iterator
    end = row_info.rend();
  
  for (; ri != end; --dst_ptr, --src_ptr, ++ri)
    {
      Assert (ri->diagonal != RowInfo::invalid_diagonal, ExcNoDiagonal());
      number s = *src_ptr;
      const unsigned int end_row = ri->start + ri->length;
      for (unsigned int i=ri->start+ri->diagonal+1;i<end_row;++i)
	s -= data[i].value * dst(data[i].column);
      
      *dst_ptr = om * s / data[ri->start + ri->diagonal].value;
    }
}



template <typename number>
template <typename somenumber>
void
SparseMatrixEZ<number>::precondition_SSOR (Vector<somenumber>       &dst,
					   const Vector<somenumber> &src,
					   const number              om) const
{
  Assert (m() == n(), ExcNoSquareMatrix());
  Assert (dst.size() == n(), ExcDimensionMismatch (dst.size(), n()));
  Assert (src.size() == n(), ExcDimensionMismatch (src.size(), n()));

  somenumber       *dst_ptr = dst.begin();
  const somenumber *src_ptr = src.begin();
  typename std::vector<RowInfo>::const_iterator ri;
  const typename std::vector<RowInfo>::const_iterator end = row_info.end();

				   // Forward
  for (ri = row_info.begin(); ri != end; ++dst_ptr, ++src_ptr, ++ri)
    {
      Assert (ri->diagonal != RowInfo::invalid_diagonal, ExcNoDiagonal());
      number s = *src_ptr;
      const unsigned int end_row = ri->start + ri->diagonal;
      for (unsigned int i=ri->start;i<end_row;++i)
	s -= om * data[i].value * dst(data[i].column);
      
      *dst_ptr = s / data[ri->start + ri->diagonal].value;
    }
				   // Diagonal
  dst_ptr = dst.begin();
  for (ri = row_info.begin(); ri != end; ++dst_ptr, ++ri)
    *dst_ptr *= (2.-om) * data[ri->start + ri->diagonal].value;

				   // Backward
  typename std::vector<RowInfo>::const_reverse_iterator rri;
   const typename std::vector<RowInfo>::const_reverse_iterator
     rend = row_info.rend();
   dst_ptr = dst.begin()+dst.size()-1;
   for (rri = row_info.rbegin(); rri != rend; --dst_ptr, ++rri)
     {
      const unsigned int end_row = rri->start + rri->length;
      for (unsigned int i=rri->start+rri->diagonal+1;i<end_row;++i)
	*dst_ptr -= om * data[i].value * dst(data[i].column);
      *dst_ptr /= data[rri->start + rri->diagonal].value;
     }
}



template <typename number>
unsigned int
SparseMatrixEZ<number>::memory_consumption() const
{
  return
    sizeof (*this)
    + sizeof(unsigned int) * row_info.capacity()
    + sizeof(typename SparseMatrixEZ<number>::Entry) * data.capacity();
}


template <typename number>
template <class STREAM>
void
SparseMatrixEZ<number>::print_statistics(STREAM& out, bool full)
{
  typename std::vector<RowInfo>::const_iterator row = row_info.begin();
  const typename std::vector<RowInfo>::const_iterator endrow = row_info.end();

				   // Add up entries actually used
  unsigned int entries_used = 0;
  unsigned int max_length = 0;
  for (; row != endrow ; ++ row)
    {
      entries_used += row->length;
      if (max_length < row->length)
	max_length = row->length;
    }
  
				   // Number of entries allocated is
				   // position of last entry used
  --row;
  unsigned int entries_alloc = row->start + row->length;

  out << "SparseMatrixEZ:used     entries:" << entries_used << std::endl
      << "SparseMatrixEZ:alloc    entries:" << entries_alloc << std::endl
      << "SparseMatrixEZ:reserved entries:" << data.capacity() << std::endl;
  
  if (full)
    {
      std::vector<unsigned int> length_used (max_length+1);
      
      for (row = row_info.begin() ; row != endrow; ++row)
	{
	  ++length_used[row->length];
	}
      for (unsigned int i=0; i< length_used.size();++i)
	if (length_used[i] != 0)
	  out << "SparseMatrixEZ:entries\t" << i
	      << "\trows\t" << length_used[i]
	      << std::endl;
    }
}


template <typename number>
void
SparseMatrixEZ<number>::block_write (std::ostream &out) const 
{
  AssertThrow (out, ExcIO());
  
                                   // first the simple objects,
                                   // bracketed in [...]
  out << '[' << row_info.size() << "]["
      << n_columns << "]["
      << data.size() << "]["
      << increment << "][";
                                   // then write out real data
  typename std::vector<RowInfo>::const_iterator r = row_info.begin();
  out.write(reinterpret_cast<const char*>(&*r),
	    sizeof(RowInfo) * row_info.size());

//   Just in case that vector entries are not stored consecutively  
//   const typename std::vector<RowInfo>::const_iterator re = row_info.end();
  
//   while (r != re)
//     {
//       out.write(reinterpret_cast<const char*>(&*r),
// 		sizeof(RowInfo));
//       ++r;
//     }

  out << "][";
  
  typename std::vector<Entry>::const_iterator d = data.begin();
  out.write(reinterpret_cast<const char*>(&*d),
	    sizeof(Entry) * data.size());

//   Just in case that vector entries are not stored consecutively  
//   const typename std::vector<Entry>::const_iterator de = data.end();

//   while (d != de)
//     {
//       out.write(reinterpret_cast<const char*>(&*d),
// 		sizeof(Entry));
//       ++d;
//     }
  
  out << ']';
  
  AssertThrow (out, ExcIO());
}


#define CHECKFOR(in,a,c) {in >> c; AssertThrow(c == a, ExcIO());}

template <typename number>
void
SparseMatrixEZ<number>::block_read (std::istream &in)
{
  AssertThrow (in, ExcIO());

  char c;
  int n;
                                   // first read in simple data
  CHECKFOR(in,'[',c);
  in >> n;
  row_info.resize(n);
  
  CHECKFOR(in,']',c);
  CHECKFOR(in,'[',c);
  in >> n_columns;
  
  CHECKFOR(in,']',c);
  CHECKFOR(in,'[',c);
  in >> n;
  data.resize(n);
  
  CHECKFOR(in,']',c);
  CHECKFOR(in,'[',c);
  in >> increment;

  CHECKFOR(in,']',c);
  CHECKFOR(in,'[',c);

                                   // then read data
  in.read(reinterpret_cast<char*>(&row_info[0]),
	  sizeof(RowInfo) * row_info.size());
  
  CHECKFOR(in,']',c);
  CHECKFOR(in,'[',c);
	    
  in.read(reinterpret_cast<char*>(&data[0]),
	  sizeof(Entry) * data.size());
  
  CHECKFOR(in,']',c);
}

#endif // __deal2__sparse_matrix_ez_templates_h
