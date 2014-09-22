// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/base/index_set.h>
#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/compressed_sparsity_pattern.h>
#include <deal.II/lac/compressed_simple_sparsity_pattern.h>
#include <deal.II/lac/compressed_set_sparsity_pattern.h>
#include <deal.II/lac/chunk_sparsity_pattern.h>
#include <deal.II/lac/full_matrix.h>
#include "testmatrix.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <list>
#include <set>
#include <cstdio>


const unsigned int N = 15;

// chunk size to be used for ChunkSparsityPattern. may be overwritten in
// main()
unsigned int chunk_size = 1;


// reinitialize sparsity patterns for 5-point star
void do_reinit (SparsityPattern &sp)
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1), 5);
}


void do_reinit (ChunkSparsityPattern &sp)
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1), 5, chunk_size);
}


void do_reinit (CompressedSparsityPattern &sp)
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1));
}

void do_reinit (CompressedSimpleSparsityPattern &sp,
                const IndexSet &index_set = IndexSet())
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1), index_set);
}

void do_reinit (CompressedSetSparsityPattern &sp)
{
  sp.reinit((N-1)*(N-1), (N-1)*(N-1));
}



template <typename SP>
void build_sparsity (SP &sparsity_pattern)
{
  // generate usual 5-point sparsity pattern
  do_reinit (sparsity_pattern);
  FDMatrix(N,N).five_point_structure (sparsity_pattern);
  sparsity_pattern.compress ();

  deallog << sparsity_pattern.n_rows() << " "
          << sparsity_pattern.n_cols() << " "
          << sparsity_pattern.bandwidth() << " "
          << sparsity_pattern.n_nonzero_elements()
          << std::endl;
}


template <typename SP>
void row_length ()
{
  SP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
    deallog << sparsity_pattern.row_length(i) << std::endl;

  deallog << "OK" << std::endl;
}


template <typename SP>
void print_gnuplot ()
{
  SP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  sparsity_pattern.print_gnuplot(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



template <typename SP>
void print ()
{
  SP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  sparsity_pattern.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}


template <typename SP>
void copy_with_offdiagonals_check (SP &sp2)
{
  deallog << sp2.n_rows() << " " << sp2.n_cols() << " "
          << sp2.bandwidth() << " " << sp2.n_nonzero_elements()
          << std::endl;
  for (unsigned int i=0; i<sp2.n_rows(); ++i)
    deallog << sp2.row_length(i) << std::endl;
  sp2.print_gnuplot(deallog.get_file_stream());
}



template <typename SP>
void copy_with_offdiagonals_1 ()
{
  SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  // generate copy of sp1 with extra
  // off-diagonals
  SP sp2(sparsity_pattern, 10, 2);
  sp2.compress ();

  copy_with_offdiagonals_check (sp2);

  deallog << "OK" << std::endl;
}



template <>
void copy_with_offdiagonals_1<ChunkSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



template <>
void copy_with_offdiagonals_1<CompressedSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}

template <>
void copy_with_offdiagonals_1<CompressedSimpleSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}

template <>
void copy_with_offdiagonals_1<CompressedSetSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



template <typename SP>
void copy_with_offdiagonals_2 ()
{
  SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  // generate copy of sp1 with
  // extra off-diagonals, add some
  // non-symmetric elements and symmetrize
  // again
  SP sp3(sparsity_pattern, (N-1)*(N-1), 2);
  for (unsigned int i=0; i<(N-1)*(N-1); ++i)
    sp3.add (0,i);
  sp3.symmetrize ();
  sp3.compress ();

  copy_with_offdiagonals_check (sp3);

  deallog << "OK" << std::endl;
}


template <>
void copy_with_offdiagonals_2<ChunkSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



template <>
void copy_with_offdiagonals_2<CompressedSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



template <>
void copy_with_offdiagonals_2<CompressedSimpleSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



template <>
void copy_with_offdiagonals_2<CompressedSetSparsityPattern> ()
{
  // this sparsity pattern doesn't have this
  // function
  deallog << "OK" << std::endl;
}



void
do_copy_from (const std::list<std::set<unsigned int,std::greater<unsigned int> > > &sparsity,
              SparsityPattern &sp4)
{
  sp4.copy_from ((N-1)*(N-1), (N-1)*(N-1),
                 sparsity.begin(), sparsity.end());
}



void
do_copy_from (const std::list<std::set<unsigned int,std::greater<unsigned int> > > &sparsity,
              ChunkSparsityPattern &sp4)
{
  sp4.copy_from ((N-1)*(N-1), (N-1)*(N-1),
                 sparsity.begin(), sparsity.end(),
                 chunk_size);
}


template <typename SP>
void
do_copy_from (const CompressedSparsityPattern &sparsity,
              SP &sp4)
{
  std::list<std::set<unsigned int,std::greater<unsigned int> > > sparsity_x;
  for (unsigned int i=0; i<sparsity.n_rows(); ++i)
    {
      sparsity_x.push_back
      (std::set<unsigned int,std::greater<unsigned int> >());

      for (unsigned int j=0; j<sparsity.n_cols(); ++j)
        if (sparsity.exists(i,j))
          sparsity_x.back().insert (j);
    }

  do_copy_from (sparsity_x, sp4);
}


template <typename SP>
void
do_copy_from (const CompressedSimpleSparsityPattern &sparsity,
              SP &sp4)
{
  std::list<std::set<unsigned int,std::greater<unsigned int> > > sparsity_x;
  for (unsigned int i=0; i<sparsity.n_rows(); ++i)
    {
      sparsity_x.push_back
      (std::set<unsigned int,std::greater<unsigned int> >());

      for (unsigned int j=0; j<sparsity.n_cols(); ++j)
        if (sparsity.exists(i,j))
          sparsity_x.back().insert (j);
    }

  do_copy_from (sparsity_x, sp4);
}


template <typename SP>
void
do_copy_from (const CompressedSetSparsityPattern &sparsity,
              SP &sp4)
{
  std::list<std::set<unsigned int,std::greater<unsigned int> > > sparsity_x;
  for (unsigned int i=0; i<sparsity.n_rows(); ++i)
    {
      sparsity_x.push_back
      (std::set<unsigned int,std::greater<unsigned int> >());

      for (unsigned int j=0; j<sparsity.n_cols(); ++j)
        if (sparsity.exists(i,j))
          sparsity_x.back().insert (j);
    }

  do_copy_from (sparsity_x, sp4);
}



template <typename SP>
void
do_copy_from (const FullMatrix<double> &sparsity,
              SP &sp4)
{
  sp4.copy_from (sparsity);
}


void
do_copy_from (const FullMatrix<double> &sparsity,
              ChunkSparsityPattern &sp4)
{
  sp4.copy_from (sparsity, chunk_size);
}




template <typename SP>
void copy_from_1 ()
{
  SparsityPattern sparsity_pattern;
  build_sparsity (sparsity_pattern);

  // now test the copy_from function. for
  // this copy over the column indices, but
  // in different order as the order should
  // not matter to that function
  std::list<std::set<unsigned int,std::greater<unsigned int> > > sparsity;
  for (unsigned int row=0; row<sparsity_pattern.n_rows(); ++row)
    {
      sparsity.push_back
      (std::set<unsigned int,std::greater<unsigned int> >());
      for (SparsityPattern::iterator p = sparsity_pattern.begin(row);
           p != sparsity_pattern.end(row);
           ++p)
        sparsity.back().insert (p->column());
    }
  SP sp4;
  do_copy_from (sparsity, sp4);

  // now check for equivalence of original
  // and copy if both only store explicitly
  // added elements. otherwise, check that
  // the presence of an element in the source
  // implies an element in the copy
  if ((sparsity_pattern.stores_only_added_elements() == true)
      &&
      (sp4.stores_only_added_elements() == true))
    Assert (sparsity_pattern.n_nonzero_elements() ==
            sp4.n_nonzero_elements(),
            ExcInternalError());
  for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
    for (unsigned int j=0; j<sparsity_pattern.n_cols(); ++j)
      if ((sparsity_pattern.stores_only_added_elements() == true)
          &&
          (sp4.stores_only_added_elements() == true))
        Assert (sparsity_pattern.exists(i,j) == sp4.exists (i,j),
                ExcInternalError())
        else if (sparsity_pattern.exists(i,j))
          Assert (sp4.exists (i,j), ExcInternalError());

  deallog << "OK" << std::endl;
}



template <typename SP, typename CSP>
void copy_from_2 ()
{
  CSP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  SP sp4;
  do_copy_from (sparsity_pattern, sp4);

  // now check for equivalence of original
  // and copy if both only store explicitly
  // added elements. otherwise, check that
  // the presence of an element in the source
  // implies an element in the copy
  if ((sparsity_pattern.stores_only_added_elements() == true)
      &&
      (sp4.stores_only_added_elements() == true))
    Assert (sparsity_pattern.n_nonzero_elements() ==
            sp4.n_nonzero_elements(),
            ExcInternalError());
  for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
    for (unsigned int j=0; j<sparsity_pattern.n_cols(); ++j)
      if ((sparsity_pattern.stores_only_added_elements() == true)
          &&
          (sp4.stores_only_added_elements() == true))
        Assert (sparsity_pattern.exists(i,j) == sp4.exists (i,j),
                ExcInternalError())
        else if (sparsity_pattern.exists(i,j))
          Assert (sp4.exists (i,j), ExcInternalError());

  deallog << "OK" << std::endl;
}



template <typename SP>
void copy_from_4 ()
{
  const unsigned int M = (N-1)*(N-1);
  FullMatrix<double> sparsity_pattern(M,M);
  for (unsigned int i=0; i<M; ++i)
    for (unsigned int j=0; j<M; ++j)
      if (std::abs((int)(i-j)) == 3)
        sparsity_pattern (i,j) = 1;

  SP sp4;
  do_reinit (sp4);
  do_copy_from (sparsity_pattern, sp4);

  // now check for equivalence of original
  // and copy
  if (sp4.stores_only_added_elements() == true)
    Assert (sp4.n_nonzero_elements() ==
            static_cast<unsigned int>(sparsity_pattern.frobenius_norm() *
                                      sparsity_pattern.frobenius_norm()
                                      + 0.5),
            ExcInternalError())
    else
      Assert (sp4.n_nonzero_elements() >=
              static_cast<unsigned int>(sparsity_pattern.frobenius_norm() *
                                        sparsity_pattern.frobenius_norm()
                                        + 0.5),
              ExcInternalError());

  for (unsigned int i=0; i<M; ++i)
    for (unsigned int j=0; j<M; ++j)
      if (std::abs((int)(i-j)) == 3)
        Assert (sp4.exists (i,j)
                == true,
                ExcInternalError());

  deallog << "OK" << std::endl;
}



template <typename SP>
void matrix_position ()
{
  SP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  // check the matrix_position
  // function. the checked
  // function should be the inverse
  // of operator()
  for (unsigned int i=0; i<sparsity_pattern.n_nonzero_elements(); ++i)
    Assert (sparsity_pattern(sparsity_pattern.matrix_position(i).first,
                             sparsity_pattern.matrix_position(i).second) == i,
            ExcInternalError());
  for (types::global_dof_index row=0; row<sparsity_pattern.n_rows(); ++row)
    for (types::global_dof_index col=0; col<sparsity_pattern.n_cols(); ++col)
      if (sparsity_pattern(row,col) != SparsityPattern::invalid_entry)
        Assert (sparsity_pattern.matrix_position(sparsity_pattern(row,col)) ==
                std::make_pair(row,col),
                ExcInternalError());

  deallog << "OK" << std::endl;
}


template <>
void matrix_position<ChunkSparsityPattern> ()
{
  // this class doesn't have that function
  deallog << "OK" << std::endl;
}


template <>
void matrix_position<CompressedSparsityPattern> ()
{
  // this class doesn't have that function
  deallog << "OK" << std::endl;
}



template <>
void matrix_position<CompressedSimpleSparsityPattern> ()
{
  // this class doesn't have that function
  deallog << "OK" << std::endl;
}



template <>
void matrix_position<CompressedSetSparsityPattern> ()
{
  // this class doesn't have that function
  deallog << "OK" << std::endl;
}



template <typename SP>
void block_read_write ()
{
  SP sparsity_pattern;
  build_sparsity (sparsity_pattern);

  // check block_write/block_read by
  // dumping a sparsity pattern and
  // checking whether the
  // read-back-in pattern is the same
  std::ostringstream tmp_write;
  sparsity_pattern.block_write (tmp_write);

  SP sp5;

  std::istringstream tmp_read(tmp_write.str());
  sp5.block_read (tmp_read);

  // now check for equivalence of
  // sparsity_pattern and sp5
  deallog << sparsity_pattern.n_rows() - sp5.n_rows() << ' '
          << sparsity_pattern.n_cols() - sp5.n_cols() << ' '
          << std::endl;

  for (unsigned int i=0; i<sparsity_pattern.n_rows(); ++i)
    for (unsigned int j=0; j<sparsity_pattern.n_cols(); ++j)
      Assert (sparsity_pattern.exists(i,j) == sp5.exists (i,j),
              ExcInternalError());

  deallog << "OK" << std::endl;
}



template <>
void block_read_write<CompressedSparsityPattern> ()
{
  // not implemented for this sparsity
  // pattern
  deallog << "OK" << std::endl;
}



template <>
void block_read_write<CompressedSimpleSparsityPattern> ()
{
  // not implemented for this sparsity
  // pattern
  deallog << "OK" << std::endl;
}



template <>
void block_read_write<CompressedSetSparsityPattern> ()
{
  // not implemented for this sparsity
  // pattern
  deallog << "OK" << std::endl;
}



template <typename SP>
void test_index_set (const bool contiguous)
{
  SP sp1, sp2;

  IndexSet index_set ((N-1)*(N-1));
  index_set.add_range (5,10);

  if (!contiguous)
    for (unsigned int i=3*(N-1); i<4*(N-1); i += 3)
      index_set.add_index (i);

  do_reinit (sp1);
  do_reinit (sp2, index_set);

  FDMatrix(N,N).five_point_structure (sp1);
  FDMatrix(N,N).five_point_structure (sp2);

  sp1.compress ();
  sp2.compress ();


  for (unsigned int i=0; i<sp1.n_rows(); ++i)
    {
      deallog << sp1.row_length(i) << ' '
              << (index_set.is_element(i) ? (int)sp2.row_length(i) : -1)
              << std::endl;
      if (index_set.is_element(i))
        Assert (sp2.row_length(i) == sp1.row_length(i), ExcInternalError());

      if (index_set.is_element(i))
        {
          deallog << "      columns=";
          for (unsigned int j=0; j<sp2.row_length(i); ++j)
            {
              deallog << sp1.column_number(i,j) << ' ';
              Assert (sp1.column_number(i,j) ==
                      sp2.column_number(i,j),
                      ExcInternalError());
            }
          deallog << std::endl;
        }
    }

  deallog << "OK" << std::endl;
}
