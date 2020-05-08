// ---------------------------------------------------------------------
//
// Copyright (C) 2001 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------



#include <deal.II/lac/sparsity_pattern.h>

#include <list>
#include <set>

#include "../tests.h"

#include "../testmatrix.h"


int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  // generate usual sparsity pattern
  const unsigned int N = 15;
  SparsityPattern    sp1((N - 1) * (N - 1), (N - 1) * (N - 1), 5);
  FDMatrix(N, N).five_point_structure(sp1);
  sp1.compress();
  deallog << sp1.n_rows() << " " << sp1.n_cols() << " " << sp1.bandwidth()
          << " " << sp1.n_nonzero_elements() << std::endl;
  for (unsigned int i = 0; i < sp1.n_rows(); ++i)
    deallog << sp1.row_length(i) << std::endl;
  sp1.print_gnuplot(deallog.get_file_stream());

  // generate copy of sp1 with extra
  // off-diagonals
  SparsityPattern sp2(sp1, 10, 2);
  sp2.compress();
  deallog << sp2.n_rows() << " " << sp2.n_cols() << " " << sp2.bandwidth()
          << " " << sp2.n_nonzero_elements() << std::endl;
  for (unsigned int i = 0; i < sp2.n_rows(); ++i)
    deallog << sp2.row_length(i) << std::endl;
  sp2.print_gnuplot(deallog.get_file_stream());


  // generate copy of sp1 with extra
  // off-diagonals, add some
  // non-symmetric elements and
  // symmetrize again
  SparsityPattern sp3(sp1, (N - 1) * (N - 1), 2);
  for (unsigned int i = 0; i < (N - 1) * (N - 1); ++i)
    sp3.add(0, i);
  sp3.symmetrize();
  sp3.compress();
  deallog << sp3.n_rows() << " " << sp3.n_cols() << " " << sp3.bandwidth()
          << " " << sp3.n_nonzero_elements() << std::endl;
  for (unsigned int i = 0; i < sp3.n_rows(); ++i)
    deallog << sp3.row_length(i) << std::endl;
  sp3.print_gnuplot(deallog.get_file_stream());


  // now test the copy_from
  // function. for this copy over the
  // column indices, but in different
  // order as the order should not
  // matter to that function
  std::list<std::set<unsigned int, std::greater<unsigned int>>> sparsity;
  for (unsigned int row = 0; row < sp3.n_rows(); ++row)
    {
      sparsity.push_back(std::set<unsigned int, std::greater<unsigned int>>());
      for (SparsityPattern::const_iterator p = sp3.begin(row);
           p != sp3.end(row);
           ++p)
        sparsity.back().insert(p->column());
    };
  SparsityPattern sp4;
  sp4.copy_from((N - 1) * (N - 1),
                (N - 1) * (N - 1),
                sparsity.begin(),
                sparsity.end());

  // now check for equivalence of sp3 and sp4
  for (unsigned int row = 0; row < sp3.n_rows(); ++row)
    {
      SparsityPattern::const_iterator p3 = sp3.begin(row), p4 = sp4.begin(row);
      for (; p3 != sp3.end(row); ++p3, ++p4)
        AssertThrow(p3->column() == p4->column(), ExcInternalError());
    };


  // check the matrix_position
  // function with sparsity patterns
  // sp1 through sp4. the checked
  // function should be the inverse
  // of operator()
  //
  // check inverseness property first
  // forward, then backward
  for (unsigned int loop = 1; loop <= 4; ++loop)
    {
      const SparsityPattern &sp =
        (loop == 1 ? sp1 : (loop == 2 ? sp2 : (loop == 3 ? sp3 : sp4)));
      for (unsigned int i = 0; i < sp.n_nonzero_elements(); ++i)
        AssertThrow(sp(sp.matrix_position(i).first,
                       sp.matrix_position(i).second) == i,
                    ExcInternalError());
      for (types::global_dof_index row = 0; row < sp.n_rows(); ++row)
        for (types::global_dof_index col = 0; col < sp.n_cols(); ++col)
          if (sp(row, col) != SparsityPattern::invalid_entry)
            AssertThrow(sp.matrix_position(sp(row, col)) ==
                          std::make_pair(row, col),
                        ExcInternalError());
    };


  // check block_write/block_read by
  // dumping a sparsity pattern and
  // checking whether the
  // read-back-in pattern is the same
  std::ofstream tmp_write("sparsity_pattern.tmp");
  sp3.block_write(tmp_write);
  tmp_write.close();

  SparsityPattern sp5;

  std::ifstream tmp_read("sparsity_pattern.tmp");
  sp5.block_read(tmp_read);
  tmp_read.close();

  // delete temporary file
  std::remove("sparsity_pattern.tmp");

  // now check for equivalence of sp3 and sp5
  deallog << sp3.n_rows() - sp5.n_rows() << ' ' << sp3.n_cols() - sp5.n_cols()
          << ' ' << (sp3.is_compressed() ^ sp5.is_compressed()) << ' '
          << (sp3.is_compressed() ^ sp5.is_compressed()) << ' ' << std::endl;

  for (unsigned int row = 0; row < sp3.n_rows(); ++row)
    {
      SparsityPattern::const_iterator p3 = sp3.begin(row), p5 = sp5.begin(row);
      for (; p3 != sp3.end(row); ++p3, ++p5)
        AssertThrow(p3->column() == p5->column(), ExcInternalError());
    }
}
