// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/lac/block_sparse_matrix_ez.h>
#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/vector.h>

#include <fstream>
#include <iomanip>

using namespace dealii;

namespace
{
  using size_type = types::global_dof_index;

  template <typename Number>
  void
  setup_blocks(BlockSparseMatrixEZ<Number>  &M,
               const std::vector<size_type> &row_block_sizes,
               const std::vector<size_type> &col_block_sizes,
               const unsigned int            max_entries_per_row)
  {
    const unsigned int n_br = row_block_sizes.size();
    const unsigned int n_bc = col_block_sizes.size();

    M.reinit(n_br, n_bc);

    for (unsigned int br = 0; br < n_br; ++br)
      for (unsigned int bc = 0; bc < n_bc; ++bc)
        M.block(br, bc).reinit(row_block_sizes[br],
                               col_block_sizes[bc],
                               max_entries_per_row);

    M.collect_sizes();
  }

  template <typename VectorType>
  void
  fill_linear(VectorType &v)
  {
    for (size_type i = 0; i < v.size(); ++i)
      v(i) = static_cast<double>(i + 1);
  }

  template <typename VectorType>
  void
  print_vector(const std::string &label, const VectorType &v)
  {
    deallog << label;
    for (size_type i = 0; i < v.size(); ++i)
      deallog << ' ' << v(i);
    deallog << std::endl;
  }

  // Fill a sparse pattern with at most two entries per row:
  // For row i, put entries at c1 and c2:
  //   c1 = i % n_cols
  //   c2 = (i+1) % n_cols   (if n_cols > 1)
  // with values:
  //   a(i,c1) = 1 + i
  //   a(i,c2) = 10 + i
  template <typename Number>
  void
  fill_two_entries_per_row(BlockSparseMatrixEZ<Number> &M)
  {
    const size_type n_rows = M.m();
    const size_type n_cols = M.n();

    for (size_type i = 0; i < n_rows; ++i)
      {
        const size_type c1 = i % n_cols;
        const Number    v1 = Number(1 + i);
        M.set(i, c1, v1);

        if (n_cols > 1)
          {
            const size_type c2 = (i + 1) % n_cols;
            const Number    v2 = Number(10 + i);
            if (c2 != c1)
              M.set(i, c2, v2);
            else
              M.add(i, c2, v2);
          }
      }
  }

  template <typename Number>
  Vector<Number>
  reference_vmult(const BlockSparseMatrixEZ<Number> &bsm_ez,
                  const Vector<Number>              &x)
  {
    const size_type n_rows = bsm_ez.m();
    const size_type n_cols = bsm_ez.n();

    AssertThrow(x.size() == n_cols, ExcInternalError());

    Vector<Number> result(n_rows);
    for (size_type i = 0; i < n_rows; ++i)
      {
        Number          sum      = Number();
        const size_type column_1 = i % n_cols;
        const Number    value_1  = Number(1 + i);
        sum += value_1 * x(column_1);

        if (n_cols > 1)
          {
            const size_type column_2 = (i + 1) % n_cols;
            const Number    value_2  = Number(10 + i);
            sum += value_2 * x(column_2);
          }

        result(i) = sum;
      }
    return result;
  }

  template <typename Number>
  Vector<Number>
  reference_Tvmult(const BlockSparseMatrixEZ<Number> &bsm_ez,
                   const Vector<Number>              &x)
  {
    const size_type n_rows = bsm_ez.m();
    const size_type n_cols = bsm_ez.n();

    AssertThrow(x.size() == n_rows, ExcInternalError());

    Vector<Number> result(n_cols);
    result = 0;

    for (size_type i = 0; i < n_rows; ++i)
      {
        const size_type column_1 = i % n_cols;
        const Number    value_1  = Number(1 + i);
        result(column_1) += value_1 * x(i);

        if (n_cols > 1)
          {
            const size_type column_2 = (i + 1) % n_cols;
            const Number    value_2  = Number(10 + i);
            result(column_2) += value_2 * x(i);
          }
      }
    return result;
  }



  void
  test_get_set_add_and_scaling()
  {
    deallog.push("get/set/add/scalar");

    BlockSparseMatrixEZ<double> bsm_ez;

    // 2x2 blocks; total size 5x5
    setup_blocks<double>(bsm_ez,
                         /*row blocks*/ {2, 3},
                         /*col blocks*/ {3, 2},
                         /*max entries/row*/ 4);

    fill_two_entries_per_row(bsm_ez);

    // Basic "get" checks (via operator()(i,j)):
    deallog << "M(0,0) " << bsm_ez(0, 0) << std::endl;
    deallog << "M(0,1) " << bsm_ez(0, 1) << std::endl;
    deallog << "M(3,3) " << bsm_ez(3, 3) << std::endl;
    deallog << "M(3,4) " << bsm_ez(3, 4) << std::endl;

    // Entry add:
    bsm_ez.add(3, 4, 5.0); // modifies an entry that is present in our pattern
    deallog << "after add M(3,4) " << bsm_ez(3, 4) << std::endl;

    // Scalar mult/div:
    bsm_ez *= 2.0;
    deallog << "after *=2 M(3,4) " << bsm_ez(3, 4) << std::endl;

    bsm_ez /= 4.0; // net scale = 1/2
    deallog << "after /=4 M(3,4) " << bsm_ez(3, 4) << std::endl;

    // Also check that scaling affects vmult consistently:
    BlockVector<double> x, result;
    x.reinit({3, 2});      // col blocks
    result.reinit({2, 3}); // row blocks
    fill_linear(x);

    bsm_ez.vmult(result, x);
    print_vector("vmult_scaled", result);

    deallog.pop();
  }

  void
  test_vmult_and_Tvmult_block_block()
  {
    deallog.push("vmult/Tvmult block-block");

    BlockSparseMatrixEZ<double> bsm_ez;
    setup_blocks<double>(bsm_ez, {2, 3}, {3, 2}, 4);
    fill_two_entries_per_row(bsm_ez);

    BlockVector<double> x, result;
    x.reinit({3, 2});      // columns
    result.reinit({2, 3}); // rows
    fill_linear(x);

    bsm_ez.vmult(result, x);

    Vector<double> x_flat(bsm_ez.n());
    for (size_type i = 0; i < x_flat.size(); ++i)
      x_flat(i) = x(i);

    const Vector<double> ref_result = reference_vmult(bsm_ez, x_flat);

    for (size_type i = 0; i < bsm_ez.m(); ++i)
      AssertThrow(result(i) == ref_result(i), ExcInternalError());

    print_vector("vmult", result);

    // Tvmult (block-block)
    BlockVector<double> xt, result_T;
    xt.reinit({2, 3});       // rows
    result_T.reinit({3, 2}); // cols
    fill_linear(xt);

    bsm_ez.Tvmult(result_T, xt);

    Vector<double> xt_flat(bsm_ez.m());
    for (size_type i = 0; i < xt_flat.size(); ++i)
      xt_flat(i) = xt(i);

    const Vector<double> ref_result_T = reference_Tvmult(bsm_ez, xt_flat);

    for (size_type j = 0; j < bsm_ez.n(); ++j)
      AssertThrow(result_T(j) == ref_result_T(j), ExcInternalError());

    print_vector("Tvmult", result_T);

    deallog.pop();
  }

  void
  test_vmult_block_nonblock_and_Tvmult_nonblock_block()
  {
    deallog.push("vmult block-Vector + Tvmult Vector-block");

    // 2x1 blocks => one block column
    BlockSparseMatrixEZ<double> bsm_ez;
    setup_blocks<double>(bsm_ez, {2, 3}, {5}, 4);
    fill_two_entries_per_row(bsm_ez);

    Vector<double> x(bsm_ez.n());
    fill_linear(x);

    BlockVector<double> result;
    result.reinit({2, 3});

    bsm_ez.vmult(result, x);

    const Vector<double> y_ref = reference_vmult(bsm_ez, x);
    for (size_type i = 0; i < bsm_ez.m(); ++i)
      AssertThrow(result(i) == y_ref(i), ExcInternalError());

    print_vector("vmult", result);

    // Tvmult(Vector&, BlockVector&) valid because: one block column
    BlockVector<double> xt;
    xt.reinit({2, 3});
    fill_linear(xt);

    Vector<double> result_T(bsm_ez.n());
    bsm_ez.Tvmult(result_T, xt);

    Vector<double> xt_flat(bsm_ez.m());
    for (size_type i = 0; i < xt_flat.size(); ++i)
      xt_flat(i) = xt(i);

    const Vector<double> ref_result_T = reference_Tvmult(bsm_ez, xt_flat);
    for (size_type j = 0; j < bsm_ez.n(); ++j)
      AssertThrow(result_T(j) == ref_result_T(j), ExcInternalError());

    print_vector("Tvmult", result_T);

    deallog.pop();
  }

  void
  test_vmult_nonblock_block_and_Tvmult_block_nonblock()
  {
    deallog.push("vmult Vector-block + Tvmult block-Vector");

    // 1x2 blocks => one block row
    BlockSparseMatrixEZ<double> bsm_ez;
    setup_blocks<double>(bsm_ez, {5}, {2, 3}, 4);
    fill_two_entries_per_row(bsm_ez);

    BlockVector<double> x;
    x.reinit({2, 3});
    fill_linear(x);

    Vector<double> result(bsm_ez.m());
    bsm_ez.vmult(result, x);

    Vector<double> x_flat(bsm_ez.n());
    for (size_type i = 0; i < x_flat.size(); ++i)
      x_flat(i) = x(i);

    const Vector<double> ref_result = reference_vmult(bsm_ez, x_flat);
    for (size_type i = 0; i < bsm_ez.m(); ++i)
      AssertThrow(result(i) == ref_result(i), ExcInternalError());

    print_vector("vmult", result);

    // Tvmult(BlockVector&, Vector&) valid because: one block row
    Vector<double> xt(bsm_ez.m());
    fill_linear(xt);

    BlockVector<double> result_T;
    result_T.reinit({2, 3});
    bsm_ez.Tvmult(result_T, xt);

    const Vector<double> ref_result_T = reference_Tvmult(bsm_ez, xt);
    for (size_type j = 0; j < bsm_ez.n(); ++j)
      AssertThrow(result_T(j) == ref_result_T(j), ExcInternalError());

    print_vector("Tvmult", result_T);

    deallog.pop();
  }

  void
  test_vmult_nonblock_nonblock_and_Tvmult_nonblock_nonblock()
  {
    deallog.push("vmult/Tvmult Vector-Vector");

    // 1x1 blocks => single block
    BlockSparseMatrixEZ<double> bsm_ez;
    setup_blocks<double>(bsm_ez, {5}, {5}, 4);
    fill_two_entries_per_row(bsm_ez);

    Vector<double> x(bsm_ez.n()), result(bsm_ez.m());
    fill_linear(x);

    bsm_ez.vmult(result, x);

    const Vector<double> y_ref = reference_vmult(bsm_ez, x);
    for (size_type i = 0; i < bsm_ez.m(); ++i)
      AssertThrow(result(i) == y_ref(i), ExcInternalError());

    print_vector("vmult", result);

    Vector<double> xt(bsm_ez.m()), result_T(bsm_ez.n());
    fill_linear(xt);

    bsm_ez.Tvmult(result_T, xt);

    const Vector<double> ref_result_T = reference_Tvmult(bsm_ez, xt);
    for (size_type j = 0; j < bsm_ez.n(); ++j)
      AssertThrow(result_T(j) == ref_result_T(j), ExcInternalError());

    print_vector("Tvmult", result_T);

    deallog.pop();
  }

  void
  test_utilities()
  {
    deallog.push("ctors/operators");

    BlockSparseMatrixEZ<double> bsm_ez_1;
    deallog << "default empty " << bsm_ez_1.empty() << std::endl;

    BlockSparseMatrixEZ<double> bsm_ez_2(2, 3);
    deallog << "ctor blocks " << bsm_ez_2.n_block_rows() << ' '
            << bsm_ez_2.n_block_cols() << std::endl;
    deallog << "ctor empty " << bsm_ez_2.empty() << std::endl;

    // Test copy constructor (only works if the rhs has empty blocks)
    BlockSparseMatrixEZ<double> bsm_ez_3(bsm_ez_2);
    deallog << "copy_ctor blocks " << bsm_ez_3.n_block_rows() << ' '
            << bsm_ez_3.n_block_cols() << std::endl;
    deallog << "copy_ctor empty " << bsm_ez_3.empty() << std::endl;

    // Test assignment operator (only works if the rhs has empty blocks)
    BlockSparseMatrixEZ<double> bsm_ez_4(2, 3);
    bsm_ez_4 = bsm_ez_3;
    deallog << "copy_assign blocks " << bsm_ez_4.n_block_rows() << ' '
            << bsm_ez_4.n_block_cols() << std::endl;
    deallog << "copy_assign empty " << bsm_ez_4.empty() << std::endl;

    BlockSparseMatrixEZ<double> bsm_ez_5;
    bsm_ez_5.reinit(1, 1);
    bsm_ez_5.collect_sizes();
    deallog << "reinit_collect_sizes m n " << bsm_ez_5.m() << ' '
            << bsm_ez_5.n() << std::endl;

    BlockSparseMatrixEZ<double> bsm_ez_6;
    setup_blocks(bsm_ez_6, {2, 3}, {3, 2}, 4);
    fill_two_entries_per_row(bsm_ez_6);

    deallog << "before eq0 M(3,4) " << bsm_ez_6(3, 4) << std::endl;
    bsm_ez_6 = 0.0;
    deallog << "after eq0 M(3,4) " << bsm_ez_6(3, 4) << std::endl;

    bsm_ez_6.clear();
    deallog << "after clear empty " << bsm_ez_6.empty() << std::endl;

    BlockSparseMatrixEZ<double> bsm_ez_7;
    setup_blocks(bsm_ez_7, {5}, {5}, 2);
    fill_two_entries_per_row(bsm_ez_7);

    std::ostringstream s0, s1;
    bsm_ez_7.print_statistics(deallog, false);
    bsm_ez_7.print_statistics(deallog, true);

    deallog.pop();
  }

} // namespace



void
test()
{
  std::ofstream logfile("output");
  deallog.attach(logfile);
  deallog << std::fixed << std::setprecision(2);

  test_get_set_add_and_scaling();
  test_utilities();
  test_vmult_and_Tvmult_block_block();
  test_vmult_block_nonblock_and_Tvmult_nonblock_block();
  test_vmult_nonblock_block_and_Tvmult_block_nonblock();
  test_vmult_nonblock_nonblock_and_Tvmult_nonblock_nonblock();
}



int
main()
{
  try
    {
      test();
    }
  catch (const std::exception &e)
    {
      std::cerr << "\n\n"
                << "----------------------------------------------------\n"
                << "Exception on processing: " << e.what() << "\n"
                << "Aborting!\n"
                << "----------------------------------------------------\n";
      return 2;
    }
  catch (...)
    {
      std::cerr << "\n\n"
                << "----------------------------------------------------\n"
                << "Unknown exception!\n"
                << "Aborting!\n"
                << "----------------------------------------------------\n";
      return 3;
    }

  return 0;
}
