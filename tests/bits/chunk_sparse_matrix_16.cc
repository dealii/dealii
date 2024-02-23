// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// compare matrix-vector product in chunk sparse matrix with usual vmult and
// extract_row_copy.

#include <deal.II/lac/chunk_sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include "../tests.h"


void
test(const unsigned int chunk_size)
{
  deallog << "Chunk size = " << chunk_size << std::endl;

  for (unsigned int n_cols = 4; n_cols < 7; ++n_cols)
    {
      deallog << "n_cols = " << n_cols << std::endl;
      ChunkSparsityPattern sp(5, n_cols, 3, chunk_size);
      for (unsigned int i = 0; i < 5; ++i)
        for (unsigned int j = 0; j < n_cols; ++j)
          if ((i + 2 * j + 1) % 3 == 0)
            sp.add(i, j);
      sp.compress();

      ChunkSparseMatrix<double> m(sp);

      // first set a few entries
      for (unsigned int i = 0; i < m.m(); ++i)
        for (unsigned int j = 0; j < m.n(); ++j)
          if ((i + 2 * j + 1) % 3 == 0)
            m.set(i, j, i * j * .5 + .5);

      // next perform a matrix-vector product using the entries as given by
      // the entries in extract_row_copy and compare it with the exact value
      Vector<double> src(m.n()), dst(m.m()), dst_ref(m.m());
      for (unsigned int i = 0; i < src.size(); ++i)
        src(i) = random_value<double>();
      std::vector<types::global_dof_index> indices(sp.max_entries_per_row());
      std::vector<double>                  values(sp.max_entries_per_row());
      for (unsigned int i = 0; i < m.m(); ++i)
        {
          types::global_dof_index n_entries = numbers::invalid_unsigned_int;
          m.extract_row_copy(
            i, values.size(), n_entries, &indices[0], &values[0]);
          double sum = 0;
          for (unsigned int j = 0; j < n_entries; ++j)
            sum += values[j] * src(indices[j]);
          dst(i) = sum;
        }
      m.vmult(dst_ref, src);
      dst -= dst_ref;
      deallog << "Error in matrix-vector product done via extract_row_copy: "
              << dst.linfty_norm() << std::endl;
    }
}



int
main()
{
  initlog();

  try
    {
      const unsigned int chunk_sizes[] = {1, 2, 4, 5, 7};
      for (unsigned int i = 0; i < sizeof(chunk_sizes) / sizeof(chunk_sizes[0]);
           ++i)
        test(chunk_sizes[i]);
    }
  catch (const std::exception &exc)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;

      return 1;
    }
  catch (...)
    {
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
