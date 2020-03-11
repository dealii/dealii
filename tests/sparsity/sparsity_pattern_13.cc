// ---------------------------------------------------------------------
//
// Copyright (C) 2018 - 2019 by the deal.II authors
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



// proof of concept for DynamicSparsityPattern-like CSR object without
// special treatment of diagonals

#include "sparsity_pattern_common.h"

class SparsityPatternStandard : public SparsityPatternBase
{
public:
  using size_type = SparsityPatternBase::size_type;

  using SparsityPatternBase::reinit;

  SparsityPatternStandard()
    : SparsityPatternBase()
  {
    reinit(0, 0, 0);
  };

  virtual void
  reinit(const size_type                      m,
         const size_type                      n,
         const ArrayView<const unsigned int> &row_lengths) override
  {
    AssertDimension(row_lengths.size(), m);

    rows = m;
    cols = n;

    // delete empty matrices
    if ((m == 0) || (n == 0))
      {
        rowstart.reset();
        colnums.reset();

        max_vec_len = max_dim = rows = cols = 0;
        // if dimension is zero: ignore max_per_row
        max_row_length = 0;
        compressed     = false;

        return;
      }

    // find out how many entries we need in the @p{colnums} array. if this
    // number is larger than @p{max_vec_len}, then we will need to reallocate
    // memory
    //
    // note that the number of elements per row is bounded by the number of
    // columns
    //
    std::size_t vec_len = 0;
    for (size_type i = 0; i < m; ++i)
      vec_len += std::min(static_cast<size_type>(row_lengths[i]), n);

    // sometimes, no entries are requested in the matrix (this most often
    // happens when blocks in a block matrix are simply zero). in that case,
    // allocate exactly one element, to have a valid pointer to some memory
    if (vec_len == 0)
      {
        vec_len     = 1;
        max_vec_len = vec_len;
        colnums     = std_cxx14::make_unique<size_type[]>(max_vec_len);
      }

    max_row_length =
      (row_lengths.size() == 0 ?
         0 :
         std::min(static_cast<size_type>(
                    *std::max_element(row_lengths.begin(), row_lengths.end())),
                  n));

    // allocate memory for the rowstart values, if necessary. even though we
    // re-set the pointers again immediately after deleting their old content,
    // set them to zero in between because the allocation might fail, in which
    // case we get an exception and the destructor of this object will be called
    // -- where we look at the non-nullness of the (now invalid) pointer again
    // and try to delete the memory a second time.
    if (rows > max_dim)
      {
        max_dim  = rows;
        rowstart = std_cxx14::make_unique<std::size_t[]>(max_dim + 1);
      }

    // allocate memory for the column numbers if necessary
    if (vec_len > max_vec_len)
      {
        max_vec_len = vec_len;
        colnums     = std_cxx14::make_unique<size_type[]>(max_vec_len);
      }

    // set the rowstart array
    rowstart[0] = 0;
    for (size_type i = 1; i <= rows; ++i)
      rowstart[i] = rowstart[i - 1] +
                    std::min(static_cast<size_type>(row_lengths[i - 1]), n);
    Assert((rowstart[rows] == vec_len) ||
             ((vec_len == 1) && (rowstart[rows] == 0)),
           ExcInternalError());

    // preset the column numbers by a value indicating it is not in use
    std::fill_n(colnums.get(), vec_len, invalid_entry);

    compressed = false;
  };

  void
  copy_from(const DynamicSparsityPattern &dsp)
  {
    std::vector<unsigned int> row_lengths(dsp.n_rows());
    for (size_type i = 0; i < dsp.n_rows(); ++i)
      row_lengths[i] = dsp.row_length(i);

    reinit(dsp.n_rows(), dsp.n_cols(), row_lengths);

    if (n_rows() != 0 && n_cols() != 0)
      for (size_type row = 0; row < dsp.n_rows(); ++row)
        {
          size_type *        cols       = &colnums[rowstart[row]];
          const unsigned int row_length = dsp.row_length(row);
          for (unsigned int index = 0; index < row_length; ++index)
            {
              const size_type col = dsp.column_number(row, index);
              *cols++             = col;
            }
        }

    compressed = true;
  };
};



int
main()
{
  initlog();
  deallog << std::setprecision(3) << std::fixed;

  DynamicSparsityPattern dsp(2);
  dsp.add(0, 1);
  dsp.add(1, 0);

  SparsityPattern sp_usual;
  sp_usual.copy_from(dsp);
  deallog << "SparsityPattern:" << std::endl;
  sp_usual.print(deallog.get_file_stream());

  SparsityPatternStandard sp;
  sp.copy_from(dsp);

  deallog << "SparsityPatternStandard:" << std::endl;
  sp.print(deallog.get_file_stream());
}
