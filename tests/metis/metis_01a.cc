// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2014 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// A test meant to identify why GridTools::partition_triangulation
// produces different output whether we're in 32- or 64-bit mode. it
// turns out that when we get METIS from PETSc and PETSc's downloaded
// version of METIS is the source of the METIS installation we use
// here, then PETSc configures METIS to use the same integer size as
// PETSc. However, then, METIS produces different output.
//
// We can not currently test for this, so the tests in this directory
// have different output depending on whether *we* (and consequently
// PETSc) use 32- or 64-bit indices, not on whether or not we use 32-
// or 64-bit METIS.

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/sparsity_pattern.h>
#include <deal.II/lac/sparsity_tools.h>

#include <metis.h>

#include "../tests.h"


void
partition(const SparsityPattern &sparsity_pattern,
          const unsigned int     n_partitions)
{
  // generate the data structures for
  // METIS. Note that this is particularly
  // simple, since METIS wants exactly our
  // compressed row storage format. we only
  // have to set up a few auxiliary arrays
  idx_t n    = static_cast<signed int>(sparsity_pattern.n_rows()),
        ncon = 1, // number of balancing constraints (should be >0)
    nparts   = static_cast<int>(n_partitions), // number of subdomains to create
    dummy; // the numbers of edges cut by the
  // resulting partition

  // use default options for METIS
  idx_t options[METIS_NOPTIONS];
  METIS_SetDefaultOptions(options);

  // one more nuisance: we have to copy our
  // own data to arrays that store signed
  // integers :-(
  std::vector<idx_t> int_rowstart(1);
  int_rowstart.reserve(sparsity_pattern.n_rows() + 1);
  std::vector<idx_t> int_colnums;
  int_colnums.reserve(sparsity_pattern.n_nonzero_elements());
  for (SparsityPattern::size_type row = 0; row < sparsity_pattern.n_rows();
       ++row)
    {
      for (SparsityPattern::iterator col = sparsity_pattern.begin(row);
           col < sparsity_pattern.end(row);
           ++col)
        int_colnums.push_back(col->column());
      int_rowstart.push_back(int_colnums.size());
    }

  std::vector<idx_t> int_partition_indices(sparsity_pattern.n_rows());

  // log the inputs to METIS
  deallog << "METIS inputs:" << std::endl;
  deallog << "IDXTYPEWIDTH=" << IDXTYPEWIDTH << std::endl;
  deallog << n << ' ' << ncon << ' ' << nparts << std::endl;
  for (unsigned int i = 0; i < int_rowstart.size(); ++i)
    deallog << int_rowstart[i] << ' ';
  deallog << std::endl;
  for (unsigned int i = 0; i < int_colnums.size(); ++i)
    deallog << int_colnums[i] << ' ';
  deallog << std::endl;
  for (unsigned int i = 0; i < METIS_NOPTIONS; ++i)
    deallog << options[i] << ' ';
  deallog << std::endl;
  deallog << sizeof(idx_t) << std::endl;


  // Make use of METIS' error code.
  int ierr;

  // Select which type of partitioning to
  // create

  ierr = METIS_PartGraphRecursive(&n,
                                  &ncon,
                                  &int_rowstart[0],
                                  &int_colnums[0],
                                  nullptr,
                                  nullptr,
                                  nullptr,
                                  &nparts,
                                  nullptr,
                                  nullptr,
                                  &options[0],
                                  &dummy,
                                  &int_partition_indices[0]);

  deallog << "METIS outputs:" << std::endl;
  deallog << dummy << std::endl;
  for (unsigned int i = 0; i < int_partition_indices.size(); ++i)
    deallog << i << ' ' << int_partition_indices[i] << std::endl;
  deallog << std::endl;
}


template <int dim>
void
test()
{
  Triangulation<dim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);

  DynamicSparsityPattern cell_connectivity;
  GridTools::get_face_connectivity_of_cells(triangulation, cell_connectivity);

  SparsityPattern sp_cell_connectivity;
  sp_cell_connectivity.copy_from(cell_connectivity);
  partition(sp_cell_connectivity, 5);
}



int
main()
{
  initlog();

  try
    {
      test<1>();
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
