/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2026 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE.md at
 * the top level directory of deal.II.
 *
 * ---------------------------------------------------------------------
 */

#include <gtest/gtest.h>

#include "step-80.cc"


TEST(Step80, GetDimensionAndSpacedimensionReadsParameterHandler)
{
  dealii::ParameterHandler prm;

  prm.declare_entry("dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.declare_entry("space dimension", "2", dealii::Patterns::Integer(2, 3));
  prm.set("dimension", "3");
  prm.set("space dimension", "3");

  const auto [dim, spacedim] = Step80::get_dimension_and_spacedimension(prm);

  EXPECT_EQ(dim, 3U);
  EXPECT_EQ(spacedim, 3U);
}


TEST(Step80, MPIIsInitialized)
{
  int initialized = 0;
  ASSERT_EQ(MPI_Initialized(&initialized), MPI_SUCCESS);
  EXPECT_NE(initialized, 0);

  EXPECT_EQ(dealii::Utilities::MPI::n_mpi_processes(MPI_COMM_WORLD), 1U);
}


int main(int argc, char *argv[])
{
  dealii::Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
