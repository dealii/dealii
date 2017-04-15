// ---------------------------------------------------------------------
//
// Copyright (C) 2003 - 2016 by the deal.II authors
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


#include <deal.II/lac/la_parallel_block_vector.h>
#include "fe_tools_extrapolate_common.h"

// check FETools::extrapolate on distributed triangulations
// for LinearAlgebra::distributed::BlockVector<double>

template <int dim>
void
check (const FiniteElement<dim> &fe1,
       const FiniteElement<dim> &fe2,
       const std::string        &name)
{
  if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD)==0)
    deallog << "Checking " << name
            << " in " << dim << "d:"
            << std::endl;

  // call main function in .cc files
  check_this_dealii<dim, LinearAlgebra::distributed::BlockVector<double> > (fe1, fe2);
}

#define CHECK(EL1,deg1,EL2,deg2,dim)\
  { FE_ ## EL1<dim> fe1(deg1);   \
    FE_ ## EL2<dim> fe2(deg2);   \
    check(fe1, fe2, #EL1 #deg1 " against " #EL2 #deg2); \
  }

#define CHECK_ALL(EL1,deg1,EL2,deg2)\
  { CHECK(EL1,deg1,EL2,deg2,2); \
    CHECK(EL1,deg1,EL2,deg2,3); \
  }


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization (argc, argv, 1);

  MPILogInitAll log;
  deallog.threshold_double(1.e-10);

  CHECK_ALL(Q,1,Q,1);
  CHECK_ALL(Q,1,Q,2);
  CHECK_ALL(Q,2,Q,1);
  CHECK_ALL(Q,2,Q,2);

  CHECK_ALL(DGQ,0,DGQ,0);
  CHECK_ALL(DGQ,0,DGQ,1);
  CHECK_ALL(DGQ,1,DGQ,0);
  CHECK_ALL(DGQ,1,DGQ,1);
}

