// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// not really a testcase for deal.II but one that verifies that the recurring
// Trilinos bug 5609 is indeed fixed, see
//   http://software.sandia.gov/bugzilla/show_bug.cgi?id=5609

#include <deal.II/base/mpi.h>

#include <Epetra_FEVector.h>
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>

#include <iostream>

#include "../tests.h"


void
test()
{
  int NumElements = 4;

  Epetra_MpiComm  Comm(MPI_COMM_WORLD);
  Epetra_Map      Map(NumElements, 0, Comm);
  Epetra_FEVector x1(Map);
  x1.PutScalar(0);

  // let all processors set global entry 0 to 1
  const int    GID   = 0;
  const double value = 1;
  x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble(Insert);
  if (Comm.MyPID() == 0)
    AssertThrow(x1[0][0] == 1, ExcInternalError());

  // copy vector
  Epetra_FEVector x2(x1);

  x2.PutScalar(0);

  // re-apply 1 to the vector, but only on the
  // owning processor. should be enough to set
  // the value (as non-local data in x1 should
  // be been eliminated after calling
  // GlobalAssemble).
  if (Comm.MyPID() == 0)
    x2.ReplaceGlobalValues(1, &GID, &value);
  x2.GlobalAssemble(Insert);

  if (Comm.MyPID() == 0)
    AssertThrow(x1[0][0] == 1, ExcInternalError());
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(
    argc, argv, testing_max_num_threads());

  unsigned int myid = Utilities::MPI::this_mpi_process(MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      initlog();
      deallog << std::setprecision(4);

      test();

      deallog << "OK" << std::endl;
    }
  else
    test();
}
