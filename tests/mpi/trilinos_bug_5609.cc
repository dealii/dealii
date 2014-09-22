// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2013 by the deal.II authors
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



// not really a testcase for deal.II but one that verifies that the recurring
// Trilinos bug 5609 is indeed fixed, see
//   http://software.sandia.gov/bugzilla/show_bug.cgi?id=5609


#include "../tests.h"
#include <Epetra_Map.h>
#include <Epetra_MpiComm.h>
#include <Epetra_FEVector.h>

#include <mpi.h>

#include <iostream>


void test ()
{
  int NumElements = 4;

  Epetra_MpiComm Comm (MPI_COMM_WORLD);
  Epetra_Map     Map(NumElements, 0, Comm);
  Epetra_FEVector x1(Map);
  x1.PutScalar (0);

  // let all processors set global entry 0 to 1
  const int GID = 0;
  const double value = 1;
  x1.ReplaceGlobalValues(1, &GID, &value);
  x1.GlobalAssemble (Insert);
  if (Comm.MyPID()==0)
    Assert (x1[0][0] == 1, ExcInternalError());

  // copy vector
  Epetra_FEVector x2 (x1);

  x2.PutScalar(0);

  // re-apply 1 to the vector, but only on the
  // owning processor. should be enough to set
  // the value (as non-local data in x1 should
  // be been eliminated after calling
  // GlobalAssemble).
  if (Comm.MyPID()==0)
    x2.ReplaceGlobalValues(1, &GID, &value);
  x2.GlobalAssemble (Insert);

  if (Comm.MyPID()==0)
    Assert (x1[0][0] == 1, ExcInternalError());
}


int main (int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile("output");
      deallog.attach(logfile);
      deallog << std::setprecision(4);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();

      deallog << "OK" << std::endl;
    }
  else
    test();
}
