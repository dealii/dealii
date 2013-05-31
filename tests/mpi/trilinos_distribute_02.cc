//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2009, 2010, 2012, 2013 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


// check ConstraintMatrix.distribute() for a trilinos vector
//
// we do this by creating a vector where each processor has 100
// elements but no ghost elements. then we add constraints on each
// processor that constrain elements within each processor's local
// range to ones outside. these have to be added on all
// processors. then call distribute() and verify that the result is
// true.
//
// we use constraints of the form x_i = x_j with sequentially growing
// x_j's so that we can verify the correctness analytically
//
// compared to the _01 test, here the ConstraintMatrix object acts on an index
// set that only includes the locally owned vector elements, without any
// overlap. this means that the owners of a source DoF in a constraint
// do not know about the constraint if the constrained DoF is owned by
// another processor. This is not allowed, and should trigger an
// exception.


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/trilinos_vector.h>
#include <deal.II/lac/constraint_matrix.h>

#include <fstream>
#include <sstream>



void test()
{
  const unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);
  const unsigned int n_processes = Utilities::MPI::n_mpi_processes (MPI_COMM_WORLD);

  // create a vector that consists of elements indexed from 0 to n
  TrilinosWrappers::MPI::Vector vec;
  {
    IndexSet is (100*n_processes);
    is.add_range (100*myid, 100*myid+100);
    vec.reinit (is, MPI_COMM_WORLD);
  }
  Assert (vec.local_size() == 100, ExcInternalError());
  Assert (vec.local_range().first == 100*myid, ExcInternalError());
  Assert (vec.local_range().second == 100*myid+100, ExcInternalError());
  for (unsigned int i=vec.local_range().first; i<vec.local_range().second; ++i)
    vec(i) = i;
  vec.compress(VectorOperation::insert);

  // verify correctness so far
  {
    double exact_l1 = 0;
    for (unsigned int i=0; i<vec.size(); ++i)
      exact_l1 += i;
    Assert (vec.l1_norm() == exact_l1, ExcInternalError());
  }


  // create a ConstraintMatrix with a range that only includes locally owned
  // DoFs
  IndexSet locally_owned_range (vec.size());
  locally_owned_range.add_range (100*myid,
				    100*myid+100);
  ConstraintMatrix cm (locally_owned_range);

  // add constraints that constrain an element in the middle of the
  // local range of each processor against an element outside, both in
  // the ghost range before and after
  //
  // note that we tell each processor about all constraints, but most
  // of them will throw away this information since it is not for a
  // DoF inside the locally relevant range
  for (unsigned int p=0; p<n_processes; ++p)
    {
      if ((p != 0) && locally_owned_range.is_element (p*100+10))
	{
	  cm.add_line (p*100+10);
	  cm.add_entry (p*100+10,
			p*100-25,
			1);
	}

      if ((p != n_processes-1) && locally_owned_range.is_element (p*100+90))
	{
	  cm.add_line (p*100+90);
	  cm.add_entry (p*100+90,
			p*100+105,
			1);
	}
    }
  cm.close ();

  // now distribute these constraints. this should trigger an exception
  bool exception_caught = false;
  deal_II_exceptions::disable_abort_on_exception();
  try
    {
      cm.distribute (vec);
    }
  catch (...)
    {
      exception_caught = true;
    }
  if (n_processes > 1)
    Assert (exception_caught == true, ExcInternalError());

  if (myid == 0)
    deallog << "OK" << std::endl;
}


int main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  unsigned int myid = Utilities::MPI::this_mpi_process (MPI_COMM_WORLD);


  deallog.push(Utilities::int_to_string(myid));

  if (myid == 0)
    {
      std::ofstream logfile(output_file_for_mpi("trilinos_distribute_02").c_str());
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-10);

      test();
    }
  else
    test();

}
