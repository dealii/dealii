//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------

#include "../tests.h"

#include <deal.II/base/exceptions.h>
#include <deal.II/lac/vector_memory.h>
#include <deal.II/lac/vector.h>

#include <fstream>

using namespace dealii;

template<class VECTOR>
void
test_leak()
{
  GrowingVectorMemory<VECTOR> mem;
  VECTOR* v = mem.alloc();
  v->reinit(5);
}


template<class VECTOR>
void
test_stat()
{
  GrowingVectorMemory<VECTOR> mem(1, true);
  VECTOR* v1 = mem.alloc();
  VECTOR* v2 = mem.alloc();
  VECTOR* v3 = mem.alloc();
  VECTOR* v4 = mem.alloc();
  VECTOR* v5 = mem.alloc();
  VECTOR* v6 = mem.alloc();
  v1->reinit(5);
  v2->reinit(5);
  v3->reinit(5);
  v4->reinit(5);
  v5->reinit(5);
  v6->reinit(5);
  mem.free(v1);
  mem.free(v2);
  mem.free(v3);
  mem.free(v4);
  mem.free(v5);
  mem.free(v6);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
  v1 = mem.alloc();
  mem.free(v1);
}


int
main()
{
  std::ofstream logfile("vector_memory/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  test_stat<Vector<double> >();
  
  try
    {
      test_leak<Vector<double> >();
      test_leak<Vector<float> >();
    }
  catch (StandardExceptions::ExcMemoryLeak e)
    {
      e.print_exc_data(logfile);
      e.print_info(logfile);
    }
}
