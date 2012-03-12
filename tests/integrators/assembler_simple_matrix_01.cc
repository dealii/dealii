//----------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 2012 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

/**
 * @file Test initialization of Assembler::MatrixSimple and
 * LocalResults
 */

#include "../tests.h"
#include <base/logstream.h>
#include <lac/full_matrix.h>
#include <lac/block_indices.h>
#include <deal.II/meshworker/local_results.h>
#include <deal.II/meshworker/simple.h>

using namespace dealii;

template <class MATRIX>
void test(MeshWorker::Assembler::MatrixSimple<MATRIX>& ass)
{
  MeshWorker::LocalResults<double> info1;
  ass.initialize_info(info1, false);
  deallog << "No faces" << std::endl;
  info1.print_debug(deallog);
  
  MeshWorker::LocalResults<double> info2;
  ass.initialize_info(info2, true);
  deallog << "With faces" << std::endl;
  info2.print_debug(deallog);
}

int main()
{
  const std::string logname = JobIdentifier::base_name(__FILE__) + std::string("/output");
  std::ofstream logfile(logname.c_str());
  deallog.attach(logfile);
  deallog.depth_console (0);
  
  MeshWorker::Assembler::MatrixSimple<FullMatrix<double> > ass1;
  deallog.push("Single block");
  test(ass1);
  BlockIndices ind;
  ind.push_back(3);
  ind.push_back(5);
  ind.push_back(1);
  ass1.initialize_local_blocks(ind);
  deallog.pop();
  deallog.push("Multiple blocks");
  test(ass1);
  ind.push_back(2);
  deallog.pop();
  deallog.push("Same blocks");
  test(ass1);
  ass1.initialize_local_blocks(ind);
  deallog.pop();
  deallog.push("More blocks");
  test(ass1);  
}
