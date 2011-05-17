//----------------------------------------------------------------------
//    $Id: block_sparsity_pattern_01.cc 15661 2008-01-24 04:59:09Z kanschat $
//    Version: $Name$ 
//
//    Copyright (C) 2006, 2008 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------


// Test BlockCompressedSimpleSparsityPattern with IndexSets


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main()
{
  std::ofstream logfile("block_compressed_simple_sparsity_pattern_indexset_01/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  IndexSet a(5);
  IndexSet b(3);
  a.add_index(0);
  a.add_index(3);
  b.add_index(0);
  std::vector<IndexSet> part;
  part.push_back(a);
  part.push_back(b);
  part.push_back(a);
  BlockCompressedSimpleSparsityPattern csp(part);
  
  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols() << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1,0).n_rows() << "x" << csp.block(1,0).n_cols() << std::endl;

  part.pop_back();
  csp.reinit(part); //also check the reinit variant.

  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols() << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1,0).n_rows() << "x" << csp.block(1,0).n_cols() << std::endl;

  for (int i=0;i<13;++i)
    {
      if (i==0 || i==3 || i==5)
	{
	  csp.add(i,0);
	  csp.add(i,i);
	}
    }  

  csp.print(logfile);  
  
  return 0;
}
