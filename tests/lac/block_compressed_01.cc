//----------------------------------------------------------------------
//    $Id: block_compressed_simple_sparsity_pattern_indexset_01.cc 25686 2012-07-04 14:18:43Z bangerth $
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


// Test BlockCompressedSimpleSparsityPattern column_number


#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>

#include <iomanip>
#include <fstream>

int main()
{
  std::ofstream logfile("block_compressed_01/output");
  logfile.setf(std::ios::fixed);
  deallog << std::setprecision(2);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  BlockCompressedSimpleSparsityPattern csp(2,2);
  csp.block(0,0).reinit(2,2);
  csp.block(1,0).reinit(1,2);
  csp.block(0,1).reinit(2,1);
  csp.block(1,1).reinit(1,1);
  csp.collect_sizes();
  for (int i=0;i<3;++i)
    csp.add(i,i);
  csp.add(1,0);
  csp.add(1,2);
  csp.compress();
  
  deallog << "blocks: " << csp.n_block_rows() << "x" << csp.n_block_cols() << std::endl;
  deallog << "size: " << csp.n_rows() << "x" << csp.n_cols() << std::endl;
  deallog << "size block(1,0):" << csp.block(1,0).n_rows() << "x" << csp.block(1,0).n_cols() << std::endl;
  
  csp.print(logfile);

  deallog << std::endl;
  for (unsigned int i=0;i<3;++i)
    {
      deallog << "row " << i << ": ";
      
      for (unsigned int j=0;j<csp.row_length(i);++j)
	deallog << csp.column_number(i,j) << " ";

      deallog << std::endl;
    }
  
  return 0;
}
