//----------------------------------------------------------------------
//    $Id$
//    Version: $Name$ 
//
//    Copyright (C) 2007 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------------

// Test if block indices are handled properly

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_indices.h>

#include <fstream>


void test (const BlockIndices& idx)
{
  const unsigned int n = idx.size();
  deallog << "blocks " << idx.size() << std::endl
	  << "elements " << idx.total_size() << std::endl
	  << "block sizes ";
  for (unsigned i=0;i<n;++i)
    deallog << '\t' << idx.block_size(i);
  deallog << std::endl << "Block start ";
  for (unsigned i=0;i<n;++i)
    deallog << '\t' << idx.block_start(i);

  deallog << std::endl;
  
  for (unsigned int i=0;i<idx.total_size();++i)
    {
      const unsigned int b = idx.global_to_local(i).first;
      const unsigned int j = idx.global_to_local(i).second;
      deallog << ' '<< i << ':' << b << ':' << j;
    }
  
  deallog << std::endl;

  for (unsigned int b=0;b<n;++b)
    for (unsigned int j=0;j<idx.block_size(b);++j)
      {
	const unsigned int i = idx.local_to_global(b,j);
	deallog << ' '<< i << ':' << b << ':' << j;
      }
  
  deallog << std::endl;
}


int main()
{
  std::ofstream logfile("block_indices/output");
  deallog.attach(logfile);
  deallog.depth_console(0);
  
  BlockIndices bi1(3);
  test(bi1);
  bi1.reinit(3,4);
  test(bi1);
  
  std::vector<unsigned int> v(4);
  for (unsigned int i=0;i<v.size();++i)
    v[i] = 4-i;
  
  BlockIndices bi2(v);
  test(bi2);
}
