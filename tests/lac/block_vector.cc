//--------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//--------------------------------------------------------------------


#include <base/logstream.h>
#include <lac/block_sparsity_pattern.h>
#include <lac/block_sparse_matrix.h>
#include <lac/block_vector.h>
#include <fstream>
#include <vector>

void abort ()
{}


void test ()
{
  deallog.push("BlockIndices");
  
  vector<unsigned int> ivector(4);
  ivector[0] = 3;
  ivector[1] = 0;
  ivector[2] = 1;
  ivector[3] = 2;
  
  BlockIndices i1(ivector);
  BlockIndices i2 = i1;
  BlockIndices i3(13);
  i3.reinit(ivector);

  deallog.push("global->local");
  
  unsigned int n = i1.total_size();
  for (unsigned int i=0;i<n;++i)
    {
      deallog << i
	      << '\t' << i1.global_to_local(i).first
	      << '\t' << i1.global_to_local(i).second
	      << '\t' << i2.global_to_local(i).first
	      << '\t' << i2.global_to_local(i).second
	      << '\t' << i3.global_to_local(i).first
	      << '\t' << i3.global_to_local(i).second
	      << endl;
    }

  deallog.pop();

  deallog.push("local->global");
  for (unsigned int i=0 ; i<i1.size() ; ++i)
    for (unsigned int j=0 ; j<ivector[i] ; ++j)
      deallog << i << '\t' << j << '\t'
	      << i1.local_to_global(i,j) << endl;
  
  deallog.pop();

  deallog.push("reinit");
  
  ivector.insert(ivector.begin(), 5);
  i1.reinit(ivector);
  n = i1.total_size();
  for (unsigned int i=0;i<n;++i)
    {
      deallog << i
	      << '\t' << i1.global_to_local(i).first
	      << '\t' << i1.global_to_local(i).second
	      << endl;
    }
  deallog << "---" << endl;
  
  ivector.erase(ivector.begin());
  ivector.erase(ivector.begin());
  ivector.erase(ivector.begin());
  i1.reinit(ivector);
  n = i1.total_size();
  for (unsigned int i=0;i<n;++i)
    {
      deallog << i
	      << '\t' << i1.global_to_local(i).first
	      << '\t' << i1.global_to_local(i).second
	      << endl;
    }
  deallog.pop();

  deallog << "Range" << endl;
  unsigned int i = i1.global_to_local(3).first;
  i = i1.local_to_global(1,2);
  i = i1.local_to_global(2,0);
}




int main ()
{
  ofstream logfile("block_vector.output");
  logfile.setf(ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  cerr = logfile;
  
  try
    {
      test ();
    }
  catch (exception &e)
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Exception on processing: " << e.what() << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 2;
    }
  catch (...) 
    {
      cerr << endl << endl
	   << "----------------------------------------------------"
	   << endl;
      cerr << "Unknown exception!" << endl
	   << "Aborting!" << endl
	   << "----------------------------------------------------"
	   << endl;
				       // abort
      return 3;
    };
  
  
  return 0;
};

