//----------------------------  $RCSfile$  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  $RCSfile$  ---------------------------


#include <base/logstream.h>
#include <lac/sparsity_pattern.h>
#include "testmatrix.h"
#include <fstream>
#include <list>
#include <set>


int
main ()
{
  ofstream logfile("sparsity_pattern.output");
  logfile.setf(ios::fixed);
  logfile.precision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);

				   // generate usual sparsity pattern
  const unsigned int N = 15;
  SparsityPattern sp1((N-1)*(N-1), (N-1)*(N-1), 5);
  FDMatrix(N,N).build_structure (sp1);
  sp1.compress ();
  deallog << sp1.n_rows() << " " << sp1.n_cols() << " "
	  << sp1.bandwidth() << " " << sp1.n_nonzero_elements()
	  << std::endl;
  for (unsigned int i=0; i<sp1.n_rows(); ++i)
    deallog << sp1.row_length(i) << std::endl;
  sp1.print_gnuplot(deallog.get_file_stream());

				   // generate copy of sp1 with extra
				   // off-diagonals
  SparsityPattern sp2(sp1, 10, 2);
  sp2.compress ();
  deallog << sp2.n_rows() << " " << sp2.n_cols() << " "
	  << sp2.bandwidth() << " " << sp2.n_nonzero_elements()
	  << std::endl;
  for (unsigned int i=0; i<sp2.n_rows(); ++i)
    deallog << sp2.row_length(i) << std::endl;
  sp2.print_gnuplot(deallog.get_file_stream());


				   // generate copy of sp1 with extra
				   // off-diagonals, add some
				   // non-symmetric elements and
				   // symmetrize again
  SparsityPattern sp3(sp1, (N-1)*(N-1), 2);
  for (unsigned int i=0; i<(N-1)*(N-1); ++i)
    sp3.add (0,i);
  sp3.symmetrize ();
  sp3.compress ();
  deallog << sp3.n_rows() << " " << sp3.n_cols() << " "
	  << sp3.bandwidth() << " " << sp3.n_nonzero_elements()
	  << std::endl;
  for (unsigned int i=0; i<sp3.n_rows(); ++i)
    deallog << sp3.row_length(i) << std::endl;
  sp3.print_gnuplot(deallog.get_file_stream());


				   // now test the copy_from
				   // function. for this copy over the
				   // column indices, but in different
				   // order as the order should not
				   // matter to that function
  std::list<std::set<unsigned int,std::greater<unsigned int> > > sparsity;
  for (unsigned int row=0; row<sp3.n_rows(); ++row)
    {
      sparsity.push_back (std::set<unsigned int,std::greater<unsigned int> >());
      for (const unsigned int *p=sp3.get_column_numbers()+sp3.get_rowstart_indices()[row];
	   p != sp3.get_column_numbers()+sp3.get_rowstart_indices()[row+1]; ++p)
	sparsity.back().insert (*p);
    };
  SparsityPattern sp4;
  sp4.copy_from ((N-1)*(N-1), (N-1)*(N-1),
		 sparsity.begin(), sparsity.end());

				   // now check for equivalence of sp3 and sp4
  for (unsigned int row=0; row<sp3.n_rows(); ++row)
    {
      const unsigned int *sp3_p=sp3.get_column_numbers()+sp3.get_rowstart_indices()[row];
      const unsigned int *sp4_p=sp4.get_column_numbers()+sp4.get_rowstart_indices()[row];
      for (; sp3_p != sp3.get_column_numbers()+sp3.get_rowstart_indices()[row+1]; ++sp3_p, ++sp4_p)
	Assert (*sp3_p == *sp4_p, ExcInternalError());
    };
};

  
  
