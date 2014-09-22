// ---------------------------------------------------------------------
//
// Copyright (C) 2000 - 2013 by the deal.II authors
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



#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/block_sparsity_pattern.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/block_vector.h>
#include <fstream>
#include <iomanip>
#include <vector>
#include <list>

void test ()
{
  deallog.push("BlockIndices");

  std::vector<types::global_dof_index> ivector(4);
  ivector[0] = 3;
  ivector[1] = 0;
  ivector[2] = 1;
  ivector[3] = 2;

  const std::vector<types::global_dof_index> vector_indices(ivector);

  BlockIndices i1(ivector);
  BlockIndices i2 = i1;
  BlockIndices i3;
  // no output expected here
  deallog.push("empty constructor");
  for (unsigned int i=0 ; i<i3.size() ; ++i)
    deallog << i << '\t' << i3.local_to_global(i,0) << std::endl;
  for (unsigned int i=0 ; i<i3.total_size() ; ++i)
    deallog << i
            << '\t' << i3.global_to_local(i).first
            << '\t' << i3.global_to_local(i).second
            << std::endl;
  deallog.pop();


  i3.reinit(ivector);

  deallog.push("global->local");

  unsigned int n = i1.total_size();
  for (unsigned int i=0; i<n; ++i)
    {
      deallog << i
              << '\t' << i1.global_to_local(i).first
              << '\t' << i1.global_to_local(i).second
              << '\t' << i2.global_to_local(i).first
              << '\t' << i2.global_to_local(i).second
              << '\t' << i3.global_to_local(i).first
              << '\t' << i3.global_to_local(i).second
              << std::endl;
    }

  deallog.pop();

  deallog.push("local->global");
  for (unsigned int i=0 ; i<i1.size() ; ++i)
    for (unsigned int j=0 ; j<ivector[i] ; ++j)
      deallog << i << '\t' << j << '\t'
              << i1.local_to_global(i,j) << std::endl;

  deallog.pop();

  deallog.push("reinit");

  ivector.insert(ivector.begin(), 5);
  i1.reinit(ivector);
  n = i1.total_size();
  for (unsigned int i=0; i<n; ++i)
    {
      deallog << i
              << '\t' << i1.global_to_local(i).first
              << '\t' << i1.global_to_local(i).second
              << std::endl;
    }
  deallog << "---" << std::endl;

  ivector.erase(ivector.begin());
  ivector.erase(ivector.begin());
  ivector.erase(ivector.begin());
  i1.reinit(ivector);
  n = i1.total_size();
  for (unsigned int i=0; i<n; ++i)
    {
      deallog << i
              << '\t' << i1.global_to_local(i).first
              << '\t' << i1.global_to_local(i).second
              << std::endl;
    }
  deallog.pop();

  deallog.pop ();
  deallog.push("BlockVector");
  // initialization by an iterator
  // range
  deallog.push ("Constructor with iterators");
  double array[] = { 0, 1, 2, 3, 4, 5 };
  BlockVector<double> v1(vector_indices, &array[0], &array[6]);
  for (unsigned int i=0; i<v1.size(); ++i)
    deallog << v1(i) << ' ';
  deallog << std::endl;

  // same test, but do not initialize
  // from double*'s, but from
  // std::list iterators.
  std::list<double> l(&array[0], &array[6]);
  BlockVector<double> v2(vector_indices, l.begin(), l.end());
  for (unsigned int i=0; i<v2.n_blocks(); ++i)
    for (unsigned int j=0; j<v2.block(i).size(); ++j)
      deallog << i << '\t' << j << '\t' <<v2.block(i)(j) << std::endl;

  deallog.pop();
  deallog.push("reinit block");
  v2.block(1).reinit(5);
  v2.collect_sizes();
  for (unsigned int i=0; i<v2.size(); ++i)
    deallog << v2(i) << ' ';
  deallog << std::endl;

  deallog.pop ();
  deallog.pop ();
}




int main ()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(3);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // do the same weird stuff as in
  // tests/base/reference.cc
#if __GNUC__ != 2
  std::basic_streambuf<char> *old_cerr_buf = std::cerr.rdbuf();
#else
  streambuf *old_cerr_buf = std::cerr.rdbuf();
#endif
  std::cerr.rdbuf(logfile.rdbuf());

  try
    {
      test ();
    }
  catch (std::exception &e)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << e.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 0;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      // abort
      return 0;
    };

  std::cerr.rdbuf(old_cerr_buf);
}

