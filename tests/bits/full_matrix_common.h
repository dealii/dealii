// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2013 by the deal.II authors
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


// common framework for the various full_matrix_*.cc tests

#include "../tests.h"
#include <deal.II/base/logstream.h>
#include <deal.II/lac/vector.h>
#include <deal.II/lac/full_matrix.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <complex>


// forward declaration of the function that must be provided in the
// .cc files
template <typename number>
void
check ();

// forward declaration of a variable with the name of the output file
extern std::string output_file_name;



template <typename number>
void make_matrix (FullMatrix<number> &m)
{
  m.reinit (5,5);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      m(i,j) = (i+1)*(j+2);
}


template <typename number>
void make_complex_matrix (FullMatrix<std::complex<number> > &m)
{
  m.reinit (5,5);
  for (unsigned int i=0; i<5; ++i)
    for (unsigned int j=0; j<5; ++j)
      m(i,j) = std::complex<number>((i+1)*(j+2), (i+3)*(j+4));
}


template <typename number>
void make_vector (Vector<number> &v)
{
  v.reinit (5);
  for (unsigned int i=0; i<5; ++i)
    v(i) = (i+1);
}



template <typename number>
void make_complex_vector (Vector<std::complex<number> > &v)
{
  v.reinit (5);
  for (unsigned int i=0; i<5; ++i)
    v(i) = std::complex<number>(i+1, i+3);
}



template <typename number>
void
print_matrix (const FullMatrix<number> &m)
{
  for (unsigned int i=0; i<m.m(); ++i)
    for (unsigned int j=0; j<m.n(); ++j)
      deallog << i << ' ' << j << ' ' << m(i,j)
              << std::endl;
}



template <typename number>
void
print_vector (const Vector<number> &v)
{
  for (unsigned int i=0; i<v.size(); ++i)
    deallog << i << ' ' << v(i)
            << std::endl;
}

template <typename number>
void
display_matrix(FullMatrix<number> M)
{
  deallog<<M.m()<<"x"<<M.n()<<" matrix"<<std::endl;
  for (unsigned int i=0; i<M.m(); i++)
    {
      for (unsigned int j=0; j<M.n(); j++)
        deallog<<M(i,j)<<" ";
      deallog<<std::endl;
    }
}


template <typename number>
void
fill_matrix(FullMatrix<number> &A)
{
  for (unsigned int i=0; i<A.m(); i++)
    for (unsigned int j=0; j<A.n(); j++)
      A(i,j)=number(i*A.n() + j+1);
}

int
main()
{
  try
    {
      std::ofstream logfile(output_file_name.c_str());
      deallog << std::setprecision (2);
      deallog.attach(logfile);
      deallog.depth_console(0);
      deallog.threshold_double(1.e-7);

      check<double> ();
      check<float> ();

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Exception on processing: " << std::endl
              << exc.what() << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    }
  catch (...)
    {
      deallog << std::endl << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}

