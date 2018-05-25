// ---------------------------------------------------------------------
//
// Copyright (C) 2007 - 2018 by the deal.II authors
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

#include <deal.II/base/logstream.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/vector.h>

#include <complex>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

#include "../tests.h"


// forward declaration of the function that must be provided in the
// .cc files
template <typename number>
void
check();

static constexpr int n_rows = 5;
static constexpr int n_cols = 4;


template <typename number>
void
make_square_matrix(FullMatrix<number> &m)
{
  m.reinit(n_rows, n_rows);
  for (unsigned int i = 0; i < n_rows; ++i)
    for (unsigned int j = 0; j < n_rows; ++j)
      m(i, j) = (i + 1) * (j + 2);
}



template <typename number>
void
make_matrix(FullMatrix<number> &m)
{
  m.reinit(n_rows, n_cols);
  for (unsigned int i = 0; i < n_rows; ++i)
    for (unsigned int j = 0; j < n_cols; ++j)
      m(i, j) = (i + 1) * (j + 2);
}


template <typename number>
void
make_complex_square_matrix(FullMatrix<std::complex<number>> &m)
{
  m.reinit(n_rows, n_rows);
  for (unsigned int i = 0; i < n_rows; ++i)
    for (unsigned int j = 0; j < n_rows; ++j)
      m(i, j) = std::complex<number>((i + 1) * (j + 2), (i + 3) * (j + 4));
}


template <typename number>
void
make_complex_matrix(FullMatrix<std::complex<number>> &m)
{
  m.reinit(n_rows, n_cols);
  for (unsigned int i = 0; i < n_rows; ++i)
    for (unsigned int j = 0; j < n_cols; ++j)
      m(i, j) = std::complex<number>((i + 1) * (j + 2), (i + 3) * (j + 4));
}


template <typename number>
void
make_domain_vector(Vector<number> &v)
{
  v.reinit(n_cols);
  for (unsigned int i = 0; i < n_cols; ++i)
    v(i) = (i + 1);
}



template <typename number>
void
make_range_vector(Vector<number> &v)
{
  v.reinit(n_rows);
  for (unsigned int i = 0; i < n_rows; ++i)
    v(i) = (i + 1);
}



template <typename number>
void
make_complex_domain_vector(Vector<std::complex<number>> &v)
{
  v.reinit(n_cols);
  for (unsigned int i = 0; i < n_cols; ++i)
    v(i) = std::complex<number>(i + 1, i + 3);
}



template <typename number>
void
make_complex_range_vector(Vector<std::complex<number>> &v)
{
  v.reinit(n_rows);
  for (unsigned int i = 0; i < n_rows; ++i)
    v(i) = std::complex<number>(i + 1, i + 3);
}



template <typename number>
void
print_matrix(const FullMatrix<number> &m)
{
  const number tolerance = 100. * std::numeric_limits<number>::epsilon();
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      deallog << i << ' ' << j << ' '
              << filter_out_small_numbers(m(i, j), tolerance) << std::endl;
}



template <typename number>
void
print_matrix(const FullMatrix<std::complex<number>> &m)
{
  const number tolerance = 100. * std::numeric_limits<number>::epsilon();
  for (unsigned int i = 0; i < m.m(); ++i)
    for (unsigned int j = 0; j < m.n(); ++j)
      deallog << i << ' ' << j << ' '
              << filter_out_small_numbers(m(i, j), tolerance) << std::endl;
}



template <typename number>
void
print_vector(const Vector<number> &v)
{
  const typename numbers::NumberTraits<number>::real_type tolerance =
    100. * std::numeric_limits<
             typename numbers::NumberTraits<number>::real_type>::epsilon();
  for (unsigned int i = 0; i < v.size(); ++i)
    deallog << i << ' ' << filter_out_small_numbers(v(i), tolerance)
            << std::endl;
}

template <typename number>
void
display_matrix(FullMatrix<number> M)
{
  const number tolerance = 100. * std::numeric_limits<number>::epsilon();
  deallog << M.m() << "x" << M.n() << " matrix" << std::endl;
  for (unsigned int i = 0; i < M.m(); i++)
    {
      for (unsigned int j = 0; j < M.n(); j++)
        deallog << filter_out_small_numbers(M(i, j), tolerance) << " ";
      deallog << std::endl;
    }
}

template <typename number>
void
display_matrix(FullMatrix<std::complex<number>> M)
{
  const number tolerance = 100. * std::numeric_limits<number>::epsilon();
  deallog << M.m() << "x" << M.n() << " matrix" << std::endl;
  for (unsigned int i = 0; i < M.m(); i++)
    {
      for (unsigned int j = 0; j < M.n(); j++)
        deallog << filter_out_small_numbers(M(i, j), tolerance) << " ";
      deallog << std::endl;
    }
}


template <typename number>
void
fill_matrix(FullMatrix<number> &A)
{
  for (unsigned int i = 0; i < A.m(); i++)
    for (unsigned int j = 0; j < A.n(); j++)
      A(i, j) = number(i * A.n() + j + 1);
}

int
main()
{
  try
    {
      deallog << std::setprecision(2);
      initlog();
      deallog.depth_console(0);

      check<double>();

      check<float>();

      return 0;
    }
  catch (std::exception &exc)
    {
      deallog << std::endl
              << std::endl
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
      deallog << std::endl
              << std::endl
              << "----------------------------------------------------"
              << std::endl;
      deallog << "Unknown exception!" << std::endl
              << "Aborting!" << std::endl
              << "----------------------------------------------------"
              << std::endl;
      return 1;
    };
}
