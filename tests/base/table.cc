// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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


// check equivalence of operator[] and operator() on table objects


#include "../tests.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <cstdlib>

#include <deal.II/base/logstream.h>
#include <deal.II/base/table.h>

const int entries[] = { 11,12,13,21,
                        22,23,31,32,
                        33,58,65,78,
                        80,83,87,91,
                        1, 2, 3, 4,
                        6, 8, 10, 14
                      };

int
main ()
{
  std::ofstream logfile("output");
  deallog << std::fixed;
  deallog << std::setprecision(0);
  deallog.attach(logfile);
  deallog.depth_console(0);
  deallog.threshold_double(1.e-10);

  // a square table
  if (true)
    {
      Table<2,double> Td(3,3);
      Table<2,int> Ti(3,3);

      for (unsigned int i=0; i<9; ++i)
        {
          Td[i/3][i%3] = entries[i];
          Ti[i/3][i%3] = entries[i];
        };

      for (unsigned int i=0; i<3; ++i)
        for (unsigned int j=0; j<3; ++j)
          {
            Assert (Td[i][j] == Td(i,j), ExcInternalError());
            Assert (Ti[i][j] == Ti(i,j), ExcInternalError());
            Assert (Ti[i][j] == Td(i,j), ExcInternalError());

            Assert (*(Td[i].begin()+j) == Td(i,j), ExcInternalError());
            Assert (*(Ti[i].begin()+j) == Ti(i,j), ExcInternalError());
            Assert (*(Ti[i].begin()+j) == Td(i,j), ExcInternalError());

            Assert (&*(Td[i].begin()+j) == &Td(i,j), ExcInternalError());
            Assert (&*(Ti[i].begin()+j) == &Ti(i,j), ExcInternalError());

            deallog << i << " " << j << " " << Td[i][j] << " ok" << std::endl;
          };
    };

  // a rectangular table
  if (true)
    {
      Table<2,double> Td(4,3);
      Table<2,int> Ti(4,3);

      for (unsigned int i=0; i<12; ++i)
        {
          Td(i/3,i%3) = entries[i];
          Ti(i/3,i%3) = entries[i];
        };

      for (unsigned int i=0; i<4; ++i)
        for (unsigned int j=0; j<3; ++j)
          {
            Assert (Td[i][j] == Td(i,j), ExcInternalError());
            Assert (Ti[i][j] == Ti(i,j), ExcInternalError());
            Assert (Ti[i][j] == Td(i,j), ExcInternalError());

            Assert (*(Td[i].begin()+j) == Td(i,j), ExcInternalError());
            Assert (*(Ti[i].begin()+j) == Ti(i,j), ExcInternalError());
            Assert (*(Ti[i].begin()+j) == Td(i,j), ExcInternalError());

            Assert (&*(Td[i].begin()+j) == &Td(i,j), ExcInternalError());
            Assert (&*(Ti[i].begin()+j) == &Ti(i,j), ExcInternalError());

            deallog << i << " " << j << " " << Td[i][j] << " ok" << std::endl;
          };

      for (unsigned int i=0; i<4; ++i)
        Assert (Td[i].end() - Td[i].begin() == 3,
                ExcInternalError());
    };


  // a transposed rectangular table
  if (true)
    {
      TransposeTable<float> Td(4,3);
      TransposeTable<int> Ti(4,3);

      // Fill the float matrix in row
      // first ordering
      for (unsigned int i=0; i<12; ++i)
        {
          Td(i/3,i%3) = entries[i];
        };

      // This fills the integer
      // matrix in column first
      // ordering
      Ti.fill(entries);
      // Output both matrices...
      // They should be different
      for (unsigned int i=0; i<4; ++i)
        {
          for (unsigned int j=0; j<3; ++j)
            {
              deallog << '\t' << Ti(i,j) << '/' << Td(i,j);
            }
          deallog << std::endl;
        }
    }


  // a 1d-table
  if (true)
    {
      const unsigned int N=10;
      Table<1,double> Td(N);
      Table<1,int> Ti(N);

      for (unsigned int i=0; i<N; ++i)
        {
          Td(i) = entries[i];
          Ti(i) = entries[i];
        };

      for (unsigned int i=0; i<N; ++i)
        {
          Assert (Td[i] == Td(i), ExcInternalError());
          Assert (Ti[i] == Ti(i), ExcInternalError());
          Assert (Ti[i] == Td(i), ExcInternalError());

          deallog << i << " " << Td[i] << " ok" << std::endl;
        };
    };

  // a 3d-table
  if (true)
    {
      const unsigned int I=4, J=3, K=2;
      Table<3,double> Td(I,J,K);
      Table<3,int> Ti(I,J,K);

      unsigned int index=0;
      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k, ++index)
            {
              Td(i,j,k) = entries[index];
              Ti(i,j,k) = entries[index];
            };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            {
              Assert (Td[i][j][k] == Td(i,j,k), ExcInternalError());
              Assert (Ti[i][j][k] == Ti(i,j,k), ExcInternalError());
              Assert (Ti[i][j][k] == Td(i,j,k), ExcInternalError());

              Assert (*(Td[i][j].begin()+k) == Td(i,j,k), ExcInternalError());
              Assert (*(Ti[i][j].begin()+k) == Ti(i,j,k), ExcInternalError());
              Assert (*(Ti[i][j].begin()+k) == Td(i,j,k), ExcInternalError());

              Assert (&*(Td[i][j].begin()+k) == &Td(i,j,k), ExcInternalError());
              Assert (&*(Ti[i][j].begin()+k) == &Ti(i,j,k), ExcInternalError());

              deallog << i << " " << j << " " << k << " "
                      << Td[i][j][k] << " ok" << std::endl;
            };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          Assert (Td[i][j].end() - Td[i][j].begin() ==
                  static_cast<signed int>(K),
                  ExcInternalError());
    };

  // a 4d-table
  if (true)
    {
      const unsigned int I=5, J=4, K=3, L=2;
      Table<4,double> Td(I,J,K,L);
      Table<4,int> Ti(I,J,K,L);

      unsigned int index=0;
      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l, ++index)
              {
                // Take a simple
                // number generator
                // rather than
                // dreaming up lots
                // of `arbitrary'
                // numbers in the
                // entries array
                Td(i,j,k,l) = 2*index+1;
                Ti(i,j,k,l) = 2*index+1;
              };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              {
                Assert (Td[i][j][k][l] == Td(i,j,k,l), ExcInternalError());
                Assert (Ti[i][j][k][l] == Ti(i,j,k,l), ExcInternalError());
                Assert (Ti[i][j][k][l] == Td(i,j,k,l), ExcInternalError());

                Assert (*(Td[i][j][k].begin()+l) == Td(i,j,k,l), ExcInternalError());
                Assert (*(Ti[i][j][k].begin()+l) == Ti(i,j,k,l), ExcInternalError());
                Assert (*(Ti[i][j][k].begin()+l) == Td(i,j,k,l), ExcInternalError());

                Assert (&*(Td[i][j][k].begin()+l) == &Td(i,j,k,l), ExcInternalError());
                Assert (&*(Ti[i][j][k].begin()+l) == &Ti(i,j,k,l), ExcInternalError());

                deallog << i << " " << j << " " << k << " " << l << " "
                        << Td[i][j][k][l] << " ok" << std::endl;
              };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            Assert (Td[i][j][k].end() - Td[i][j][k].begin() ==
                    static_cast<signed int>(L),
                    ExcDimensionMismatch(Td[i][j][k].end() - Td[i][j][k].begin(),
                                         static_cast<signed int>(L)));
    };

  // 5d-table
  if (true)
    {
      const unsigned int I=6, J=2, K=3, L=4, M=5;
      Table<5,double> Td(I,J,K,L,M);
      Table<5,int> Ti(I,J,K,L,M);

      unsigned int index=0;
      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              for (unsigned int m=0; m<M; ++m, ++index)
                {
                  Td(i,j,k,l,m) = 2*index+1;
                  Ti(i,j,k,l,m) = 2*index+1;
                };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              for (unsigned int m=0; m<M; ++m)
                {
                  Assert (Td[i][j][k][l][m] == Td(i,j,k,l,m), ExcInternalError());
                  Assert (Ti[i][j][k][l][m] == Ti(i,j,k,l,m), ExcInternalError());
                  Assert (Ti[i][j][k][l][m] == Td(i,j,k,l,m), ExcInternalError());

                  Assert (*(Td[i][j][k][l].begin()+m) == Td(i,j,k,l,m), ExcInternalError());
                  Assert (*(Ti[i][j][k][l].begin()+m) == Ti(i,j,k,l,m), ExcInternalError());
                  Assert (*(Ti[i][j][k][l].begin()+m) == Td(i,j,k,l,m), ExcInternalError());

                  Assert (&*(Td[i][j][k][l].begin()+m) == &Td(i,j,k,l,m), ExcInternalError());
                  Assert (&*(Ti[i][j][k][l].begin()+m) == &Ti(i,j,k,l,m), ExcInternalError());

                  deallog << i << " " << j << " " << k << " " << l << " " << m << " "
                          << Td[i][j][k][l][m] << " ok" << std::endl;
                };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              Assert (Td[i][j][k][l].end() - Td[i][j][k][l].begin() ==
                      static_cast<signed int>(M),
                      ExcDimensionMismatch(Td[i][j][k][l].end() - Td[i][j][k][l].begin(),
                                           static_cast<signed int>(M)));
    };

  // a 6d-table
  if (true)
    {
      const unsigned int I=6, J=2, K=4, L=3, M=5, N=7;
      Table<6,double> Td(I,J,K,L,M,N);
      Table<6,int> Ti(I,J,K,L,M,N);

      unsigned int index=0;
      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              for (unsigned int m=0; m<M; ++m)
                for (unsigned int n=0; n<N; ++n, ++index)
                  {
                    Td(i,j,k,l,m,n) = 2*index+1;
                    Ti(i,j,k,l,m,n) = 2*index+1;
                  };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              for (unsigned int m=0; m<M; ++m)
                for (unsigned int n=0; n<N; ++n)
                  {
                    Assert (Td[i][j][k][l][m][n] == Td(i,j,k,l,m,n), ExcInternalError());
                    Assert (Ti[i][j][k][l][m][n] == Ti(i,j,k,l,m,n), ExcInternalError());
                    Assert (Ti[i][j][k][l][m][n] == Td(i,j,k,l,m,n), ExcInternalError());

                    Assert (*(Td[i][j][k][l][m].begin()+n) == Td(i,j,k,l,m,n), ExcInternalError());
                    Assert (*(Ti[i][j][k][l][m].begin()+n) == Ti(i,j,k,l,m,n), ExcInternalError());
                    Assert (*(Ti[i][j][k][l][m].begin()+n) == Td(i,j,k,l,m,n), ExcInternalError());

                    Assert (&*(Td[i][j][k][l][m].begin()+n) == &Td(i,j,k,l,m,n), ExcInternalError());
                    Assert (&*(Ti[i][j][k][l][m].begin()+n) == &Ti(i,j,k,l,m,n), ExcInternalError());

                    deallog << i << " " << j << " " << k << " " << l << " " << m << " " << n << " "
                            << Td[i][j][k][l][m][n] << " ok" << std::endl;
                  };

      for (unsigned int i=0; i<I; ++i)
        for (unsigned int j=0; j<J; ++j)
          for (unsigned int k=0; k<K; ++k)
            for (unsigned int l=0; l<L; ++l)
              for (unsigned int m=0; m<M; ++m)
                Assert (Td[i][j][k][l][m].end() - Td[i][j][k][l][m].begin() ==
                        static_cast<signed int>(N),
                        ExcDimensionMismatch(Td[i][j][k][l][m].end() - Td[i][j][k][l][m].begin(),
                                             static_cast<signed int>(N)));
    };
}



