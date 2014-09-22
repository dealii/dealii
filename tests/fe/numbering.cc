// ---------------------------------------------------------------------
//
// Copyright (C) 2013 by the deal.II authors
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
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_tools.h>
#include <vector>
#include <fstream>


std::ofstream logfile ("output");


template <int dim>
void check (const FE_Q<dim> &fe)
{
  Assert (fe.n_components() == 1, ExcInternalError());

  std::vector<unsigned int> hierarchic_to_lexicographic_numbering (fe.dofs_per_cell);

  const unsigned int dofs_per_cell = fe.dofs_per_cell;
  // polynomial degree
  const unsigned int degree = fe.dofs_per_line+1;
  // number of grid points in each
  // direction
  const unsigned int n = degree+1;

  // the following lines of code are
  // somewhat odd, due to the way the
  // hierarchic numbering is
  // organized. if someone would
  // really want to understand these
  // lines, you better draw some
  // pictures where you indicate the
  // indices and orders of vertices,
  // lines, etc, along with the
  // numbers of the degrees of
  // freedom in hierarchical and
  // lexicographical order
  switch (dim)
    {
    case 1:
    {
      hierarchic_to_lexicographic_numbering[0] = 0;
      hierarchic_to_lexicographic_numbering[1] = dofs_per_cell-1;
      for (unsigned int i=2; i<dofs_per_cell; ++i)
        hierarchic_to_lexicographic_numbering[i] = i-1;

      break;
    }

    case 2:
    {
      unsigned int next_index = 0;
      // first the four vertices
      hierarchic_to_lexicographic_numbering[next_index++] = 0;
      hierarchic_to_lexicographic_numbering[next_index++] = n-1;
      hierarchic_to_lexicographic_numbering[next_index++] = n*(n-1);
      hierarchic_to_lexicographic_numbering[next_index++] = n*n-1;

      // left   line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (1+i)*n;

      // right  line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (2+i)*n-1;

      // bottom line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = 1+i;

      // top    line
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n*(n-1)+i+1;

      // inside quad
      Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
              ExcInternalError());
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = n*(i+1)+j+1;

      Assert (next_index == fe.dofs_per_cell, ExcInternalError());

      break;
    }

    case 3:
    {
      unsigned int next_index = 0;
      // first the eight vertices
      hierarchic_to_lexicographic_numbering[next_index++] = 0;                 // 0
      hierarchic_to_lexicographic_numbering[next_index++] = (      1)*degree;  // 1
      hierarchic_to_lexicographic_numbering[next_index++] = (    n  )*degree;  // 2
      hierarchic_to_lexicographic_numbering[next_index++] = (    n+1)*degree;  // 3
      hierarchic_to_lexicographic_numbering[next_index++] = (n*n    )*degree;  // 4
      hierarchic_to_lexicographic_numbering[next_index++] = (n*n  +1)*degree;  // 5
      hierarchic_to_lexicographic_numbering[next_index++] = (n*n+n  )*degree;  // 6
      hierarchic_to_lexicographic_numbering[next_index++] = (n*n+n+1)*degree;  // 7

      // line 0
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n;
      // line 1
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n;
      // line 2
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = 1+i;
      // line 3
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = 1+i+n*(n-1);

      // line 4
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*n*n+(i+1)*n;
      // line 5
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*(n*n+1)+(i+1)*n;
      // line 6
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n*n*(n-1)+i+1;
      // line 7
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n*n*(n-1)+i+1+n*(n-1);

      // line 8
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n;
      // line 9
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n*n;
      // line 10
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n*(n-1);
      // line 11
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        hierarchic_to_lexicographic_numbering[next_index++] = n-1+(i+1)*n*n+n*(n-1);


      // inside quads
      Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
              ExcInternalError());
      // face 0
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n*(j+1);
      // face 1
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = (i+1)*n*n+n-1+n*(j+1);
      // face 2, note the orientation!
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = (j+1)*n*n+i+1;
      // face 3, note the orientation!
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = (j+1)*n*n+n*(n-1)+i+1;
      // face 4
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = n*(i+1)+j+1;
      // face 5
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          hierarchic_to_lexicographic_numbering[next_index++] = (n-1)*n*n+n*(i+1)+j+1;

      // inside hex
      Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
              ExcInternalError());
      for (unsigned int i=0; i<fe.dofs_per_line; ++i)
        for (unsigned int j=0; j<fe.dofs_per_line; ++j)
          for (unsigned int k=0; k<fe.dofs_per_line; ++k)
            hierarchic_to_lexicographic_numbering[next_index++] = n*n*(i+1)+n*(j+1)+k+1;

      Assert (next_index == fe.dofs_per_cell, ExcInternalError());

      break;
    }

    default:
      Assert (false, ExcNotImplemented());
    }

  // now check with data from the
  // lib: there we have the mapping
  // the other way round, and we
  // check that the concatenation of
  // the two mappings is the
  // identity. output the two maps to
  // generate some output for
  // automatic comparison
  std::vector<unsigned int> l2h (fe.dofs_per_cell);
  FETools::lexicographic_to_hierarchic_numbering (fe, l2h);
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    {
      Assert (l2h[hierarchic_to_lexicographic_numbering[i]] == i,
              ExcInternalError());
      logfile << dim << "d, degree=" << degree << ": "
              << l2h[i]
              << ' '
              << hierarchic_to_lexicographic_numbering[i]
              << std::endl;
    };

  // finally, we also have the
  // forward map in the lib, so check
  // for equality
  std::vector<unsigned int> h2l (fe.dofs_per_cell);
  FETools::hierarchic_to_lexicographic_numbering (fe, h2l);
  Assert (hierarchic_to_lexicographic_numbering == h2l,
          ExcInternalError());
}



template <int dim>
void check_dim ()
{
  for (unsigned int degree=1; degree<6; ++degree)
    {
      FE_Q<dim> fe(degree);
      check (fe);
    };
}


int main ()
{
  deallog << std::setprecision(2);
  deallog << std::fixed;
  deallog.attach(logfile);
  deallog.depth_console (0);

  check_dim<1> ();
  check_dim<2> ();
  check_dim<3> ();
}


