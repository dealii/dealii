//----------------------------  testmatrix.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  testmatrix.cc  ---------------------------


#include "testmatrix.h"
#include <lac/sparse_matrix.h>
#include <lac/vector.h>

FDMatrix::FDMatrix(unsigned int nx, unsigned int ny)
		:
		nx(nx), ny(ny)
{}

void
FDMatrix::build_structure(SparsityPattern& structure) const
{
  for(unsigned int i=0;i<=ny-2;i++)
    {
      for(unsigned int j=0;j<=nx-2; j++)
	{
					   // Number of the row to be entered
	  unsigned int row = j+(nx-1)*i;
	  structure.add(row, row);
	  if (j>0)
	    {
	      structure.add(row-1, row);
	      structure.add(row, row-1);
	    }
	  if (j<nx-2)
	    {
	      structure.add(row+1, row);
	      structure.add(row, row+1);
	    }
	  if (i>0)
	    {
	      structure.add(row-(nx-1), row);
	      structure.add(row, row-(nx-1));
	    }
	  if (i<ny-2)
	    {
	      structure.add(row+(nx-1), row);
	      structure.add(row, row+(nx-1));
	    }
	}
    }
}

template<typename number>
void
FDMatrix::laplacian(SparseMatrix<number>& A) const
{
   for(unsigned int i=0;i<=ny-2;i++)
    {
      for(unsigned int j=0;j<=nx-2; j++)
	{
					   // Number of the row to be entered
	  unsigned int row = j+(nx-1)*i;
	  
	  A.set(row, row, 4.);
	  if (j>0)
	    {
	      A.set(row-1, row, -1.);
	      A.set(row, row-1, -1.);
	    }
	  if (j<nx-2)
	    {
	      A.set(row+1, row, -1.);
	      A.set(row, row+1, -1.);
	    }
	  if (i>0)
	    {
	      A.set(row-(nx-1), row, -1.);
	      A.set(row, row-(nx-1), -1.);
	    }
	  if (i<ny-2)
	    {
	      A.set(row+(nx-1), row, -1.);
	      A.set(row, row+(nx-1), -1.);
	    }
	}
    } 
}

template<typename number>
void
FDMatrix::gnuplot_print(std::ostream& s, const Vector<number>& V) const
{
   for(unsigned int i=0;i<=ny-2;i++)
    {
      for(unsigned int j=0;j<=nx-2; j++)
	{
					   // Number of the row to be entered
	  unsigned int row = j+(nx-1)*i;
	  s << (j+1)/(float)nx << '\t' << (i+1)/(float)ny << '\t' << V(row) << std::endl;
	}
      s << std::endl;
    } 
   s << std::endl;
}



template void FDMatrix::laplacian(SparseMatrix<long double>&) const;
template void FDMatrix::laplacian(SparseMatrix<double>&) const;
template void FDMatrix::laplacian(SparseMatrix<float>&) const;
template void FDMatrix::gnuplot_print(std::ostream& s, const Vector<long double>& V) const;
template void FDMatrix::gnuplot_print(std::ostream& s, const Vector<double>& V) const;
template void FDMatrix::gnuplot_print(std::ostream& s, const Vector<float>& V) const;

