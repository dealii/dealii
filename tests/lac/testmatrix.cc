// $Id$

#include "testmatrix.h"
#include <lac/sparsematrix.h>

FDMatrix::FDMatrix(unsigned int nx, unsigned int ny)
		:
		nx(nx), ny(ny)
{}

void
FDMatrix::build_structure(SparseMatrixStruct& structure) const
{
  for(unsigned int i=0;i<=ny;i++)
    {
      for(unsigned int j=0;j<=nx; j++)
	{
					   // Number of the row to be entered
	  unsigned int row = j+(nx+1)*i;
	  structure.add(row, row);
	  if (j>0)
	    {
	      structure.add(row-1, row);
	      structure.add(row, row-1);
	    }
	  if (j<nx)
	    {
	      structure.add(row+1, row);
	      structure.add(row, row+1);
	    }
	  if (i>0)
	    {
	      structure.add(row-nx, row);
	      structure.add(row, row-nx);
	    }
	  if (i<ny)
	    {
	      structure.add(row+nx, row);
	      structure.add(row, row+nx);
	    }
	}
    }
}

template<typename number>
void
FDMatrix::laplacian(SparseMatrix<number>& A) const
{
   for(unsigned int i=0;i<=ny;i++)
    {
      for(unsigned int j=0;j<=nx; j++)
	{
					   // Number of the row to be entered
	  unsigned int row = j+(nx+1)*i;
	  
	  A.set(row, row, 4.);
	  if (j>0)
	    {
	      A.set(row-1, row, -1.);
	      A.set(row, row-1, -1.);
	    }
	  if (j<nx)
	    {
	      A.set(row+1, row, -1.);
	      A.set(row, row+1, -1.);
	    }
	  if (i>0)
	    {
	      A.set(row-nx, row, -1.);
	      A.set(row, row-nx, -1.);
	    }
	  if (i<ny)
	    {
	      A.set(row+nx, row, -1.);
	      A.set(row, row+nx, -1.);
	    }
	}
    } 
}

template void FDMatrix::laplacian(SparseMatrix<double>&) const;
template void FDMatrix::laplacian(SparseMatrix<float>&) const;
