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

FDMGTransfer::FDMGTransfer(unsigned int nx, unsigned int ny,
			   unsigned int nlevels)
		:
		structures(nlevels-1), matrices(nlevels-1)
{
  unsigned int power = 1 << nlevels;
  
  Assert ((nx%power)==0, ExcDivide(nx,power));
  Assert ((ny%power)==0, ExcDivide(ny,power));
  
  nx /= power;
  ny /= power;
  
  for (unsigned int level = 0; level < nlevels-1; ++ level)
    {
      build_matrix(nx,ny,structures[level],matrices[level]);
      nx *= 2;
      ny *= 2;
    }
}

void
FDMGTransfer::build_matrix(unsigned int nx, unsigned int ny,
			   SparseMatrixStruct& structure, SparseMatrix<float>& matrix)
{
  structure.reinit((nx+1)*(ny+1),(2*nx+1)*(2*ny+1),9);
  
				   // Traverse all points of coarse grid
  for (unsigned int i=0 ; i<=ny; ++i)
    for (unsigned int j=0 ; j<=nx; ++j)
      {
					 // coarse grid point number
	unsigned int ncoarse =j+(nx+1)*i;
					 // same point on fine grid
	unsigned int nfine = 2*j+(4*nx+2)*i;
	
	structure.add(ncoarse,nfine);
					 // left
	if (j>0)
	  {
	    structure.add(ncoarse,nfine-1);
					     // lower left
	    if (i>0)
	      structure.add(ncoarse,nfine-(2*nx+1)-1);
					     // upper left
	    if (i<ny)
	      structure.add(ncoarse,nfine+(2*nx+1)-1);
	  }
					 // right
	if (j<nx)
	  {
	    structure.add(ncoarse, nfine+1);
					     // lower right
	    if (i>0)
	      structure.add(ncoarse,nfine-(2*nx+1)+1);
					     // upper right
	    if (i<ny)
	      structure.add(ncoarse,nfine+(2*nx+1)+1);
	  }

					     // lower
	    if (i>0)
	      structure.add(ncoarse,nfine-(2*nx+1));
					     // upper
	    if (i<ny)
	      structure.add(ncoarse,nfine+(2*nx+1));
      }

  structure.compress();
  matrix.reinit(structure);
  
  for (unsigned int i=0 ; i<=ny; ++i)
    for (unsigned int j=0 ; j<=nx; ++j)
      {
					 // coarse grid point number
	unsigned int ncoarse =j+(nx+1)*i;
					 // same point on fine grid
	unsigned int nfine = 2*j+(4*nx+2)*i;
	
	matrix.set(ncoarse,nfine,1.);
					 // left
	if (j>0)
	  {
	    matrix.set(ncoarse,nfine-1,.5);
					     // lower left
	    if (i>0)
	      matrix.set(ncoarse,nfine-(2*nx+1)-1,.25);
					     // upper left
	    if (i<ny)
	      matrix.set(ncoarse,nfine+(2*nx+1)-1,.25);
	  }
					 // right
	if (j<nx)
	  {
	    matrix.set(ncoarse,nfine+1,.5);
					     // lower right
	    if (i>0)
	      matrix.set(ncoarse,nfine-(2*nx+1)+1,.25);
					     // upper right
	    if (i<ny)
	      matrix.set(ncoarse,nfine+(2*nx+1)+1,.25);
	  }

					     // lower
	    if (i>0)
	      matrix.set(ncoarse,nfine-(2*nx+1),.5);
					     // upper
	    if (i<ny)
	      matrix.set(ncoarse,nfine+(2*nx+1),.5);
      }
  
}

void
FDMGTransfer::prolongate (const unsigned int   to_level,
			  Vector<float>       &dst,
			  const Vector<float> &src) const
{
  Assert((to_level>0) && (to_level<=matrices.size()),
	 ExcIndexRange(to_level, 0, matrices.size()+1));
  
  matrices[to_level-1].Tvmult(dst,src);
}


void
FDMGTransfer::restrict (const unsigned int   from_level,
	 		Vector<float>       &dst,
			const Vector<float> &src) const
{
  Assert((from_level>0) && (from_level<=matrices.size()),
	 ExcIndexRange(from_level, 0, matrices.size()+1));

  matrices[from_level-1].vmult(dst,src);
}


template void FDMatrix::laplacian(SparseMatrix<double>&) const;
template void FDMatrix::laplacian(SparseMatrix<float>&) const;
