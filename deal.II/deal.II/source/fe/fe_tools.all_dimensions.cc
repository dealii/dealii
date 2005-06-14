//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <fe/fe.h>
#include <fe/fe_tools.h>



template <int dim>
void
FETools::hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe,
						std::vector<unsigned int> &h2l)
{
  Assert (fe.n_components() == 1, ExcInvalidFE());
  Assert (h2l.size() == fe.dofs_per_cell,
	  ExcDimensionMismatch (h2l.size(), fe.dofs_per_cell));

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
	h2l[0] = 0;
	h2l[1] = dofs_per_cell-1;
	for (unsigned int i=2; i<dofs_per_cell; ++i)
	  h2l[i] = i-1;

	break;
      }

      case 2:
      {
	unsigned int next_index = 0;
					 // first the four vertices
	h2l[next_index++] = 0;
	h2l[next_index++] = n-1;
	h2l[next_index++] = n*n-1;
	h2l[next_index++] = n*(n-1);
					 // first line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = 1+i;
					 // second line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n-1;
					 // third line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n*(n-1)+i+1;
					 // fourth line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (1+i)*n;
					 // inside quad
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = n*(i+1)+j+1;

	Assert (next_index == fe.dofs_per_cell, ExcInternalError());

	break;
      }

      case 3:
      {
	unsigned int next_index = 0;
					 // first the eight vertices
	h2l[next_index++] = 0;
	h2l[next_index++] = n-1;
	h2l[next_index++] = (n-1)*(n*n+1);
	h2l[next_index++] = (n-1)*n*n;
	h2l[next_index++] = n*(n-1);
	h2l[next_index++] = n*n-1;
	h2l[next_index++] = n*n*n-1;
	h2l[next_index++] = (n-1)*(n*n+n);

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = 1+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n-1+(i+1)*n*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n*n*(n-1)+i+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (i+1)*n*n;

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = 1+i+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n-1+(i+1)*n*n+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n*n*(n-1)+i+1+n*(n-1);
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (i+1)*n*n+n*(n-1);

	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n-1+(i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (n-1)*(n*n+1)+(i+1)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (n-1)*n*n+(i+1)*n;

					 // inside quads
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
					 // quad 1
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (i+1)*n*n+j+1;
					 // quad 2
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (i+1)*n*n+n*(n-1)+j+1;
					 // quad 3
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = n*(i+1)+j+1;
					 // quad 4
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (i+1)*n*n+n-1+n*(j+1);
					 // quad 5
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (n-1)*n*n+n*(i+1)+j+1;
					 // quad 6
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (i+1)*n*n+n*(j+1);

					 // inside hex
	Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    for (unsigned int k=0; k<fe.dofs_per_line; ++k)
	      h2l[next_index++]	= n*n*(i+1)+n*(j+1)+k+1;

	Assert (next_index == fe.dofs_per_cell, ExcInternalError());
	
	break;
      }       

      default:
	    Assert (false, ExcNotImplemented());
    }
}




template <int dim>
void
FETools::lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe,
						std::vector<unsigned int>    &l2h)
{
  std::vector<unsigned int> tmp(l2h.size());
  FETools::hierarchic_to_lexicographic_numbering(fe,tmp);
  for (unsigned int i=0; i<l2h.size(); ++i)
    l2h[tmp[i]]=i;  
}




template
void
FETools::hierarchic_to_lexicographic_numbering<1>
(const FiniteElementData<1> &fe,
 std::vector<unsigned int> &h2l);
template
void
FETools::hierarchic_to_lexicographic_numbering<2>
(const FiniteElementData<2> &fe,
 std::vector<unsigned int> &h2l);
template
void
FETools::hierarchic_to_lexicographic_numbering<3>
(const FiniteElementData<3> &fe,
 std::vector<unsigned int> &h2l);


template
void
FETools::lexicographic_to_hierarchic_numbering<1>
(const FiniteElementData<1> &fe,
 std::vector<unsigned int> &l2h);
template
void
FETools::lexicographic_to_hierarchic_numbering<2>
(const FiniteElementData<2> &fe,
 std::vector<unsigned int> &l2h);
template
void
FETools::lexicographic_to_hierarchic_numbering<3>
(const FiniteElementData<3> &fe,
 std::vector<unsigned int> &l2h);
