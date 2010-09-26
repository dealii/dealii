//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <base/utilities.h>
#include <fe/fe.h>
#include <fe/fe_tools.h>

DEAL_II_NAMESPACE_OPEN

namespace FETools
{

  template <int dim>
  void
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe,
					 std::vector<unsigned int> &h2l)
  {
    Assert (h2l.size() == fe.dofs_per_cell,
	    ExcDimensionMismatch (h2l.size(), fe.dofs_per_cell));
    h2l = hierarchic_to_lexicographic_numbering (fe);
  }



  template <int dim>
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering (const FiniteElementData<dim> &fe)
  {
    Assert (fe.n_components() == 1, ExcInvalidFE());

    std::vector<unsigned int> h2l (fe.dofs_per_cell);

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
	  h2l[next_index++] = n*(n-1);
	  h2l[next_index++] = n*n-1;

					   // left   line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (1+i)*n;

					   // right  line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (2+i)*n-1;

					   // bottom line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i;

					   // top    line
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*(n-1)+i+1;

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
	  h2l[next_index++] = 0;                 // 0
	  h2l[next_index++] = (      1)*degree;  // 1
	  h2l[next_index++] = (    n  )*degree;  // 2
	  h2l[next_index++] = (    n+1)*degree;  // 3
	  h2l[next_index++] = (n*n    )*degree;  // 4
	  h2l[next_index++] = (n*n  +1)*degree;  // 5
	  h2l[next_index++] = (n*n+n  )*degree;  // 6
	  h2l[next_index++] = (n*n+n+1)*degree;  // 7

					   // line 0
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n;
					   // line 1
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n;
					   // line 2
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i;
					   // line 3
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = 1+i+n*(n-1);

					   // line 4
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (n-1)*n*n+(i+1)*n;
					   // line 5
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (n-1)*(n*n+1)+(i+1)*n;
					   // line 6
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*n*(n-1)+i+1;
					   // line 7
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n*n*(n-1)+i+1+n*(n-1);

					   // line 8
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n*n;
					   // line 9
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n*n;
					   // line 10
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = (i+1)*n*n+n*(n-1);
					   // line 11
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    h2l[next_index++] = n-1+(i+1)*n*n+n*(n-1);


					   // inside quads
	  Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		  ExcInternalError());
					   // face 0
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (i+1)*n*n+n*(j+1);
					   // face 1
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (i+1)*n*n+n-1+n*(j+1);
					   // face 2, note the orientation!
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (j+1)*n*n+i+1;
					   // face 3, note the orientation!
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (j+1)*n*n+n*(n-1)+i+1;
					   // face 4
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = n*(i+1)+j+1;
					   // face 5
	  for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	    for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	      h2l[next_index++] = (n-1)*n*n+n*(i+1)+j+1;

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

    return h2l;
  }



  template <int dim>
  void
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe,
					 std::vector<unsigned int>    &l2h)
  {
    l2h = lexicographic_to_hierarchic_numbering (fe);
  }



  template <int dim>
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe)
  {
    return Utilities::invert_permutation(hierarchic_to_lexicographic_numbering (fe));
  }

}


// explicit instantiations
namespace FETools
{
  template
  void
  hierarchic_to_lexicographic_numbering<1>
  (const FiniteElementData<1> &fe,
   std::vector<unsigned int> &h2l);
  template
  void
  hierarchic_to_lexicographic_numbering<2>
  (const FiniteElementData<2> &fe,
   std::vector<unsigned int> &h2l);
  template
  void
  hierarchic_to_lexicographic_numbering<3>
  (const FiniteElementData<3> &fe,
   std::vector<unsigned int> &h2l);


  template
  void
  lexicographic_to_hierarchic_numbering<1>
  (const FiniteElementData<1> &fe,
   std::vector<unsigned int> &l2h);
  template
  void
  lexicographic_to_hierarchic_numbering<2>
  (const FiniteElementData<2> &fe,
   std::vector<unsigned int> &l2h);
  template
  void
  lexicographic_to_hierarchic_numbering<3>
  (const FiniteElementData<3> &fe,
   std::vector<unsigned int> &l2h);



  template
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering<1>
  (const FiniteElementData<1> &fe);
  template
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering<2>
  (const FiniteElementData<2> &fe);
  template
  std::vector<unsigned int>
  hierarchic_to_lexicographic_numbering<3>
  (const FiniteElementData<3> &fe);


  template
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering<1>
  (const FiniteElementData<1> &fe);
  template
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering<2>
  (const FiniteElementData<2> &fe);
  template
  std::vector<unsigned int>
  lexicographic_to_hierarchic_numbering<3>
  (const FiniteElementData<3> &fe);
}


DEAL_II_NAMESPACE_CLOSE
