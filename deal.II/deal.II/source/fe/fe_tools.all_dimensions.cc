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
FETools::lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe,
						std::vector<unsigned int>    &l2h)
{
				   // note: this function does the
				   // reverse operation of the
				   // previous one. nevertheless, they
				   // have been written independently
				   // from each other. the test
				   // "fe/numbering" checks that the
				   // output of the two functions is
				   // indeed the reverse of each other
				   // by checking that the
				   // concatenation of the two maps is
				   // the identity operation
				   //
				   // The experienced code reader will
				   // note that this function was not
				   // written by the same author than
				   // the previous one (although the
				   // author of the previous function
				   // cleaned up this if-block a
				   // little bit by introducing the
				   // arrays of numbers). Therefore,
				   // both authors have experienced
				   // the downsides of the hierarchic
				   // numbering of degrees of freedom
				   // in deal.II. Just to also provide
				   // some fun while reading code,
				   // here is the rant of the author
				   // of this function about the
				   // author of the previous one:
				   //
				   // "Unfortunately, somebody
				   // switched the upper corner points
				   // of a quad. The same person
				   // decided to find a very creative
				   // numbering of the vertices of a
				   // hexahedron. Therefore, this code
				   // looks quite sophisticated."
				   //
				   // NB: The "accused" same person
				   // claims to have had good reasons
				   // then, but seems to have
				   // forgotten about them. At least,
				   // the numbering was discussed with
				   // the complaining person back then
				   // when all began :-)
  Assert (fe.n_components() == 1, ExcInvalidFE());
  Assert (l2h.size() == fe.dofs_per_cell,
	  ExcDimensionMismatch (l2h.size(), fe.dofs_per_cell));
				   // polynomial degree
  const unsigned int degree = fe.dofs_per_line+1;
				   // number of grid points in each
				   // direction
  const unsigned int n = degree+1;

  if (degree > 0)
    {
      Assert (fe.dofs_per_vertex == 1, ExcInternalError());
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
	{
	  unsigned int index = 0;
					   // Find indices of vertices.
	  switch (dim)
	    {
	      case 1:
	      {
		const unsigned int values[GeometryInfo<1>::vertices_per_cell]
		  = { 0, degree };
		index = values[i];
		break;
	      };
	     
	      case 2:
	      {
		const unsigned int values[GeometryInfo<2>::vertices_per_cell]
		  = { 0, degree, n*degree+degree, n*degree };
		index = values[i];
		break;
	      };
	     
	      case 3:
	      {
		const unsigned int values[GeometryInfo<3>::vertices_per_cell]
		  = { 0, degree,
		      n*n*degree + degree, n*n*degree,
		      n*degree, n*degree+degree,
		      n*n*degree + n*degree+degree, n*n*degree + n*degree};
		index = values[i];
		break;
	      };
	     
	      default:
		    Assert(false, ExcNotImplemented());
	    }
	
	  l2h[index] = i;
	}
    };
  
				   // for degree 2 and higher: Lines,
				   // quads, hexes etc also carry
				   // degrees of freedom
  if (degree > 1)
    {
      Assert (fe.dofs_per_line == degree-1, ExcInternalError());
      Assert ((fe.dofs_per_quad == (degree-1)*(degree-1)) ||
	      (dim < 2), ExcInternalError());
      Assert ((fe.dofs_per_hex == (degree-1)*(degree-1)*(degree-1)) ||
	      (dim < 3), ExcInternalError());
	    
      for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::lines_per_cell); ++i)
	{
	  unsigned int index = fe.first_line_index + i*fe.dofs_per_line;
	  unsigned int incr = 0;
	  unsigned int tensorstart = 0;
					   // This again looks quite
					   // strange because of the odd
					   // numbering scheme.
	  switch (i+100*dim)
	    {
					       // lines in x-direction
	      case 100:
	      case 200: case 202:
	      case 300: case 302: case 304: case 306:
		    incr = 1;
		    break;
						     // lines in y-direction
	      case 201: case 203:
	      case 308: case 309: case 310: case 311:
		    incr = n;
		    break;
						     // lines in z-direction
	      case 301: case 303: case 305: case 307:
		    incr = n*n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());
	    }
	  switch (i+100*dim)
	    {
					       // x=y=z=0
	      case 100:
	      case 200: case 203:
	      case 300: case 303: case 308:
		    tensorstart = 0;
		    break;
						     // x=1 y=z=0
	      case 201:
	      case 301: case 309:
		    tensorstart = degree;
		    break;
						     // y=1 x=z=0
	      case 202:
	      case 304: case 307:
		    tensorstart = n*degree;
		    break;
						     // x=z=1 y=0
	      case 310:
		    tensorstart = n*n*degree+degree;
		    break;
						     // z=1 x=y=0
	      case 302: case 311:
		    tensorstart = n*n*degree;
		    break;
						     // x=y=1 z=0
	      case 305:
		    tensorstart = n*degree+degree;
		    break;
						     // y=z=1 x=0
	      case 306:
		    tensorstart = n*n*n-n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());	      
	    }
	  
	  for (unsigned int jx = 1; jx<degree ;++jx)
	    {
	      unsigned int tensorindex = tensorstart + jx * incr;
	      l2h[tensorindex] = index++;
	    }
	}

      for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::quads_per_cell); ++i)
	{
	  unsigned int index = fe.first_quad_index+i*fe.dofs_per_quad;
	  unsigned int tensorstart = 0;
	  unsigned int incx = 0;
	  unsigned int incy = 0;
	  switch (i)
	    {
	      case 0:
		    tensorstart = 0; incx = 1;
		    if (dim==2)
		      incy = n;
		    else
		      incy = n*n;
		    break;
	      case 1:
		    tensorstart = n*degree; incx = 1; incy = n*n;
		    break;
	      case 2:
		    tensorstart = 0; incx = 1; incy = n;
		    break;
	      case 3:
		    tensorstart = degree; incx = n; incy = n*n;
		    break;
	      case 4:
		    tensorstart = n*n*degree; incx = 1; incy = n;
		    break;
	      case 5:
		    tensorstart = 0; incx = n; incy = n*n;
		    break;
	      default:
		    Assert(false, ExcNotImplemented());	      
	    }
	  
	  for (unsigned int jy = 1; jy<degree; jy++)
	    for (unsigned int jx = 1; jx<degree ;++jx)
	      {
		unsigned int tensorindex = tensorstart
					   + jx * incx + jy * incy;
		l2h[tensorindex] = index++;
	      }
	}

      if (GeometryInfo<dim>::hexes_per_cell > 0)
	for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::hexes_per_cell); ++i)
	  {
	    unsigned int index = fe.first_hex_index;
	    
	    for (unsigned int jz = 1; jz<degree; jz++)
	      for (unsigned int jy = 1; jy<degree; jy++)
		for (unsigned int jx = 1; jx<degree; jx++)
		  {
		    const unsigned int tensorindex = jx + jy*n + jz*n*n;
		    l2h[tensorindex]=index++;
		  }  
	  } 
    }
}


template
void
FETools::lexicographic_to_hierarchic_numbering<1>
(const FiniteElementData<1> &fe,
 std::vector<unsigned int> &h2l);
template
void
FETools::lexicographic_to_hierarchic_numbering<2>
(const FiniteElementData<2> &fe,
 std::vector<unsigned int> &h2l);
template
void
FETools::lexicographic_to_hierarchic_numbering<3>
(const FiniteElementData<3> &fe,
 std::vector<unsigned int> &h2l);
