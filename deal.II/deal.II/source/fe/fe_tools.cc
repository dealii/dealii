//----------------------------  fe_tools.cc  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  fe_tools.cc  ---------------------------


#include <base/quadrature.h>
#include <lac/full_matrix.h>
#include <lac/vector.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <fe/fe_tools.h>
#include <fe/fe.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>
#include <fe/mapping_q1.h>
#include <dofs/dof_handler.h>
#include <dofs/dof_accessor.h>

// if necessary try to work around a bug in the IBM xlC compiler
#ifdef XLC_WORK_AROUND_STD_BUG
using namespace std;
#endif



template <int dim, typename number>
void FETools::get_interpolation_matrix(const FiniteElement<dim> &fe1,
				       const FiniteElement<dim> &fe2,
				       FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe2.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell,
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe2.dofs_per_cell,
				    fe1.dofs_per_cell));

				   // Initialize FEValues for fe1 at
				   // the unit support points of the
				   // fe2 element.
  const typename std::vector<Point<dim> > &
    fe2_support_points = fe2.get_unit_support_points ();

  Assert(fe2_support_points.size()==fe2.dofs_per_cell,
	 typename FiniteElementBase<dim>::ExcFEHasNoSupportPoints());

  for (unsigned int i=0; i<fe2.dofs_per_cell; ++i)	
    {
      const unsigned int i1 = fe2.system_to_component_index(i).first;
      for (unsigned int j=0; j<fe1.dofs_per_cell; ++j)
	{
	  const unsigned int j1 = fe1.system_to_component_index(j).first;
	  if (i1==j1)
	    interpolation_matrix(i,j) = fe1.shape_value (j,fe2_support_points[i]);
	  else
	    interpolation_matrix(i,j)=0.;
	}  
    }
}



template <int dim, typename number>
void FETools::get_back_interpolation_matrix(const FiniteElement<dim> &fe1,
					    const FiniteElement<dim> &fe2,
					    FullMatrix<number> &interpolation_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(interpolation_matrix.m()==fe1.dofs_per_cell &&
	 interpolation_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(interpolation_matrix.m(),
				    interpolation_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
    
  FullMatrix<number> first_matrix (fe2.dofs_per_cell, fe1.dofs_per_cell);
  FullMatrix<number> second_matrix(fe1.dofs_per_cell, fe2.dofs_per_cell);
  
  get_interpolation_matrix(fe1, fe2, first_matrix);
  get_interpolation_matrix(fe2, fe1, second_matrix);

				   // int_matrix=second_matrix*first_matrix
  second_matrix.mmult(interpolation_matrix, first_matrix);
}



template <int dim, typename number>
void FETools::get_interpolation_difference_matrix(const FiniteElement<dim> &fe1,
						  const FiniteElement<dim> &fe2,
						  FullMatrix<number> &difference_matrix)
{
  Assert (fe1.n_components() == fe2.n_components(),
	  ExcDimensionMismatch(fe1.n_components(), fe2.n_components()));
  Assert(difference_matrix.m()==fe1.dofs_per_cell &&
	 difference_matrix.n()==fe1.dofs_per_cell, 
	 ExcMatrixDimensionMismatch(difference_matrix.m(),
				    difference_matrix.n(),
				    fe1.dofs_per_cell,
				    fe1.dofs_per_cell));
   
  FullMatrix<number> interpolation_matrix(fe1.dofs_per_cell);
  get_back_interpolation_matrix(fe1, fe2, interpolation_matrix);
  
  for (unsigned int i=0; i<fe1.dofs_per_cell; ++i)
    difference_matrix(i,i) = 1.;
  
				   // compute difference
  difference_matrix.add (-1, interpolation_matrix);
}



template <int dim, typename number>
void FETools::interpolate(const DoFHandler<dim> &dof1,
			  const Vector<number> &u1,
			  const DoFHandler<dim> &dof2,
			  Vector<number> &u2)
{
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;
  const unsigned int dofs_per_cell2=dof2.get_fe().dofs_per_cell;
  
  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u2_local(dofs_per_cell2);

  FullMatrix<number> interpolation_matrix(dofs_per_cell2,
					  dofs_per_cell1);
  FETools::get_interpolation_matrix(dof1.get_fe(), dof2.get_fe(),
				    interpolation_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end(),
						 cell2 = dof2.begin_active(),
						 endc2 = dof2.end();

  std::vector<unsigned int> index_multiplicity(dof2.n_dofs(),0);
  std::vector<unsigned int> dofs (dofs_per_cell2);
  u2.clear ();
  
  for (; cell1!=endc1, cell2!=endc2; ++cell1, ++cell2) 
    {
      cell1->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u2_local, u1_local);
      cell2->get_dof_indices(dofs);
      for (unsigned int i=0; i<dofs_per_cell2; ++i)
	{
	  u2(dofs[i])+=u2_local(i);
	  ++index_multiplicity[dofs[i]];
	}
    }
  for (unsigned int i=0; i<dof2.n_dofs(); ++i)
    {
      Assert(index_multiplicity[i]!=0, ExcInternalError());
      u2(i) /= index_multiplicity[i];
    }
}



template <int dim, typename number>
void FETools::back_interpolate(const DoFHandler<dim> &dof1,
			       const Vector<number> &u1,
			       const FiniteElement<dim> &fe2,
			       Vector<number> &u1_interpolated)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_interpolated.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_interpolated.size(), dof1.n_dofs()));  

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u1_int_local(dofs_per_cell1);
  
  FullMatrix<number> interpolation_matrix(dofs_per_cell1, dofs_per_cell1);
  FETools::get_back_interpolation_matrix(dof1.get_fe(), fe2,
					 interpolation_matrix);

  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end();
  
  for (; cell1!=endc1; ++cell1) 
    {
      cell1->get_dof_values(u1, u1_local);
      interpolation_matrix.vmult(u1_int_local, u1_local);
      cell1->set_dof_values(u1_int_local, u1_interpolated);
    }
}



template <int dim, typename number>
void FETools::interpolation_difference(const DoFHandler<dim> &dof1,
				       const Vector<number> &u1,
				       const FiniteElement<dim> &fe2,
				       Vector<number> &u1_difference)
{
  Assert(dof1.get_fe().n_components() == fe2.n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), fe2.n_components()));
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));  
  Assert(u1_difference.size()==dof1.n_dofs(),
	 ExcDimensionMismatch(u1_difference.size(), dof1.n_dofs()));  

  const unsigned int dofs_per_cell1=dof1.get_fe().dofs_per_cell;

  Vector<number> u1_local(dofs_per_cell1);
  Vector<number> u1_diff_local(dofs_per_cell1);
  
  FullMatrix<number> difference_matrix(dofs_per_cell1, dofs_per_cell1);
  FETools::get_interpolation_difference_matrix(dof1.get_fe(), fe2,
					       difference_matrix);
  
  typename DoFHandler<dim>::active_cell_iterator cell1 = dof1.begin_active(),
						 endc1 = dof1.end();
  
  for (; cell1!=endc1; ++cell1) 
    {
      cell1->get_dof_values(u1, u1_local);
      difference_matrix.vmult(u1_diff_local, u1_local);
      cell1->set_dof_values(u1_diff_local, u1_difference);
    }
}



template <int dim, typename number>
void FETools::extrapolate(const DoFHandler<dim> &dof1,
			  const Vector<number> &u1,
			  const DoFHandler<dim> &dof2,
			  Vector<number> &u2)
{
  Assert(dof1.get_fe().n_components() == dof2.get_fe().n_components(),
	 ExcDimensionMismatch(dof1.get_fe().n_components(), dof2.get_fe().n_components()));
  Assert(&dof1.get_tria()==&dof2.get_tria(), ExcTriangulationMismatch());
  Assert(u1.size()==dof1.n_dofs(), ExcDimensionMismatch(u1.size(), dof1.n_dofs()));
  Assert(u2.size()==dof2.n_dofs(), ExcDimensionMismatch(u2.size(), dof2.n_dofs()));

  Vector<number> u3(dof2.n_dofs());
  interpolate(dof1, u1, dof2, u3);

  const unsigned int dofs_per_cell  = dof2.get_fe().dofs_per_cell;
  const Triangulation<dim> &tria=dof1.get_tria();
  Vector<number> dof_values(dofs_per_cell);
  
  for (unsigned int level=0; level<tria.n_levels()-1; ++level)
    {
      typename DoFHandler<dim>::cell_iterator cell=dof2.begin(level),
					      endc=dof2.end(level);

      for (; cell!=endc; ++cell)
	if (!cell->active())
	  {
					     // check whether this
					     // cell has active
					     // children
	    bool active_children=false;
	    for (unsigned int child_n=0;
		 child_n<GeometryInfo<dim>::children_per_cell; ++child_n)
	      if (cell->child(child_n)->active())
		{
		  active_children=true;
		  break;
		}

					     // if there are active
					     // children, the we have
					     // to work on this
					     // cell. get the data
					     // from the one vector
					     // and set it on the
					     // other
	    if (active_children)
	      {
		cell->get_interpolated_dof_values(u3, dof_values);
		cell->set_dof_values_by_interpolation(dof_values, u2);
	      }
	  }
    }
}



template <int dim>
void
FETools::hierarchic_to_lexicographic_numbering (const FE_Q<dim>           &fe,
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
      };

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
      };

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
      };       

      default:
	    Assert (false, ExcNotImplemented());
    };
};



template <int dim>
void
FETools::lexicographic_to_hierarchic_numbering (const FE_Q<dim>           &fe,
						std::vector<unsigned int> &l2h)
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
};




/*-------------- Explicit Instantiations -------------------------------*/

template
void FETools::get_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
				       const FiniteElement<deal_II_dimension> &,
				       FullMatrix<double> &);
template
void FETools::get_back_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
					    const FiniteElement<deal_II_dimension> &,
					    FullMatrix<double> &);
template
void FETools::get_interpolation_difference_matrix(const FiniteElement<deal_II_dimension> &,
						  const FiniteElement<deal_II_dimension> &,
						  FullMatrix<double> &);
template
void FETools::interpolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<double> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<double> &);
template
void FETools::back_interpolate(const DoFHandler<deal_II_dimension> &,
			       const Vector<double> &,
			       const FiniteElement<deal_II_dimension> &,
			       Vector<double> &);
template
void FETools::interpolation_difference(const DoFHandler<deal_II_dimension> &,
				       const Vector<double> &,
				       const FiniteElement<deal_II_dimension> &,
				       Vector<double> &);
template
void FETools::extrapolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<double> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<double> &);


template
void FETools::get_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
				       const FiniteElement<deal_II_dimension> &,
				       FullMatrix<float> &);
template
void FETools::get_back_interpolation_matrix(const FiniteElement<deal_II_dimension> &,
					    const FiniteElement<deal_II_dimension> &,
					    FullMatrix<float> &);
template
void FETools::get_interpolation_difference_matrix(const FiniteElement<deal_II_dimension> &,
						  const FiniteElement<deal_II_dimension> &,
						  FullMatrix<float> &);
template
void FETools::interpolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<float> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<float> &);
template
void FETools::back_interpolate(const DoFHandler<deal_II_dimension> &,
			       const Vector<float> &,
			       const FiniteElement<deal_II_dimension> &,
			       Vector<float> &);
template
void FETools::interpolation_difference(const DoFHandler<deal_II_dimension> &,
				       const Vector<float> &,
				       const FiniteElement<deal_II_dimension> &,
				       Vector<float> &);
template
void FETools::extrapolate(const DoFHandler<deal_II_dimension> &,
			  const Vector<float> &,
			  const DoFHandler<deal_II_dimension> &,
			  Vector<float> &);
template
void
FETools::hierarchic_to_lexicographic_numbering (const FE_Q<deal_II_dimension> &fe,
						std::vector<unsigned int>     &h2l);
template
void
FETools::lexicographic_to_hierarchic_numbering (const FE_Q<deal_II_dimension> &fe,
						std::vector<unsigned int>     &h2l);


/*----------------------------   fe_tools.cc     ---------------------------*/
