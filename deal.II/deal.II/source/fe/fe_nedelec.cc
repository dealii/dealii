//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <base/quadrature.h>
#include <base/qprojector.h>
#include <base/table.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_nedelec.h>
#include <fe/fe_values.h>

#include <sstream>

//TODO: implement the adjust_line_dof_index_for_line_orientation_table field,
//and write a test similar to bits/face_orientation_and_fe_q_02

DEAL_II_NAMESPACE_OPEN


//TODO: Remove doubled member variable 'degree'
template <int dim, int spacedim>
FE_Nedelec<dim,spacedim>::FE_Nedelec (const unsigned int degree)
		:
		FiniteElement<dim,spacedim> (
		  FiniteElementData<dim>(get_dpo_vector(degree), dim,
					 degree+1, FiniteElementData<dim>::Hcurl, 1),
				    std::vector<bool> (
				      FiniteElementData<dim>(get_dpo_vector(degree), dim,
							     degree+1).dofs_per_cell,false),
				    std::vector<std::vector<bool> >(
				      FiniteElementData<dim>(get_dpo_vector(degree), dim,
							     degree+1).dofs_per_cell,
				      std::vector<bool>(dim,true))),
		degree(degree)
{
  Assert (dim >= 2, ExcImpossibleInDim(dim));
  
				   // copy constraint and embedding
				   // matrices if they are
				   // defined. otherwise leave them at
				   // invalid size
  initialize_constraints ();
  initialize_embedding ();
  initialize_restriction ();

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();

				   // finite element classes need to
				   // initialize the
				   // adjust_quad_dof_index... table. however,
				   // for the current element, there are no
				   // dofs on quads in 3d (i.e. in the
				   // interior of a face), so there is nothing
				   // to do
  if (dim == 3)
    {
      Assert (this->dofs_per_quad == 0,
	      ExcInternalError());
      Assert (this->adjust_quad_dof_index_for_face_orientation_table.size(0)==
	      this->dofs_per_quad,
	      ExcInternalError());
    }
}



template <int dim, int spacedim>
std::string
FE_Nedelec<dim,spacedim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

  std::ostringstream namebuf;  
  namebuf << "FE_Nedelec<" << dim << ">(" << degree << ")";

  return namebuf.str();
}



template <int dim, int spacedim>
FiniteElement<dim,spacedim> *
FE_Nedelec<dim,spacedim>::clone() const
{
  return new FE_Nedelec<dim,spacedim>(*this);
}


#if deal_II_dimension == 1

template <>
double
FE_Nedelec<1>::shape_value_component (const unsigned int ,
				     const Point<1>    &,
				     const unsigned int ) const
{
  Assert (false, ExcNotImplemented());
  return 0.;
}

#endif

#if deal_II_dimension == 2

template <>
double
FE_Nedelec<2>::shape_value_component (const unsigned int i,
				      const Point<2>    &p,
				      const unsigned int component) const
{
  const unsigned int dim = 2;
  
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));
  
  switch (degree)
    {
				       // first order Nedelec elements
      case 1:
      {
	switch (i)
	  {
						   // (0, 1-x)
	    case 0: return (component == 0 ? 0 : 1-p(0));
						   // (0,x)
	    case 1: return (component == 0 ? 0 : p(0));
						   // (1-y, 0)
	    case 2: return (component == 0 ? 1-p(1) : 0);
						   // (y, 0)
	    case 3: return (component == 0 ? p(1) : 0);
                        
						   // there are only
						   // four shape
						   // functions!?
	    default:
		  Assert (false, ExcInternalError());
		  return 0;
	  };
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return 0;
}

#endif

#if deal_II_dimension == 3

template <>
double
FE_Nedelec<3>::shape_value_component (const unsigned int i,
				      const Point<3>    &p,
				      const unsigned int component) const
{
  const unsigned int dim = 3;
  
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));
  
  switch (degree)
    {
				       // first order Nedelec
				       // elements
      case 1:
      {
					 // note that the degrees of
					 // freedom on opposite faces
					 // have a common vector
					 // direction, so simplify that
					 // a little. these directions
					 // are:
					 //
					 // for lines 2, 3, 6, 7:
					 //    (1,0,0)
					 // for lines 0, 1, 4, 5:
					 //    (0,1,0)
					 // for lines 8, 9, 10, 11:
					 //    (0,0,1)
					 //
					 // thus, sort out all those
					 // cases where the component
					 // is zero anyway, and only
					 // otherwise compute the
					 // spatially dependent part
					 // which is then also the
					 // return value
	if (((i<8) && (((i%4>=2) && (component!=0)) ||
		       ((i%4<2) && (component!=1)))) ||
	    ((i>=8) && (component != 2)))
	  return 0;

					       // now we know that the
					       // only non-zero
					       // component is
					       // requested:
	const double x = p(0),
		     y = p(1),
		     z = p(2);
	switch (i)
	  {
	    case  0: return (1-x)*(1-z);
	    case  1: return     x*(1-z);
	    case  2: return (1-y)*(1-z);
	    case  3: return     y*(1-z);

	    case  4: return (1-x)*z;
	    case  5: return     x*z;
	    case  6: return (1-y)*z;
	    case  7: return     y*z;

	    case  8: return (1-x)*(1-y);
	    case  9: return     x*(1-y);
	    case 10: return (1-x)*y;
	    case 11: return     x*y;
			
	    default:
		  Assert (false, ExcInternalError());
		  return 0;
	  };
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return 0;
}

#endif

#if deal_II_dimension == 1

template <>
Tensor<1,1>
FE_Nedelec<1>::shape_grad_component (const unsigned int ,
				     const Point<1>    &,
				     const unsigned int ) const
{
  Assert (false, ExcNotImplemented());
  return Tensor<1,1>();
}

#endif

#if deal_II_dimension == 2

template <>
Tensor<1,2>
FE_Nedelec<2>::shape_grad_component (const unsigned int i,
				     const Point<2>    &,
				     const unsigned int component) const
{
  const unsigned int dim = 2;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

  switch (degree)
    {
				       // first order Nedelec elements
      case 1:
      {
					 // on the unit cell, the
					 // gradients of these shape
					 // functions are constant, so
					 // we pack them into a table
					 // for simpler lookup
					 //
					 // the format is: first
					 // index=shape function
					 // number; second
					 // index=vector component,
					 // third index=component
					 // within gradient
	static const double unit_gradients[4][2][2]
	  = { { {0., 0.}, {-1.,0.} },
	      { {0., 0.}, {+1.,0.} },
	      { {0.,-1.}, { 0.,0.} },
	      { {0.,+1.}, { 0.,0.} } };
	return Tensor<1,dim>(unit_gradients[i][component]);
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return Tensor<1,dim>();
}

#endif

#if deal_II_dimension == 3

template <>
Tensor<1,3>
FE_Nedelec<3>::shape_grad_component (const unsigned int i,
				     const Point<3>    &p,
				     const unsigned int component) const
{
  const unsigned int dim = 3;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

  switch (degree)
    {
				       // first order Nedelec elements
      case 1:
      {
					 // on the unit cell, the
					 // gradients of these shape
					 // functions are linear. we
					 // pack them into an array,
					 // knowing that it may be
					 // expensive to recompute the
					 // whole array each
					 // time. maybe some clever
					 // compiler can optimize this
					 // out, seeing that except
					 // for one element all the
					 // other ones are dead
					 // stores...
					 //
					 // the format is: first
					 // index=shape function
					 // number; second
					 // index=vector component,
					 // third index=component
					 // within gradient
	const double x = p(0),
		     y = p(1),
		     z = p(2);
	const double unit_gradients[12][3][3]
	  = { { {0,      0,      0}, {-(1-z), 0, -(1-x)}, {0, 0, 0} },
	      { {0,      0,      0}, { (1-z), 0,     -x}, {0, 0, 0} },
	      { {0, -(1-z), -(1-y)}, {0,      0,      0}, {0, 0, 0} },
	      { {0,  (1-z),     -y}, {0,      0,      0}, {0, 0, 0} },
	      
	      { {0,  0,     0}, {-z, 0, (1-x)}, {0, 0, 0} },
	      { {0,  0,     0}, { z, 0,     x}, {0, 0, 0} },
	      { {0, -z, (1-y)}, { 0, 0,     0}, {0, 0, 0} },
	      { {0,  z,     y}, { 0, 0,     0}, {0, 0, 0} },

	      { {0, 0, 0}, {0, 0, 0}, {-(1-y), -(1-x), 0} },
	      { {0, 0, 0}, {0, 0, 0}, { (1-y),     -x, 0} },
	      { {0, 0, 0}, {0, 0, 0}, {    -y,  (1-x), 0} },
	      { {0, 0, 0}, {0, 0, 0}, {     y,      x, 0} } };
					 // note: simple check whether
					 // this can at all be: build
					 // the sum over all these
					 // tensors. since the sum of
					 // the shape functions is a
					 // constant, the gradient
					 // must necessarily be
					 // zero. this is in fact the
					 // case here, so test
					 // successfull
	return Tensor<1,dim>(unit_gradients[i][component]);
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };
  
  return Tensor<1,dim>();
}

#endif


#if deal_II_dimension == 1

template <>
Tensor<2,1>
FE_Nedelec<1>::shape_grad_grad_component (const unsigned int ,
					  const Point<1>    &,
					  const unsigned int ) const
{
  Assert (false, ExcNotImplemented());
  return Tensor<2,1>();
}

#endif


#if deal_II_dimension == 2

template <>
Tensor<2,2>
FE_Nedelec<2>::shape_grad_grad_component (const unsigned int i,
					  const Point<2> &/*p*/,
					  const unsigned int component) const
{
  const unsigned int dim = 2;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

  switch (degree)
    {
				       // first order Nedelec
				       // elements. their second
				       // derivatives on the unit cell
				       // are zero
      case 1:
      {
	return Tensor<2,dim>();
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };

  return Tensor<2,dim>();
}

#endif

#if deal_II_dimension == 3

template <>
Tensor<2,3>
FE_Nedelec<3>::shape_grad_grad_component (const unsigned int i,
					  const Point<3>    &/*p*/,
					  const unsigned int component) const
{
  const unsigned int dim = 3;
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component < dim, ExcIndexRange (component, 0, dim));

  switch (degree)
    {
				       // first order Nedelec
				       // elements. their second
				       // derivatives on the unit cell
				       // are constant, but non-zero
      case 1:
      {
					 // the format is: first
					 // index=shape function
					 // number; second
					 // index=vector component,
					 // third and fourth
					 // index=component within
					 // second derivative
	static const double unit_grad_grads[12][3][3][3]
	  = {
		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 1}, {0, 0, 0}, {1, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0,-1}, {0, 0, 0}, {-1, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },
		
		{ { {0, 0, 0}, {0, 0, 1}, {0, 1, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0,-1}, {0,-1, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0,-1}, {0, 0, 0}, {-1, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 1}, {0, 0, 0}, {1, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0,-1}, {0,-1, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 1}, {0, 1, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 1, 0}, {1, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0,-1, 0}, {-1, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0,-1, 0}, {-1, 0, 0}, {0, 0, 0} } },

		{ { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 0, 0}, {0, 0, 0}, {0, 0, 0} },
		  { {0, 1, 0}, {1, 0, 0}, {0, 0, 0} } }
	  };

	return Tensor<2,dim>(unit_grad_grads[i][component]);
      };

					// no other degrees
					// implemented
      default:
	    Assert (false, ExcNotImplemented());
    };

  return Tensor<2,dim>();
}

#endif

//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------



template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::initialize_constraints ()
{
				   // copy constraint matrices if they
				   // are defined. otherwise leave
				   // them at zero size
  if (degree<Matrices::n_constraint_matrices+1)
    {
      this->interface_constraints.
        TableBase<2,double>::reinit (this->interface_constraints_size());
      this->interface_constraints.fill (Matrices::constraint_matrices[degree-1]);
    };
}



template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::initialize_embedding ()
{
  unsigned int iso=RefinementCase<dim>::isotropic_refinement-1;
  if ((degree < Matrices::n_embedding_matrices+1) &&
      (Matrices::embedding[degree-1][0] != 0))
    for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
      {
                                         // copy
        this->prolongation[iso][c].reinit (this->dofs_per_cell,
                                      this->dofs_per_cell);
        this->prolongation[iso][c].fill (Matrices::embedding[degree-1][c]);
                                         // and make sure that the row
                                         // sum is 0.5 (for usual
                                         // elements, the row sum must
                                         // be 1, but here the shape
                                         // function is multiplied by
                                         // the inverse of the
                                         // Jacobian, which introduces
                                         // a factor of 1/2 when going
                                         // from mother to child)
        for (unsigned int row=0; row<this->dofs_per_cell; ++row)
          {
            double sum = 0;
            for (unsigned int col=0; col<this->dofs_per_cell; ++col)
              sum += this->prolongation[iso][c](row,col);
            Assert (std::fabs(sum-.5) < 1e-14,
                    ExcInternalError());
          };
      };
}



template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::initialize_restriction ()
{
  unsigned int iso=RefinementCase<dim>::isotropic_refinement-1;
  switch (dim)
    {
      case 2:   // 2d
      {
	switch (degree)
	  {
	    case 1:
	    {
                                               // this is a strange
                                               // element, since it is
                                               // both additive and
                                               // then it is also
                                               // not. ideally, we
                                               // would like to have
                                               // the value of the
                                               // shape function on
                                               // the coarse line to
                                               // be the mean value of
                                               // that on the two
                                               // child ones. thus,
                                               // one should make it
                                               // additive. however,
                                               // additivity only
                                               // works if an element
                                               // does not have any
                                               // continuity
                                               // requirements, since
                                               // otherwise degrees of
                                               // freedom are shared
                                               // between adjacent
                                               // elements, and when
                                               // we make the element
                                               // additive, that would
                                               // mean that we end up
                                               // adding up
                                               // contributions not
                                               // only from the child
                                               // cells of this cell,
                                               // but also from the
                                               // child cells of the
                                               // neighbor, and since
                                               // we cannot know
                                               // whether there even
                                               // exists a neighbor we
                                               // cannot simply make
                                               // the element
                                               // additive.
					       //
                                               // so, until someone
                                               // comes along with a
                                               // better alternative,
                                               // we do the following:
                                               // make the element
                                               // non-additive, and
                                               // simply pick the
                                               // value of one of the
                                               // child lines for the
                                               // value of the mother
                                               // line (note that we
                                               // have to multiply by
                                               // two, since the shape
                                               // functions scale with
                                               // the inverse
                                               // Jacobian). we thus
                                               // throw away the
                                               // information of one
                                               // of the child lines,
                                               // but there seems to
                                               // be no other way than
                                               // that...
                                               //
                                               // note: to make things
                                               // consistent, and
                                               // restriction
                                               // independent of the
                                               // order in which we
                                               // travel across the
                                               // cells of the coarse
                                               // grid, we have to
                                               // make sure that we
                                               // take the same small
                                               // line when visiting
                                               // its two neighbors,
                                               // to get the value for
                                               // the mother line. we
                                               // take the first line
                                               // always, in the
                                               // canonical direction
                                               // of lines
              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                this->restriction[iso][c].reinit (this->dofs_per_cell,
						  this->dofs_per_cell);
	      
	      this->restriction[iso][0](0,0) = 2.;
	      this->restriction[iso][1](1,1) = 2.;
	      this->restriction[iso][0](2,2) = 2.;
	      this->restriction[iso][2](3,3) = 2.;

	      break;
	    };
	    
	    default:
	    {
					       // in case we don't
					       // have the matrices
					       // (yet), leave them
					       // empty. this does not
					       // prevent the use of
					       // this FE, but will
					       // prevent the use of
					       // these matrices
              break;
	    };
	  };
	
	break;
      };


      case 3:   // 3d
      {
	switch (degree)
	  {
	    case 1:
	    {
					       // same principle as in
					       // 2d, take one child
					       // cell to get at the
					       // values of each of
					       // the 12 lines
              for (unsigned int c=0; c<GeometryInfo<dim>::max_children_per_cell; ++c)
                this->restriction[iso][c].reinit (this->dofs_per_cell,
						  this->dofs_per_cell);
	      this->restriction[iso][0](0,0) = 2.;
	      this->restriction[iso][1](1,1) = 2.;
	      this->restriction[iso][0](2,2) = 2.;
	      this->restriction[iso][2](3,3) = 2.;
              
	      this->restriction[iso][4](4,4) = 2.;
	      this->restriction[iso][5](5,5) = 2.;
	      this->restriction[iso][4](6,6) = 2.;
	      this->restriction[iso][6](7,7) = 2.;
              
	      this->restriction[iso][0](8,8) = 2.;
	      this->restriction[iso][1](9,9) = 2.;
	      this->restriction[iso][2](10,10) = 2.;
	      this->restriction[iso][3](11,11) = 2.;
              
	      break;
	    };
	    
	    default:
	    {
					       // in case we don't
					       // have the matrices
					       // (yet), leave them
					       // empty. this does not
					       // prevent the use of
					       // this FE, but will
					       // prevent the use of
					       // these matrices
              break;
	    };
	  };
	
	break;
      };
      
      default:
	    Assert (false,ExcNotImplemented());
    }
}



template <int dim, int spacedim>
void FE_Nedelec<dim,spacedim>::initialize_unit_support_points ()
{
  switch (degree)
    {
      case 1:
      {
					 // all degrees of freedom are
					 // on edges, and their order
					 // is the same as the edges
					 // themselves
	this->unit_support_points.resize(GeometryInfo<dim>::lines_per_cell);
	for (unsigned int line=0; line<GeometryInfo<dim>::lines_per_cell; ++line)
	  {
	    const unsigned int
	      vertex_index_0 = GeometryInfo<dim>::line_to_cell_vertices(line,0),
	      vertex_index_1 = GeometryInfo<dim>::line_to_cell_vertices(line,1);
	    
	    const Point<dim>
	      vertex_0 = GeometryInfo<dim>::unit_cell_vertex(vertex_index_0),
	      vertex_1 = GeometryInfo<dim>::unit_cell_vertex(vertex_index_1);
	    
					     // place dofs right
					     // between the vertices
					     // of each line
	    this->unit_support_points[line] = (vertex_0 + vertex_1) / 2;
	  };
	    
	break;
      };

      default:
					     // no higher order
					     // elements implemented
					     // right now
	    Assert (false, ExcNotImplemented());
    };
}


#if deal_II_dimension == 1

template <>
void FE_Nedelec<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

#endif


template <int dim, int spacedim>
void FE_Nedelec<dim,spacedim>::initialize_unit_face_support_points ()
{
  switch (degree)
    {
      case 1:
      {
					 // do this the same as above, but
					 // for one dimension less
	this->unit_face_support_points.resize(GeometryInfo<dim-1>::lines_per_cell);
	for (unsigned int line=0; line<GeometryInfo<dim-1>::lines_per_cell; ++line)
	  {
	    const unsigned int
	      vertex_index_0 = GeometryInfo<dim-1>::line_to_cell_vertices(line,0),
	      vertex_index_1 = GeometryInfo<dim-1>::line_to_cell_vertices(line,1);
      
	    const Point<dim-1>
	      vertex_0 = GeometryInfo<dim-1>::unit_cell_vertex(vertex_index_0),
	      vertex_1 = GeometryInfo<dim-1>::unit_cell_vertex(vertex_index_1);

					     // place dofs right
					     // between the vertices of each
					     // line
	      this->unit_face_support_points[line] = (vertex_0 + vertex_1) / 2;
	  };
	break;
      };

      default:
					     // no higher order
					     // elements implemented
					     // right now
	    Assert (false, ExcNotImplemented());
    };	    
}



template <int dim, int spacedim>
std::vector<unsigned int>
FE_Nedelec<dim,spacedim>::get_dpo_vector(const unsigned int degree)
{
  Assert (degree == 1, ExcNotImplemented());

				   // for degree==1, put all degrees
				   // of freedom on the lines, and in
				   // particular @p{degree} DoFs per
				   // line:
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[1] = degree;

  return dpo;
}



template <int dim, int spacedim>
UpdateFlags
FE_Nedelec<dim,spacedim>::update_once (const UpdateFlags) const
{
				   // even the values have to be
				   // computed on the real cell, so
				   // nothing can be done in advance
  return update_default;
}



template <int dim, int spacedim>
UpdateFlags
FE_Nedelec<dim,spacedim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values             | update_covariant_transformation;
  if (flags & update_gradients)
    out |= update_gradients          | update_covariant_transformation;
  if (flags & update_hessians)
    out |= update_hessians | update_covariant_transformation;

  return out;
}



//---------------------------------------------------------------------------
// Data field initialization
//---------------------------------------------------------------------------

template <int dim, int spacedim>
typename Mapping<dim,spacedim>::InternalDataBase *
FE_Nedelec<dim,spacedim>::get_data (const UpdateFlags      update_flags,
			   const Mapping<dim,spacedim>    &mapping,
			   const Quadrature<dim> &quadrature) const
{
 				   // generate a new data object and
 				   // initialize some fields
   InternalData* data = new InternalData;

 				   // check what needs to be
 				   // initialized only once and what
 				   // on every cell/face/subface we
 				   // visit
   data->update_once = update_once(update_flags);
   data->update_each = update_each(update_flags);
   data->update_flags = data->update_once | data->update_each;

   const UpdateFlags flags(data->update_flags);
   const unsigned int n_q_points = quadrature.size();

 				   // initialize fields only if really
 				   // necessary. otherwise, don't
 				   // allocate memory
   if (flags & update_values)
     data->shape_values.resize (this->dofs_per_cell,
                                std::vector<Tensor<1,dim> > (n_q_points));

   if (flags & update_gradients)
     data->shape_gradients.resize (this->dofs_per_cell,
                                   std::vector<Tensor<2,dim> > (n_q_points));

 				   // if second derivatives through
 				   // finite differencing is required,
 				   // then initialize some objects for
 				   // that
   if (flags & update_hessians)
     data->initialize_2nd (this, mapping, quadrature);

 				   // next already fill those fields
 				   // of which we have information by
 				   // now. note that the shape values
 				   // and gradients are only those on
 				   // the unit cell, and need to be
 				   // transformed when visiting an
 				   // actual cell
   for (unsigned int i=0; i<this->dofs_per_cell; ++i)
     for (unsigned int q=0; q<n_q_points; ++q)
       {
	 if (flags & update_values)
	   for (unsigned int c=0; c<dim; ++c)
	     data->shape_values[i][q][c]
	       = shape_value_component(i,quadrature.point(q),c);
	
	 if (flags & update_gradients)
	   for (unsigned int c=0; c<dim; ++c)
	     data->shape_gradients[i][q][c]
	       = shape_grad_component(i,quadrature.point(q),c);
       }
   
   return data;
}




//---------------------------------------------------------------------------
// Fill data of FEValues
//---------------------------------------------------------------------------

template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::fill_fe_values 
  (const Mapping<dim,spacedim>                      &mapping,
   const typename Triangulation<dim,spacedim>::cell_iterator &cell,
   const Quadrature<dim>                            &quadrature,
   typename Mapping<dim,spacedim>::InternalDataBase &mapping_data,
   typename Mapping<dim,spacedim>::InternalDataBase &fedata,
   FEValuesData<dim,spacedim>                       &data,
   enum CellSimilarity::Similarity                  &/*cell_similarity*/) const
{
 				   // convert data object to internal
 				   // data for this class. fails with
 				   // an exception if that is not
 				   // possible
  Assert (dynamic_cast<InternalData *> (&fedata) != 0,
	  ExcInternalError());
  InternalData &fe_data = static_cast<InternalData &> (fedata);

				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.size();
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  Assert (fe_data.shape_values[k].size() == n_q_points,
		  ExcInternalError());
	  mapping.transform(fe_data.shape_values[k], shape_values,
			    mapping_data, mapping_covariant);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
	  Assert (fe_data.shape_gradients[k].size() == n_q_points,
		  ExcInternalError());
                                           // do first transformation
	  mapping.transform(fe_data.shape_gradients[k], shape_grads1,
			    mapping_data, mapping_covariant);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform(shape_grads2, shape_grads1,
			    mapping_data, mapping_covariant);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_hessians)
    this->compute_2nd (mapping, cell,
                       QProjector<dim>::DataSetDescriptor::cell(),
                       mapping_data, fe_data, data);
}



template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::fill_fe_face_values (const Mapping<dim,spacedim>                   &mapping,
				      const typename Triangulation<dim,spacedim>::cell_iterator &cell,
				      const unsigned int                    face,
				      const Quadrature<dim-1>              &quadrature,
				      typename Mapping<dim,spacedim>::InternalDataBase       &mapping_data,
				      typename Mapping<dim,spacedim>::InternalDataBase       &fedata,
				      FEValuesData<dim,spacedim>                    &data) const
{
 				   // convert data object to internal
 				   // data for this class. fails with
 				   // an exception if that is not
 				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

                                   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset
    = (QProjector<dim>::DataSetDescriptor::
       face (face,
	     cell->face_orientation(face),
	     cell->face_flip(face),
	     cell->face_rotation(face),
             quadrature.size()));

  				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.size();
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
                                       // check size of array. in 3d,
                                       // we have faces oriented both
                                       // ways
      Assert (fe_data.shape_values[0].size() ==
              GeometryInfo<dim>::faces_per_cell * n_q_points *
              (dim == 3 ? 8 : 1),
              ExcInternalError());
      
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  mapping.transform(make_slice(fe_data.shape_values[k], offset, n_q_points),
			    shape_values, mapping_data, mapping_covariant);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
                                       // check size of array. in 3d,
                                       // we have faces oriented both
                                       // ways
      Assert (fe_data.shape_gradients[0].size() ==
              GeometryInfo<dim>::faces_per_cell * n_q_points *
              (dim == 3 ? 8 : 1),
              ExcInternalError());

      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
                                           // 
                                           // do first transformation
	  mapping.transform(make_slice(fe_data.shape_gradients[k], offset, n_q_points),
			    shape_grads1, mapping_data, mapping_covariant);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform(shape_grads2, shape_grads1,
			    mapping_data, mapping_covariant);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_hessians)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim, int spacedim>
void
FE_Nedelec<dim,spacedim>::fill_fe_subface_values (const Mapping<dim,spacedim>                   &mapping,
					 const typename Triangulation<dim,spacedim>::cell_iterator &cell,
					 const unsigned int                    face,
					 const unsigned int                    subface,
					 const Quadrature<dim-1>              &quadrature,
					 typename Mapping<dim,spacedim>::InternalDataBase       &mapping_data,
					 typename Mapping<dim,spacedim>::InternalDataBase       &fedata,
					 FEValuesData<dim,spacedim>                    &data) const
{
 				   // convert data object to internal
 				   // data for this class. fails with
 				   // an exception if that is not
 				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

                                   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  const typename QProjector<dim>::DataSetDescriptor offset
    = (QProjector<dim>::DataSetDescriptor::
       subface (face, subface,
		cell->face_orientation(face),
		cell->face_flip(face),
		cell->face_rotation(face),
		quadrature.size(),
		cell->subface_case(face)));

  				   // get the flags indicating the
				   // fields that have to be filled
  const UpdateFlags flags(fe_data.current_update_flags());

  const unsigned int n_q_points = quadrature.size();
				  
				   // fill shape function
				   // values. these are vector-valued,
				   // so we have to transform
				   // them. since the output format
				   // (in data.shape_values) is a
				   // sequence of doubles (one for
				   // each non-zero shape function
				   // value, and for each quadrature
				   // point, rather than a sequence of
				   // small vectors, we have to use a
				   // number of conversions
  if (flags & update_values)
    {
      Assert (fe_data.shape_values[0].size() ==
 	      GeometryInfo<dim>::max_children_per_face *
              GeometryInfo<dim>::faces_per_cell *
	      n_q_points,
             ExcInternalError());
      
      std::vector<Tensor<1,dim> > shape_values (n_q_points);

      Assert (data.shape_values.n_rows() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_values.n_cols() == n_q_points,
	      ExcInternalError());
      
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
					   // first transform shape
					   // values...
	  mapping.transform(make_slice(fe_data.shape_values[k], offset, n_q_points),
			    shape_values, mapping_data, mapping_covariant);

					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_values[k*dim+d][q] = shape_values[q][d];
	};
    };
  
      
  if (flags & update_gradients)
    {
      Assert (fe_data.shape_gradients.size() ==
              GeometryInfo<dim>::faces_per_cell *
	      GeometryInfo<dim>::max_children_per_face *
	      n_q_points,
              ExcInternalError());

      std::vector<Tensor<2,dim> > shape_grads1 (n_q_points);
      std::vector<Tensor<2,dim> > shape_grads2 (n_q_points);

      Assert (data.shape_gradients.size() == this->dofs_per_cell * dim,
	      ExcInternalError());
      Assert (data.shape_gradients[0].size() == n_q_points,
	      ExcInternalError());

                                       // loop over all shape
                                       // functions, and treat the
                                       // gradients of each shape
                                       // function at all quadrature
                                       // points
      for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	{
                                           // treat the gradients of
                                           // this particular shape
                                           // function at all
                                           // q-points. if Dv is the
                                           // gradient of the shape
                                           // function on the unit
                                           // cell, then
                                           // (J^-T)Dv(J^-1) is the
                                           // value we want to have on
                                           // the real cell. so, we
                                           // will have to apply a
                                           // covariant transformation
                                           // to Dv twice. since the
                                           // interface only allows
                                           // multiplication with
                                           // (J^-1) from the right,
                                           // we have to trick a
                                           // little in between
                                           // 
                                           // do first transformation
	  mapping.transform(make_slice(fe_data.shape_gradients[k], offset, n_q_points),
			    shape_grads1, mapping_data, mapping_covariant);
                                           // transpose matrix
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
                                           // do second transformation
	  mapping.transform(shape_grads2, shape_grads1,
			    mapping_data, mapping_covariant);
                                           // transpose back
          for (unsigned int q=0; q<n_q_points; ++q)
            shape_grads2[q] = transpose(shape_grads1[q]);
          
					   // then copy over to target:
	  for (unsigned int q=0; q<n_q_points; ++q)
	    for (unsigned int d=0; d<dim; ++d)
	      data.shape_gradients[k*dim+d][q] = shape_grads2[q][d];
	};
    }

  if (flags & update_hessians)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim, int spacedim>
unsigned int
FE_Nedelec<dim,spacedim>::n_base_elements () const
{
  return 1;
}



template <int dim, int spacedim>
const FiniteElement<dim,spacedim> &
FE_Nedelec<dim,spacedim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
}



template <int dim, int spacedim>
unsigned int
FE_Nedelec<dim,spacedim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
}



template <int dim, int spacedim>
bool
FE_Nedelec<dim,spacedim>::has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));

  switch (degree)
    {
      case 1:
      {
        switch (dim)
          {
            case 2:
            {
                                               // only on the one
                                               // non-adjacent face
                                               // are the values
                                               // actually zero. list
                                               // these in a table
              const unsigned int
                opposite_faces[GeometryInfo<2>::faces_per_cell]
                = { 1, 0, 3, 2};
              
              return (face_index != opposite_faces[shape_index]);
            };
            
            case 3:
            {
                                               // the shape functions
                                               // are zero on the two
                                               // faces opposite the
                                               // two faces adjacent
                                               // to the line the
                                               // shape function is
                                               // defined on
              const unsigned int
                opposite_faces[GeometryInfo<3>::lines_per_cell][2]
                = { {1,5}, {0,5}, {3,5}, {2,5},
		    {1,4}, {0,4}, {3,4}, {2,4},
		    {1,3}, {0,3}, {1,2}, {0,2}};
              
              return ((face_index != opposite_faces[shape_index][0])
                      &&
                      (face_index != opposite_faces[shape_index][1]));
            };
            
            default: Assert (false, ExcNotImplemented());
          };
      };
      
      default:  // other degree
            Assert (false, ExcNotImplemented());
    };
  
  return true;
}



template <int dim, int spacedim>
unsigned int
FE_Nedelec<dim,spacedim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim, int spacedim>
unsigned int
FE_Nedelec<dim,spacedim>::get_degree () const
{
  return degree;
}


template class FE_Nedelec<deal_II_dimension>;

DEAL_II_NAMESPACE_CLOSE

