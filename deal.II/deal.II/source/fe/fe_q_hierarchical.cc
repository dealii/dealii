//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------------------------------------------

#include <base/quadrature.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <grid/tria.h>
#include <grid/tria_iterator.h>
#include <dofs/dof_accessor.h>
#include <fe/fe.h>
#include <fe/mapping.h>
#include <fe/fe_q_hierarchical.h>
#include <fe/fe_values.h>

#include <cmath>



template <int dim>
FE_Q_Hierarchical<dim>::FE_Q_Hierarchical (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1),
				    std::vector<bool> (FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
                                                       false),
				    std::vector<std::vector<bool> >(FiniteElementData<dim>(get_dpo_vector(degree),1).dofs_per_cell,
								    std::vector<bool>(1,true))),
		degree(degree),
		renumber(this->dofs_per_cell, 0),
		renumber_inverse(this->dofs_per_cell, 0),
	        face_renumber(this->dofs_per_face, 0),
		polynomial_space(Polynomials::Hierarchical<double>::generate_complete_basis(degree))
{
				   // do some internal book-keeping on
				   // cells and faces. if in 1d, the
				   // face function is empty
  lexicographic_to_hierarchic_numbering (*this, degree, renumber);
  face_lexicographic_to_hierarchic_numbering (degree, face_renumber);
  
				   // build inverse of renumbering
				   // vector
  for (unsigned int i=0; i<this->dofs_per_cell; ++i)
    renumber_inverse[renumber[i]]=i;

                                  // build the dofs_subcell, dofs_cell 
                                  // matrices for use in building prolongation,                                   // restriction, and constraint matrices.
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;

  for (unsigned int c=0; c<GeometryInfo<1>::children_per_cell; ++c)
  {
    dofs_cell.push_back (FullMatrix<double> (dofs_1d,dofs_1d) );
    dofs_subcell.push_back (FullMatrix<double> (dofs_1d,dofs_1d) );

    for (unsigned int j=0; j<dofs_1d; ++j)
    {
      for (unsigned int k=0; k<dofs_1d; ++k)
      {	
	                          // upper diagonal block
	if ((j<=1) && (k<=1))
	{
	  if (((c==0) && (j==0) && (k==0)) || 
	      ((c==1) && (j==1) && (k==1)))
	    dofs_cell[c](j,k) = 1.;
	  else 
	    dofs_cell[c](j,k) = 0.;
	  
	  if      (((c==0) && (j==1)) || ((c==1) && (j==0)))
	    dofs_subcell[c](j,k) = .5;
	  else if (((c==0) && (k==0)) || ((c==1) && (k==1)))
	    dofs_subcell[c](j,k) = 1.;
	  else
	    dofs_subcell[c](j,k) = 0.;
	}
	                         // upper right block
	else if ((j<=1) && (k>=2))
	{
	  if (((c==0) && (j==1) && ((k % 2)==0)) ||
	      ((c==1) && (j==0) && ((k % 2)==0)))
	    dofs_subcell[c](j,k) = -1.;
	}
	                        // lower diagonal block
	else if ((j>=2) && (k>=2) && (j<=k))
	{
	  double factor = 1.;
	  for (unsigned int i=1; i<=j;++i)
	    factor *= ((double) (k-i+1))/((double) i);
	  if (c==0)
	  {
	    dofs_subcell[c](j,k) = ((k+j) % 2 == 0) ? 
                                   std::pow(.5,static_cast<double>(k))*factor :
                                   -std::pow(.5,static_cast<double>(k))*factor;
	    dofs_cell[c](j,k) = std::pow(2.,static_cast<double>(j))*factor;
	  }
	  else
	  {
	    dofs_subcell[c](j,k) = std::pow(.5,static_cast<double>(k))*factor;
	    dofs_cell[c](j,k) = ((k+j) % 2 == 0) ? 
                                std::pow(2.,static_cast<double>(j))*factor :
                                -std::pow(2.,static_cast<double>(j))*factor;
	  }
	}
      }
    }
  }
				   // fill constraint matrices
  if (dim==2 || dim==3)
  {
    this->interface_constraints.reinit ( (dim==2) ? 1 + 2*(degree-1) : 
                         5 + 12*(degree-1) + 4*(degree-1)*(degree-1),
					 (dim==2) ? (degree+1) : 
					 (degree+1)*(degree+1) );
    switch (dim)
    {
      case 2:
	// vertex node
	for (unsigned int i=0; i<dofs_1d; ++i)
	  this->interface_constraints(0,i) = dofs_subcell[0](1,i); 
	// edge nodes
	for (unsigned int c=0; c<GeometryInfo<1>::children_per_cell; ++c)
	  for (unsigned int i=0; i<dofs_1d; ++i)
	    for (unsigned int j=2; j<dofs_1d; ++j)
	      this->interface_constraints(1 + c*(degree-1) + j - 2,i) = 
                                       dofs_subcell[c](j,i);
	break;
      case 3:
	for (unsigned int i=0; i<dofs_1d * dofs_1d; i++)
	{
	  // center vertex node	  
	  this->interface_constraints(0,face_renumber[i]) = 
	    dofs_subcell[0](1,i % dofs_1d) * 
	    dofs_subcell[0](1,(i - (i % dofs_1d)) / dofs_1d);

	  // boundary vertex nodes
	  this->interface_constraints(1,face_renumber[i]) = 
	    dofs_subcell[0](1, i % dofs_1d) * 
	    dofs_subcell[0](0, (i - (i % dofs_1d)) / dofs_1d);
	  this->interface_constraints(2,face_renumber[i]) = 
	    dofs_subcell[1](1, i % dofs_1d) * 
	    dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
	  this->interface_constraints(3,face_renumber[i]) = 
	    dofs_subcell[1](0, i % dofs_1d) * 
	    dofs_subcell[1](1, (i - (i % dofs_1d)) / dofs_1d);
	  this->interface_constraints(4,face_renumber[i]) = 
	    dofs_subcell[0](0, i % dofs_1d) * 
	    dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
	  
	  // interior edges
	  for (unsigned int j=0; j<degree-1; j++)
	  {
	    this->interface_constraints(5 + j,face_renumber[i]) = 
	      dofs_subcell[0](1, i % dofs_1d) * 
	      dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + (degree-1) + j,face_renumber[i]) = 
	      dofs_subcell[1](2 + j, i % dofs_1d) * 
	      dofs_subcell[0](1, (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 2*(degree-1) + j,face_renumber[i]) = 
	      dofs_subcell[0](1,i % dofs_1d) * 
	      dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 3*(degree-1) + j,face_renumber[i]) = 
	      dofs_subcell[0](2 + j,i % dofs_1d) * 
	      dofs_subcell[1](0, (i - (i % dofs_1d)) / dofs_1d);
	  }

	  // boundary edges
	  for (unsigned int j=0; j<degree-1; j++)
	  {
	    // bottom edge 
	    this->interface_constraints(5 + 4*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[0](2 + j, i % dofs_1d) * 
	      dofs_subcell[0](0,     (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 4*(degree-1) + (degree-1) + j,face_renumber[i]) =
	      dofs_subcell[1](2 + j, i % dofs_1d) * 
	      dofs_subcell[0](0,     (i - (i % dofs_1d)) / dofs_1d);
	    // right edge
	    this->interface_constraints(5 + 4*(degree-1) + 2*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[1](1,     i % dofs_1d) * 
	      dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 4*(degree-1) + 3*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[1](1,     i % dofs_1d) * 
	      dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	    // top edge
	    this->interface_constraints(5 + 4*(degree-1) + 4*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[0](2 + j, i % dofs_1d) * 
	      dofs_subcell[1](1,     (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 4*(degree-1) + 5*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[1](2 + j, i % dofs_1d) * 
	      dofs_subcell[1](1,     (i - (i % dofs_1d)) / dofs_1d);
	    // left edge
	    this->interface_constraints(5 + 4*(degree-1) + 6*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[0](0,     i % dofs_1d) * 
	      dofs_subcell[0](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	    this->interface_constraints(5 + 4*(degree-1) + 7*(degree-1) + j,face_renumber[i]) =
	      dofs_subcell[0](0,     i % dofs_1d) * 
	      dofs_subcell[1](2 + j, (i - (i % dofs_1d)) / dofs_1d);
	  }

	  // interior faces
	  for (unsigned int j=0; j<degree-1; j++)
	    for (unsigned int k=0; k<degree-1; k++)
	    {
	      // subcell 0
	      this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1),face_renumber[i]) =
		dofs_subcell[0](2 + j, i % dofs_1d) * 
		dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
	      // subcell 1
	      this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + (degree-1)*(degree-1),face_renumber[i]) =
		dofs_subcell[1](2 + j, i % dofs_1d) * 
		dofs_subcell[0](2 + k, (i - (i % dofs_1d)) / dofs_1d);
	      // subcell 2
	      this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + 2*(degree-1)*(degree-1),face_renumber[i]) =
		dofs_subcell[1](2 + j, i % dofs_1d) * 
		dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
	      // subcell 3
	      this->interface_constraints(5 + 12*(degree-1) + j + k*(degree-1) + 3*(degree-1)*(degree-1),face_renumber[i]) =
		dofs_subcell[0](2 + j, i % dofs_1d) * 
		dofs_subcell[1](2 + k, (i - (i % dofs_1d)) / dofs_1d);
	    }
	}
	break;
    }
  };

				   // fill prolongation and restriction 
                                   // matrices
  if (dim==1)
  {
    for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      this->prolongation[c].reinit (this->dofs_per_cell,this->dofs_per_cell);
      this->prolongation[c].fill (dofs_subcell[c]);

      this->restriction[c].reinit (this->dofs_per_cell,this->dofs_per_cell);
      this->restriction[c].fill (dofs_cell[c]);
    }
  }
  else if (dim==2 || dim==3)
  {
    for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      this->prolongation[c].reinit (this->dofs_per_cell,this->dofs_per_cell);
      this->restriction[c].reinit (this->dofs_per_cell,this->dofs_per_cell);
    }
                                       // j loops over dofs in the subcell. 
                                       // These are the rows in the 
                                       // embedding matrix.
    for (unsigned int j=0; j<this->dofs_per_cell; ++j)
    {
                                       // i loops over the dofs in the master 
                                       // cell. These are the columns in 
                                       // the embedding matrix.
      for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      {
	switch (dim)
	{
	  case 2:
	    for (unsigned int c=0; c<GeometryInfo<2>::children_per_cell; ++c)
	    {
	      unsigned int c0 = ((c==1) || (c==2)) ? 1 : 0;
	      unsigned int c1 = ((c==2) || (c==3)) ? 1 : 0;

	      this->prolongation[c](j,i) = 
		dofs_subcell[c0](renumber_inverse[j] % dofs_1d,
				 renumber_inverse[i] % dofs_1d) *
		dofs_subcell[c1]((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d,
				 (renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d);

	      this->restriction[c](j,i) = 
		dofs_cell[c0](renumber_inverse[j] % dofs_1d,
			      renumber_inverse[i] % dofs_1d) *
		dofs_cell[c1]((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d,
			      (renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d);
	    }
	    break;

	  case 3:
	    for (unsigned int c=0; c<GeometryInfo<3>::children_per_cell; ++c)
	    {
	      unsigned int c0 = ((c==1) || (c==2) || (c==5) || (c==6)) ? 1 : 0;
	      unsigned int c1 = ((c==4) || (c==5) || (c==6) || (c==7)) ? 1 : 0;
              unsigned int c2 = ((c==2) || (c==3) || (c==6) || (c==7)) ? 1 : 0;

	      this->prolongation[c](j,i) = 
		dofs_subcell[c0](renumber_inverse[j] % dofs_1d,
				 renumber_inverse[i] % dofs_1d) *
		dofs_subcell[c1](((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d) % dofs_1d,
				 ((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
		dofs_subcell[c2](((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d - (((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
				 ((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d - (((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);

	      this->restriction[c](j,i) = 
		dofs_cell[c0](renumber_inverse[j] % dofs_1d,
			      renumber_inverse[i] % dofs_1d) *
		dofs_cell[c1](((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d) % dofs_1d,
			      ((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
		dofs_cell[c2](((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d - (((renumber_inverse[j] - (renumber_inverse[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
			      ((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d - (((renumber_inverse[i] - (renumber_inverse[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);
	    }
	    break;
	}
      }
    }
  }
  else
    Assert (false, ExcNotImplemented());
				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
};



template <int dim>
FiniteElement<dim> *
FE_Q_Hierarchical<dim>::clone() const
{
  return new FE_Q_Hierarchical<dim>(degree);
}



template <int dim>
double
FE_Q_Hierarchical<dim>::shape_value (const unsigned int i,
  			             const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_value(renumber_inverse[i], p);
}



template <int dim>
double
FE_Q_Hierarchical<dim>::shape_value_component (const unsigned int i,
					       const Point<dim> &p,
					       const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_value(renumber_inverse[i], p);
}



template <int dim>
Tensor<1,dim>
FE_Q_Hierarchical<dim>::shape_grad (const unsigned int i,
				    const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<1,dim>
FE_Q_Hierarchical<dim>::shape_grad_component (const unsigned int i,
					      const Point<dim> &p,
					      const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<2,dim>
FE_Q_Hierarchical<dim>::shape_grad_grad (const unsigned int i,
			      const Point<dim> &p) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  return polynomial_space.compute_grad_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<2,dim>
FE_Q_Hierarchical<dim>::shape_grad_grad_component (const unsigned int i,
					const Point<dim> &p,
					const unsigned int component) const
{
  Assert (i<this->dofs_per_cell, ExcIndexRange(i,0,this->dofs_per_cell));
  Assert (component == 0, ExcIndexRange (component, 0, 1));
  return polynomial_space.compute_grad_grad(renumber_inverse[i], p);
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------



template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_support_points ()
{
				   // number of points: (degree+1)^dim
  unsigned int n = degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= degree+1;
  
  this->unit_support_points.resize(n);
  
  Point<dim> p;
                                   // the method of numbering allows
                                   // each dof to be associated with a
                                   // support point. There is
                                   // only one support point per
                                   // vertex, line, quad, hex, etc.
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((dim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=degree; ++ix)
      {
	if (ix==0)
	  p(0) =  0.;
	else if (ix==1)
	  p(0) =  1.;
	else
	  p(0) = .5;
	if (dim>1)
	{
	  if (iy==0)
	    p(1) =  0.;
	  else if (iy==1)
	    p(1) =  1.;
	  else
	    p(1) = .5;
	}
	if (dim>2)
	{
	  if (iz==0)
	    p(2) =  0.;
	  else if (iz==1)
	    p(2) =  1.;
	  else
	    p(2) = .5;
	}
	this->unit_support_points[renumber[k++]] = p;
      };
};


#if deal_II_dimension == 1

template <>
void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
};

#endif


template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_face_support_points ()
{
  const unsigned int codim = dim-1;
  
				   // number of points: (degree+1)^codim
  unsigned int n = degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= degree+1;
  
  this->unit_face_support_points.resize(n);
  
  Point<codim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=degree; ++ix)
      {
	if (ix==0)
	  p(0) =  0.;
	else if (ix==1)
	  p(0) =  1.;
	else
	  p(0) = .5;
	if (codim>1)
	{
	  if (iy==0)
	    p(1) =  0.;
	  else if (iy==1)
	    p(1) =  1.;
	  else
	    p(1) = .5;
	}
	if (codim>2)
	{
	  if (iz==0)
	    p(2) =  0.;
	  else if (iz==1)
	    p(2) =  1.;
	  else
	    p(2) = .5;
	}
	this->unit_face_support_points[face_renumber[k++]] = p;
      };
};


                                           // we use same dpo_vector as FE_Q
template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, static_cast<unsigned int>(1));
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim>
void
FE_Q_Hierarchical<dim>::lexicographic_to_hierarchic_numbering (
			     const FiniteElementData<dim> &fe_data,
                             const unsigned int            degree,
			     std::vector<unsigned int>    &renumber)
{
  const unsigned int n = degree+1;


  if (degree == 0)
  {
    Assert ( (fe_data.dofs_per_vertex == 0) &&
	     ((fe_data.dofs_per_line == 0) || (dim == 1)) &&
	     ((fe_data.dofs_per_quad == 0) || (dim == 2)) &&
	     ((fe_data.dofs_per_hex == 0)  || (dim == 3)),
	     ExcInternalError() );
    renumber[0] = 0;
  };

  if (degree > 0)
  {
    Assert (fe_data.dofs_per_vertex == 1, ExcInternalError());
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
    {
      unsigned int index = 0;
					   // Find indices of vertices.
					   // Unfortunately, somebody
					   // switched the upper corner
					   // points of a quad. The same
					   // person decided to find a very
					   // creative numbering of the
					   // vertices of a hexahedron.
					   // Therefore, this looks quite
					   // sophisticated.
					   //
					   // NB: This same person
					   // claims to have had good
					   // reasons then, but seems to
					   // have forgotten about
					   // them. At least, the
					   // numbering was discussed
					   // with the complaining
					   // person back then when all
					   // began :-)
      switch (dim)
      {
        case 1:
	{
	  const unsigned int values[GeometryInfo<1>::vertices_per_cell]
	    = { 0, 1 };
	  index = values[i];
	  break;
	};
	     
        case 2:
	{
	  const unsigned int values[GeometryInfo<2>::vertices_per_cell]
	    = { 0, 1, n + 1, n };
	  index = values[i];
	  break;
	};
	     
        case 3:
	{
	  const unsigned int values[GeometryInfo<3>::vertices_per_cell]
	    = { 0,             1,
		n * n + 1,     n * n,
		n,             n + 1,
		n * n + n + 1, n * n + n};
	  index = values[i];
	  break;
	};
	
        default:
	  Assert(false, ExcNotImplemented());
      }
      
      renumber[index] = i;
    }
  };
				   // for degree 2 and higher: Lines,
				   // quads, hexes etc also carry
				   // degrees of freedom
  if (degree > 1)
  {
    Assert (fe_data.dofs_per_line == degree-1, ExcInternalError());
    Assert ((fe_data.dofs_per_quad == (degree-1)*(degree-1)) ||
	    (dim < 2), ExcInternalError());
    Assert ((fe_data.dofs_per_hex == (degree-1)*(degree-1)*(degree-1)) ||
	    (dim < 3), ExcInternalError());
    
    for (unsigned int i=0; i<GeometryInfo<dim>::lines_per_cell; ++i)
    {
      unsigned int index = fe_data.first_line_index
	               + i*fe_data.dofs_per_line;
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
	  incr = n * n;
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
	  tensorstart = 1;
	  break;
						 // y=1 x=z=0
        case 202:
        case 304: case 307:
	  tensorstart = n;
	  break;
						 // x=z=1 y=0
        case 310:
	  tensorstart = n * n + 1;
	  break;
						 // z=1 x=y=0
        case 302: case 311:
	  tensorstart = n * n;
	  break;
						 // x=y=1 z=0
        case 305:
	  tensorstart = n + 1;
	  break;
						 // y=z=1 x=0
        case 306:
	  tensorstart = n * n + n;
	  break;
        default:
	  Assert(false, ExcNotImplemented());	      
      }
	  
      for (unsigned int jx = 2; jx<=degree ;++jx)
      {
	unsigned int tensorindex = tensorstart + jx * incr;
	renumber[tensorindex] = index++;
      }
    }
    
    for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::quads_per_cell); ++i)
    {
      unsigned int index = fe_data.first_quad_index+i*fe_data.dofs_per_quad;
      unsigned int tensorstart = 0;
      unsigned int incx = 0;
      unsigned int incy = 0;
      switch (i)
      {
    	                                         // z=0 (dim==2), y=0 (dim==3)
        case 0:
	  tensorstart = 0; incx = 1;
	  if (dim==2)
	    incy = n;
	  else
	    incy = n * n;
	  break;
 	                                        // y=1
        case 1:
	  tensorstart = n; incx = 1; incy = n * n;
	  break;
	                                        // z=0
        case 2:
	  tensorstart = 0; incx = 1; incy = n;
	  break;
	                                        // x=1
        case 3:
	  tensorstart = 1; incx = n; incy = n * n;
	  break;
	                                        // z=1
        case 4:
	  tensorstart = n * n; incx = 1; incy = n;
	  break;
	                                        // x=0
        case 5:
	  tensorstart = 0; incx = n; incy = n * n;
	  break;
        default:
	  Assert(false, ExcNotImplemented());	      
      }
	  
      for (unsigned int jy = 2; jy<=degree; jy++)
	for (unsigned int jx = 2; jx<=degree ;++jx)
	{
	  unsigned int tensorindex = tensorstart
	                           + jx * incx + jy * incy;
	  renumber[tensorindex] = index++;
	}
    }
    
    if (GeometryInfo<dim>::hexes_per_cell > 0)
      for (int i=0; i<static_cast<signed int>(GeometryInfo<dim>::hexes_per_cell); ++i)
      {
	unsigned int index = fe_data.first_hex_index;
	    
	for (unsigned int jz = 2; jz<=degree; jz++)
	  for (unsigned int jy = 2; jy<=degree; jy++)
	    for (unsigned int jx = 2; jx<=degree; jx++)
	    {
	      const unsigned int tensorindex = jx + jy * n + jz * n * n;
	      renumber[tensorindex]=index++;
	    }  
      } 
  }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::face_lexicographic_to_hierarchic_numbering (
		         const unsigned int         degree,
			 std::vector<unsigned int> &numbering)
{
  FiniteElementData<dim-1> fe_data(FE_Q_Hierarchical<dim-1>::get_dpo_vector(degree),1);
  FE_Q_Hierarchical<dim-1>::lexicographic_to_hierarchic_numbering (fe_data,
							degree,
							numbering); 
}


#if (deal_II_dimension == 1)

template <>
void
FE_Q_Hierarchical<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int,
						 std::vector<unsigned int>&)
{}

#endif


template <int dim>
UpdateFlags
FE_Q_Hierarchical<dim>::update_once (const UpdateFlags flags) const
{
				   // for this kind of elements, only
				   // the values can be precomputed
				   // once and for all. set this flag
				   // if the values are requested at
				   // all
  return (update_default | (flags & update_values));
}



template <int dim>
UpdateFlags
FE_Q_Hierarchical<dim>::update_each (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_gradients)
    out |= update_gradients | update_covariant_transformation;
  if (flags & update_second_derivatives)
    out |= update_second_derivatives | update_covariant_transformation;

  return out;
}



//----------------------------------------------------------------------
// Data field initialization
//----------------------------------------------------------------------

template <int dim>
typename Mapping<dim>::InternalDataBase *
FE_Q_Hierarchical<dim>::get_data (const UpdateFlags      update_flags,
				  const Mapping<dim>    &mapping,
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
  const unsigned int n_q_points = quadrature.n_quadrature_points;

				   // some scratch arrays
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);

				   // initialize fields only if really
				   // necessary. otherwise, don't
				   // allocate memory
  if (flags & update_values)
    {
      values.resize (this->dofs_per_cell);
      data->shape_values.reinit (this->dofs_per_cell,
				 n_q_points);
    }

  if (flags & update_gradients)
    {
      grads.resize (this->dofs_per_cell);
      data->shape_gradients.reinit (this->dofs_per_cell,
				    n_q_points);
    }

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);

				   // next already fill those fields
				   // of which we have information by
				   // now. note that the shape
				   // gradients are only those on the
				   // unit cell, and need to be
				   // transformed when visiting an
				   // actual cell
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
	polynomial_space.compute(quadrature.point(i),
				 values, grads, grad_grads);
	
	if (flags & update_values)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_values[renumber[k]][i] = values[k];
	
	if (flags & update_gradients)
	  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
	    data->shape_gradients[renumber[k]][i] = grads[k];
      }
  return data;
}




//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_Q_Hierarchical<dim>::fill_fe_values (
		 const Mapping<dim>                            &mapping,
		 const typename DoFHandler<dim>::cell_iterator &cell,
		 const Quadrature<dim>                         &quadrature,
		 typename Mapping<dim>::InternalDataBase       &mapping_data,
		 typename Mapping<dim>::InternalDataBase       &fedata,
		 FEValuesData<dim>                             &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
	for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());	  
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin(),
				      mapping_data);
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, 0, mapping_data, fe_data, data);
}



template <int dim>
void
FE_Q_Hierarchical<dim>::fill_fe_face_values (
		  const Mapping<dim>                            &mapping,
		  const typename DoFHandler<dim>::cell_iterator &cell,
		  const unsigned int                             face,
		  const Quadrature<dim-1>                       &quadrature,
	      	  typename Mapping<dim>::InternalDataBase       &mapping_data,
		  typename Mapping<dim>::InternalDataBase       &fedata,
		  FEValuesData<dim>                             &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // faces are stored contiguously)
  const unsigned int offset = face * quadrature.n_quadrature_points;
  
  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin()+offset,
				      mapping_data);
	};
    }

  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
void
FE_Q_Hierarchical<dim>::fill_fe_subface_values (
		   const Mapping<dim>                            &mapping,
	           const typename DoFHandler<dim>::cell_iterator &cell,
		   const unsigned int                             face,
		   const unsigned int                             subface,
		   const Quadrature<dim-1>                       &quadrature,
		   typename Mapping<dim>::InternalDataBase       &mapping_data,
		   typename Mapping<dim>::InternalDataBase       &fedata,
		   FEValuesData<dim>                             &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // sub-faces are stored contiguously)
  const unsigned int offset = ((face * GeometryInfo<dim>::subfaces_per_face + subface)
                               * quadrature.n_quadrature_points);

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<this->dofs_per_cell; ++k)
    {
      if (flags & update_values)
        for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	{
	  Assert (data.shape_gradients[k].size() + offset <=
		  fe_data.shape_gradients[k].size(),
		  ExcInternalError());
	  mapping.transform_covariant(data.shape_gradients[k].begin(),
				      data.shape_gradients[k].end(),
				      fe_data.shape_gradients[k].begin()+offset,
				      mapping_data);
	};
    }
  
  if (flags & update_second_derivatives)
    this->compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
}



template <int dim>
unsigned int
FE_Q_Hierarchical<dim>::n_base_elements () const
{
  return 1;
};



template <int dim>
const FiniteElement<dim> &
FE_Q_Hierarchical<dim>::base_element (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return *this;
};



template <int dim>
unsigned int
FE_Q_Hierarchical<dim>::element_multiplicity (const unsigned int index) const
{
  Assert (index==0, ExcIndexRange(index, 0, 1));
  return 1;
};



template <int dim>
bool
FE_Q_Hierarchical<dim>::has_support_on_face (const unsigned int shape_index,
				             const unsigned int face_index) const
{
  Assert (shape_index < this->dofs_per_cell,
	  ExcIndexRange (shape_index, 0, this->dofs_per_cell));
  Assert (face_index < GeometryInfo<dim>::faces_per_cell,
	  ExcIndexRange (face_index, 0, GeometryInfo<dim>::faces_per_cell));


				   // in 1d, things are simple. since
				   // there is only one degree of
				   // freedom per vertex in this
				   // class, the first is on vertex 0
				   // (==face 0 in some sense), the
				   // second on face 1:
  if (dim==1)
    return (((shape_index == 0) && (face_index == 0)) ||
	    ((shape_index == 1) && (face_index == 1)));
  else
				     // more dimensions
    {
                                       // first, special-case interior
                                       // shape functions, since they
                                       // have no support no-where on
                                       // the boundary
      if (((dim==2) && (shape_index>=this->first_quad_index))
          ||
          ((dim==3) && (shape_index>=this->first_hex_index)))
        return false;
                                       
                                       // let's see whether this is a
                                       // vertex
      if (shape_index < this->first_line_index) 
        {
                                           // for Q elements, there is
                                           // one dof per vertex, so
                                           // shape_index==vertex_number. check
                                           // whether this vertex is
                                           // on the given face. thus,
                                           // for each face, give a
                                           // list of vertices
          const unsigned int vertex_no = shape_index;
          Assert (vertex_no < GeometryInfo<dim>::vertices_per_cell,
                  ExcInternalError());
          switch (dim)
            {
              case 2:
              {
                static const unsigned int face_vertices[4][2] =
                  { {0,1},{1,2},{2,3},{0,3} };
                return ((face_vertices[face_index][0] == vertex_no)
                        ||
                        (face_vertices[face_index][1] == vertex_no));
              };

              case 3:
              {
                static const unsigned int face_vertices[6][4] =
                  { {0,1,2,3},{4,5,6,7},{0,1,5,4},
                    {1,5,6,2},{3,2,6,7},{0,4,7,3} };
                return ((face_vertices[face_index][0] == vertex_no)||
                        (face_vertices[face_index][1] == vertex_no)||
                        (face_vertices[face_index][2] == vertex_no)||
                        (face_vertices[face_index][3] == vertex_no));
              };

              default:
                    Assert (false, ExcNotImplemented());
            };
        }
      else if (shape_index < this->first_quad_index)
                                         // ok, dof is on a line
        {
          const unsigned int line_index
            = (shape_index - this->first_line_index) / this->dofs_per_line;
          Assert (line_index < GeometryInfo<dim>::lines_per_cell,
                  ExcInternalError());

                                           // in 2d, the line is the
                                           // face, so get the line
                                           // index
          if (dim == 2)
            return (line_index == face_index);
          else if (dim == 3)
            {
                                               // see whether the
                                               // given line is on the
                                               // given face. use
                                               // table technique
                                               // again
              static const unsigned int face_lines[6][4] =
                { {0,1,2,3},{4,5,6,7},{0,8,9,4},
                  {1,9,5,10},{2,10,6,11},{3,8,7,11} };
              return ((face_lines[face_index][0] == line_index)||
                      (face_lines[face_index][1] == line_index)||
                      (face_lines[face_index][2] == line_index)||
                      (face_lines[face_index][3] == line_index));
            }
          else
            Assert (false, ExcNotImplemented());
        }
      else if (shape_index < this->first_hex_index)
                                         // dof is on a quad
        {
          const unsigned int quad_index 
            = (shape_index - this->first_quad_index) / this->dofs_per_quad;
          Assert (static_cast<signed int>(quad_index) <
                  static_cast<signed int>(GeometryInfo<dim>::quads_per_cell),
                  ExcInternalError());
          
                                           // in 2d, cell bubble are
                                           // zero on all faces. but
                                           // we have treated this
                                           // case above already
          Assert (dim != 2, ExcInternalError());

                                           // in 3d,
                                           // quad_index=face_index
          if (dim == 3)
            return (quad_index == face_index);
          else
            Assert (false, ExcNotImplemented());
        }
      else
                                         // dof on hex
        {
                                           // can only happen in 3d,
                                           // but this case has
                                           // already been covered
                                           // above
          Assert (false, ExcNotImplemented());
          return false;
        };
    };

                                   // we should not have gotten here
  Assert (false, ExcInternalError());
  return false;

}



template <int dim>
std::vector<unsigned int> 
FE_Q_Hierarchical<dim>::get_embedding_dofs (const unsigned int sub_degree) const
{
  Assert ((sub_degree>0) && (sub_degree<=degree),
	  ExcIndexRange(sub_degree, 1, degree));

  if (dim==1)
  {
    std::vector<unsigned int> embedding_dofs (sub_degree+1);
    for (unsigned int i=0; i<(sub_degree+1); ++i)
      embedding_dofs[i] = i;

    return embedding_dofs;
  }
  
  if (sub_degree==1)
  {
    std::vector<unsigned int> embedding_dofs (GeometryInfo<dim>::vertices_per_cell);
    for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
      embedding_dofs[i] = i;

    return embedding_dofs;
  }
  else if (sub_degree==degree)
  {
    std::vector<unsigned int> embedding_dofs (this->dofs_per_cell);
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      embedding_dofs[i] = i;

    return embedding_dofs;
  }

  if ((dim==2) || (dim==3))
  {
    std::vector<unsigned int> embedding_dofs ( (dim==2) ? 
					       (sub_degree+1) * (sub_degree+1) :
					       (sub_degree+1) * (sub_degree+1) * (sub_degree+1) );

    for (unsigned int i=0; i<( (dim==2) ? 
			       (sub_degree+1) * (sub_degree+1) :
			       (sub_degree+1) * (sub_degree+1) * (sub_degree+1) ); ++i)
    {
      // vertex
      if (i<GeometryInfo<dim>::vertices_per_cell)
        embedding_dofs[i] = i;
      // line
      else if (i<(GeometryInfo<dim>::vertices_per_cell + 
	    	  GeometryInfo<dim>::lines_per_cell * (sub_degree-1)))
      {
        const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell) %
	                       (sub_degree-1);
        const unsigned int line = (i - GeometryInfo<dim>::vertices_per_cell - j) / (sub_degree-1);

        embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell + 
	                    line * (degree-1) + j;
      }
      // quad
      else if (i<(GeometryInfo<dim>::vertices_per_cell + 
	  	  GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) + 
	          GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1))
      {
        const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell -
                                GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) % (sub_degree-1);
        const unsigned int k = ( (i - GeometryInfo<dim>::vertices_per_cell -
                                  GeometryInfo<dim>::lines_per_cell * (sub_degree-1) - j) / (sub_degree-1) ) % (sub_degree-1);
        const unsigned int face = (i - GeometryInfo<dim>::vertices_per_cell - 
                                   GeometryInfo<dim>::lines_per_cell * (sub_degree-1) - k * (sub_degree-1) - j) / ( (sub_degree-1) * (sub_degree-1) );

        embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell + 
                            GeometryInfo<dim>::lines_per_cell * (degree-1) +
	                    face * (degree-1) * (degree-1) +
                            k * (degree-1) + j;
      }
      // hex
      else if (i<(GeometryInfo<dim>::vertices_per_cell + 
	  	  GeometryInfo<dim>::lines_per_cell * (sub_degree-1)) + 
	          GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) +
	          GeometryInfo<dim>::hexes_per_cell * (sub_degree-1) * (sub_degree-1) * (sub_degree-1))
      {
	const unsigned int j = (i - GeometryInfo<dim>::vertices_per_cell -
                                GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
                                GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) ) % (sub_degree-1);
        const unsigned int k = ( (i - GeometryInfo<dim>::vertices_per_cell -
                                  GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
                                  GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) - j) / (sub_degree-1) ) % (sub_degree-1);
        const unsigned int l = (i - GeometryInfo<dim>::vertices_per_cell -
                                GeometryInfo<dim>::lines_per_cell * (sub_degree-1) -
	                        GeometryInfo<dim>::quads_per_cell * (sub_degree-1) * (sub_degree-1) - j - k * (sub_degree-1)) / ( (sub_degree-1) * (sub_degree-1) );
        
        embedding_dofs[i] = GeometryInfo<dim>::vertices_per_cell + 
                            GeometryInfo<dim>::lines_per_cell * (degree-1) +
                            GeometryInfo<dim>::quads_per_cell * (degree-1) * (degree-1) +
                            l * (degree-1) * (degree-1) + k * (degree-1) + j;
      }
    }

  return embedding_dofs;
  }
  else
  {
    Assert(false, ExcNotImplemented ());
    return std::vector<unsigned int> ();
  }
}



template <int dim>
unsigned int
FE_Q_Hierarchical<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_Q_Hierarchical<dim>::get_degree () const
{
  return degree;
};



template class FE_Q_Hierarchical<deal_II_dimension>;
