//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <fe/fe_q_hierarchical.h>

#include <cmath>
#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


namespace 
{
  inline
  std::vector<unsigned int>
  invert_numbering (const std::vector<unsigned int> &in)
  {
    std::vector<unsigned int> out (in.size());
    for (unsigned int i=0; i<in.size(); ++i)
      out[in[i]]=i;
    return out;
  }
}



template <int dim>
FE_Q_Hierarchical<dim>::FE_Q_Hierarchical (const unsigned int degree)
		:
		FE_Poly<TensorProductPolynomials<dim>, dim> (
		  Polynomials::Hierarchical::generate_complete_basis(degree),
		  FiniteElementData<dim>(get_dpo_vector(degree),1, degree, FiniteElementData<dim>::H1),
		  std::vector<bool> (FiniteElementData<dim>(
		    get_dpo_vector(degree),1, degree).dofs_per_cell, false),
		  std::vector<std::vector<bool> >(FiniteElementData<dim>(
		    get_dpo_vector(degree),1, degree).dofs_per_cell, std::vector<bool>(1,true))),
						    face_renumber(face_fe_q_hierarchical_to_hierarchic_numbering (degree))
{
  this->poly_space.set_numbering(
    hierarchic_to_fe_q_hierarchical_numbering(*this));
  
				   // The matrix @p{dofs_cell} contains the 
				   // values of the linear functionals of 
				   // the master 1d cell applied to the 
				   // shape functions of the two 1d subcells.
				   // The matrix @p{dofs_subcell} constains
				   // the values of the linear functionals 
				   // on each 1d subcell applied to the 
				   // shape functions on the master 1d 
				   // subcell. 
				   // We use @p{dofs_cell} and 
				   // @p{dofs_subcell} to compute the 
				   // @p{prolongation}, @p{restriction} and 
				   // @p{interface_constraints} matrices 
				   // for all dimensions.
  std::vector<FullMatrix<double> >
    dofs_cell (GeometryInfo<1>::children_per_cell,
	       FullMatrix<double> (2*this->dofs_per_vertex + this->dofs_per_line,
				   2*this->dofs_per_vertex + this->dofs_per_line));
  std::vector<FullMatrix<double> >
    dofs_subcell (GeometryInfo<1>::children_per_cell,
		  FullMatrix<double> (2*this->dofs_per_vertex + this->dofs_per_line,
				      2*this->dofs_per_vertex + this->dofs_per_line));
				   // build these fields, as they are
				   // needed as auxiliary fields later
				   // on
  build_dofs_cell (dofs_cell, dofs_subcell);

				   // then use them to initialize
				   // other fields
  initialize_constraints (dofs_subcell);
  initialize_embedding_and_restriction (dofs_cell, dofs_subcell);

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
}



template <int dim>
std::string
FE_Q_Hierarchical<dim>::get_name () const
{
				   // note that the
				   // FETools::get_fe_from_name
				   // function depends on the
				   // particular format of the string
				   // this function returns, so they
				   // have to be kept in synch

#ifdef HAVE_STD_STRINGSTREAM
  std::ostringstream namebuf;
#else
  std::ostrstream namebuf;
#endif
  
  namebuf << "FE_Q_Hierarchical<" << dim << ">(" << this->degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_Q_Hierarchical<dim>::clone() const
{
  return new FE_Q_Hierarchical<dim>(this->degree);
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
void
FE_Q_Hierarchical<dim>::build_dofs_cell (std::vector<FullMatrix<double> > &dofs_cell,
					 std::vector<FullMatrix<double> > &dofs_subcell) const
{
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;

  for (unsigned int c=0; c<GeometryInfo<1>::children_per_cell; ++c)
    for (unsigned int j=0; j<dofs_1d; ++j)
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



template <int dim>
void
FE_Q_Hierarchical<dim>::
initialize_constraints (const std::vector<FullMatrix<double> > &dofs_subcell)
{
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;
  const unsigned int degree=this->degree;

  this->interface_constraints
    .TableBase<2,double>::reinit (this->interface_constraints_size());

  switch (dim)
    {
      case 1:
      {
					 // no constraints in 1d
	break;
      }
       
      case 2:
      {
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
      }
       
      case 3:
      {
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
	    for (unsigned int j=0; j<(degree-1); j++)
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
	    for (unsigned int j=0; j<(degree-1); j++)
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
	    for (unsigned int j=0; j<(degree-1); j++)
	      for (unsigned int k=0; k<(degree-1); k++)
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
       
      default:
	    Assert (false, ExcNotImplemented());
    }
}



template <int dim>
void
FE_Q_Hierarchical<dim>::
initialize_embedding_and_restriction (const std::vector<FullMatrix<double> > &dofs_cell,
				      const std::vector<FullMatrix<double> > &dofs_subcell)
{
  const unsigned int dofs_1d = 2*this->dofs_per_vertex + this->dofs_per_line;
  const std::vector<unsigned int> &renumber=
    this->poly_space.get_numbering();
  
  for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
    {
      this->prolongation[c].reinit (this->dofs_per_cell, this->dofs_per_cell);
      this->restriction[c].reinit (this->dofs_per_cell, this->dofs_per_cell);
    }

				   // the 1d case is particularly
				   // simple, so special case it:
  if (dim==1)
    {
      for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
	{
	  this->prolongation[c].fill (dofs_subcell[c]);
	  this->restriction[c].fill (dofs_cell[c]);
	}
      return;
    }

				   // for higher dimensions, things
				   // are a little more tricky:
  
				   // j loops over dofs in the
				   // subcell.  These are the rows in
				   // the embedding matrix.
				   //
				   // i loops over the dofs in the
				   // master cell. These are the
				   // columns in the embedding matrix.
  for (unsigned int j=0; j<this->dofs_per_cell; ++j)
    for (unsigned int i=0; i<this->dofs_per_cell; ++i)
      switch (dim)
	{
	  case 2:
	  {
	    for (unsigned int c=0; c<GeometryInfo<2>::children_per_cell; ++c)
	      {
		unsigned int c0 = ((c==1) || (c==2)) ? 1 : 0;
		unsigned int c1 = ((c==2) || (c==3)) ? 1 : 0;

		this->prolongation[c](j,i) = 
		  dofs_subcell[c0](renumber[j] % dofs_1d,
				   renumber[i] % dofs_1d) *
		  dofs_subcell[c1]((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d,
				   (renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d);

		this->restriction[c](j,i) = 
		  dofs_cell[c0](renumber[j] % dofs_1d,
				renumber[i] % dofs_1d) *
		  dofs_cell[c1]((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d,
				(renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d);
	      }
	    break;
	  }
	     
	  case 3:
	  {
	    for (unsigned int c=0; c<GeometryInfo<3>::children_per_cell; ++c)
	      {
		unsigned int c0 = ((c==1) || (c==2) || (c==5) || (c==6)) ? 1 : 0;
		unsigned int c1 = ((c==4) || (c==5) || (c==6) || (c==7)) ? 1 : 0;
		unsigned int c2 = ((c==2) || (c==3) || (c==6) || (c==7)) ? 1 : 0;

		this->prolongation[c](j,i) = 
		  dofs_subcell[c0](renumber[j] % dofs_1d,
				   renumber[i] % dofs_1d) *
		  dofs_subcell[c1](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) % dofs_1d,
				   ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
		  dofs_subcell[c2](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d - (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
				   ((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d - (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);

		this->restriction[c](j,i) = 
		  dofs_cell[c0](renumber[j] % dofs_1d,
				renumber[i] % dofs_1d) *
		  dofs_cell[c1](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d) % dofs_1d,
				((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d) % dofs_1d) *
		  dofs_cell[c2](((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d - (((renumber[j] - (renumber[j] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d,
				((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d - (((renumber[i] - (renumber[i] % dofs_1d)) / dofs_1d ) % dofs_1d)) / dofs_1d);
	      }
	    break;
	  }

	  default:
		Assert (false, ExcNotImplemented());
	}
}



template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_support_points ()
{
				   // number of points: (degree+1)^dim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= this->degree+1;
  
  this->unit_support_points.resize(n);

  const std::vector<unsigned int> &index_map_inverse=
    this->poly_space.get_numbering_inverse();
  
  Point<dim> p;
                                   // the method of numbering allows
                                   // each dof to be associated with a
                                   // support point. There is
                                   // only one support point per
                                   // vertex, line, quad, hex, etc.
                                   //
                                   // note, however, that the support
                                   // points thus associated with
                                   // shape functions are not unique:
                                   // the linear shape functions are
                                   // associated with the vertices,
                                   // but all others are associated
                                   // with either line, quad, or hex
                                   // midpoints, and there may be
                                   // multiple shape functions
                                   // associated with them. there
                                   // really is no other useful
                                   // numbering, since the
                                   // hierarchical shape functions do
                                   // not vanish at all-but-one
                                   // interpolation points (like the
                                   // Lagrange functions used in
                                   // FE_Q), so there's not much we
                                   // can do here.
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((dim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
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
	  this->unit_support_points[index_map_inverse[k++]] = p;
	};
}


#if deal_II_dimension == 1

template <>
void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
}

#endif


template <int dim>
void FE_Q_Hierarchical<dim>::initialize_unit_face_support_points ()
{
  const unsigned int codim = dim-1;
  
				   // number of points: (degree+1)^codim
  unsigned int n = this->degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= this->degree+1;
  
  this->unit_face_support_points.resize(n);
  
  Point<codim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? this->degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? this->degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=this->degree; ++ix)
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
}


				 // we use same dpo_vector as FE_Q
template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 1U);
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::
hierarchic_to_fe_q_hierarchical_numbering (const FiniteElementData<dim> &fe)
{
  Assert (fe.n_components() == 1, ExcInternalError());
  std::vector<unsigned int> h2l(fe.dofs_per_cell);

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
	for (unsigned int i=0; i<fe.dofs_per_cell; ++i)
	  h2l[i] = i;

	break;
      }

      case 2:
      {
					 // Example: degree=3
					 //
					 // hierarchical numbering:
					 //  3  8  9  2
					 // 11 14 15  7
					 // 10 12 13  6
					 //  0  4  5  1
					 //
					 // fe_q_hierarchical numbering:
					 //  4  6  7  5
					 // 12 14 15 13
					 //  8 10 11  9
					 //  0  2  3  1
	unsigned int next_index = 0;
					 // first the four vertices
	h2l[next_index++] = 0;
	h2l[next_index++] = 1;
	h2l[next_index++] = n+1;
	h2l[next_index++] = n;
					 // bottom line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = 2+i;
					 // right line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n+1;
					 // top line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n+2+i;
					 // left line
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n;
					 // inside quad
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n+2+j;

	Assert (next_index == fe.dofs_per_cell, ExcInternalError());

	break;
      }

      case 3:
      {
	unsigned int next_index = 0;
	const unsigned int n2=n*n;
					 // first the eight vertices
					 // front face, counterclock wise
	h2l[next_index++] = 0;
	h2l[next_index++] = 1;
	h2l[next_index++] = n2+1;
	h2l[next_index++] = n2;
					 // back face, counterclock wise
	h2l[next_index++] = n;
	h2l[next_index++] = n+1;
	h2l[next_index++] = n2+n+1;
	h2l[next_index++] = n2+n;

					 // now the lines
					 // front face, counterclock wise
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = 2+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n2+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n2+2+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n2;
					 // back face, counterclock wise
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n+2+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n2+n+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n2+n+2+i;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n2+n;
					 // lines in y-direction,
					 // counterclock wise
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = (2+i)*n+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n2+(2+i)*n+1;
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  h2l[next_index++] = n2+(2+i)*n;

					 // inside quads
	Assert (fe.dofs_per_quad == fe.dofs_per_line*fe.dofs_per_line,
		ExcInternalError());
					 // front face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n2+2+j;
					 // back face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n2+n+2+j;
					 // bottom face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n+2+j;
					 // right face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n2+(2+j)*n+1;
					 // top face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = n2+(2+i)*n+2+j;
					 // left face
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    h2l[next_index++] = (2+i)*n2+(2+j)*n;

					 // inside hex
	Assert (fe.dofs_per_hex == fe.dofs_per_quad*fe.dofs_per_line,
		ExcInternalError());
	for (unsigned int i=0; i<fe.dofs_per_line; ++i)
	  for (unsigned int j=0; j<fe.dofs_per_line; ++j)
	    for (unsigned int k=0; k<fe.dofs_per_line; ++k)
	      h2l[next_index++]	= (2+i)*n2+(2+j)*n+2+k;

	Assert (next_index == fe.dofs_per_cell, ExcInternalError());
	
	break;
      }       

      default:
	    Assert (false, ExcNotImplemented());
    }
  return h2l;
}


template <int dim>
std::vector<unsigned int>
FE_Q_Hierarchical<dim>::
face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int degree)
{
  FiniteElementData<dim-1> fe_data(FE_Q_Hierarchical<dim-1>::get_dpo_vector(degree),1);
  return invert_numbering(FE_Q_Hierarchical<dim-1>::
			  hierarchic_to_fe_q_hierarchical_numbering (fe_data));
}


#if (deal_II_dimension == 1)

template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_fe_q_hierarchical_to_hierarchic_numbering (const unsigned int)
{
  return std::vector<unsigned int> ();
}

#endif



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
  Assert ((sub_degree>0) && (sub_degree<=this->degree),
	  ExcIndexRange(sub_degree, 1, this->degree));

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
  else if (sub_degree==this->degree)
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
				  line * (this->degree-1) + j;
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
				  GeometryInfo<dim>::lines_per_cell * (this->degree-1) +
				  face * (this->degree-1) * (this->degree-1) +
				  k * (this->degree-1) + j;
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
				  GeometryInfo<dim>::lines_per_cell * (this->degree-1) +
				  GeometryInfo<dim>::quads_per_cell * (this->degree-1) * (this->degree-1) +
				  l * (this->degree-1) * (this->degree-1) + k * (this->degree-1) + j;
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



template class FE_Q_Hierarchical<deal_II_dimension>;
