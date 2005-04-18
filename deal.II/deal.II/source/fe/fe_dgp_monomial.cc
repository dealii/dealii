//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------

#include <fe/fe_dgp_monomial.h>
#include <fe/fe_tools.h>

#ifdef HAVE_STD_STRINGSTREAM
#  include <sstream>
#else
#  include <strstream>
#endif


// namespace for some functions that are used in this file.
namespace 
{
				   // storage of hand-chosen support
				   // points
				   //
				   // For dim=2, dofs_per_cell of
				   // FE_DGPMonomial(k) is given by
				   // 0.5(k+1)(k+2), i.e.
				   // 
				   // k    0  1  2  3  4  5  6  7
				   // dofs 1  3  6 10 15 21 28 36
				   // 
				   // indirect access of unit points:
				   // the points for degree k are
				   // located at
				   //
				   // points[start_index[k]..start_index[k+1]-1]
  const unsigned int start_index2d[6]={0,1,4,10,20,35};
  const double points2d[35][2]=
  {{0,0},
   {0,0},{1,0},{0,1},
   {0,0},{1,0},{0,1},{1,1},{0.5,0},{0,0.5},
   {0,0},{1,0},{0,1},{1,1},{1./3.,0},{2./3.,0},{0,1./3.},{0,2./3.},{0.5,1},{1,0.5},
   {0,0},{1,0},{0,1},{1,1},{0.25,0},{0.5,0},{0.75,0},{0,0.25},{0,0.5},{0,0.75},{1./3.,1},{2./3.,1},{1,1./3.},{1,2./3.},{0.5,0.5}
  };

				   //
				   // For dim=3, dofs_per_cell of
				   // FE_DGPMonomial(k) is given by
				   // 1./6.(k+1)(k+2)(k+3), i.e.
				   // 
				   // k    0  1  2  3  4  5  6   7
				   // dofs 1  4 10 20 35 56 84 120
  const unsigned int start_index3d[6]={0,1,5,15/*,35*/};
  const double points3d[35][3]=
  {{0,0,0},
   {0,0,0},{1,0,0},{0,1,0},{0,0,1},
   {0,0,0},{1,0,0},{0,1,0},{0,0,1},{0.5,0,0},{0,0.5,0},{0,0,0.5},{1,1,0},{1,0,1},{0,1,1}
  };

  
  template<int dim>
  void generate_unit_points (const unsigned int,
			     std::vector<Point<dim> > &);

  template <>
  void generate_unit_points (const unsigned int k,
			     std::vector<Point<1> > &p)
  {
    Assert(p.size()==k+1, ExcDimensionMismatch(p.size(), k+1));
    const double h = 1./k;
    for (unsigned int i=0; i<p.size(); ++i)
      p[i](0)=i*h;
  }
  
  
  template <>
  void generate_unit_points (const unsigned int k,
			     std::vector<Point<2> > &p)
  {
    Assert(k<=4, ExcNotImplemented());
    Assert(p.size()==start_index2d[k+1]-start_index2d[k], ExcInternalError());
    for (unsigned int i=0; i<p.size(); ++i)
      {
	p[i](0)=points2d[start_index2d[k]+i][0];
	p[i](1)=points2d[start_index2d[k]+i][1];
      }
  }
  
  template <>
  void generate_unit_points (const unsigned int k,
			     std::vector<Point<3> > &p)
  {
    Assert(k<=2, ExcNotImplemented());
    Assert(p.size()==start_index3d[k+1]-start_index3d[k], ExcInternalError());
    for (unsigned int i=0; i<p.size(); ++i)
      {
	p[i](0)=points3d[start_index3d[k]+i][0];
	p[i](1)=points3d[start_index3d[k]+i][1];
	p[i](2)=points3d[start_index3d[k]+i][2];
      }
  }
}



template <int dim>
FE_DGPMonomial<dim>::FE_DGPMonomial (const unsigned int degree)
		:
		FE_Poly<PolynomialsP<dim>, dim> (
		  degree,
		  PolynomialsP<dim>(degree),
		  FiniteElementData<dim>(get_dpo_vector(degree), 1, degree),
		  std::vector<bool>(FiniteElementData<dim>(get_dpo_vector(degree), 1, degree).dofs_per_cell,true),
		  std::vector<std::vector<bool> >(FiniteElementData<dim>(
		    get_dpo_vector(degree), 1, degree).dofs_per_cell, std::vector<bool>(1,true)))
{
  Assert(this->poly_space.n()==this->dofs_per_cell, ExcInternalError());
  Assert(this->poly_space.degree()==this->degree, ExcInternalError());
  
				   // DG doesn't have constraints, so
				   // leave them empty

				   // initialize the interpolation
				   // matrices
  for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell; ++i)
    this->prolongation[i].reinit (this->dofs_per_cell,
				  this->dofs_per_cell);
  FETools::compute_embedding_matrices (*this, &this->prolongation[0]);
//  initialize_restriction ();

                                   // note, that these elements have
                                   // neither support nor face-support
                                   // points, so leave these fields
                                   // empty
}



template <int dim>
std::string
FE_DGPMonomial<dim>::get_name () const
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
  
  namebuf << "FE_DGPMonomial<" << dim << ">(" << this->degree << ")";

#ifndef HAVE_STD_STRINGSTREAM
  namebuf << std::ends;
#endif
  return namebuf.str();
}



template <int dim>
FiniteElement<dim> *
FE_DGPMonomial<dim>::clone() const
{
  return new FE_DGPMonomial<dim>(this->degree);
}



template <int dim>
void
FE_DGPMonomial<dim>::
get_interpolation_matrix (const FiniteElementBase<dim> &source_fe,
			  FullMatrix<double>           &interpolation_matrix) const
{  
  const FE_DGPMonomial<dim> *source_dgp_monomial
    = dynamic_cast<const FE_DGPMonomial<dim> *>(&source_fe);

  if (source_dgp_monomial)
    {
      				   // ok, source_fe is a DGP_Monomial
      				   // element. Then, the interpolation
      				   // matrix is simple
      const unsigned int m=interpolation_matrix.m();
      const unsigned int n=interpolation_matrix.n();
      Assert (m == this->dofs_per_cell, ExcDimensionMismatch (m, this->dofs_per_cell));
      Assert (n == source_dgp_monomial->dofs_per_cell,
	      ExcDimensionMismatch (n, source_dgp_monomial->dofs_per_cell));
      
      const unsigned int min_mn=
	interpolation_matrix.m()<interpolation_matrix.n() ?
	interpolation_matrix.m() : interpolation_matrix.n();
      
      for (unsigned int i=0; i<min_mn; ++i)
	interpolation_matrix(i,i)=1.;
    }
  else
    {
      std::vector<Point<dim> > unit_points(this->dofs_per_cell);
      generate_unit_points(this->degree, unit_points);

      FullMatrix<double> source_fe_matrix(unit_points.size(), source_fe.dofs_per_cell);
      for (unsigned int j=0; j<source_fe.dofs_per_cell; ++j)
	for (unsigned int k=0; k<unit_points.size(); ++k)
	  source_fe_matrix(k,j)=source_fe.shape_value(j, unit_points[k]);

      FullMatrix<double> this_matrix(this->dofs_per_cell, this->dofs_per_cell);
      for (unsigned int j=0; j<this->dofs_per_cell; ++j)
	for (unsigned int k=0; k<unit_points.size(); ++k)
	  this_matrix(k,j)=this->poly_space.compute_value (j, unit_points[k]);

      this_matrix.gauss_jordan();
      
      this_matrix.mmult(interpolation_matrix, source_fe_matrix);
    }
}



template <int dim>
void
FE_DGPMonomial<dim>::initialize_restriction ()
{
  Assert(false, ExcNotImplemented());
}


//---------------------------------------------------------------------------
// Auxiliary functions
//---------------------------------------------------------------------------


template <int dim>
std::vector<unsigned int>
FE_DGPMonomial<dim>::get_dpo_vector(unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, 0U);
  dpo[dim] = ++deg;
  for (unsigned int i=1;i<dim;++i)
    {
      dpo[dim] *= deg+i;
      dpo[dim] /= i+1;
    }
  return dpo;
}


#if deal_II_dimension == 1

template <>
bool
FE_DGPMonomial<1>::has_support_on_face (const unsigned int,
					const unsigned int face_index) const
{
  return face_index==1 || (face_index==0 && this->degree==0);
}

#endif

#if deal_II_dimension == 2

template <>
bool
FE_DGPMonomial<2>::has_support_on_face (const unsigned int shape_index,
					const unsigned int face_index) const
{
  bool support_on_face=false;
  if (face_index==1 || face_index==2)
    support_on_face=true;
  else
    {
      unsigned int degrees[2];
      this->poly_space.directional_degrees(shape_index, degrees);
      if ((face_index==0 && degrees[1]==0) ||
	  (face_index==3 && degrees[0]==0))
	support_on_face=true;
    }
  return support_on_face;
}

#endif

#if deal_II_dimension == 3

template <>
bool
FE_DGPMonomial<3>::has_support_on_face (const unsigned int shape_index,
					const unsigned int face_index) const
{
  bool support_on_face=false;
  if (face_index==1 || face_index==3 || face_index==4)
    support_on_face=true;
  else
    {
      unsigned int degrees[3];
      this->poly_space.directional_degrees(shape_index, degrees);
      if ((face_index==0 && degrees[1]==0) ||
	  (face_index==2 && degrees[2]==0) ||
	  (face_index==5 && degrees[0]==0))
	support_on_face=true;
    }
  return support_on_face;
}

#endif


template <int dim>
unsigned int
FE_DGPMonomial<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template class FE_DGPMonomial<deal_II_dimension>;
