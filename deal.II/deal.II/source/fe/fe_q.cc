//----------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001 by the deal.II authors
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
#include <fe/mapping_q1.h>
#include <fe/fe_q.h>
#include <fe/fe_values.h>



//TODO:[RH] move build_renumbering to mapping class

template <int dim>
FE_Q<dim>::FE_Q (const unsigned int degree)
		:
		FiniteElement<dim> (FiniteElementData<dim>(get_dpo_vector(degree),1),
				    std::vector<bool> (1,false)),
		degree(degree),
		renumber(dofs_per_cell, 0),
		renumber_inverse(dofs_per_cell, 0),
		face_renumber(dofs_per_face, 0),
		poly(0)
{
				   // Q0 elements cannot be
				   // continuous, use FE_DGQ<dim>(0)
				   // instead
  Assert (degree>0, ExcNotImplemented());
  std::vector<LagrangeEquidistant> v;
  for (unsigned int i=0;i<=degree;++i)
    v.push_back(LagrangeEquidistant(degree,i));
  
  poly = new TensorProductPolynomials<dim> (v);

				   // do some internal book-keeping on
				   // cells and faces. if in 1d, the
				   // face function is empty
  build_renumbering (*this, degree, renumber);
  build_face_renumbering (degree, face_renumber);
  
				   // build inverse of renumbering
				   // vector
  for (unsigned int i=0; i<dofs_per_cell; ++i)
    renumber_inverse[renumber[i]]=i;

				   // copy constraint matrices if they
				   // are defined. otherwise set them
				   // to invalid size
  if ((dim>1) && (degree<Matrices::n_constraint_matrices+1))
    interface_constraints.fill (Matrices::constraint_matrices[degree-1]);
  else
    interface_constraints.reinit(0,0);

				   // next copy over embedding
				   // matrices if they are defined
  if ((degree < Matrices::n_embedding_matrices+1) &&
      (Matrices::embedding[degree-1][0] != 0))
    for (unsigned int c=0; c<GeometryInfo<dim>::children_per_cell; ++c)
      prolongation[c].fill (Matrices::embedding[degree-1][c]);
  else
    for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell;++i)
      prolongation[i].reinit(0,0);

				   // then fill restriction
				   // matrices. they are hardcoded for
				   // the first few elements
  switch (dim)
    {
      case 1:
	    switch (degree)
	      {
		case 1:
		      restriction[0](0,0) = 1;
		      restriction[1](1,1) = 1;
		      break;
		case 2:
		      restriction[0](0,0) = 1;
		      restriction[0](2,1) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](1,1) = 1;
		      break;
		case 3:
		      restriction[0](0,0) = 1;
		      restriction[0](2,3) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](3,2) = 1;
		      break;
		case 4:
		      restriction[0](0,0) = 1;
		      restriction[0](2,3) = 1;
		      restriction[0](3,1) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](3,0) = 1;
		      restriction[1](4,3) = 1;
		      break;
		default:
		      for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell;++i)
			restriction[i].reinit(0,0);
	      }
	    break;
      case 2:
	    switch (degree)
	      {
		case 1:
		      restriction[0](0,0) = 1;
		      restriction[1](1,1) = 1;
		      restriction[2](2,2) = 1;
		      restriction[3](3,3) = 1;
		      break;
		case 2:
		      restriction[0](0,0) = 1;
		      restriction[0](4,1) = 1;
		      restriction[0](7,3) = 1;
		      restriction[0](8,2) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](4,0) = 1;
		      restriction[1](5,2) = 1;
		      restriction[1](8,3) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](5,1) = 1;
		      restriction[2](6,3) = 1;
		      restriction[2](8,0) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](6,2) = 1;
		      restriction[3](7,0) = 1;
		      restriction[3](8,1) = 1;
		      break;
		case 3:
		      restriction[0](0,0) = 1;
		      restriction[0](4,5) = 1;
		      restriction[0](10,11) = 1;
		      restriction[0](12,15) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](5,4) = 1;
		      restriction[1](6,7) = 1;
		      restriction[1](13,14) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](7,6) = 1;
		      restriction[2](9,8) = 1;
		      restriction[2](15,12) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](8,9) = 1;
		      restriction[3](11,10) = 1;
		      restriction[3](14,13) = 1;
		      break;
		case 4:
		      restriction[0](0,0) = 1;
		      restriction[0](4,5) = 1;
		      restriction[0](5,1) = 1;
		      restriction[0](13,14) = 1;
		      restriction[0](14,3) = 1;
		      restriction[0](16,20) = 1;
		      restriction[0](17,8) = 1;
		      restriction[0](19,11) = 1;
		      restriction[0](20,2) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](5,0) = 1;
		      restriction[1](6,5) = 1;
		      restriction[1](7,8) = 1;
		      restriction[1](8,2) = 1;
		      restriction[1](17,14) = 1;
		      restriction[1](18,20) = 1;
		      restriction[1](20,3) = 1;
		      restriction[1](21,11) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](8,1) = 1;
		      restriction[2](9,8) = 1;
		      restriction[2](11,3) = 1;
		      restriction[2](12,11) = 1;
		      restriction[2](20,0) = 1;
		      restriction[2](21,5) = 1;
		      restriction[2](23,14) = 1;
		      restriction[2](24,20) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](10,11) = 1;
		      restriction[3](11,2) = 1;
		      restriction[3](14,0) = 1;
		      restriction[3](15,14) = 1;
		      restriction[3](19,5) = 1;
		      restriction[3](20,1) = 1;
		      restriction[3](22,20) = 1;
		      restriction[3](23,8) = 1;
		      break;
		default:
		      for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell;++i)
			restriction[i].reinit(0,0);
	      }
	    break;
      case 3:
	    switch (degree)
	      {
		case 1:
		      restriction[0](0,0) = 1;
		      restriction[1](1,1) = 1;
		      restriction[2](2,2) = 1;
		      restriction[3](3,3) = 1;
		      restriction[4](4,4) = 1;
		      restriction[5](5,5) = 1;
		      restriction[6](6,6) = 1;
		      restriction[7](7,7) = 1;
		      break;
		case 2:
		      restriction[0](0,0) = 1;
		      restriction[0](8,1) = 1;
		      restriction[0](11,3) = 1;
		      restriction[0](16,4) = 1;
		      restriction[0](20,2) = 1;
		      restriction[0](22,5) = 1;
		      restriction[0](25,7) = 1;
		      restriction[0](26,6) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](8,0) = 1;
		      restriction[1](9,2) = 1;
		      restriction[1](17,5) = 1;
		      restriction[1](20,3) = 1;
		      restriction[1](22,4) = 1;
		      restriction[1](23,6) = 1;
		      restriction[1](26,7) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](9,1) = 1;
		      restriction[2](10,3) = 1;
		      restriction[2](18,6) = 1;
		      restriction[2](20,0) = 1;
		      restriction[2](23,5) = 1;
		      restriction[2](24,7) = 1;
		      restriction[2](26,4) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](10,2) = 1;
		      restriction[3](11,0) = 1;
		      restriction[3](19,7) = 1;
		      restriction[3](20,1) = 1;
		      restriction[3](24,6) = 1;
		      restriction[3](25,4) = 1;
		      restriction[3](26,5) = 1;
		      restriction[4](4,4) = 1;
		      restriction[4](12,5) = 1;
		      restriction[4](15,7) = 1;
		      restriction[4](16,0) = 1;
		      restriction[4](21,6) = 1;
		      restriction[4](22,1) = 1;
		      restriction[4](25,3) = 1;
		      restriction[4](26,2) = 1;
		      restriction[5](5,5) = 1;
		      restriction[5](12,4) = 1;
		      restriction[5](13,6) = 1;
		      restriction[5](17,1) = 1;
		      restriction[5](21,7) = 1;
		      restriction[5](22,0) = 1;
		      restriction[5](23,2) = 1;
		      restriction[5](26,3) = 1;
		      restriction[6](6,6) = 1;
		      restriction[6](13,5) = 1;
		      restriction[6](14,7) = 1;
		      restriction[6](18,2) = 1;
		      restriction[6](21,4) = 1;
		      restriction[6](23,1) = 1;
		      restriction[6](24,3) = 1;
		      restriction[6](26,0) = 1;
		      restriction[7](7,7) = 1;
		      restriction[7](14,6) = 1;
		      restriction[7](15,4) = 1;
		      restriction[7](19,3) = 1;
		      restriction[7](21,5) = 1;
		      restriction[7](24,2) = 1;
		      restriction[7](25,0) = 1;
		      restriction[7](26,1) = 1;
		      break;
		case 3:
		      restriction[0](0,0) = 1;
		      restriction[0](8,9) = 1;
		      restriction[0](14,15) = 1;
		      restriction[0](24,25) = 1;
		      restriction[0](32,35) = 1;
		      restriction[0](40,43) = 1;
		      restriction[0](52,55) = 1;
		      restriction[0](56,63) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](9,8) = 1;
		      restriction[1](10,11) = 1;
		      restriction[1](26,27) = 1;
		      restriction[1](33,34) = 1;
		      restriction[1](41,42) = 1;
		      restriction[1](44,47) = 1;
		      restriction[1](57,62) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](11,10) = 1;
		      restriction[2](13,12) = 1;
		      restriction[2](28,29) = 1;
		      restriction[2](35,32) = 1;
		      restriction[2](46,45) = 1;
		      restriction[2](49,50) = 1;
		      restriction[2](61,58) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](12,13) = 1;
		      restriction[3](15,14) = 1;
		      restriction[3](30,31) = 1;
		      restriction[3](34,33) = 1;
		      restriction[3](48,51) = 1;
		      restriction[3](54,53) = 1;
		      restriction[3](60,59) = 1;
		      restriction[4](4,4) = 1;
		      restriction[4](16,17) = 1;
		      restriction[4](22,23) = 1;
		      restriction[4](25,24) = 1;
		      restriction[4](36,39) = 1;
		      restriction[4](42,41) = 1;
		      restriction[4](53,54) = 1;
		      restriction[4](58,61) = 1;
		      restriction[5](5,5) = 1;
		      restriction[5](17,16) = 1;
		      restriction[5](18,19) = 1;
		      restriction[5](27,26) = 1;
		      restriction[5](37,38) = 1;
		      restriction[5](43,40) = 1;
		      restriction[5](45,46) = 1;
		      restriction[5](59,60) = 1;
		      restriction[6](6,6) = 1;
		      restriction[6](19,18) = 1;
		      restriction[6](21,20) = 1;
		      restriction[6](29,28) = 1;
		      restriction[6](39,36) = 1;
		      restriction[6](47,44) = 1;
		      restriction[6](51,48) = 1;
		      restriction[6](63,56) = 1;
		      restriction[7](7,7) = 1;
		      restriction[7](20,21) = 1;
		      restriction[7](23,22) = 1;
		      restriction[7](31,30) = 1;
		      restriction[7](38,37) = 1;
		      restriction[7](50,49) = 1;
		      restriction[7](55,52) = 1;
		      restriction[7](62,57) = 1;
		      break;
		case 4:
		      restriction[0](0,0) = 1;
		      restriction[0](8,9) = 1;
		      restriction[0](9,1) = 1;
		      restriction[0](17,18) = 1;
		      restriction[0](18,3) = 1;
		      restriction[0](32,33) = 1;
		      restriction[0](33,4) = 1;
		      restriction[0](44,48) = 1;
		      restriction[0](45,12) = 1;
		      restriction[0](47,15) = 1;
		      restriction[0](48,2) = 1;
		      restriction[0](62,66) = 1;
		      restriction[0](63,36) = 1;
		      restriction[0](65,21) = 1;
		      restriction[0](66,5) = 1;
		      restriction[0](89,93) = 1;
		      restriction[0](90,30) = 1;
		      restriction[0](92,42) = 1;
		      restriction[0](93,7) = 1;
		      restriction[0](98,111) = 1;
		      restriction[0](99,75) = 1;
		      restriction[0](101,57) = 1;
		      restriction[0](102,24) = 1;
		      restriction[0](107,84) = 1;
		      restriction[0](108,39) = 1;
		      restriction[0](110,27) = 1;
		      restriction[0](111,6) = 1;
		      restriction[1](1,1) = 1;
		      restriction[1](9,0) = 1;
		      restriction[1](10,9) = 1;
		      restriction[1](11,12) = 1;
		      restriction[1](12,2) = 1;
		      restriction[1](35,36) = 1;
		      restriction[1](36,5) = 1;
		      restriction[1](45,18) = 1;
		      restriction[1](46,48) = 1;
		      restriction[1](48,3) = 1;
		      restriction[1](49,15) = 1;
		      restriction[1](63,33) = 1;
		      restriction[1](64,66) = 1;
		      restriction[1](66,4) = 1;
		      restriction[1](67,21) = 1;
		      restriction[1](71,75) = 1;
		      restriction[1](72,24) = 1;
		      restriction[1](74,39) = 1;
		      restriction[1](75,6) = 1;
		      restriction[1](99,93) = 1;
		      restriction[1](100,111) = 1;
		      restriction[1](102,30) = 1;
		      restriction[1](103,57) = 1;
		      restriction[1](108,42) = 1;
		      restriction[1](109,84) = 1;
		      restriction[1](111,7) = 1;
		      restriction[1](112,27) = 1;
		      restriction[2](2,2) = 1;
		      restriction[2](12,1) = 1;
		      restriction[2](13,12) = 1;
		      restriction[2](15,3) = 1;
		      restriction[2](16,15) = 1;
		      restriction[2](38,39) = 1;
		      restriction[2](39,6) = 1;
		      restriction[2](48,0) = 1;
		      restriction[2](49,9) = 1;
		      restriction[2](51,18) = 1;
		      restriction[2](52,48) = 1;
		      restriction[2](74,36) = 1;
		      restriction[2](75,5) = 1;
		      restriction[2](77,75) = 1;
		      restriction[2](78,24) = 1;
		      restriction[2](81,42) = 1;
		      restriction[2](82,84) = 1;
		      restriction[2](84,7) = 1;
		      restriction[2](85,27) = 1;
		      restriction[2](108,33) = 1;
		      restriction[2](109,66) = 1;
		      restriction[2](111,4) = 1;
		      restriction[2](112,21) = 1;
		      restriction[2](117,93) = 1;
		      restriction[2](118,111) = 1;
		      restriction[2](120,30) = 1;
		      restriction[2](121,57) = 1;
		      restriction[3](3,3) = 1;
		      restriction[3](14,15) = 1;
		      restriction[3](15,2) = 1;
		      restriction[3](18,0) = 1;
		      restriction[3](19,18) = 1;
		      restriction[3](41,42) = 1;
		      restriction[3](42,7) = 1;
		      restriction[3](47,9) = 1;
		      restriction[3](48,1) = 1;
		      restriction[3](50,48) = 1;
		      restriction[3](51,12) = 1;
		      restriction[3](80,84) = 1;
		      restriction[3](81,39) = 1;
		      restriction[3](83,27) = 1;
		      restriction[3](84,6) = 1;
		      restriction[3](92,33) = 1;
		      restriction[3](93,4) = 1;
		      restriction[3](95,93) = 1;
		      restriction[3](96,30) = 1;
		      restriction[3](107,66) = 1;
		      restriction[3](108,36) = 1;
		      restriction[3](110,21) = 1;
		      restriction[3](111,5) = 1;
		      restriction[3](116,111) = 1;
		      restriction[3](117,75) = 1;
		      restriction[3](119,57) = 1;
		      restriction[3](120,24) = 1;
		      restriction[4](4,4) = 1;
		      restriction[4](20,21) = 1;
		      restriction[4](21,5) = 1;
		      restriction[4](29,30) = 1;
		      restriction[4](30,7) = 1;
		      restriction[4](33,0) = 1;
		      restriction[4](34,33) = 1;
		      restriction[4](53,57) = 1;
		      restriction[4](54,24) = 1;
		      restriction[4](56,27) = 1;
		      restriction[4](57,6) = 1;
		      restriction[4](65,9) = 1;
		      restriction[4](66,1) = 1;
		      restriction[4](68,66) = 1;
		      restriction[4](69,36) = 1;
		      restriction[4](90,18) = 1;
		      restriction[4](91,93) = 1;
		      restriction[4](93,3) = 1;
		      restriction[4](94,42) = 1;
		      restriction[4](101,48) = 1;
		      restriction[4](102,12) = 1;
		      restriction[4](104,111) = 1;
		      restriction[4](105,75) = 1;
		      restriction[4](110,15) = 1;
		      restriction[4](111,2) = 1;
		      restriction[4](113,84) = 1;
		      restriction[4](114,39) = 1;
		      restriction[5](5,5) = 1;
		      restriction[5](21,4) = 1;
		      restriction[5](22,21) = 1;
		      restriction[5](23,24) = 1;
		      restriction[5](24,6) = 1;
		      restriction[5](36,1) = 1;
		      restriction[5](37,36) = 1;
		      restriction[5](54,30) = 1;
		      restriction[5](55,57) = 1;
		      restriction[5](57,7) = 1;
		      restriction[5](58,27) = 1;
		      restriction[5](66,0) = 1;
		      restriction[5](67,9) = 1;
		      restriction[5](69,33) = 1;
		      restriction[5](70,66) = 1;
		      restriction[5](72,12) = 1;
		      restriction[5](73,75) = 1;
		      restriction[5](75,2) = 1;
		      restriction[5](76,39) = 1;
		      restriction[5](102,18) = 1;
		      restriction[5](103,48) = 1;
		      restriction[5](105,93) = 1;
		      restriction[5](106,111) = 1;
		      restriction[5](111,3) = 1;
		      restriction[5](112,15) = 1;
		      restriction[5](114,42) = 1;
		      restriction[5](115,84) = 1;
		      restriction[6](6,6) = 1;
		      restriction[6](24,5) = 1;
		      restriction[6](25,24) = 1;
		      restriction[6](27,7) = 1;
		      restriction[6](28,27) = 1;
		      restriction[6](39,2) = 1;
		      restriction[6](40,39) = 1;
		      restriction[6](57,4) = 1;
		      restriction[6](58,21) = 1;
		      restriction[6](60,30) = 1;
		      restriction[6](61,57) = 1;
		      restriction[6](75,1) = 1;
		      restriction[6](76,36) = 1;
		      restriction[6](78,12) = 1;
		      restriction[6](79,75) = 1;
		      restriction[6](84,3) = 1;
		      restriction[6](85,15) = 1;
		      restriction[6](87,42) = 1;
		      restriction[6](88,84) = 1;
		      restriction[6](111,0) = 1;
		      restriction[6](112,9) = 1;
		      restriction[6](114,33) = 1;
		      restriction[6](115,66) = 1;
		      restriction[6](120,18) = 1;
		      restriction[6](121,48) = 1;
		      restriction[6](123,93) = 1;
		      restriction[6](124,111) = 1;
		      restriction[7](7,7) = 1;
		      restriction[7](26,27) = 1;
		      restriction[7](27,6) = 1;
		      restriction[7](30,4) = 1;
		      restriction[7](31,30) = 1;
		      restriction[7](42,3) = 1;
		      restriction[7](43,42) = 1;
		      restriction[7](56,21) = 1;
		      restriction[7](57,5) = 1;
		      restriction[7](59,57) = 1;
		      restriction[7](60,24) = 1;
		      restriction[7](83,15) = 1;
		      restriction[7](84,2) = 1;
		      restriction[7](86,84) = 1;
		      restriction[7](87,39) = 1;
		      restriction[7](93,0) = 1;
		      restriction[7](94,33) = 1;
		      restriction[7](96,18) = 1;
		      restriction[7](97,93) = 1;
		      restriction[7](110,9) = 1;
		      restriction[7](111,1) = 1;
		      restriction[7](113,66) = 1;
		      restriction[7](114,36) = 1;
		      restriction[7](119,48) = 1;
		      restriction[7](120,12) = 1;
		      restriction[7](122,111) = 1;
		      restriction[7](123,75) = 1;
		      break;
		default:
		      for (unsigned int i=0; i<GeometryInfo<dim>::children_per_cell;++i)
			restriction[i].reinit(0,0);
	      }
	    break;
      default:
	    Assert (false,ExcNotImplemented());
    }

				   // finally fill in support points
				   // on cell and face
  initialize_unit_support_points ();
  initialize_unit_face_support_points ();
};



template <int dim>
FE_Q<dim>::~FE_Q ()
{
				   // delete poly member and set it to
				   // zero to prevent accidental use
  delete poly;
  poly = 0;
}



template <int dim>
FiniteElement<dim> *
FE_Q<dim>::clone() const
{
  return new FE_Q<dim>(degree);
}



template <int dim>
double
FE_Q<dim>::shape_value (const unsigned int i,
			const Point<dim> &p) const
{
  return poly->compute_value(renumber_inverse[i], p);
}



template <int dim>
Tensor<1,dim>
FE_Q<dim>::shape_grad (const unsigned int i,
		       const Point<dim> &p) const
{
  return poly->compute_grad(renumber_inverse[i], p);
}



template <int dim>
Tensor<2,dim>
FE_Q<dim>::shape_grad_grad (const unsigned int i,
			    const Point<dim> &p) const
{
  return poly->compute_grad_grad(renumber_inverse[i], p);
}


//----------------------------------------------------------------------
// Auxiliary functions
//----------------------------------------------------------------------



template <int dim>
void FE_Q<dim>::initialize_unit_support_points ()
{
				   // number of points: (degree+1)^dim
  unsigned int n = degree+1;
  for (unsigned int i=1; i<dim; ++i)
    n *= degree+1;
  
  unit_support_points.resize(n);
  
  const double step = 1./degree;
  Point<dim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((dim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((dim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=degree; ++ix)
	{
	  p(0) = ix * step;
	  if (dim>1)
	    p(1) = iy * step;
	  if (dim>2)
	    p(2) = iz * step;
	  
	  unit_support_points[renumber[k++]] = p;
	};
};


#if deal_II_dimension == 1

template <>
void FE_Q<1>::initialize_unit_face_support_points ()
{
				   // no faces in 1d, so nothing to do
};

#endif


template <int dim>
void FE_Q<dim>::initialize_unit_face_support_points ()
{
  const unsigned int codim = dim-1;
  
				   // number of points: (degree+1)^codim
  unsigned int n = degree+1;
  for (unsigned int i=1; i<codim; ++i)
    n *= degree+1;
  
  unit_face_support_points.resize(n);
  
  const double step = 1./degree;
  Point<codim> p;
  
  unsigned int k=0;
  for (unsigned int iz=0; iz <= ((codim>2) ? degree : 0) ; ++iz)
    for (unsigned int iy=0; iy <= ((codim>1) ? degree : 0) ; ++iy)
      for (unsigned int ix=0; ix<=degree; ++ix)
	{
	  p(0) = ix * step;
	  if (codim>1)
	    p(1) = iy * step;
	  if (codim>2)
	    p(2) = iz * step;
	  
	  unit_face_support_points[face_renumber[k++]] = p;
	};
};



template <int dim>
std::vector<unsigned int>
FE_Q<dim>::get_dpo_vector(const unsigned int deg)
{
  std::vector<unsigned int> dpo(dim+1, static_cast<unsigned int>(1));
  for (unsigned int i=1; i<dpo.size(); ++i)
    dpo[i]=dpo[i-1]*(deg-1);
  return dpo;
}



template <int dim>
void
FE_Q<dim>::build_renumbering (const FiniteElementData<dim> &fe_data,
			      const unsigned int            degree,
			      std::vector<unsigned int>    &renumber)
{
  const unsigned int n = degree+1;
  
  if (degree > 0)
    for (unsigned int i=0;i<GeometryInfo<dim>::vertices_per_cell;++i)
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
	switch (dim)
	  {
	    case 1:
		  if (i==1)
		    index += degree;
		  break;
	    case 2:
		  switch (i)
		    {
		      case 1:
			    index += degree;
			    break;
		      case 3:
			    index += n*degree;
			    break;
		      case 2:
			    index += n*degree+degree;
			    break;
		    }
		  break;
	    case 3:
		  switch (i)
		    {
		      case 1:
			    index += degree;
			    break;
		      case 4:
			    index += n*degree;
			    break;
		      case 5:
			    index += n*degree+degree;
			    break;
		      case 3:
			    index += n*n*degree;
			    break;
		      case 2:
			    index += n*n*degree + degree;
			    break;
		      case 7:
			    index += n*n*degree + n*degree;
			    break;
		      case 6:
			    index += n*n*degree + n*degree+degree;
			    break;
		    }
		  break;
		  
	    default:
		  Assert(false, ExcNotImplemented());
	  }
	
	renumber[index] = i;
      }
  else
				     // degree == 0
    renumber[0] = 0;
  
				   // Lines and higher
  if (degree > 1)
    {
      for (int i=0;i< (int) GeometryInfo<dim>::lines_per_cell;++i)
	{
	  unsigned int index = fe_data.first_line_index + i*fe_data.dofs_per_line;
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
	      renumber[tensorindex] = index++;
	    }
	}

      for (int i=0;i< (int) GeometryInfo<dim>::quads_per_cell;++i)
	{
	  unsigned int index = fe_data.first_quad_index+i*fe_data.dofs_per_quad;
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
		renumber[tensorindex] = index++;
	      }
	}

      for (int i=0;i< (int) GeometryInfo<dim>::hexes_per_cell;++i)
	{
	  unsigned int index = fe_data.first_hex_index;

	  for (unsigned int jz = 1; jz<degree; jz++)
	    for (unsigned int jy = 1; jy<degree; jy++)
	      for (unsigned int jx = 1; jx<degree; jx++)
		{
		  const unsigned int tensorindex = jx + jy*n + jz*n*n;
		  renumber[tensorindex]=index++;
		}  
	}
      
    }
}



template <int dim>
void
FE_Q<dim>::build_face_renumbering (const unsigned int              degree,
				   std::vector<unsigned int>      &numbering)
{
  FiniteElementData<dim-1> fe_data(FE_Q<dim-1>::get_dpo_vector(degree),1);
  FE_Q<dim-1>::build_renumbering (fe_data, degree, numbering); 
}


#if (deal_II_dimension == 1)

template <>
void
FE_Q<1>::build_face_renumbering (const unsigned int,
				 std::vector<unsigned int>&)
{}

#endif


template <int dim>
UpdateFlags
FE_Q<dim>::update_once (const UpdateFlags flags) const
{
  UpdateFlags out = update_default;

  if (flags & update_values)
    out |= update_values;

  return out;
}



template <int dim>
UpdateFlags
FE_Q<dim>::update_each (const UpdateFlags flags) const
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
FE_Q<dim>::get_data (const UpdateFlags      update_flags,
		     const Mapping<dim>    &mapping,
		     const Quadrature<dim> &quadrature) const
{
  InternalData* data = new InternalData;
  std::vector<double> values(0);
  std::vector<Tensor<1,dim> > grads(0);
  std::vector<Tensor<2,dim> > grad_grads(0);

				   // check what needs to be
				   // initialized only once and what
				   // on every cell/face/subface we
				   // visit
  data->update_once = update_once(update_flags);
  data->update_each = update_each(update_flags);
  data->update_flags = data->update_once | data->update_each;

  const UpdateFlags flags(data->update_flags);
  const unsigned int n_q_points = quadrature.n_quadrature_points;
  
  if (flags & update_values)
    {
      values.resize (dofs_per_cell);
      data->shape_values.resize(dofs_per_cell,
				std::vector<double>(n_q_points));
    }

  if (flags & update_gradients)
    {
      grads.resize (dofs_per_cell);
      data->shape_gradients.resize(dofs_per_cell,
				   std::vector<Tensor<1,dim> >(n_q_points));
    }

				   // if second derivatives through
				   // finite differencing is required,
				   // then initialize some objects for
				   // that
  if (flags & update_second_derivatives)
    data->initialize_2nd (this, mapping, quadrature);
  
  if (flags & (update_values | update_gradients))
    for (unsigned int i=0; i<n_q_points; ++i)
      {
	poly->compute(quadrature.point(i), values, grads, grad_grads);
	for (unsigned int k=0; k<dofs_per_cell; ++k)
	  {
	    if (flags & update_values)
	      data->shape_values[renumber[k]][i] = values[k];
	    if (flags & update_gradients)
	      data->shape_gradients[renumber[k]][i] = grads[k];
	  }
      }
  return data;
}




//----------------------------------------------------------------------
// Fill data of FEValues
//----------------------------------------------------------------------

template <int dim>
void
FE_Q<dim>::fill_fe_values (const Mapping<dim>                   &mapping,
			   const typename DoFHandler<dim>::cell_iterator &cell,
			   const Quadrature<dim>                &quadrature,
			   typename Mapping<dim>::InternalDataBase &mapping_data,
			   typename Mapping<dim>::InternalDataBase &fedata,
			   FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);
  
  const UpdateFlags flags(fe_data.current_update_flags());

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0; i<quadrature.n_quadrature_points; ++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, 0);
    }

  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, 0, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_Q<dim>::fill_fe_face_values (const Mapping<dim>                   &mapping,
				const typename DoFHandler<dim>::cell_iterator &cell,
				const unsigned int                    face,
				const Quadrature<dim-1>              &quadrature,
				typename Mapping<dim>::InternalDataBase       &mapping_data,
				typename Mapping<dim>::InternalDataBase       &fedata,
				FEValuesData<dim>                    &data) const
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

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
    }

  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
void
FE_Q<dim>::fill_fe_subface_values (const Mapping<dim>                   &mapping,
				   const typename DoFHandler<dim>::cell_iterator &cell,
				   const unsigned int                    face,
				   const unsigned int                    subface,
				   const Quadrature<dim-1>              &quadrature,
				   typename Mapping<dim>::InternalDataBase       &mapping_data,
				   typename Mapping<dim>::InternalDataBase       &fedata,
				   FEValuesData<dim>                    &data) const
{
				   // convert data object to internal
				   // data for this class. fails with
				   // an exception if that is not
				   // possible
  InternalData &fe_data = dynamic_cast<InternalData &> (fedata);

				   // offset determines which data set
				   // to take (all data sets for all
				   // sub-faces are stored contiguously)
  const unsigned int offset = (face * GeometryInfo<dim>::subfaces_per_face + subface)
			      * quadrature.n_quadrature_points;

  const UpdateFlags flags(fe_data.update_once | fe_data.update_each);

  for (unsigned int k=0; k<dofs_per_cell; ++k)
    {
      for (unsigned int i=0;i<quadrature.n_quadrature_points;++i)
	if (flags & update_values)
	  data.shape_values(k,i) = fe_data.shape_values[k][i+offset];
      
      if (flags & update_gradients)
	mapping.transform_covariant(data.shape_gradients[k],
				    fe_data.shape_gradients[k],
				    mapping_data, offset);
    }
  
  if (flags & update_second_derivatives)
    compute_2nd (mapping, cell, offset, mapping_data, fe_data, data);
  
  fe_data.first_cell = false;
}



template <int dim>
unsigned int
FE_Q<dim>::memory_consumption () const
{
  Assert (false, ExcNotImplemented ());
  return 0;
}



template <int dim>
unsigned int
FE_Q<dim>::get_degree () const
{
  return degree;
};



template class FE_Q<deal_II_dimension>;
