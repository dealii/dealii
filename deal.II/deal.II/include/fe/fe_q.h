//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_continuous_h
#define __deal2__fe_continuous_h

#include <base/polynomial.h>

template <int dim> class TensorProductPolynomials;

/**
 * Tensor product elements based on equitdistant support points.
 */
template <int dim>
class FE_Q : public FiniteElement<dim>
{
public:
				   /**
				    * Constructor for tensor product
				    * polynomials of degree @p{k}.
				    */
  FE_Q (unsigned int k);
				   /**
				    * Destructor.
				    */
  ~FE_Q ();
  
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual FEValuesBase<dim>::InternalData*
    get_data (const UpdateFlags,
	      const Quadrature<dim>& quadrature) const ;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of faces and fill in values
				      * independent of the cell.
				      */
    virtual FEValuesBase<dim>::InternalData*
    get_face_data (const UpdateFlags flags,
		   const Quadrature<dim>& quadrature) const ;

				     /**
				      * Prepare internal data
				      * structure for transformation
				      * of children of faces and fill
				      * in values independent of the
				      * cell.
				      */
    virtual FEValuesBase<dim>::InternalData*
    get_subface_data (const UpdateFlags flags,
		       const Quadrature<dim>& quadrature) const;

    virtual void
    fill_fe_values (const Mapping<dim> &mapping,
		   const DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    FEValuesBase<dim>::InternalData      &mapping_internal,
		    FEValuesBase<dim>::InternalData      &fe_internal,
		    FEValuesData<dim>& data) const;
    
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim>                &quadrature,
			 FEValuesBase<dim>::InternalData      &mapping_internal,
			 FEValuesBase<dim>::InternalData      &fe_internal,
			 FEValuesData<dim>& data) const ;
    
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
		   const DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim>                &quadrature,
			    FEValuesBase<dim>::InternalData      &mapping_internal,
			    FEValuesBase<dim>::InternalData      &fe_internal,
			    FEValuesData<dim>& data) const ;

private:
				   /**
				    * Map tensor product data to shape
				    * function numbering.
				    *
				    * The node values are ordered such
				    * that vertices are first,
				    * followed by lines,
				    * quadrilaterals and
				    * hexahedra. Furthermore, the
				    * ordering inside each group may
				    * be confused, too. Therefore,
				    * this function computes a mapping
				    * from lexicographic ordering
				    * (x,y,z) to the shape function
				    * structure.
				    */
  void build_renumbering (unsigned int degree,
			  vector<unsigned int>& numbering);

				   /**
				    * Compute flags for initial update only.
				    */
  static UpdateFlags update_once (UpdateFlags flags);
  
				   /**
				    * Compute flags for update on each cell.
				    */
  static UpdateFlags update_each (UpdateFlags flags);
  
				   /**
				    * Degree of the polynomials.
				    */  
  const unsigned int degree;
				   /**
				    * Mapping from lexicographic to
				    * shape function numbering.
				    */
  vector<unsigned int> renumber;
				   /**
				    * Fields of cell-independent data.
				    */
  class InternalData : public FEValuesBase<dim>::InternalData
  {
  public:
    vector<double> shape_values;
    vector<Tensor<1,dim> > shape_grads;
  };
  
				   /**
				    * Vector of one-dimensional
				    * polynomials used.
				    */
  vector<LagrangeEquidistant> polynomials;

				   /**
				    * Implementation of the tensor
				    * product of polynomials.
				    */
  TensorProductPolynomials<dim>* poly;
};



#endif
