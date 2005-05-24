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
#ifndef __deal2__fe_poly_tensor_h
#define __deal2__fe_poly_tensor_h


#include <lac/full_matrix.h>
#include <fe/fe.h>

/*!@addtogroup fe */
/*@{*/

/**
 * This class gives a unified framework for the implementation of
 * FiniteElement classes based on Tensor valued polynomial spaces like
 * PolynomialsBDM and PolynomialsRaviartThomas.
 *
 * Every class that implements following function can be used as
 * template parameter POLY.
 *
 * @code
 * void compute (const Point<dim>            &unit_point,
 *               std::vector<Tensor<1,dim> > &values,
 *               std::vector<Tensor<2,dim> > &grads,
 *               std::vector<Tensor<3,dim> > &grad_grads) const;
 * @endcode
 *
 * The polynomial spaces are usually described as direct sum of
 * simpler spaces. In such a case, the usual basis of node functionals
 * is not dual to the basis of the polynomial space. Therefore, the
 * matrix #inverse_node_matrix can be filled by the constructor of a
 * derived class such that the usual interpolation condition
 * <i>N<sub>i</sub>(v<sub>j</sub>)</i> holds on the reference cell.
 *
 * In many cases, the node functionals depend on the shape of the mesh
 * cell, since they evaluate normal or tangential components on the
 * faces. In order to allow for a set of transformations, the variable
 * #mapping_type has been introduced. It should also be set in the
 * constructor of a derived class.
 *
 * This class is not a fully implemented FiniteElement class, but
 * implements some common features of vector valued elements based on
 * vector valued polynomial classes. What's missing here in particular
 * is information on the topological location of the node values.
 *
 * @see PolynomialsBDM, PolynomialsRaviartThomas
 *
 * @author Guido Kanschat, 2005
 **/
template <class POLY, int dim>
class FE_PolyTensor : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor.
				      */
    FE_PolyTensor (unsigned int degree,
		   const FiniteElementData<dim> &fe_data,
		   const std::vector<bool> &restriction_is_additive_flags,
		   const std::vector<std::vector<bool> > &nonzero_components);

				     /**
				      * Since these elements are
				      * vector valued, an exception is
				      * thrown.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;

				     /**
				      * Since these elements are
				      * vector valued, an exception is
				      * thrown.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Since these elements are
				      * vector valued, an exception is
				      * thrown.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
						     const Point<dim> &p,
						     const unsigned int component) const;

                                     /**
				      * Number of base elements in a
				      * mixed discretization. Since
				      * this is not a composed element,
				      * return one.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects. Since this element is
				      * not composed of several elements,
				      * <tt>base_element(0)</tt> is
				      * <tt>this</tt>, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * Multiplicity of base element
                                      * <tt>index</tt>. Since this is
                                      * not a composed element,
                                      * <tt>element_multiplicity(0)</tt>
                                      * returns one, and all other
                                      * indices will throw an error.
                                      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;

    
  protected:
				     /**
				      * Different options for
				      * transforming the basis
				      * functions from the reference
				      * cell to the actual mesh cell.
				      *
				      * Most vector valued elements
				      * either transform shape
				      * functions to keep node values
				      * on edges meaningful. Still, in
				      * special cases, it may be
				      * possible to avoid the mapping.
				      */
    enum MappingType {
					   /// Shape functions do not depend on actual mesh cell
	  independent,
					   /// Shape functions are transformed covariant.
	  covariant,
					   /// Shape functions are transformed contravariant.
	  contravariant
    };

				     /**
				      * The mapping type to be used to
				      * map shape functions from the
				      * reference cell to the msh
				      * cell.
				      */
    MappingType mapping_type;
    
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags,
	      const Mapping<dim>& mapping,
	      const Quadrature<dim>& quadrature) const ;

    virtual void
    fill_fe_values (const Mapping<dim> &mapping,
		    const typename Triangulation<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase      &mapping_internal,
		    typename Mapping<dim>::InternalDataBase      &fe_internal,
		    FEValuesData<dim>& data) const;
    
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const typename Triangulation<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>                &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_internal,
			 typename Mapping<dim>::InternalDataBase      &fe_internal,
			 FEValuesData<dim>& data) const ;
    
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
			    const typename Triangulation<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>                &quadrature,
			    typename Mapping<dim>::InternalDataBase      &mapping_internal,
			    typename Mapping<dim>::InternalDataBase      &fe_internal,
			    FEValuesData<dim>& data) const ;

				     /**
				      * Fields of cell-independent
				      * data for FE_PolyTensor. Stores
				      * the values of the shape
				      * functions and their
				      * derivatives on the reference
				      * cell for later use.
				      *
				      * All tables are organized in a
				      * way, that the value for shape
				      * function <i>i</i> at
				      * quadrature point <i>k</i> is
				      * accessed by indices
				      * <i>(i,k)</i>.
				      */
    class InternalData : public FiniteElementBase<dim>::InternalDataBase
    {
      public:
					 /**
					  * Array with shape function
					  * values in quadrature
					  * points. There is one
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  */
	std::vector<std::vector<Tensor<1,dim> > > shape_values;

					 /**
					  * Array with shape function
					  * gradients in quadrature
					  * points. There is one
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  */      
	std::vector<std::vector<Tensor<2,dim> > > shape_grads;
    };
    
				     /**
				      * Degree of the polynomials.
				      */  
    unsigned int degree;

                                     /**
                                      * The polynomial space. Its type
                                      * is given by the template
                                      * parameter POLY.
                                      */    
    POLY poly_space;
				     /**
				      * The inverse of the matrix
				      * <i>a<sub>ij</sub></i> of node
				      * values <i>N<sub>i</sub></i>
				      * applied to polynomial
				      * <i>p<sub>j</sub></i>. This
				      * matrix is used to convert
				      * polynomials in the "raw" basis
				      * provided in #poly_space to the
				      * basis dual to the node
				      * functionals on the reference cell.
				      *
				      * This object is not filled by
				      * FE_PolyTensor, but is a chance
				      * for a derived class to allow
				      * for reorganization of the
				      * basis functions. If it is left
				      * empty, the basis in
				      * #poly_space is used.
				      */
    FullMatrix<double> inverse_node_matrix;

				     /**
				      * If a shape function is
				      * computed at a single point, we
				      * must compute all of them to
				      * apply #inverse_node_matrix. In
				      * order to avoid too much
				      * overhead, we cache the point
				      * and the function values for
				      * the next evaluation.
				      */ 
    mutable Point<dim> cached_point;
    
				     /**
				      * Cached shape function values after
				      * call to
				      * shape_value_component().
				      */
    mutable std::vector<Tensor<1,dim> > cached_values;
    
				     /**
				      * Cached shape function gradients after
				      * call to
				      * shape_grad_component().
				      */
    mutable std::vector<Tensor<2,dim> > cached_grads;
    
				     /**
				      * Cached second derivatives of
				      * shape functions after call to
				      * shape_grad_grad_component().
				      */
    mutable std::vector<Tensor<3,dim> > cached_grad_grads;
};

/*@}*/

#endif
