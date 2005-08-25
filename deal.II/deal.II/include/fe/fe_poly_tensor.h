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

/*!@addtogroup febase */
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
 * <h3>Deriving classes</h3>
 *
 * Any derived class must decide on the polynomial space to use.  This
 * polynomial space should be implemented simply as a set of vector
 * valued polynomials like PolynomialsBDM and
 * PolynomialsRaviartThomas.  In order to facilitate this
 * implementation, the basis of this space may be arbitrary.
 *
 * <h4>Determining the correct basis</h4>
 *
 * In most cases, the set of desired node values <i>N<sub>i</sub></i>
 * and the basis functions <i>v<sub>j</sub></i> will not fulfil the
 * interpolation condition <i>N<sub>i</sub>(v<sub>j</sub>) =
 * &delta;<sub>ij</sub></i>.
 *
 * The use of the membaer data #inverse_node_matrix allows to compute
 * the basis <i>v<sub>j</sub></i> automatically, after the node values
 * for ech original basis function of the polynomial space have been
 * computed.
 *
 * Therefore, the constructor of a derived class should have a
 * structure like this (example for interpolation in support points):
 *
 * @verbatim
 *  fill_support_points();
 *
 *  const unsigned int n_dofs = this->dofs_per_cell;
 *  FullMatrix<double> N(n_dofs, n_dofs);
 *
 *  for (unsigned int i=0;i<n_dofs;++i)
 *    {
 *      const Point<dim>& p = this->unit_support_point[i];
 *
 *      for (unsigned int j=0;j<n_dofs;++j)		
 *  	  for (unsigned int d=0;d<dim;++d)		
 *	    N(i,j) += node_vector[i][d]
 *                  * this->shape_value_component(j, p, d);
 *    }
 *
 *  this->inverse_node_matrix.reinit(n_dofs, n_dofs);
 *  this->inverse_node_matrix.invert(N);
 * @endverbatim
 *
 * @note The matrix #inverse_node_matrix should have dimensions zero
 * before this piece of code is executed. Only then,
 * shape_value_component() will return the raw bolynomial <i>j</i> as
 * defined in the polynomial space POLY.
 *
 * <h4>Setting the transformation</h4>
 *
 * In most cases, vector valued basis functions must be transformed
 * when mapped from the reference cell to the actual grid cell. These
 * transformations can be selected from the set MappingType and stored
 * in #mapping_type. Therefore, each constructor should contain a line
 * like:
 * @verbatim
 * this->mapping_type = this->independent_on_cartesian;
 * @endverbatim
 *
 *@see PolynomialsBDM, PolynomialsRaviartThomas
 *
 * @author Guido Kanschat, 2005
 **/
template <class POLY, int dim>
class FE_PolyTensor : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor.
				      *
				      * @arg @c degree: constructor
				      * argument for poly. May be
				      * different from @p
				      * fe_data.degree.
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

				     /**
				      * Given <tt>flags</tt>,
				      * determines the values which
				      * must be computed only for the
				      * reference cell. Make sure,
				      * that #mapping_type is set by
				      * the derived class, such that
				      * this function can operate
				      * correctly.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
				     /**
				      * Given <tt>flags</tt>,
				      * determines the values which
				      * must be computed in each cell
				      * cell. Make sure, that
				      * #mapping_type is set by the
				      * derived class, such that this
				      * function can operate
				      * correctly.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
    
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
					   /**
					    * No mapping has been
					    * selected, throw an error
					    * if needed.
					    */
	  no_mapping,
					   /**
					    * Shape functions do not
					    * depend on actual mesh
					    * cell
					    */
	  independent,
					   /**
					    * Shape functions do not
					    * depend on actual mesh
					    * cell. The mapping class
					    * must be
					    * MappingCartesian.
					    */
	  independent_on_cartesian,
					   /**
					    * Shape functions are
					    * transformed covariant.
					    */ 
	  covariant,
					   /**
					    * Shape functions are
					    * transformed
					    * contravariant.
					    */
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
    class InternalData : public FiniteElement<dim>::InternalDataBase
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
