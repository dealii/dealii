//---------------------------------------------------------------
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
//---------------------------------------------------------------
#ifndef __deal2__fe_nedelec_h
#define __deal2__fe_nedelec_h

#include <base/config.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe.h>

template <int dim> class TensorProductPolynomials;
template <int dim> class MappingQ;



/**
 * Implementation of continuous Nedelec elements for the space H_curl.
 *
 * The constructor of this class takes the degree @p{p} of this finite
 * element.
 *
 *
 * @sect3{Numbering of the degrees of freedom (DoFs)}
 *
 * Nedelec elements have their degrees of freedom on edges, with shape
 * functions being vector valued and pointing in tangential
 * direction. We use the standard enumeration and direction of edges
 * in deal.II, yielding the following shape functions in 2d:
 *
 *   @begin{verbatim}
 *          2
 *      *--->---*
 *      |       |
 *     3^       ^1
 *      |       |
 *      *--->---*
 *          0
 *   @end{verbatim}
 *
 * For the 3d case, the ordering follows the same scheme: the lines
 * are numbered as described in the documentation of the
 * @ref{Triangulation} class, i.e.
 *   @begin{verbatim}
 *         *---6---*        *---6---*
 *        /|       |       /       /|
 *      11 |       5      11     10 5
 *      /  7       |     /       /  |
 *     *   |       |    *---2---*   |
 *     |   *---4---*    |       |   *
 *     |  /       /     |       1  /
 *     3 8       9      3       | 9
 *     |/       /       |       |/
 *     *---0---*        *---0---*
 *   @end{verbatim}
 * and their directions are as follows:
 *   @begin{verbatim}
 *         *--->---*        *--->---*
 *        /|       |       /       /|
 *       ^ |       ^      ^       ^ ^
 *      /  ^       |     /       /  |
 *     *   |       |    *--->---*   |
 *     |   *--->---*    |       |   *
 *     |  /       /     |       ^  /
 *     ^ ^       ^      ^       | ^
 *     |/       /       |       |/
 *     *--->---*        *--->---*
 *   @end{verbatim}
 *
 * The element does not make much sense in 1d, so it is not
 * implemented there.
 *
 *
 * @author Anna Schneebeli, Wolfgang Bangerth, 2002
 */
template <int dim>
class FE_Nedelec : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor for the Nedelec
				      * element of degree @p{p}.
				      */
    FE_Nedelec (const unsigned int p);
    
				     /**
				      * Return the value of the
				      * @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;

				     /**
				      * Return the gradient of the
				      * @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Return the second derivative
				      * of the @p{component}th vector
				      * component of the @p{i}th shape
				      * function at the point
				      * @p{p}. See the
				      * @ref{FiniteElementBase} base
				      * class for more information
				      * about the semantics of this
				      * function.
				      */
    virtual Tensor<2,dim> shape_grad_grad_component (const unsigned int i,
						     const Point<dim> &p,
						     const unsigned int component) const;

				     /**
				      * Return the polynomial degree
				      * of this finite element,
				      * i.e. the value passed to the
				      * constructor.
				      */
    unsigned int get_degree () const;
    
				     /**
				      * Number of base elements in a
				      * mixed discretization. Here,
				      * this is of course equal to
				      * one.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects. Since this element is
				      * atomic, @p{base_element(0)} is
				      * @p{this}, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> & base_element (const unsigned int index) const;
    
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p{true}, if the shape
				      * function @p{shape_index} has
				      * non-zero values on the face
				      * @p{face_index}.
				      *
				      * Implementation of the
				      * interface in
				      * @ref{FiniteElement}
				      */
    virtual bool has_support_on_face (const unsigned int shape_index,
				      const unsigned int face_index) const;

				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      *
				      * This function is made virtual,
				      * since finite element objects
				      * are usually accessed through
				      * pointers to their base class,
				      * rather than the class itself.
				      */
    virtual unsigned int memory_consumption () const;


				     /**
				      * Declare a nested class which
				      * will hold static definitions of
				      * various matrices such as
				      * constraint and embedding
				      * matrices. The definition of
				      * the various static fields are
				      * in the files @p{fe_nedelec_[23]d.cc}
				      * in the source directory.
				      */
    struct Matrices
    {
					 /**
					  * Embedding matrices. For
					  * each element type (the
					  * first index) there are as
					  * many embedding matrices as
					  * there are children per
					  * cell. The first index
					  * starts with linear
					  * elements and goes up in
					  * polynomial degree. The
					  * array may grow in the
					  * future with the number of
					  * elements for which these
					  * matrices have been
					  * computed. If for some
					  * element, the matrices have
					  * not been computed then you
					  * may use the element
					  * nevertheless but can not
					  * access the respective
					  * fields.
					  */
	static const double * const
	embedding[][GeometryInfo<dim>::children_per_cell];

					 /**
					  * Number of elements (first
					  * index) the above field
					  * has. Equals the highest
					  * polynomial degree for
					  * which the embedding
					  * matrices have been
					  * computed.
					  */
	static const unsigned int n_embedding_matrices;

					 /**
					  * As the
					  * @p{embedding_matrices}
					  * field, but for the
					  * interface constraints. One
					  * for each element for which
					  * it has been computed.
					  */
	static const double * const constraint_matrices[];

					 /**
					  * Like
					  * @p{n_embedding_matrices},
					  * but for the number of
					  * interface constraint
					  * matrices.
					  */
	static const unsigned int n_constraint_matrices;
    };
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotUsefulInThisDimension);
    
  protected:    
				     /**
				      * @p{clone} function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p{FESystem}.
				      */
    virtual FiniteElement<dim> * clone() const;
  
				     /**
				      * Prepare internal data
				      * structures and fill in values
				      * independent of the cell.
				      */
    virtual
    typename Mapping<dim>::InternalDataBase *
    get_data (const UpdateFlags,
	      const Mapping<dim>& mapping,
	      const Quadrature<dim>& quadrature) const ;

				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_values (const Mapping<dim> &mapping,
		    const typename DoFHandler<dim>::cell_iterator &cell,
		    const Quadrature<dim>                &quadrature,
		    typename Mapping<dim>::InternalDataBase      &mapping_internal,
		    typename Mapping<dim>::InternalDataBase      &fe_internal,
		    FEValuesData<dim>& data) const;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                    face_no,
			 const Quadrature<dim-1>                &quadrature,
			 typename Mapping<dim>::InternalDataBase      &mapping_internal,
			 typename Mapping<dim>::InternalDataBase      &fe_internal,
			 FEValuesData<dim>& data) const ;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
			    const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int                    face_no,
			    const unsigned int                    sub_no,
			    const Quadrature<dim-1>                &quadrature,
			    typename Mapping<dim>::InternalDataBase      &mapping_internal,
			    typename Mapping<dim>::InternalDataBase      &fe_internal,
			    FEValuesData<dim>& data) const ;

  private:
    
				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p{get_dofs_per_object_vector}
				      * function and it creates the
				      * @p{dofs_per_object} vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p{FiniteElementData}.
				      */
    static std::vector<unsigned int> get_dpo_vector(const unsigned int degree);

				     /**
				      * Initialize the
				      * @p{unit_support_points} field
				      * of the @ref{FiniteElementBase}
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_support_points ();

				     /**
				      * Initialize the
				      * @p{unit_face_support_points} field
				      * of the @ref{FiniteElementBase}
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_face_support_points ();
    
				     /**
				      * Given a set of flags indicating
				      * what quantities are requested
				      * from a @p{FEValues} object,
				      * return which of these can be
				      * precomputed once and for
				      * all. Often, the values of
				      * shape function at quadrature
				      * points can be precomputed, for
				      * example, in which case the
				      * return value of this function
				      * would be the logical and of
				      * the input @p{flags} and
				      * @p{update_values}.
				      *
				      * For the present kind of finite
				      * element, this is exactly the
				      * case.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * This is the opposite to the
				      * above function: given a set of
				      * flags indicating what we want
				      * to know, return which of these
				      * need to be computed each time
				      * we visit a new cell.
				      *
				      * If for the computation of one
				      * quantity something else is
				      * also required (for example, we
				      * often need the covariant
				      * transformation when gradients
				      * need to be computed), include
				      * this in the result as well.
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
    
				     /**
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;

				     /**
				      * Fields of cell-independent data.
				      *
				      * For information about the
				      * general purpose of this class,
				      * see the documentation of the
				      * base class.
				      */
    class InternalData : public FiniteElementBase<dim>::InternalDataBase
    {
      public:
					 /**
					  * Array with shape function
					  * values in quadrature
					  * points. There is one
					  * vector for each shape
					  * function, containing
					  * values for each quadrature
					  * point. Since the shape
					  * functions are
					  * vector-valued (with as
					  * many components as there
					  * are space dimensions), the
					  * value is a tensor.
					  *
					  * In this array, we store
					  * the values of the shape
					  * function in the quadrature
					  * points on the unit
					  * cell. The transformation
					  * to the real space cell is
					  * then simply done by
					  * multiplication with the
					  * Jacobian of the mapping.
					  */
	std::vector<std::vector<Tensor<1,dim> > > shape_values;

					 /**
					  * Array with shape function
					  * gradients in quadrature
					  * points. There is one
					  * vector for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  *
					  * We store the gradients in
					  * the quadrature points on
					  * the unit cell. We then
					  * only have to apply the
					  * transformation (which is a
					  * matrix-vector
					  * multiplication) when
					  * visiting an actual cell.
					  */
	typename std::vector<typename std::vector<Tensor<2,dim> > > shape_gradients;
    };
    
				     /**
				      * Allow access from other
				      * dimensions.
				      */
    template <int dim1> friend class FE_Nedelec;
};


/* -------------- declaration of explicit specializations ------------- */

template <> void FE_Nedelec<1>::initialize_unit_face_support_points ();

// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <> 
const double * const 
FE_Nedelec<1>::Matrices::embedding[][GeometryInfo<1>::children_per_cell];

template <>
const unsigned int FE_Nedelec<1>::Matrices::n_embedding_matrices;

template <>
const double * const FE_Nedelec<1>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Nedelec<1>::Matrices::n_constraint_matrices;

template <> 
const double * const 
FE_Nedelec<2>::Matrices::embedding[][GeometryInfo<2>::children_per_cell];

template <>
const unsigned int FE_Nedelec<2>::Matrices::n_embedding_matrices;

template <>
const double * const FE_Nedelec<2>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Nedelec<2>::Matrices::n_constraint_matrices;

template <> 
const double * const 
FE_Nedelec<3>::Matrices::embedding[][GeometryInfo<3>::children_per_cell];

template <>
const unsigned int FE_Nedelec<3>::Matrices::n_embedding_matrices;

template <>
const double * const FE_Nedelec<3>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Nedelec<3>::Matrices::n_constraint_matrices;

#endif

#endif
