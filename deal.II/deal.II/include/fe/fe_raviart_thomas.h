//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_raviart_thomas_h
#define __deal2__fe_raviart_thomas_h

#include <base/config.h>
#include <base/tensor_product_polynomials.h>
#include <grid/geometry_info.h>
#include <fe/fe.h>

template <int dim> class MappingQ;



/**
 * Implementation of continuous Raviart-Thomas elements for the space
 * H_div. Note, however, that continuity only concerns the normal
 * component of the vector field.
 *
 * The constructor of this class takes the degree @p{p} of this finite
 * element. The numbering of the degree of this element in the
 * literature is somewhat funny: the degree is defined not as the
 * polynomial degree of the finite element space, but as that of the
 * normal component of the traces onto the boundary. Thus, the lowest
 * order, zero, has linear shape functions, but on the faces, the
 * traces of the normal component of these elements is constant on
 * each face.
 * 
 * 
 * @sect3{Interpolation to finer and coarser meshes}
 *
 * Each finite element class in deal.II provides matrices that are
 * used to interpolate from coarser to finer meshes and the other way
 * round. Interpolation from a mother cell to its children is usually
 * trivial, since finite element spaces are normally nested and this
 * kind of interpolation is therefore exact. On the other hand, when
 * we interpolate from child cells to the mother cell, we usually have
 * to throw away some information.
 *
 * For continuous elements, this transfer usually happens by
 * interpolating the values on the child cells at the support points
 * of the shape functions of the mother cell. However, for
 * discontinuous elements, we often use a projection from the child
 * cells to the mother cell. The projection approach is only possible
 * for discontinuous elements, since it cannot be guaranteed that the
 * values of the projected functions on one cell and its neighbor
 * match. In this case, only an interpolation can be
 * used. (Internally, whether the values of a shape function are
 * interpolated or projected, or better: whether the matrices the
 * finite element provides are to be treated with the properties of a
 * projection or of an interpolation, is controlled by the
 * @p{restriction_is_additive} flag. See there for more information.)
 *
 * Here, things are not so simple: since the element has some
 * continuity requirements across faces, we can only resort to some
 * kind of interpolation. On the other hand, for the lowest order
 * elements, the values of generating functionals are the (constant)
 * tangential values of the shape functions. We would therefore really
 * like to take the mean value of the tangential values of the child
 * faces, and make this the value of the mother face. Then, however,
 * taking a mean value of two piecewise constant function is not an
 * interpolation, but a restriction. Since this is not possible, we
 * cannot use this.
 *
 * To make a long story somewhat shorter, when interpolating from
 * refined edges to a coarse one, we do not take the mean value, but
 * pick only one (the one from the first child edge). While this is
 * not optimal, it is certainly a valid choice (using an interpolation
 * point that is not in the middle of the cell, but shifted to one
 * side), and it also preserves the order of the interpolation.
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
 *      *---^---*
 *      |       |
 *     3>       >1
 *      |       |
 *      *---^---*
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
 * @author Wolfgang Bangerth, 2003
 */
template <int dim>
class FE_RaviartThomas : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor for the Nedelec
				      * element of degree @p{p}.
				      */
    FE_RaviartThomas (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * @p{FE_RaviartThomas<dim>(degree)}, with
				      * @p{dim} and @p{degree}
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;

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
				      * Return the matrix
				      * interpolating from the given
				      * finite element to the present
				      * one. The size of the matrix is
				      * then @p{dofs_per_cell} times
				      * @p{source.dofs_per_cell}.
				      *
				      * These matrices are only
				      * available if the source
				      * element is also a Raviart
				      * Thomas element. Otherwise, an
				      * exception of type
				      * @ref{FiniteElementBase<dim>::ExcInterpolationNotImplemented}
				      * is thrown.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElementBase<dim> &source,
			      FullMatrix<double>           &matrix) const;

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
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;

                                     /**
                                      * Multiplicity of base element
                                      * @p{index}. Since this is an
                                      * atomic element,
                                      * @p{element_multiplicity(0)}
                                      * returns one, and all other
                                      * indices will throw an error.
                                      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;
    
				     /**
				      * This function returns
				      * @p{true}, if the shape
				      * function @p{shape_index} has
				      * non-zero values on the face
				      * @p{face_index}. For the lowest
				      * order Nedelec elements, this
				      * is actually the case for the
				      * one on which the shape
				      * function is defined and all
				      * neighboring ones.
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
			 FEValuesData<dim>& data) const;
    
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
			    FEValuesData<dim>& data) const;

  private:
				     /**
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;

                                     /**
                                      * Spaces describing the
                                      * anisotropic polynomial spaces
                                      * for each vector component,
                                      * i.e. there are @p{dim}
                                      * elements of this field. The
                                      * values for this member are
                                      * created in
                                      * @ref{create_polynomials}.
                                      */
    const std::vector<AnisotropicPolynomials<dim> > polynomials;

                                     /**
                                      * For each shape function, store
                                      * to which vector component (on
                                      * the unit cell, they are mixed
                                      * on the real cell by the
                                      * transformation) they belong,
                                      * and which index they have
                                      * within the anisotropic tensor
                                      * product polynomial space
                                      * describing this vector
                                      * component.
                                      *
                                      * These values are computed by
                                      * the @ref{compute_renumber}
                                      * function.
                                      */
    const std::vector<std::pair<unsigned int, unsigned int> > renumber;
    
    
                                     /**
                                      * Generate the polynomial spaces
                                      * for the @ref{polynomials}
                                      * member.
                                      */
    static std::vector<AnisotropicPolynomials<dim> >
    create_polynomials (const unsigned int degree);
    
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
    static std::vector<unsigned int>
    get_dpo_vector (const unsigned int degree);

				     /**
				      * Compute the vector used for
				      * the
				      * @p{restriction_is_additive}
				      * field passed to the base
				      * class's constructor.
				      */
    static std::vector<bool>
    get_ria_vector (const unsigned int degree);

                                     /**
                                      * Compute the values of the
                                      * @p{renumber} field.
                                      */
    static std::vector<std::pair<unsigned int, unsigned int> >
    compute_renumber (const unsigned int);

				     /**
				      * Initialize the hanging node
				      * constraints matrices. Called
				      * from the constructor.
				      */
    void initialize_constraints ();

				     /**
				      * Initialize the embedding
				      * matrices. Called from the
				      * constructor.
				      */
    void initialize_embedding ();

				     /**
				      * Initialize the restriction
				      * matrices. Called from the
				      * constructor.
				      */
    void initialize_restriction ();
    
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
					  * points. There is one row
					  * for each shape function,
					  * containing values for each
					  * quadrature point. Since
					  * the shape functions are
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
	Table<2,Tensor<1,dim> > shape_values;

					 /**
					  * Array with shape function
					  * gradients in quadrature
					  * points. There is one
					  * row for each shape
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
	Table<2,Tensor<2,dim> > shape_gradients;
    };
    
				     /**
				      * Allow access from other
				      * dimensions.
				      */
    template <int dim1> friend class FE_RaviartThomas;
};


/* -------------- declaration of explicit specializations ------------- */

template <>
void FE_RaviartThomas<1>::initialize_unit_face_support_points ();

template <>
std::vector<unsigned int> FE_RaviartThomas<1>::get_dpo_vector (const unsigned int);

template <>
std::vector<AnisotropicPolynomials<1> >
FE_RaviartThomas<1>::create_polynomials (const unsigned int);

template <>
std::vector<AnisotropicPolynomials<2> >
FE_RaviartThomas<2>::create_polynomials (const unsigned int);

template <>
std::vector<AnisotropicPolynomials<3> >
FE_RaviartThomas<3>::create_polynomials (const unsigned int);

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<1>::compute_renumber (const unsigned int);

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<2>::compute_renumber (const unsigned int);

template <>
std::vector<std::pair<unsigned int, unsigned int> >
FE_RaviartThomas<3>::compute_renumber (const unsigned int);

template <>
void
FE_RaviartThomas<1>::initialize_constraints ();

template <>
void
FE_RaviartThomas<2>::initialize_constraints ();

template <>
void
FE_RaviartThomas<3>::initialize_constraints ();

template <>
void
FE_RaviartThomas<1>::initialize_embedding ();

template <>
void
FE_RaviartThomas<1>::initialize_restriction ();

template <>
void
FE_RaviartThomas<2>::initialize_restriction ();

template <>
void
FE_RaviartThomas<3>::initialize_restriction ();


#endif
