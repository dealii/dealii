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
#ifndef __deal2__fe_q_hierarchical_h
#define __deal2__fe_q_hierarchical_h

#include <base/config.h>
#include <base/polynomial.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe.h>
#include <lac/full_matrix.h>

template <int dim> class TensorProductPolynomials;
template <int dim> class MappingQ;



/**
 * Implementation of Hierarchical finite elements @p{Qp} that yield the
 * finite element space of continuous, piecewise polynomials of degree
 * @p{p}. This class is realized using tensor product polynomials
 * based on a hierarchical basis @p{Hierarchical} of the interval 
 * @p{[0,1]} which is suitable for building an @p{hp} tensor product 
 * finite element, if we assume that each element has a single degree.
 * 
 * There are not many differences between @p{FE_Q_Hierarchical} and 
 * @p{FE_Q}, except that we add a function @p{embedding_dofs} that takes 
 * a given integer @p{q}, between @p{1} and @p{p}, and 
 * returns the numbering of basis functions of the element of order 
 * @p{q} in basis of order @p{p}.  This function is 
 * useful if one wants to make calculations using the hierarchical
 * nature of these shape functions.
 *
 * The unit support points now are reduced to @p{0}, @p{1}, and @p{0.5} in 
 * one dimension, and tensor products in higher dimensions. Thus, various 
 * interpolation functions will only work correctly for the linear case. 
 * Future work will involve writing projection--interpolation operators
 * that can interpolate onto the higher order bubble functions.
 *
 * The various constraint, prolongation, and restriction matrices are 
 * now available in all dimensions for all degrees @p{p}, currently up to 
 * order 19.
 *
 * The constructor of this class takes the degree @p{p} of this finite
 * element.
 *
 * @sect3{Implementation}
 *
 * The constructor creates a @ref{TensorProductPolynomials} object
 * that includes the tensor product of @p{Hierarchical}
 * polynomials of degree @p{p}. This @p{TensorProductPolynomials}
 * object provides all values and derivatives of the shape functions.
 *
 * @sect3{Numbering of the degrees of freedom (DoFs)}
 *
 * The original ordering of the shape functions represented by the
 * @ref{TensorProductPolynomials} is a tensor product
 * numbering. However, the shape functions on a cell are renumbered
 * beginning with the shape functions whose support points are at the
 * vertices, then on the line, on the quads, and finally (for 3d) on
 * the hexes. To be explicit, these numberings are listed in the
 * following:
 *
 * @sect4{Q1 elements}
 * @begin{itemize}
 * @item 1D case:
 *   @begin{verbatim}
 *      0-------1
 *   @end{verbatim}
 *
 * @item 2D case:
 *   @begin{verbatim}
 *      3-------2
 *      |       |
 *      |       |
 *      |       |
 *      0-------1
 *   @end{verbatim}
 *
 * @item 3D case:
 *   @begin{verbatim}
 *         7-------6        7-------6
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     3   |       |    3-------2   |
 *     |   4-------5    |       |   5
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     0-------1        0-------1
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   @begin{itemize}
 *   @item Index 0: @p{[0, 0, 0]};
 *   @item Index 1: @p{[1, 0, 0]};
 *   @item Index 2: @p{[1, 0, 1]};
 *   @item Index 3: @p{[0, 0, 1]};
 *   @item Index 4: @p{[0, 1, 0]};
 *   @item Index 5: @p{[1, 1, 0]};
 *   @item Index 6: @p{[1, 1, 1]};
 *   @item Index 7: @p{[0, 1, 1]};
 *   @end{itemize}
 * @end{itemize}
 * @sect4{Q2 elements}
 * @begin{itemize}
 * @item 1D case:
 *   @begin{verbatim}
 *      0---2---1
 *   @end{verbatim}
 *
 * @item 2D case:
 *   @begin{verbatim}
 *      3---6---2
 *      |       |
 *      7   8   5
 *      |       |
 *      0---4---1
 *   @end{verbatim}
 *
 * @item 3D case:
 *   @begin{verbatim}
 *         7--14---6        7--14---6
 *        /|       |       /       /|
 *      19 |       13     19      1813
 *      /  15      |     /       /  |
 *     3   |       |    3---10--2   |
 *     |   4--12---5    |       |   5
 *     |  /       /     |       9  /
 *    11 16      17     11      | 17
 *     |/       /       |       |/
 *     0---8---1        0---8---1
 *
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |  21   |      /  24   / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     |25 *-------*    |       |23 *
 *     |  /       /     |   20  |  /
 *     | /  22   /      |       | /
 *     |/       /       |       |/
 *     *-------*        *-------* 
 *   @end{verbatim}
 *   The center vertex has number 26.
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   @begin{itemize}
 *   @item Index 0: @p{[0, 0, 0]};
 *   @item Index 1: @p{[1, 0, 0]};
 *   @item Index 2: @p{[1, 0, 1]};
 *   @item Index 3: @p{[0, 0, 1]};
 *   @item Index 4: @p{[0, 1, 0]};
 *   @item Index 5: @p{[1, 1, 0]};
 *   @item Index 6: @p{[1, 1, 1]};
 *   @item Index 7: @p{[0, 1, 1]};
 *   @item Index 8: @p{[1/2, 0, 0]};
 *   @item Index 9: @p{[1, 0, 1/2]};
 *   @item Index 10: @p{[1/2, 0, 1]};
 *   @item Index 11: @p{[0, 0, 1/2]};
 *   @item Index 12: @p{[1/2, 1, 0]};
 *   @item Index 13: @p{[1, 1, 1/2]};
 *   @item Index 14: @p{[1/2, 1, 1]};
 *   @item Index 15: @p{[0, 1, 1/2]};
 *   @item Index 16: @p{[0, 1/2, 0]};
 *   @item Index 17: @p{[1, 1/2, 0]};
 *   @item Index 18: @p{[1, 1/2, 1]};
 *   @item Index 19: @p{[0, 1/2, 1]};
 *   @item Index 20: @p{[1/2, 0, 1/2]};
 *   @item Index 21: @p{[1/2, 1, 1/2]};
 *   @item Index 22: @p{[1/2, 1/2, 0]};
 *   @item Index 23: @p{[1, 1/2, 1/2]};
 *   @item Index 24: @p{[1/2, 1/2, 1]};
 *   @item Index 25: @p{[0, 1/2, 1/2]};
 *   @item Index 26: @p{[1/2, 1/2, 1/2]}; 
 *   @end{itemize}
 * @end{itemize}
 * @sect4{Q3 elements}
 * @begin{itemize}
 * @item 1D case:
 *   @begin{verbatim}
 *      0--2--3--1
 *   @end{verbatim}
 *
 * @item 2D case:
 *   @begin{verbatim}
 *      3--8--9--2
 *      |        |
 *      11 14 15 7
 *      |        |
 *      10 12 13 6
 *      |        |
 *      0--4--5--1
 *   @end{verbatim}
 *   Note the reverse ordering of degrees of freedom on the left and
 *   upper line.
 * @end{itemize}
 * @sect4{Q4 elements}
 * @begin{itemize}
 * @item 1D case:
 *   @begin{verbatim}
 *      0--2--3--4--1
 *   @end{verbatim}
 *
 * @item 2D case:
 *   @begin{verbatim}
 *      3--10-11-12-2
 *      |           |
 *      15 22 23 24 9
 *      |           |
 *      14 19 20 21 8
 *      |           |
 *      13 16 17 18 7
 *      |           |
 *      0--4--5--6--1
 *   @end{verbatim}
 * @end{itemize}
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line.
 *
 * @author Brian Carnes, 2002
 */
template <int dim>
class FE_Q_Hierarchical : public FiniteElement<dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p{p}.
				      */
    FE_Q_Hierarchical (const unsigned int p);
    
				     /**
				      * Return the value of the
				      * @p{i}th shape function at the
				      * point @p{p}.  @p{p} is a point
				      * on the reference element.
				      */
    virtual double shape_value (const unsigned int i,
			        const Point<dim> &p) const;
    
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
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
				      * specified component is zero.
				      */
    virtual double shape_value_component (const unsigned int i,
					  const Point<dim> &p,
					  const unsigned int component) const;

				     /**
				      * Return the gradient of the
				      * @p{i}th shape function at the
				      * point @p{p}. @p{p} is a point
				      * on the reference element, and
				      * likewise the gradient is the
				      * gradient on the unit cell with
				      * respect to unit cell
				      * coordinates.
				      */
    virtual Tensor<1,dim> shape_grad (const unsigned int  i,
				      const Point<dim>   &p) const;

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
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
				      * specified component is zero.
				      */
    virtual Tensor<1,dim> shape_grad_component (const unsigned int i,
						const Point<dim> &p,
						const unsigned int component) const;

				     /**
				      * Return the tensor of second
				      * derivatives of the @p{i}th
				      * shape function at point @p{p}
				      * on the unit cell. The
				      * derivatives are derivatives on
				      * the unit cell with respect to
				      * unit cell coordinates.
				      */
    virtual Tensor<2,dim> shape_grad_grad (const unsigned int  i,
					   const Point<dim> &p) const;

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
				      *
				      * Since this element is scalar,
				      * the returned value is the same
				      * as if the function without the
				      * @p{_component} suffix were
				      * called, provided that the
				      * specified component is zero.
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
				      * mixed discretization. Since
				      * this is a scalar element,
				      * return one.
				      */
    virtual unsigned int n_base_elements () const;
    
				     /**
				      * Access to base element
				      * objects. Since this element is
				      * scalar, @p{base_element(0)} is
				      * @p{this}, and all other
				      * indices throw an error.
				      */
    virtual const FiniteElement<dim> &
    base_element (const unsigned int index) const;
    
                                     /**
                                      * Multiplicity of base element
                                      * @p{index}. Since this is a
                                      * scalar element,
                                      * @p{element_multiplicity(0)}
                                      * returns one, and all other
                                      * indices will throw an error.
                                      */
    virtual unsigned int element_multiplicity (const unsigned int index) const;
    
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
				      * For a finite element of degree
				      * @p{sub_degree} < @p{degree}, we 
				      * return a vector which maps the 
				      * numbering on an FE
				      * of degree @p{sub_degree} into the 
				      * numbering on this element.
				      */
    std::vector<unsigned int> get_embedding_dofs (const unsigned int sub_degree) const;

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
		    const Quadrature<dim>                         &quadrature,
		    typename Mapping<dim>::InternalDataBase       &mapping_internal,
		    typename Mapping<dim>::InternalDataBase       &fe_internal,
		    FEValuesData<dim>& data) const;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_face_values (const Mapping<dim> &mapping,
			 const typename DoFHandler<dim>::cell_iterator &cell,
			 const unsigned int                            face_no,
			 const Quadrature<dim-1>                       &quadrature,
			 typename Mapping<dim>::InternalDataBase       &mapping_internal,
			 typename Mapping<dim>::InternalDataBase       &fe_internal,
			 FEValuesData<dim>& data) const ;
    
				     /**
				      * Implementation of the same
				      * function in
				      * @ref{FiniteElement}.
				      */
    virtual void
    fill_fe_subface_values (const Mapping<dim> &mapping,
			    const typename DoFHandler<dim>::cell_iterator &cell,
			    const unsigned int                            face_no,
			    const unsigned int                            sub_no,
			    const Quadrature<dim-1>                       &quadrature,
			    typename Mapping<dim>::InternalDataBase       &mapping_internal,
			    typename Mapping<dim>::InternalDataBase       &fe_internal,
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
				      * Map tensor product data to
				      * shape function numbering. This
				      * function is actually an alike
				      * replica of the respective
				      * function in the @ref{FETools}
				      * class, but is kept for three
				      * reasons:
				      *
				      * 1. It only operates on a
				      * @ref{FiniteElementData}
				      * structure. This is ok in the
				      * present context, since we can
				      * control which types of
				      * arguments it is called with
				      * because this is a private
				      * function. However, the
				      * publicly visible function in
				      * the @ref{FETools} class needs
				      * to make sure that the
				      * @ref{FiniteElementData} object
				      * it works on actually
				      * represents a continuous finite
				      * element, which we found too
				      * difficult if we do not pass an
				      * object of type @ref{FE_Q}
				      * directly.
				      *
				      * 2. If we would call the
				      * publicly available version of
				      * this function instead of this
				      * one, we would have to pass a
				      * finite element
				      * object. However, since the
				      * construction of an entire
				      * finite element object can be
				      * costly, we rather chose to
				      * retain this function.
				      *
				      * 3. Third reason is that we
				      * want to call this function for
				      * faces as well, by just calling
				      * this function for the finite
				      * element of one dimension
				      * less. If we would call the
				      * global function instead, this
				      * would require us to construct
				      * a second finite element object
				      * of one dimension less, just to
				      * call this function. Since that
				      * function does not make use of
				      * hanging nodes constraints,
				      * interpolation and restriction
				      * matrices, etc, this would have
				      * been a waste. Furthermore, it
				      * would have posed problems with
				      * template instantiations.
				      *
				      * To sum up, the existence of
				      * this function is a compromise
				      * between simplicity and proper
				      * library design, where we have
				      * chosen to weigh the simplicity
				      * aspect a little more than
				      * proper design.
				      */
    static
    void
    lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
					   const unsigned int            degree,
					   std::vector<unsigned int>    &numbering);

				     /**
				      * This is an analogon to the
				      * previous function, but working
				      * on faces.
				      */
    static
    void
    face_lexicographic_to_hierarchic_numbering (const unsigned int         degree,
						std::vector<unsigned int> &numbering);


  // not sure if needed
				     /**
				      * Initialize the
				      * @p{unit_support_points} field
				      * of the @ref{FiniteElementBase}
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_support_points ();

  // not sure if needed
				     /**
				      * Initialize the
				      * @p{unit_face_support_points} field
				      * of the @ref{FiniteElementBase}
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_face_support_points ();
    
				     /**
				      * Determine the values that need
				      * to be computed on the unit
				      * cell to be able to compute all
				      * values required by @p{flags}.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * @p{FiniteElement}.
				      *
				      * The effect in this element is
				      * as follows: if
				      * @p{update_values} is set in
				      * @p{flags}, copy it to the
				      * result. All other flags of the
				      * result are cleared, since
				      * everything else must be
				      * computed for each cell.
				      */
    virtual UpdateFlags update_once (const UpdateFlags flags) const;
  
				     /**
				      * Determine the values that need
				      * to be computed on every
				      * cell to be able to compute all
				      * values required by @p{flags}.
				      *
				      * For the purpuse of this
				      * function, refer to the
				      * documentation in
				      * @p{FiniteElement}.
				      *
				      * The effect in this element is
				      * as follows:
				      * @begin{itemize}
				      * @item if @p{update_gradients}
				      * is set, the result will
				      * contain @p{update_gradients}
				      * and
				      * @p{update_covariant_transformation}.
				      * The latter is required to
				      * transform the gradient on the
				      * unit cell to the real
				      * cell. Remark, that the action
				      * required by
				      * @p{update_covariant_transformation}
				      * is actually performed by the
				      * @p{Mapping} object used in
				      * conjunction with this finite
				      * element.
				      * @item if
				      * @p{update_second_derivatives}
				      * is set, the result will
				      * contain
				      * @p{update_second_derivatives}
				      * and
				      * @p{update_covariant_transformation}.
				      * The rationale is the same as
				      * above and no higher
				      * derivatives of the
				      * transformation are required,
				      * since we use difference
				      * quotients for the actual
				      * computation.
				      * @end{itemize}
				      */
    virtual UpdateFlags update_each (const UpdateFlags flags) const;
    
				     /**
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;

				     /**
				      * Mapping from lexicographic to
				      * shape function numbering.
				      */
    std::vector<unsigned int> renumber;

				     /**
				      * Inverse renumber
				      * vector. i.e. mapping from
				      * shape function numbering to
				      * lexicographic numbering.
				      */
    std::vector<unsigned int> renumber_inverse;
             
				     /**
				      * Mapping from lexicographic to
				      * shape function numbering on first face.
				      */
    std::vector<unsigned int> face_renumber;

                                     /**
                                      * The matrix @p{dofs_cell} contains the 
				      * values of the linear functionals of 
                                      * the master 1d cell applied to the 
				      * shape functions of the two 1d subcells.
				      * The matrix @p{dofs_subcell} constains
				      * the values of the linear functionals 
				      * on each 1d subcell applied to the 
				      * shape functions on the master 1d 
				      * subcell. 
				      * We use @p{dofs_cell} and 
				      * @p{dofs_subcell} to compute the 
				      * @p{prolongation}, @p{restriction} and 
				      * @p{interface_constraints} matrices 
				      * for all dimensions.
				      */
    std::vector<FullMatrix<double> > dofs_cell;
    std::vector<FullMatrix<double> > dofs_subcell;

				     /**
				      * Pointer to the tensor
				      * product polynomials.
				      */
    const TensorProductPolynomials<dim> polynomial_space;

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
					  * row for each shape
					  * function, containing
					  * values for each quadrature
					  * point.
					  *
					  * In this array, we store
					  * the values of the shape
					  * function in the quadrature
					  * points on the unit
					  * cell. Since these values
					  * do not change under
					  * transformation to the real
					  * cell, we only need to copy
					  * them over when visiting a
					  * concrete cell.
					  */
	Table<2,double> shape_values;

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
	Table<2,Tensor<1,dim> > shape_gradients;
    };
    
				     /**
				      * Allow access from other
				      * dimensions. We need this since
				      * we want to call the functions
				      * @p{get_dpo_vector} and
				      * @p{lexicographic_to_hierarchic_numbering}
				      * for the faces of the finite
				      * element of dimension dim+1.
				      */
    template <int dim1> friend class FE_Q_Hierarchical;
};


/* -------------- declaration of explicit specializations ------------- */

template <> void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ();
template <> void FE_Q_Hierarchical<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int,
								                   std::vector<unsigned int>&);

#endif
