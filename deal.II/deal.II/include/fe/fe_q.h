//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002, 2003, 2004 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------
#ifndef __deal2__fe_q_h
#define __deal2__fe_q_h

#include <base/config.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe_poly.h>


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of Lagrange finite elements @p Qp that yield the
 * finite element space of continuous, piecewise polynomials of degree
 * @p p. This class is realized using tensor product polynomials
 * based on equidistant support points.
 *
 * The constructor of this class takes the degree @p p of this finite
 * element.
 *
 * @sect3{Implementation}
 *
 * The constructor creates a TensorProductPolynomials object
 * that includes the tensor product of @p LagrangeEquidistant
 * polynomials of degree @p p. This @p TensorProductPolynomials
 * object provides all values and derivatives of the shape functions.
 *
 * Furthermore the constructor filles the @p interface_constraints,
 * the @p prolongation (embedding) and the @p restriction
 * matrices. These are implemented only up to a certain degree, that
 * is listed in the following:
 *
 * <ul>
 * <li> <tt>dim==1</tt>
 *   <ul>
 *   <li> the @p interface_constraints are not needed
 *   <li> the @p prolongation matrices up to degree 4, and
 *   <li> the @p restriction matrices up to degree 4.
 *   </ul>
 * <li> <tt>dim==2</tt>
 *   <ul>
 *   <li> the @p interface_constraints up to degree 4,
 *   <li> the @p prolongation matrices up to degree 3, and
 *   <li> the @p restriction matrices up to degree 4.
 *   </ul>
 * <li> <tt>dim==3</tt>
 *   <ul>
 *   <li> the @p interface_constraints up to degree 2,
 *   <li> the @p prolongation matrices up to degree 2, and
 *   <li> the @p restriction matrices up to degree 4.
 *   </ul>
 * </ul>
 *
 * @sect3{Numbering of the degrees of freedom (DoFs)}
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomials is a tensor product
 * numbering. However, the shape functions on a cell are renumbered
 * beginning with the shape functions whose support points are at the
 * vertices, then on the line, on the quads, and finally (for 3d) on
 * the hexes. To be explicit, these numberings are listed in the
 * following:
 *
 * @sect4{Q1 elements}
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0-------1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      3-------2
 *      |       |
 *      |       |
 *      |       |
 *      0-------1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
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
 *   @endverbatim
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   <ul>
 *   <li> Index 0: <tt>[0, 0, 0]</tt>;
 *   <li> Index 1: <tt>[1, 0, 0]</tt>;
 *   <li> Index 2: <tt>[1, 0, 1]</tt>;
 *   <li> Index 3: <tt>[0, 0, 1]</tt>;
 *   <li> Index 4: <tt>[0, 1, 0]</tt>;
 *   <li> Index 5: <tt>[1, 1, 0]</tt>;
 *   <li> Index 6: <tt>[1, 1, 1]</tt>;
 *   <li> Index 7: <tt>[0, 1, 1]</tt>;
 *   </ul>
 * </ul>
 * @sect4{Q2 elements}
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0---2---1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      3---6---2
 *      |       |
 *      7   8   5
 *      |       |
 *      0---4---1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
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
 *   @endverbatim
 *   The center vertex has number 26.
 *
 *   The respective coordinate values of the support points of the degrees
 *   of freedom are as follows:
 *   <ul>
 *   <li> Index 0: <tt>[0, 0, 0]</tt>;
 *   <li> Index 1: <tt>[1, 0, 0]</tt>;
 *   <li> Index 2: <tt>[1, 0, 1]</tt>;
 *   <li> Index 3: <tt>[0, 0, 1]</tt>;
 *   <li> Index 4: <tt>[0, 1, 0]</tt>;
 *   <li> Index 5: <tt>[1, 1, 0]</tt>;
 *   <li> Index 6: <tt>[1, 1, 1]</tt>;
 *   <li> Index 7: <tt>[0, 1, 1]</tt>;
 *   <li> Index 8: <tt>[1/2, 0, 0]</tt>;
 *   <li> Index 9: <tt>[1, 0, 1/2]</tt>;
 *   <li> Index 10: <tt>[1/2, 0, 1]</tt>;
 *   <li> Index 11: <tt>[0, 0, 1/2]</tt>;
 *   <li> Index 12: <tt>[1/2, 1, 0]</tt>;
 *   <li> Index 13: <tt>[1, 1, 1/2]</tt>;
 *   <li> Index 14: <tt>[1/2, 1, 1]</tt>;
 *   <li> Index 15: <tt>[0, 1, 1/2]</tt>;
 *   <li> Index 16: <tt>[0, 1/2, 0]</tt>;
 *   <li> Index 17: <tt>[1, 1/2, 0]</tt>;
 *   <li> Index 18: <tt>[1, 1/2, 1]</tt>;
 *   <li> Index 19: <tt>[0, 1/2, 1]</tt>;
 *   <li> Index 20: <tt>[1/2, 0, 1/2]</tt>;
 *   <li> Index 21: <tt>[1/2, 1, 1/2]</tt>;
 *   <li> Index 22: <tt>[1/2, 1/2, 0]</tt>;
 *   <li> Index 23: <tt>[1, 1/2, 1/2]</tt>;
 *   <li> Index 24: <tt>[1/2, 1/2, 1]</tt>;
 *   <li> Index 25: <tt>[0, 1/2, 1/2]</tt>;
 *   <li> Index 26: <tt>[1/2, 1/2, 1/2]</tt>; 
 *   </ul>
 * </ul>
 * @sect4{Q3 elements}
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      3--8--9--2
 *      |        |
 *      11 14 15 7
 *      |        |
 *      10 12 13 6
 *      |        |
 *      0--4--5--1
 *   @endverbatim
 *   Note the reverse ordering of degrees of freedom on the left and
 *   upper line.
 * </ul>
 * @sect4{Q4 elements}
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--4--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      3--10-11-12-2
 *      |           |
 *      15 22 23 24 9
 *      |           |
 *      14 19 20 21 8
 *      |           |
 *      13 16 17 18 7
 *      |           |
 *      0--4--5--6--1
 *   @endverbatim
 * </ul>
 * Note the reverse ordering of degrees of freedom on the left and upper
 * line.
 *
 * @author Wolfgang Bangerth, 1998, Guido Kanschat, 2001, Ralf Hartmann, 2001, 2004
 */
template <int dim>
class FE_Q : public FE_Poly<TensorProductPolynomials<dim>,dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p p.
				      */
    FE_Q (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * <tt>FE_Q<dim>(degree)</tt>, with
				      * @p dim and @p degree
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;
    
				     /**
				      * Return the matrix
				      * interpolating from the given
				      * finite element to the present
				      * one. The size of the matrix is
				      * then @p dofs_per_cell times
				      * <tt>source.dofs_per_cell</tt>.
				      *
				      * These matrices are only
				      * available if the source
				      * element is also a @p FE_Q
				      * element. Otherwise, an
				      * exception of type
				      * @ref{FiniteElementBase<dim>::ExcInterpolationNotImplemented}
				      * is thrown.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElementBase<dim> &source,
			      FullMatrix<double>           &matrix) const;
    
				     /**
				      * Check for non-zero values on a face.
				      *
				      * This function returns
				      * @p true, if the shape
				      * function @p shape_index has
				      * non-zero values on the face
				      * @p face_index.
				      *
				      * Implementation of the
				      * interface in
				      * FiniteElement
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
				      * in the files <tt>fe_q_[123]d.cc</tt>
				      * in the source directory.
				      */
    struct Matrices
    {
					 /**
					  * As the
					  * @p embedding_matrices
					  * field, but for the
					  * interface constraints. One
					  * for each element for which
					  * it has been computed.
					  */
	static const double * const constraint_matrices[];

					 /**
					  * Like
					  * @p n_embedding_matrices,
					  * but for the number of
					  * interface constraint
					  * matrices.
					  */
	static const unsigned int n_constraint_matrices;
    };

  protected:    
				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p FESystem.
				      */
    virtual FiniteElement<dim> * clone() const;

  private:
    
				     /**
				      * Only for internal use. Its
				      * full name is
				      * @p get_dofs_per_object_vector
				      * function and it creates the
				      * @p dofs_per_object vector that is
				      * needed within the constructor to
				      * be passed to the constructor of
				      * @p FiniteElementData.
				      */
    static std::vector<unsigned int> get_dpo_vector(const unsigned int degree);
    
				     /**
				      * Map tensor product data to
				      * shape function numbering. This
				      * function is actually an alike
				      * replica of the respective
				      * function in the FETools
				      * class, but is kept for three
				      * reasons:
				      *
				      * 1. It only operates on a
				      * FiniteElementData
				      * structure. This is ok in the
				      * present context, since we can
				      * control which types of
				      * arguments it is called with
				      * because this is a private
				      * function. However, the
				      * publicly visible function in
				      * the FETools class needs
				      * to make sure that the
				      * FiniteElementData object
				      * it works on actually
				      * represents a continuous finite
				      * element, which we found too
				      * difficult if we do not pass an
				      * object of type FE_Q()
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
				      *
				      * This function is called from
				      * the constructor.
				      */
    static
    std::vector<unsigned int>
    lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
					   const unsigned int            degree);

				     /**
				      * This is an analogon to the
				      * previous function, but working
				      * on faces. Called from the
				      * constructor.
				      */
    static
    std::vector<unsigned int>
    face_lexicographic_to_hierarchic_numbering (const unsigned int degree);

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
				      * @p unit_support_points field
				      * of the FiniteElementBase
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_support_points ();

				     /**
				      * Initialize the
				      * @p unit_face_support_points field
				      * of the FiniteElementBase
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_face_support_points ();
             
				     /**
				      * Mapping from hierarchic to
				      * lexicographic numbering on
				      * first face. Hierarchic is the
				      * numbering of the shape
				      * functions.
				      */
    const std::vector<unsigned int> face_index_map;
    
				     /**
				      * Allow access from other
				      * dimensions. We need this since
				      * we want to call the functions
				      * @p get_dpo_vector and
				      * @p lexicographic_to_hierarchic_numbering
				      * for the faces of the finite
				      * element of dimension dim+1.
				      */
    template <int dim1> friend class FE_Q;
};



/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <>
void FE_Q<1>::initialize_unit_face_support_points ();

template <>
std::vector<unsigned int>
FE_Q<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int);

template <>
void FE_Q<1>::initialize_constraints ();

template <>
void FE_Q<2>::initialize_constraints ();

template <>
void FE_Q<3>::initialize_constraints ();


// declaration of explicit specializations of member variables, if the
// compiler allows us to do that (the standard says we must)
#ifndef DEAL_II_MEMBER_VAR_SPECIALIZATION_BUG
template <>
const double * const FE_Q<1>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Q<1>::Matrices::n_constraint_matrices;

template <>
const double * const FE_Q<2>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Q<2>::Matrices::n_constraint_matrices;

template <>
const double * const FE_Q<3>::Matrices::constraint_matrices[];

template <>
const unsigned int FE_Q<3>::Matrices::n_constraint_matrices;

#endif

#endif
