//---------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2002, 2003, 2004 by the deal.II authors
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
#include <base/tensor_product_polynomials.h>
#include <fe/fe_poly.h>
#include <fe/fe.h>
//#include <lac/full_matrix.h>

template <int dim> class TensorProductPolynomials;
template <int dim> class MappingQ;


/*!@addtogroup fe */
/*@{*/

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
 * @author Brian Carnes, 2002, Ralf Hartmann 2004
 */
template <int dim>
class FE_Q_Hierarchical : public FE_Poly<TensorProductPolynomials<dim>,dim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p p.
				      */
    FE_Q_Hierarchical (const unsigned int p);
    
				     /**
				      * Return a string that uniquely
				      * identifies a finite
				      * element. This class returns
				      * @p{FE_Q_Hierarchical<dim>(degree)},
				      * with @p{dim} and @p{degree}
				      * replaced by appropriate
				      * values.
				      */
    virtual std::string get_name () const;

				     /**
				      * Return the polynomial degree
				      * of this finite element,
				      * i.e. the value passed to the
				      * constructor.
				      */
    unsigned int get_degree () const;
    
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
    std::vector<unsigned int>
    lexicographic_to_hierarchic_numbering (const FiniteElementData<dim> &fe_data,
					   const unsigned int            degree);

				     /**
				      * This is an analogon to the
				      * previous function, but working
				      * on faces.
				      */
    static
    std::vector<unsigned int>
    face_lexicographic_to_hierarchic_numbering (const unsigned int degree);

				     /**
				      * Initialize two auxiliary
				      * fields that will be used in
				      * setting up the various
				      * matrices in the constructor.
				      */
    void build_dofs_cell (std::vector<FullMatrix<double> > &dofs_cell,
			  std::vector<FullMatrix<double> > &dofs_subcell) const;
    
				     /**
				      * Initialize the hanging node
				      * constraints matrices. Called
				      * from the constructor.
				      */
    void initialize_constraints (const std::vector<FullMatrix<double> > &dofs_subcell);

				     /**
				      * Initialize the embedding
				      * matrices. Called from the
				      * constructor.
				      */
    void initialize_embedding_and_restriction (const std::vector<FullMatrix<double> > &dofs_cell,
					       const std::vector<FullMatrix<double> > &dofs_subcell);

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
				      * Degree of the polynomials.
				      */  
    const unsigned int degree;
             
				     /**
				      * Mapping from lexicographic to
				      * shape function numbering on first face.
				      */
    const std::vector<unsigned int> face_renumber;
    
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

/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <>
void FE_Q_Hierarchical<1>::initialize_unit_face_support_points ();

template <>
std::vector<unsigned int>
FE_Q_Hierarchical<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int);

#endif
