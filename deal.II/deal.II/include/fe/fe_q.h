//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__fe_q_h
#define __deal2__fe_q_h

#include <base/config.h>
#include <base/tensor_product_polynomials.h>
#include <fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN


/*!@addtogroup fe */
/*@{*/

/**
 * Implementation of a scalar Lagrange finite element @p Qp that yields the
 * finite element space of continuous, piecewise polynomials of degree @p p in
 * each coordinate direction. This class is realized using tensor product
 * polynomials based on equidistant or given support points.
 *
 * The standard constructor of this class takes the degree @p p of this finite
 * element. Alternatively, it can take a quadrature formula @p points defining
 * the support points of the Lagrange interpolation in one coordinate direction.
 *
 * For more information about the <tt>spacedim</tt> template parameter
 * check the documentation of FiniteElement or the one of
 * Triangulation.
 *
 * <h3>Implementation</h3>
 *
 * The constructor creates a TensorProductPolynomials object that includes the
 * tensor product of @p LagrangeEquidistant polynomials of degree @p p. This
 * @p TensorProductPolynomials object provides all values and derivatives of
 * the shape functions.  In case a quadrature rule is given, the constructure
 * creates a TensorProductPolynomials object that includes the tensor product
 * of @p Lagrange polynomials with the support points from @p points.
 *
 * Furthermore the constructor fills the @p interface_constraints, the
 * @p prolongation (embedding) and the @p restriction matrices. These
 * are implemented only up to a certain degree and may not be
 * available for very high polynomial degree.
 *
 *
 * <h3>Numbering of the degrees of freedom (DoFs)</h3>
 *
 * The original ordering of the shape functions represented by the
 * TensorProductPolynomials is a tensor product
 * numbering. However, the shape functions on a cell are renumbered
 * beginning with the shape functions whose support points are at the
 * vertices, then on the line, on the quads, and finally (for 3d) on
 * the hexes. To be explicit, these numberings are listed in the
 * following:
 *
 * <h4>Q1 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0-------1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2-------3
 *      |       |
 *      |       |
 *      |       |
 *      0-------1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
 *         6-------7        6-------7
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     4   |       |    4-------5   |
 *     |   2-------3    |       |   3
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
 *   <li> Index 2: <tt>[0, 1, 0]</tt>;
 *   <li> Index 3: <tt>[1, 1, 0]</tt>;
 *   <li> Index 4: <tt>[0, 0, 1]</tt>;
 *   <li> Index 5: <tt>[1, 0, 1]</tt>;
 *   <li> Index 6: <tt>[0, 1, 1]</tt>;
 *   <li> Index 7: <tt>[1, 1, 1]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q2 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0---2---1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2---7---3
 *      |       |
 *      4   8   5
 *      |       |
 *      0---6---1
 *   @endverbatim
 *
 * <li> 3D case:
 *   @verbatim
 *         6--15---7        6--15---7
 *        /|       |       /       /|
 *      12 |       19     12      1319
 *      /  18      |     /       /  |
 *     4   |       |    4---14--5   |
 *     |   2---11--3    |       |   3
 *     |  /       /     |      17  /
 *    16 8       9     16       | 9
 *     |/       /       |       |/
 *     0---10--1        0---8---1
 *
 *         *-------*        *-------*
 *        /|       |       /       /|
 *       / |  23   |      /  25   / |
 *      /  |       |     /       /  |
 *     *   |       |    *-------*   |
 *     |20 *-------*    |       |21 *
 *     |  /       /     |   22  |  /
 *     | /  24   /      |       | /
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
 *   <li> Index 2: <tt>[0, 1, 0]</tt>;
 *   <li> Index 3: <tt>[1, 1, 0]</tt>;
 *   <li> Index 4: <tt>[0, 0, 1]</tt>;
 *   <li> Index 5: <tt>[1, 0, 1]</tt>;
 *   <li> Index 6: <tt>[0, 1, 1]</tt>;
 *   <li> Index 7: <tt>[1, 1, 1]</tt>;
 *   <li> Index 8: <tt>[0, 1/2, 0]</tt>;
 *   <li> Index 9: <tt>[1, 1/2, 0]</tt>;
 *   <li> Index 10: <tt>[1/2, 0, 0]</tt>;
 *   <li> Index 11: <tt>[1/2, 1, 0]</tt>;
 *   <li> Index 12: <tt>[0, 1/2, 1]</tt>;
 *   <li> Index 13: <tt>[1, 1/2, 1]</tt>;
 *   <li> Index 14: <tt>[1/2, 0, 1]</tt>;
 *   <li> Index 15: <tt>[1/2, 1, 1]</tt>;
 *   <li> Index 16: <tt>[0, 0, 1/2]</tt>;
 *   <li> Index 17: <tt>[1, 0, 1/2]</tt>;
 *   <li> Index 18: <tt>[0, 1, 1/2]</tt>;
 *   <li> Index 19: <tt>[1, 1, 1/2]</tt>;
 *   <li> Index 20: <tt>[0, 1/2, 1/2]</tt>;
 *   <li> Index 21: <tt>[1, 1/2, 1/2]</tt>;
 *   <li> Index 22: <tt>[1/2, 0, 1/2]</tt>;
 *   <li> Index 23: <tt>[1/2, 1, 1/2]</tt>;
 *   <li> Index 24: <tt>[1/2, 1/2, 0]</tt>;
 *   <li> Index 25: <tt>[1/2, 1/2, 1]</tt>;
 *   <li> Index 26: <tt>[1/2, 1/2, 1/2]</tt>;
 *   </ul>
 * </ul>
 * <h4>Q3 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2--10-11-3
 *      |        |
 *      5  14 15 7
 *      |        |
 *      4  12 13 6
 *      |        |
 *      0--8--9--1
 *   @endverbatim
 * </ul>
 * <h4>Q4 elements</h4>
 * <ul>
 * <li> 1D case:
 *   @verbatim
 *      0--2--3--4--1
 *   @endverbatim
 *
 * <li> 2D case:
 *   @verbatim
 *      2--13-14-15-3
 *      |           |
 *      6  22 23 24 9
 *      |           |
 *      5  19 20 21 8
 *      |           |
 *      4  16 17 18 7
 *      |           |
 *      0--10-11-12-1
 *   @endverbatim
 * </ul>
 *
 * @author Wolfgang Bangerth, 1998, 2003; Guido Kanschat, 2001; Ralf Hartmann, 2001, 2004, 2005; Oliver Kayser-Herold, 2004; Katharina Kormann, 2008; Martin Kronbichler, 2008
 */
template <int dim, int spacedim=dim>
class FE_Q : public FE_Poly<TensorProductPolynomials<dim>,dim,spacedim>
{
  public:
				     /**
				      * Constructor for tensor product
				      * polynomials of degree @p p.
				      */
    FE_Q (const unsigned int p);

				     /**
				      * Constructor for tensor product
				      * polynomials with support points @p
				      * points based on a one-dimensional
				      * quadrature formula. The degree of the
				      * finite element is
				      * <tt>points.size()-1</tt>.  Note that
				      * the first point has to be 0 and the
				      * last one 1.
				      */

    FE_Q (const Quadrature<1> &points);
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
				      * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented
				      * is thrown.
				      */
    virtual void
    get_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
			      FullMatrix<double>       &matrix) const;


				     /**
				      * Return the matrix
				      * interpolating from a face of
				      * of one element to the face of
				      * the neighboring element.
				      * The size of the matrix is
				      * then <tt>source.dofs_per_face</tt> times
				      * <tt>this->dofs_per_face</tt>.
				      *
				      * Derived elements will have to
				      * implement this function. They
				      * may only provide interpolation
				      * matrices for certain source
				      * finite elements, for example
				      * those from the same family. If
				      * they don't implement
				      * interpolation from a given
				      * element, then they must throw
				      * an exception of type
				      * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
				      */
    virtual void
    get_face_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
				   FullMatrix<double>       &matrix) const;

				     /**
				      * Return the matrix
				      * interpolating from a face of
				      * of one element to the face of
				      * the neighboring element.
				      * The size of the matrix is
				      * then <tt>source.dofs_per_face</tt> times
				      * <tt>this->dofs_per_face</tt>.
				      *
				      * Derived elements will have to
				      * implement this function. They
				      * may only provide interpolation
				      * matrices for certain source
				      * finite elements, for example
				      * those from the same family. If
				      * they don't implement
				      * interpolation from a given
				      * element, then they must throw
				      * an exception of type
				      * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented.
				      */
    virtual void
    get_subface_interpolation_matrix (const FiniteElement<dim,spacedim> &source,
				      const unsigned int        subface,
				      FullMatrix<double>       &matrix) const;

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
				      * @name Functions to support hp
				      * @{
				      */

                                     /**
                                      * Return whether this element
                                      * implements its hanging node
                                      * constraints in the new way,
				      * which has to be used to make
				      * elements "hp compatible".
                                      *
				      * For the FE_Q class the result is
				      * always true (independent of the degree
				      * of the element), as it implements the
				      * complete set of functions necessary
				      * for hp capability.
                                      */
    virtual bool hp_constraints_are_implemented () const;

				     /**
				      * If, on a vertex, several
				      * finite elements are active,
				      * the hp code first assigns the
				      * degrees of freedom of each of
				      * these FEs different global
				      * indices. It then calls this
				      * function to find out which of
				      * them should get identical
				      * values, and consequently can
				      * receive the same global DoF
				      * index. This function therefore
				      * returns a list of identities
				      * between DoFs of the present
				      * finite element object with the
				      * DoFs of @p fe_other, which is
				      * a reference to a finite
				      * element object representing
				      * one of the other finite
				      * elements active on this
				      * particular vertex. The
				      * function computes which of the
				      * degrees of freedom of the two
				      * finite element objects are
				      * equivalent, and returns a list
				      * of pairs of global dof indices
				      * in @p identities. The first
				      * index of each pair denotes one
				      * of the vertex dofs of the
				      * present element, whereas the
				      * second is the corresponding
				      * index of the other finite
				      * element.
				      */
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_vertex_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

				     /**
				      * Same as
				      * hp_vertex_dof_indices(),
				      * except that the function
				      * treats degrees of freedom on
				      * lines.
				      */
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_line_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

				     /**
				      * Same as
				      * hp_vertex_dof_indices(),
				      * except that the function
				      * treats degrees of freedom on
				      * quads.
				      */
    virtual
    std::vector<std::pair<unsigned int, unsigned int> >
    hp_quad_dof_identities (const FiniteElement<dim,spacedim> &fe_other) const;

				     /**
				      * Return whether this element dominates
				      * the one given as argument when they
				      * meet at a common face,
				      * whether it is the other way around,
				      * whether neither dominates, or if
				      * either could dominate.
				      *
				      * For a definition of domination, see
				      * FiniteElementBase::Domination and in
				      * particular the @ref hp_paper "hp paper".
				      */
    virtual
    FiniteElementDomination::Domination
    compare_for_face_domination (const FiniteElement<dim,spacedim> &fe_other) const;
				     //@}

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

  protected:
				     /**
				      * @p clone function instead of
				      * a copy constructor.
				      *
				      * This function is needed by the
				      * constructors of @p FESystem.
				      */
    virtual FiniteElement<dim,spacedim> * clone() const;

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
				      * This is an analogon to the
				      * FETools::lexicographic_to_hierarchic_numbering
				      * function, but working on
				      * faces. Called from the
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
				      * of the FiniteElement
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_support_points ();

					     /**
				      * Initialize the @p unit_support_points
				      * field of the FiniteElement
				      * class. Called from the constructor in
				      * case the finite element is based on
				      * quadrature points.
				      */
    void initialize_unit_support_points (const Quadrature<1> &points);

				     /**
				      * Initialize the
				      * @p unit_face_support_points field
				      * of the FiniteElement
				      * class. Called from the
				      * constructor.
				      */
    void initialize_unit_face_support_points ();

					     /**
				      * Initialize the @p
				      * unit_face_support_points field of the
				      * FiniteElement class. Called from the
				      * constructor in case the finite element
				      * is based on quadrature points.
				      */
    void initialize_unit_face_support_points (const Quadrature<1> &points);

				     /**
				      * Initialize the
				      * @p adjust_quad_dof_index_for_face_orientation_table field
				      * of the FiniteElement
				      * class. Called from the
				      * constructor.
				      */
    void initialize_quad_dof_index_permutation ();

				     /**
				      * Mapping from hierarchic to
				      * lexicographic numbering on
				      * first face. Hierarchic is the
				      * numbering of the shape
				      * functions.
				      */
    const std::vector<unsigned int> face_index_map;


				     /**
				      * Forward declaration of a class
				      * into which we put significant
				      * parts of the implementation.
				      *
				      * See the .cc file for more
				      * information.
				      */
    struct Implementation;

				     /**
				      * Allow access from other
				      * dimensions. We need this since
				      * we want to call the functions
				      * @p get_dpo_vector and
				      * @p lexicographic_to_hierarchic_numbering
				      * for the faces of the finite
				      * element of dimension dim+1.
				      */
    template <int, int> friend class FE_Q;

    friend class FE_Q<dim,spacedim>::Implementation;
};



/*@}*/

/* -------------- declaration of explicit specializations ------------- */

template <>
void FE_Q<1>::initialize_unit_face_support_points ();

template <>
std::vector<unsigned int>
FE_Q<1>::face_lexicographic_to_hierarchic_numbering (const unsigned int);


#if deal_II_dimension != 3

template <>
void FE_Q<1,2>::initialize_unit_face_support_points ();

template <>
std::vector<unsigned int>
FE_Q<1,2>::face_lexicographic_to_hierarchic_numbering (const unsigned int);

#endif


DEAL_II_NAMESPACE_CLOSE

#endif
