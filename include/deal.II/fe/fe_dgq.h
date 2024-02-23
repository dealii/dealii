// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_dgq_h
#define dealii_fe_dgq_h

#include <deal.II/base/config.h>

#include <deal.II/base/mutex.h>
#include <deal.II/base/tensor_product_polynomials.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
class MappingQ;
template <int dim>
class Quadrature;
#endif

/**
 * @addtogroup fe
 * @{
 */

/**
 * Implementation of scalar, discontinuous tensor product elements based on
 * equidistant support points.
 *
 * This is a discontinuous finite element based on tensor products of
 * Lagrangian polynomials. The shape functions are Lagrangian interpolants of
 * an equidistant grid of points on the unit cell. The points are numbered in
 * lexicographical order, with <i>x</i> running fastest, then <i>y</i>, then
 * <i>z</i> (if these coordinates are present for a given space dimension at
 * all). For example, these are the node orderings for <tt>FE_DGQ(1)</tt> in
 * 3d:
 *  @verbatim
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
 *  @endverbatim
 * and <tt>FE_DGQ(2)</tt>:
 *  @verbatim
 *         24--25--26       24--25--26
 *        /|       |       /       /|
 *      21 |       |     21  22  23 |
 *      /  15  16  17    /       /  17
 *    18   |       |   18--19--20   |
 *     |12 6---7---8    |       |14 8
 *     9  /       /     9  10  11  /
 *     | 3   4   5      |       | 5
 *     |/       /       |       |/
 *     0---1---2        0---1---2
 *  @endverbatim
 * with node 13 being placed in the interior of the hex.
 *
 * Note, however, that these are just the Lagrange interpolation points of the
 * shape functions. Even though they may physically be on the boundary of the
 * cell, they are logically in the interior since there are no continuity
 * requirements for these shape functions across cell boundaries. While
 * discontinuous, when restricted to a single cell the shape functions of this
 * element are exactly the same as those of the FE_Q element where they are
 * shown visually.
 *
 * <h3>Unit support point distribution and conditioning of interpolation</h3>
 *
 * When constructing an FE_DGQ element at polynomial degrees one or two,
 * equidistant support points at 0 and 1 (linear case) or 0, 0.5, and 1
 * (quadratic case) are used. The unit support or nodal points
 * <i>x<sub>i</sub></i> are those points where the <i>j</i>th Lagrange
 * polynomial satisfies the $\delta_{ij}$ property, i.e., where one polynomial
 * is one and all the others are zero.  For higher polynomial degrees, the
 * support points are non-equidistant by default, and chosen to be the support
 * points of the <tt>(degree+1)</tt>-order Gauss-Lobatto quadrature rule. This
 * point distribution yields well-conditioned Lagrange interpolation at
 * arbitrary polynomial degrees. By contrast, polynomials based on equidistant
 * points get increasingly ill-conditioned as the polynomial degree
 * increases. In interpolation, this effect is known as the Runge
 * phenomenon. For Galerkin methods, the Runge phenomenon is typically not
 * visible in the solution quality but rather in the condition number of the
 * associated system matrices. For example, the elemental @ref GlossMassMatrix "mass matrix" of
 * equidistant points at degree 10 has condition number 2.6e6, whereas the
 * condition number for Gauss-Lobatto points is around 400.
 *
 * The Gauss-Lobatto points in 1d include the end points 0 and +1 of the unit
 * interval. The interior points are shifted towards the end points, which
 * gives a denser point distribution close to the element boundary.
 */
template <int dim, int spacedim = dim>
class FE_DGQ : public FE_Poly<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials of degree <tt>p</tt>. The
   * shape functions created using this constructor correspond to Lagrange
   * interpolation polynomials for Gauss-Lobatto support (node) points in each
   * coordinate direction.
   */
  FE_DGQ(const unsigned int p);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGQ<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. The size of the matrix is then @p dofs_per_cell times
   * <tt>source.n_dofs_per_cell()</tt>.
   *
   * These matrices are only available if the source element is also a @p
   * FE_DGQ element. Otherwise, an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double>                 &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element. The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>.
   *
   * Derived elements will have to implement this function. They may only
   * provide interpolation matrices for certain source finite elements, for
   * example those from the same family. If they don't implement interpolation
   * from a given element, then they must throw an exception of type
   * FiniteElement<dim>::ExcInterpolationNotImplemented.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double>                 &matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * Projection from a fine grid space onto a coarse grid space. Overrides the
   * respective method in FiniteElement, implementing lazy evaluation
   * (initialize when requested).
   *
   * If this projection operator is associated with a matrix @p P, then the
   * restriction of this matrix @p P_i to a single child cell is returned
   * here.
   *
   * The matrix @p P is the concatenation or the sum of the cell matrices @p
   * P_i, depending on the #restriction_is_additive_flags. This distinguishes
   * interpolation (concatenation) and projection with respect to scalar
   * products (summation).
   *
   * Row and column indices are related to coarse grid and fine grid spaces,
   * respectively, consistent with the definition of the associated operator.
   */
  virtual const FullMatrix<double> &
  get_restriction_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Embedding matrix between grids. Overrides the respective method in
   * FiniteElement, implementing lazy evaluation (initialize when queried).
   *
   * The identity operator from a coarse grid space into a fine grid space is
   * associated with a matrix @p P. The restriction of this matrix @p P_i to a
   * single child cell is returned here.
   *
   * The matrix @p P is the concatenation, not the sum of the cell matrices @p
   * P_i. That is, if the same non-zero entry <tt>j,k</tt> exists in two
   * different child matrices @p P_i, the value should be the same in both
   * matrices and it is copied into the matrix @p P only once.
   *
   * Row and column indices are related to fine grid and coarse grid spaces,
   * respectively, consistent with the definition of the associated operator.
   *
   * These matrices are used by routines assembling the prolongation matrix
   * for multi-level methods.  Upon assembling the transfer matrix between
   * cells using this matrix array, zero elements in the prolongation matrix
   * are discarded and will not fill up the transfer matrix.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * If, on a vertex, several finite elements are active, the hp-code first
   * assigns the degrees of freedom of each of these FEs different global
   * indices. It then calls this function to find out which of them should get
   * identical values, and consequently can receive the same global DoF index.
   * This function therefore returns a list of identities between DoFs of the
   * present finite element object with the DoFs of @p fe_other, which is a
   * reference to a finite element object representing one of the other finite
   * elements active on this particular vertex. The function computes which of
   * the degrees of freedom of the two finite element objects are equivalent,
   * both numbered between zero and the corresponding value of
   * n_dofs_per_vertex() of the two finite elements. The first index of each
   * pair denotes one of the vertex dofs of the present element, whereas the
   * second is the corresponding index of the other finite element.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   *
   * This being a discontinuous element, the set of such constraints is of
   * course empty.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   *
   * For the FE_DGQ class the result is always true (independent of the degree
   * of the element), as it has no hanging nodes (being a discontinuous
   * element).
   */
  virtual bool
  hp_constraints_are_implemented() const override;

  /**
   * @copydoc FiniteElement::compare_for_domination()
   */
  virtual FiniteElementDomination::Domination
  compare_for_domination(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int codim = 0) const override final;

  /**
   * @}
   */

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

  /**
   * Return a list of constant modes of the element. For this element, it
   * simply returns one row with all entries set to true.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;

protected:
  /**
   * Constructor for tensor product polynomials based on an arbitrary vector
   * of polynomials. This constructor is used in derived classes to construct
   * e.g. elements with arbitrary nodes or elements based on Legendre
   * polynomials.
   *
   * The degree of these polynomials is <tt>polynomials.size()-1</tt>.
   */
  FE_DGQ(const std::vector<Polynomials::Polynomial<double>> &polynomials);

private:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Compute renumbering for rotation of degrees of freedom.
   *
   * This function rotates a tensor product numbering of degrees of
   * freedom by 90 degrees.  It is used to compute the transfer
   * matrices of the children by using only the matrix for the first
   * child.
   *
   * The direction parameter determines the type of rotation. It is one
   * character of @p xXyYzZ. The character determines the axis of rotation,
   * case determines the direction. Lower case is counter-clockwise seen in
   * direction of the axis.
   *
   * Since rotation around the y-axis is not used, it is not implemented
   * either.
   */
  void
  rotate_indices(std::vector<unsigned int> &indices,
                 const char                 direction) const;

  /**
   * Mutex variables used for protecting the initialization of restriction
   * and embedding matrices.
   */
  mutable Threads::Mutex restriction_matrix_mutex;
  mutable Threads::Mutex prolongation_matrix_mutex;

  // Allow access from other dimensions.
  template <int dim1, int spacedim1>
  friend class FE_DGQ;

  // Allow @p MappingQ class to access to build_renumbering function.
  template <int dim1, int spacedim1>
  friend class MappingQ;
};



/**
 * Implementation of scalar, discontinuous tensor product elements based on
 * Lagrange polynomials with arbitrary nodes. The primary purpose of this
 * class is to provide an element for which the @ref GlossMassMatrix "mass matrix" can be made
 * diagonal by choosing basis functions that are not either zero or one at the
 * vertices of the cell, but instead are zero or one at a given set of
 * quadrature points. If this set of quadrature points is then also used in
 * integrating the mass matrix, then it will be diagonal. The number of
 * quadrature points automatically determines the polynomial degree chosen for
 * this element. The typical applications are the Gauss quadrature or the
 * Gauss-Lobatto quadrature (provided through the base class).
 *
 * See the base class documentation in FE_DGQ for details.
 */
template <int dim, int spacedim = dim>
class FE_DGQArbitraryNodes : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials based on Polynomials::Lagrange
   * interpolation of the support points in the quadrature rule
   * <tt>points</tt>. The degree of these polynomials is
   * <tt>points.size()-1</tt>.
   */
  FE_DGQArbitraryNodes(const Quadrature<1> &points);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGQArbitraryNodes<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by appropriate values.
   */
  virtual std::string
  get_name() const override;

  /**
   * Implementation of the corresponding function in the FiniteElement
   * class.  Since the current element is interpolatory, the nodal
   * values are exactly the support point values. Furthermore, since
   * the current element is scalar, the support point values need to
   * be vectors of length 1.
   */
  virtual void
  convert_generalized_support_point_values_to_dof_values(
    const std::vector<Vector<double>> &support_point_values,
    std::vector<double>               &nodal_values) const override;
  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};



/**
 * Implementation of scalar, discontinuous tensor product elements based on
 * Legendre polynomials, described by the tensor product of the polynomial
 * space Polynomials::Legendre. The tensor product is achieved using
 * TensorProductPolynomials and the ordering of shape functions, like in
 * TensorProductPolynomials, is lexicographic. For instance, the ordering in 2d
 * is $P_0(x)P_0(y),\ P_1(x)P_0(y),\ \ldots,\ P_n(x)P_0(y),\ P_0(x)P_1(y),
 * \ \ldots,\ P_n(x)P_1(y),\ \ldots,\ P_0(x)P_n(y),\ \ldots,\ P_n(x)P_n(y)$
 * when <tt>degree=n</tt> where $\{P_i\}_{i=0}^{n}$ are the one-dimensional
 * Legendre polynomials defined on $[0,1]$. As opposed to the basic FE_DGQ
 * element, these elements are not interpolatory and no support points are
 * defined.
 *
 * See the base class documentation in FE_DGQ for details.
 */
template <int dim, int spacedim = dim>
class FE_DGQLegendre : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials based on Polynomials::Legendre
   * interpolation.
   */
  FE_DGQLegendre(const unsigned int degree);

  /**
   * Return a list of constant modes of the element. For the Legendre basis,
   * it returns one row where the first element (corresponding to the constant
   * mode) is set to true and all other elements are set to false.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGQLegendre<dim>(degree)</tt> with <tt>dim</tt> and
   * <tt>degree</tt> replaced by the values given by the template parameter
   * and the argument passed to the constructor, respectively.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};



/**
 * Implementation of scalar, discontinuous tensor product elements based on
 * Hermite-like polynomials, described by the polynomial space
 * Polynomials::HermiteLikeInterpolation. As opposed to the basic FE_DGQ
 * element, these elements are not interpolatory and no support points are
 * defined.
 *
 * Note that Hermite polynomials are only available for degrees larger or
 * equal to three, and thus the beneficial properties of
 * Polynomials::HermiteLikeInterpolation with only two basis functions having
 * a non-trivial value or derivative on a face per dimension is only present
 * for higher degrees. To facilitate usage also for degrees zero to two, a
 * usual Lagrange basis is constructed by this class.
 *
 * See the base class documentation in FE_DGQ for details.
 */
template <int dim, int spacedim = dim>
class FE_DGQHermite : public FE_DGQ<dim, spacedim>
{
public:
  /**
   * Constructor for tensor product polynomials based on
   * Polynomials::HermiteLikeInterpolation.
   */
  FE_DGQHermite(const unsigned int degree);

  /**
   * Return a string that uniquely identifies a finite element. This class
   * returns <tt>FE_DGQHermite<dim>(degree)</tt>, with <tt>dim</tt> and
   * <tt>degree</tt> replaced by the values given by the template parameter
   * and the argument passed to the constructor, respectively.
   */
  virtual std::string
  get_name() const override;

  virtual std::unique_ptr<FiniteElement<dim, spacedim>>
  clone() const override;
};


/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
