// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_fe_q_base_h
#define dealii_fe_q_base_h

#include <deal.II/base/config.h>

#include <deal.II/base/mutex.h>

#include <deal.II/fe/fe_poly.h>

DEAL_II_NAMESPACE_OPEN


/**
 * @addtogroup fe
 * @{
 */

/**
 * This class collects the basic methods used in FE_Q, FE_Q_DG0 and
 * FE_Q_Bubbles. There is no public constructor for this class as it is not
 * functional as a stand-alone. The completion of definitions is left to the
 * derived classes.
 */
template <int dim, int spacedim = dim>
class FE_Q_Base : public FE_Poly<dim, spacedim>
{
public:
  /**
   * Constructor.
   */
  FE_Q_Base(const ScalarPolynomialsBase<dim> &poly_space,
            const FiniteElementData<dim>     &fe_data,
            const std::vector<bool>          &restriction_is_additive_flags);

  /**
   * Return the matrix interpolating from the given finite element to the
   * present one. The size of the matrix is then @p dofs_per_cell times
   * <tt>source.n_dofs_per_cell()</tt>.
   *
   * These matrices are only available if the source element is also a @p FE_Q
   * element. Otherwise, an exception of type
   * FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is thrown.
   */
  virtual void
  get_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                           FullMatrix<double> &matrix) const override;


  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. The
   * FE_Q element family only provides interpolation matrices for elements of
   * the same type and FE_Nothing. For all other elements, an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is
   * thrown.
   */
  virtual void
  get_face_interpolation_matrix(const FiniteElement<dim, spacedim> &source,
                                FullMatrix<double>                 &matrix,
                                const unsigned int face_no = 0) const override;

  /**
   * Return the matrix interpolating from a face of one element to the face
   * of the neighboring element.  The size of the matrix is then
   * <tt>source.dofs_per_face</tt> times <tt>this->dofs_per_face</tt>. The
   * FE_Q element family only provides interpolation matrices for elements of
   * the same type and FE_Nothing. For all other elements, an exception of
   * type FiniteElement<dim,spacedim>::ExcInterpolationNotImplemented is
   * thrown.
   */
  virtual void
  get_subface_interpolation_matrix(
    const FiniteElement<dim, spacedim> &source,
    const unsigned int                  subface,
    FullMatrix<double>                 &matrix,
    const unsigned int                  face_no = 0) const override;

  /**
   * This function returns @p true, if the shape function @p shape_index has
   * non-zero function values somewhere on the face @p face_index.
   */
  virtual bool
  has_support_on_face(const unsigned int shape_index,
                      const unsigned int face_index) const override;

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
   *
   * If projection matrices are not implemented in the derived finite element
   * class, this function aborts with ExcProjectionVoid. You can check whether
   * this is the case by calling the restriction_is_implemented() or the
   * isotropic_restriction_is_implemented() function.
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
   *
   * If projection matrices are not implemented in the derived finite element
   * class, this function aborts with ExcEmbeddingVoid. You can check whether
   * this is the case by calling the prolongation_is_implemented() or the
   * isotropic_prolongation_is_implemented() function.
   */
  virtual const FullMatrix<double> &
  get_prolongation_matrix(
    const unsigned int         child,
    const RefinementCase<dim> &refinement_case =
      RefinementCase<dim>::isotropic_refinement) const override;

  /**
   * Given an index in the natural ordering of indices on a face, return the
   * index of the same degree of freedom on the cell.
   *
   * To explain the concept, consider the case where we would like to know
   * whether a degree of freedom on a face, for example as part of an FESystem
   * element, is primitive. Unfortunately, the is_primitive() function in the
   * FiniteElement class takes a cell index, so we would need to find the cell
   * index of the shape function that corresponds to the present face index.
   * This function does that.
   *
   * Code implementing this would then look like this:
   * @code
   * for (i=0; i<dofs_per_face; ++i)
   *  if (fe.is_primitive(fe.face_to_cell_index(i, some_face_no)))
   *   ... do whatever
   * @endcode
   * The function takes additional arguments that account for the fact that
   * actual faces can be in their standard ordering with respect to the cell
   * under consideration, or can be flipped, oriented, etc.
   *
   * @param face_dof_index The index of the degree of freedom on a face. This
   * index must be between zero and dofs_per_face.
   * @param face The number of the face this degree of freedom lives on. This
   * number must be between zero and GeometryInfo::faces_per_cell.
   * @param combined_orientation The combined orientation flag containing the
   * orientation, rotation, and flip of the face. See
   * @ref GlossCombinedOrientation.
   * @return The index of this degree of freedom within the set of degrees of
   * freedom on the entire cell. The returned value will be between zero and
   * dofs_per_cell.
   */
  virtual unsigned int
  face_to_cell_index(const unsigned int                 face_dof_index,
                     const unsigned int                 face,
                     const types::geometric_orientation combined_orientation =
                       numbers::default_geometric_orientation) const override;

  /**
   * Return a list of constant modes of the element. For this element, the
   * list consists of true arguments for all components.
   */
  virtual std::pair<Table<2, bool>, std::vector<unsigned int>>
  get_constant_modes() const override;

  /**
   * @name Functions to support hp
   * @{
   */

  /**
   * Return whether this element implements its hanging node constraints in
   * the new way, which has to be used to make elements "hp-compatible".
   *
   * For the FE_Q class the result is always true (independent of the degree
   * of the element), as it implements the complete set of functions necessary
   * for hp-capability.
   */
  virtual bool
  hp_constraints_are_implemented() const override;

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
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_vertex_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on lines.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_line_dof_identities(
    const FiniteElement<dim, spacedim> &fe_other) const override;

  /**
   * Same as hp_vertex_dof_indices(), except that the function treats degrees
   * of freedom on quads.
   */
  virtual std::vector<std::pair<unsigned int, unsigned int>>
  hp_quad_dof_identities(const FiniteElement<dim, spacedim> &fe_other,
                         const unsigned int face_no = 0) const override;

  /** @} */

  /**
   * Attempt to construct an FE_Q object of degree 0
   *
   * @ingroup Exceptions
   */
  DeclExceptionMsg(ExcFEQCannotHaveDegree0,
                   "FE_Q can only be used for polynomial degrees "
                   "greater than zero. If you want an element of polynomial "
                   "degree zero, then it cannot be continuous and you "
                   "will want to use FE_DGQ<dim>(0).");

protected:
  /**
   * Only for internal use. Its full name is @p get_dofs_per_object_vector
   * function and it creates the @p dofs_per_object vector that is needed
   * within the constructor to be passed to the constructor of @p
   * FiniteElementData.
   */
  static std::vector<unsigned int>
  get_dpo_vector(const unsigned int degree);

  /**
   * Perform the initialization of the element based on 1d support points,
   * i.e., sets renumbering, initializes unit support points, initializes
   * constraints as well as restriction and prolongation matrices.
   */
  void
  initialize(const std::vector<Point<1>> &support_points_1d);

  /**
   * Initialize the hanging node constraints matrices. Called from
   * initialize().
   */
  void
  initialize_constraints(const std::vector<Point<1>> &points);

  /**
   * Initialize the @p unit_support_points field of the FiniteElement class.
   * Called from initialize().
   */
  void
  initialize_unit_support_points(const std::vector<Point<1>> &points);

  /**
   * Initialize the @p unit_face_support_points field of the FiniteElement
   * class. Called from initialize().
   */
  void
  initialize_unit_face_support_points(const std::vector<Point<1>> &points);

  /**
   * Initialize the @p adjust_quad_dof_index_for_face_orientation_table and
   * adjust_line_dof_index_for_line_orientation_table tables of the
   * FiniteElement class. Called from initialize().
   */
  void
  initialize_dof_index_permutations();

  /**
   * Forward declaration of a class into which we put significant parts of the
   * implementation.
   *
   * See the .cc file for more information.
   */
  struct Implementation;

  // Declare implementation friend.
  friend struct FE_Q_Base<dim, spacedim>::Implementation;

private:
  /**
   * Mutex variables used for protecting the initialization of restriction
   * and embedding matrices.
   */
  mutable Threads::Mutex restriction_matrix_mutex;
  mutable Threads::Mutex prolongation_matrix_mutex;

  /**
   * The highest polynomial degree of the underlying tensor product space
   * without any enrichment. For FE_Q*(p) this is p. Note that enrichments
   * may lead to a difference to degree.
   */
  const unsigned int q_degree;
};


/** @} */

DEAL_II_NAMESPACE_CLOSE

#endif
