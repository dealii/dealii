// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2012 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#ifndef dealii_matrix_free_shape_info_h
#define dealii_matrix_free_shape_info_h


#include <deal.II/base/config.h>

#include <deal.II/base/aligned_vector.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/quadrature.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vectorization.h>


DEAL_II_NAMESPACE_OPEN


// forward declaration
template <int dim, int spacedim>
class FiniteElement;


namespace internal
{
  namespace MatrixFreeFunctions
  {
    /**
     * An enum that encodes the type of element detected during
     * initialization. FEEvaluation will select the most efficient algorithm
     * based on the given element type.
     *
     * There is an implied ordering in the type ElementType::tensor_symmetric
     * in the sense that both ElementType::tensor_symmetric_collocation and
     * ElementType::tensor_symmetric_hermite are also of type
     * ElementType::tensor_symmetric. Likewise, a configuration of type
     * ElementType::tensor_symmetric is also of type
     * ElementType::tensor_general. As a consequence, we support `<=`
     * operations between the types with this sorting, but not against the
     * even higher indexed types such as ElementType::truncated_tensor.
     *
     * @ingroup matrixfree
     */
    enum ElementType
    {
      /**
       * Tensor product shape function where the shape value evaluation in the
       * quadrature point corresponds to the identity operation and no
       * interpolation needs to be performed (collocation approach, also
       * called spectral evaluation). This is for example the case for an
       * element with nodes in the Gauss-Lobatto support points and
       * integration in the Gauss-Lobatto quadrature points of the same order.
       */
      tensor_symmetric_collocation = 0,

      /**
       * Symmetric tensor product shape functions fulfilling a Hermite
       * identity with values and first derivatives zero at the element end
       * points in 1d.
       */
      tensor_symmetric_hermite = 1,

      /**
       * Usual tensor product shape functions whose shape values and
       * quadrature points are symmetric about the midpoint of the unit
       * interval 0.5
       */
      tensor_symmetric = 2,

      /**
       * For function spaces that are not equivalent to a polynom of degree p,
       * which the assumption of collocation is.
       */
      tensor_symmetric_no_collocation = 3,

      /**
       * Tensor product shape functions without further particular properties
       */
      tensor_general = 4,

      /**
       * Polynomials of complete degree rather than tensor degree which can be
       * described by a truncated tensor product
       */
      truncated_tensor = 5,

      /**
       * Tensor product shape functions that are symmetric about the midpoint
       * of the unit interval 0.5 that additionally add a constant shape
       * function according to FE_Q_DG0.
       */
      tensor_symmetric_plus_dg0 = 6,

      /**
       * Special case of the FE_RaviartThomasNodal element with anisotropic
       * tensor product shape functions, i.e degree (k + 1) in normal direction,
       * and k in tangential direction.
       */
      tensor_raviart_thomas = 7,

      /**
       * Shape functions without a tensor product properties.
       */
      tensor_none = 8


    };



    /**
     * This struct stores the shape functions, their gradients and Hessians
     * evaluated for a one-dimensional section of a tensor product finite
     * element and tensor product quadrature formula in reference
     * coordinates. This data structure also includes the evaluation of
     * quantities at the cell boundary and on the sub-interval $(0, 0.5)$ and
     * $(0.5, 1)$ for face integrals.
     */
    template <typename Number>
    struct UnivariateShapeData
    {
      /**
       * Empty constructor. Sets default configuration.
       */
      UnivariateShapeData();

      /**
       * Return the memory consumption of this class in bytes.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Evaluate the finite element shape functions at the points of the
       * given quadrature formula, filling the fields
       * shape_[values,gradients,hessians] and related information.
       *
       * The two last arguments 'lexicographic' and 'direction' are used to
       * describe the unknowns along a single dimension, and the respective
       * direction of derivatives.
       */
      template <int dim, int spacedim>
      void
      evaluate_shape_functions(const FiniteElement<dim, spacedim> &fe,
                               const Quadrature<1>                &quad,
                               const std::vector<unsigned int> &lexicographic,
                               const unsigned int               direction);

      /**
       * Evaluate the auxiliary polynomial space associated with the Lagrange
       * polynomials in points of the given quadrature formula, filling the
       * fields shape_[gradients,hessians]_collocation and related
       * information.
       */
      template <int dim, int spacedim>
      void
      evaluate_collocation_space(const FiniteElement<dim, spacedim> &fe,
                                 const Quadrature<1>                &quad,
                                 const std::vector<unsigned int> &lexicographic,
                                 const unsigned int               direction);

      /**
       * Check whether we have symmetries in the shape values. In that case,
       * also fill the shape_???_eo fields.
       */
      bool
      check_and_set_shapes_symmetric();

      /**
       * Check whether symmetric 1d basis functions are such that the shape
       * values form a diagonal matrix, i.e., the nodal points are collocated
       * with the quadrature points. This allows for specialized algorithms
       * that save some operations in the evaluation.
       */
      bool
      check_shapes_collocation() const;

      /**
       * Encodes the type of element detected at construction. FEEvaluation
       * will select the most efficient algorithm based on the given element
       * type.
       */
      ElementType element_type;

      /**
       * Stores the shape values of the 1d finite element evaluated at all 1d
       * quadrature points. The length of
       * this array is <tt>n_dofs_1d * n_q_points_1d</tt> and quadrature
       * points are the index running fastest.
       */
      AlignedVector<Number> shape_values;

      /**
       * Stores the shape gradients of the 1d finite element evaluated at all
       * 1d quadrature points. The length of
       * this array is <tt>n_dofs_1d * n_q_points_1d</tt> and quadrature
       * points are the index running fastest.
       */
      AlignedVector<Number> shape_gradients;

      /**
       * Stores the shape Hessians of the 1d finite element evaluated at all
       * 1d quadrature points. The length of
       * this array is <tt>n_dofs_1d * n_q_points_1d</tt> and quadrature
       * points are the index running fastest.
       */
      AlignedVector<Number> shape_hessians;

      /**
       * Stores the shape gradients of the shape function space associated to
       * the quadrature (collocation), given by FE_DGQ<1>(Quadrature<1>).
       */
      AlignedVector<Number> shape_gradients_collocation;

      /**
       * Stores the shape hessians of the shape function space associated to
       * the quadrature (collocation), given by FE_DGQ<1>(Quadrature<1>).
       */
      AlignedVector<Number> shape_hessians_collocation;

      /**
       * Stores the shape values in a different format, namely the so-called
       * even-odd scheme where the symmetries in shape_values are used for
       * faster evaluation.
       */
      AlignedVector<Number> shape_values_eo;

      /**
       * Stores the shape gradients in a different format, namely the
       * so-called even-odd scheme where the symmetries in shape_gradients are
       * used for faster evaluation.
       */
      AlignedVector<Number> shape_gradients_eo;

      /**
       * Stores the shape second derivatives in a different format, namely the
       * so-called even-odd scheme where the symmetries in shape_hessians are
       * used for faster evaluation.
       */
      AlignedVector<Number> shape_hessians_eo;

      /**
       * Stores the shape gradients of the shape function space associated to
       * the quadrature (collocation), given by FE_DGQ<1>(Quadrature<1>). This
       * array provides an alternative representation of the
       * shape_gradients_collocation field in the even-odd format.
       */
      AlignedVector<Number> shape_gradients_collocation_eo;

      /**
       * Stores the shape hessians of the shape function space associated to
       * the quadrature (collocation), given by FE_DGQ<1>(Quadrature<1>). This
       * array provides an alternative representation of the
       * shape_hessians_collocation field in the even-odd format.
       */
      AlignedVector<Number> shape_hessians_collocation_eo;

      /**
       * Stores the inverse transformation from the data at quadrature points
       * to the basis defined by the shape_values fields. The data at
       * quadrature points is interpreted either implicitly by its polynomial
       * interpolation, or explicitly in terms of separate polynomials such as
       * with the `_collocation` fields. The size of the array equals the
       * layout of the `shape_values` array, and it is combined with the shape
       * values array such that this matrix is the pseudo inverse of
       * shape_values. In case the number of 1d quadrature points equals the
       * size of the basis, this array is exactly the inverse of the
       * shape_values array. The length of this array is <tt>n_dofs_1d *
       * n_q_points_1d</tt> and quadrature points are the index running
       * fastest.
       */
      AlignedVector<Number> inverse_shape_values;

      /**
       * Stores the even-odd variant of the `inverse_shape_values` field.
       */
      AlignedVector<Number> inverse_shape_values_eo;

      /**
       * Collects all data of 1d shape values evaluated at the point 0 and 1
       * (the vertices) in one data structure. The sorting of data is to
       * start with the values, then gradients, then second derivatives.
       */
      std::array<AlignedVector<Number>, 2> shape_data_on_face;

      /**
       * Collects all data of 1d nodal shape values (defined by the Lagrange
       * polynomials in the points of the quadrature rule) evaluated at the
       * point 0 and 1 (the vertices) in one data structure.
       *
       * This data structure can be used to interpolate from the cell to the
       * face quadrature points. The sorting of data is to start with the
       * values, then gradients, then second derivatives.
       */
      std::array<AlignedVector<Number>, 2> quadrature_data_on_face;

      /**
       * Stores one-dimensional values of shape functions on subface. Since
       * there are two subfaces, store two variants.
       */
      std::array<AlignedVector<Number>, 2> values_within_subface;

      /**
       * Stores one-dimensional gradients of shape functions on subface. Since
       * there are two subfaces, store two variants.
       */
      std::array<AlignedVector<Number>, 2> gradients_within_subface;

      /**
       * Stores one-dimensional gradients of shape functions on subface. Since
       * there are two subfaces, store two variants.
       */
      std::array<AlignedVector<Number>, 2> hessians_within_subface;

      /**
       * A 1d subface interpolation matrices to the first and second quadrant.
       * This data structure is only set up for FE_Q for dim > 1.
       */
      std::array<AlignedVector<Number>, 2> subface_interpolation_matrices;

      /**
       * Same as above but stored in a scalar format independent of the type of
       * Number
       */
      std::array<AlignedVector<typename dealii::internal::VectorizedArrayTrait<
                   Number>::value_type>,
                 2>
        subface_interpolation_matrices_scalar;

      /**
       * We store a copy of the one-dimensional quadrature formula
       * used for initialization.
       */
      Quadrature<1> quadrature;

      /**
       * Stores the degree of the element.
       */
      unsigned int fe_degree;

      /**
       * Stores the number of quadrature points per dimension.
       */
      unsigned int n_q_points_1d;

      /**
       * Indicates whether the basis functions are nodal in 0 and 1, i.e., the
       * end points of the unit cell.
       */
      bool nodal_at_cell_boundaries;

      /**
       * Stores the shape values of the finite element evaluated at all
       * quadrature points for all faces and orientations (no tensor-product
       * structure exploited).
       */
      Table<3, Number> shape_values_face;

      /**
       * Stores the shape gradients of the finite element evaluated at all
       * quadrature points for all faces, orientations, and directions
       * (no tensor-product structure exploited).
       */
      Table<3, Number> shape_gradients_face;
    };



    /**
     * This struct stores a tensor (Kronecker) product view of the finite
     * element and quadrature formula used for evaluation. It is based on a
     * single or a collection of UnivariateShapeData object(s) that describe
     * one-dimensional ingredients, plus some additional information about how
     * these are combined and how indices are laid out in the multi-dimensional
     * case such as the hierarchical -> lexicographic ordering of FE_Q.
     *
     * @ingroup matrixfree
     */
    template <typename Number>
    struct ShapeInfo
    {
      /**
       * Encodes the type of element detected at construction. FEEvaluation
       * will select the most efficient algorithm based on the given element
       * type.
       */
      ElementType element_type;

      /**
       * Empty constructor. Does nothing.
       */
      ShapeInfo();

      /**
       * Constructor that initializes the data fields using the reinit method.
       */
      template <int dim, int spacedim, int dim_q>
      ShapeInfo(const Quadrature<dim_q>            &quad,
                const FiniteElement<dim, spacedim> &fe,
                const unsigned int                  base_element = 0);

      /**
       * Initializes the data fields. Takes a one-dimensional quadrature
       * formula and a finite element as arguments and evaluates the shape
       * functions, gradients and Hessians on the one-dimensional unit cell.
       * This function assumes that the finite element is derived from a
       * one-dimensional element by a tensor product and that the zeroth shape
       * function in zero evaluates to one.
       */
      template <int dim, int spacedim, int dim_q>
      void
      reinit(const Quadrature<dim_q>            &quad,
             const FiniteElement<dim, spacedim> &fe_dim,
             const unsigned int                  base_element = 0);

      /**
       * Return whether an element is supported by MatrixFree.
       *
       * The following scalar elements are supported:
       * - FE_Nothing, FE_DGP, and FE_Q_DG0
       * - polynomial tensor-product elements based on
       *   Polynomials::Polynomial (FE_Q, FE_DGQ, FE_DGQArbitraryNodes,
       *   FE_DGQHermite, FE_DGQLegendre) or Polynomials::PiecewisePolynomial
       *   (FE_Q_iso_Q1).
       * - elements for simplex, pyramids, and wedges (FE_SimplexP,
       *   FE_SimplexDGP, FE_PyramidP, FE_PyramidDGP, FE_WedgeP, FE_WedgeDGP)
       *
       * In the case of vectorial elements, FE_RaviartThomasNodal
       * and FESystem with base elements from the scalar elements
       * listed above are supported.
       */
      template <int dim, int spacedim>
      static bool
      is_supported(const FiniteElement<dim, spacedim> &fe);

      /**
       * Compute a table with numbers of re-orientation for all versions of
       * face flips, orientation, and rotation (relating only to 3d elements).
       */
      static Table<2, unsigned int>
      compute_orientation_table(const unsigned int n_points_per_dim);

      /**
       * Return data of univariate shape functions which defines the
       * dimension @p dimension of tensor product shape functions
       * regarding vector component @p component of the underlying
       * finite element.
       */
      const UnivariateShapeData<Number> &
      get_shape_data(const unsigned int dimension = 0,
                     const unsigned int component = 0) const;

      /**
       * Return the memory consumption of this class in bytes.
       */
      std::size_t
      memory_consumption() const;

      /**
       * Renumbering from deal.II's numbering of cell degrees of freedom to
       * lexicographic numbering used inside the FEEvaluation schemes of the
       * underlying element in the DoFHandler. For vector-valued elements, the
       * renumbering starts with a lexicographic numbering of the first
       * component, then everything of the second component, and so on.
       */
      std::vector<unsigned int> lexicographic_numbering;

      /**
       * Stores data of univariate shape functions defining the
       * underlying tensor product finite element.
       */
      std::vector<UnivariateShapeData<Number>> data;

      /**
       * Grants access to univariate shape function data of given
       * dimension and vector component. Rows identify dimensions and
       * columns identify vector components.
       */
      dealii::Table<2, UnivariateShapeData<Number> *> data_access;

      /**
       * Stores the number of space dimensions.
       */
      unsigned int n_dimensions;

      /**
       * Stores the number of vector components of the underlying
       * vector-valued finite element.
       */
      unsigned int n_components;

      /**
       * Stores the number of quadrature points in @p dim dimensions for a
       * cell.
       */
      unsigned int n_q_points;

      /**
       * Stores the number of DoFs per cell of the scalar element in @p dim
       * dimensions.
       */
      unsigned int dofs_per_component_on_cell;

      /**
       * Stores the number of quadrature points per face in @p dim dimensions.
       */
      unsigned int n_q_points_face;

      /**
       * Stores the number of quadrature points of a face in @p dim dimensions
       * for simplex, wedge and pyramid reference cells.
       */
      std::vector<unsigned int> n_q_points_faces;

      /**
       * Stores the number of DoFs per face in @p dim dimensions.
       */
      unsigned int dofs_per_component_on_face;

      /**
       * For nodal basis functions with nodes located at the boundary of the
       * unit cell, face integrals that involve only the values of the shape
       * functions (approximations of first derivatives in DG) do not need to
       * load all degrees of freedom of the cell but rather only the degrees
       * of freedom located on the face. While it would also be possible to
       * compute these indices on the fly, we choose to simplify the code and
       * store the indirect addressing in this class.
       *
       * The first table index runs through the faces of a cell, and the
       * second runs through the nodal degrees of freedom of the face, using
       * @p dofs_per_face entries.
       *
       * The indices stored in this member variable are as follows. Consider
       * for example a 2d element of degree 3 with the following degrees of
       * freedom in lexicographic numbering:
       * @code
       * 12   13   14   15
       * 8    9    10   11
       * 4    5     6    7
       * 0    1     2    3
       * @endcode
       *
       * The first row stores the indices on the face with index 0, i.e., the
       * numbers <code>0, 4, 8, 12</code>, the second row holds the indices
       * <code>3, 7, 11, 15</code> for face 1, the third row holds the indices
       * <code>0, 1, 2, 3</code> for face 2, and the last (fourth) row holds
       * the indices <code>12, 13, 14, 15</code>. Similarly, the indices are
       * stored in 3d. (Note that the y faces in 3d use indices reversed in
       * terms of the lexicographic numbers due to the orientation of the
       * coordinate system.)
       *
       * @note This object is only filled in case @p nodal_at_cell_boundaries
       * evaluates to @p true.
       */
      dealii::Table<2, unsigned int> face_to_cell_index_nodal;

      /**
       * The @p face_to_cell_index_nodal provides a shortcut for the
       * evaluation of values on the faces. For Hermite-type basis functions,
       * one can go one step further and also use shortcuts to get derivatives
       * more cheaply where only two layers of degrees of freedom contribute
       * to the derivative on the face. In the lexicographic ordering, the
       * respective indices is in the next "layer" of degrees of freedom as
       * compared to the nodal values. This array stores the indirect
       * addressing of both the values and the gradient.
       *
       * The first table index runs through the faces of a cell, and the
       * second runs through the pairs of the nodal degrees of freedom of the
       * face and the derivatives, using <code>2*dofs_per_face</code> entries.
       *
       * The indices stored in this member variable are as follows. Consider
       * for example a 2d element of degree 3 with the following degrees of
       * freedom in lexicographic numbering:
       * @code
       * 20   21   22   23   24
       * 15   16   17   18   19
       * 10   11   12   13   14
       * 5    6     7    8    9
       * 0    1     2    3    4
       * @endcode
       *
       * The first row stores the indices for values and gradients on the face
       * with index 0, i.e., the numbers <code>0, 1, 5, 6, 10, 11, 15, 16, 20,
       * 21</code>, the second row holds the indices <code>4, 3, 9, 8, 14, 13,
       * 19, 18, 24, 23</code> for face 1, the third row holds the indices
       * <code>0, 5, 1, 6, 2, 7, 3, 8, 4, 9</code> for face 2, and the last
       * (fourth) row holds the indices <code>20, 15, 21, 16, 22, 17, 23, 18,
       * 24, 19</code>. Similarly, the indices are stored in 3d. (Note that
       * the y faces in 3d use indices reversed in terms of the lexicographic
       * numbers due to the orientation of the coordinate system.)
       *
       * @note This object is only filled in case @p element_type evaluates to
       * @p tensor_symmetric_hermite.
       */
      dealii::Table<2, unsigned int> face_to_cell_index_hermite;

      /**
       * For unknowns located on faces, the basis functions are not
       * in the correct order if a face is not in the standard orientation
       * to a given element. This data structure is used to re-order the
       * basis functions to represent the correct order.
       */
      dealii::Table<2, unsigned int> face_orientations_dofs;

      /**
       * For interpretation of values at quadrature points, the order of
       * points is not correct if a face is not in the standard orientation to
       * a given element. This data structure is used to re-order the
       * quadrature points to represent the correct order.
       */
      dealii::Table<2, unsigned int> face_orientations_quad;
    };



    // ------------------------------------------ inline functions

    template <typename Number>
    inline const UnivariateShapeData<Number> &
    ShapeInfo<Number>::get_shape_data(const unsigned int dimension,
                                      const unsigned int component) const
    {
      AssertDimension(n_dimensions, data_access.size(0));
      AssertDimension(n_components, data_access.size(1));
      AssertIndexRange(dimension, n_dimensions);
      AssertIndexRange(component, n_components);
      return *(data_access(dimension, component));
    }

  } // end of namespace MatrixFreeFunctions

} // end of namespace internal

DEAL_II_NAMESPACE_CLOSE

#endif
