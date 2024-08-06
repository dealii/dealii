// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2001 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mapping_q_h
#define dealii_mapping_q_h


#include <deal.II/base/config.h>

#include <deal.II/base/derivative_form.h>
#include <deal.II/base/polynomial.h>
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/table.h>
#include <deal.II/base/vectorization.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/matrix_free/shape_info.h>

#include <array>
#include <cmath>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
template <int, int>
class MappingQCache;
#endif

/**
 * @addtogroup mapping
 * @{
 */

/**
 * This class implements the functionality for polynomial mappings $Q_p$ of
 * polynomial degree $p$ that will be used on all cells of the mesh.  In order
 * to get a genuine higher-order mapping for all cells, it is important to
 * provide information about how interior edges and faces of the mesh should
 * be curved. This is typically done by associating a Manifold with interior
 * cells and edges. A simple example of this is discussed in the "Results"
 * section of step-6; a full discussion of manifolds is provided in
 * step-53. If manifolds are only attached to the boundaries of a domain, the
 * current class with higher polynomial degrees will provide the same
 * information as a mere MappingQ1 object. If you are working on meshes that
 * describe a (curved) manifold embedded in higher space dimensions, i.e., if
 * dim!=spacedim, then every cell is at the boundary of the domain you will
 * likely already have attached a manifold object to all cells that can then
 * also be used by the mapping classes for higher order mappings.
 *
 * <h4>Behavior along curved boundaries and with different manifolds</h4>
 *
 * For a number of applications, one only knows a manifold description of a
 * surface but not the interior of the computational domain. In such a case, a
 * FlatManifold object will be assigned to the interior entities that
 * describes a usual planar coordinate system where the additional points for
 * the higher order mapping are placed exactly according to a bi-/trilinear
 * mapping. When combined with a non-flat manifold on the boundary, for
 * example a circle bulging into the interior of a square cell, the two
 * manifold descriptions are in general incompatible. For example, a
 * FlatManifold defined solely through the cell's vertices would put an
 * interior point located at some small distance epsilon away from the
 * boundary along a straight line and thus in general outside the concave part
 * of a circle. If the polynomial degree of MappingQ is sufficiently high, the
 * transformation from the reference cell to such a cell would in general
 * contain inverted regions close to the boundary.
 *
 * In order to avoid this situation, this class applies an algorithm for
 * making this transition smooth using a so-called transfinite interpolation
 * that is essentially a linear blend between the descriptions along the
 * surrounding entities. In the algorithm that computes additional points, the
 * compute_mapping_support_points() method, all the entities of the cells are
 * passed through hierarchically, starting from the lines to the quads and
 * finally hexes. Points on objects higher up in the hierarchy are obtained
 * from the manifold associated with that object, taking into account all the
 * points previously computed by the manifolds associated with the
 * lower-dimensional objects, not just the vertices. If only a line is
 * assigned a curved boundary but the adjacent quad is on a flat manifold, the
 * flat manifold on the quad will take the points on the deformed line into
 * account when interpolating the position of the additional points inside the
 * quad and thus always result in a well-defined transformation.
 *
 * The interpolation scheme used in this class makes sure that curved
 * descriptions can go over to flat descriptions within a single layer of
 * elements, maintaining the overall optimal convergence rates of the finite
 * element interpolation. However, this only helps as long as opposite faces
 * of a cell are far enough away from each other: If a curved part is indeed
 * curved to the extent that it would come close or even intersect some of the
 * other faces, as is often the case with long and sliver cells, the current
 * approach still leads to bad mesh quality. Therefore, the recommended way is
 * to spread the transition between curved boundaries and flat interior
 * domains over a larger range as the mesh is refined. This is provided by the
 * special manifold TransfiniteInterpolationManifold.
 */
template <int dim, int spacedim = dim>
class MappingQ : public Mapping<dim, spacedim>
{
public:
  /**
   * Constructor.  @p polynomial_degree denotes the polynomial degree of the
   * polynomials that are used to map cells from the reference to the real
   * cell.
   */
  MappingQ(const unsigned int polynomial_degree);

  /**
   * Copy constructor.
   */
  MappingQ(const MappingQ<dim, spacedim> &mapping);

  // for documentation, see the Mapping base class
  virtual std::unique_ptr<Mapping<dim, spacedim>>
  clone() const override;

  /**
   * Return the degree of the mapping, i.e. the value which was passed to the
   * constructor.
   */
  unsigned int
  get_degree() const;

  /**
   * Always returns @p true because the default implementation of functions in
   * this class preserves vertex locations.
   */
  virtual bool
  preserves_vertex_locations() const override;

  // for documentation, see the Mapping base class
  virtual BoundingBox<spacedim>
  get_bounding_box(const typename Triangulation<dim, spacedim>::cell_iterator
                     &cell) const override;

  virtual bool
  is_compatible_with(const ReferenceCell &reference_cell) const override;

  /**
   * @name Mapping points between reference and real cells
   * @{
   */

  // for documentation, see the Mapping base class
  virtual Point<spacedim>
  transform_unit_to_real_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<dim> &p) const override;

  // for documentation, see the Mapping base class
  virtual Point<dim>
  transform_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim> &p) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform_points_real_to_unit_cell(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<spacedim>>                     &real_points,
    const ArrayView<Point<dim>> &unit_points) const override;

  /**
   * @}
   */

  /**
   * @name Functions to transform tensors from reference to real coordinates
   * @{
   */

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<1, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<1, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<1, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<2, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<2, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const DerivativeForm<2, dim, spacedim>> &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  // for documentation, see the Mapping base class
  virtual void
  transform(const ArrayView<const Tensor<3, dim>>                   &input,
            const MappingKind                                        kind,
            const typename Mapping<dim, spacedim>::InternalDataBase &internal,
            const ArrayView<Tensor<3, spacedim>> &output) const override;

  /**
   * @}
   */

  /**
   * As opposed to the other fill_fe_values() and fill_fe_face_values()
   * functions that rely on pre-computed information of InternalDataBase, this
   * function chooses the flexible evaluation path on the cell and points
   * passed in to the current function.
   *
   * @param[in] cell The cell where to evaluate the mapping
   *
   * @param[in] unit_points The points in reference coordinates where the
   * transformation (Jacobians, positions) should be computed.
   *
   * @param[in] update_flags The kind of information that should be computed.
   *
   * @param[out] output_data A struct containing the evaluated quantities such
   * as the Jacobian resulting from application of the mapping on the given
   * cell with its underlying manifolds.
   */
  void
  fill_mapping_data_for_generic_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const ArrayView<const Point<dim>>                          &unit_points,
    const UpdateFlags                                           update_flags,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * As opposed to the fill_fe_face_values()
   * function that relies on pre-computed information of InternalDataBase, this
   * function chooses the flexible evaluation path on the cell and points
   * passed in to the current function.
   *
   * @param[in] cell The cell where to evaluate the mapping.
   *
   * @param[in] face_number The face number where to evaluate the mapping.
   *
   * @param[in] face_quadrature The quadrature points where the
   * transformation (Jacobians, positions) should be computed.
   *
   * @param[in] internal_data A reference to an object previously created
   * that may be used to store information the mapping can compute once on the
   * reference cell. See the documentation of the Mapping::InternalDataBase
   * class for an extensive description of the purpose of these objects.
   *
   * @param[out] output_data A struct containing the evaluated quantities such
   * as the Jacobian resulting from application of the mapping on the given
   * cell with its underlying manifolds.
   */
  void
  fill_mapping_data_for_face_quadrature(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_number,
    const Quadrature<dim - 1>                                  &face_quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const;

  /**
   * @name Interface with FEValues and friends
   * @{
   */

  /**
   * Storage for internal data of polynomial mappings. See
   * Mapping::InternalDataBase for an extensive description.
   *
   * For the current class, the InternalData class stores data that is
   * computed once when the object is created (in get_data()) as well as data
   * the class wants to store from between the call to fill_fe_values(),
   * fill_fe_face_values(), or fill_fe_subface_values() until possible later
   * calls from the finite element to functions such as transform(). The
   * latter class of member variables are marked as 'mutable'.
   */
  class InternalData : public Mapping<dim, spacedim>::InternalDataBase
  {
  public:
    /**
     * Constructor. The argument denotes the polynomial degree of the mapping
     * to which this object will correspond.
     */
    InternalData(const unsigned int polynomial_degree);

    // Documentation see Mapping::InternalDataBase.
    virtual void
    reinit(const UpdateFlags      update_flags,
           const Quadrature<dim> &quadrature) override;

    /**
     * Initialize the object's member variables related to cell and face data
     * based on the given arguments. In order to initialize cell data, this
     * function calls initialize().
     */
    void
    initialize_face(const UpdateFlags      update_flags,
                    const Quadrature<dim> &quadrature,
                    const unsigned int     n_original_q_points);

    /**
     * Return an estimate (in bytes) for the memory consumption of this object.
     */
    virtual std::size_t
    memory_consumption() const override;

    /**
     * Location of quadrature points of faces or subfaces in 3d with all
     * possible orientations. Can be accessed with the correct offset provided
     * via QProjector::DataSetDescriptor. Not needed/used for cells.
     */
    std::vector<Point<dim>> quadrature_points;

    /**
     * Unit tangential vectors. Used for the computation of boundary forms and
     * normal vectors.
     *
     * This array has `(dim-1) * GeometryInfo::faces_per_cell` entries. The
     * first GeometryInfo::faces_per_cell contain the vectors in the first
     * tangential direction for each face; the second set of
     * GeometryInfo::faces_per_cell entries contain the vectors in the second
     * tangential direction (only in 3d, since there we have 2 tangential
     * directions per face), etc.
     *
     * Filled once.
     */
    std::array<std::vector<Tensor<1, dim>>,
               GeometryInfo<dim>::faces_per_cell *(dim - 1)>
      unit_tangentials;

    /**
     * The polynomial degree of the mapping. Since the objects here are also
     * used (with minor adjustments) by MappingQ, we need to store this.
     */
    const unsigned int polynomial_degree;

    /**
     * Number of shape functions. If this is a Q1 mapping, then it is simply
     * the number of vertices per cell. However, since also derived classes
     * use this class (e.g. the Mapping_Q() class), the number of shape
     * functions may also be different.
     *
     * In general, it is $(p+1)^\text{dim}$, where $p$ is the polynomial
     * degree of the mapping.
     */
    const unsigned int n_shape_functions;

    /*
     * The default line support points. Is used in when the shape function
     * values are computed.
     *
     * The number of quadrature points depends on the degree of this
     * class, and it matches the number of degrees of freedom of an
     * FE_Q<1>(this->degree).
     */
    QGaussLobatto<1> line_support_points;

    /**
     * For the fast tensor-product path of the MappingQ class, we choose SIMD
     * vectors that are as wide as possible to minimize the number of
     * arithmetic operations. However, we do not want to choose it wider than
     * necessary, e.g., we avoid something like 8-wide AVX-512 when we only
     * compute 3 components of a 3d computation. This is because the
     * additional lanes would not do useful work, but a few operations on very
     * wide vectors can already lead to a lower clock frequency of processors
     * over long time spans (thousands of clock cycles). Hence, we choose
     * 2-wide SIMD for 1D and 2d and 4-wide SIMD for 3d. Note that we do not
     * immediately fall back to no SIMD for 1d because all architectures that
     * support SIMD also support 128-bit vectors (and none is reported to
     * reduce clock frequency for 128-bit SIMD).
     */
    using VectorizedArrayType =
      VectorizedArray<double,
                      std::min<std::size_t>(VectorizedArray<double>::size(),
                                            (dim <= 2 ? 2 : 4))>;

    /**
     * In case the quadrature rule given represents a tensor product
     * we need to store the evaluations of the 1d polynomials at
     * the 1d quadrature points. That is what this variable is for.
     */
    internal::MatrixFreeFunctions::ShapeInfo<double> shape_info;

    /**
     * In case the quadrature rule given represents a tensor product
     * we need to store temporary data in this object.
     */
    mutable AlignedVector<VectorizedArrayType> scratch;

    /**
     * Indicates whether the given Quadrature object is a tensor product.
     */
    bool tensor_product_quadrature;

    /**
     * Auxiliary vectors for internal use.
     */
    mutable std::vector<AlignedVector<Tensor<1, spacedim>>> aux;

    /**
     * Stores the support points of the mapping shape functions on the @p
     * cell_of_current_support_points.
     */
    mutable std::vector<Point<spacedim>> mapping_support_points;

    /**
     * Stores the cell of which the @p mapping_support_points are stored.
     */
    mutable typename Triangulation<dim, spacedim>::cell_iterator
      cell_of_current_support_points;

    /**
     * The determinant of the Jacobian in each quadrature point. Filled if
     * #update_volume_elements.
     */
    mutable AlignedVector<double> volume_elements;

    /**
     * Pointer to the mapping output data that holds most of the arrays,
     * including the Jacobians representing the covariant and contravariant
     * transformations.
     */
    mutable internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      *output_data;
  };

protected:
  // documentation can be found in Mapping::requires_update_flags()
  virtual UpdateFlags
  requires_update_flags(const UpdateFlags update_flags) const override;

  // documentation can be found in Mapping::get_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_data(const UpdateFlags, const Quadrature<dim> &quadrature) const override;

  using Mapping<dim, spacedim>::get_face_data;

  // documentation can be found in Mapping::get_face_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_face_data(const UpdateFlags               flags,
                const hp::QCollection<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::get_subface_data()
  virtual std::unique_ptr<typename Mapping<dim, spacedim>::InternalDataBase>
  get_subface_data(const UpdateFlags          flags,
                   const Quadrature<dim - 1> &quadrature) const override;

  // documentation can be found in Mapping::fill_fe_values()
  virtual CellSimilarity::Similarity
  fill_fe_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const CellSimilarity::Similarity                            cell_similarity,
    const Quadrature<dim>                                      &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  using Mapping<dim, spacedim>::fill_fe_face_values;

  // documentation can be found in Mapping::fill_fe_face_values()
  virtual void
  fill_fe_face_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const hp::QCollection<dim - 1>                             &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_subface_values()
  virtual void
  fill_fe_subface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const unsigned int                                          face_no,
    const unsigned int                                          subface_no,
    const Quadrature<dim - 1>                                  &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  // documentation can be found in Mapping::fill_fe_immersed_surface_values()
  virtual void
  fill_fe_immersed_surface_values(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const NonMatching::ImmersedSurfaceQuadrature<dim>          &quadrature,
    const typename Mapping<dim, spacedim>::InternalDataBase    &internal_data,
    internal::FEValuesImplementation::MappingRelatedData<dim, spacedim>
      &output_data) const override;

  /**
   * @}
   */

  /**
   * The degree of the polynomials used as shape functions for the mapping of
   * cells.
   */
  const unsigned int polynomial_degree;

  /*
   * The default line support points. These are used when computing the
   * location in real space of the support points on lines and quads, which
   * are needed by the Manifold<dim,spacedim> class.
   *
   * The number of points depends on the degree of this class, and it matches
   * the number of degrees of freedom of an FE_Q<1>(this->degree).
   */
  const std::vector<Point<1>> line_support_points;

  /*
   * The one-dimensional polynomials defined as Lagrange polynomials from the
   * line support points. These are used for point evaluations and match the
   * polynomial space of an FE_Q<1>(this->degree).
   */
  const std::vector<Polynomials::Polynomial<double>> polynomials_1d;

  /*
   * The numbering from the lexicographic to the hierarchical ordering used
   * when expanding the tensor product with the mapping support points (which
   * come in hierarchical numbers).
   */
  const std::vector<unsigned int> renumber_lexicographic_to_hierarchic;

  /*
   * The support points in reference coordinates. These are used for
   * constructing approximations of the output of
   * compute_mapping_support_points() when evaluating the mapping on the fly,
   * rather than going through the FEValues interface provided by
   * InternalData.
   *
   * The number of points depends on the degree of this class, and it matches
   * the number of degrees of freedom of an FE_Q<dim>(this->degree).
   */
  const std::vector<Point<dim>> unit_cell_support_points;

  /**
   * A vector of tables of weights by which we multiply the locations of the
   * support points on the perimeter of an object (line, quad, hex) to get the
   * location of interior support points.
   *
   * Access into this table is by @p [structdim-1], i.e., use 0 to access the
   * support point weights on a line (i.e., the interior points of the
   * GaussLobatto quadrature), use 1 to access the support point weights from
   * to perimeter to the interior of a quad, and use 2 to access the support
   * point weights from the perimeter to the interior of a hex.
   *
   * The table itself contains as many columns as there are surrounding points
   * to a particular object (2 for a line, <code>4 + 4*(degree-1)</code> for
   * a quad, <code>8 + 12*(degree-1) + 6*(degree-1)*(degree-1)</code> for a
   * hex) and as many rows as there are strictly interior points.
   *
   * For the definition of this table see equation (8) of the `mapping'
   * report.
   */
  const std::vector<Table<2, double>>
    support_point_weights_perimeter_to_interior;

  /**
   * A table of weights by which we multiply the locations of the vertex
   * points of the cell to get the location of all additional support points,
   * both on lines, quads, and hexes (as appropriate). This data structure is
   * used when we fill all support points at once, which is the case if the
   * same manifold is attached to all sub-entities of a cell. This way, we can
   * avoid some of the overhead in transforming data for mappings.
   *
   * The table has as many rows as there are vertices to the cell (2 in 1d, 4
   * in 2d, 8 in 3d), and as many rows as there are additional support points
   * in the mapping, i.e., <code>(degree+1)^dim - 2^dim</code>.
   */
  const Table<2, double> support_point_weights_cell;

  /**
   * Return the locations of support points for the mapping. For example, for
   * $Q_1$ mappings these are the vertices, and for higher order polynomial
   * mappings they are the vertices plus interior points on edges, faces, and
   * the cell interior that are placed in consultation with the Manifold
   * description of the domain and its boundary. However, other classes may
   * override this function differently. In particular, the MappingQ1Eulerian
   * class does exactly this by not computing the support points from the
   * geometry of the current cell but instead evaluating an externally given
   * displacement field in addition to the geometry of the cell.
   *
   * The default implementation of this function is appropriate for most
   * cases. It takes the locations of support points on the boundary of the
   * cell from the underlying manifold. Interior support points (ie. support
   * points in quads for 2d, in hexes for 3d) are then computed using an
   * interpolation from the lower-dimensional entities (lines, quads) in order
   * to make the transformation as smooth as possible without introducing
   * additional boundary layers within the cells due to the placement of
   * support points.
   *
   * The function works its way from the vertices (which it takes from the
   * given cell) via the support points on the line (for which it calls the
   * add_line_support_points() function) and the support points on the quad
   * faces (in 3d, for which it calls the add_quad_support_points() function).
   * It then adds interior support points that are either computed by
   * interpolation from the surrounding points using weights for transfinite
   * interpolation, or if dim<spacedim, it asks the underlying manifold for
   * the locations of interior points.
   */
  virtual std::vector<Point<spacedim>>
  compute_mapping_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell) const;

  /**
   * Transform the point @p p on the real cell to the corresponding point on
   * the unit cell @p cell by a Newton iteration.
   */
  Point<dim>
  transform_real_to_unit_cell_internal(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    const Point<spacedim>                                      &p,
    const Point<dim> &initial_p_unit) const;

  /**
   * Append the support points of all shape functions located on bounding
   * lines of the given cell to the vector @p a. Points located on the
   * vertices of a line are not included.
   *
   * This function uses the underlying manifold object of the line (or, if
   * none is set, of the cell) for the location of the requested points. This
   * function is usually called by compute_mapping_support_points() function.
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual void
  add_line_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>>                               &a) const;

  /**
   * Append the support points of all shape functions located on bounding
   * faces (quads in 3d) of the given cell to the vector @p a. This function
   * is only defined for <tt>dim=3</tt>. Points located on the vertices or
   * lines of a quad are not included.
   *
   * This function uses the underlying manifold object of the quad (or, if
   * none is set, of the cell) for the location of the requested points. This
   * function is usually called by compute_mapping_support_points().
   *
   * This function is made virtual in order to allow derived classes to choose
   * shape function support points differently than the present class, which
   * chooses the points as interpolation points on the boundary.
   */
  virtual void
  add_quad_support_points(
    const typename Triangulation<dim, spacedim>::cell_iterator &cell,
    std::vector<Point<spacedim>>                               &a) const;

  // Make MappingQCache a friend since it needs to call the
  // compute_mapping_support_points() function.
  template <int, int>
  friend class MappingQCache;
};



/**
 * A class that implements a polynomial mapping $Q_p$ of degree $p$ on all
 * cells. This class is completely equivalent to the MappingQ class and there
 * for backward compatibility.
 */
template <int dim, int spacedim = dim>
using MappingQGeneric = MappingQ<dim, spacedim>;

/** @} */


/*----------------------------------------------------------------------*/

#ifndef DOXYGEN

template <int dim, int spacedim>
inline bool
MappingQ<dim, spacedim>::preserves_vertex_locations() const
{
  return true;
}

#endif // DOXYGEN

/* -------------- declaration of explicit specializations ------------- */


DEAL_II_NAMESPACE_CLOSE

#endif
