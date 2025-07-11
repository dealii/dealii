// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2005 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_qprojector_h
#define dealii_qprojector_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/hp/q_collection.h>

DEAL_II_NAMESPACE_OPEN

#ifndef DOXYGEN
class ReferenceCell;
#endif

/**
 * @addtogroup Quadrature
 * @{
 */

/**
 * @brief Class which transforms `dim - 1`-dimensional quadrature rules to
 * `dim`-dimensional face quadratures.
 *
 * The majority of the finite element infrastructure, such as FE_Q and
 * FE_SimplexP, uses polynomials defined on a reference cell: for example,
 * FE_Q<3> defines a polynomial space whose domain is the unit hexahedron (i.e.,
 * ReferenceCells::Hexahedron). Hence, computing quadratures using shape
 * functions on a face of a reference cell requires converting a
 * lower-dimensional Quadrature into one defined on the boundary of the
 * higher-dimensional object, e.g., converting a Quadrature defined on a
 * quadrilateral into one defined on one face of a hexahedron.
 *
 * QProjector computes the locations of quadrature points on faces or subfaces
 * of reference cells from a Quadrature of one dimension less than that of the
 * cell, face (and possibly also subface) number, and orientation. For example,
 * calling QProjector::project_to_face() with QSimpson<1>, face number 1, and
 * numbers::default_geometric_orientation returns a Quadrature with points
 * (1,0), (1,0.5), and (1,1) with weights equal to the 1d case. Similarly, if we
 * instead use face 3 and numbers::reverse_line_orientation we obtain points
 * (1,1), (0.5,1) and (0,1). Projection to subfaces works in the same way.
 *
 * In practice, computing face integrals (e.g., via FEFaceValues or
 * FESubfaceValues) requires quadrature rules for all possible permutations of
 * face number, subface number, and combined orientation. This class provides
 * several functions for doing just that, such as
 * QProjector::project_to_all_faces(). Furthermore, the DataSetDescriptor class
 * implements indexing for converting between face number, subface number, and
 * combined orientation to the index of the associated Quadrature rule.
 */
template <int dim>
class QProjector
{
public:
  /**
   * Define an alias for a quadrature that acts on an object of one dimension
   * less. For cells, this would then be a face quadrature.
   */
  using SubQuadrature = Quadrature<dim - 1>;

  /**
   * Compute the quadrature points on the cell if the given quadrature formula
   * is used on face <tt>face_no</tt>. For further details, see the general
   * doc for this class.
   *
   * @deprecated Use the version of this function which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of this function which takes a combined_orientation "
    "argument instead.")
  static void
  project_to_face(const ReferenceCell     &reference_cell,
                  const SubQuadrature     &quadrature,
                  const unsigned int       face_no,
                  std::vector<Point<dim>> &q_points);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on face <tt>face_no</tt>. For further details, see
   * the general doc for this class.
   *
   * @deprecated Use the version of this function which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of this function which takes a combined_orientation "
    "argument instead.")
  static Quadrature<dim>
  project_to_face(const ReferenceCell &reference_cell,
                  const SubQuadrature &quadrature,
                  const unsigned int   face_no);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on face <tt>face_no</tt> taking into account the
   * orientation of the face. For further details, see the general doc for this
   * class.
   *
   * @deprecated Use the version of project_to_face() which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of project_to_face() which takes a combined_orientation "
    "argument instead.")
  static Quadrature<dim>
  project_to_oriented_face(const ReferenceCell &reference_cell,
                           const SubQuadrature &quadrature,
                           const unsigned int   face_no,
                           const bool           face_orientation,
                           const bool           face_flip,
                           const bool           face_rotation);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on face <tt>face_no</tt>. For further details, see
   * the general doc for this class.
   */
  static Quadrature<dim>
  project_to_face(const ReferenceCell               &reference_cell,
                  const SubQuadrature               &quadrature,
                  const unsigned int                 face_no,
                  const types::geometric_orientation combined_orientation);

  /**
   * Compute the quadrature points on the cell if the given quadrature formula
   * is used on face <tt>face_no</tt>, subface number <tt>subface_no</tt>
   * corresponding to RefineCase::Type <tt>ref_case</tt>. The last argument is
   * only used in 3d.
   *
   * @note Only the points are transformed. The quadrature weights are the
   * same as those of the original rule.
   *
   * @deprecated Use the version of project_to_subface() which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of project_to_subface() which takes a "
    "combined_orientation argument instead.")
  static void
  project_to_subface(const ReferenceCell           &reference_cell,
                     const SubQuadrature           &quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     std::vector<Point<dim>>       &q_points,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on subface <tt>subface_no</tt> of face
   * <tt>face_no</tt> with RefinementCase<dim-1> <tt>ref_case</tt>. The last
   * argument is only used in 3d.
   *
   * @note Only the points are transformed. The quadrature weights are the
   * same as those of the original rule.
   *
   * @note This function is deprecated since it makes an implicit assumption
   * that the cell is a line (1D), a quad (2d), or a hex (3d). Use the other
   * version of this function that takes the reference cell type instead.
   *
   * @deprecated Use the version of project_to_subface() which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of project_to_subface() which takes a "
    "combined_orientation argument instead.")
  static Quadrature<dim>
  project_to_subface(const ReferenceCell           &reference_cell,
                     const SubQuadrature           &quadrature,
                     const unsigned int             face_no,
                     const unsigned int             subface_no,
                     const RefinementCase<dim - 1> &ref_case =
                       RefinementCase<dim - 1>::isotropic_refinement);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on subface <tt>subface_no</tt> of face
   * <tt>face_no</tt> with SubfaceCase<dim> <tt>ref_case</tt>. The last
   * argument is only used in 3d.
   *
   * @note Only the points are transformed. The quadrature weights are the
   * same as those of the original rule.
   *
   * @deprecated Use the version of project_to_subface() which takes a
   * combined_orientation argument instead.
   */
  DEAL_II_DEPRECATED_EARLY_WITH_COMMENT(
    "Use the version of project_to_subface() which takes a "
    "combined_orientation argument instead.")
  static Quadrature<dim>
  project_to_oriented_subface(const ReferenceCell             &reference_cell,
                              const SubQuadrature             &quadrature,
                              const unsigned int               face_no,
                              const unsigned int               subface_no,
                              const bool                       face_orientation,
                              const bool                       face_flip,
                              const bool                       face_rotation,
                              const internal::SubfaceCase<dim> ref_case);

  /**
   * Compute the cell quadrature formula corresponding to using
   * <tt>quadrature</tt> on subface <tt>subface_no</tt> of face
   * <tt>face_no</tt> with RefinementCase<dim-1> <tt>ref_case</tt>. The last
   * argument is only used in 3d.
   *
   * @note Only the points are transformed. The quadrature weights are the
   * same as those of the original rule.
   */
  static Quadrature<dim>
  project_to_subface(const ReferenceCell               &reference_cell,
                     const SubQuadrature               &quadrature,
                     const unsigned int                 face_no,
                     const unsigned int                 subface_no,
                     const types::geometric_orientation combined_orientation,
                     const RefinementCase<dim - 1>     &ref_case);

  /**
   * Take a collection of face quadrature formulas and generate a cell
   * quadrature formula from it where the quadrature points of the given
   * argument are projected on all faces.
   *
   * The weights of the new rule are replications of the original weights.
   * Thus, the sum of the weights is not one, but the number of faces, which
   * is the surface of the reference cell.
   *
   * This in particular allows us to extract a subset of points corresponding
   * to a single face and use it as a quadrature on this face, as is done in
   * FEFaceValues.
   *
   * @note This function creates ReferenceCell::n_face_orientations() sets of
   * quadrature points for each face which are indexed (by combined orientation
   * and face number) by a DataSetDescriptor.
   */
  static Quadrature<dim>
  project_to_all_faces(const ReferenceCell            &reference_cell,
                       const hp::QCollection<dim - 1> &quadrature);

  /**
   * Like the above function, applying the same face quadrature
   * formula on all faces.
   */
  static Quadrature<dim>
  project_to_all_faces(const ReferenceCell       &reference_cell,
                       const Quadrature<dim - 1> &quadrature);

  /**
   * Take a face quadrature formula and generate a cell quadrature formula
   * from it where the quadrature points of the given argument are projected
   * on all subfaces.
   *
   * Like in project_to_all_faces(), the weights of the new rule sum up to the
   * number of faces (not subfaces), which is the surface of the reference
   * cell.
   *
   * This in particular allows us to extract a subset of points corresponding
   * to a single subface and use it as a quadrature on this face, as is done
   * in FESubfaceValues.
   */
  static Quadrature<dim>
  project_to_all_subfaces(const ReferenceCell &reference_cell,
                          const SubQuadrature &quadrature);

  /**
   * Project a given quadrature formula to a child of a cell. You may want to
   * use this function in case you want to extend an integral only over the
   * area which a potential child would occupy. The child numbering is the
   * same as the children would be numbered upon refinement of the cell.
   *
   * As integration using this quadrature formula now only extends over a
   * fraction of the cell, the weights of the resulting object are divided by
   * GeometryInfo<dim>::children_per_cell.
   *
   * @warning This function is only implemented for hypercube elements.
   */
  static Quadrature<dim>
  project_to_child(const ReferenceCell   &reference_cell,
                   const Quadrature<dim> &quadrature,
                   const unsigned int     child_no);

  /**
   * Project a quadrature rule to all children of a cell. Similarly to
   * project_to_all_subfaces(), this function replicates the formula generated
   * by project_to_child() for all children, such that the weights sum up to
   * one, the volume of the total cell again.
   *
   * The child numbering is the same as the children would be numbered upon
   * refinement of the cell.
   *
   * @warning This function is only implemented for hypercube elements.
   */
  static Quadrature<dim>
  project_to_all_children(const ReferenceCell   &reference_cell,
                          const Quadrature<dim> &quadrature);

  /**
   * Project the one dimensional rule <tt>quadrature</tt> to the straight line
   * connecting the points <tt>p1</tt> and <tt>p2</tt>.
   */
  static Quadrature<dim>
  project_to_line(const ReferenceCell &reference_cell,
                  const Quadrature<1> &quadrature,
                  const Point<dim>    &p1,
                  const Point<dim>    &p2);

  /**
   * @brief Class storing the offset index into a Quadrature rule created by
   * project_to_all_faces() or project_to_all_subfaces().
   *
   * The functions QProjector::project_to_all_faces() and
   * QProjector::project_to_all_subfaces() each combine all quadrature rules
   * (i.e., all possible combinations of face, subface, and combined
   * orientation) into a single Quadrature object. DataSetDescriptor implements
   * the correct indexing for extracting from that Quadrature rule the correct
   * index for those values.
   */
  class DataSetDescriptor
  {
  public:
    /**
     * Default constructor. This doesn't do much except generating an invalid
     * index, since you didn't give a valid descriptor of the cell, face, or
     * subface you wanted.
     */
    DataSetDescriptor();

    /**
     * Static function to generate the offset of a cell. Since we only have
     * one cell per quadrature object, this offset is of course zero, but we
     * carry this function around for consistency with the other static
     * functions.
     */
    static DataSetDescriptor
    cell();

    /**
     * Static function to generate an offset object for a given face of a cell
     * with the given face orientation, flip and rotation. This function of
     * course is only allowed if <tt>dim>=2</tt>, and the face orientation,
     * flip and rotation are ignored if the space dimension equals 2.
     *
     * The last argument denotes the number of quadrature points the
     * lower-dimensional face quadrature formula (the one that has been
     * projected onto the faces) has.
     *
     * @deprecated Use the version of this function which takes a
     * combined_orientation argument instead.
     */
    DEAL_II_DEPRECATED_WITH_COMMENT(
      "Use the version of this function which takes a combined_orientation "
      "argument instead.")
    static DataSetDescriptor
    face(const ReferenceCell &reference_cell,
         const unsigned int   face_no,
         const bool           face_orientation,
         const bool           face_flip,
         const bool           face_rotation,
         const unsigned int   n_quadrature_points);

    /**
     * Static function to generate an offset object for a given face of a cell
     * with the given combined face orientation.
     *
     * @p n_quadrature_points is the number of quadrature points the
     * lower-dimensional face quadrature formula (the one that has been
     * projected onto the faces) has.
     */
    static DataSetDescriptor
    face(const ReferenceCell               &reference_cell,
         const unsigned int                 face_no,
         const types::geometric_orientation combined_orientation,
         const unsigned int                 n_quadrature_points);

    /**
     * Compute an offset object for the given face number and orientation,
     * taking into account the possibility of different quadrature rules being
     * used on each face.
     *
     * @deprecated Use the version of this function which takes a
     * combined_orientation argument instead.
     */
    DEAL_II_DEPRECATED_WITH_COMMENT(
      "Use the version of this function which takes a combined_orientation "
      "argument instead.")
    static DataSetDescriptor
    face(const ReferenceCell            &reference_cell,
         const unsigned int              face_no,
         const bool                      face_orientation,
         const bool                      face_flip,
         const bool                      face_rotation,
         const hp::QCollection<dim - 1> &quadrature);

    /**
     * Compute an offset object for the given face number and orientation,
     * taking into account the possibility of different quadrature rules being
     * used on each face.
     *
     */
    static DataSetDescriptor
    face(const ReferenceCell               &reference_cell,
         const unsigned int                 face_no,
         const types::geometric_orientation combined_orientation,
         const hp::QCollection<dim - 1>    &quadrature);

    /**
     * Static function to generate an offset object for a given subface of a
     * cell with the given face orientation, flip and rotation. This function
     * of course is only allowed if <tt>dim>=2</tt>, and the face orientation,
     * flip and rotation are ignored if the space dimension equals 2.
     *
     * The last but one argument denotes the number of quadrature points the
     * lower-dimensional face quadrature formula (the one that has been
     * projected onto the faces) has.
     *
     * Through the last argument anisotropic refinement can be respected.
     *
     * @deprecated Use the version of this function which takes a
     * combined_orientation argument instead.
     */
    DEAL_II_DEPRECATED_WITH_COMMENT(
      "Use the version of this function which takes a combined_orientation "
      "argument instead.")
    static DataSetDescriptor
    subface(const ReferenceCell             &reference_cell,
            const unsigned int               face_no,
            const unsigned int               subface_no,
            const bool                       face_orientation,
            const bool                       face_flip,
            const bool                       face_rotation,
            const unsigned int               n_quadrature_points,
            const internal::SubfaceCase<dim> ref_case =
              internal::SubfaceCase<dim>::case_isotropic);

    /**
     * Static function to generate an offset object for a given subface of a
     * cell with the given combined face orientation. This function of course is
     * only allowed if <tt>dim>=2</tt>, and the orientation is ignored if the
     * space dimension equals 2.
     *
     * @p n_quadrature_points denotes the number of quadrature points the
     * lower-dimensional face quadrature formula (the one that has been
     * projected onto the faces) has.
     *
     * Through the last argument anisotropic refinement can be respected.
     */
    static DataSetDescriptor
    subface(const ReferenceCell               &reference_cell,
            const unsigned int                 face_no,
            const unsigned int                 subface_no,
            const types::geometric_orientation combined_orientation,
            const unsigned int                 n_quadrature_points,
            const internal::SubfaceCase<dim>   ref_case =
              internal::SubfaceCase<dim>::case_isotropic);

    /**
     * Conversion operator to an integer denoting the offset of the first
     * element of this dataset in the set of quadrature formulas all projected
     * onto faces and subfaces. This conversion operator allows us to use
     * offset descriptor objects in place of integer offsets.
     */
    operator unsigned int() const;

  private:
    /**
     * Store the integer offset for a given cell, face, or subface.
     */
    const unsigned int dataset_offset;

    /**
     * This is the real constructor, but it is private and thus only available
     * to the static member functions above.
     */
    DataSetDescriptor(const unsigned int dataset_offset);
  };
};

/** @} */


// -------------------  inline and template functions ----------------



template <int dim>
inline QProjector<dim>::DataSetDescriptor::DataSetDescriptor(
  const unsigned int dataset_offset)
  : dataset_offset(dataset_offset)
{}


template <int dim>
inline QProjector<dim>::DataSetDescriptor::DataSetDescriptor()
  : dataset_offset(numbers::invalid_unsigned_int)
{}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::cell()
{
  return 0;
}



template <int dim>
inline QProjector<dim>::DataSetDescriptor::operator unsigned int() const
{
  return dataset_offset;
}



template <int dim>
Quadrature<dim> inline QProjector<dim>::project_to_all_faces(
  const ReferenceCell       &reference_cell,
  const Quadrature<dim - 1> &quadrature)
{
  return project_to_all_faces(reference_cell,
                              hp::QCollection<dim - 1>(quadrature));
}


/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN


template <>
void
QProjector<1>::project_to_face(const ReferenceCell &reference_cell,
                               const Quadrature<0> &,
                               const unsigned int,
                               std::vector<Point<1>> &);
template <>
void
QProjector<2>::project_to_face(const ReferenceCell   &reference_cell,
                               const Quadrature<1>   &quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<2>> &q_points);
template <>
void
QProjector<3>::project_to_face(const ReferenceCell   &reference_cell,
                               const Quadrature<2>   &quadrature,
                               const unsigned int     face_no,
                               std::vector<Point<3>> &q_points);

template <>
void
QProjector<1>::project_to_subface(const ReferenceCell &reference_cell,
                                  const Quadrature<0> &,
                                  const unsigned int,
                                  const unsigned int,
                                  std::vector<Point<1>> &,
                                  const RefinementCase<0> &);
template <>
void
QProjector<2>::project_to_subface(const ReferenceCell   &reference_cell,
                                  const Quadrature<1>   &quadrature,
                                  const unsigned int     face_no,
                                  const unsigned int     subface_no,
                                  std::vector<Point<2>> &q_points,
                                  const RefinementCase<1> &);
template <>
void
QProjector<3>::project_to_subface(const ReferenceCell     &reference_cell,
                                  const Quadrature<2>     &quadrature,
                                  const unsigned int       face_no,
                                  const unsigned int       subface_no,
                                  std::vector<Point<3>>   &q_points,
                                  const RefinementCase<2> &face_ref_case);

template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces(const ReferenceCell &reference_cell,
                                       const Quadrature<0> &quadrature);


#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
