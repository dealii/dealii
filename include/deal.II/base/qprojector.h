// ---------------------------------------------------------------------
//
// Copyright (C) 2005 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__qprojector_h
#define __deal2__qprojector_h


#include <deal.II/base/quadrature.h>
#include <deal.II/base/geometry_info.h>

DEAL_II_NAMESPACE_OPEN
/*!@addtogroup Quadrature */
/*@{*/


/**
 *  This class is a helper class to facilitate the usage of quadrature formulae
 *  on faces or subfaces of cells. It computes the locations of quadrature
 *  points on the unit cell from a quadrature object for a manifold of
 *  one dimension less than that of the cell and the number of the face.
 *  For example, giving the Simpson rule in one dimension and using the
 *  project_to_face() function with face number 1, the returned points will
 *  be (1,0), (1,0.5) and (1,1). Note that faces have an orientation,
 *  so when projecting to face 3, you will get (0,0), (0,0.5) and (0,1),
 *  which is in clockwise sense, while for face 1 the points were in
 *  counterclockwise sense.
 *
 *  For the projection to subfaces (i.e. to the children of a face of the
 *  unit cell), the same applies as above. Note the order in which the
 *  children of a face are numbered, which in two dimensions coincides
 *  with the orientation of the face.
 *
 *  The second set of functions generates a quadrature formula by
 *  projecting a given quadrature rule on <b>all</b> faces and
 *  subfaces. This is used in the FEFaceValues and
 *  FESubfaceValues classes. Since we now have the quadrature
 *  points of all faces and subfaces in one array, we need to have a
 *  way to find the starting index of the points and weights
 *  corresponding to one face or subface within this array. This is
 *  done through the DataSetDescriptor member class.
 *
 *  The different functions are grouped into a common class to avoid
 *  putting them into global namespace. However, since they have no
 *  local data, all functions are declared <tt>static</tt> and can be
 *  called without creating an object of this class.
 *
 *  For the 3d case, you should note that the orientation of faces is
 *  even more intricate than for two dimensions. Quadrature formulae
 *  are projected upon the faces in their standard orientation, not to
 *  the inside or outside of the hexahedron. To make things more
 *  complicated, in 3d we allow faces in two orientations (which can
 *  be identified using <tt>cell->face_orientation(face)</tt>), so we
 *  have to project quadrature formula onto faces and subfaces in two
 *  orientations. (Refer to the documentation of the Triangulation
 *  class for a description of the orientation of the different faces,
 *  as well as to
 *  @ref GlossFaceOrientation "the glossary entry on face orientation"
 *  for more information on this.) The
 *  DataSetDescriptor member class is used to identify where each
 *  dataset starts.
 *
 *  @author Wolfgang Bangerth, Guido Kanschat, 1998, 1999, 2003, 2005
 */
template <int dim>
class QProjector
{
public:
  /**
   * Define a typedef for a
   * quadrature that acts on an
   * object of one dimension
   * less. For cells, this would
   * then be a face quadrature.
   */
  typedef Quadrature<dim-1> SubQuadrature;

  /**
   * Compute the quadrature points
   * on the cell if the given
   * quadrature formula is used on
   * face <tt>face_no</tt>. For further
   * details, see the general doc
   * for this class.
   */
  static void project_to_face (const SubQuadrature &quadrature,
                               const unsigned int      face_no,
                               std::vector<Point<dim> > &q_points);

  /**
   * Compute the cell quadrature
   * formula corresponding to using
   * <tt>quadrature</tt> on face
   * <tt>face_no</tt>. For further
   * details, see the general doc
   * for this class.
   */
  static Quadrature<dim>
  project_to_face (const SubQuadrature &quadrature,
                   const unsigned int      face_no);

  /**
   * Compute the quadrature points on the
   * cell if the given quadrature formula is
   * used on face <tt>face_no</tt>, subface
   * number <tt>subface_no</tt> corresponding
   * to RefineCase::Type
   * <tt>ref_case</tt>. The last argument is
   * only used in 3D.
   *
   * @note Only the points are
   * transformed. The quadrature
   * weights are the same as those
   * of the original rule.
   */
  static void project_to_subface (const SubQuadrature       &quadrature,
                                  const unsigned int         face_no,
                                  const unsigned int         subface_no,
                                  std::vector<Point<dim> >  &q_points,
                                  const RefinementCase<dim-1> &ref_case=RefinementCase<dim-1>::isotropic_refinement);

  /**
   * Compute the cell quadrature formula
   * corresponding to using
   * <tt>quadrature</tt> on subface
   * <tt>subface_no</tt> of face
   * <tt>face_no</tt> with
   * RefinementCase<dim-1>
   * <tt>ref_case</tt>. The last argument is
   * only used in 3D.
   *
   * @note Only the points are
   * transformed. The quadrature
   * weights are the same as those
   * of the original rule.
   */
  static Quadrature<dim>
  project_to_subface (const SubQuadrature       &quadrature,
                      const unsigned int         face_no,
                      const unsigned int         subface_no,
                      const RefinementCase<dim-1> &ref_case=RefinementCase<dim-1>::isotropic_refinement);

  /**
   * Take a face quadrature formula
   * and generate a cell quadrature
   * formula from it where the
   * quadrature points of the given
   * argument are projected on all
   * faces.
   *
   * The weights of the new rule
   * are replications of the
   * original weights. Thus, the
   * sum of the weights is not one,
   * but the number of faces, which
   * is the surface of the
   * reference cell.
   *
   * This in particular allows us
   * to extract a subset of points
   * corresponding to a single face
   * and use it as a quadrature on
   * this face, as is done in
   * FEFaceValues.
   *
   * @note In 3D, this function
   * produces eight sets of
   * quadrature points for each
   * face, in order to cope
   * possibly different
   * orientations of the mesh.
   */
  static Quadrature<dim>
  project_to_all_faces (const SubQuadrature &quadrature);

  /**
   * Take a face quadrature formula
   * and generate a cell quadrature
   * formula from it where the
   * quadrature points of the given
   * argument are projected on all
   * subfaces.
   *
   * Like in project_to_all_faces(),
   * the weights of the new rule
   * sum up to the number of faces
   * (not subfaces), which
   * is the surface of the
   * reference cell.
   *
   * This in particular allows us
   * to extract a subset of points
   * corresponding to a single subface
   * and use it as a quadrature on
   * this face, as is done in
   * FESubfaceValues.
   */
  static Quadrature<dim>
  project_to_all_subfaces (const SubQuadrature &quadrature);

  /**
   * Project a given quadrature
   * formula to a child of a
   * cell. You may want to use this
   * function in case you want to
   * extend an integral only over
   * the area which a potential
   * child would occupy. The child
   * numbering is the same as the
   * children would be numbered
   * upon refinement of the cell.
   *
   * As integration using this
   * quadrature formula now only
   * extends over a fraction of the
   * cell, the weights of the
   * resulting object are divided by
   * GeometryInfo<dim>::children_per_cell.
   */
  static
  Quadrature<dim>
  project_to_child (const Quadrature<dim>  &quadrature,
                    const unsigned int      child_no);

  /**
   * Project a quadrature rule to
   * all children of a
   * cell. Similarly to
   * project_to_all_subfaces(),
   * this function replicates the
   * formula generated by
   * project_to_child() for all
   * children, such that the
   * weights sum up to one, the
   * volume of the total cell
   * again.
   *
   * The child
   * numbering is the same as the
   * children would be numbered
   * upon refinement of the cell.
   */
  static
  Quadrature<dim>
  project_to_all_children (const Quadrature<dim>  &quadrature);

  /**
   * Project the onedimensional
   * rule <tt>quadrature</tt> to
   * the straight line connecting
   * the points <tt>p1</tt> and
   * <tt>p2</tt>.
   */
  static
  Quadrature<dim>
  project_to_line(const Quadrature<1> &quadrature,
                  const Point<dim> &p1,
                  const Point<dim> &p2);

  /**
   * Since the
   * project_to_all_faces() and
   * project_to_all_subfaces()
   * functions chain together the
   * quadrature points and weights
   * of all projections of a face
   * quadrature formula to the
   * faces or subfaces of a cell,
   * we need a way to identify
   * where the starting index of
   * the points and weights for a
   * particular face or subface
   * is. This class provides this:
   * there are static member
   * functions that generate
   * objects of this type, given
   * face or subface indices, and
   * you can then use the generated
   * object in place of an integer
   * that denotes the offset of a
   * given dataset.
   *
   * @author Wolfgang Bangerth, 2003
   */
  class DataSetDescriptor
  {
  public:
    /**
     * Default constructor. This
     * doesn't do much except
     * generating an invalid
     * index, since you didn't
     * give a valid descriptor of
     * the cell, face, or subface
     * you wanted.
     */
    DataSetDescriptor ();

    /**
     * Static function to
     * generate the offset of a
     * cell. Since we only have
     * one cell per quadrature
     * object, this offset is of
     * course zero, but we carry
     * this function around for
     * consistency with the other
     * static functions.
     */
    static DataSetDescriptor cell ();

    /**
     * Static function to generate an
     * offset object for a given face of a
     * cell with the given face
     * orientation, flip and rotation. This
     * function of course is only allowed
     * if <tt>dim>=2</tt>, and the face
     * orientation, flip and rotation are
     * ignored if the space dimension
     * equals 2.
     *
     * The last argument denotes
     * the number of quadrature
     * points the
     * lower-dimensional face
     * quadrature formula (the
     * one that has been
     * projected onto the faces)
     * has.
     */
    static
    DataSetDescriptor
    face (const unsigned int face_no,
          const bool         face_orientation,
          const bool         face_flip,
          const bool         face_rotation,
          const unsigned int n_quadrature_points);

    /**
     * Static function to generate an
     * offset object for a given subface of
     * a cell with the given face
     * orientation, flip and rotation. This
     * function of course is only allowed
     * if <tt>dim>=2</tt>, and the face
     * orientation, flip and rotation are
     * ignored if the space dimension
     * equals 2.
     *
     * The last but one argument denotes
     * the number of quadrature
     * points the
     * lower-dimensional face
     * quadrature formula (the
     * one that has been
     * projected onto the faces)
     * has.
     *
     * Through the last argument
     * anisotropic refinement can be
     * respected.
     */
    static
    DataSetDescriptor
    subface (const unsigned int face_no,
             const unsigned int subface_no,
             const bool         face_orientation,
             const bool         face_flip,
             const bool         face_rotation,
             const unsigned int n_quadrature_points,
             const internal::SubfaceCase<dim> ref_case=internal::SubfaceCase<dim>::case_isotropic);

    /**
     * Conversion operator to an
     * integer denoting the
     * offset of the first
     * element of this dataset in
     * the set of quadrature
     * formulas all projected
     * onto faces and
     * subfaces. This conversion
     * operator allows us to use
     * offset descriptor objects
     * in place of integer
     * offsets.
     */
    operator unsigned int () const;

  private:
    /**
     * Store the integer offset
     * for a given cell, face, or
     * subface.
     */
    const unsigned int dataset_offset;

    /**
     * This is the real
     * constructor, but it is
     * private and thus only
     * available to the static
     * member functions above.
     */
    DataSetDescriptor (const unsigned int dataset_offset);
  };

private:
  /**
   * Given a quadrature object in
   * 2d, reflect all quadrature
   * points at the main diagonal
   * and return them with their
   * original weights.
   *
   * This function is necessary for
   * projecting a 2d quadrature
   * rule onto the faces of a 3d
   * cube, since there we need both
   * orientations.
   */
  static Quadrature<2> reflect (const Quadrature<2> &q);

  /**
   * Given a quadrature object in
   * 2d, rotate all quadrature
   * points by @p n_times * 90 degrees
   * counterclockwise
   * and return them with their
   * original weights.
   *
   * This function is necessary for
   * projecting a 2d quadrature
   * rule onto the faces of a 3d
   * cube, since there we need all
   * rotations to account for
   * face_flip and face_rotation
   * of non-standard faces.
   */
  static Quadrature<2> rotate (const Quadrature<2> &q,
                               const unsigned int n_times);
};

/*@}*/


// -------------------  inline and template functions ----------------



template <int dim>
inline
QProjector<dim>::DataSetDescriptor::
DataSetDescriptor (const unsigned int dataset_offset)
  :
  dataset_offset (dataset_offset)
{}


template <int dim>
inline
QProjector<dim>::DataSetDescriptor::
DataSetDescriptor ()
  :
  dataset_offset (numbers::invalid_unsigned_int)
{}



template <int dim>
typename QProjector<dim>::DataSetDescriptor
QProjector<dim>::DataSetDescriptor::cell ()
{
  return 0;
}



template <int dim>
inline
QProjector<dim>::DataSetDescriptor::operator unsigned int () const
{
  return dataset_offset;
}


/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN


template <>
void
QProjector<1>::project_to_face (const Quadrature<0> &,
                                const unsigned int,
                                std::vector<Point<1> > &);
template <>
void
QProjector<2>::project_to_face (const Quadrature<1>      &quadrature,
                                const unsigned int        face_no,
                                std::vector<Point<2> >   &q_points);
template <>
void
QProjector<3>::project_to_face (const Quadrature<2>    &quadrature,
                                const unsigned int      face_no,
                                std::vector<Point<3> > &q_points);

template <>
Quadrature<1>
QProjector<1>::project_to_all_faces (const Quadrature<0> &quadrature);


template <>
void
QProjector<1>::project_to_subface (const Quadrature<0> &,
                                   const unsigned int,
                                   const unsigned int,
                                   std::vector<Point<1> > &,
                                   const RefinementCase<0> &);
template <>
void
QProjector<2>::project_to_subface (const Quadrature<1>    &quadrature,
                                   const unsigned int      face_no,
                                   const unsigned int      subface_no,
                                   std::vector<Point<2> > &q_points,
                                   const RefinementCase<1> &);
template <>
void
QProjector<3>::project_to_subface (const Quadrature<2>       &quadrature,
                                   const unsigned int         face_no,
                                   const unsigned int         subface_no,
                                   std::vector<Point<3> >    &q_points,
                                   const RefinementCase<2> &face_ref_case);

template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces (const Quadrature<0> &quadrature);


template <>
bool
QIterated<1>::uses_both_endpoints (const Quadrature<1> &base_quadrature);

template <>
QIterated<1>::QIterated (const Quadrature<1> &base_quadrature,
                         const unsigned int   n_copies);


#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
