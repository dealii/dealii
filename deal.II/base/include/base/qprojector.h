//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__qprojector_h
#define __deal2__qprojector_h


#include <base/quadrature.h>

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
 *  the inside or outside of the hexahedron. Refer to the
 *  documentation of the <tt>Triangulation</tt> class for a description of
 *  the orientation of the different faces. To make things more
 *  complicated, in 3d we allow faces in two orientations (which can
 *  be identified using <tt>cell->face_orientation(face)</tt>), so we have
 *  to project quadrature formula onto faces and subfaces in two
 *  orientations. The DataSetDescriptor member class is used to
 *  identify where each dataset starts.
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
				      * Compute the quadrature points
				      * on the cell if the given
				      * quadrature formula is used on
				      * face <tt>face_no</tt>, subface
				      * number <tt>subface_no</tt>. For
				      * further details, see the
				      * general doc for this class.
				      */
    static void project_to_subface (const SubQuadrature &quadrature,
				    const unsigned int   face_no,
				    const unsigned int   subface_no,
				    std::vector<Point<dim> > &q_points);

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
				      * original weights. This is not
				      * a proper handling, in that the
				      * sum of weights does not equal
				      * one, but it is consistent with
				      * the use of this function,
				      * namely to generate sets of
				      * face quadrature points on a
				      * cell, one set of which will
				      * then be selected at each
				      * time. This is used in the
				      * FEFaceValues class,
				      * where we initialize the
				      * values, derivatives, etc on
				      * all faces at once, while
				      * selecting the data of one
				      * particular face only happens
				      * later.
				      *
				      * @note In 3D, this function
				      * produces two sets of
				      * quadrature points for each
				      * face, in order to cope
				      * possibly different
				      * orientations of the mesh. 
				      */
    static Quadrature<dim>
    project_to_all_faces (const SubQuadrature &quadrature);

				     /**
				      * This function is alike the
				      * previous one, but projects the
				      * given face quadrature formula
				      * to the subfaces of a cell,
				      * i.e. to the children of the
				      * faces of the unit cell.
				      */
    static Quadrature<dim>
    project_to_all_subfaces (const SubQuadrature &quadrature);

				     /**
				      * Project a give quadrature
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
				      * GeometryInfo@<dim@>::children_per_cell.
				      */
    static
    Quadrature<dim>
    project_to_child (const Quadrature<dim>  &quadrature,
		      const unsigned int      child_no);

				     /**
				      * Project the onedimensional
				      * rule <tt>quadrature</tt> to
				      * the straight line connecting
				      * the points <tt>p1</tt> and
				      * <tt>p2</tt>.
				      */
    static
    Quadrature<dim>
    project_to_line(const Quadrature<1>& quadrature,
		    const Point<dim>& p1,
		    const Point<dim>& p2);
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
                                          * Static function to
                                          * generate an offset object
                                          * for a given face of a cell
                                          * with the given face
                                          * orientation. This function
                                          * of course is only allowed
                                          * if <tt>dim>=2</tt>, and the
                                          * face orientation is
                                          * ignored if the space
                                          * dimension equals 2.
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
              const unsigned int n_quadrature_points);

                                         /**
                                          * Static function to
                                          * generate an offset object
                                          * for a given subface of a
                                          * cell with the given face
                                          * orientation. This function
                                          * of course is only allowed
                                          * if <tt>dim>=2</tt>, and the
                                          * face orientation is
                                          * ignored if the space
                                          * dimension equals 2.
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
        sub_face (const unsigned int face_no,
                  const unsigned int subface_no,
                  const bool         face_orientation,
                  const unsigned int n_quadrature_points);

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
};

/*@}*/

/// @if NoDoc

/* -------------- declaration of explicit specializations ------------- */

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
                                   std::vector<Point<1> > &);
template <>
void
QProjector<2>::project_to_subface (const Quadrature<1>    &quadrature,
                                   const unsigned int      face_no,
                                   const unsigned int      subface_no,
                                   std::vector<Point<2> > &q_points);
template <>
void
QProjector<3>::project_to_subface (const Quadrature<2>    &quadrature,
                                   const unsigned int      face_no,
                                   const unsigned int      subface_no,
                                   std::vector<Point<3> > &q_points);

template <>
Quadrature<1>
QProjector<1>::project_to_all_subfaces (const Quadrature<0> &quadrature);


template <>
bool
QIterated<1>::uses_both_endpoints (const Quadrature<1> &base_quadrature);

template <>
QIterated<1>::QIterated (const Quadrature<1> &base_quadrature,
			 const unsigned int   n_copies);


/// @endif

#endif
