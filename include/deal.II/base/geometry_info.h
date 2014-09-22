// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
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

#ifndef __deal2__geometry_info_h
#define __deal2__geometry_info_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

DEAL_II_NAMESPACE_OPEN



/**
 * A class that provides possible choices for isotropic and
 * anisotropic refinement flags in the current space dimension.
 *
 * This general template is unused except in some weird template
 * constructs. Actual is made, however, of the specializations
 * <code>RefinementPossibilities@<1@></code>,
 * <code>RefinementPossibilities@<2@></code>, and
 * <code>RefinementPossibilities@<3@></code>.
 *
 * @ingroup aniso
 * @author Ralf Hartmann, 2005, Wolfgang Bangerth, 2007
 */
template <int dim>
struct RefinementPossibilities
{
  /**
   * Possible values for refinement
   * cases in the current
   * dimension.
   *
   * Note the construction of the
   * values: the lowest bit
   * describes a cut of the x-axis,
   * the second to lowest bit
   * corresponds to a cut of the
   * y-axis and the third to lowest
   * bit corresponds to a cut of
   * the z-axis. Thus, the
   * following relations hold
   * (among others):
   *
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   *
   * Only those cuts that are
   * reasonable in a given space
   * dimension are offered, of
   * course.
   *
   * In addition, the tag
   * <code>isotropic_refinement</code>
   * denotes isotropic refinement
   * in the space dimension
   * selected by the template
   * argument of this class.
   */
  enum Possibilities
  {
    no_refinement= 0,

    isotropic_refinement = static_cast<unsigned char>(-1)
  };
};



/**
 * A class that provides possible choices for isotropic and
 * anisotropic refinement flags in the current space dimension.
 *
 * This specialization is used for <code>dim=1</code>, where it offers
 * refinement in x-direction.
 *
 * @ingroup aniso
 * @author Ralf Hartmann, 2005, Wolfgang Bangerth, 2007
 */
template <>
struct RefinementPossibilities<1>
{
  /**
   * Possible values for refinement
   * cases in the current
   * dimension.
   *
   * Note the construction of the
   * values: the lowest bit
   * describes a cut of the x-axis,
   * the second to lowest bit
   * corresponds to a cut of the
   * y-axis and the third to lowest
   * bit corresponds to a cut of
   * the z-axis. Thus, the
   * following relations hold
   * (among others):
   *
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   *
   * Only those cuts that are
   * reasonable in a given space
   * dimension are offered, of
   * course.
   *
   * In addition, the tag
   * <code>isotropic_refinement</code>
   * denotes isotropic refinement
   * in the space dimension
   * selected by the template
   * argument of this class.
   */
  enum Possibilities
  {
    no_refinement= 0,
    cut_x        = 1,

    isotropic_refinement = cut_x
  };
};



/**
 * A class that provides possible choices for isotropic and
 * anisotropic refinement flags in the current space dimension.
 *
 * This specialization is used for <code>dim=2</code>, where it offers
 * refinement in x- and y-direction separately, as well as isotropic
 * refinement in both directions at the same time.
 *
 * @ingroup aniso
 * @author Ralf Hartmann, 2005, Wolfgang Bangerth, 2007
 */
template <>
struct RefinementPossibilities<2>
{
  /**
   * Possible values for refinement
   * cases in the current
   * dimension.
   *
   * Note the construction of the
   * values: the lowest bit
   * describes a cut of the x-axis,
   * the second to lowest bit
   * corresponds to a cut of the
   * y-axis and the third to lowest
   * bit corresponds to a cut of
   * the z-axis. Thus, the
   * following relations hold
   * (among others):
   *
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   *
   * Only those cuts that are
   * reasonable in a given space
   * dimension are offered, of
   * course.
   *
   * In addition, the tag
   * <code>isotropic_refinement</code>
   * denotes isotropic refinement
   * in the space dimension
   * selected by the template
   * argument of this class.
   */
  enum Possibilities
  {
    no_refinement= 0,
    cut_x        = 1,
    cut_y        = 2,
    cut_xy       = cut_x | cut_y,

    isotropic_refinement = cut_xy
  };
};



/**
 * A class that provides possible choices for isotropic and
 * anisotropic refinement flags in the current space dimension.
 *
 * This specialization is used for <code>dim=3</code>, where it offers
 * refinement in x-, y- and z-direction separately, as well as
 * combinations of these and isotropic refinement in all directions at
 * the same time.
 *
 * @ingroup aniso
 * @author Ralf Hartmann, 2005, Wolfgang Bangerth, 2007
 */
template <>
struct RefinementPossibilities<3>
{
  /**
   * Possible values for refinement
   * cases in the current
   * dimension.
   *
   * Note the construction of the
   * values: the lowest bit
   * describes a cut of the x-axis,
   * the second to lowest bit
   * corresponds to a cut of the
   * y-axis and the third to lowest
   * bit corresponds to a cut of
   * the z-axis. Thus, the
   * following relations hold
   * (among others):
   *
   * @code
   * cut_xy  == cut_x  | cut_y
   * cut_xyz == cut_xy | cut_xz
   * cut_x   == cut_xy & cut_xz
   * @endcode
   *
   * Only those cuts that are
   * reasonable in a given space
   * dimension are offered, of
   * course.
   *
   * In addition, the tag
   * <code>isotropic_refinement</code>
   * denotes isotropic refinement
   * in the space dimension
   * selected by the template
   * argument of this class.
   */
  enum Possibilities
  {
    no_refinement= 0,
    cut_x        = 1,
    cut_y        = 2,
    cut_xy       = cut_x | cut_y,
    cut_z        = 4,
    cut_xz       = cut_x | cut_z,
    cut_yz       = cut_y | cut_z,
    cut_xyz      = cut_x | cut_y | cut_z,

    isotropic_refinement = cut_xyz
  };
};



/**
 * A class storing the possible anisotropic and isotropic refinement
 * cases of an object with <code>dim</code> dimensions (for example,
 * for a line <code>dim=1</code> in whatever space dimension we are,
 * for a quad <code>dim=2</code>, etc.). Possible values of this class
 * are the ones listed in the enumeration declared within the class.
 *
 * @ingroup aniso
 * @author Ralf Hartmann, 2005, Wolfgang Bangerth, 2007
 */
template <int dim>
class RefinementCase : public RefinementPossibilities<dim>
{
public:
  /**
   * Default constructor. Initialize the
   * refinement case with no_refinement.
   */
  RefinementCase ();

  /**
   * Constructor. Take and store a
   * value indicating a particular
   * refinement from the list of
   * possible refinements specified
   * in the base class.
   */
  RefinementCase (const typename RefinementPossibilities<dim>::Possibilities refinement_case);

  /**
   * Constructor. Take and store a
   * value indicating a particular
   * refinement as a bit field. To
   * avoid implicit conversions to
   * and from integral values, this
   * constructor is marked as
   * explicit.
   */
  explicit RefinementCase (const unsigned char refinement_case);

  /**
   * Return the numeric value
   * stored by this class. While
   * the presence of this operator
   * might seem dangerous, it is
   * useful in cases where one
   * would like to have code like
   * <tt>switch
   * (refinement_flag)... case
   * RefinementCase<dim>::cut_x:
   * ... </tt>, which can be
   * written as <code>switch
   * (static_cast@<unsigned
   * char@>(refinement_flag)</code>. Another
   * application is to use an
   * object of the current type as
   * an index into an array;
   * however, this use is
   * deprecated as it assumes a
   * certain mapping from the
   * symbolic flags defined in the
   * RefinementPossibilities base
   * class to actual numerical
   * values (the array indices).
   */
  operator unsigned char () const;

  /**
   * Return the union of the
   * refinement flags represented
   * by the current object and the
   * one given as argument.
   */
  RefinementCase operator | (const RefinementCase &r) const;

  /**
   * Return the intersection of the
   * refinement flags represented
   * by the current object and the
   * one given as argument.
   */
  RefinementCase operator & (const RefinementCase &r) const;

  /**
   * Return the negation of the
   * refinement flags represented
   * by the current object. For
   * example, in 2d, if the current
   * object holds the flag
   * <code>cut_x</code>, then the
   * returned value will be
   * <code>cut_y</code>; if the
   * current value is
   * <code>isotropic_refinement</code>
   * then the result will be
   * <code>no_refinement</code>;
   * etc.
   */
  RefinementCase operator ~ () const;


  /**
   * Return the flag that
   * corresponds to cutting a cell
   * along the axis given as
   * argument. For example, if
   * <code>i=0</code> then the
   * returned value is
   * <tt>RefinementPossibilities<dim>::cut_x</tt>.
   */
  static
  RefinementCase cut_axis (const unsigned int i);

  /**
   * Return the amount of memory
   * occupied by an object of this
   * type.
   */
  static std::size_t memory_consumption ();

  /**
   * Read or write the data of this object to or
   * from a stream for the purpose of serialization
   */
  template <class Archive>
  void serialize(Archive &ar,
                 const unsigned int version);

  /**
   * Exception.
   */
  DeclException1 (ExcInvalidRefinementCase,
                  int,
                  << "The refinement flags given (" << arg1 << ") contain set bits that do not "
                  << "make sense for the space dimension of the object to which they are applied.");

private:
  /**
   * Store the refinement case as a
   * bit field with as many bits as
   * are necessary in any given
   * dimension.
   */
unsigned char value :
  (dim > 0 ? dim : 1);
};


namespace internal
{


  /**
   * A class that provides all possible situations a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces. For <code>dim=1</code> and <code>dim=2</code> they
   * correspond to the cases given in
   * <code>RefinementPossibilities@<dim-1@></code>. However,
   * <code>SubfacePossibilities@<3@></code> includes the refinement
   * cases of <code>RefinementPossibilities@<2@></code>, but
   * additionally some subface possibilities a face might be subdivided
   * into which occur through repeated anisotropic refinement steps
   * performed on one of two neighboring cells.
   *
   * This general template is unused except in some weird template
   * constructs. Actual is made, however, of the specializations
   * <code>SubfacePossibilities@<1@></code>,
   * <code>SubfacePossibilities@<2@></code> and
   * <code>SubfacePossibilities@<3@></code>.
   *
   * @ingroup aniso
   * @author Tobias Leicht 2007, Ralf Hartmann, 2008
   */
  template <int dim>
  struct SubfacePossibilities
  {
    /**
     * Possible cases of faces
     * being subdivided into
     * subface.
     */
    enum Possibilities
    {
      case_none = 0,

      case_isotropic = static_cast<unsigned char>(-1)
    };
  };


  /**
   * A class that provides all possible situations a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces.
   *
   * For <code>dim=0</code> we provide a dummy implementation only.
   *
   * @ingroup aniso
   * @author Ralf Hartmann, 2008
   */
  template <>
  struct SubfacePossibilities<0>
  {
    /**
     * Possible cases of faces
     * being subdivided into
     * subface.
     *
     * Dummy implementation.
     */
    enum Possibilities
    {
      case_none = 0,

      case_isotropic = case_none
    };
  };



  /**
   * A class that provides all possible situations a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces.
   *
   * For <code>dim=1</code> there are no faces. Thereby, there are no
   * subface possibilities.
   *
   * @ingroup aniso
   * @author Ralf Hartmann, 2008
   */
  template <>
  struct SubfacePossibilities<1>
  {
    /**
     * Possible cases of faces
     * being subdivided into
     * subface.
     *
     * In 1d there are no faces,
     * thus no subface
     * possibilities.
     */
    enum Possibilities
    {
      case_none = 0,

      case_isotropic = case_none
    };
  };



  /**
   * A class that provides all possible situations a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces.
   *
   * This specialization is used for <code>dim=2</code>, where it offers
   * the following possibilities: a face (line) being refined
   * (<code>case_x</code>) or not refined (<code>case_no</code>).
   *
   * @ingroup aniso
   * @author Ralf Hartmann, 2008
   */
  template <>
  struct SubfacePossibilities<2>
  {
    /**
     * Possible cases of faces
     * being subdivided into
     * subface.
     *
     * In 2d there are following
     * possibilities: a face (line)
     * being refined *
     * (<code>case_x</code>) or not
     * refined
     * (<code>case_no</code>).
    */
    enum Possibilities
    {
      case_none = 0,
      case_x    = 1,

      case_isotropic = case_x
    };
  };



  /**
   * A class that provides all possible situations a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces.
   *
   * This specialization is used for dim=3, where it offers following
   * possibilities: a face (quad) being refined in x- or y-direction (in
   * the face-intern coordinate system) separately, (<code>case_x</code>
   * or (<code>case_y</code>), and in both directions
   * (<code>case_x</code> which corresponds to
   * (<code>case_isotropic</code>). Additionally, it offers the
   * possibilities a face can have through repeated anisotropic
   * refinement steps performed on one of the two neighboring cells.  It
   * might be possible for example, that a face (quad) is refined with
   * <code>cut_x</code> and afterwards the left child is again refined
   * with <code>cut_y</code>, so that there are three active
   * subfaces. Note, however, that only refinement cases are allowed
   * such that each line on a face between two hexes has not more than
   * one hanging node. Furthermore, it is not allowed that two
   * neighboring hexes are refined such that one of the hexes refines
   * the common face with <code>cut_x</code> and the other hex refines
   * that face with <code>cut_y</code>. In fact,
   * Triangulation::prepare_coarsening_and_refinement takes care of this
   * situation and ensures that each face of a refined cell is
   * completely contained in a single face of neighboring cells.
   *
   * The following drawings explain the SubfacePossibilities and give
   * the corresponding subface numbers:
   * @code

   *-------*
   |       |
   |   0   |    case_none
   |       |
   *-------*

   *---*---*
   |   |   |
   | 0 | 1 |    case_x
   |   |   |
   *---*---*

   *---*---*
   | 1 |   |
   *---* 2 |    case_x1y
   | 0 |   |
   *---*---*

   *---*---*
   |   | 2 |
   | 0 *---*    case_x2y
   |   | 1 |
   *---*---*

   *---*---*
   | 1 | 3 |
   *---*---*    case_x1y2y   (successive refinement: first cut_x, then cut_y for both children)
   | 0 | 2 |
   *---*---*

   *-------*
   |   1   |
   *-------*    case_y
   |   0   |
   *-------*

   *-------*
   |   2   |
   *---*---*    case_y1x
   | 0 | 1 |
   *---*---*

   *---*---*
   | 1 | 2 |
   *---*---*    case_y2x
   |   0   |
   *-------*

   *---*---*
   | 2 | 3 |
   *---*---*    case_y1x2x   (successive refinement: first cut_y, then cut_x for both children)
   | 0 | 1 |
   *---+---*

   *---*---*
   | 2 | 3 |
   *---*---*    case_xy      (one isotropic refinement step)
   | 0 | 1 |
   *---*---*

   * @endcode
   *
   * @ingroup aniso
   * @author Tobias Leicht 2007, Ralf Hartmann, 2008
   */
  template <>
  struct SubfacePossibilities<3>
  {
    /**
     * Possible cases of faces
     * being subdivided into
     * subface.
     *
     * See documentation to the
     * SubfacePossibilities<3> for
     * more details on the subface
     * possibilities.
     */
    enum Possibilities
    {
      case_none  = 0,
      case_x     = 1,
      case_x1y   = 2,
      case_x2y   = 3,
      case_x1y2y = 4,
      case_y     = 5,
      case_y1x   = 6,
      case_y2x   = 7,
      case_y1x2x = 8,
      case_xy    = 9,

      case_isotropic = case_xy
    };
  };




  /**
   * A class that provides all possible cases a face (in the
   * current space dimension @p dim) might be subdivided into
   * subfaces.
   *
   * @ingroup aniso
   * @author Ralf Hartmann, 2008
   */
  template <int dim>
  class SubfaceCase : public SubfacePossibilities<dim>
  {
  public:
    /**
     * Constructor. Take and store
     * a value indicating a
     * particular subface
     * possibility in the list of
     * possible situations
     * specified in the base class.
     */
    SubfaceCase (const typename SubfacePossibilities<dim>::Possibilities subface_possibility);

    /**
     * Return the numeric value
     * stored by this class. While
     * the presence of this operator
     * might seem dangerous, it is
     * useful in cases where one
     * would like to have code like
     * <code>switch
     * (subface_case)... case
     * SubfaceCase::case_x:
     * ... </code>, which can be
     * written as <code>switch
     * (static_cast@<unsigned
     * char@>(subface_case)</code>. Another
     * application is to use an
     * object of the current type as
     * an index into an array;
     * however, this use is
     * deprecated as it assumes a
     * certain mapping from the
     * symbolic flags defined in the
     * SubfacePossibilities
     * base class to actual numerical
     * values (the array indices).
     */
    operator unsigned char () const;

    /**
     * Return the amount of memory
     * occupied by an object of this
     * type.
     */
    static std::size_t memory_consumption ();

    /**
     * Exception.
     */
    DeclException1 (ExcInvalidSubfaceCase,
                    int,
                    << "The subface case given (" << arg1 << ") does not make sense "
                    << "for the space dimension of the object to which they are applied.");

  private:
    /**
     * Store the refinement case as a
     * bit field with as many bits as
     * are necessary in any given
     * dimension.
     */
  unsigned char value :
    (dim == 3 ? 4 : 1);
  };

} // namespace internal



template <int dim> struct GeometryInfo;




/**
 * Topological description of zero dimensional cells,
 * i.e. points. This class might not look too useful but often is if
 * in a certain dimension we would like to enquire information about
 * objects with dimension one lower than the present, e.g. about
 * faces.
 *
 * This class contains as static members information on vertices and
 * faces of a @p dim-dimensional grid cell. The interface is the same
 * for all dimensions. If a value is of no use in a low dimensional
 * cell, it is (correctly) set to zero, e.g. #max_children_per_face in
 * 1d.
 *
 * This information should always replace hard-coded numbers of
 * vertices, neighbors and so on, since it can be used dimension
 * independently.
 *
 * @ingroup grid geomprimitives aniso
 * @author Wolfgang Bangerth, 1998
 */
template <>
struct GeometryInfo<0>
{

  /**
   * Maximum number of children of
   * a cell, i.e. the number of
   * children of an isotropically
   * refined cell.
   *
   * If a cell is refined
   * anisotropically, the actual
   * number of children may be less
   * than the value given here.
   */
  static const unsigned int max_children_per_cell = 1;

  /**
   * Number of faces a cell has.
   */
  static const unsigned int faces_per_cell    = 0;

  /**
   * Maximum number of children of
   * a refined face, i.e. the
   * number of children of an
   * isotropically refined face.
   *
   * If a cell is refined
   * anisotropically, the actual
   * number of children may be less
   * than the value given here.
   */
  static const unsigned int max_children_per_face = 0;

  /**
   * Return the number of children
   * of a cell (or face) refined
   * with <tt>ref_case</tt>. Since
   * we are concerned here with
   * points, the number of children
   * is equal to one.
   */
  static unsigned int n_children(const RefinementCase<0> &refinement_case);

  /**
   * Number of vertices a cell has.
   */
  static const unsigned int vertices_per_cell = 1;

  /**
   * Number of vertices each face has.
   * Since this is not useful in one
   * dimension, we provide a useless
   * number (in the hope that a compiler
   * may warn when it sees constructs like
   * <tt>for (i=0; i<vertices_per_face; ++i)</tt>,
   * at least if @p i is an <tt>unsigned int</tt>.
   */
  static const unsigned int vertices_per_face = 0;

  /**
   * Number of lines each face has.
   */
  static const unsigned int lines_per_face    = 0;

  /**
   * Number of quads on each face.
   */
  static const unsigned int quads_per_face    = 0;

  /**
   * Number of lines of a cell.
   */
  static const unsigned int lines_per_cell    = 0;

  /**
   * Number of quadrilaterals of a
   * cell.
   */
  static const unsigned int quads_per_cell    = 0;

  /**
   * Number of hexahedra of a
   * cell.
   */
  static const unsigned int hexes_per_cell    = 0;
};





/**
 * This class provides dimension independent information to all topological
 * structures that make up the unit, or
 * @ref GlossReferenceCell "reference cell".
 *
 * It is the one central point in the library where information about the
 * numbering of vertices, lines, or faces of the reference cell is
 * collected. Consequently, the information of this class is used extensively
 * in the geometric description of Triangulation objects, as well as in
 * various other parts of the code. In particular, it also serves as the focus
 * of writing code in a dimension independent way; for example, instead of
 * writing a loop over vertices 0<=v<4 in 2d, one would write it as
 * 0<=v<GeometryInfo<dim>::vertices_per_cell, thus allowing the code to work
 * in 3d as well without changes.
 *
 * The most frequently used parts of the class are its static members like
 * vertices_per_cell, faces_per_cell, etc. However, the class also offers
 * information about more abstract questions like the orientation of faces,
 * etc. The following documentation gives a textual description of many of
 * these concepts.
 *
 *
 * <h3>Implementation conventions for two spatial dimensions</h3>
 *
 * From version 5.2 onwards deal.II is based on a numbering scheme
 * that uses a lexicographic ordering (with x running fastest)
 * wherever possible, hence trying to adopt a kind of 'canonical'
 * ordering.
 *
 * The ordering of vertices and faces (lines) in 2d is defined by
 *
 * - Vertices are numbered in lexicographic ordering
 *
 * - Faces (lines in 2d): first the two faces with normals in x-
 *   and then y-direction. For each two faces: first the face with
 *   normal in negative coordinate direction, then the one with normal
 *   in positive direction, i.e. the faces are ordered according to
 *   their normals pointing in -x, x, -y, y direction.
 *
 * - The direction of a line is represented by the direction of
 *   point 0 towards point 1 and is always in one of the coordinate
 *   directions
 *
 * - Face lines in 3d are ordered, such that the induced 2d local
 *   coordinate system (x,y) implies (right hand rule) a normal in
 *   face normal direction, see N2/.
 *
 * The resulting numbering of vertices and faces (lines) in 2d as
 * well as the directions of lines is shown in the following.
 * @verbatim
 *       3
 *    2-->--3
 *    |     |
 *   0^     ^1
 *    |     |
 *    0-->--1
 *        2
 * @endverbatim
 *
 * Note that the orientation of lines has to be correct upon construction of a
 * grid; however, it is automatically preserved upon refinement.
 *
 * Further we define that child lines have the same direction as their parent,
 * i.e. that <tt>line->child(0)->vertex(0)==line->vertex(0)</tt> and
 * <tt>line->child(1)->vertex(1)==line->vertex(1)</tt>. This also implies,
 * that the first sub-line (<tt>line->child(0)</tt>) is the one at vertex(0)
 * of the old line.
 *
 * Similarly we define, that the four children of a quad are adjacent to the
 * vertex with the same number of the old quad.
 *
 * Note that information about several of these conventions can be
 * extracted at run- or compile-time from the member functions and
 * variables of the present class.
 *
 *
 * <h4>Coordinate systems</h4>
 *
 * When explicit coordinates are required for points in a cell (e.g for
 * quadrature formulae or the point of definition of trial functions), we
 * define the following coordinate system for the unit cell:
 * @verbatim
 *  y^   2-----3
 *   |   |     |
 *   |   |     |
 *   |   |     |
 *   |   0-----1
 *   *------------>x
 * @endverbatim
 *
 * Here, vertex 0 is the origin of the coordinate system, vertex 1 has
 * coordinates <tt>(1,0)</tt>, vertex 2 at <tt>(0,1)</tt> and vertex 3 at
 * <tt>(1,1)</tt>. The GeometryInfo<dim>::unit_cell_vertex() function can be
 * used to query this information at run-time.
 *
 *
 * <h3>Implementation conventions for three spatial dimensions</h3>
 *
 * By convention, we will use the following numbering conventions
 * for vertices, lines and faces of hexahedra in three space
 * dimensions. Before giving these conventions we declare the
 * following sketch to be the standard way of drawing 3d pictures of
 * hexahedra:
 * @verbatim
 *                       *-------*        *-------*
 *                      /|       |       /       /|
 *                     / |       |      /       / |
 *  z                 /  |       |     /       /  |
 *  ^                *   |       |    *-------*   |
 *  |   ^y           |   *-------*    |       |   *
 *  |  /             |  /       /     |       |  /
 *  | /              | /       /      |       | /
 *  |/               |/       /       |       |/
 *  *------>x        *-------*        *-------*
 * @endverbatim
 * The left part of the picture shows the left, bottom and back face of the
 * cube, while the right one shall be the top, right and front face. You may
 * recover the whole cube by moving the two parts together into one.
 *
 * Note again that information about several of the following
 * conventions can be extracted at run- or compile-time from the
 * member functions and variables of the present class.
 *
 * <h4>Vertices</h4>
 *
 * The ordering of vertices in 3d is defined by the same rules as in
 * the 2d case. In particular, the following is still true:
 *
 * - Vertices are numbered in lexicographic ordering.
 *
 * Hence, the vertices are numbered as follows
 * @verbatim
 *       6-------7        6-------7
 *      /|       |       /       /|
 *     / |       |      /       / |
 *    /  |       |     /       /  |
 *   4   |       |    4-------5   |
 *   |   2-------3    |       |   3
 *   |  /       /     |       |  /
 *   | /       /      |       | /
 *   |/       /       |       |/
 *   0-------1        0-------1
 * @endverbatim
 *
 * We note, that first the vertices on the bottom face (z=0) are numbered
 * exactly the same way as are the vertices on a quadrilateral. Then the
 * vertices on the top face (z=1) are numbered similarly by moving the bottom
 * face to the top. Again, the GeometryInfo<dim>::unit_cell_vertex() function
 * can be used to query this information at run-time.
 *
 *
 * <h4>Lines</h4>
 *
 * Here, the same holds as for the vertices:
 *
 * - Line ordering in 3d:
 *   <ul>
 *   <li>first the lines of face (z=0) in 2d line ordering,
 *   <li>then the lines of face (z=1) in 2d line ordering,
 *   <li>finally the lines in z direction in lexicographic ordering
 *   </ul>
 * @verbatim
 *       *---7---*        *---7---*
 *      /|       |       /       /|
 *     4 |       11     4       5 11
 *    /  10      |     /       /  |
 *   *   |       |    *---6---*   |
 *   |   *---3---*    |       |   *
 *   |  /       /     |       9  /
 *   8 0       1      8       | 1
 *   |/       /       |       |/
 *   *---2---*        *---2---*
 * @endverbatim
 * As in 2d lines are directed in coordinate directions.
 * @verbatim
 *       *--->---*        *--->---*
 *      /|       |       /       /|
 *     ^ |       ^      ^       ^ ^
 *    /  ^       |     /       /  |
 *   *   |       |    *--->---*   |
 *   |   *--->---*    |       |   *
 *   |  /       /     |       ^  /
 *   ^ ^       ^      ^       | ^
 *   |/       /       |       |/
 *   *--->---*        *--->---*
 * @endverbatim
 *
 * The fact that edges (just as vertices and faces) are entities that
 * are stored in their own right rather than constructed from cells
 * each time they are needed, means that adjacent cells actually have
 * pointers to edges that are thus shared between them. This implies
 * that the convention that sets of parallel edges have parallel
 * directions is not only a local condition. Before a list of cells is
 * passed to an object of the Triangulation class for creation of a
 * triangulation, you therefore have to make sure that cells are
 * oriented in a compatible fashion, so that edge directions are
 * globally according to above convention. However, the GridReordering
 * class can do this for you, by reorienting cells and edges of an
 * arbitrary list of input cells that need not be already sorted.
 *
 * <h4>Faces</h4>
 *
 * The numbering of faces in 3d is defined by a rule analogous to 2d:
 *
 * - Faces (quads in 3d): first the two faces with normals in x-,
 * then y- and z-direction. For each two faces: first the face with
 * normal in negative coordinate direction, then the one with normal
 * in positive direction, i.e. the faces are ordered according to
 * their normals pointing in -x, x, -y, y, -z, z direction.
 *
 * Therefore, the faces are numbered in the ordering: left, right,
 * front, back, bottom and top face:
 * @verbatim
 *       *-------*        *-------*
 *      /|       |       /       /|
 *     / |   3   |      /   5   / |
 *    /  |       |     /       /  |
 *   *   |       |    *-------*   |
 *   | 0 *-------*    |       | 1 *
 *   |  /       /     |       |  /
 *   | /   4   /      |   2   | /
 *   |/       /       |       |/
 *   *-------*        *-------*
 * @endverbatim
 *
 * The <em>standard</em> direction of the faces is such, that the
 * induced 2d local coordinate system (x,y) implies (right hand
 * rule) a normal in face normal direction, see N2a).  In the
 * following we show the local coordinate system and the numbering
 * of face lines:
 * <ul>
 * <li> Faces 0 and 1:
 *  @verbatim
 *          Face 0           Face 1
 *        *-------*        *-------*
 *       /|       |       /       /|
 *      3 1       |      /       3 1
 *    y/  |       |     /      y/  |
 *    *   |x      |    *-------*   |x
 *    |   *-------*    |       |   *
 *    0  /       /     |       0  /
 *    | 2       /      |       | 2
 *    |/       /       |       |/
 *    *-------*        *-------*
 *  @endverbatim
 *
 * <li> Faces 2 and 3:
 *  @verbatim
 *        x Face 3           Face 2
 *        *---1---*        *-------*
 *       /|       |       /       /|
 *      / |       3      /       / |
 *     /  2       |    x/       /  |
 *    *   |       |    *---1---*   |
 *    |   *---0---*y   |       |   *
 *    |  /       /     |       3  /
 *    | /       /      2       | /
 *    |/       /       |       |/
 *    *-------*        *---0---*y
 *  @endverbatim
 *
 * <li> Faces 4 and 5:
 *  @verbatim
 *          Face 4         y Face 5
 *        *-------*        *---3---*
 *       /|       |       /       /|
 *      / |       |      0       1 |
 *     /  |       |     /       /  |
 *    *   |y      |    *---2---* x |
 *    |   *---3---*    |       |   *
 *    |  /       /     |       |  /
 *    | 0       1      |       | /
 *    |/       /       |       |/
 *    *---2---* x      *-------*
 *  @endverbatim
 * </ul>
 *
 * The face line numbers (0,1,2,3) correspond to following cell line
 * numbers.
 * <ul>
 * <li> Face 0: lines 8, 10, 0, 4;
 * <li> Face 1: lines 9, 11, 1, 5;
 * <li> Face 2: lines 2, 6, 8, 9;
 * <li> Face 3: lines 3, 7, 10, 11;
 * <li> Face 4: lines 0, 1, 2, 3;
 * <li> Face 5: lines 4, 5, 6, 7;
 * </ul>
 * You can get these numbers using the
 * GeometryInfo<3>::face_to_cell_lines() function.
 *
 * The face normals can be deduced from the face orientation by
 * applying the right hand side rule (x,y -> normal).  We note, that
 * in the standard orientation of faces in 2d, faces 0 and 2 have
 * normals that point into the cell, and faces 1 and 3 have normals
 * pointing outward. In 3d, faces 0, 2, and 4
 * have normals that point into the cell, while the normals of faces
 * 1, 3, and 5 point outward. This information, again, can be queried from
 * GeometryInfo<dim>::unit_normal_orientation.
 *
 * However, it turns out that a significant number of 3d meshes cannot
 * satisfy this convention. This is due to the fact that the face
 * convention for one cell already implies something for the
 * neighbor, since they share a common face and fixing it for the
 * first cell also fixes the normal vectors of the opposite faces of
 * both cells. It is easy to construct cases of loops of cells for
 * which this leads to cases where we cannot find orientations for
 * all faces that are consistent with this convention.
 *
 * For this reason, above convention is only what we call the <em>standard
 * orientation</em>. deal.II actually allows faces in 3d to have either the
 * standard direction, or its opposite, in which case the lines that make up a
 * cell would have reverted orders, and the above line equivalences would not
 * hold any more. You can ask a cell whether a given face has standard
 * orientation by calling <tt>cell->face_orientation(face_no)</tt>: if the
 * result is @p true, then the face has standard orientation, otherwise its
 * normal vector is pointing the other direction. There are not very many
 * places in application programs where you need this information actually,
 * but a few places in the library make use of this. Note that in 2d, the
 * result is always @p true. More information on the topic can be found in this
 * @ref GlossFaceOrientation "glossary" article.
 *
 * In order to allow all kinds of meshes in 3d, including
 * <em>Moebius</em>-loops, a face might even be rotated looking
 * from one cell, whereas it is according to the standard looking at it from the
 * neighboring cell sharing that particular face. In order to cope with this,
 * two flags <tt>face_flip</tt> and <tt>face_rotation</tt> are available, to
 * represent rotations by 180 and 90 degree, respectively. Setting both flags
 * amounts to a rotation of 270 degrees (all counterclockwise). You can ask
 * the cell for these flags like for the <tt>face_orientation</tt>. In order to
 * enable rotated faces, even lines can deviate from their standard direction in
 * 3d. This information is available as the <tt>line_orientation</tt> flag for
 * cells and faces in 3d. Again, this is something that should be internal to
 * the library and application program will probably never have to bother about
 * it. For more information on this see also
 * @ref GlossFaceOrientation "this glossary entry" .
 *
 *
 * <h4>Children</h4>
 *
 * The eight children of an isotropically refined cell are numbered according to
 * the vertices they are adjacent to:
 * @verbatim
 *       *----*----*        *----*----*
 *      /| 6  |  7 |       / 6  /  7 /|
 *     *6|    |    |      *----*----*7|
 *    /| *----*----*     / 4  /  5 /| *
 *   * |/|    |    |    *----*----* |/|
 *   |4* | 2  |  3 |    | 4  |  5 |5*3|
 *   |/|2*----*----*    |    |    |/| *
 *   * |/ 2  /  3 /     *----*----* |/
 *   |0*----*----*      |    |    |1*
 *   |/0   /  1 /       | 0  |  1 |/
 *   *----*----*        *----*----*
 * @endverbatim
 *
 * Taking into account the orientation of the faces, the following
 * children are adjacent to the respective faces:
 * <ul>
 * <li> Face 0: children 0, 2, 4, 6;
 * <li> Face 1: children 1, 3, 5, 7;
 * <li> Face 2: children 0, 4, 1, 5;
 * <li> Face 3: children 2, 6, 3, 7;
 * <li> Face 4: children 0, 1, 2, 3;
 * <li> Face 5: children 4, 5, 6, 7.
 * </ul>
 * You can get these numbers using the
 * GeometryInfo<3>::child_cell_on_face() function. As each child is
 * adjacent to the vertex with the same number these numbers are
 * also given by the GeometryInfo<3>::face_to_cell_vertices()
 * function.
 *
 * Note that, again, the above list only holds for faces in their
 * standard orientation. If a face is not in standard orientation,
 * then the children at positions 1 and 2 (counting from 0 to 3)
 * would be swapped. In fact, this is what the child_cell_on_face
 * and the face_to_cell_vertices functions of GeometryInfo<3> do,
 * when invoked with a <tt>face_orientation=false</tt> argument.
 *
 * The information which child cell is at which position of which face
 * is most often used when computing jump terms across faces with
 * hanging nodes, using objects of type FESubfaceValues. Sitting on
 * one cell, you would look at a face and figure out which child of
 * the neighbor is sitting on a given subface between the present and
 * the neighboring cell. To avoid having to query the standard
 * orientation of the faces of the two cells every time in such cases,
 * you should use a function call like
 * <tt>cell->neighbor_child_on_subface(face_no,subface_no)</tt>, which
 * returns the correct result both in 2d (where face orientations are
 * immaterial) and 3d (where it is necessary to use the face
 * orientation as additional argument to
 * <tt>GeometryInfo<3>::child_cell_on_face</tt>).
 *
 * For anisotropic refinement, the child cells can not be numbered according to
 * adjacent vertices, thus the following conventions are used:
 * @verbatim
 *            RefinementCase<3>::cut_x
 *
 *       *----*----*        *----*----*
 *      /|    |    |       /    /    /|
 *     / |    |    |      / 0  /  1 / |
 *    /  | 0  |  1 |     /    /    /  |
 *   *   |    |    |    *----*----*   |
 *   | 0 |    |    |    |    |    | 1 |
 *   |   *----*----*    |    |    |   *
 *   |  /    /    /     | 0  | 1  |  /
 *   | / 0  /  1 /      |    |    | /
 *   |/    /    /       |    |    |/
 *   *----*----*        *----*----*
 * @endverbatim
 *
 * @verbatim
 *            RefinementCase<3>::cut_y
 *
 *       *---------*        *---------*
 *      /|         |       /    1    /|
 *     * |         |      *---------* |
 *    /| |    1    |     /    0    /| |
 *   * |1|         |    *---------* |1|
 *   | | |         |    |         | | |
 *   |0| *---------*    |         |0| *
 *   | |/    1    /     |    0    | |/
 *   | *---------*      |         | *
 *   |/    0    /       |         |/
 *   *---------*        *---------*
 * @endverbatim
 *
 * @verbatim
 *            RefinementCase<3>::cut_z
 *
 *       *---------*        *---------*
 *      /|    1    |       /         /|
 *     / |         |      /    1    / |
 *    /  *---------*     /         /  *
 *   * 1/|         |    *---------* 1/|
 *   | / |    0    |    |    1    | / |
 *   |/  *---------*    |         |/  *
 *   * 0/         /     *---------* 0/
 *   | /    0    /      |         | /
 *   |/         /       |    0    |/
 *   *---------*        *---------*
 * @endverbatim
 *
 * @verbatim
 *            RefinementCase<3>::cut_xy
 *
 *       *----*----*        *----*----*
 *      /|    |    |       / 2  /  3 /|
 *     * |    |    |      *----*----* |
 *    /| | 2  |  3 |     / 0  /  1 /| |
 *   * |2|    |    |    *----*----* |3|
 *   | | |    |    |    |    |    | | |
 *   |0| *----*----*    |    |    |1| *
 *   | |/ 2  /  3 /     | 0  |  1 | |/
 *   | *----*----*      |    |    | *
 *   |/ 0  /  1 /       |    |    |/
 *   *----*----*        *----*----*
 * @endverbatim
 *
 * @verbatim
 *            RefinementCase<3>::cut_xz
 *
 *       *----*----*        *----*----*
 *      /| 1  |  3 |       /    /    /|
 *     / |    |    |      / 1  /  3 / |
 *    /  *----*----*     /    /    /  *
 *   * 1/|    |    |    *----*----* 3/|
 *   | / | 0  |  2 |    | 1  |  3 | / |
 *   |/  *----*----*    |    |    |/  *
 *   * 0/    /    /     *----*----* 2/
 *   | / 0  /  2 /      |    |    | /
 *   |/    /    /       | 0  |  2 |/
 *   *----*----*        *----*----*
 * @endverbatim
 *
 * @verbatim
 *            RefinementCase<3>::cut_yz
 *
 *       *---------*        *---------*
 *      /|    3    |       /    3    /|
 *     * |         |      *---------* |
 *    /|3*---------*     /    2    /|3*
 *   * |/|         |    *---------* |/|
 *   |2* |    1    |    |    2    |2* |
 *   |/|1*---------*    |         |/|1*
 *   * |/    1    /     *---------* |/
 *   |0*---------*      |         |0*
 *   |/    0    /       |    0    |/
 *   *---------*        *---------*
 * @endverbatim
 *
 * This information can also be obtained by the
 * <tt>GeometryInfo<3>::child_cell_on_face</tt> function.
 *
 * <h4>Coordinate systems</h4>
 *
 * We define the following coordinate system for the explicit coordinates of
 * the vertices of the unit cell:
 * @verbatim
 *                       6-------7        6-------7
 *                      /|       |       /       /|
 *                     / |       |      /       / |
 *  z                 /  |       |     /       /  |
 *  ^                4   |       |    4-------5   |
 *  |   ^y           |   2-------3    |       |   3
 *  |  /             |  /       /     |       |  /
 *  | /              | /       /      |       | /
 *  |/               |/       /       |       |/
 *  *------>x        0-------1        0-------1
 * @endverbatim
 *
 * By the convention laid down as above, the vertices have the following
 * coordinates (lexicographic, with x running fastest):
 * <ul>
 *    <li> Vertex 0: <tt>(0,0,0)</tt>;
 *    <li> Vertex 1: <tt>(1,0,0)</tt>;
 *    <li> Vertex 2: <tt>(0,1,0)</tt>;
 *    <li> Vertex 3: <tt>(1,1,0)</tt>;
 *    <li> Vertex 4: <tt>(0,0,1)</tt>;
 *    <li> Vertex 5: <tt>(1,0,1)</tt>;
 *    <li> Vertex 6: <tt>(0,1,1)</tt>;
 *    <li> Vertex 7: <tt>(1,1,1)</tt>.
 * </ul>
 *
 *
 *
 * @note Instantiations for this template are provided for dimensions 1,2,3,4,
 * and there is a specialization for dim=0 (see the section on @ref
 * Instantiations in the manual).
 *
 * @ingroup grid geomprimitives aniso
 * @author Wolfgang Bangerth, 1998, Ralf Hartmann, 2005, Tobias Leicht, 2007
 */
template <int dim>
struct GeometryInfo
{

  /**
   * Maximum number of children of
   * a refined cell, i.e. the
   * number of children of an
   * isotropically refined cell.
   *
   * If a cell is refined
   * anisotropically, the actual
   * number of children may be less
   * than the value given here.
   */
  static const unsigned int max_children_per_cell = 1 << dim;

  /**
   * Number of faces of a cell.
   */
  static const unsigned int faces_per_cell = 2 * dim;

  /**
   * Maximum number of children of
   * a refined face, i.e. the
   * number of children of an
   * isotropically refined face.
   *
   * If a cell is refined
   * anisotropically, the actual
   * number of children may be less
   * than the value given here.
   */
  static const unsigned int max_children_per_face = GeometryInfo<dim-1>::max_children_per_cell;

  /**
   * Number of vertices of a cell.
   */
  static const unsigned int vertices_per_cell = 1 << dim;

  /**
   * Number of vertices on each
   * face.
   */
  static const unsigned int vertices_per_face = GeometryInfo<dim-1>::vertices_per_cell;

  /**
   * Number of lines on each face.
   */
  static const unsigned int lines_per_face
    = GeometryInfo<dim-1>::lines_per_cell;

  /**
   * Number of quads on each face.
   */
  static const unsigned int quads_per_face
    = GeometryInfo<dim-1>::quads_per_cell;

  /**
   * Number of lines of a cell.
   *
   * The formula to compute this makes use
   * of the fact that when going from one
   * dimension to the next, the object of
   * the lower dimension is copied once
   * (thus twice the old number of lines)
   * and then a new line is inserted
   * between each vertex of the old object
   * and the corresponding one in the copy.
   */
  static const unsigned int lines_per_cell
    = (2*GeometryInfo<dim-1>::lines_per_cell +
       GeometryInfo<dim-1>::vertices_per_cell);

  /**
   * Number of quadrilaterals of a
   * cell.
   *
   * This number is computed recursively
   * just as the previous one, with the
   * exception that new quads result from
   * connecting an original line and its
   * copy.
   */
  static const unsigned int quads_per_cell
    = (2*GeometryInfo<dim-1>::quads_per_cell +
       GeometryInfo<dim-1>::lines_per_cell);

  /**
   * Number of hexahedra of a
   * cell.
   */
  static const unsigned int hexes_per_cell
    = (2*GeometryInfo<dim-1>::hexes_per_cell +
       GeometryInfo<dim-1>::quads_per_cell);

  /**
   * Rearrange vertices for UCD
   * output.  For a cell being
   * written in UCD format, each
   * entry in this field contains
   * the number of a vertex in
   * <tt>deal.II</tt> that corresponds
   * to the UCD numbering at this
   * location.
   *
   * Typical example: write a cell
   * and arrange the vertices, such
   * that UCD understands them.
   *
   * @code
   * for (i=0; i< n_vertices; ++i)
   *   out << cell->vertex(ucd_to_deal[i]);
   * @endcode
   *
   * As the vertex numbering in
   * deal.II versions <= 5.1
   * happened to coincide with the
   * UCD numbering, this field can
   * also be used like a
   * old_to_lexicographic mapping.
   */
  static const unsigned int ucd_to_deal[vertices_per_cell];

  /**
   * Rearrange vertices for OpenDX
   * output.  For a cell being
   * written in OpenDX format, each
   * entry in this field contains
   * the number of a vertex in
   * <tt>deal.II</tt> that corresponds
   * to the DX numbering at this
   * location.
   *
   * Typical example: write a cell
   * and arrange the vertices, such
   * that OpenDX understands them.
   *
   * @code
   * for (i=0; i< n_vertices; ++i)
   *   out << cell->vertex(dx_to_deal[i]);
   * @endcode
   */
  static const unsigned int dx_to_deal[vertices_per_cell];

  /**
   * This field stores for each vertex
   * to which faces it belongs. In any
   * given dimension, the number of
   * faces is equal to the dimension.
   * The first index in this 2D-array
   * runs over all vertices, the second
   * index over @p dim faces to which
   * the vertex belongs.
  *
  * The order of the faces for
  * each vertex is such that the
  * first listed face bounds the
  * reference cell in <i>x</i>
  * direction, the second in
  * <i>y</i> direction, and so on.
   */
  static const unsigned int vertex_to_face[vertices_per_cell][dim];

  /**
   * Return the number of children
   * of a cell (or face) refined
   * with <tt>ref_case</tt>.
   */
  static
  unsigned int
  n_children(const RefinementCase<dim> &refinement_case);

  /**
   * Return the number of subfaces
   * of a face refined according to
   * internal::SubfaceCase
   * @p face_ref_case.
   */
  static
  unsigned int
  n_subfaces(const internal::SubfaceCase<dim> &subface_case);

  /**
   * Given a face on the reference
   * element with a
   * <code>internal::SubfaceCase@<dim@></code>
   * @p face_refinement_case this
   * function returns the ratio
   * between the area of the @p
   * subface_no th subface and the
   * area(=1) of the face.
   *
   * E.g. for
   * internal::SubfaceCase::cut_xy
   * the ratio is 1/4 for each of
   * the subfaces.
   */
  static
  double
  subface_ratio(const internal::SubfaceCase<dim> &subface_case,
                const unsigned int subface_no);

  /**
   * Given a cell refined with the
   * <code>RefinementCase</code>
   * @p cell_refinement_case
   * return the
   * <code>SubfaceCase</code> of
   * the @p face_no th face.
   */
  static
  RefinementCase<dim-1>
  face_refinement_case (const RefinementCase<dim> &cell_refinement_case,
                        const unsigned int face_no,
                        const bool face_orientation = true,
                        const bool face_flip        = false,
                        const bool face_rotation    = false);

  /**
   * Given the SubfaceCase @p
   * face_refinement_case of the @p
   * face_no th face, return the
   * smallest RefinementCase of the
   * cell, which corresponds to
   * that refinement of the face.
   */
  static
  RefinementCase<dim>
  min_cell_refinement_case_for_face_refinement
  (const RefinementCase<dim-1> &face_refinement_case,
   const unsigned int face_no,
   const bool face_orientation = true,
   const bool face_flip        = false,
   const bool face_rotation    = false);

  /**
   * Given a cell refined with the
   * RefinementCase @p
   * cell_refinement_case return
   * the RefinementCase of the @p
   * line_no th face.
   */
  static
  RefinementCase<1>
  line_refinement_case(const RefinementCase<dim> &cell_refinement_case,
                       const unsigned int line_no);

  /**
   * Return the minimal / smallest
   * RefinementCase of the cell, which
   * ensures refinement of line
   * @p line_no.
   */
  static
  RefinementCase<dim>
  min_cell_refinement_case_for_line_refinement(const unsigned int line_no);

  /**
   * This field stores which child
   * cells are adjacent to a
   * certain face of the mother
   * cell.
   *
   * For example, in 2D the layout of
   * a cell is as follows:
   * @verbatim
   * .      3
   * .   2-->--3
   * .   |     |
   * . 0 ^     ^ 1
   * .   |     |
   * .   0-->--1
   * .      2
   * @endverbatim
   * Vertices and faces are indicated
   * with their numbers, faces also with
   * their directions.
   *
   * Now, when refined, the layout is
   * like this:
   * @verbatim
   * *--*--*
   * | 2|3 |
   * *--*--*
   * | 0|1 |
   * *--*--*
   * @endverbatim
   *
   * Thus, the child cells on face
   * 0 are (ordered in the
   * direction of the face) 0 and
   * 2, on face 3 they are 2 and 3,
   * etc.
   *
   * For three spatial dimensions, the exact
   * order of the children is laid down in
   * the general documentation of this
   * class.
   *
   * Through the <tt>face_orientation</tt>,
   * <tt>face_flip</tt> and
   * <tt>face_rotation</tt> arguments this
   * function handles faces oriented in the
   * standard and non-standard orientation.
   * <tt>face_orientation</tt> defaults to
   * <tt>true</tt>, <tt>face_flip</tt> and
   * <tt>face_rotation</tt> default to
   * <tt>false</tt> (standard orientation)
   * and has no effect in 2d. The concept of
   * face orientations is explained in this
   * @ref GlossFaceOrientation "glossary"
   * entry.
   *
   * In the case of anisotropically refined
   * cells and faces, the @p RefinementCase of
   * the face, <tt>face_ref_case</tt>,
   * might have an influence on
   * which child is behind which given
   * subface, thus this is an additional
   * argument, defaulting to isotropic
   * refinement of the face.
   */
  static
  unsigned int
  child_cell_on_face (const RefinementCase<dim> &ref_case,
                      const unsigned int face,
                      const unsigned int subface,
                      const bool face_orientation = true,
                      const bool face_flip        = false,
                      const bool face_rotation    = false,
                      const RefinementCase<dim-1> &face_refinement_case
                      = RefinementCase<dim-1>::isotropic_refinement);

  /**
   * Map line vertex number to cell
   * vertex number, i.e. give the
   * cell vertex number of the
   * <tt>vertex</tt>th vertex of
   * line <tt>line</tt>, e.g.
   * <tt>GeometryInfo<2>::line_to_cell_vertices(3,0)=2</tt>.
   *
   * The order of the lines, as well as
   * their direction (which in turn
   * determines which is the first and
   * which the second vertex on a line) is
   * the canonical one in deal.II, as
   * described in the general documentation
   * of this class.
   *
   * For <tt>dim=2</tt> this call
   * is simply passed down to the
   * face_to_cell_vertices()
   * function.
   */
  static
  unsigned int
  line_to_cell_vertices (const unsigned int line,
                         const unsigned int vertex);

  /**
   * Map face vertex number to cell
   * vertex number, i.e. give the
   * cell vertex number of the
   * <tt>vertex</tt>th vertex of
   * face <tt>face</tt>, e.g.
   * <tt>GeometryInfo<2>::face_to_cell_vertices(3,0)=2</tt>,
   * see the image under point N4
   * in the 2d section of this
   * class's documentation.
   *
   * Through the <tt>face_orientation</tt>,
   * <tt>face_flip</tt> and
   * <tt>face_rotation</tt> arguments this
   * function handles faces oriented in the
   * standard and non-standard orientation.
   * <tt>face_orientation</tt> defaults to
   * <tt>true</tt>, <tt>face_flip</tt> and
   * <tt>face_rotation</tt> default to
   * <tt>false</tt> (standard orientation).
   * In 2d only <tt>face_flip</tt> is considered.
   * See this @ref GlossFaceOrientation "glossary"
   * article for more information.
   *
   * As the children of a cell are
   * ordered according to the
   * vertices of the cell, this
   * call is passed down to the
   * child_cell_on_face() function.
   * Hence this function is simply
   * a wrapper of
   * child_cell_on_face() giving it
   * a suggestive name.
   */
  static
  unsigned int
  face_to_cell_vertices (const unsigned int face,
                         const unsigned int vertex,
                         const bool face_orientation = true,
                         const bool face_flip        = false,
                         const bool face_rotation    = false);

  /**
   * Map face line number to cell
   * line number, i.e. give the
   * cell line number of the
   * <tt>line</tt>th line of face
   * <tt>face</tt>, e.g.
   * <tt>GeometryInfo<3>::face_to_cell_lines(5,0)=4</tt>.
   *
   * Through the <tt>face_orientation</tt>,
   * <tt>face_flip</tt> and
   * <tt>face_rotation</tt> arguments this
   * function handles faces oriented in the
   * standard and non-standard orientation.
   * <tt>face_orientation</tt> defaults to
   * <tt>true</tt>, <tt>face_flip</tt> and
   * <tt>face_rotation</tt> default to
   * <tt>false</tt> (standard orientation)
   * and has no effect in 2d.
   */
  static
  unsigned int
  face_to_cell_lines (const unsigned int face,
                      const unsigned int line,
                      const bool face_orientation = true,
                      const bool face_flip        = false,
                      const bool face_rotation    = false);

  /**
   * Map the vertex index @p vertex of a face
   * in standard orientation to one of a face
   * with arbitrary @p face_orientation, @p
   * face_flip and @p face_rotation. The
   * values of these three flags default to
   * <tt>true</tt>, <tt>false</tt> and
   * <tt>false</tt>, respectively. this
   * combination describes a face in standard
   * orientation.
   *
   * This function is only implemented in 3D.
   */
  static
  unsigned int
  standard_to_real_face_vertex (const unsigned int vertex,
                                const bool face_orientation = true,
                                const bool face_flip        = false,
                                const bool face_rotation    = false);

  /**
   * Map the vertex index @p vertex of a face
   * with arbitrary @p face_orientation, @p
   * face_flip and @p face_rotation to a face
   * in standard orientation. The values of
   * these three flags default to
   * <tt>true</tt>, <tt>false</tt> and
   * <tt>false</tt>, respectively. this
   * combination describes a face in standard
   * orientation.
   *
   * This function is only implemented in 3D.
   */
  static
  unsigned int
  real_to_standard_face_vertex (const unsigned int vertex,
                                const bool face_orientation = true,
                                const bool face_flip        = false,
                                const bool face_rotation    = false);

  /**
   * Map the line index @p line of a face
   * in standard orientation to one of a face
   * with arbitrary @p face_orientation, @p
   * face_flip and @p face_rotation. The
   * values of these three flags default to
   * <tt>true</tt>, <tt>false</tt> and
   * <tt>false</tt>, respectively. this
   * combination describes a face in standard
   * orientation.
   *
   * This function is only implemented in 3D.
   */
  static
  unsigned int
  standard_to_real_face_line (const unsigned int line,
                              const bool face_orientation = true,
                              const bool face_flip        = false,
                              const bool face_rotation    = false);

  /**
   * Map the line index @p line of a face
   * with arbitrary @p face_orientation, @p
   * face_flip and @p face_rotation to a face
   * in standard orientation. The values of
   * these three flags default to
   * <tt>true</tt>, <tt>false</tt> and
   * <tt>false</tt>, respectively. this
   * combination describes a face in standard
   * orientation.
   *
   * This function is only implemented in 3D.
   */
  static
  unsigned int
  real_to_standard_face_line (const unsigned int line,
                              const bool face_orientation = true,
                              const bool face_flip        = false,
                              const bool face_rotation    = false);

  /**
   * Return the position of the @p ith
   * vertex on the unit cell. The order of
   * vertices is the canonical one in
   * deal.II, as described in the general
   * documentation of this class.
   */
  static
  Point<dim>
  unit_cell_vertex (const unsigned int vertex);

  /**
   * Given a point @p p in unit
   * coordinates, return the number
   * of the child cell in which it
   * would lie in. If the point
   * lies on the interface of two
   * children, return any one of
   * their indices. The result is
   * always less than
   * GeometryInfo<dimension>::max_children_per_cell.
   *
   * The order of child cells is described
   * the general documentation of this
   * class.
   */
  static
  unsigned int
  child_cell_from_point (const Point<dim> &p);

  /**
   * Given coordinates @p p on the
   * unit cell, return the values
   * of the coordinates of this
   * point in the coordinate system
   * of the given child. Neither
   * original nor returned
   * coordinates need actually be
   * inside the cell, we simply
   * perform a scale-and-shift
   * operation with a shift that
   * depends on the number of the
   * child.
   */
  static
  Point<dim>
  cell_to_child_coordinates (const Point<dim>          &p,
                             const unsigned int         child_index,
                             const RefinementCase<dim>  refine_case
                             = RefinementCase<dim>::isotropic_refinement);

  /**
   * The reverse function to the
   * one above: take a point in the
   * coordinate system of the
   * child, and transform it to the
   * coordinate system of the
   * mother cell.
   */
  static
  Point<dim>
  child_to_cell_coordinates (const Point<dim>          &p,
                             const unsigned int         child_index,
                             const RefinementCase<dim>  refine_case
                             = RefinementCase<dim>::isotropic_refinement);

  /**
   * Return true if the given point
   * is inside the unit cell of the
   * present space dimension.
   */
  static
  bool
  is_inside_unit_cell (const Point<dim> &p);

  /**
   * Return true if the given point
   * is inside the unit cell of the
   * present space dimension. This
   * * function accepts an
   * additional * parameter which
   * specifies how * much the point
   * position may * actually be
   * outside the true * unit
   * cell. This is useful because
   * in practice we may often not
   * be able to compute the
   * coordinates of a point in
   * reference coordinates exactly,
   * but only up to numerical
   * roundoff.
   *
   * The tolerance parameter may be
   * less than zero, indicating
   * that the point should be
   * safely inside the cell.
   */
  static
  bool
  is_inside_unit_cell (const Point<dim> &p,
                       const double eps);

  /**
   * Projects a given point onto the
   * unit cell, i.e. each coordinate
   * outside [0..1] is modified
   * to lie within that interval.
   */
  static
  Point<dim>
  project_to_unit_cell (const Point<dim> &p);

  /**
   * Returns the infinity norm of
   * the vector between a given point @p p
   * outside the unit cell to the closest
   * unit cell boundary.
   * For points inside the cell, this is
   * defined as zero.
   */
  static
  double
  distance_to_unit_cell (const Point<dim> &p);

  /**
   * Compute the value of the $i$-th
   * $d$-linear (i.e. (bi-,tri-)linear)
   * shape function at location $\xi$.
   */
  static
  double
  d_linear_shape_function (const Point<dim> &xi,
                           const unsigned int i);

  /**
   * Compute the gradient of the $i$-th
   * $d$-linear (i.e. (bi-,tri-)linear)
   * shape function at location $\xi$.
   */
  static
  Tensor<1,dim>
  d_linear_shape_function_gradient (const Point<dim> &xi,
                                    const unsigned int i);

  /**
   * For a (bi-, tri-)linear
   * mapping from the reference
   * cell, face, or edge to the
   * object specified by the given
   * vertices, compute the
   * alternating form of the
   * transformed unit vectors
   * vertices. For an object of
   * dimensionality @p dim, there
   * are @p dim vectors with @p
   * spacedim components each, and
   * the alternating form is a
   * tensor of rank spacedim-dim
   * that corresponds to the wedge
   * product of the @p dim unit
   * vectors, and it corresponds to
   * the volume and normal vectors
   * of the mapping from reference
   * element to the element
   * described by the vertices.
   *
   * For example, if dim==spacedim==2, then
   * the alternating form is a scalar
   * (because spacedim-dim=0) and its value
   * equals $\mathbf v_1\wedge \mathbf
   * v_2=\mathbf v_1^\perp \cdot\mathbf
   * v_2$, where $\mathbf v_1^\perp$ is a
   * vector that is rotated to the right by
   * 90 degrees from $\mathbf v_1$. If
   * dim==spacedim==3, then the result is
   * again a scalar with value $\mathbf
   * v_1\wedge \mathbf v_2 \wedge \mathbf
   * v_3 = (\mathbf v_1\times \mathbf
   * v_2)\cdot \mathbf v_3$, where $\mathbf
   * v_1, \mathbf v_2, \mathbf v_3$ are the
   * images of the unit vectors at a vertex
   * of the unit dim-dimensional cell under
   * transformation to the dim-dimensional
   * cell in spacedim-dimensional space. In
   * both cases, i.e. for dim==2 or 3, the
   * result happens to equal the
   * determinant of the Jacobian of the
   * mapping from reference cell to cell in
   * real space. Note that it is the actual
   * determinant, not its absolute value as
   * often used in transforming integrals
   * from one coordinate system to
   * another. In particular, if the object
   * specified by the vertices is a
   * parallelogram (i.e. a linear
   * transformation of the reference cell)
   * then the computed values are the same
   * at all vertices and equal the (signed)
   * area of the cell; similarly, for
   * parallel-epipeds, it is the volume of
   * the cell.
   *
   * Likewise, if we have dim==spacedim-1
   * (e.g. we have a quad in 3d space, or a
   * line in 2d), then the alternating
   * product denotes the normal vector
   * (i.e. a rank-1 tensor, since
   * spacedim-dim=1) to the object at each
   * vertex, where the normal vector's
   * magnitude denotes the area element of
   * the transformation from the reference
   * object to the object given by the
   * vertices. In particular, if again the
   * mapping from reference object to the
   * object under consideration here is
   * linear (not bi- or trilinear), then
   * the returned vectors are all
   * %parallel, perpendicular to the mapped
   * object described by the vertices, and
   * have a magnitude equal to the
   * area/volume of the mapped object. If
   * dim=1, spacedim=2, then the returned
   * value is $\mathbf v_1^\perp$, where
   * $\mathbf v_1$ is the image of the sole
   * unit vector of a line mapped to the
   * line in 2d given by the vertices; if
   * dim=2, spacedim=3, then the returned
   * values are $\mathbf v_1 \wedge \mathbf
   * v_2=\mathbf v_1 \times \mathbf v_2$
   * where $\mathbf v_1,\mathbf v_2$ are
   * the two three-dimensional vectors that
   * are tangential to the quad mapped into
   * three-dimensional space.
   *
   * This function is used in order to
   * determine how distorted a cell is (see
   * the entry on
   * @ref GlossDistorted "distorted cells"
   * in the glossary).
   */
  template <int spacedim>
  static void
  alternating_form_at_vertices
#ifndef DEAL_II_CONSTEXPR_BUG
  (const Point<spacedim> (&vertices)[vertices_per_cell],
   Tensor<spacedim-dim,spacedim> (&forms)[vertices_per_cell]);
#else
  (const Point<spacedim> *vertices,
   Tensor<spacedim-dim,spacedim> *forms);
#endif

  /**
   * For each face of the reference
   * cell, this field stores the
   * coordinate direction in which
   * its normal vector points. In
   * <tt>dim</tt> dimension these
   * are the <tt>2*dim</tt> first
   * entries of
   * <tt>{0,0,1,1,2,2,3,3}</tt>.
   *
   * Note that this is only the
   * coordinate number. The actual
   * direction of the normal vector
   * is obtained by multiplying the
   * unit vector in this direction
   * with #unit_normal_orientation.
   */
  static const unsigned int unit_normal_direction[faces_per_cell];

  /**
   * Orientation of the unit normal
   * vector of a face of the
   * reference cell. In
   * <tt>dim</tt> dimension these
   * are the <tt>2*dim</tt> first
   * entries of
   * <tt>{-1,1,-1,1,-1,1,-1,1}</tt>.
   *
   * Each value is either
   * <tt>1</tt> or <tt>-1</tt>,
   * corresponding to a normal
   * vector pointing in the
   * positive or negative
   * coordinate direction,
   * respectively.
   *
   * Note that this is only the
   * <em>standard orientation</em>
   * of faces. At least in 3d,
   * actual faces of cells in a
   * triangulation can also have
   * the opposite orientation,
   * depending on a flag that one
   * can query from the cell it
   * belongs to. For more
   * information, see the
   * @ref GlossFaceOrientation "glossary"
   * entry on
   * face orientation.
   */
  static const int unit_normal_orientation[faces_per_cell];

  /**
   * List of numbers which denotes which
   * face is opposite to a given face. Its
   * entries are the first <tt>2*dim</tt>
   * entries of
   * <tt>{ 1, 0, 3, 2, 5, 4, 7, 6}</tt>.
   */
  static const unsigned int opposite_face[faces_per_cell];


  /**
   * Exception
   */
  DeclException1 (ExcInvalidCoordinate,
                  double,
                  << "The coordinates must satisfy 0 <= x_i <= 1, "
                  << "but here we have x_i=" << arg1);

  /**
   * Exception
   */
  DeclException3 (ExcInvalidSubface,
                  int, int, int,
                  << "RefinementCase<dim> " << arg1 << ": face " << arg2
                  << " has no subface " << arg3);
};




#ifndef DOXYGEN


/* -------------- declaration of explicit specializations ------------- */

#ifndef DEAL_II_MEMBER_ARRAY_SPECIALIZATION_BUG
template <>
const unsigned int GeometryInfo<1>::unit_normal_direction[faces_per_cell];
template <>
const unsigned int GeometryInfo<2>::unit_normal_direction[faces_per_cell];
template <>
const unsigned int GeometryInfo<3>::unit_normal_direction[faces_per_cell];
template <>
const unsigned int GeometryInfo<4>::unit_normal_direction[faces_per_cell];

template <>
const int GeometryInfo<1>::unit_normal_orientation[faces_per_cell];
template <>
const int GeometryInfo<2>::unit_normal_orientation[faces_per_cell];
template <>
const int GeometryInfo<3>::unit_normal_orientation[faces_per_cell];
template <>
const int GeometryInfo<4>::unit_normal_orientation[faces_per_cell];

template <>
const unsigned int GeometryInfo<1>::opposite_face[faces_per_cell];
template <>
const unsigned int GeometryInfo<2>::opposite_face[faces_per_cell];
template <>
const unsigned int GeometryInfo<3>::opposite_face[faces_per_cell];
template <>
const unsigned int GeometryInfo<4>::opposite_face[faces_per_cell];
#endif


template <>
Tensor<1,1>
GeometryInfo<1>::
d_linear_shape_function_gradient (const Point<1> &xi,
                                  const unsigned int i);
template <>
Tensor<1,2>
GeometryInfo<2>::
d_linear_shape_function_gradient (const Point<2> &xi,
                                  const unsigned int i);
template <>
Tensor<1,3>
GeometryInfo<3>::
d_linear_shape_function_gradient (const Point<3> &xi,
                                  const unsigned int i);




/* -------------- inline functions ------------- */

namespace internal
{

  template <int dim>
  inline
  SubfaceCase<dim>::SubfaceCase (const typename SubfacePossibilities<dim>::Possibilities subface_possibility)
    :
    value (subface_possibility)
  {}


  template <int dim>
  inline
  SubfaceCase<dim>::operator unsigned char () const
  {
    return value;
  }


} // namespace internal


template <int dim>
inline
RefinementCase<dim>
RefinementCase<dim>::cut_axis (const unsigned int)
{
  Assert (false, ExcInternalError());
  return static_cast<unsigned char>(-1);
}


template <>
inline
RefinementCase<1>
RefinementCase<1>::cut_axis (const unsigned int i)
{
  const unsigned int dim = 1;
  Assert (i < dim, ExcIndexRange(i, 0, dim));

  static const RefinementCase options[dim] = { RefinementPossibilities<1>::cut_x };
  return options[i];
}



template <>
inline
RefinementCase<2>
RefinementCase<2>::cut_axis (const unsigned int i)
{
  const unsigned int dim = 2;
  Assert (i < dim, ExcIndexRange(i, 0, dim));

  static const RefinementCase options[dim] = { RefinementPossibilities<2>::cut_x,
                                               RefinementPossibilities<2>::cut_y
                                             };
  return options[i];
}



template <>
inline
RefinementCase<3>
RefinementCase<3>::cut_axis (const unsigned int i)
{
  const unsigned int dim = 3;
  Assert (i < dim, ExcIndexRange(i, 0, dim));

  static const RefinementCase options[dim] = { RefinementPossibilities<3>::cut_x,
                                               RefinementPossibilities<3>::cut_y,
                                               RefinementPossibilities<3>::cut_z
                                             };
  return options[i];
}



template <int dim>
inline
RefinementCase<dim>::RefinementCase ()
  :
  value (RefinementPossibilities<dim>::no_refinement)
{}



template <int dim>
inline
RefinementCase<dim>::
RefinementCase (const typename RefinementPossibilities<dim>::Possibilities refinement_case)
  :
  value (refinement_case)
{
  // check that only those bits of
  // the given argument are set that
  // make sense for a given space
  // dimension
  Assert ((refinement_case & RefinementPossibilities<dim>::isotropic_refinement) ==
          refinement_case,
          ExcInvalidRefinementCase (refinement_case));
}



template <int dim>
inline
RefinementCase<dim>::RefinementCase (const unsigned char refinement_case)
  :
  value (refinement_case)
{
  // check that only those bits of
  // the given argument are set that
  // make sense for a given space
  // dimension
  Assert ((refinement_case & RefinementPossibilities<dim>::isotropic_refinement) ==
          refinement_case,
          ExcInvalidRefinementCase (refinement_case));
}



template <int dim>
inline
RefinementCase<dim>::operator unsigned char () const
{
  return value;
}



template <int dim>
inline
RefinementCase<dim>
RefinementCase<dim>::operator | (const RefinementCase<dim> &r) const
{
  return RefinementCase<dim>(static_cast<unsigned char> (value | r.value));
}



template <int dim>
inline
RefinementCase<dim>
RefinementCase<dim>::operator & (const RefinementCase<dim> &r) const
{
  return RefinementCase<dim>(static_cast<unsigned char> (value & r.value));
}



template <int dim>
inline
RefinementCase<dim>
RefinementCase<dim>::operator ~ () const
{
  return RefinementCase<dim>(static_cast<unsigned char> (
                               (~value) & RefinementPossibilities<dim>::isotropic_refinement));
}




template <int dim>
inline
std::size_t
RefinementCase<dim>::memory_consumption ()
{
  return sizeof(RefinementCase<dim>);
}



template <int dim>
template <class Archive>
void RefinementCase<dim>::serialize (Archive &ar,
                                     const unsigned int)
{
  // serialization can't deal with bitfields, so copy from/to a full sized
  // unsigned char
  unsigned char uchar_value = value;
  ar &uchar_value;
  value = uchar_value;
}




template <>
inline
Point<1>
GeometryInfo<1>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
          ExcIndexRange (vertex, 0, vertices_per_cell));

  return Point<1>(static_cast<double>(vertex));
}



template <>
inline
Point<2>
GeometryInfo<2>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
          ExcIndexRange (vertex, 0, vertices_per_cell));

  return Point<2>(vertex%2, vertex/2);
}



template <>
inline
Point<3>
GeometryInfo<3>::unit_cell_vertex (const unsigned int vertex)
{
  Assert (vertex < vertices_per_cell,
          ExcIndexRange (vertex, 0, vertices_per_cell));

  return Point<3>(vertex%2, vertex/2%2, vertex/4);
}



template <int dim>
inline
Point<dim>
GeometryInfo<dim>::unit_cell_vertex (const unsigned int)
{
  Assert(false, ExcNotImplemented());

  return Point<dim> ();
}



template <>
inline
unsigned int
GeometryInfo<1>::child_cell_from_point (const Point<1> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));

  return (p[0] <= 0.5 ? 0 : 1);
}



template <>
inline
unsigned int
GeometryInfo<2>::child_cell_from_point (const Point<2> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert ((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));

  return (p[0] <= 0.5 ?
          (p[1] <= 0.5 ? 0 : 2) :
          (p[1] <= 0.5 ? 1 : 3));
}



template <>
inline
unsigned int
GeometryInfo<3>::child_cell_from_point (const Point<3> &p)
{
  Assert ((p[0] >= 0) && (p[0] <= 1), ExcInvalidCoordinate(p[0]));
  Assert ((p[1] >= 0) && (p[1] <= 1), ExcInvalidCoordinate(p[1]));
  Assert ((p[2] >= 0) && (p[2] <= 1), ExcInvalidCoordinate(p[2]));

  return (p[0] <= 0.5 ?
          (p[1] <= 0.5 ?
           (p[2] <= 0.5 ? 0 : 4) :
           (p[2] <= 0.5 ? 2 : 6)) :
          (p[1] <= 0.5 ?
           (p[2] <= 0.5 ? 1 : 5) :
           (p[2] <= 0.5 ? 3 : 7)));
}


template <int dim>
inline
unsigned int
GeometryInfo<dim>::child_cell_from_point (const Point<dim> &)
{
  Assert(false, ExcNotImplemented());

  return 0;
}



template <>
inline
Point<1>
GeometryInfo<1>::cell_to_child_coordinates (const Point<1>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<1> refine_case)

{
  Assert (child_index < 2,
          ExcIndexRange (child_index, 0, 2));
  Assert (refine_case==RefinementCase<1>::cut_x,
          ExcInternalError());
  (void)refine_case; // removes -Wunused-parameter warning in optimized mode

  return p*2.0-unit_cell_vertex(child_index);
}



template <>
inline
Point<2>
GeometryInfo<2>::cell_to_child_coordinates (const Point<2>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<2> refine_case)

{
  Assert (child_index < GeometryInfo<2>::n_children(refine_case),
          ExcIndexRange (child_index, 0, GeometryInfo<2>::n_children(refine_case)));

  Point<2> point=p;
  switch (refine_case)
    {
    case RefinementCase<2>::cut_x:
      point[0]*=2.0;
      if (child_index==1)
        point[0]-=1.0;
      break;
    case RefinementCase<2>::cut_y:
      point[1]*=2.0;
      if (child_index==1)
        point[1]-=1.0;
      break;
    case RefinementCase<2>::cut_xy:
      point*=2.0;
      point-=unit_cell_vertex(child_index);
      break;
    default:
      Assert(false, ExcInternalError());
    }

  return point;
}



template <>
inline
Point<3>
GeometryInfo<3>::cell_to_child_coordinates (const Point<3>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<3> refine_case)

{
  Assert (child_index < GeometryInfo<3>::n_children(refine_case),
          ExcIndexRange (child_index, 0, GeometryInfo<3>::n_children(refine_case)));

  Point<3> point=p;
  // there might be a cleverer way to do
  // this, but since this function is called
  // in very few places for initialization
  // purposes only, I don't care at the
  // moment
  switch (refine_case)
    {
    case RefinementCase<3>::cut_x:
      point[0]*=2.0;
      if (child_index==1)
        point[0]-=1.0;
      break;
    case RefinementCase<3>::cut_y:
      point[1]*=2.0;
      if (child_index==1)
        point[1]-=1.0;
      break;
    case RefinementCase<3>::cut_z:
      point[2]*=2.0;
      if (child_index==1)
        point[2]-=1.0;
      break;
    case RefinementCase<3>::cut_xy:
      point[0]*=2.0;
      point[1]*=2.0;
      if (child_index%2==1)
        point[0]-=1.0;
      if (child_index/2==1)
        point[1]-=1.0;
      break;
    case RefinementCase<3>::cut_xz:
      // careful, this is slightly
      // different from xy and yz due to
      // different internal numbering of
      // children!
      point[0]*=2.0;
      point[2]*=2.0;
      if (child_index/2==1)
        point[0]-=1.0;
      if (child_index%2==1)
        point[2]-=1.0;
      break;
    case RefinementCase<3>::cut_yz:
      point[1]*=2.0;
      point[2]*=2.0;
      if (child_index%2==1)
        point[1]-=1.0;
      if (child_index/2==1)
        point[2]-=1.0;
      break;
    case RefinementCase<3>::cut_xyz:
      point*=2.0;
      point-=unit_cell_vertex(child_index);
      break;
    default:
      Assert(false, ExcInternalError());
    }

  return point;
}



template <int dim>
inline
Point<dim>
GeometryInfo<dim>::cell_to_child_coordinates (const Point<dim>         &/*p*/,
                                              const unsigned int        /*child_index*/,
                                              const RefinementCase<dim> /*refine_case*/)

{
  AssertThrow (false, ExcNotImplemented());
  return Point<dim>();
}



template <>
inline
Point<1>
GeometryInfo<1>::child_to_cell_coordinates (const Point<1>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<1> refine_case)

{
  Assert (child_index < 2,
          ExcIndexRange (child_index, 0, 2));
  Assert (refine_case==RefinementCase<1>::cut_x,
          ExcInternalError());
  (void)refine_case; // removes -Wunused-parameter warning in optimized mode

  return (p+unit_cell_vertex(child_index))*0.5;
}



template <>
inline
Point<3>
GeometryInfo<3>::child_to_cell_coordinates (const Point<3>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<3> refine_case)

{
  Assert (child_index < GeometryInfo<3>::n_children(refine_case),
          ExcIndexRange (child_index, 0, GeometryInfo<3>::n_children(refine_case)));

  Point<3> point=p;
  // there might be a cleverer way to do
  // this, but since this function is called
  // in very few places for initialization
  // purposes only, I don't care at the
  // moment
  switch (refine_case)
    {
    case RefinementCase<3>::cut_x:
      if (child_index==1)
        point[0]+=1.0;
      point[0]*=0.5;
      break;
    case RefinementCase<3>::cut_y:
      if (child_index==1)
        point[1]+=1.0;
      point[1]*=0.5;
      break;
    case RefinementCase<3>::cut_z:
      if (child_index==1)
        point[2]+=1.0;
      point[2]*=0.5;
      break;
    case RefinementCase<3>::cut_xy:
      if (child_index%2==1)
        point[0]+=1.0;
      if (child_index/2==1)
        point[1]+=1.0;
      point[0]*=0.5;
      point[1]*=0.5;
      break;
    case RefinementCase<3>::cut_xz:
      // careful, this is slightly
      // different from xy and yz due to
      // different internal numbering of
      // children!
      if (child_index/2==1)
        point[0]+=1.0;
      if (child_index%2==1)
        point[2]+=1.0;
      point[0]*=0.5;
      point[2]*=0.5;
      break;
    case RefinementCase<3>::cut_yz:
      if (child_index%2==1)
        point[1]+=1.0;
      if (child_index/2==1)
        point[2]+=1.0;
      point[1]*=0.5;
      point[2]*=0.5;
      break;
    case RefinementCase<3>::cut_xyz:
      point+=unit_cell_vertex(child_index);
      point*=0.5;
      break;
    default:
      Assert(false, ExcInternalError());
    }

  return point;
}



template <>
inline
Point<2>
GeometryInfo<2>::child_to_cell_coordinates (const Point<2>         &p,
                                            const unsigned int      child_index,
                                            const RefinementCase<2> refine_case)
{
  Assert (child_index < GeometryInfo<2>::n_children(refine_case),
          ExcIndexRange (child_index, 0, GeometryInfo<2>::n_children(refine_case)));

  Point<2> point=p;
  switch (refine_case)
    {
    case RefinementCase<2>::cut_x:
      if (child_index==1)
        point[0]+=1.0;
      point[0]*=0.5;
      break;
    case RefinementCase<2>::cut_y:
      if (child_index==1)
        point[1]+=1.0;
      point[1]*=0.5;
      break;
    case RefinementCase<2>::cut_xy:
      point+=unit_cell_vertex(child_index);
      point*=0.5;
      break;
    default:
      Assert(false, ExcInternalError());
    }

  return point;
}



template <int dim>
inline
Point<dim>
GeometryInfo<dim>::child_to_cell_coordinates (const Point<dim>         &/*p*/,
                                              const unsigned int        /*child_index*/,
                                              const RefinementCase<dim> /*refine_case*/)
{
  AssertThrow (false, ExcNotImplemented());
  return Point<dim>();
}



template <>
inline
bool
GeometryInfo<1>::is_inside_unit_cell (const Point<1> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.);
}



template <>
inline
bool
GeometryInfo<2>::is_inside_unit_cell (const Point<2> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) &&
         (p[1] >= 0.) && (p[1] <= 1.);
}



template <>
inline
bool
GeometryInfo<3>::is_inside_unit_cell (const Point<3> &p)
{
  return (p[0] >= 0.) && (p[0] <= 1.) &&
         (p[1] >= 0.) && (p[1] <= 1.) &&
         (p[2] >= 0.) && (p[2] <= 1.);
}

template <>
inline
bool
GeometryInfo<1>::is_inside_unit_cell (const Point<1> &p,
                                      const double eps)
{
  return (p[0] >= -eps) && (p[0] <= 1.+eps);
}



template <>
inline
bool
GeometryInfo<2>::is_inside_unit_cell (const Point<2> &p,
                                      const double eps)
{
  const double l = -eps, u = 1+eps;
  return (p[0] >= l) && (p[0] <= u) &&
         (p[1] >= l) && (p[1] <= u);
}



template <>
inline
bool
GeometryInfo<3>::is_inside_unit_cell (const Point<3> &p,
                                      const double eps)
{
  const double l = -eps, u = 1.0+eps;
  return (p[0] >= l) && (p[0] <= u) &&
         (p[1] >= l) && (p[1] <= u) &&
         (p[2] >= l) && (p[2] <= u);
}


#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
