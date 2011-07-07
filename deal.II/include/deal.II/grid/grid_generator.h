//---------------------------------------------------------------------------
//    $Id$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__grid_generator_h
#define __deal2__grid_generator_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/grid/tria.h>
#include <map>

DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim> class Triangulation;
template <typename number> class Vector;
template <typename number> class SparseMatrix;


/**
 * This class provides a collection of functions for generating basic
 * triangulations. Below, we try to provide some pictures in order to
 * illustrate at least the more complex ones.
 *
 * Some of these functions receive a flag @p colorize. If this is
 * set, parts of the boundary receive different boundary numbers,
 * allowing them to be distinguished by application programs. See the
 * documentation of the functions for details.
 *
 * Additionally this class provides a function
 * (@p laplace_transformation) that smoothly transforms a grid
 * according to given new boundary points. This can be used to
 * transform (simple-shaped) grids to a more complicated ones, like a
 * shell onto a grid of an airfoil, for example.
 *
 * No meshes for the codimension one case are provided at the moment.
 *
 *
 * @ingroup grid
 * @author Wolfgang Bangerth, Ralf Hartmann, Guido Kanschat, Stefan
 * Nauber, Joerg Weimar, Yaqi Wang, Luca Heltai, 1998, 1999, 2000, 2001, 2002,
 * 2003, 2006, 2007, 2008, 2009, 2010, 2011.
 */
class GridGenerator
{
  public:
    				     /**
				      * Initialize the given triangulation
				      * with a hypercube (line in 1D, square
				      * in 2D, etc) consisting of exactly one
				      * cell. The hypercube volume is the
				      * tensor product interval
				      * <i>[left,right]<sup>dim</sup></i> in
				      * the present number of dimensions,
				      * where the limits are given as
				      * arguments. They default to zero and
				      * unity, then producing the unit
				      * hypercube. All boundary indicators are
				      * set to zero ("not colorized") for 2d
				      * and 3d. In 1d the indicators are
				      * colorized, see hyper_rectangle().
				      *
				      * @image html hyper_cubes.png
				      *
				      * See also
				      * subdivided_hyper_cube() for a
				      * coarse mesh consisting of
				      * several cells. See
				      * hyper_rectangle(), if
				      * different lengths in different
				      * ordinate directions are
				      * required.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim, int spacedim>
    static void hyper_cube (Triangulation<dim,spacedim>  &tria,
			    const double        left = 0.,
			    const double        right= 1.);

    				     /**
				      * Same as hyper_cube(), but
				      * with the difference that not
				      * only one cell is created but
				      * each coordinate direction is
				      * subdivided into
				      * @p repetitions cells. Thus,
				      * the number of cells filling
				      * the given volume is
				      * <tt>repetitions<sup>dim</sup></tt>.
				      *
				      * If spacedim=dim+1 the same
				      * mesh as in the case
				      * spacedim=dim is created, but
				      * the vertices have an
				      * additional coordinate =0. So,
				      * if dim=1 one obtains line
				      * along the x axis in the xy
				      * plane, and if dim=3 one
				      * obtains a square in lying in
				      * the xy plane in 3d space.
				      *
				      * @note The triangulation needs
				      * to be void upon calling this
				      * function.
				      */
    template <int dim>
    static void subdivided_hyper_cube (Triangulation<dim>  &tria,
                                       const unsigned int  repetitions,
                                       const double        left = 0.,
                                       const double        right= 1.);

    				     /**
				      * Create a coordinate-parallel
				      * brick from the two
				      * diagonally opposite corner
				      * points @p p1 and @p p2.
				      *
				      * If the @p colorize flag is
				      * set, the
				      * @p boundary_indicators of the
				      * surfaces are assigned, such
				      * that the lower one in
				      * @p x-direction is 0, the
				      * upper one is 1. The indicators
				      * for the surfaces in
				      * @p y-direction are 2 and 3,
				      * the ones for @p z are 4 and
				      * 5. Additionally, material ids
				      * are assigned to the cells
				      * according to the octant their
				      * center is in: being in the right half
				      * plane for any coordinate
				      * direction <i>x<sub>i</sub></i>
				      * adds 2<sup>i</sup>. For
				      * instance, the center point
				      * (1,-1,1) yields a material id 5.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim, int spacedim>
    static void hyper_rectangle (Triangulation<dim,spacedim> &tria,
				 const Point<spacedim>       &p1,
				 const Point<spacedim>       &p2,
				 const bool                  colorize = false);

				     /**
				      * Create a coordinate-parallel
				      * parallelepiped from the two
				      * diagonally opposite corner
				      * points @p p1 and @p p2. In
				      * dimension @p i,
				      * <tt>repetitions[i]</tt> cells are
				      * generated.
				      *
				      * To get cells with an aspect
				      * ratio different from that of
				      * the domain, use different
				      * numbers of subdivisions in
				      * different coordinate
				      * directions. The minimum number
				      * of subdivisions in each
				      * direction is
				      * 1. @p repetitions is a list
				      * of integers denoting the
				      * number of subdivisions in each
				      * coordinate direction.
				      *
				      * If the @p colorize flag is
				      * set, the
				      * @p boundary_indicators of the
				      * surfaces are assigned, such
				      * that the lower one in
				      * @p x-direction is 0, the
				      * upper one is 1. The indicators
				      * for the surfaces in
				      * @p y-direction are 2 and 3,
				      * the ones for @p z are 4 and
				      * 5.  Additionally, material ids
				      * are assigned to the cells
				      * according to the octant their
				      * center is in: being in the right half
				      * plane for any coordinate
				      * direction <i>x<sub>i</sub></i>
				      * adds 2<sup>i</sup>. For
				      * instance, the center point
				      * (1,-1,1) yields a material id 5.
				      *
				      * Note that the @p colorize flag is
				      * ignored in 1d and is assumed to always
				      * be true. That means the boundary
				      * indicator is 0 on the left and 1 on
				      * the right.  See step-15 for details.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      *
				      * @note For an example of the
				      * use of this function see the
				      * step-28
				      * tutorial program.
				      */
    template <int dim>
    static
    void
    subdivided_hyper_rectangle (Triangulation<dim>              &tria,
				const std::vector<unsigned int> &repetitions,
				const Point<dim>                &p1,
				const Point<dim>                &p2,
				const bool                      colorize=false);

				     /**
				      * Like the previous
				      * function. However, here the
				      * second argument does not
				      * denote the number of
				      * subdivisions in each
				      * coordinate direction, but a
				      * sequence of step sizes for
				      * each coordinate direction. The
				      * domain will therefore be
				      * subdivided into
				      * <code>step_sizes[i].size()</code>
				      * cells in coordinate direction
				      * <code>i</code>, with widths
				      * <code>step_sizes[i][j]</code>
				      * for the <code>j</code>th cell.
				      *
				      * This function is therefore the
				      * right one to generate graded
				      * meshes where cells are
				      * concentrated in certain areas,
				      * rather than a uniformly
				      * subdivided mesh as the
				      * previous function generates.
				      *
				      * The step sizes have to add up
				      * to the dimensions of the hyper
				      * rectangle specified by the
				      * points @p p1 and @p p2.
				      */
    template <int dim>
    static
    void
    subdivided_hyper_rectangle(Triangulation<dim>                      &tria,
			       const std::vector<std::vector<double> > &step_sizes,
			       const Point<dim>                        &p_1,
			       const Point<dim>                        &p_2,
			       const bool                              colorize);

				     /**
				      * Like the previous function, but with
				      * the following twist: the @p
				      * material_id argument is a
				      * dim-dimensional array that, for each
				      * cell, indicates which material_id
				      * should be set. In addition, and this
				      * is the major new functionality, if the
				      * material_id of a cell is <tt>(unsigned
				      * char)(-1)</tt>, then that cell is
				      * deleted from the triangulation,
				      * i.e. the domain will have a void
				      * there.
				      */
    template <int dim>
    static
    void
    subdivided_hyper_rectangle (Triangulation<dim>                       &tria,
				const std::vector< std::vector<double> > &spacing,
				const Point<dim>                         &p,
				const Table<dim,unsigned char>           &material_id,
				const bool                               colorize=false);

				     /**
				      * A parallelogram. The first
				      * corner point is the
				      * origin. The <tt>dim</tt>
				      * adjacent points are the
				      * one-dimensional subtensors of
				      * the tensor provided and
				      * additional points will be sums
				      * of these two vectors.
				      * Colorizing is done according
				      * to hyper_rectangle().
				      *
				      * @note This function is
				      * implemented in 2d only.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void
    parallelogram(Triangulation<dim>&  tria,
		  const Tensor<2,dim>& corners,
		  const bool           colorize=false);


				     /**
				      * Hypercube with a layer of
				      * hypercubes around it. The
				      * first two parameters give the
				      * lower and upper bound of the
				      * inner hypercube in all
				      * coordinate directions.
				      * @p thickness marks the size of
				      * the layer cells.
				      *
				      * If the flag colorize is set,
				      * the outer cells get material
				      * id's according to the
				      * following scheme: extending
				      * over the inner cube in
				      * (+/-) x-direction: 1/2. In y-direction
				      * 4/8, in z-direction 16/32. The cells
				      * at corners and edges (3d) get
				      * these values bitwise or'd.
				      *
				      * Presently only available in 2d
				      * and 3d.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void enclosed_hyper_cube (Triangulation<dim> &tria,
	 			     const double      left = 0.,
				     const double      right= 1.,
				     const double      thickness = 1.,
				     const bool        colorize = false);

				     /**
				      * Initialize the given
				      * triangulation with a
				      * hyperball, i.e. a circle or a
				      * ball around <tt>center</tt>
				      * with given <tt>radius</tt>.
				      *
				      * In order to avoid degenerate
				      * cells at the boundaries, the
				      * circle is triangulated by five
				      * cells, the ball by seven
				      * cells. The diameter of the
				      * center cell is chosen so that
				      * the aspect ratio of the
				      * boundary cells after one
				      * refinement is optimized.
				      *
				      * This function is declared to
				      * exist for triangulations of
				      * all space dimensions, but
				      * throws an error if called in
				      * 1d.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void hyper_ball (Triangulation<dim> &tria,
			    const Point<dim>   &center = Point<dim>(),
			    const double      radius = 1.);

				     /**
				      * This class produces a half
				      * hyper-ball around
				      * <tt>center</tt>, which
				      * contains four elements in 2d
				      * and 6 in 3d. The cut plane is
				      * perpendicular to the
				      * <i>x</i>-axis.
				      *
				      * The boundary indicators for the final
				      * triangulation are 0 for the curved boundary and
				      * 1 for the cut plane.
				      *
				      * The appropriate
				      * boundary class is
				      * HalfHyperBallBoundary, or HyperBallBoundary.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void half_hyper_ball (Triangulation<dim> &tria,
				 const Point<dim>   &center = Point<dim>(),
				 const double      radius = 1.);

				     /**
				      * Create a cylinder around the
				      * x-axis.  The cylinder extends
				      * from <tt>x=-half_length</tt> to
				      * <tt>x=+half_length</tt> and its
				      * projection into the
				      * @p yz-plane is a circle of
				      * radius @p radius.
				      *
				      * In two dimensions, the
				      * cylinder is a rectangle from
				      * <tt>x=-half_length</tt> to
				      * <tt>x=+half_length</tt> and
				      * from <tt>y=-radius</tt> to
				      * <tt>y=radius</tt>.
				      *
				      * The boundaries are colored
				      * according to the following
				      * scheme: 0 for the hull of the
				      * cylinder, 1 for the left hand
				      * face and 2 for the right hand
				      * face.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void cylinder (Triangulation<dim> &tria,
			  const double      radius = 1.,
			  const double      half_length = 1.);

                                     /**
                                      * Create a cutted cone around
                                      * the x-axis.  The cone extends
                                      * from <tt>x=-half_length</tt>
                                      * to <tt>x=half_length</tt> and
                                      * its projection into the @p
                                      * yz-plane is a circle of radius
                                      * @p radius_0 at
                                      * <tt>x=-half_length</tt> and a
                                      * circle of radius @p radius_1
                                      * at <tt>x=+half_length</tt>.
                                      * In between the radius is
                                      * linearly decreasing.
                                      *
                                      * In two dimensions, the cone is
                                      * a trapezoid from
                                      * <tt>x=-half_length</tt> to
                                      * <tt>x=+half_length</tt> and
                                      * from <tt>y=-radius_0</tt> to
                                      * <tt>y=radius_0</tt> at
                                      * <tt>x=-half_length</tt> and
                                      * from <tt>y=-radius_1</tt> to
                                      * <tt>y=radius_1</tt> at
                                      * <tt>x=+half_length</tt>.  In
                                      * between the range of
                                      * <tt>y</tt> is linearly
                                      * decreasing.
				      *
                                      * The boundaries are colored
                                      * according to the following
                                      * scheme: 0 for the hull of the
                                      * cone, 1 for the left hand
                                      * face and 2 for the right hand
                                      * face.
                                      *
				      * An example of use can be found in the
				      * documentation of the ConeBoundary
				      * class, with which you probably want to
				      * associate boundary indicator 0 (the
				      * hull of the cone).
				      *
                                      * @note The triangulation needs to be
                                      * void upon calling this
                                      * function.
				      *
				      * @author Markus B&uuml;rg, 2009
                                      */
    template <int dim>
    static void
    truncated_cone (Triangulation<dim> &tria,
		    const double radius_0 = 1.0,
		    const double radius_1 = 0.5,
		    const double half_length = 1.0);

				     /**
				      * Initialize the given
				      * triangulation with a hyper-L
				      * consisting of exactly
				      * <tt>2^dim-1</tt> cells. It
				      * produces the hypercube with
				      * the interval [<i>left,right</i>] without
				      * the hypercube made out of the
				      * interval [<i>(a+b)/2,b</i>].
				      *
				      * @image html hyper_l.png
				      *
				      * The triangulation needs to be
				      * void upon calling this
				      * function.
				      *
				      * This function is declared to
				      * exist for triangulations of
				      * all space dimensions, but
				      * throws an error if called in
				      * 1d.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void hyper_L (Triangulation<dim> &tria,
			 const double      left = -1.,
			 const double      right= 1.);

                                     /**
				      * Initialize the given
				      * Triangulation with a hypercube
				      * with a slit. In each
				      * coordinate direction, the
				      * hypercube extends from @p left
				      * to @p right.
				      *
				      * In 2d, the split goes in
				      * vertical direction from
				      * <tt>x=(left+right)/2,
				      * y=left</tt> to the center of
				      * the square at
				      * <tt>x=y=(left+right)/2</tt>.
				      *
				      * In 3d, the 2d domain is just
				      * extended in the
				      * <i>z</i>-direction, such that
				      * a plane cuts the lower half of
				      * a rectangle in two.

				      * This function is declared to
				      * exist for triangulations of
				      * all space dimensions, but
				      * throws an error if called in
				      * 1d.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void hyper_cube_slit (Triangulation<dim> &tria,
				 const double      left = 0.,
				 const double      right= 1.,
				 const bool colorize = false);

				     /**
				      * Produce a hyper-shell,
				      * the region between two
				      * spheres around <tt>center</tt>,
				      * with given
				      * <tt>inner_radius</tt> and
				      * <tt>outer_radius</tt>.
				      *
				      * If the flag @p colorize is @p
				      * true, then the outer boundary
				      * will have the id 1, while the
				      * inner boundary has id zero. If
				      * the flag is @p false, both
				      * have id zero.
				      *
				      * In 2D, the number
				      * <tt>n_cells</tt> of elements
				      * for this initial triangulation
				      * can be chosen arbitrarily. If
				      * the number of initial cells is
				      * zero (as is the default), then
				      * it is computed adaptively such
				      * that the resulting elements
				      * have the least aspect ratio.
				      *
				      * In 3D, only two different numbers are
				      * meaningful, 6 for a surface based on a
				      * hexahedron (i.e. 6 panels on the inner
				      * sphere extruded in radial direction to
				      * form 6 cells) and 12 for the rhombic
				      * dodecahedron. These give rise to the
				      * following meshes upon one refinement:
				      *
				      * @image html hypershell3d-6.png
				      * @image html hypershell3d-12.png
				      *
				      * Neither of these meshes is
				      * particularly good since one ends up
				      * with poorly shaped cells at the inner
				      * edge upon refinement. For example,
				      * this is the middle plane of the mesh
				      * for the <code>n_cells=6</code>:
				      *
				      * @image html hyper_shell_6_cross_plane.png
				      *
				      * The mesh generated with
				      * <code>n_cells=6</code> is better but
				      * still not good. As a consequence, you
				      * may also specify
				      * <code>n_cells=96</code> as a third
				      * option. The mesh generated in this way
				      * is based on a once refined version of
				      * the one with <code>n_cells=12</code>,
				      * where all internal nodes are re-placed
				      * along a shell somewhere between the
				      * inner and outer boundary of the
				      * domain. The following two images
				      * compare half of the hyper shell for
				      * <code>n_cells=12</code> and
				      * <code>n_cells=96</code> (note that the
				      * doubled radial lines on the cross
				      * section are artifacts of the
				      * visualization):
				      *
				      * @image html hyper_shell_12_cut.png
				      * @image html hyper_shell_96_cut.png
				      *
				      * @note This function is declared to
				      * exist for triangulations of
				      * all space dimensions, but
				      * throws an error if called in
				      * 1d.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void hyper_shell (Triangulation<dim>   &tria,
			     const Point<dim>     &center,
			     const double        inner_radius,
			     const double        outer_radius,
			     const unsigned int  n_cells = 0,
			     bool colorize = false);

				     /**
				      * Produce a half hyper-shell,
				      * i.e. the space between two
				      * circles in two space
				      * dimensions and the region
				      * between two spheres in 3d,
				      * with given inner and outer
				      * radius and a given number of
				      * elements for this initial
				      * triangulation.  However,
				      * opposed to the previous
				      * function, it does not produce
				      * a whole shell, but only one
				      * half of it, namely that part
				      * for which the first component
				      * is restricted to non-negative
				      * values. The purpose of this
				      * class is to enable
				      * computations for solutions
				      * which have rotational
				      * symmetry, in which case the
				      * half shell in 2d represents a
				      * shell in 3d.
				      *
				      * If the number of
				      * initial cells is zero (as is
				      * the default), then it is
				      * computed adaptively such that
				      * the resulting elements have
				      * the least aspect ratio.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function.
				      */
    template <int dim>
    static void half_hyper_shell (Triangulation<dim>   &tria,
				  const Point<dim>     &center,
				  const double        inner_radius,
				  const double        outer_radius,
				  const unsigned int  n_cells = 0);

    
				     /**
				      * Produce a quarter hyper-shell,
				      * i.e. the space between two
				      * circles in two space
				      * dimensions and the region
				      * between two spheres in 3d,
				      * with given inner and outer
				      * radius and a given number of
				      * elements for this initial
				      * triangulation. All components are
				      * restricted to be positive, so the
				      * opening angle is 90 degrees.
				      *
				      * If the number of
				      * initial cells is zero (as is
				      * the default), then it is
				      * computed adaptively such that
				      * the resulting elements have
				      * the least aspect ratio.
				      *
				      * @note The triangulation needs to be
				      * void upon calling this
				      * function. Only implemented in 2d so far.
				      */
    template <int dim>
    static void quarter_hyper_shell (Triangulation<dim>   &tria,
				  const Point<dim>     &center,
				  const double        inner_radius,
				  const double        outer_radius,
				  const unsigned int  n_cells = 0,
				  const bool colorize = false);
   
				     /**
				      * Produce a domain that is the space
				      * between two cylinders in 3d, with
				      * given length, inner and outer radius
				      * and a given number of elements for
				      * this initial triangulation. If @p
				      * n_radial_cells is zero (as is the
				      * default), then it is computed
				      * adaptively such that the resulting
				      * elements have the least aspect
				      * ratio. The same holds for @p
				      * n_axial_cells.
				      *
				      * @note Although this function
				      * is declared as a template, it
				      * does not make sense in 1D and
				      * 2D.
				      *
				      * @note The triangulation needs
				      * to be void upon calling this
				      * function.
				      */
    template <int dim>
    static void cylinder_shell (Triangulation<dim>   &tria,
                                const double        length,
                                const double        inner_radius,
                                const double        outer_radius,
                                const unsigned int  n_radial_cells = 0,
                                const unsigned int  n_axial_cells = 0);

				     /**
				      * This class produces a square
				      * on the <i>xy</i>-plane with a
				      * circular hole in the middle,
				      * times the interval [0.L]
				      * (only in 3d).
				      *
				      *	 @image html cubes_hole.png
				      *
				      * It is implemented in 2d and
				      * 3d, and takes the following
				      * arguments:
				      *
				      * @arg @p inner_radius: size of the
                                      *    internal hole
				      * @arg @p  outer_radius: size of the
                                      *    biggest enclosed cylinder
				      * @arg @p L: extension on the @p z-direction
				      * @arg @p repetitions: number of subdivisions
				      *      along the @p z-direction
				      * @arg @p colorize: wether to assign different
				      *     boundary indicators to different faces.
				      *    The colors are given in lexicographic
				      *    ordering for the flat faces (0 to 3 in 2d,
				      *    0 to 5 in 3d) plus the curved hole
				      *    (4 in 2d, and 6 in 3d).
				      *    If @p colorize is set to false, then flat faces
				      *    get the number 0 and the hole gets number 1.
				      */
    template<int dim>
    static void hyper_cube_with_cylindrical_hole (Triangulation<dim> &triangulation,
						const double inner_radius = .25,
						const double outer_radius = .5,
						const double L = .5,
						const unsigned int repetition = 1,
						const bool colorize = false);

				     /**
				      * Produce a ring of cells in 3D that is
				      * cut open, twisted and glued together
				      * again. This results in a kind of
				      * moebius-loop.
				      *
				      * @param tria        The triangulation to be worked on.
				      * @param n_cells     The number of cells in the loop. Must be greater than 4.
				      * @param n_rotations The number of rotations (Pi/2 each) to be performed before glueing the loop together.
				      * @param R           The radius of the circle, which forms the middle line of the torus containing the loop of cells. Must be greater than @p r.
				      * @param r           The radius of the cylinder bend together as loop.
				      */
    static void moebius (Triangulation<3,3>&  tria,
			 const unsigned int   n_cells,
			 const unsigned int   n_rotations,
			 const double         R,
			 const double         r);

				     /**
				      * Given the two triangulations
				      * specified as the first two
				      * arguments, create the
				      * triangulation that contains
				      * the cells of both
				      * triangulation and store it in
				      * the third parameter. Previous
				      * content of @p result will be
				      * deleted.
				      * 
				      * This function is most often used 
				      * to compose meshes for more
				      * complicated geometries if the
				      * geometry can be composed of
				      * simpler parts for which functions
				      * exist to generate coarse meshes.
				      * For example, the channel mesh used
				      * in step-35 could in principle be
				      * created using a mesh created by the 
				      * GridGenerator::hyper_cube_with_cylindrical_hole
				      * function and several rectangles,
				      * and merging them using the current
				      * function. The rectangles will
				      * have to be translated to the
				      * right for this, a task that can
				      * be done using the GridTools::shift
				      * function (other tools to transform
				      * individual mesh building blocks are
				      * GridTools::transform, GridTools::rotate,
				      * and GridTools::scale).
				      * 
				      * @note The two input triangulations
				      * must be coarse meshes that have
				      * no refined cells. 
				      * 
				      * @note The function copies the material ids
				      * of the cells of the two input
				      * triangulations into the output
				      * triangulation but it currently makes
				      * no attempt to do the same for boundary
				      * ids. In other words, if the two
				      * coarse meshes have anything but
				      * the default boundary indicators,
				      * then you will currently have to set
				      * boundary indicators again by hand
				      * in the output triangulation.
				      * 
				      * @note For a related operation
				      * on refined meshes when both
				      * meshes are derived from the
				      * same coarse mesh, see
				      * GridTools::create_union_triangulation .
				      */
    template <int dim, int spacedim>
    static
    void
    merge_triangulations (const Triangulation<dim, spacedim> &triangulation_1,
			  const Triangulation<dim, spacedim> &triangulation_2,
			  Triangulation<dim, spacedim>       &result);

                                     /**
				      * This function transformes the
				      * @p Triangulation @p tria
				      * smoothly to a domain that is
				      * described by the boundary
				      * points in the map
				      * @p new_points. This map maps
				      * the point indices to the
				      * boundary points in the
				      * transformed domain.
				      *
				      * Note, that the
				      * @p Triangulation is changed
				      * in-place, therefore you don't
				      * need to keep two
				      * triangulations, but the given
				      * triangulation is changed
				      * (overwritten).
				      *
				      * In 1d, this function is not
				      * currently implemented.
				      */
    template <int dim>
    static void laplace_transformation (Triangulation<dim> &tria,
					const std::map<unsigned int,Point<dim> > &new_points);

				     /**
				      * Exception
				      */
    DeclException0 (ExcInvalidRadii);
                                     /**
                                      * Exception
                                      */
    DeclException1 (ExcInvalidRepetitions,
                    int,
                    << "The number of repetitions " << arg1
                    << " must be >=1.");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidRepetitionsDimension,
                    int,
                    << "The vector of repetitions  must have "
		    << arg1 <<" elements.");

  private:
				     /**
				      * Perform the action specified
				      * by the @p colorize flag of
				      * the hyper_rectangle()
				      * function of this class.
				      */
    template <int dim, int spacedim>
    static
    void
    colorize_hyper_rectangle (Triangulation<dim,spacedim> &tria);

				     /**
				      * Perform the action specified
				      * by the @p colorize flag of
				      * the
				      * subdivided_hyper_rectangle()
				      * function of this class. This
				      * function is singled out
				      * because it is dimension
				      * specific.
				      */
    template <int dim>
    static
    void
    colorize_subdivided_hyper_rectangle (Triangulation<dim> &tria,
					 const Point<dim>   &p1,
					 const Point<dim>   &p2,
					 const double        epsilon);

				     /**
				      * Assign boundary number zero to
				      * the inner shell boundary and 1
				      * to the outer.
				      */
    template<int dim>
    static
    void
    colorize_hyper_shell (Triangulation<dim>& tria,
			  const Point<dim>& center,
			  const double inner_radius,
			  const double outer_radius);

				     /**
				      * Solve the Laplace equation for
				      * @p laplace_transformation
				      * function for one of the
				      * @p dim space
				      * dimensions. Externalized into
				      * a function of its own in order
				      * to allow parallel execution.
				      */
    static
    void
    laplace_solve (const SparseMatrix<double>          &S,
		   const std::map<unsigned int,double> &m,
		   Vector<double>                      &u);
};



DEAL_II_NAMESPACE_CLOSE

#endif
