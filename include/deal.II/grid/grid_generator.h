// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__grid_generator_h
#define dealii__grid_generator_h


#include <deal.II/base/config.h>
#include <deal.II/base/std_cxx11/array.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>
#include <deal.II/base/function.h>
#include <deal.II/grid/tria.h>
#include <map>

DEAL_II_NAMESPACE_OPEN

/**
 * This namespace provides a collection of functions for generating
 * triangulations for some basic geometries.
 *
 * Some of these functions receive a flag @p colorize. If this is set, parts
 * of the boundary receive different
 * @ref GlossBoundaryIndicator "boundary indicators"),
 * allowing them to be distinguished for the purpose of attaching geometry
 * objects and evaluating different boundary conditions.
 *
 * @ingroup grid
 */
namespace GridGenerator
{
  /**
   * @name Creating meshes for basic geometries
   */
  ///@{

  /**
   * Initialize the given triangulation with a hypercube (line in 1D, square
   * in 2D, etc) consisting of exactly one cell. The hypercube volume is the
   * tensor product interval $[left,right]^{\text{dim}}$ in the present number
   * of dimensions, where the limits are given as arguments. They default to
   * zero and unity, then producing the unit hypercube. If the argument @p
   * colorize is false, all boundary indicators are set to zero ("not
   * colorized") for 2d and 3d. If it is true, the boundary is colorized as in
   * hyper_rectangle(). In 1d the indicators are always colorized, see
   * hyper_rectangle().
   *
   * @image html hyper_cubes.png
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero. For example, a
   * <tt>Triangulation@<2,3@></tt> will be a square in the xy plane with z=0.
   *
   * See also subdivided_hyper_cube() for a coarse mesh consisting of several
   * cells. See hyper_rectangle(), if different lengths in different ordinate
   * directions are required.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim, int spacedim>
  void hyper_cube (Triangulation<dim,spacedim>  &tria,
                   const double                  left = 0.,
                   const double                  right= 1.,
                   const bool                    colorize= false);

  /**
   * \brief %Triangulation of a d-simplex with (d+1) vertices and mesh cells.
   *
   * The @p vertices argument contains a vector with all d+1 vertices of the
   * simplex. They must be given in an order such that the vectors from the
   * first vertex to each of the others form a right-handed system. And I am
   * not happy about the discrimination involved here.
   *
   * The meshes generated in two and three dimensions are
   *
   * @image html simplex_2d.png
   * @image html simplex_3d.png
   *
   * @param tria The Triangulation to create. It needs to be empty upon
   * calling this function.
   *
   * @param vertices The dim+1 corners of the simplex.
   *
   * @note Implemented for <tt>Triangulation@<2,2@></tt>,
   * <tt>Triangulation@<3,3@></tt>.
   *
   * @author Guido Kanschat
   * @date 2015
   */
  template <int dim>
  void simplex(Triangulation<dim, dim> &tria,
               const std::vector<Point<dim> > &vertices);

  /**
   * Same as hyper_cube(), but with the difference that not only one cell is
   * created but each coordinate direction is subdivided into @p repetitions
   * cells. Thus, the number of cells filling the given volume is
   * <tt>repetitions<sup>dim</sup></tt>.
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero. For example, a
   * <tt>Triangulation@<2,3@></tt> will be a square in the xy plane with z=0.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim, int spacedim>
  void subdivided_hyper_cube (Triangulation<dim,spacedim>  &tria,
                              const unsigned int   repetitions,
                              const double         left = 0.,
                              const double         right= 1.);

  /**
   * Create a coordinate-parallel brick from the two diagonally opposite
   * corner points @p p1 and @p p2.
   *
   * If the @p colorize flag is set, the @p boundary_ids of the surfaces are
   * assigned, such that the lower one in @p x-direction is 0, the upper one
   * is 1. The indicators for the surfaces in @p y-direction are 2 and 3, the
   * ones for @p z are 4 and 5. Additionally, material ids are assigned to the
   * cells according to the octant their center is in: being in the right half
   * plane for any coordinate direction <i>x<sub>i</sub></i> adds
   * 2<sup>i</sup>. For instance, the center point (1,-1,1) yields a material
   * id 5.
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero. For example, a
   * <tt>Triangulation@<2,3@></tt> will be a rectangle in the xy plane with
   * z=0, defined by the two opposing corners @p p1 and @p p2.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim, int spacedim>
  void hyper_rectangle (Triangulation<dim,spacedim> &tria,
                        const Point<dim>            &p1,
                        const Point<dim>            &p2,
                        const bool                  colorize = false);

  /**
   * Create a coordinate-parallel parallelepiped from the two diagonally
   * opposite corner points @p p1 and @p p2. In direction @p i,
   * <tt>repetitions[i]</tt> cells are generated.
   *
   * To get cells with an aspect ratio different from that of the domain, use
   * different numbers of subdivisions in different coordinate directions. The
   * minimum number of subdivisions in each direction is 1. @p repetitions is
   * a list of integers denoting the number of subdivisions in each coordinate
   * direction.
   *
   * If the @p colorize flag is set, the @p boundary_ids of the surfaces are
   * assigned, such that the lower one in @p x-direction is 0, the upper one
   * is 1 (the left and the right vertical face). The indicators for the
   * surfaces in @p y-direction are 2 and 3, the ones for @p z are 4 and 5.
   * Additionally, material ids are assigned to the cells according to the
   * octant their center is in: being in the right half plane for any
   * coordinate direction <i>x<sub>i</sub></i> adds 2<sup>i</sup>. For
   * instance, the center point (1,-1,1) yields a material id 5 (this means
   * that in 2d only material ids 0,1,2,3 are assigned independent from the
   * number of repetitions).
   *
   * Note that the @p colorize flag is ignored in 1d and is assumed to always
   * be true. That means the boundary indicator is 0 on the left and 1 on the
   * right.  See step-15 for details.
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero. For example, a
   * <tt>Triangulation@<2,3@></tt> will be a rectangle in the xy plane with
   * z=0, defined by the two opposing corners @p p1 and @p p2.
   *
   * @note For an example of the use of this function see the step-28 tutorial
   * program.
   *
   * @param tria The Triangulation to create. It needs to be empty upon
   * calling this function.
   *
   * @param repetitions A vector of dim positive values denoting the number of
   * cells to generate in that direction.
   *
   * @param p1 First corner point.
   *
   * @param p2 Second corner opposite to @p p1.
   *
   * @param colorize Assign different boundary ids if set to true.
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle (Triangulation<dim,spacedim>     &tria,
                              const std::vector<unsigned int> &repetitions,
                              const Point<dim>                &p1,
                              const Point<dim>                &p2,
                              const bool                      colorize=false);

  /**
   * Like the previous function. However, here the second argument does not
   * denote the number of subdivisions in each coordinate direction, but a
   * sequence of step sizes for each coordinate direction. The domain will
   * therefore be subdivided into <code>step_sizes[i].size()</code> cells in
   * coordinate direction <code>i</code>, with widths
   * <code>step_sizes[i][j]</code> for the <code>j</code>th cell.
   *
   * This function is therefore the right one to generate graded meshes where
   * cells are concentrated in certain areas, rather than a uniformly
   * subdivided mesh as the previous function generates.
   *
   * The step sizes have to add up to the dimensions of the hyper rectangle
   * specified by the points @p p1 and @p p2.
   */
  template <int dim>
  void
  subdivided_hyper_rectangle (Triangulation<dim>                      &tria,
                              const std::vector<std::vector<double> > &step_sizes,
                              const Point<dim>                        &p_1,
                              const Point<dim>                        &p_2,
                              const bool                              colorize);

  /**
   * Like the previous function, but with the following twist: the @p
   * material_id argument is a dim-dimensional array that, for each cell,
   * indicates which material_id should be set. In addition, and this is the
   * major new functionality, if the material_id of a cell is <tt>(unsigned
   * char)(-1)</tt>, then that cell is deleted from the triangulation, i.e.
   * the domain will have a void there.
   *
   * @note If you need a lot of holes, you may consider cheese().
   */
  template <int dim>
  void
  subdivided_hyper_rectangle (Triangulation<dim>                       &tria,
                              const std::vector< std::vector<double> > &spacing,
                              const Point<dim>                         &p,
                              const Table<dim,types::material_id>      &material_id,
                              const bool                                colorize=false);

  /**
   * \brief Rectangular domain with rectangular pattern of holes
   *
   * The domain itself is rectangular, very much as if it had been generated
   * by subdivided_hyper_rectangle(). The argument <code>holes</code>
   * specifies how many square holes the domain should have in each coordinate
   * direction. The total number of mesh cells in that direction is then twice
   * this number plus one.
   *
   * The number of holes in one direction must be at least one.
   *
   * An example with two by three holes is
   *
   * @image html cheese_2d.png
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero.
   *
   * @param tria The Triangulation to create. It needs to be empty upon
   * calling this function.
   *
   * @param holes Positive number of holes in each of the dim directions.
   * @author Guido Kanschat
   * @date 2015
   */
  template <int dim, int spacedim>
  void
  cheese (Triangulation<dim, spacedim> &tria,
          const std::vector<unsigned int> &holes);

  /**
   * A parallelogram. The first corner point is the origin. The @p dim
   * adjacent points are the ones given in the second argument and the fourth
   * point will be the sum of these two vectors.  Colorizing is done in the
   * same way as in hyper_rectangle().
   *
   * @note This function is implemented in 2d only.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void
  parallelogram (Triangulation<dim> &tria,
                 const Point<dim>  (&corners)[dim],
                 const bool          colorize=false);

  /**
   * A parallelepiped. The first corner point is the origin. The @p dim
   * adjacent points are vectors describing the edges of the parallelepiped
   * with respect to the origin. Additional points are sums of these dim
   * vectors. Colorizing is done according to hyper_rectangle().
   *
   * @note This function silently reorders the vertices on the cells to
   * lexicographic ordering (see <code>GridReordering::reorder_grid</code>).
   * In other words, if reordering of the vertices does occur, the ordering of
   * vertices in the array of <code>corners</code> will no longer refer to the
   * same triangulation.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void
  parallelepiped (Triangulation<dim> &tria,
                  const Point<dim>  (&corners) [dim],
                  const bool          colorize = false);

  /**
   * A subdivided parallelepiped. The first corner point is the origin. The @p
   * dim adjacent points are vectors describing the edges of the
   * parallelepiped with respect to the origin. Additional points are sums of
   * these dim vectors. The variable @p n_subdivisions designates the number
   * of subdivisions in each of the @p dim directions. Colorizing is done
   * according to hyper_rectangle().
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void
  subdivided_parallelepiped (Triangulation<dim>  &tria,
                             const unsigned int   n_subdivisions,
                             const Point<dim>   (&corners) [dim],
                             const bool           colorize = false);

  /**
   * A subdivided parallelepiped, i.e., the same as above, but where the
   * number of subdivisions in each of the @p dim directions may vary.
   * Colorizing is done according to hyper_rectangle().
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void
  subdivided_parallelepiped (Triangulation<dim>  &tria,
#ifndef _MSC_VER
                             const unsigned int(&n_subdivisions)[dim],
#else
                             const unsigned int *n_subdivisions,
#endif
                             const Point<dim>   (&corners) [dim],
                             const bool           colorize = false);

  /**
   * A subdivided parallelepiped.
   *
   * @param tria The Triangulation to create. It needs to be empty upon
   * calling this function.
   *
   * @param origin First corner of the parallelepiped.
   *
   * @param edges An array of @p dim tensors describing the length and
   * direction of the edges from @p origin.
   *
   * @param subdivisions Number of subdivisions in each of the dim directions.
   * Each entry must be positive. An empty vector is equivalent to one
   * subdivision in each direction.
   *
   * @param colorize Assign different boundary ids if set to true.
   *
   * @note Implemented for all combinations of @p dim and @p spacedim.
   *
   * @note You likely need to help the compiler by explicitly specifying the
   * two template parameters when calling this function.
   */
  template <int dim, int spacedim>
  void
  subdivided_parallelepiped (Triangulation<dim, spacedim>  &tria,
                             const Point<spacedim> &origin,
                             const std_cxx11::array<Tensor<1,spacedim>,dim> &edges,
                             const std::vector<unsigned int> &subdivisions = std::vector<unsigned int>(),
                             const bool colorize = false);

  /**
   * Hypercube with a layer of hypercubes around it. The first two parameters
   * give the lower and upper bound of the inner hypercube in all coordinate
   * directions.  @p thickness marks the size of the layer cells.
   *
   * If the flag @p colorize is set, the outer cells get material id's
   * according to the following scheme: extending over the inner cube in (+/-)
   * x-direction: 1/2. In y-direction 4/8, in z-direction 16/32. The cells at
   * corners and edges (3d) get these values bitwise or'd.
   *
   * Presently only available in 2d and 3d.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void enclosed_hyper_cube (Triangulation<dim> &tria,
                            const double        left = 0.,
                            const double        right= 1.,
                            const double        thickness = 1.,
                            const bool          colorize = false);

  /**
   * Initialize the given triangulation with a hyperball, i.e. a circle or a
   * ball around @p center with given @p radius.
   *
   * In order to avoid degenerate cells at the boundaries, the circle is
   * triangulated by five cells, the ball by seven cells. The diameter of the
   * center cell is chosen so that the aspect ratio of the boundary cells
   * after one refinement is optimized.
   *
   * This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * You should attach a SphericalManifold to the cells and faces for correct
   * placement of vertices upon refinement and to be able to use higher order
   * mappings.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void hyper_ball (Triangulation<dim> &tria,
                   const Point<dim>   &center = Point<dim>(),
                   const double        radius = 1.);

  /**
   * Creates a hyper sphere, i.e., a surface of a ball in @p spacedim
   * dimensions. This function only exists for dim+1=spacedim in 2 and 3 space
   * dimensions.
   *
   * You should attach a SphericalManifold to the cells and faces for correct
   * placement of vertices upon refinement and to be able to use higher order
   * mappings.
   *
   * The following pictures are generated with:
   * @code
   * Triangulation<2,3>   triangulation;
   *
   * static SphericalManifold<2,3> surface_description;
   *
   * GridGenerator::hyper_sphere(triangulation);
   *
   * triangulation.set_all_manifold_ids(0);
   * triangulation.set_manifold (0, surface_description);
   * triangulation.refine_global(3);
   * @endcode
   *
   * See the
   * @ref manifold "documentation module on manifolds"
   * for more details.
   *
   * @image html sphere.png
   * @image html sphere_section.png
   *
   * @note The triangulation needs to be void upon calling this function.
   */

  template <int dim, int spacedim>
  void hyper_sphere (Triangulation<dim,spacedim> &tria,
                     const Point<spacedim>   &center = Point<spacedim>(),
                     const double        radius = 1.);

  /**
   * This class produces a half hyper-ball around @p center, which contains
   * four elements in 2d and 6 in 3d. The cut plane is perpendicular to the
   * <i>x</i>-axis.
   *
   * The boundary indicators for the final triangulation are 0 for the curved
   * boundary and 1 for the cut plane.
   *
   * The appropriate boundary class is HalfHyperBallBoundary, or
   * HyperBallBoundary.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void half_hyper_ball (Triangulation<dim> &tria,
                        const Point<dim>   &center = Point<dim>(),
                        const double        radius = 1.);

  /**
   * Create a cylinder around the $x$-axis.  The cylinder extends from
   * <tt>x=-half_length</tt> to <tt>x=+half_length</tt> and its projection
   * into the @p yz-plane is a circle of radius @p radius.
   *
   * In two dimensions, the cylinder is a rectangle from
   * <tt>x=-half_length</tt> to <tt>x=+half_length</tt> and from
   * <tt>y=-radius</tt> to <tt>y=radius</tt>.
   *
   * The boundaries are colored according to the following scheme: 0 for the
   * hull of the cylinder, 1 for the left hand face and 2 for the right hand
   * face.
   *
   * If you want the cylinder to revolve around a different axis than the
   * $x$-axis, then simply rotate the mesh generated by this function using
   * the GridTools::transform() function using a rotation operator as
   * argument.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void cylinder (Triangulation<dim> &tria,
                 const double        radius = 1.,
                 const double        half_length = 1.);

  /**
   * Create a cut cone around the x-axis.  The cone extends from
   * <tt>x=-half_length</tt> to <tt>x=half_length</tt> and its projection into
   * the @p yz-plane is a circle of radius @p radius_0 at
   * <tt>x=-half_length</tt> and a circle of radius @p radius_1 at
   * <tt>x=+half_length</tt>.  In between the radius is linearly decreasing.
   *
   * In two dimensions, the cone is a trapezoid from <tt>x=-half_length</tt>
   * to <tt>x=+half_length</tt> and from <tt>y=-radius_0</tt> to
   * <tt>y=radius_0</tt> at <tt>x=-half_length</tt> and from
   * <tt>y=-radius_1</tt> to <tt>y=radius_1</tt> at <tt>x=+half_length</tt>.
   * In between the range of <tt>y</tt> is linearly decreasing.
   *
   * The boundaries are colored according to the following scheme: 0 for the
   * hull of the cone, 1 for the left hand face and 2 for the right hand face.
   *
   * An example of use can be found in the documentation of the ConeBoundary
   * class, with which you probably want to associate boundary indicator 0
   * (the hull of the cone).
   *
   * @note The triangulation needs to be void upon calling this function.
   *
   * @author Markus B&uuml;rg, 2009
   */
  template <int dim>
  void
  truncated_cone (Triangulation<dim> &tria,
                  const double        radius_0 = 1.0,
                  const double        radius_1 = 0.5,
                  const double        half_length = 1.0);

  /**
   * \brief A center cell with stacks of cell protruding from each surface.
   *
   * Each of the square mesh cells is Cartesian and has size one in each
   * coordinate direction. The center of cell number zero is the origin.
   *
   * @param tria A Triangulation object which has to be empty.
   *
   * @param sizes A vector of integers of dimension
   * GeometryInfo<dim>::faces_per_cell with the following meaning: the legs of
   * the cross are stacked on the faces of the center cell, in the usual order
   * of deal.II cells, namely first $-x$, then $x$, then $-y$ and so on. The
   * corresponding entries in <code>sizes</code> name the number of cells
   * stacked on this face. All numbers may be zero, thus L- and T-shaped
   * domains are specializations of this domain.
   *
   * @param colorize_cells If colorization is chosen, then the material id of
   * a cells corresponds to the leg it is in. The id of the center cell is
   * zero, and then the legs are numbered starting at one.
   *
   * Examples in two and three dimensions are
   *
   * @image html hyper_cross_2d.png
   * @image html hyper_cross_3d.png
   *
   * @author Guido Kanschat
   * @date 2015
   */
  template <int dim, int spacedim>
  void hyper_cross(Triangulation<dim, spacedim> &tria,
                   const std::vector<unsigned int> &sizes,
                   const bool colorize_cells = false);

  /**
   * Initialize the given triangulation with a hyper-L (in 2d or 3d)
   * consisting of exactly <tt>2^dim-1</tt> cells. It produces the hypercube
   * with the interval [<i>left,right</i>] without the hypercube made out of
   * the interval [<i>(left+right)/2,right</i>] for each coordinate. If the
   * @p colorize flag is set, the @p boundary_ids of the surfaces are
   * assigned, such that the left boundary is 0, and the others are set with
   * growing number accordingly to the counterclockwise. Colorize option works
   * only with 2-dimensional problem. This function will create the classical
   * L-shape in 2d and it will look like the following in 3d:
   *
   * @image html hyper_l.png
   *
   * This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void hyper_L (Triangulation<dim> &tria,
                const double        left = -1.,
                const double        right= 1.,
                const bool          colorize = false);

  /**
   * Initialize the given Triangulation with a hypercube with a slit. In each
   * coordinate direction, the hypercube extends from @p left to @p right.
   *
   * In 2d, the split goes in vertical direction from <tt>x=(left+right)/2,
   * y=left</tt> to the center of the square at <tt>x=y=(left+right)/2</tt>.
   *
   * In 3d, the 2d domain is just extended in the <i>z</i>-direction, such
   * that a plane cuts the lower half of a rectangle in two.  This function is
   * declared to exist for triangulations of all space dimensions, but throws
   * an error if called in 1d.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void hyper_cube_slit (Triangulation<dim> &tria,
                        const double        left = 0.,
                        const double        right = 1.,
                        const bool          colorize = false);

  /**
   * Produce a hyper-shell, the region between two spheres around
   * <tt>center</tt>, with given <tt>inner_radius</tt> and
   * <tt>outer_radius</tt>. The number <tt>n_cells</tt> indicates the number
   * of cells of the resulting triangulation, i.e., how many cells form the
   * ring (in 2d) or the shell (in 3d).
   *
   * If the flag @p colorize is @p true, then the outer boundary will have the
   * indicator 1, while the inner boundary has id zero. In 3d, this applies to
   * both the faces and the edges of these boundaries. If the flag is @p
   * false, both have indicator zero.
   *
   * You should attach a SphericalManifold to the cells and faces for correct
   * placement of vertices upon refinement and to be able to use higher order
   * mappings. Alternatively, it is also possible to attach a
   * HyperShellBoundary to the inner and outer boundary. This will create
   * inferior meshes as described below.
   *
   * In 2d, the number <tt>n_cells</tt> of elements for this initial
   * triangulation can be chosen arbitrarily. If the number of initial cells
   * is zero (as is the default), then it is computed adaptively such that the
   * resulting elements have the least aspect ratio.
   *
   * In 3d, only certain numbers are allowed, 6 (or the default 0) for a
   * surface based on a hexahedron (i.e. 6 panels on the inner sphere extruded
   * in radial direction to form 6 cells), 12 for the rhombic dodecahedron,
   * and 96 (see below).
   *
   * While the SphericalManifold, that is demonstrated in the documentation of
   * the
   * @ref manifold "documentation module on manifolds",
   * creates reasonable meshes for any number of @p n_cells if attached to all
   * cells and boundaries, the situation is less than ideal when only
   * attaching a HyperShellBoundary. Then, only vertices on the boundaries are
   * placed at the correct distance from the center. As an example, the 3d
   * meshes give rise to the following meshes upon one refinement:
   *
   * @image html hypershell3d-6.png
   * @image html hypershell3d-12.png
   *
   * Neither of these meshes is particularly good since one ends up with
   * poorly shaped cells at the inner edge upon refinement. For example, this
   * is the middle plane of the mesh for the <code>n_cells=6</code>:
   *
   * @image html hyper_shell_6_cross_plane.png
   *
   * The mesh generated with <code>n_cells=12</code> is better but still not
   * good. As a consequence, you may also specify <code>n_cells=96</code> as a
   * third option. The mesh generated in this way is based on a once refined
   * version of the one with <code>n_cells=12</code>, where all internal nodes
   * are re-placed along a shell somewhere between the inner and outer
   * boundary of the domain. The following two images compare half of the
   * hyper shell for <code>n_cells=12</code> and <code>n_cells=96</code> (note
   * that the doubled radial lines on the cross section are artifacts of the
   * visualization):
   *
   * @image html hyper_shell_12_cut.png
   * @image html hyper_shell_96_cut.png
   *
   * @note This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void hyper_shell (Triangulation<dim> &tria,
                    const Point<dim>   &center,
                    const double        inner_radius,
                    const double        outer_radius,
                    const unsigned int  n_cells = 0,
                    bool                colorize = false);

  /**
   * Produce a half hyper-shell, i.e. the space between two circles in two
   * space dimensions and the region between two spheres in 3d, with given
   * inner and outer radius and a given number of elements for this initial
   * triangulation.  However, opposed to the previous function, it does not
   * produce a whole shell, but only one half of it, namely that part for
   * which the first component is restricted to non-negative values. The
   * purpose of this class is to enable computations for solutions which have
   * rotational symmetry, in which case the half shell in 2d represents a
   * shell in 3d.
   *
   * If the number of initial cells is zero (as is the default), then it is
   * computed adaptively such that the resulting elements have the least
   * aspect ratio.
   *
   * If colorize is set to true, the inner, outer, and the part of the
   * boundary where $x=0$, get indicator 0, 1, and 2, respectively. Otherwise
   * all indicators are set to 0.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void half_hyper_shell (Triangulation<dim> &tria,
                         const Point<dim>   &center,
                         const double        inner_radius,
                         const double        outer_radius,
                         const unsigned int  n_cells = 0,
                         const bool          colorize = false);


  /**
   * Produce a domain that is the intersection between a hyper-shell with
   * given inner and outer radius, i.e. the space between two circles in two
   * space dimensions and the region between two spheres in 3d, and the
   * positive quadrant (in 2d) or octant (in 3d). In 2d, this is indeed a
   * quarter of the full annulus, while the function is a misnomer in 3d
   * because there the domain is not a quarter but one eighth of the full
   * shell.
   *
   * If the number of initial cells is zero (as is the default), then it is
   * computed adaptively such that the resulting elements have the least
   * aspect ratio in 2d.
   *
   * If @p colorize is set to true, the inner, outer, left, and right boundary
   * get indicator 0, 1, 2, and 3 in 2d, respectively. Otherwise all
   * indicators are set to 0. In 3d indicator 2 is at the face x=0, 3 at y=0,
   * 4 at z=0.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void quarter_hyper_shell (Triangulation<dim> &tria,
                            const Point<dim>   &center,
                            const double        inner_radius,
                            const double        outer_radius,
                            const unsigned int  n_cells = 0,
                            const bool          colorize = false);

  /**
   * Produce a domain that is the space between two cylinders in 3d, with
   * given length, inner and outer radius and a given number of elements for
   * this initial triangulation. If @p n_radial_cells is zero (as is the
   * default), then it is computed adaptively such that the resulting elements
   * have the least aspect ratio. The same holds for @p n_axial_cells.
   *
   * @note Although this function is declared as a template, it does not make
   * sense in 1D and 2D.
   *
   * @note The triangulation needs to be void upon calling this function.
   */
  template <int dim>
  void cylinder_shell (Triangulation<dim> &tria,
                       const double        length,
                       const double        inner_radius,
                       const double        outer_radius,
                       const unsigned int  n_radial_cells = 0,
                       const unsigned int  n_axial_cells = 0);



  /**
   * Produce the volume or surface mesh of a torus. The axis of the torus is
   * the $y$-axis while the plane of the torus is the $x$-$z$ plane.
   *
   * If @p dim is 3, the mesh will be the volume of the torus. By default,
   * the boundary faces will have manifold id 0 and you should attach a
   * TorusManifold to it. The cells will have manifold id 1 and you should
   * attach a SphericalManifold to it.
   *
   * If @p dim is 2, the mesh will describe the surface of the torus. All
   * cells and faces will have manifold id 0 and you should attach a
   * TorusManifold to it.
   *
   * @param tria The triangulation to be filled.
   *
   * @param R The radius of the circle, which forms the middle line of the
   * torus containing the loop of cells. Must be greater than @p r.
   *
   * @param r The inner radius of the torus.
   *
   * @note Implemented for Triangulation<2,3> and Triangulation<3,3>.
   */
  template <int dim, int spacedim>
  void torus (Triangulation<dim,spacedim> &tria,
              const double R,
              const double r);



  /**
   * This class produces a square in the <i>xy</i>-plane with a circular hole
   * in the middle. Square and circle are centered at the origin. In 3d, this
   * geometry is extruded in $z$ direction to the interval $[0,L]$.
   *
   * @image html cubes_hole.png
   *
   * It is implemented in 2d and 3d, and takes the following arguments:
   *
   * @param triangulation The triangulation to be filled.
   * @param inner_radius  Radius of the internal hole.
   * @param outer_radius Half of the edge length of the square.
   * @param L  Extension in @p z-direction (only used in 3d).
   * @param repetitions Number of subdivisions along the @p z-direction.
   * @param colorize Whether to assign different boundary indicators to
   * different faces. The colors are given in lexicographic ordering for the
   * flat faces (0 to 3 in 2d, 0 to 5 in 3d) plus the curved hole (4 in 2d,
   * and 6 in 3d). If @p colorize is set to false, then flat faces get the
   * number 0 and the hole gets number 1.
   */
  template<int dim>
  void hyper_cube_with_cylindrical_hole (
    Triangulation<dim> &triangulation,
    const double        inner_radius = .25,
    const double        outer_radius = .5,
    const double        L = .5,
    const unsigned int  repetitions = 1,
    const bool          colorize = false);

  /**
   * Produce a ring of cells in 3d that is cut open, twisted and glued
   * together again. This results in a kind of moebius-loop.
   *
   * @param tria        The triangulation to be worked on.
   * @param n_cells     The number of cells in the loop. Must be greater than
   * 4.
   * @param n_rotations The number of rotations (Pi/2 each) to be performed
   * before gluing the loop together.
   * @param R           The radius of the circle, which forms the middle line
   * of the torus containing the loop of cells. Must be greater than @p r.
   * @param r           The radius of the cylinder bend together as loop.
   */
  void moebius (Triangulation<3,3> &tria,
                const unsigned int  n_cells,
                const unsigned int  n_rotations,
                const double        R,
                const double        r);

  ///@}

  /**
   * @name Creating meshes from other meshes
   */
  ///@{

  /**
   * Given the two triangulations specified as the first two arguments, create
   * the triangulation that contains the cells of both triangulation and store
   * it in the third parameter. Previous content of @p result will be deleted.
   *
   * This function is most often used to compose meshes for more complicated
   * geometries if the geometry can be composed of simpler parts for which
   * functions exist to generate coarse meshes.  For example, the channel mesh
   * used in step-35 could in principle be created using a mesh created by the
   * GridGenerator::hyper_cube_with_cylindrical_hole function and several
   * rectangles, and merging them using the current function. The rectangles
   * will have to be translated to the right for this, a task that can be done
   * using the GridTools::shift function (other tools to transform individual
   * mesh building blocks are GridTools::transform, GridTools::rotate, and
   * GridTools::scale).
   *
   * @note The two input triangulations must be coarse meshes that have no
   * refined cells.
   *
   * @note The function copies the material ids of the cells of the two input
   * triangulations into the output triangulation but it currently makes no
   * attempt to do the same for boundary ids. In other words, if the two
   * coarse meshes have anything but the default boundary indicators, then you
   * will currently have to set boundary indicators again by hand in the
   * output triangulation.
   *
   * @note For a related operation on refined meshes when both meshes are
   * derived from the same coarse mesh, see
   * GridGenerator::create_union_triangulation().
   */
  template <int dim, int spacedim>
  void
  merge_triangulations (const Triangulation<dim, spacedim> &triangulation_1,
                        const Triangulation<dim, spacedim> &triangulation_2,
                        Triangulation<dim, spacedim>       &result);

  /**
   * Given the two triangulations specified as the first two arguments, create
   * the triangulation that contains the finest cells of both triangulation
   * and store it in the third parameter. Previous content of @p result will
   * be deleted.
   *
   * @note This function is intended to create an adaptively refined
   * triangulation that contains the <i>most refined cells</i> from two input
   * triangulations that were derived from the <i>same</i> coarse grid by
   * adaptive refinement. This is an operation sometimes needed when one
   * solves for two variables of a coupled problem on separately refined
   * meshes on the same domain (for example because these variables have
   * boundary layers in different places) but then needs to compute something
   * that involves both variables or wants to output the result into a single
   * file. In both cases, in order not to lose information, the two solutions
   * can not be interpolated onto the respectively other mesh because that may
   * be coarser than the ones on which the variable was computed. Rather, one
   * needs to have a mesh for the domain that is at least as fine as each of
   * the two initial meshes. This function computes such a mesh.
   *
   * @note If you want to create a mesh that is the merger of two other coarse
   * meshes, for example in order to compose a mesh for a complicated geometry
   * from meshes for simpler geometries, then this is not the function for
   * you. Instead, consider GridGenerator::merge_triangulations().
   *
   * @pre Both of the source conditions need to be available entirely locally.
   * In other words, they can not be objects of type
   * parallel::distributed::Triangulation.
   */
  template <int dim, int spacedim>
  void
  create_union_triangulation (const Triangulation<dim, spacedim> &triangulation_1,
                              const Triangulation<dim, spacedim> &triangulation_2,
                              Triangulation<dim, spacedim>       &result);

  /**
   * This function creates a triangulation that consists of the same cells as
   * are present in the first argument, except those cells that are listed in
   * the second argument. The purpose of the function is to generate
   * geometries <i>subtractively</i> from the geometry described by an
   * existing triangulation. A prototypical case is a 2d domain with
   * rectangular holes. This can be achieved by first meshing the entire
   * domain and then using this function to get rid of the cells that are
   * located at the holes. Likewise, you could create the mesh that
   * GridGenerator::hyper_L() produces by starting with a
   * GridGenerator::hyper_cube(), refining it once, and then calling the
   * current function with a single cell in the second argument.
   *
   * @param[in] input_triangulation The original triangulation that serves as
   * the template from which the new one is to be created.
   * @param[in] cells_to_remove A list of cells of the triangulation provided
   * as first argument that should be removed (i.e., that should not show up
   * in the result.
   * @param[out] result The resulting triangulation that consists of the same
   * cells as are in @p input_triangulation, with the exception of the cells
   * listed in @p cells_to_remove.
   *
   * @pre Because we cannot create triangulations de novo that contain
   * adaptively refined cells, the input triangulation needs to have all of
   * its cells on the same level. Oftentimes, this will in fact be the
   * coarsest level, but it is allowed to pass in a triangulation that has
   * been refined <i>globally</i> a number of times. The output triangulation
   * will in that case simply be a mesh with only one level that consists of
   * the active cells of the input minus the ones listed in the second
   * argument. However, the input triangulation must not have been
   * <i>adaptively</i> refined.
   */
  template <int dim, int spacedim>
  void
  create_triangulation_with_removed_cells (const Triangulation<dim, spacedim> &input_triangulation,
                                           const std::set<typename Triangulation<dim, spacedim>::active_cell_iterator> &cells_to_remove,
                                           Triangulation<dim, spacedim>       &result);


  /**
   * Take a 2d Triangulation that is being extruded in z direction by the
   * total height of @p height using @p n_slices slices (minimum is 2). The
   * boundary indicators of the faces of @p input are going to be assigned to
   * the corresponding side walls in z direction. The bottom and top get the
   * next two free boundary indicators.
   *
   * @note The 2d input triangulation @p input must be a coarse mesh that has
   * no refined cells.
   */
  void
  extrude_triangulation (const Triangulation<2, 2> &input,
                         const unsigned int         n_slices,
                         const double               height,
                         Triangulation<3,3>        &result);

  /**
   * Given an input triangulation @p in_tria, this function makes a new flat
   * triangulation @p out_tria which contains a single level with all active
   * cells of the input triangulation. If @p spacedim1 and @p spacedim2 are
   * different, only the smallest spacedim components of the vertices are
   * copied over. This is useful to create a Triangulation<2,3> out of a
   * Triangulation<2,2>, or to project a Triangulation<2,3> into a
   * Triangulation<2,2>, by neglecting the z components of the vertices.
   *
   * No internal checks are performed on the vertices, which are assumed to
   * make sense topologically in the target @p spacedim2 dimensional space. If
   * this is not the case, you will encounter problems when using the
   * triangulation later on.
   *
   * All information about cell manifold_ids and material ids are copied from
   * one triangulation to the other, and only the boundary manifold_ids and
   * boundary_ids are copied over from the faces of @p in_tria to the faces of
   * @p out_tria. If you need to specify manifold ids on interior faces, they
   * have to be specified manually after the triangulation is created.
   *
   * This function will fail if the input Triangulation is of type
   * parallel::distributed::Triangulation, as well as when the input
   * Triangulation contains hanging nodes.
   *
   * @author Luca Heltai, 2014
   */
  template <int dim, int spacedim1, int spacedim2>
  void flatten_triangulation(const Triangulation<dim,spacedim1> &in_tria,
                             Triangulation<dim,spacedim2> &out_tria);

  ///@}

  /**
   * @name Creating lower-dimensional meshes from parts of higher-dimensional
   * meshes
   */
  ///@{

#ifdef _MSC_VER
  // Microsoft's VC++ has a bug where it doesn't want to recognize that
  // an implementation (definition) of the extract_boundary_mesh function
  // matches a declaration. This can apparently only be avoided by
  // doing some contortion with the return type using the following
  // intermediate type. This is only used when using MS VC++ and uses
  // the direct way of doing it otherwise
  template <template <int,int> class MeshType, int dim, int spacedim>
  struct ExtractBoundaryMesh
  {
    typedef
    std::map<typename MeshType<dim-1,spacedim>::cell_iterator,
        typename MeshType<dim,spacedim>::face_iterator>
        return_type;
  };
#endif

  /**
   * This function implements a boundary subgrid extraction.  Given a
   * <dim,spacedim>-Triangulation (the "volume mesh") the function extracts a
   * subset of its boundary (the "surface mesh").  The boundary to be
   * extracted is specified by a list of boundary_ids.  If none is specified
   * the whole boundary will be extracted. The function is used in step-38.
   *
   * The function also builds a mapping linking the cells on the surface mesh
   * to the corresponding faces on the volume one. This mapping is the return
   * value of the function.
   *
   * @note The function builds the surface mesh by creating a coarse mesh from
   * the selected faces of the coarse cells of the volume mesh. It copies the
   * boundary indicators of these faces to the cells of the coarse surface
   * mesh. The surface mesh is then refined in the same way as the faces of
   * the volume mesh are. In order to ensure that the surface mesh has the
   * same vertices as the volume mesh, it is therefore important that you
   * assign appropriate boundary objects through Triangulation::set_boundary()
   * to the surface mesh object before calling this function. If you don't,
   * the refinement will happen under the assumption that all faces are
   * straight (i.e using the StraightBoundary class) rather than any curved
   * boundary object you may want to use to determine the location of new
   * vertices.
   *
   * @tparam MeshType A type that satisfies the requirements of the
   * @ref ConceptMeshType "MeshType concept".
   * The map that is returned will be between cell iterators pointing into the
   * container describing the surface mesh and face iterators of the volume
   * mesh container. If MeshType is DoFHandler or hp::DoFHandler, then the
   * function will re-build the triangulation underlying the second argument
   * and return a map between appropriate iterators into the MeshType
   * arguments. However, the function will not actually distribute degrees of
   * freedom on this newly created surface mesh.
   *
   * @tparam dim The dimension of the cells of the volume mesh. For example,
   * if dim==2, then the cells are quadrilaterals that either live in the
   * plane, or form a surface in a higher-dimensional space. The dimension of
   * the cells of the surface mesh is consequently dim-1.
   * @tparam spacedim The dimension of the space in which both the volume and
   * the surface mesh live.
   *
   * @param[in] volume_mesh A container of cells that define the volume mesh.
   * @param[out] surface_mesh A container whose associated triangulation will
   * be built to consist of the cells that correspond to the (selected portion
   * of) the boundary of the volume mesh.
   * @param[in] boundary_ids A list of boundary indicators denoting that
   * subset of faces of volume cells for which this function should extract
   * the surface mesh. If left at its default, i.e., if the set is empty, then
   * the function operates on <i>all</i> boundary faces.
   *
   * @return A map that for each cell of the surface mesh (key) returns an
   * iterator to the corresponding face of a cell of the volume mesh (value).
   * The keys include both active and non-active cells of the surface mesh.
   * For dim=2 (i.e., where volume cells are quadrilaterals and surface cells
   * are lines), the order of vertices of surface cells and the corresponding
   * volume faces match. For dim=3 (i.e., where volume cells are hexahedra and
   * surface cells are quadrilaterals), the order of vertices may not match in
   * order to ensure that each surface cell has a right-handed coordinate
   * system when viewed from one of the two sides of the surface connecting
   * the cells of the surface mesh.
   *
   * @note The algorithm outlined above assumes that all faces on higher
   * refinement levels always have exactly the same boundary indicator as
   * their parent face. Consequently, we can start with coarse level faces and
   * build the surface mesh based on that. It would not be very difficult to
   * extend the function to also copy boundary indicators from finer level
   * faces to their corresponding surface mesh cells, for example to
   * accommodate different geometry descriptions in the case of curved
   * boundaries (but this is not currently implemented).
   */
  template <template <int,int> class MeshType, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename MeshType<dim-1,spacedim>::cell_iterator,
      typename MeshType<dim,spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<MeshType,dim,spacedim>::return_type
#endif
      extract_boundary_mesh (const MeshType<dim,spacedim>       &volume_mesh,
                             MeshType<dim-1,spacedim>           &surface_mesh,
                             const std::set<types::boundary_id> &boundary_ids
                             = std::set<types::boundary_id>());

  ///@}


  /**
   * @name Exceptions
   */
  ///@{


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

  ///@}
}



DEAL_II_NAMESPACE_CLOSE

#endif
