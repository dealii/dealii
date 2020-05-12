// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_grid_generator_h
#define dealii_grid_generator_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/function.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/grid/tria.h>

#include <array>
#include <map>

DEAL_II_NAMESPACE_OPEN

/**
 * This namespace provides a collection of functions for generating
 * triangulations for some basic geometries.
 *
 * Some of these functions receive a flag @p colorize (see
 * @ref GlossColorization "the glossary entry on colorization").
 * If this is set, parts of the boundary receive different
 * @ref GlossBoundaryIndicator "boundary indicators"
 * allowing them to be distinguished for the purpose of evaluating
 * different boundary conditions.
 *
 * If the domain is curved, each of the domain parts that should be
 * refined by following an appropriate Manifold description will
 * receive a different
 * @ref GlossManifoldIndicator "manifold indicator",
 * and the correct Manifold descriptor will be attached to
 * the Triangulation. Notice that if you later transform the
 * triangulation, you have to make sure you attach the correct new Manifold
 * to the triangulation.
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
   * zero and unity, then producing the unit hypercube.
   *
   * If the argument @p colorize is false, then all boundary indicators are
   * set to zero (the default boundary indicator) for 2d and 3d. If it is
   * true, the boundary is
   * @ref GlossColorization "colorized"
   * as in hyper_rectangle(). In 1d the
   * indicators are always colorized, see hyper_rectangle().
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
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim, int spacedim>
  void
  hyper_cube(Triangulation<dim, spacedim> &tria,
             const double                  left     = 0.,
             const double                  right    = 1.,
             const bool                    colorize = false);

  /**
   * Create a $d$-<a href="https://en.wikipedia.org/wiki/Simplex">simplex</a>
   * (i.e., a triangle in 2d, or a tetrahedron in 3d) with
   * $d+1$ corners. Since deal.II does not support triangular and
   * tetrahedral cells, the simplex described by the input arguments
   * is subdivided into quadrilaterals and hexahedra by adding edge,
   * face, and simplex midpoints, resulting in a mesh that consists of
   * $d+1$ quadrilateral or hexahedral cells.
   *
   * The @p vertices argument contains a vector with all d+1 vertices defining
   * the corners of the simplex. They must be given in an order such that the
   * vectors from the first vertex to each of the others form a right-handed
   * system.
   *
   * The meshes generated in two and three dimensions are:
   *
   * @image html simplex_2d.png
   * @image html simplex_3d.png
   *
   * @param tria The triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * @param vertices The dim+1 corners of the simplex.
   *
   * @note Implemented for <tt>Triangulation@<2,2@></tt>,
   * <tt>Triangulation@<3,3@></tt>.
   *
   * @author Guido Kanschat
   */
  template <int dim>
  void
  simplex(Triangulation<dim, dim> &      tria,
          const std::vector<Point<dim>> &vertices);

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
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   *
   * @param tria The triangulation to create. It needs to be empty upon
   * calling this function.
   *
   * @param repetitions A vector of @p dim positive values denoting the number
   * of cells to generate in that direction. For example, the first component of
   * the vector specifies the number of cells in the x-direction, the second
   * component specifies that for the y-direction for dim=2.
   *
   * @param left Lower bound for the interval used to create the hyper cube.
   *
   * @param right Upper bound for the interval used to create the hyper cube.
   *
   * @param colorize Assign different boundary ids if set to true.
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_cube(Triangulation<dim, spacedim> &tria,
                        const unsigned int            repetitions,
                        const double                  left     = 0.,
                        const double                  right    = 1.,
                        const bool                    colorize = false);

  /**
   * Create a coordinate-parallel brick from the two diagonally opposite
   * corner points @p p1 and @p p2.
   *
   * If the @p colorize flag is <code>true</code>, then the @p boundary_ids of
   * the boundary faces are assigned, such that the lower one in @p
   * x-direction is 0, the upper one is 1. The indicators for the surfaces in
   * @p y-direction are 2 and 3, the ones for @p z are 4 and 5. This
   * corresponds to the numbers of faces of the unit square of cube as laid
   * out in the documentation of the GeometryInfo class; see also
   * @ref GlossColorization "the glossary entry on colorization".
   * Importantly,
   * however, in 3d
   * @ref GlossColorization "colorization"
   * does not set @p
   * boundary_ids of <i>edges</i>, but only of <i>faces</i>, because each
   * boundary edge is shared between two faces and it is not clear how the
   * boundary id of an edge should be set in that case.
   *
   * Additionally, if @p colorize is @p true, material ids are assigned to the
   * cells according to the octant their center is in: being in the right half
   * space for any coordinate direction <i>x<sub>i</sub></i> adds
   * 2<sup>i</sup>. For instance, a cell with center point (1,-1,1) yields a
   * material id 5, assuming that the center of the hyper rectangle lies at
   * the origin. No manifold id is set for the cells.
   *
   * If @p dim < @p spacedim, this will create a @p dim dimensional object in
   * the first @p dim coordinate directions embedded into the @p spacedim
   * dimensional space with the remaining entries set to zero. For example, a
   * <tt>Triangulation@<2,3@></tt> will be a rectangle in the xy plane with
   * z=0, defined by the two opposing corners @p p1 and @p p2.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim, int spacedim>
  void
  hyper_rectangle(Triangulation<dim, spacedim> &tria,
                  const Point<dim> &            p1,
                  const Point<dim> &            p2,
                  const bool                    colorize = false);

  /**
   * Create a coordinate-parallel brick from the two diagonally opposite
   * corner points @p p1 and @p p2. The number of cells in coordinate
   * direction @p i is given by the integer <tt>repetitions[i]</tt>.
   *
   * To get cells with an aspect ratio different from that of the domain, use
   * different numbers of subdivisions, given by @p repetitions, in different
   * coordinate directions. The minimum number of subdivisions in each
   * direction is 1.
   *
   * If the @p colorize flag is <code>true</code>, then the @p boundary_ids of
   * the surfaces are assigned, such that the lower one in @p x-direction is
   * 0, the upper one is 1 (the left and the right vertical face). The
   * indicators for the surfaces in @p y-direction are 2 and 3, the ones for
   * @p z are 4 and 5.  Additionally, material ids are assigned to the cells
   * according to the octant their center is in: being in the right half plane
   * for any coordinate direction <i>x<sub>i</sub></i> adds 2<sup>i</sup> (see
   * @ref GlossColorization "the glossary entry on colorization").
   * For
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
   * @param tria The triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * @param repetitions A vector of @p dim positive values denoting the number
   * of cells to generate in that direction.
   *
   * @param p1 First corner point.
   *
   * @param p2 Second corner opposite to @p p1.
   *
   * @param colorize Assign different boundary ids if set to true. The same
   * comments apply as for the hyper_rectangle() function.
   *
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_rectangle(Triangulation<dim, spacedim> &   tria,
                             const std::vector<unsigned int> &repetitions,
                             const Point<dim> &               p1,
                             const Point<dim> &               p2,
                             const bool                       colorize = false);

  /**
   * Like the previous function. However, here the second argument does not
   * denote the number of subdivisions in each coordinate direction, but a
   * sequence of step sizes for each coordinate direction. The domain will
   * therefore be subdivided into <code>step_sizes[i].size()</code> cells in
   * coordinate direction <code>i</code>, with width
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
  subdivided_hyper_rectangle(Triangulation<dim> &                    tria,
                             const std::vector<std::vector<double>> &step_sizes,
                             const Point<dim> &                      p_1,
                             const Point<dim> &                      p_2,
                             const bool colorize = false);

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
  subdivided_hyper_rectangle(Triangulation<dim> &                    tria,
                             const std::vector<std::vector<double>> &spacing,
                             const Point<dim> &                      p,
                             const Table<dim, types::material_id> &material_id,
                             const bool colorize = false);

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
   * @param tria The triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * @param holes Positive number of holes in each of the dim directions.
   * @author Guido Kanschat
   */
  template <int dim, int spacedim>
  void
  cheese(Triangulation<dim, spacedim> &   tria,
         const std::vector<unsigned int> &holes);

  /**
   * \brief Rectangular plate with an (offset) cylindrical hole.
   *
   * Generate a rectangular plate with an (offset) cylindrical hole. The
   * geometry consists of 2 regions:
   * The first is a square region with length @p outer_radius and a hole of radius @p inner_radius .
   * Cells in this region will have TransfiniteInterpolationManifold with
   * manifold id @p tfi_manifold_id attached to them. Additionally, the boundary
   * faces of the hole will be associated with a PolarManifold (in 2D) or
   * CylindricalManifold (in 3D). The center of this
   * region can be prescribed via @p center , namely the axis of the hole will
   * be located at @p center .
   * The second region describes the remainder of the bulk material. It is
   * specified via padding
   * parameters @p pad_bottom, @p padding_top, @p padding_left and @p padding_right.
   * All cells in this region will have a FlatManifold attached to them.
   * The final width of the plate will be <code>padding_left + 2*outer_radius +
   * padding_right</code>, while its length is <code>padding_top +
   * 2*outer_radius + padding_bottom</code>.
   *
   * Here is the non-symmetric grid (after one global refinement, colored
   * according to manifold id) in 2D and 3D, respectively:
   *
   * \htmlonly <style>div.image
   * img[src="plate_with_a_hole.png"]{width:25%;}</style> \endhtmlonly
   * @image html plate_with_a_hole.png
   * \htmlonly <style>div.image
   * img[src="plate_with_a_hole_3D.png"]{width:25%;}</style> \endhtmlonly
   * @image html plate_with_a_hole_3D.png
   *
   * In 3D, triangulation will be extruded in the z-direction by the total
   * height of @p L using @p n_slices slices (minimum is 2).

   * If the @p colorize flag is <code>true</code>, the boundary_ids of the
   * boundary faces are assigned such that the lower one in the x-direction is
   * 0, and the upper one is 1 (see
   * @ref GlossColorization "the glossary entry on colorization").
   * The indicators for the surfaces in the y-direction are 2 and 3, and the
   * ones for the z-direction are 5 and 6. The hole boundary has indicator 4.
   *
   * @p tria is the triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * @author Denis Davydov, 2018
   */
  template <int dim>
  void
  plate_with_a_hole(Triangulation<dim> &     tria,
                    const double             inner_radius      = 0.4,
                    const double             outer_radius      = 1.,
                    const double             pad_bottom        = 2.,
                    const double             pad_top           = 2.,
                    const double             pad_left          = 1.,
                    const double             pad_right         = 1.,
                    const Point<dim>         center            = Point<dim>(),
                    const types::manifold_id polar_manifold_id = 0,
                    const types::manifold_id tfi_manifold_id   = 1,
                    const double             L                 = 1.,
                    const unsigned int       n_slices          = 2,
                    const bool               colorize          = false);

  /**
   * Generate a grid consisting of a channel with a cylinder. This is a common
   * benchmark for Navier-Stokes solvers. The geometry consists of a channel
   * of size $[0, 2.2] \times [0, 0.41] \times [0, 0.41] $ (where the $z$
   * dimension is omitted in 2D) with a cylinder, parallel to the $z$ axis
   * with diameter $0.1$, centered at $(0.2, 0.2, 0)$. The channel has three
   * distinct regions:
   * <ol>
   *   <li>If @p n_shells is greater than zero, then there are that many shells
   *   centered around the cylinder,</li>
   *   <li>a blending region between the shells and the rest of the
   *   triangulation, and</li>
   *   <li>a bulk region consisting of Cartesian cells.</li>
   * </ol>
   * Since the cylinder is slightly offset from the center of the channel,
   * this geometry results in vortex shedding at moderate Reynolds
   * numbers. Here is the grid (without additional global refinement) in 2D:
   *
   * @image html channel_with_cylinder_2d.png
   *
   * and in 3D:
   *
   * @image html channel_with_cylinder_3d.png
   *
   * The resulting Triangulation uses three manifolds: a PolarManifold (in 2D)
   * or CylindricalManifold (in 3D) with manifold id $0$, a
   * TransfiniteInterpolationManifold with manifold id $1$, and a FlatManifold
   * everywhere else. For more information on this topic see
   * @ref GlossManifoldIndicator "the glossary entry on manifold indicators".
   * The
   * cell faces on the cylinder and surrounding shells have manifold ids of
   * $0$, while the cell volumes adjacent to the shells (or, if they do not
   * exist, the cylinder) have a manifold id of $1$. Put another way: this
   * grid uses TransfiniteInterpolationManifold to smoothly transition from
   * the shells (generated with GridGenerator::concentric_hyper_shells) to the
   * bulk region. All other cell volumes and faces have manifold id
   * numbers::flat_manifold_id and use FlatManifold. All cells with id
   * numbers::flat_manifold_id are rectangular prisms aligned with the
   * coordinate axes.
   *
   * The picture below shows part of the 2D grid (using all default arguments
   * to this function) after two global refinements. The cells with manifold
   * id $0$ are orange (the polar manifold id), cells with manifold id $1$ are
   * yellow (the transfinite interpolation manifold id), and the ones with
   * manifold id numbers::flat_manifold_id are cyan:
   *
   * @image html channel_with_cylinder_2d_manifolds.png
   *
   * @param tria Triangulation to be created. Must be empty upon calling this
   * function.
   *
   * @param shell_region_width Width of the layer of shells around the cylinder.
   * This value should be between $0$ and $0.05$; the default value is $0.03$.
   *
   * @param n_shells Number of shells to use in the shell layer.
   *
   * @param skewness Parameter controlling how close the shells are
   * to the cylinder: see the mathematical definition given in
   * GridGenerator::concentric_hyper_shells.
   *
   * @param colorize Assign different boundary ids if set to true. For more
   * information on boundary indicators see
   * @ref GlossBoundaryIndicator "this glossary entry".
   * The left boundary (at $x = 0$) is assigned an id of $0$, the right
   * boundary (at $x = 2.2$) is assigned an id of $1$, the cylinder boundary
   * is assigned an id of $2$, and the channel walls are assigned an id of
   * $3$.
   *
   * See the original paper for more information:
   * @code{.bib}
   * @inbook{schafer1996,
   * author    = {Sch{\"a}fer, M. and Turek, S. and Durst, F. and Krause, E.
   *              and Rannacher, R.},
   * title     = {Benchmark Computations of Laminar Flow Around a Cylinder},
   * bookTitle = {Flow Simulation with High-Performance Computers II: DFG
   *              Priority Research Programme Results 1993--1995},
   * year      = {1996},
   * publisher = {Vieweg+Teubner Verlag},
   * address   = {Wiesbaden},
   * pages     = {547--566},
   * isbn      = {978-3-322-89849-4},
   * doi       = {10.1007/978-3-322-89849-4_39},
   * url       = {https://doi.org/10.1007/978-3-322-89849-4_39}
   * }
   * @endcode
   */
  template <int dim>
  void
  channel_with_cylinder(Triangulation<dim> &tria,
                        const double        shell_region_width = 0.03,
                        const unsigned int  n_shells           = 2,
                        const double        skewness           = 2.0,
                        const bool          colorize           = false);

  /**
   * A general @p dim -dimensional cell (a segment if dim is 1, a quadrilateral
   * if @p dim is 2, or a hexahedron if @p dim is 3) immersed in a
   * @p spacedim -dimensional space. It is the responsibility of the user to
   * provide the vertices in the right order (see the documentation of the
   * GeometryInfo class) because the vertices are stored in the same order as
   * they are given. It is also important to make sure that the volume of the
   * cell is positive.
   *
   * If the argument @p colorize is false, then all boundary indicators are
   * set to zero for 2d and 3d. If it is true, the boundary is colorized as in
   * hyper_rectangle() (see
   * @ref GlossColorization "the glossary entry on colorization").
   * In 1d the
   * indicators are always colorized, see hyper_rectangle().
   *
   * @param tria The triangulation that will be created
   * @param vertices The 2^dim vertices of the cell
   * @param colorize If true, set different boundary ids.
   *
   * @author Bruno Turcksin
   */
  template <int dim, int spacedim>
  void
  general_cell(Triangulation<dim, spacedim> &      tria,
               const std::vector<Point<spacedim>> &vertices,
               const bool                          colorize = false);

  /**
   * A parallelogram. The first corner point is the origin. The next @p dim
   * vertices are the ones given in the second argument and the last vertex
   * will be the sum of the two vectors connecting the origin to those
   * points. Colorizing is done in the same way as in hyper_rectangle().
   *
   * @note This function is implemented in 2d only.
   *
   * @param tria The triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * @param corners Second and third vertices of the parallelogram.
   *
   * @param colorize Assign different boundary ids if true. (see
   * @ref GlossColorization "the glossary entry on colorization").
   */
  template <int dim>
  void
  parallelogram(Triangulation<dim> &tria,
                const Point<dim> (&corners)[dim],
                const bool colorize = false);

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
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  parallelepiped(Triangulation<dim> &tria,
                 const Point<dim> (&corners)[dim],
                 const bool colorize = false);

  /**
   * A subdivided parallelepiped. The first corner point is the origin. The @p
   * dim adjacent points are vectors describing the edges of the
   * parallelepiped with respect to the origin. Additional points are sums of
   * these dim vectors. The variable @p n_subdivisions designates the number
   * of subdivisions in each of the @p dim directions. Colorizing is done
   * according to hyper_rectangle().
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  subdivided_parallelepiped(Triangulation<dim> &tria,
                            const unsigned int  n_subdivisions,
                            const Point<dim> (&corners)[dim],
                            const bool colorize = false);

  /**
   * A subdivided parallelepiped, i.e., the same as above, but where the
   * number of subdivisions in each of the @p dim directions may vary.
   * Colorizing is done according to hyper_rectangle().
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  subdivided_parallelepiped(Triangulation<dim> &tria,
#ifndef _MSC_VER
                            const unsigned int (&n_subdivisions)[dim],
#else
                            const unsigned int *n_subdivisions,
#endif
                            const Point<dim> (&corners)[dim],
                            const bool colorize = false);

  /**
   * A subdivided parallelepiped.
   *
   * @param tria The triangulation to be created. It needs to be empty upon
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
   * @param colorize Assign different boundary ids if set to true (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * @note Implemented for all combinations of @p dim and @p spacedim.
   *
   * @note You likely need to help the compiler by explicitly specifying the
   * two template parameters when calling this function.
   */
  template <int dim, int spacedim>
  void
  subdivided_parallelepiped(Triangulation<dim, spacedim> &              tria,
                            const Point<spacedim> &                     origin,
                            const std::array<Tensor<1, spacedim>, dim> &edges,
                            const std::vector<unsigned int> &subdivisions = {},
                            const bool                       colorize = false);

  /**
   * Hypercube with a layer of hypercubes around it. Parameters @p left and
   * @p right give the lower and upper bound of the inner hypercube in all
   * coordinate directions.  @p thickness marks the size of the layer cells.
   *
   * If the flag @p colorize is set, the outer cells get material ids
   * according to the following scheme: extending over the inner cube in (+/-)
   * x-direction 1/2, y-direction 4/8, z-direction 16/32. A bitwise OR operation
   * is used to get these values at the corners and edges (3d), (see also
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * Presently only available in 2d and 3d.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  enclosed_hyper_cube(Triangulation<dim> &tria,
                      const double        left      = 0.,
                      const double        right     = 1.,
                      const double        thickness = 1.,
                      const bool          colorize  = false);

  /**
   * Initialize the given triangulation with several
   * @ref GlossCoarseMesh "coarse mesh cells"
   * that cover a hyperball, i.e. a circle in 2d or a
   * ball in 3d, around @p center with given @p radius. The function is
   * used in step-6.
   *
   * In order to avoid degenerate cells at the boundaries, the circle is
   * triangulated by five cells, whereas in 3d the ball is subdivided by
   * seven cells. Specifically, these
   * cells are one cell in the center plus one "cap" cell on each of the faces
   * of this center cell. This ensures that under repeated refinement, none
   * of the cells at the outer boundary will degenerate to have an interior
   * angle approaching 180 degrees, as opposed to the case where one might
   * start with just one square (or cube) to approximate the domain.
   * The diameter of the
   * center cell is chosen so that the aspect ratio of the boundary cells
   * after one refinement is optimized.
   *
   * This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * By default, the manifold_id is set to 0 on the boundary faces, 1 on the
   * boundary cells, and types::flat_manifold_id on the central cell and on
   * internal faces.
   *
   * A SphericalManifold is attached by default to the boundary faces for
   * correct placement of boundary vertices upon refinement and to be able to
   * use higher order mappings. However, it turns out that this strategy may
   * not be the optimal one to create a good a mesh for a hyperball. The
   * "Possibilities for extensions" section of step-6 has an extensive
   * discussion of how one would construct better meshes and what one needs to
   * do for it. Setting the argument
   * `attach_spherical_manifold_on_boundary_cells` to true attaches a
   * SphericalManifold manifold also to the cells adjacent to the boundary, and
   * not only to the boundary faces.
   *
   * @note Since this is likely one of the earliest functions users typically
   *   consider to create meshes with curved boundaries, let us also comment
   *   on one aspect that is often confusing: Namely, that what one sees is not
   *   always what is actually happening. Specifically, if you output the coarse
   *   mesh with a function such as GridOut::write_vtk() using default options,
   *   then one doesn't generally get to see curved faces at the boundary.
   *   That's because most file formats by default only store vertex locations,
   *   with the implicit understanding that cells are composed from these
   *   vertices and bounded by straight edges. At the same time, the fact
   *   that this function attaches a SphericalManifold object to the boundary
   *   faces means that at least *internally*, edges really are curved. If
   *   you want to see them that way, you need to make sure that the function
   *   you use to output the mesh actually plots boundary faces as curved
   *   lines rather than straight lines characterized by only the locations
   *   of the two end points. For example, GridOut::write_gnuplot() can do
   *   that if you set the corresponding flag in the GridOutFlags::Gnuplot
   *   structure. It is, however, an entirely separate consideration whether
   *   you are actually *computing* on curved cells. In typical finite
   *   element computations, one has to compute integrals and these are
   *   computed by transforming back actual cells using a mapping to the
   *   reference cell. What mapping is used determines what shape the
   *   cells have for these internal computations: For example, with the
   *   widely used $Q_1$ mapping (implicitly used in step-6), integration
   *   always happens on cells that are assumed to have straight boundaries
   *   described by only the vertex locations. In other words, if such a
   *   mapping is used, then the cells of the domain really do have
   *   straight edges, regardless of the manifold description attached
   *   to these edges and regardless of the flags given when generating
   *   output. As a consequence of all of this, it is important to
   *   distinguish three things: (i) the manifold description attached to an
   *   object in the mesh; (ii) the mapping used in integration; and (iii) the
   *   style used in outputting graphical information about the mesh. All of
   *   these can be chosen more or less independently of each other, and
   *   what you see visualized is not necessarily exactly what is
   *   happening.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  hyper_ball(Triangulation<dim> &tria,
             const Point<dim> &  center = Point<dim>(),
             const double        radius = 1.,
             const bool attach_spherical_manifold_on_boundary_cells = false);

  /**
   * This is an alternative to hyper_ball with 12 cells in 2d and 32 cells in
   * 3d, which provides a better balance between the size of the cells around
   * the outer curved boundaries and the cell in the interior. The mesh is
   * based on the cells used by GridGenerator::quarter_hyper_ball() with
   * appropriate copies and rotations to fill the whole ball.
   *
   * The following pictures show the resulting mesh in 2D (left) and 3D:
   * <table align="center" class="doxtable">
   *   <tr>
   *     <td>
   *       \htmlonly <style>div.image
   *         img[src="hyper_ball_balanced_2d.png"]{width:40%}</style>
   *       \endhtmlonly
   *       @image html hyper_ball_balanced_2d.png
   *     </td>
   *     <td>
   *       \htmlonly <style>div.image
   *         img[src="hyper_ball_balanced_3d.png"]{width:40%}</style>
   *       \endhtmlonly
   *       @image html hyper_ball_balanced_3d.png
   *     </td>
   *   </tr>
   * </table>
   *
   * By default, the manifold_id is set to 0 on the boundary faces, 1 on the
   * boundary cells, and types::flat_manifold_id on the central cell and on
   * internal faces.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  hyper_ball_balanced(Triangulation<dim> &tria,
                      const Point<dim> &  center = Point<dim>(),
                      const double        radius = 1.);

  /**
   * Creates a hyper sphere, i.e., a surface of a ball in @p spacedim
   * dimensions. This function only exists for dim+1=spacedim in 2 and 3 space
   * dimensions. (To create a mesh of a ball, use GridGenerator::hyper_ball().)
   *
   * By default, all manifold ids of the triangulation are set to zero, and a
   * SphericalManifold is attached to the grid.
   *
   * The following pictures are generated with:
   * @code
   * Triangulation<2,3>   triangulation;
   * GridGenerator::hyper_sphere(triangulation);
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
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */

  template <int spacedim>
  void hyper_sphere(Triangulation<spacedim - 1, spacedim> &tria,
                    const Point<spacedim> &center = Point<spacedim>(),
                    const double           radius = 1.);

  /**
   * This class produces a hyper-ball intersected with the positive orthant
   * relative to @p center, which contains three elements in 2d and four in
   * 3d. The interior points of the mesh are chosen to balance the minimal
   * singular value of the Jacobian of the mapping from reference to real
   * coordinates among the cells around the interior point, which corresponds
   * to a high mesh quality.
   *
   * The boundary indicators for the final triangulation are 0 for the curved
   * boundary and 1 for the cut plane. The manifold id for the curved boundary
   * is set to zero, and a SphericalManifold is attached to it.
   *
   * The resulting grid in 2D and 3D looks as follows:
   * <table align="center" class="doxtable">
   *   <tr>
   *     <td>
   *       \htmlonly <style>div.image
   *         img[src="quarter_hyper_ball_2d.png"]{width:50%}</style>
   *       \endhtmlonly
   *       @image html quarter_hyper_ball_2d.png
   *     </td>
   *     <td>
   *       \htmlonly <style>div.image
   *         img[src="quarter_hyper_ball_3d.png"]{width:46%}</style>
   *       \endhtmlonly
   *       @image html quarter_hyper_ball_3d.png
   *     </td>
   *   </tr>
   * </table>
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  quarter_hyper_ball(Triangulation<dim> &tria,
                     const Point<dim> &  center = Point<dim>(),
                     const double        radius = 1.);

  /**
   * This class produces a half hyper-ball around @p center, which contains
   * four elements in 2d and 6 in 3d. The cut plane is perpendicular to the
   * <i>x</i>-axis.
   *
   * The boundary indicators for the final triangulation are 0 for the curved
   * boundary and 1 for the cut plane. The manifold id for the curved boundary
   * is set to zero, and a SphericalManifold is attached to it.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  half_hyper_ball(Triangulation<dim> &tria,
                  const Point<dim> &  center = Point<dim>(),
                  const double        radius = 1.);

  /**
   * Create a @p dim dimensional cylinder where the $x$-axis serves as
   * the axis of the cylinder. For the purposes of this function, a
   * cylinder is defined as a (@p dim - 1) dimensional disk of given
   * @p radius, extruded along the axis of the cylinder (which is the
   * first coordinate direction). Consequently, in three dimensions,
   * the cylinder extends from `x=-half_length` to `x=+half_length`
   * and its projection into the @p yz-plane is a circle of radius @p
   * radius. In two dimensions, the cylinder is a rectangle from
   * `x=-half_length` to `x=+half_length` and from `y=-radius` to
   * `y=radius`.
   *
   * The boundaries are colored according to the following scheme: 0 for the
   * hull of the cylinder, 1 for the left hand face and 2 for the right hand
   * face (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * The manifold id for the hull of the cylinder is set to zero, and a
   * CylindricalManifold is attached to it.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  cylinder(Triangulation<dim> &tria,
           const double        radius      = 1.,
           const double        half_length = 1.);

  /**
   * Create a cut cone around the x-axis.  The cone extends from
   * <tt>x=-half_length</tt> to <tt>x=half_length</tt> and its projection into
   * the @p yz-plane is a circle of radius @p radius_0 at
   * <tt>x=-half_length</tt> and a circle of radius @p radius_1 at
   * <tt>x=+half_length</tt>. In between the radius is linearly decreasing.
   *
   * In two dimensions, the cone is a trapezoid from <tt>x=-half_length</tt>
   * to <tt>x=+half_length</tt> and from <tt>y=-radius_0</tt> to
   * <tt>y=radius_0</tt> at <tt>x=-half_length</tt> and from
   * <tt>y=-radius_1</tt> to <tt>y=radius_1</tt> at <tt>x=+half_length</tt>.
   * In between the range of <tt>y</tt> is linearly decreasing.
   *
   * The boundaries are colored according to the following scheme: 0 for the
   * hull of the cone, 1 for the left hand face, and 2 for the right hand face
   * (see
   * @ref GlossColorization "the glossary entry on colorization").
   * Both the boundary indicators and the manifold indicators are set.
   *
   * In three dimensions, the manifold id of the hull is set to zero, and a
   * CylindricalManifold is attached to it.
   *
   * Here are the grids in 2D and 3D after two mesh refinements:
   *
   * @image html truncated_cone_2d.png
   * @image html truncated_cone_3d.png
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   *
   * @author Markus B&uuml;rg, 2009
   */
  template <int dim>
  void
  truncated_cone(Triangulation<dim> &tria,
                 const double        radius_0    = 1.0,
                 const double        radius_1    = 0.5,
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
   * @param colorize_cells If colorization is enabled, then the material id of
   * a cells corresponds to the leg it is in. The id of the center cell is
   * zero, and then the legs are numbered starting at one (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * Examples in two and three dimensions are:
   *
   * @image html hyper_cross_2d.png
   * @image html hyper_cross_3d.png
   *
   * @author Guido Kanschat
   */
  template <int dim, int spacedim>
  void
  hyper_cross(Triangulation<dim, spacedim> &   tria,
              const std::vector<unsigned int> &sizes,
              const bool                       colorize_cells = false);

  /**
   * Initialize the given triangulation with a hyper-L (in 2d or 3d)
   * consisting of exactly <tt>2^dim-1</tt> cells. It produces the
   * hypercube with the interval [<i>left,right</i>] without the
   * hypercube made out of the interval [<i>(left+right)/2,right</i>]
   * for each coordinate. Because the domain is about the simplest one
   * with a reentrant (i.e., non-convex) corner, solutions of many
   * partial differential equations have singularities at this
   * corner. That is, at the corner, the gradient or a higher
   * derivative (depending on the boundary conditions chosen) does not
   * remain bounded. As a consequence, this domain is often used to
   * test convergence of schemes when the solution lacks regularity.
   *
   * If the @p colorize flag is <code>true</code>, the @p boundary_ids of the
   * surfaces are assigned such that the left boundary is 0 and the others are
   * assigned counterclockwise in ascending order (see
   * @ref GlossColorization "the glossary entry on colorization"). The @p
   * colorize option only works in two dimensions.
   *
   * This function will create the classical L-shape in 2d
   * and it will look like the following in 3d:
   *
   * @image html hyper_l.png
   *
   * @note The 3d domain is also often referred to as the "Fichera corner",
   * named after Gaetano Fichera (1922-1996) who first computed an
   * approximation of the corner singularity exponent of the lowest
   * eigenfunction of the domain.
   *
   * This function exists for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  hyper_L(Triangulation<dim> &tria,
          const double        left     = -1.,
          const double        right    = 1.,
          const bool          colorize = false);

  /**
   * Initialize the given triangulation in 2D or 3D with a generalized
   * subdivided hyper-L.
   *
   * This function produces a subdivided hyper rectangle with dimensions given
   * by @p bottom_left and @p top_right, with the given number of
   * subdivisions in each direction given in the vector @p repetitions,
   * and with a number of cells removed, given in the vector @p n_cells_to_remove.
   * Note that @p n_cells_to_remove contains integers, meaning that its entries
   * can be both positive and negative. A positive number denotes
   * cutting away cells in the 'positive' orientation, for example
   * left to right in the x-direction, bottom to top in
   * the y-direction, and front to back in the z-direction. A negative number
   * denotes cutting away cells in the reverse direction, so right to left,
   * top to bottom, and back to front.
   *
   * This function may be used to generate a mesh for a backward
   * facing step, a useful domain for benchmark problems in fluid dynamics.
   * The first image is a backward facing step in 3D, generated by
   * removing all cells in the z-direction, and 2 cells in the
   * positive x- and y-directions:
   * @image html subdivided_hyper_L_3d.png
   * And in 2D, we can cut away 1 cell in the negative x-direction, and 2 cells
   * in the negative y-direction:
   * @image html subdivided_hyper_L_2d.png
   *
   * @note This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1D.
   *
   * @author Mae Markowski, 2019
   */
  template <int dim, int spacedim>
  void
  subdivided_hyper_L(Triangulation<dim, spacedim> &   tria,
                     const std::vector<unsigned int> &repetitions,
                     const Point<dim> &               bottom_left,
                     const Point<dim> &               top_right,
                     const std::vector<int> &         n_cells_to_remove);

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
   * If @p colorize is set to @p true, the faces forming the slit are marked
   * with boundary id 1 and 2, respectively (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  hyper_cube_slit(Triangulation<dim> &tria,
                  const double        left     = 0.,
                  const double        right    = 1.,
                  const bool          colorize = false);

  /**
   * Produce a hyper-shell, the region between two spheres around
   * <tt>center</tt>, with given <tt>inner_radius</tt> and
   * <tt>outer_radius</tt>. The number <tt>n_cells</tt> indicates the number
   * of cells of the resulting triangulation, i.e., how many cells form the
   * ring (in 2d) or the shell (in 3d).
   *
   * If the flag @p colorize is <code>true</code>, then the outer boundary
   * will have the indicator 1 while the inner boundary has id zero. In 3d,
   * this applies to both the faces and the edges of these boundaries. If the
   * flag is @p false, both have indicator zero (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * All manifold ids are set to zero, and a SphericalManifold is attached to
   * every cell and face of the triangulation.
   *
   * In 2d, the number <tt>n_cells</tt> of elements for this initial
   * triangulation can be chosen arbitrarily. If the number of initial cells
   * is zero (as is the default), then it is computed adaptively such that the
   * resulting elements have the least aspect ratio.
   *
   * In 3d, only certain numbers are allowed:
   * <ul>
   * <li> 6 (or the default 0) for a surface based on a hexahedron (i.e. 6
   *      panels on the inner sphere extruded in radial direction to form 6
   *      cells),
   * <li> 12 for the rhombic dodecahedron,
   * <li> 24 for the hexahedron-based surface refined once in the azimuthal
   *      directions but not in the radial direction,
   * <li> 48 for the rhombic dodecahedron refined once in the azimuthal
   *      directions but not in the radial direction,
   * <li> 96 for the rhombic dodecahedron refined once. This choice dates from
   *      an older version of deal.II before the Manifold classes were
   *      implemented: today this choce is equivalent to the rhombic
   *      dodecahedron after performing one global refinement.
   * <li> Numbers of the kind $192\times 2^m$ with $m\geq 0$ integer. This
   *      choice is similar to the 24 and 48 cell cases, but provides
   *      additional refinements in azimuthal direction combined with a single
   *      layer in radial direction. The base mesh is either the 6 or 12 cell
   *      version, depending on whether $m$ in the power is odd or even,
   *      respectively.
   * </ul>
   * The versions with 24, 48, and $2^m 192$ cells are useful if the shell is
   * thin and the radial lengths should be made more similar to the
   * circumferential lengths.
   *
   * The 3d grids with 12 and 96 cells are plotted below:
   *
   * @image html hypershell3d-12.png
   * @image html hypershell3d-96.png
   *
   * @note This function is declared to exist for triangulations of all space
   * dimensions, but throws an error if called in 1d.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  hyper_shell(Triangulation<dim> &tria,
              const Point<dim> &  center,
              const double        inner_radius,
              const double        outer_radius,
              const unsigned int  n_cells  = 0,
              bool                colorize = false);

  /**
   * Produce an eccentric hyper-shell, the region between two spheres centered
   * on two distinct center points. One has to specify the <tt>inner_center</tt>
   * and <tt>outer_center</tt>, with given <tt>inner_radius</tt> and
   * <tt>outer_radius</tt>. The number <tt>n_cells</tt> indicates the number of
   * cells of the resulting triangulation, i.e., how many cells form the ring
   * (in 2d) or the shell (in 3d).
   *
   * By default, the outer boundary has the indicator 1 while the inner boundary
   * has id 0. In 3d, this applies to both the faces and the edges of these
   * boundaries.
   *
   * A SphericalManifold is attached to the outer boundary with an id of 1 while
   * another SphericalManifold is attached to the inner boundary with an id of
   * 0. A TransfiniteInterpolationManifold is attached to all other cells and
   * faces of the triangulation with an id of 2.
   *
   * Here, the number <tt>n_cells</tt> of elements has the same meaning as in
   * GridGenerator::hyper_shell.
   *
   * The grids with a 30% offset of the inner shell in the x direction, 12
   * initial cells and 3 levels of global refinement are plotted below:
   *
   * @image html eccentric_hyper_shell_2D.png
   * @image html eccentric_hyper_shell_3D.png
   *
   * @note Because it uses the definition of the hyper shell, this function is
   * declared to exist for triangulations of all space dimensions, but throws an
   * error if called in 1d.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  eccentric_hyper_shell(Triangulation<dim> &triangulation,
                        const Point<dim> &  inner_center,
                        const Point<dim> &  outer_center,
                        const double        inner_radius,
                        const double        outer_radius,
                        const unsigned int  n_cells);

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
   * If the number of initial cells @p n_cells is zero in 2d (as is the
   * default), then it is computed adaptively such that the resulting elements
   * have the least aspect ratio. The argument is ignored in 3d, where the
   * coarse mesh always has 5 cells.
   *
   * If colorize is set to <code>true</code>, the inner, outer, and the part
   * of the boundary where $x=0$, get indicator 0, 1, and 2,
   * respectively. Additionally, in 2d, the boundary indicator 3 is given to
   * the vertical edge below the x-axis. Otherwise, if colorize is set to
   * <code>false</code> all indicators are set to 0 (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * All manifold ids are set to zero, and a SphericalManifold is attached
   * to the triangulation.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  half_hyper_shell(Triangulation<dim> &tria,
                   const Point<dim> &  center,
                   const double        inner_radius,
                   const double        outer_radius,
                   const unsigned int  n_cells  = 0,
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
   * If @p colorize is set to <code>true</code>, the inner, outer, left, and
   * right boundary get indicator 0, 1, 2, and 3 in 2d,
   * respectively. Otherwise all indicators are set to 0. In 3d indicator 2 is
   * at the face $x=0$, 3 at $y=0$, 4 at $z=0$ (see
   * @ref GlossColorization "the glossary entry on colorization").
   *
   * All manifold ids are set to zero, and a SphericalManifold is attached
   * to the triangulation.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   */
  template <int dim>
  void
  quarter_hyper_shell(Triangulation<dim> &tria,
                      const Point<dim> &  center,
                      const double        inner_radius,
                      const double        outer_radius,
                      const unsigned int  n_cells  = 0,
                      const bool          colorize = false);

  /**
   * Produce a domain that is the space between two cylinders in 3d, with
   * given length, inner and outer radius and a given number of elements. The
   * cylinder shell is built around the $z$-axis with the two faces located
   * at $z = 0$ and $z = $ @p length.
   *
   * If @p n_radial_cells is zero (as is the
   * default), then it is computed adaptively such that the resulting elements
   * have the least aspect ratio. The same holds for @p n_axial_cells.
   *
   * @note Although this function is declared as a template, it does not make
   * sense in 1D and 2D. Also keep in mind that this object is rotated
   * and positioned differently than the one created by cylinder().
   *
   * All manifold ids are set to zero, and a CylindricalManifold is attached
   * to the triangulation.
   *
   * @pre The triangulation passed as argument needs to be empty when calling
   * this function.
   *
   * @image html cylinder_shell.png
   *
   * In this picture, a cylinder shell of length 2, inner radius 0.5, outer
   * radius 1 is shown. The default argument for n_radial_cells and
   * n_axial_cells are used and a single global refinement is carried out.
   *
   */
  template <int dim>
  void
  cylinder_shell(Triangulation<dim> &tria,
                 const double        length,
                 const double        inner_radius,
                 const double        outer_radius,
                 const unsigned int  n_radial_cells = 0,
                 const unsigned int  n_axial_cells  = 0);

  /**
   * Produce the volume or surface mesh of a torus. The axis of the torus is
   * the $y$-axis while the plane of the torus is the $x$-$z$ plane.
   *
   * If @p dim is 3, the mesh will be the volume of the torus, using a mesh
   * equivalent to the circle in the poloidal coordinates with 5 cells on the
   * cross section. This function attaches a TorusManifold to all boundary
   * faces which are marked with a manifold id of 1, a CylindricalManifold to
   * the interior cells and all their faces which are marked with a manifold
   * id of 2 (representing a flat state within the poloidal coordinates), and
   * a TransfiniteInterpolationManifold to the cells between the TorusManifold
   * on the surface and the ToroidalManifold in the center, with cells marked
   * with manifold id 0.
   *
   * An example for the case if @p dim is 3 with a cut through the domain at
   * $z=0$, 6 toroidal cells, $R=2$ and $r=0.5$ without any global refinement
   * is given here:
   *
   * @image html torus_manifold_ids.png
   *
   * In this picture, the light gray shade represents the manifold id 0 of the
   * transfinite interpolation, which is applied to smoothly add new points
   * between the toroidal shape on the domain boundary and the inner rim where
   * a cylindrical description around the y-axis is prescribed. The inner rim
   * with the manifold id 2 is shown in red shade.
   *
   * If @p dim is 2, the mesh will describe the surface of the torus and this
   * function attaches a TorusManifold to all cells and faces (which are
   * marked with a manifold id of 0).
   *
   * @param tria The triangulation to be filled.
   *
   * @param R The radius of the circle, which forms the middle line of the
   * torus containing the loop of cells. Must be greater than @p r.
   *
   * @param r The inner radius of the torus.
   *
   * @param n_cells_toroidal Optional argument to set the number of cell
   * layers in toroidal direction. The default is 6 cell layers.
   *
   * @param phi Optional argument to generate an open torus with angle
   * $0 < \varphi <= 2 \pi$. The default value is $2 \pi$,
   * in which case a closed torus is generated. If the torus is open,
   * the torus is cut at two planes perpendicular to the torus centerline.
   * The center of these two planes are located at $(x_1, y_1, z_1) = (R, 0, 0)$
   * and $(x_2, y_2, z_2) = (R \cos(\varphi), 0, R \sin(\varphi))$.
   *
   * @note Implemented for Triangulation<2,3> and Triangulation<3,3>.
   */
  template <int dim, int spacedim>
  void
  torus(Triangulation<dim, spacedim> &tria,
        const double                  R,
        const double                  r,
        const unsigned int            n_cells_toroidal = 6,
        const double                  phi              = 2.0 * numbers::PI);

  /**
   * This function produces a square in the <i>xy</i>-plane with a cylindrical
   * hole in the middle. The square and the circle are centered at the
   * origin. In 3d, this geometry is extruded in $z$ direction to the interval
   * $[0,L]$.
   *
   * The inner boundary has a manifold id of $0$ and a boundary id of
   * $6$. This function attaches a PolarManifold or CylindricalManifold to the
   * interior boundary in 2d and 3d respectively. The other faces have
   * boundary ids of $0, 1, 2, 3, 4$, or $5$ given in the standard order of
   * faces in 2d or 3d.
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
   * different faces (see
   * @ref GlossColorization "the glossary entry on colorization").
   * The colors are given in lexicographic ordering for the
   * flat faces (0 to 3 in 2d, 0 to 5 in 3d) plus the curved hole (4 in 2d,
   * and 6 in 3d). If @p colorize is set to false, then flat faces get the
   * number 0 and the hole gets number 1.
   */
  template <int dim>
  void
  hyper_cube_with_cylindrical_hole(Triangulation<dim> &triangulation,
                                   const double        inner_radius = .25,
                                   const double        outer_radius = .5,
                                   const double        L            = .5,
                                   const unsigned int  repetitions  = 1,
                                   const bool          colorize     = false);

  /**
   * Produce a grid consisting of concentric shells. The primary difference
   * between this function and GridGenerator::hyper_shell() is that this
   * function permits unevenly spaced (in the radial direction)
   * @ref GlossCoarseMesh "coarse level cells".
   *
   * The parameters @p center, @p inner_radius, and @p outer_radius behave in
   * the same way as the first three arguments to
   * GridGenerator::hyper_shell. @p n_shells gives the total number of shells
   * to use (i.e., the number of cells in the radial direction). The outer
   * radius of the $k$th shell is given by
   *
   * @f[
   *     r = r_{\text{inner}} + (r_\text{outer} - r_\text{inner})
   *     \frac{1 - \tanh(\text{skewness}(1 - k/\text{n\_shells}))}
   *          {\tanh(\text{skewness})}
   * @f]
   *
   * where @p skewness is a parameter controlling the shell spacing in the
   * radial direction: values of @p skewness close to zero correspond to even
   * spacing, while larger values of @p skewness (such as $2$ or $3$)
   * correspond to shells biased to the inner radius.
   *
   * @p n_cells_per_shell is the same as in GridGenerator::hyper_shell: in 2d
   * the default choice of zero will result in 8 cells per shell (and 12 in
   * 3d). The only valid values in 3d are 6 (the default), 12, and 96 cells:
   * see the documentation of GridGenerator::hyper_shell for more information.
   *
   * If @p colorize is <code>true</code> then the outer boundary of the merged
   * shells has a boundary id of $1$ and the inner boundary has a boundary id
   * of $0$.
   *
   * Example: The following code (see, e.g., step-10 for instructions on how
   * to visualize GNUPLOT output)
   *
   * @code
   * #include <deal.II/fe/mapping_q_generic.h>
   *
   * #include <deal.II/grid/grid_generator.h>
   * #include <deal.II/grid/grid_out.h>
   * #include <deal.II/grid/tria.h>
   *
   * #include <fstream>
   *
   * int main()
   * {
   *   using namespace dealii;
   *
   *   Triangulation<2> triangulation;
   *   GridGenerator::concentric_hyper_shells(triangulation,
   *                                          Point<2>(),
   *                                          1.0,
   *                                          2.0,
   *                                          5u,
   *                                          2.0);
   *
   *   GridOut grid_out;
   *   GridOutFlags::Gnuplot gnuplot_flags(false, 10, true);
   *   grid_out.set_flags(gnuplot_flags);
   *
   *   const MappingQGeneric<2> mapping(3);
   *   std::ofstream out("out.gpl");
   *   grid_out.write_gnuplot(triangulation, out, &mapping);
   * }
   * @endcode
   *
   * generates the following output:
   *
   * @image html concentric_hyper_shells_2d.svg
   *
   */
  template <int dim>
  void
  concentric_hyper_shells(Triangulation<dim> &triangulation,
                          const Point<dim> &  center,
                          const double        inner_radius      = 0.125,
                          const double        outer_radius      = 0.25,
                          const unsigned int  n_shells          = 1,
                          const double        skewness          = 0.1,
                          const unsigned int  n_cells_per_shell = 0,
                          const bool          colorize          = false);

  /**
   * Produce a ring of cells in 3d that is cut open, twisted and glued
   * together again. This results in a kind of moebius-loop.
   *
   * @param tria        The triangulation to be worked on.
   * @param n_cells     The number of cells in the loop. Must be greater than
   * 4.
   * @param n_rotations The number of rotations ($\pi/2$ each) to be performed
   * before gluing the loop together.
   * @param R           The radius of the circle, which forms the middle line
   * of the torus containing the loop of cells. Must be greater than @p r.
   * @param r           The radius of the cylinder bent together as a loop.
   */
  void moebius(Triangulation<3, 3> &tria,
               const unsigned int   n_cells,
               const unsigned int   n_rotations,
               const double         R,
               const double         r);

  /**
   * Call one of the other GridGenerator functions, parsing the name of the
   * function to call from the string @p grid_generator_function_name, and
   * the arguments to the function from the string
   * @p grid_generator_function_arguments.
   *
   * The string that supplies the arguments is passed to the function
   * Patterns::Tools::Convert<TupleTyple>::to_value(), where `TupleType` here is
   * a tuple containing **all** the arguments of the GridGenerator function,
   * including all optional arguments.
   *
   * An example usage of this function is given by:
   * @code
   * GridGenerator::generate_from_name_and_arguments(
   *   tria,
   *   "hyper_ball",
   *   "0.0, 0.0 : 1 : false");
   * @endcode
   * Here, the colon separates the function arguments, and the comma separates
   * the coordinates of a Point<2> argument.
   *
   * According to the arity of the `TupleType`, the arguments of the function
   * may be separated by different separators (see the documentation of
   * Patterns::Tuple for the details of how the conversion is
   * performed). If a wrong format is used, an exception is thrown, and the
   * expected format is output as an error message.
   *
   * All GridGenerator functions are supported. If you find some that are
   * missing, please open an issue on GitHub.
   *
   * @param tria                              The triangulation to be worked on
   * @param grid_generator_function_name      The name of the function to call
   * @param grid_generator_function_arguments The arguments of the function, in
   * the format of a tuple-convertible string
   *
   * @author Luca Heltai, 2019.
   */
  template <int dim, int spacedim>
  void
  generate_from_name_and_arguments(
    Triangulation<dim, spacedim> &tria,
    const std::string &           grid_generator_function_name,
    const std::string &           grid_generator_function_arguments);
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
   * functions exist to generate
   * @ref GlossCoarseMesh "coarse meshes". For example, the channel mesh
   * used in step-35 could in principle be created using a mesh created by the
   * GridGenerator::hyper_cube_with_cylindrical_hole function and several
   * rectangles, and merging them using the current function. The rectangles
   * will have to be translated to the right for this, a task that can be done
   * using the GridTools::shift function (other tools to transform individual
   * mesh building blocks are GridTools::transform, GridTools::rotate, and
   * GridTools::scale).
   *
   * Vertices that are less than @p duplicated_vertex_tolerance apart will be merged
   * together. It is usually necessary to set this value to something that
   * depends on the input triangulations in some way. One reasonable choice is
   * to use the minimum distance between all adjacent vertices of the input
   * mesh divided by some constant:
   *
   * @code
   * auto min_line_length = [](const Triangulation<dim> &tria) -> double
   * {
   *   double length = std::numeric_limits<double>::max();
   *   for (const auto &cell : tria.active_cell_iterators())
   *     for (unsigned int n = 0; n < GeometryInfo<dim>::lines_per_cell; ++n)
   *       length = std::min(length, (cell->line(n)->vertex(0) -
   *                                  cell->line(n)->vertex(1)).norm());
   *   return length;
   * };
   *
   * const double tolerance = std::min(min_line_length(triangulation_1),
   *                                   min_line_length(triangulation_2)) / 2.0;
   * @endcode
   *
   * This will merge any vertices that are closer than any pair of vertices on
   * the input meshes.
   *
   * @note The two input triangulations must be
   * @ref GlossCoarseMesh "coarse meshes", i.e., they can not have any
   * refined cells.
   *
   * @note The function copies the material ids of the cells of the two input
   * triangulations into the output triangulation. If @p copy_manifold_ids is
   * set to @p true, manifold ids will be copied. Boundary indicators are never
   * copied. In other words, if the two coarse meshes have anything but the
   * default boundary indicators, then you will have to set boundary indicators
   * again by hand in the output triangulation.
   *
   * @note This function does not attach any manifolds to @p result, nor does
   * it set any manifold ids. In particular, manifolds attached to the two
   * input triangulations will be lost in the @p result triangulation.
   *
   * @note For a related operation on refined meshes when both meshes are
   * derived from the same coarse mesh, see
   * GridGenerator::create_union_triangulation().
   */
  template <int dim, int spacedim>
  void
  merge_triangulations(const Triangulation<dim, spacedim> &triangulation_1,
                       const Triangulation<dim, spacedim> &triangulation_2,
                       Triangulation<dim, spacedim> &      result,
                       const double duplicated_vertex_tolerance = 1.0e-12,
                       const bool   copy_manifold_ids           = false);

  /**
   * Same as above but allows to merge more than two triangulations at once.
   * The following gives an example of how to use this function:
   * @code
   *   Triangulation<2> tria_1, tria_2, tria_3;
   *   // initialize tria_1, tria_2 and tria_3
   *   ...
   *   Triangulation<2> merged_triangulation;
   *   GridGenerator::merge_triangulations({&tria_1, &tria_2, &tria_3},
   *                                       merged_triangulation,
   *                                       1.0e-10,
   *                                       false);
   * @endcode
   */
  template <int dim, int spacedim>
  void
  merge_triangulations(
    const std::vector<const Triangulation<dim, spacedim> *> &triangulations,
    Triangulation<dim, spacedim> &                           result,
    const double duplicated_vertex_tolerance = 1.0e-12,
    const bool   copy_manifold_ids           = false);

  /**
   * \brief Replicate a given triangulation in multiple coordinate axes
   *
   * @param input The triangulation which will be replicated along the
   * coordinate axes.
   *
   * @param extents A vector with <tt>dim</tt> entries specifying how many
   * copies of a triangulation should be present along each coordinate axis.
   *
   * @param result The triangulation to be created. It needs to be empty upon
   * calling this function.
   *
   * This function creates a new Triangulation equal to a
   * <tt>dim</tt>-dimensional array of copies of @p input. Copies of @p input
   * are created by translating @p input along the coordinate axes. Boundary
   * ids of faces (but not lines in 3D) and all manifold ids are copied but
   * Manifold objects are not since most Manifold objects do not work
   * correctly when a Triangulation has been translated.
   *
   * To see how this works, consider the following code:
   * @code
   * Triangulation<2> input;
   * GridGenerator::hyper_cube_with_cylindrical_hole(input);
   * Triangulation<2> output;
   * GridGenerator::replicate_triangulation(input, {3, 2}, output);
   * @endcode
   * results in
   *
   * @image html replicated_tria_2d.png
   *
   * And, similarly, in 3D:
   * @code
   * Triangulation<3> input;
   * GridGenerator::hyper_cross(1, 1, 1, 2, 1, 2);
   * Triangulation<3> output;
   * GridGenerator::replicate_triangulation(input, {3, 2, 1}, output);
   * @endcode
   * results in
   *
   * @image html replicated_tria_3d.png
   *
   * @note This function determines the spacing of the copies of @p input
   * based on the BoundingBox of @p input. If the boundary faces of @p input
   * are not aligned with the coordinate axes then the copies might not share
   * common faces; i.e., this function is intended for simple geometries with
   * boundary faces aligned along the coordinate axes.
   */
  template <int dim, int spacedim = dim>
  void
  replicate_triangulation(const Triangulation<dim, spacedim> &input,
                          const std::vector<unsigned int> &   extents,
                          Triangulation<dim, spacedim> &      result);

  /**
   * Given the two triangulations specified as the first two arguments, create
   * the triangulation that contains the finest cells of both triangulation
   * and store it in the third parameter. Previous content of @p result will
   * be deleted.
   *
   * @note This function is intended to create an adaptively refined
   * triangulation that contains the <i>most refined cells</i> from two input
   * triangulations that were derived from the <i>same</i>
   * @ref GlossCoarseMesh "coarse mesh" by
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
   * @note If you want to create a mesh that is the merger of two other
   * @ref GlossCoarseMesh "coarse meshes",
   * for example in order to compose
   * a mesh for a complicated geometry
   * from meshes for simpler geometries, then this is not the function for
   * you. Instead, consider GridGenerator::merge_triangulations().
   *
   * @note This function assumes that both @p triangulation_1 and @p
   * triangulation_2 have the same manifold descriptions. The output
   * Triangulation @p has the same manifold ids as these two triangulations.
   *
   * @note Both of the source conditions need to be available entirely locally.
   * In other words, they can not be objects of type
   * parallel::distributed::Triangulation.
   */
  template <int dim, int spacedim>
  void
  create_union_triangulation(
    const Triangulation<dim, spacedim> &triangulation_1,
    const Triangulation<dim, spacedim> &triangulation_2,
    Triangulation<dim, spacedim> &      result);

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
   * @note Unlike most GridGenerator functions, this function does not attach
   * any manifolds to @p result, nor does it set any manifold ids.
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
  create_triangulation_with_removed_cells(
    const Triangulation<dim, spacedim> &input_triangulation,
    const std::set<typename Triangulation<dim, spacedim>::active_cell_iterator>
      &                           cells_to_remove,
    Triangulation<dim, spacedim> &result);

  /**
   * Extrude @p input in the $z$ direction from $z = 0$ to $z =
   * \text{height}$. The number of <em>slices</em>, or layers of cells
   * perpendicular to the $z = 0$ plane, will be @p n_slices slices (minimum is
   * 2). The boundary indicators of the faces of @p input will be assigned to
   * the corresponding side walls in $z$ direction. The bottom and top get the
   * next two free boundary indicators: i.e., if @p input has boundary ids of
   * $0$, $1$, and $42$, then the $z = 0$ boundary id of @p result will be $43$
   * and the $z = \text{height}$ boundary id will be $44$.
   *
   * This function does not, by default, copy manifold ids. The reason for
   * this is that there is no way to set the manifold ids on the lines of the
   * resulting Triangulation without more information: for example, if two
   * faces of @p input with different manifold ids meet at a shared vertex then
   * there is no <em>a priori</em> reason to pick one manifold id or another
   * for the lines created in @p result that are parallel to the $z$-axis and
   * pass through that point. If @p copy_manifold_ids is <code>true</code>
   * then this function sets line manifold ids by picking the one that appears
   * <em>first</em> in @p manifold_priorities. For example: if @p
   * manifold_priorities is <code>{0, 42, numbers::flat_manifold_id}</code>
   * and the line under consideration is adjacent to faces with manifold ids of
   * <code>0</code> and <code>42</code>, then that line will have a manifold id
   * of <code>0</code>. The correct ordering is almost always
   * <ol>
   *   <li>manifold ids set on the boundary,</li>
   *   <li>manifold ids that describe most of the cells in the Triangulation
   *   (e.g., numbers::flat_manifold_id), and</li>
   *   <li>any manifold ids corresponding to TransfiniteInterpolationManifold
   *   manifolds.</li>
   * </ol>
   *
   * In particular, since TransfiniteInterpolationManifold interpolates
   * between surrounding manifolds, its manifold id should usually not be set
   * on lines or faces that are adjacent to cells with different manifold
   * ids. The default value for @p manifold_priorities follows this ranking
   * (where each category is sorted in ascending order):
   * <ol>
   *   <li>manifold ids associated with manifolds that are not
   *   TransfiniteInterpolationManifold, and</li>
   *   <li>manifold ids associated with any TransfiniteInterpolationManifold
   *   objects.</li>
   * </ol>
   * Note that numbers::flat_manifold_id (should it be a manifold id of @p
   * input) will always be the last entry in the first category.
   *
   * @pre The 2d input triangulation @p input must be a
   * @ref GlossCoarseMesh "coarse mesh",
   * i.e., it cannot have any
   * refined cells.
   *
   * @note Since @p input and @p output have different spatial dimensions, no
   * manifold objects are copied by this function regardless of the value of
   * @p copy_manifold_ids.
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const unsigned int                     n_slices,
    const double                           height,
    Triangulation<3, 3> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});

  /**
   * Overload of extrude_triangulation() to allow dimension independent
   * code to compile. This function throws an error when called, as
   * extrude_triangulation() is only implemented to extrude a dim=2 to a dim=3
   * Triangulation.
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const unsigned int                     n_slices,
    const double                           height,
    Triangulation<2, 2> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});


  /**
   * Overload of the previous function. Take a 2d Triangulation that is being
   * extruded. Differing from the previous function taking height and number of
   * slices for uniform extrusion, this function takes z-axis values
   * @p slice_coordinates where the slicing will happen. The boundary indicators
   * of the faces of @p input are going to be assigned to the corresponding side
   * walls in z direction. The bottom and top get the next two free boundary
   * indicators.
   *
   * @pre The 2d input triangulation @p input must be a
   * @ref GlossCoarseMesh "coarse mesh",
   * i.e., it cannot have any
   * refined cells.
   *
   * @note Since @p input and @p output have different spatial dimensions no
   * manifold objects are copied (nor are any manifold ids set) by this
   * function.
   *
   * @author Weixiong Zheng, 2018
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const std::vector<double> &            slice_coordinates,
    Triangulation<3, 3> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});

  /**
   * Overload of extrude_triangulation() to allow dimension independent
   * code to compile. This function throws an error when called, as
   * extrude_triangulation() is only implemented to extrude a dim=2 to a dim=3
   * Triangulation.
   */
  void
  extrude_triangulation(
    const Triangulation<2, 2> &            input,
    const std::vector<double> &            slice_coordinates,
    Triangulation<2, 2> &                  result,
    const bool                             copy_manifold_ids   = false,
    const std::vector<types::manifold_id> &manifold_priorities = {});



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
   * @param[in] in_tria The base input for a new flat triangulation.
   * @param[out] out_tria The desired flattened triangulation constructed from
   * the in_tria.
   *
   * @note Since @p input and @p output have different spatial dimensions no
   * manifold objects are copied by this function: you must attach new
   * manifold objects to @p out_tria.
   *
   * @author Luca Heltai, 2014
   */
  template <int dim, int spacedim1, int spacedim2>
  void
  flatten_triangulation(const Triangulation<dim, spacedim1> &in_tria,
                        Triangulation<dim, spacedim2> &      out_tria);


  /**
   * Namespace Airfoil contains classes and functions in order to create a
   * C-type mesh for the (flow) field around Joukowski or NACA airfoils.
   *
   * @ingroup GridGenerator
   */
  namespace Airfoil
  {
    /**
     * AdditionalData collects all settings that are required to generate a
     * airfoil triangulation with the functions Airfoil::create_triangulation().
     */
    struct AdditionalData
    {
      /**
       * Type of the airfoil: either "NACA" or "Joukowksi" to choose airfoil
       * geometry among NACA and Joukowski airfoil.
       */
      std::string airfoil_type;

      /**
       * NACA serial number defining the airfoil shape.
       *
       * @note Currently serial numbers with length 4 are supported.
       * A good overview of NACA serial numbers is presented in Wikipedia
       * (https://en.wikipedia.org/wiki/NACA_airfoil)
       */
      std::string naca_id;

      /**
       * Center of Joukowski circle.
       *
       * @note A center on the x-axis leads to a symmetrical airfoil.
       */
      Point<2, double> joukowski_center;

      /**
       * Chord length of the airfoil, i.e. distance from leading to trailing
       * edge.
       */
      double airfoil_length;

      /**
       * Vertical distance from airfoil chord to upper boundary of the mesh
       * i.e. half of the total mesh height.
       */
      double height;

      /**
       * Length of mesh from the airfoil trailing edge to outflow boundary.
       */
      double length_b2;

      /**
       * Factor defining the inclination HG of the coarse grid
       * The figure shows the upper coarse grid with two different inclinations
       * - incline_factor = 0   --> face HG
       * - incline_factor = 0.5 --> face HG'
       * Coordinate of point G' is defined by incline_factor after interpolation
       * G'(0) = G(0) + incline_factor * (K(0) - G(0))
       * with incline_factor in [0,1).
       *
       *              o-----G---G'--K
       *           /  |     |  /    |
       *         /    o     | /     |
       *       /    /    \  |/      |
       *     o----o         H-------o
       */
      double incline_factor;

      /**
       * Factor to receive a finer mesh around the airfoil by increasing
       * bias_factor b.
       * Bias function: f(x) = tanh(bx) / tanh(x) with x in [0,1], leads to a
       * compression of values close to x = 1.
       */
      double bias_factor;

      /**
       * Number of global refinements.
       */
      unsigned int refinements;

      /**
       * Number of subdivisions along the airfoil in left block.
       */
      unsigned int n_subdivision_x_0;

      /**
       * Number of subdivisions along the airfoil in middle block.
       */
      unsigned int n_subdivision_x_1;

      /**
       * Number of subdivisions in block right of the airfoil.
       */
      unsigned int n_subdivision_x_2;

      /**
       * Number of subdivisions normal to the airfoil contour.
       */
      unsigned int n_subdivision_y;

      /**
       * Factor to enhance the approximation of the airfoil geometry that
       * happens when interpolating provided nonequidistant airfoil points to
       * equidistant airfoil points. When generating the required vector
       * consisting the equidistant airfoil points, it is interpolated between
       * nonequidistand airfoil points.
       * Increasing the provided nonequidistant airfoil points leads to
       * a better approximation of the airfoil geometry. Parameter
       * "airfoil_sampling_factor" thereby defines the relation of
       * provided_nonequidistant_points to required_equidistant_points.
       */
      unsigned int airfoil_sampling_factor;

      /**
       * Constructor.
       */
      AdditionalData();

      /**
       * This function adds the ParameterHandler entries.
       *
       * @param[in] prm Parameter handler.
       */
      void
      add_parameters(ParameterHandler &prm);
    };



    /**
     * Initialize the given triangulation with a flow field around an airfoil,
     * i.e., a mesh of C-Type approximating Joukowski or NACA (4 digit)
     * airfoils.
     *
     * The user can specify the airfoil geometry and the mesh setup by providing
     * input parameters for the struct AdditionalData.
     * Thereby, the user can choose among different types of Joukowski or NACA
     * airfoils with variable chord length, far field size and mesh density.
     *
     * @note This function creates a refined mesh (number of global refinements
     *       can be specified by the user). No manifold is attached. The
     *       vertices in the final mesh are moved by this function to the
     *       right position.
     *
     * @note This function is currently only implemented for 2D but the mesh
     *       can of course be extruded into the third dimension using
     *       GridGenerator::extrude().
     *
     * @param[out] tria The triangulation to be created. It needs to be empty
     *             upon calling this function.
     * @param[in] additional_data Configuration of the mesh.
     *
     \htmlonly <style>div.image
     *
     img[src="https://www.dealii.org/images/grids/airfoils_naca_joukowski.png"]{width:50%;}</style>
     \endhtmlonly
     * @image html https://www.dealii.org/images/grids/airfoils_naca_joukowski.png
     */
    template <int dim>
    void
    create_triangulation(
      Triangulation<dim, dim> &tria,
      const AdditionalData &   additional_data = AdditionalData());



    /**
     * The same as above but periodic boundary conditions on the
     * upper and lower faces of the far field are applied.
     *
     * @note This function is currently only implemented for 2D.
     *
     * @param[out] tria The triangulation to be created. It needs to be empty
     * upon calling this function.
     * @param[out] periodic_faces Periodic faces at upper and lower horizontal
     *                       boundaries.
     * @param[in] additional_data Configuration of the mesh.
     */
    template <int dim>
    void
    create_triangulation(
      Triangulation<dim, dim> &                            tria,
      std::vector<GridTools::PeriodicFacePair<
        typename Triangulation<dim, dim>::cell_iterator>> &periodic_faces,
      const AdditionalData &additional_data = AdditionalData());

  } // namespace Airfoil

  ///@}

  /**
   * @name Creating lower-dimensional meshes
   *
   * Created from parts of higher-dimensional meshes.
   */
  ///@{

#ifdef _MSC_VER
  // Microsoft's VC++ has a bug where it doesn't want to recognize that
  // an implementation (definition) of the extract_boundary_mesh function
  // matches a declaration. This can apparently only be avoided by
  // doing some contortion with the return type using the following
  // intermediate type. This is only used when using MS VC++ and uses
  // the direct way of doing it otherwise
  template <template <int, int> class MeshType, int dim, int spacedim>
  struct ExtractBoundaryMesh
  {
    using return_type =
      std::map<typename MeshType<dim - 1, spacedim>::cell_iterator,
               typename MeshType<dim, spacedim>::face_iterator>;
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
   * @note The function builds the surface mesh by creating a
   * @ref GlossCoarseMesh "coarse mesh"
   * from
   * the selected faces of the coarse cells of the volume mesh. It copies the
   * boundary indicators of these faces to the cells of the coarse surface
   * mesh. The surface mesh is then refined in the same way as the faces of
   * the volume mesh are. In order to ensure that the surface mesh has the
   * same vertices as the volume mesh, it is therefore important that you
   * assign appropriate boundary descriptions through
   * Triangulation::set_manifold() to the surface mesh object before calling
   * this function. If you don't, the refinement will happen under the
   * assumption that all faces are straight (i.e using the FlatManifold class)
   * rather than utilizing the Manifold object you may want to use to determine
   * the location of new vertices.
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
   * The order of vertices of surface cells and the corresponding
   * volume faces may not match in order to ensure that each surface cell is
   * associated with an outward facing normal.
   * As a consequence, if you want to match quantities on the faces of the
   * domain cells and on the cells of the surface mesh, you may have to
   * translate between vertex locations or quadrature points.
   *
   * @note The algorithm outlined above assumes that all faces on higher
   * refinement levels always have exactly the same boundary indicator as
   * their parent face. Consequently, we can start with coarse level faces and
   * build the surface mesh based on that. It would not be very difficult to
   * extend the function to also copy boundary indicators from finer level
   * faces to their corresponding surface mesh cells, for example to
   * accommodate different geometry descriptions in the case of curved
   * boundaries (but this is not currently implemented).
   *
   * @note Since @p volume_mesh and @p surface_mesh have different spatial
   * dimensions no manifold objects are copied by this function: you must
   * attach new manifold objects to @p surface_mesh.
   */
  template <template <int, int> class MeshType, int dim, int spacedim>
#ifndef _MSC_VER
  std::map<typename MeshType<dim - 1, spacedim>::cell_iterator,
           typename MeshType<dim, spacedim>::face_iterator>
#else
  typename ExtractBoundaryMesh<MeshType, dim, spacedim>::return_type
#endif
  extract_boundary_mesh(const MeshType<dim, spacedim> &     volume_mesh,
                        MeshType<dim - 1, spacedim> &       surface_mesh,
                        const std::set<types::boundary_id> &boundary_ids =
                          std::set<types::boundary_id>());

  ///@}


  /**
   * @name Exceptions
   */
  ///@{


  /**
   * Exception
   */
  DeclException0(ExcInvalidRadii);
  /**
   * Exception
   */
  DeclException1(ExcInvalidRepetitions,
                 int,
                 << "The number of repetitions " << arg1 << " must be >=1.");
  /**
   * Exception
   */
  DeclException1(ExcInvalidRepetitionsDimension,
                 int,
                 << "The vector of repetitions  must have " << arg1
                 << " elements.");

  /**
   * Exception for input that is not properly oriented.
   */
  DeclExceptionMsg(ExcInvalidInputOrientation,
                   "The input to this function is oriented in a way that will"
                   " cause all cells to have negative measure.");
  ///@}

#ifndef DOXYGEN
  // These functions are only implemented with specializations; declare them
  // here
  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<1> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<2> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void hyper_cube_with_cylindrical_hole(Triangulation<3> &,
                                        const double,
                                        const double,
                                        const double,
                                        const unsigned int,
                                        const bool);

  template <>
  void channel_with_cylinder(Triangulation<1> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);

  template <>
  void channel_with_cylinder(Triangulation<2> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);

  template <>
  void channel_with_cylinder(Triangulation<3> &,
                             const double,
                             const unsigned int,
                             const double,
                             const bool);
#endif
} // namespace GridGenerator



DEAL_II_NAMESPACE_CLOSE

#endif
