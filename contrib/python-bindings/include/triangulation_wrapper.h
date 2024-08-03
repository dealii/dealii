// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_triangulation_wrapper_h
#define dealii_triangulation_wrapper_h

#include <deal.II/base/config.h>

#include <boost/python.hpp>

#include <manifold_wrapper.h>
#include <mapping_wrapper.h>
#include <point_wrapper.h>
#include <quadrature_wrapper.h>

#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace python
{
  class CellAccessorWrapper;

  class TriangulationWrapper
  {
  public:
    typedef std::vector<CellAccessorWrapper>::iterator iterator;

    /**
     * Declare some symbolic names for mesh smoothing. These values are copied
     * from the Triangulation class, please find the documentation for their
     * meanings.
     */
    enum MeshSmoothing
    {
      none                               = 0x0,
      limit_level_difference_at_vertices = 0x1,
      eliminate_unrefined_islands        = 0x2,
      patch_level_1                      = 0x4,
      coarsest_level_1                   = 0x8,
      allow_anisotropic_smoothing        = 0x10,
      eliminate_refined_inner_islands    = 0x100,
      eliminate_refined_boundary_islands = 0x200,
      do_not_produce_unrefined_islands   = 0x400,
      smoothing_on_refinement =
        (limit_level_difference_at_vertices | eliminate_unrefined_islands),
      smoothing_on_coarsening =
        (eliminate_refined_inner_islands | eliminate_refined_boundary_islands |
         do_not_produce_unrefined_islands),
      maximum_smoothing = 0xffff ^ allow_anisotropic_smoothing
    };

    /**
     * Constructor. Takes a string @p dim with one of the following values
     * "2D", "2d", "3D", or "3d". The optional @p mesh_smoothing determines
     * the level of smoothness of the mesh size function that should be enforced
     * upon mesh refinement. The optional @p check_for_distorted_cells
     * determines whether the triangulation should check whether any of the
     * cells are distorted.
     */
    TriangulationWrapper(const std::string &dim,
                         const int          mesh_smoothing            = none,
                         const bool         check_for_distorted_cells = false);

    /**
     * Constructor. Takes a string @p dim with one of the following values
     * "2D", "2d", "3D", or "3d" and a string @p spacedim with one of the
     * following values "2D", "2d", "3D", or "3d". The dimension of @p spacedim
     * must be larger than the dimension of @p dim. The optional @p mesh_smoothing
     * determines the level of smoothness of the mesh size function that should
     * be enforced upon mesh refinement. The optional @p check_for_distorted_cells
     * determines whether the triangulation should check whether any of the
     * cells are distorted.
     */
    TriangulationWrapper(const std::string &dim,
                         const std::string &spacedim,
                         const int          mesh_smoothing            = none,
                         const bool         check_for_distorted_cells = false);

    /**
     * Destructor.
     */
    ~TriangulationWrapper();

    /**
     * Create a triangulation from a list of vertices and a list of indices,
     * each of the latter being a list of 1<<dim vertex indices.
     */
    void
    create_triangulation(const boost::python::list &vertices,
                         const boost::python::list &cells_vertices);

    /**
     * Return the number of active cells.
     */
    unsigned int
    n_active_cells() const;

    /**
     * Return the number of cells.
     */
    unsigned int
    n_cells() const;

    /*! @copydoc GridGenerator::hyper_cube
     */
    void
    generate_hyper_cube(const double left     = 0.,
                        const double right    = 1.,
                        const bool   colorize = false);

    /*! @copydoc GridGenerator::simplex
     */
    void
    generate_simplex(boost::python::list &vertices);

    /*! @copydoc GridGenerator::subdivided_hyper_cube
     */
    void
    generate_subdivided_hyper_cube(const unsigned int repetitions,
                                   const double       left  = 0.,
                                   const double       right = 1.);

    /*! @copydoc GridGenerator::hyper_rectangle
     */
    void
    generate_hyper_rectangle(PointWrapper &p1,
                             PointWrapper &p2,
                             const bool    colorize = false);

    /*! @copydoc GridGenerator::subdivided_hyper_rectangle
     */
    void
    generate_subdivided_hyper_rectangle(boost::python::list &repetitions,
                                        PointWrapper        &p1,
                                        PointWrapper        &p2,
                                        const bool           colorize = false);

    /**
     * Like the previous function. However, here the first argument does not
     * denote the number of subdivisions in each coordinate direction, but a
     * sequence of step sizes for each coordinate direction. This function is
     * therefore the right one to generate graded meshes where cells are
     * concentrated in certains areas, rather than a uniformly subdidived mesh
     * as the previous function generates.
     */
    void
    generate_subdivided_steps_hyper_rectangle(boost::python::list &step_sizes,
                                              PointWrapper        &p1,
                                              PointWrapper        &p2,
                                              const bool colorize = false);

    /**
     * Like the previous function, but with the following twist: the @p
     * material_id argument is a dim-dimensional array that, for each cell,
     * indicates which material_id should be set. In addition, and this is the
     * major new functionality, if the material_id of a cell is (-1), then that
     * cell is deleted from the triangulation, i.e. the domain will have a void
     * there.
     */
    void
    generate_subdivided_material_hyper_rectangle(
      boost::python::list &spacing,
      PointWrapper        &p,
      boost::python::list &material_id,
      const bool           colorize = false);

    /*! @copydoc GridGenerator::hyper_cube_with_cylindrical_hole
     */
    void
    generate_hyper_cube_with_cylindrical_hole(
      const double       inner_radius = .25,
      const double       outer_radius = .5,
      const double       L            = .5,
      const unsigned int repetitions  = 1,
      const bool         colorize     = false);

    /*! @copydoc GridGenerator::cheese
     */
    void
    generate_cheese(boost::python::list &holes);

    /*! @copydoc GridGenerator::plate_with_a_hole
     */
    void
    generate_plate_with_a_hole(const double        inner_radius = 0.4,
                               const double        outer_radius = 1.,
                               const double        pad_bottom   = 2.,
                               const double        pad_top      = 2.,
                               const double        pad_left     = 1.,
                               const double        pad_right    = 1.,
                               const PointWrapper &center = PointWrapper(),
                               const int           polar_manifold_id = 0,
                               const int           tfi_manifold_id   = 1,
                               const double        L                 = 1.,
                               const unsigned int  n_slices          = 2,
                               const bool          colorize          = false);

    /*! @copydoc GridGenerator::generate_channel_with_cylinder
     */
    void
    generate_channel_with_cylinder(const double       shell_region_width = 0.03,
                                   const unsigned int n_shells           = 2,
                                   const double       skewness           = 2.,
                                   const bool         colorize = false);

    /*! @copydoc GridGenerator::general_cell
     */
    void
    generate_general_cell(boost::python::list &vertices,
                          const bool           colorize = false);

    /*! @copydoc GridGenerator::parallelogram
     */
    void
    generate_parallelogram(boost::python::list &corners,
                           const bool           colorize = false);

    /*! @copydoc GridGenerator::parallelepiped
     */
    void
    generate_parallelepiped(boost::python::list &corners,
                            const bool           colorize = false);

    /*! @copydoc GridGenerator::subdivided_parallelepiped
     */
    void
    generate_fixed_subdivided_parallelepiped(const unsigned int n_subdivisions,
                                             boost::python::list &corners,
                                             const bool colorize = false);

    /**
     * A subdivided parallelepided, i.e., the same as above, but where the
     * number of subdivisions in each of the @tparam dim directions may vary.
     * Colorizing is done according to hyper_rectangle().
     */
    void
    generate_varying_subdivided_parallelepiped(
      boost::python::list &n_subdivisions,
      boost::python::list &corners,
      const bool           colorize = false);

    /*! @copydoc GridGenerator::enclosed_hyper_cube
     */
    void
    generate_enclosed_hyper_cube(const double left      = 0.,
                                 const double right     = 1.,
                                 const double thickness = 1.,
                                 const bool   colorize  = false);

    /*! @copydoc GridGenerator::hyper_ball
     */
    void
    generate_hyper_ball(PointWrapper &center, const double radius = 1.);

    /*! @copydoc GridGenerator::hyper_ball_balanced
     */
    void
    generate_hyper_ball_balanced(const PointWrapper &center = PointWrapper(),
                                 const double        radius = 1.);

    /*! @copydoc GridGenerator::hyper_sphere
     */
    void
    generate_hyper_sphere(PointWrapper &center, const double radius = 1.);

    /*! @copydoc GridGenerator::quarter_hyper_ball
     */
    void
    generate_quarter_hyper_ball(PointWrapper &center, const double radius = 1.);

    /*! @copydoc GridGenerator::half_hyper_ball
     */
    void
    generate_half_hyper_ball(PointWrapper &center, const double radius = 1.);

    /*! @copydoc GridGenerator::cylinder
     */
    void
    generate_cylinder(const double radius = 1., const double half_length = 1.);

    /*! @copydoc GridGenerator::subdivided_cylinder
     */
    void
    generate_subdivided_cylinder(const unsigned int x_subdivisions,
                                 const double       radius      = 1.,
                                 const double       half_length = 1.);

    /*! @copydoc GridGenerator::truncated_cone
     */
    void
    generate_truncated_cone(const double radius_0    = 1.,
                            const double radius_1    = 0.5,
                            const double half_length = 1.0);

    /*! @copydoc GridGenerator::hyper_shell
     */
    void
    generate_hyper_shell(PointWrapper  &center,
                         const double   inner_radius,
                         const double   outer_radius,
                         const unsigned n_cells  = 0,
                         bool           colorize = false);

    /*! @copydoc GridTools::shift
     */
    void
    shift(boost::python::list &shift_list);

    /*! @copydoc GridTools::scale
     */
    void
    scale(const double scaling_factor);

    /*! @copydoc GridGenerator::merge_triangulations
     */
    void
    merge_triangulations(boost::python::list &triangulations,
                         const double duplicated_vertex_tolerance = 1.0e-12,
                         const bool   copy_manifold_ids           = false);

    /*! @copydoc GridGenerator::replicate_triangulation
     */
    void
    replicate_triangulation(TriangulationWrapper &tria_in,
                            boost::python::list  &extents);

    /*! @copydoc GridGenerator::flatten_triangulation
     */
    void
    flatten_triangulation(TriangulationWrapper &tria_out);

    /*! @copydoc GridGenerator::extrude_triangulation
     */
    void
    extrude_triangulation(const unsigned int    n_slices,
                          const double          height,
                          TriangulationWrapper &tria_out);

    /*! @copydoc GridTools::distort_random
     */
    void
    distort_random(const double factor, const bool keep_boundary = true);

    /*! @copydoc GridTools::transform
     */
    void
    transform(boost::python::object &transformation);

    /*! @copydoc GridGenerator::convert_hypercube_to_simplex_mesh
     */
    void
    convert_hypercube_to_simplex_mesh(TriangulationWrapper &tria_out);

    /*! @copydoc GridTools::find_active_cell_around_point
     */
    CellAccessorWrapper
    find_active_cell_around_point(PointWrapper   &p,
                                  MappingQWrapper mapping = MappingQWrapper());

    /*! @copydoc GridTools::find_cells_adjacent_to_vertex
     */
    boost::python::list
    find_cells_adjacent_to_vertex(const unsigned int vertex_index);

    /*! @copydoc Triangulation::set_manifold
     */
    void
    set_manifold(const int number, ManifoldWrapper &manifold);

    /*! @copydoc Triangulation::reset_manifold
     */
    void
    reset_manifold(const int number);

    /*! @copydoc Triangulation::get_mesh_smoothing
     */
    int
    get_mesh_smoothing();

    /*! @copydoc Triangulation::set_mesh_smoothing
     */
    void
    set_mesh_smoothing(const int mesh_smoothing);

    /*! @copydoc Triangulation::refine_global
     */
    void
    refine_global(const unsigned int n);

    /*! @copydoc Triangulation::execute_coarsening_and_refinement
     */
    void
    execute_coarsening_and_refinement();

    /**
     * Return the list of active cell accessors associated to the underlying
     * Triangulation.
     */
    boost::python::list
    active_cells();

    /**
     * Return the list of cell accessors associated to the underlying
     * Triangulation.
     */
    boost::python::list
    cells();

    /*! @copydoc GridTools::minimal_cell_diameter
     */
    double
    minimal_cell_diameter() const;

    /*! @copydoc GridTools::maximal_cell_diameter
     */
    double
    maximal_cell_diameter() const;

    /*! @copydoc GridTools::compute_aspect_ratio_of_cells
     */
    boost::python::list
    compute_aspect_ratio_of_cells(const MappingQWrapper   &mapping,
                                  const QuadratureWrapper &quadrature);

    /**
     * Write mesh to the output file @filename according to the given data
     * format.
     */
    void
    write(const std::string &filename, const std::string format) const;

    /**
     * Read mesh from the file @filename using the given data
     * format.
     */
    void
    read(const std::string &filename, const std::string format) const;

    /**
     * Write the Triangulation in file.
     */
    void
    save(const std::string &filename) const;

    /**
     * Load the Triangulation from a file.
     */
    void
    load(const std::string &filename);

    /**
     * Return the dimension of the underlying Triangulation object.
     */
    int
    get_dim() const;

    /**
     * Return the space dimension of the underlying Triangulation object.
     */
    int
    get_spacedim() const;

    /**
     * Return a pointer that can be casted to the underlying Triangulation
     * object.
     */
    void *
    get_triangulation() const;

  private:
    /**
     * Helper function for the constructors.
     */
    void
    setup(const std::string &dimension,
          const std::string &spacedimension,
          const int          mesh_smoothing,
          const bool         check_for_distorted_cells);

    /**
     * Dimension of the underlying Triangulation object.
     */
    int dim;

    /**
     * Space dimension of the underlying Triangulation object.
     */
    int spacedim;

    /**
     * Pointer that can be casted to the underlying Triangulation object.
     */
    void *triangulation;
  };


  //-------------------- inline functions -----------------------//



  inline int
  TriangulationWrapper::get_dim() const
  {
    return dim;
  }



  inline int
  TriangulationWrapper::get_spacedim() const
  {
    return spacedim;
  }



  inline void *
  TriangulationWrapper::get_triangulation() const
  {
    return triangulation;
  }
} // namespace python

DEAL_II_NAMESPACE_CLOSE

#endif
