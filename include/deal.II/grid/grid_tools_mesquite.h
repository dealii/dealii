/* ---------------------------------------------------------------------
 *
 * Copyright (C) 2018 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * The deal.II library is free software; you can use it, redistribute
 * it, and/or modify it under the terms of the GNU Lesser General
 * Public License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * The full text of the license can be found in the file LICENSE at
 * the top level of the deal.II distribution.
 *
 * ---------------------------------------------------------------------
 */

#ifndef dealii_grid_tools_mesquite_h
#define dealii_grid_tools_mesquite_h


#include <deal.II/base/config.h>

#ifdef DEAL_II_WITH_MESQUITE

#  include <deal.II/base/exceptions.h>
#  include <deal.II/base/index_set.h>
#  include <deal.II/base/mpi.h>
#  include <deal.II/base/signaling_nan.h>

#  include <deal.II/distributed/shared_tria.h>
#  include <deal.II/distributed/tria.h>

#  include <deal.II/dofs/dof_handler.h>
#  include <deal.II/dofs/dof_tools.h>

#  include <deal.II/fe/fe_q.h>
#  include <deal.II/fe/fe_system.h>
#  include <deal.II/fe/mapping.h>
#  include <deal.II/fe/mapping_q1.h>
#  include <deal.II/fe/mapping_q_eulerian.h>

#  include <deal.II/grid/tria.h>

#  include <deal.II/hp/dof_handler.h>
#  include <deal.II/hp/mapping_collection.h>

#  include <deal.II/lac/block_vector.h>
#  include <deal.II/lac/vector.h>

#  include <deal.II/numerics/vector_tools.h>

#  include <boost/functional/hash.hpp>

#  include <Mesquite_all_headers.hpp>

#  ifdef DEAL_II_WITH_VERDICT
#    include <verdict.h>
#  endif

#  include <map>
#  include <memory>
#  include <utility>
#  include <vector>


DEAL_II_NAMESPACE_OPEN

namespace internal
{
  namespace Mesquite
  {
    // The meanings to the error codes can be found in the file
    // <path_to_mesquite_source>/src/Misc/Mesquite_MsqError.hpp
    // or at
    // https://github.com/trilinos/mesquite/blob/master/src/Misc/Mesquite_MsqError.hpp#L109
    std::string
    translate_error_code(const Mesquite2::MsqError::ErrorCode &code)
    {
      switch (code)
        {
          case (Mesquite2::MsqError::NO_ERROR):
            return "NO_ERROR: No error";
            break;
          case (Mesquite2::MsqError::UNKNOWN_ERROR):
            return "UNKNOWN_ERROR: Unknown error occured";
            break;
          case (Mesquite2::MsqError::OUT_OF_MEMORY):
            return "OUT_OF_MEMORY: Unable to allocate the necessary memory";
            break;
          case (Mesquite2::MsqError::INVALID_ARG):
            return "INVALID_ARG: Invalid function argument passed";
            break;
          case (Mesquite2::MsqError::NOT_INITIALIZED):
            return "NOT_INITIALIZED: Object not initialized";
            break;
          case (Mesquite2::MsqError::INVALID_STATE):
            return "INVALID_STATE: Object is in an invalid state";
            break;
          case (Mesquite2::MsqError::FILE_ACCESS):
            return "FILE_ACCESS: File cannot be opened/created";
            break;
          case (Mesquite2::MsqError::FILE_FORMAT):
            return "FILE_FORMAT: Wrong file type";
            break;
          case (Mesquite2::MsqError::PARSE_ERROR):
            return "PARSE_ERROR: Error parsing input (or input file)";
            break;
          case (Mesquite2::MsqError::IO_ERROR):
            return "IO_ERROR: An I/O error occured (e.g. read from file failed.)";
            break;
          case (Mesquite2::MsqError::INVALID_MESH):
            return "INVALID_MESH: The mesh is invalid";
            break;
          case (Mesquite2::MsqError::NO_PD_STORAGE_MODE):
            return "NO_PD_STORAGE_MODE: No storage mode chosen within PatchData";
            break;
          case (Mesquite2::MsqError::NOT_IMPLEMENTED):
            return "NOT_IMPLEMENTED: Requested functionality is not (yet) implemented";
            break;
          case (Mesquite2::MsqError::INTERNAL_ERROR):
            return "INTERNAL_ERROR: A bug in Mesquite";
            break;
          case (Mesquite2::MsqError::INTERRUPTED):
            return "INTERRUPTED: Application or user interrupted operation";
            break;
          case (Mesquite2::MsqError::TAG_ALREADY_EXISTS):
            return "TAG_ALREADY_EXISTS: Attempt to create tag that already exists";
            break;
          case (Mesquite2::MsqError::TAG_NOT_FOUND):
            return "TAG_NOT_FOUND: Specified tag does not exist";
            break;
          case (Mesquite2::MsqError::UNSUPPORTED_ELEMENT):
            return "UNSUPPORTED_ELEMENT: The element type is not supported.";
            break;
          case (Mesquite2::MsqError::PARALLEL_ERROR):
            return "PARALLEL_ERROR: An error occurred in parallel";
            break;
          case (Mesquite2::MsqError::BARRIER_VIOLATED):
            return "BARRIER_VIOLATED: Barrier violated when processing barrier Target Metric";
            break;
          default:
            AssertThrow(code < Mesquite2::MsqError::LAST_ERROR_CODE,
                        ExcMessage("Invalid error code"));
            break;
        }
      return "";
    }

    std::string
    translate_error_code(const Mesquite2::MsqError &err)
    {
      return translate_error_code(err.error_code());
    }
  } // namespace Mesquite
} // namespace internal


/**
 * An error message for Mesquite functions.
 *
 * @ingroup Exceptions
 */
DeclException1(ExcMesquiteError,
               Mesquite2::MsqError,
               << "An error with error number " << arg1.error_code()
               << " occurred while calling a Mesquite function."
               << "\n    " << internal::Mesquite::translate_error_code(arg1));


/**
 * An error message for Mesquite functions.
 *
 * @ingroup Exceptions
 */
DeclException2(ExcFENotEnoughComponents,
               unsigned int,
               unsigned int,
               << "The finite element has " << arg1
               << " components, but at least " << arg2 << " were expected.");


/**
 * An error message used to indicate that the triangulation used
 * to initialize a Mesquite mesh doesn't match the one to be
 * updated.
 *
 * @ingroup Exceptions
 */
DeclException0(ExcTriangulationMismatch);


namespace GridTools
{
  // Forward declaration
  namespace internal
  {
    template <int dim, typename RangeNumberType = double>
    class ReferencePosition;
  }

  /**
   * A function to transform from the vertex ordering used by deal.II to
   * a flattened vertex component map.
   */
  template <int dim>
  unsigned int
  vertex_map_index_from_global_vertex_index(const unsigned int vertex_index,
                                            const unsigned int component);

  /**
   * Returns the global vertex index and component for a given flattened
   * vertex component map index.
   * This is the reverse of vertex_map_index_from_global_vertex_index()
   */
  template <int dim>
  std::pair<unsigned int, unsigned int>
  global_vertex_index_component_from_vertex_map_index(
    const unsigned int vertex_map_index);


  /**
   * Get the position of triangulation vertices, as displaced by the
   * input @p mapping.
   */
  template <typename Number = double,
            int dim,
            int spacedim,
            template <int, int> class MeshType>
  Vector<Number>
  get_vertex_positions(const Mapping<dim, spacedim> & mapping,
                       const MeshType<dim, spacedim> &mesh);


  /**
   * Get the position of triangulation vertices, as displaced by the
   * input @p mapping.
   */
  template <typename Number = double,
            int dim,
            int spacedim,
            template <int, int> class MeshType>
  Vector<Number>
  get_vertex_positions(const hp::MappingCollection<dim, spacedim> &mapping,
                       const MeshType<dim, spacedim> &             mesh);


  /**
   * Get the indices of vertices that are on the boundary of the @p mesh.
   */
  template <int dim, int spacedim, template <int, int> class MeshType>
  IndexSet
  get_boundary_vertex_indices(const MeshType<dim, spacedim> &mesh);


  /**
   * Get the indices of vertices corresponding to hanging nodes of
   * the @p mesh.
   */
  template <int dim, int spacedim, template <int, int> class MeshType>
  IndexSet
  get_hanging_vertex_indices(const MeshType<dim, spacedim> &mesh);


  /**
   * Get the indices of vertices corresponding to hanging nodes of
   * the @p mesh.
   */
  template <int dim, int spacedim, template <int, int> class MeshType>
  std::vector<std::pair<unsigned int, unsigned int>>
  get_periodic_vertex_indices(const MeshType<dim, spacedim> &mesh);


  /**
   * @name Lagrangian (reference / material) configuration
   */
  //@{

  /**
   * Returns the Lagrangian vertex map for the given @p triangulation.
   * This, in effect, is simply a listing of the @p mesh coordinates in
   * a specific order.
   */
  template <int dim,
            int spacedim,
            typename Number = double,
            template <int, int> class MeshType>
  Vector<Number>
  get_lagrangian_vertex_positions(const MeshType<dim, spacedim> &mesh);

  //@}

  /**
   * @name Displacement fields (moving from the Lagrangian to Eulerian setting)
   */
  //@{


  /**
   * Returns a vertex displacement map for the given @p mesh based
   * on the specified grid @p transformation function. This vector maps vertex
   * coordinates from the reference triangulation to one in an Eulerian frame.
   */
  template <int dim,
            int spacedim,
            typename Number = double,
            template <int, int> class MeshType,
            template <int, typename> class Transformation>
  Vector<Number>
  get_vertex_displacement_map(
    const MeshType<dim, spacedim> &         mesh,
    const Transformation<spacedim, Number> &transformation,
    const unsigned int                      first_displacement_component);


  /**
   * Incrementatlly update the grid vertex positions using the
   * @p vertex_displacement_map.
   */
  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType>
  void
  displace_triangulation_vertices(
    MeshType<dim, spacedim> &mesh,
    const Vector<Number> &   vertex_displacement_map);

  //@}

  /**
   * @name Eulerian (current / spatial) configuration
   */
  //@{


  /**
   * Returns an Eulerian vertex map for the given @p mesh based on
   * the specified grid @p transformation function.
   */
  template <int dim,
            int spacedim,
            typename Number = double,
            template <int, int> class MeshType,
            template <int, typename> class Transformation>
  Vector<Number>
  get_eulerian_vertex_positions(
    const MeshType<dim, spacedim> &         mesh,
    const Transformation<spacedim, Number> &transformation,
    const unsigned int                      first_displacement_component);


  /**
   * Move the grid vertex positions from their original coordinates to
   *  a new set of coordinates using the @p eulerian_vertex_positions.
   */
  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType>
  void
  move_triangulation_vertices(MeshType<dim, spacedim> &mesh,
                              const Vector<Number> &eulerian_vertex_positions);

  //@}


  /* -------------- inline functions and template definitions -------------- */


  namespace internal
  {
    bool
    all_vertex_indices_touched(const std::vector<bool> &vertex_touched,
                               const unsigned int       n_local_vertices)
    {
      const unsigned int n_touched =
        std::accumulate(vertex_touched.begin(), vertex_touched.end(), 0);
      return n_touched == n_local_vertices;
    }
  } // namespace internal



#  ifndef DOXYGEN


  template <int dim>
  unsigned int
  vertex_map_index_from_global_vertex_index(const unsigned int vertex_index,
                                            const unsigned int component)
  {
    Assert(component < dim, ExcIndexRange(component, 0, dim));
    return dim * vertex_index + component;
  }


  template <int dim>
  std::pair<unsigned int, unsigned int>
  global_vertex_index_component_from_vertex_map_index(
    const unsigned int vertex_map_index)
  {
    const unsigned int component = vertex_map_index % dim;
    Assert(component < dim, ExcIndexRange(component, 0, dim));
    const unsigned int vertex_index = (vertex_map_index - component) / dim;
    return std::make_pair(vertex_index, component);
  }



  template <typename Number,
            int dim,
            int spacedim,
            template <int, int> class MeshType>
  Vector<Number>
  get_vertex_positions(const Mapping<dim, spacedim> & mapping,
                       const MeshType<dim, spacedim> &mesh)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const unsigned int n_vertices = tria.n_vertices();
    Vector<Number>     vertex_positions(dim * n_vertices);
    std::vector<bool>  vertex_touched(n_vertices, false);

    for (typename Triangulation<dim>::active_cell_iterator cell =
           tria.begin_active();
         cell != tria.end();
         ++cell)
      {
        const std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
          vertices = mapping.get_vertices(cell);

        // Extract the Eulerian vertex map from the mapped triangulation
        // vertices
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
#    ifdef DEBUG
            if (&mapping == &StaticMappingQ1<dim, spacedim>::mapping)
              {
                Assert((cell->vertex(v) - vertices[v]).norm() < 1e-12,
                       ExcMessage("Mapping inaccurate"));
              }
#    endif

            if (vertex_touched[cell->vertex_index(v)] == false)
              {
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    const unsigned int idx_vtx_map =
                      GridTools::vertex_map_index_from_global_vertex_index<dim>(
                        cell->vertex_index(v), d);
                    vertex_positions[idx_vtx_map] = vertices[v][d];
                  }

                vertex_touched[cell->vertex_index(v)] = true;
              }
          }
      }

    Assert(internal::all_vertex_indices_touched(vertex_touched, n_vertices),
           ExcMessage("Not all vertices have been touched."));

    return vertex_positions;
  }



  template <typename Number,
            int dim,
            int spacedim,
            template <int, int> class MeshType>
  Vector<Number>
  get_vertex_positions(const hp::MappingCollection<dim, spacedim> &mapping,
                       const MeshType<dim, spacedim> &             mesh)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    const unsigned int n_vertices = tria.n_vertices();
    Vector<Number>     vertex_positions(dim * n_vertices);
    std::vector<bool>  vertex_touched(n_vertices, false);

    for (typename Triangulation<dim>::active_cell_iterator cell =
           tria.begin_active();
         cell != tria.end();
         ++cell)
      {
        const std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell>
          vertices = mapping[cell->active_cell_index()].get_vertices(cell);

        // Extract the Eulerian vertex map from the mapped triangulation
        // vertices
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const unsigned int idx_vtx_map =
                    GridTools::vertex_map_index_from_global_vertex_index<dim>(
                      cell->vertex_index(v), d);
                  vertex_positions[idx_vtx_map] = vertices[v][d];
                }

              vertex_touched[cell->vertex_index(v)] = true;
            }
      }

    Assert(internal::all_vertex_indices_touched(vertex_touched, n_vertices),
           ExcMessage("Not all vertices have been touched."));

    return vertex_positions;
  }



  template <int dim, int spacedim, template <int, int> class MeshType>
  IndexSet
  get_boundary_vertex_indices(const MeshType<dim, spacedim> &mesh)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    IndexSet boundary_vertex_indices(tria.n_vertices());
    for (auto cell : tria.active_cell_iterators())
      {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (cell->face(f)->at_boundary() == false)
              continue;
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
                 ++v)
              {
                const unsigned int c = cell->face(f)->vertex_index(v);
                Assert(c < tria.n_vertices(),
                       ExcIndexRange(c, 0, tria.n_vertices()));
                boundary_vertex_indices.add_index(c);
              }
          }
      }

    Assert(boundary_vertex_indices.n_elements() > 0, ExcInternalError());
    return boundary_vertex_indices;
  }



  template <int dim, int spacedim, template <int, int> class MeshType>
  IndexSet
  get_hanging_vertex_indices(const MeshType<dim, spacedim> &mesh)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();
    if (!tria.has_hanging_nodes())
      return IndexSet();

    const unsigned int n_vertices = tria.n_vertices();
    IndexSet           hanging_vertex_indices(n_vertices);
    for (auto cell : tria.active_cell_iterators())
      {
        for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
          {
            if (cell->face(f)->at_boundary() == true)
              continue;
            if (cell->neighbor_is_coarser(f) == false)
              continue;
            const unsigned int fn = cell->neighbor_face_no(f);

            IndexSet vertices_1(n_vertices);
            IndexSet vertices_2(n_vertices);

            // We're now on a cell that has a coarser neighbour
            // Now find the common vertices; those that aren't
            // common are the hanging nodes
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
                 ++v)
              {
                const unsigned int c1 = cell->face(f)->vertex_index(v);
                const unsigned int c2 =
                  cell->neighbor(f)->face(fn)->vertex_index(v);
                Assert(c1 < tria.n_vertices(),
                       ExcIndexRange(c1, 0, tria.n_vertices()));
                Assert(c2 < tria.n_vertices(),
                       ExcIndexRange(c2, 0, tria.n_vertices()));
                vertices_1.add_index(c1);
                vertices_2.add_index(c2);
              }

            Assert(vertices_1.n_elements() ==
                     GeometryInfo<dim>::vertices_per_face,
                   ExcDimensionMismatch(vertices_1.n_elements(),
                                        GeometryInfo<dim>::vertices_per_face));
            Assert(vertices_2.n_elements() ==
                     GeometryInfo<dim>::vertices_per_face,
                   ExcDimensionMismatch(vertices_2.n_elements(),
                                        GeometryInfo<dim>::vertices_per_face));

            // Find the common vertices
            IndexSet vertices_diff(vertices_1);
            vertices_diff.subtract_set(vertices_2);
            hanging_vertex_indices.add_indices(vertices_diff);
          }
      }

    Assert(hanging_vertex_indices.n_elements() > 0, ExcInternalError());
    return hanging_vertex_indices;
  }



  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType>
  Vector<Number>
  get_lagrangian_vertex_positions(const MeshType<dim, spacedim> &mesh)
  {
#    ifdef DEBUG
    Vector<Number> diff =
      GridTools::get_vertex_positions(StaticMappingQ1<spacedim>::mapping, mesh);
    diff -=
      get_eulerian_vertex_positions(mesh,
                                    internal::ReferencePosition<spacedim>(),
                                    0 /*first_displacement_component*/);
    Assert(
      diff.l2_norm() < 1e-12,
      ExcMessage(
        "There is an inconsitency in the way Lagrangian vertex position maps are calculated."));
#    endif

    return GridTools::get_vertex_positions(StaticMappingQ1<spacedim>::mapping,
                                           mesh);
  }



  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType,
            template <int, typename> class Transformation>
  Vector<Number>
  get_vertex_displacement_map(
    const MeshType<dim, spacedim> &         mesh,
    const Transformation<spacedim, Number> &transformation,
    const unsigned int                      first_displacement_component)
  {
    Vector<Number> vertex_displacement_map(get_eulerian_vertex_positions(
      mesh, transformation, first_displacement_component));
    vertex_displacement_map -= get_lagrangian_vertex_positions(mesh);
    return vertex_displacement_map;
  }



  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType>
  void
  displace_triangulation_vertices(MeshType<dim, spacedim> &mesh,
                                  const Vector<Number> &vertex_displacement_map)
  {
    Triangulation<dim, spacedim> &triangulation = mesh.get_triangulation();
    Assert(vertex_displacement_map.size() == dim * triangulation.n_vertices(),
           ExcDimensionMismatch(vertex_displacement_map.size(),
                                dim * triangulation.n_vertices()));
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

    parallel::distributed::Triangulation<dim, spacedim> *pd_triangulation =
      (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
        &triangulation));
    const bool is_parallel_distributed = (pd_triangulation != nullptr);

    for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      {
        // Serial and parallel::shared Triangulations have a view of all of the
        // cells, while parallel::distributed does not.
        if (is_parallel_distributed && cell->is_locally_owned() == false)
          continue;

        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              Point<dim> vertex_displacement;
              for (unsigned int d = 0; d < dim; ++d)
                {
                  // Assume that the entries for the Eulerian
                  // vertex displacement map are sequential
                  const unsigned int idx_vtx_map =
                    GridTools::vertex_map_index_from_global_vertex_index<dim>(
                      cell->vertex_index(v), d);
                  vertex_displacement[d] = vertex_displacement_map(idx_vtx_map);
                }

              cell->vertex(v) += vertex_displacement;
              vertex_touched[cell->vertex_index(v)] = true;
            }
      }

    // Update the non-owned vertex locations
    if (is_parallel_distributed)
      pd_triangulation->communicate_locally_moved_vertices(vertex_touched);
  }



  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType,
            template <int, typename> class Transformation>
  Vector<Number>
  get_eulerian_vertex_positions(
    const MeshType<dim, spacedim> &         mesh,
    const Transformation<spacedim, Number> &transformation,
    const unsigned int                      first_displacement_component)
  {
    const unsigned int n_vertices = mesh.get_triangulation().n_vertices();
    Vector<Number>     eulerian_vertex_positions(spacedim * n_vertices);

    std::vector<bool> vertex_touched(n_vertices, false);
    for (typename MeshType<dim, spacedim>::active_cell_iterator cell =
           mesh.begin_active();
         cell != mesh.end();
         ++cell)
      {
        for (unsigned int v = 0; v < GeometryInfo<spacedim>::vertices_per_cell;
             ++v)
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              const Point<spacedim> &pt = cell->vertex(v);
              Vector<Number> new_vertex_position(transformation.n_components);
              transformation.vector_value(pt, new_vertex_position);

              for (unsigned int d = 0; d < spacedim; ++d)
                {
                  // Ensure that the entries for the Eulerian
                  // vertex displacement map are sequential
                  const unsigned int idx_vtx_map =
                    vertex_map_index_from_global_vertex_index<spacedim>(
                      cell->vertex_index(v), d);
                  eulerian_vertex_positions[idx_vtx_map] =
                    new_vertex_position(d + first_displacement_component);
                }

              vertex_touched[cell->vertex_index(v)] = true;
            }
      }

    return eulerian_vertex_positions;
  }



  template <int dim,
            int spacedim,
            typename Number,
            template <int, int> class MeshType>
  void
  move_triangulation_vertices(MeshType<dim, spacedim> &mesh,
                              const Vector<Number> &eulerian_vertex_positions)
  {
    Triangulation<dim, spacedim> &triangulation = mesh.get_triangulation();
    Assert(eulerian_vertex_positions.size() == dim * triangulation.n_vertices(),
           ExcDimensionMismatch(eulerian_vertex_positions.size(),
                                dim * triangulation.n_vertices()));
    std::vector<bool> vertex_touched(triangulation.n_vertices(), false);

    parallel::distributed::Triangulation<dim, spacedim> *pd_triangulation =
      (dynamic_cast<parallel::distributed::Triangulation<dim, spacedim> *>(
        &triangulation));
    const bool is_parallel_distributed = (pd_triangulation != nullptr);

    for (typename Triangulation<dim, spacedim>::active_cell_iterator cell =
           triangulation.begin_active();
         cell != triangulation.end();
         ++cell)
      {
        // Serial and parallel::shared Triangulations have a view of all of the
        // cells, while parallel::distributed does not.
        if (is_parallel_distributed && cell->is_locally_owned() == false)
          continue;

        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          if (vertex_touched[cell->vertex_index(v)] == false)
            {
              Point<dim> vertex_new_position;
              for (unsigned int d = 0; d < dim; ++d)
                {
                  // Assume that the entries for the Eulerian
                  // vertex displacement map are sequential
                  const unsigned int idx_vtx_map =
                    GridTools::vertex_map_index_from_global_vertex_index<dim>(
                      cell->vertex_index(v), d);
                  vertex_new_position[d] =
                    eulerian_vertex_positions(idx_vtx_map);
                }

              cell->vertex(v)                       = vertex_new_position;
              vertex_touched[cell->vertex_index(v)] = true;
            }
      }

    // Update the non-owned vertex locations
    if (is_parallel_distributed)
      pd_triangulation->communicate_locally_moved_vertices(vertex_touched);

#    ifdef DEBUG
    Vector<Number> diff(eulerian_vertex_positions);
    diff -= GridTools::get_lagrangian_vertex_positions(mesh);
    Assert(diff.l2_norm() < 1e-9,
           ExcMessage("Mesh movement was not successful"));
#    endif
  }


#  endif // DOXYGEN



  //==================



  namespace internal
  {
    // Mesqite expects counter-clockwise vertex ordering
    // (in 3d, on the bottom face first then the top face).
    // See the GeometryInfo documentation to understand how
    // these conversions between orderings are developed.
    template <int dim>
    struct MesquiteGeometryInfoHelper;

    template <>
    struct MesquiteGeometryInfoHelper<1>
    {
      static constexpr std::array<unsigned int,
                                  GeometryInfo<1>::vertices_per_cell>
        vertex_ordering{{0}};
    };

    template <>
    struct MesquiteGeometryInfoHelper<2>
    {
      static constexpr std::array<unsigned int,
                                  GeometryInfo<2>::vertices_per_cell>
        vertex_ordering{{0, 1, 3, 2}};
    };

    template <>
    struct MesquiteGeometryInfoHelper<3>
    {
      static constexpr std::array<unsigned int,
                                  GeometryInfo<3>::vertices_per_cell>
        vertex_ordering{{0, 1, 3, 2, 4, 5, 7, 6}};
    };


    /**
     * Returns the component index for the mesh motion DoF.
     */
    template <int dim, int spacedim>
    unsigned int
    get_mesh_motion_dof_index(const FiniteElement<dim, spacedim> &fe,
                              const unsigned int                  I,
                              const unsigned int first_mm_component_index)
    {
      Assert(I < fe.dofs_per_cell, ExcIndexRange(I, 0, fe.dofs_per_cell));
      if (fe.is_primitive(I))
        {
          const unsigned int component_I =
            fe.system_to_component_index(I).first;
          if (component_I >= first_mm_component_index &&
              component_I < first_mm_component_index + spacedim)
            return component_I;
        }
      else
        {
          for (unsigned int b = 0; b < fe.n_base_elements(); ++b)
            {
              try
                {
                  for (unsigned int d = 0; d < dim; ++d)
                    {
                      const unsigned int component_mm_I =
                        first_mm_component_index + d;
                      const unsigned int system_index =
                        fe.component_to_system_index(component_mm_I, b);
                      if (system_index == I)
                        return component_mm_I;
                    }
                }
              catch (...)
                {
                  // Component does not exist in that base element...
                }
            }
        }

      return std::numeric_limits<unsigned int>::max();
    }



    /**
     * Returns whether or not a local DoF (with index @p I) is one related to
     * the mesh motion field or not.
     */
    template <int dim, int spacedim>
    bool
    is_mesh_motion_dof_index(const FiniteElement<dim, spacedim> &fe,
                             const unsigned int                  I,
                             const unsigned int first_mm_component_index)
    {
      Assert(I < fe.dofs_per_cell, ExcIndexRange(I, 0, fe.dofs_per_cell));
      return (get_mesh_motion_dof_index(fe, I, first_mm_component_index) <
              fe.dofs_per_cell);
    }



    /**
     * A struct that will help determine to which type distributed vectors
     * convert to upon localization.
     */
    template <typename VectorType, typename = void>
    struct LocalVector;


    /**
     * A struct that will help determine to which type distributed vectors
     * convert to upon localization. This specialization is for
     * non-block vector types.
     */
    template <typename VectorType>
    struct LocalVector<
      VectorType,
      typename std::enable_if<!IsBlockVector<VectorType>::value>::type>
    {
      typedef Vector<typename VectorType::value_type> type;
    };


    /**
     * A struct that will help determine to which type distributed vectors
     * convert to upon localization. This specialization is for
     * block vector types.
     */
    template <typename VectorType>
    struct LocalVector<
      VectorType,
      typename std::enable_if<IsBlockVector<VectorType>::value>::type>
    {
      typedef BlockVector<typename VectorType::value_type> type;
    };



    /**
     * A intermidiary data structure to help translate between an displacement
     * fields discretized by an arbitrary finite element and that of the
     * triangulation vertices.
     *
     * @author Jean-Paul Pelteret, 2018
     */
    template <int dim, int spacedim = dim>
    struct MeshMotionData
    {
      template <template <int, int> class MeshType>
      MeshMotionData(const MeshType<dim, spacedim> &mesh)
        : degree(1)
        , fe(FE_Q<dim, spacedim>(degree), dim)
        , dof_handler(mesh.get_triangulation())
        , first_mm_dof(0)
      {
        dof_handler.distribute_dofs(fe);
        reset();
      }


      ~MeshMotionData()
      {
        cache_cell_to_local_u_fe.clear();
        vertex_touched.clear();
        dof_handler.clear();
      }


      void
      reset() const
      {
        vertex_touched = std::vector<bool>(n_vertices(), false);
      }



      unsigned int
      n_vertices() const
      {
        return dof_handler.get_triangulation().n_vertices();
      }



      unsigned int
      n_vertex_map_components() const
      {
        return dim * n_vertices();
      }



      template <typename Number = double>
      Vector<Number>
      extract_lagrangian_vertex_positions() const
      {
        return GridTools::get_lagrangian_vertex_positions(this->dof_handler);
      }


      /**
       * @name Manipulating vertex maps with an input DoFHandler
       */
      //@{

      template <template <int, int> class DoFHandlerType, class VectorType>
      Vector<typename VectorType::value_type>
      extract_vertex_displacement_map(
        const DoFHandlerType<dim, spacedim> &dof_handler_soln,
        const VectorType &                   solution_with_mesh_motion_field,
        const unsigned int                   first_mesh_motion_component) const
      {
        reset();
        const unsigned int n_vertices =
          dof_handler_soln.get_triangulation().n_vertices();
        Vector<typename VectorType::value_type> vertex_displacement_map(
          dim * n_vertices);

        // In the case of a distributed vector, we need to make a local copy
        // of the vector
        typedef
          typename internal::LocalVector<VectorType>::type LocalVectorType;
        const LocalVectorType local_solution_with_mesh_motion_field(
          solution_with_mesh_motion_field);

        const Triangulation<dim, spacedim> &tria =
          this->dof_handler.get_triangulation();
        Assert(&dof_handler_soln.get_triangulation() == &tria,
               ExcTriangulationMismatch());
        for (typename Triangulation<dim>::active_cell_iterator cell_tria =
               tria.begin_active();
             cell_tria != tria.end();
             ++cell_tria)
          {
            // When using an hp::DoFHanderl, the active cell iterator
            // type will be incompatible with that of the mesh motion
            // DoFHandler. So we use this conversion to get around that.
            typename DoFHandlerType<dim, spacedim>::active_cell_iterator
                               cell_soln(&tria,
                        cell_tria->level(),
                        cell_tria->index(),
                        &dof_handler_soln);
            const unsigned int cell_soln_active_fe_index =
              cell_soln->active_fe_index();

            // Get the equivalent cell for the mesh motion
            // DoF handler
            typename DoFHandler<dim, spacedim>::active_cell_iterator cell_mm(
              &tria,
              cell_tria->level(),
              cell_tria->index(),
              &this->dof_handler);

            // Check that the solution FE has the minimal number of required
            // components
            Assert(cell_soln->get_fe().n_components() >=
                     first_mesh_motion_component + dim,
                   ExcFENotEnoughComponents(cell_soln->get_fe().n_components(),
                                            first_mesh_motion_component + dim));

            // If the FE discretisation for the input solution exactly matches
            // that of the mesh motion field, then we can exploit it (namely the
            // continuous nature of the discretisation and the location of the
            // support points) bypass a lot of additional computations.
            if (fe_compatible_with_vertex_value_extraction(
                  cell_soln->get_fe(), first_mesh_motion_component))
              {
                // This offset helps deal with two thing:
                // 1. The input FE is a (complex) FESystem with some additional
                // components
                //    not related to the mesh motion field.
                // 2. The case where not all DoFs have support points at the
                // vertices.
                const unsigned int fe_soln_component_offset =
                  get_mesh_motion_vertex_dof_offset(
                    cell_soln->get_fe(), first_mesh_motion_component);
                Assert(fe_soln_component_offset <=
                         (cell_soln->get_fe().dofs_per_vertex - dim),
                       ExcIndexRange(fe_soln_component_offset,
                                     0,
                                     (cell_soln->get_fe().dofs_per_vertex -
                                      dim)));

                // Extract the Eulerian vertex map from the current cell
                for (unsigned int v = 0;
                     v < GeometryInfo<dim>::vertices_per_cell;
                     ++v)
                  {
                    Assert(cell_mm->vertex_index(v) < vertex_touched.size(),
                           ExcIndexRange(cell_mm->vertex_index(v),
                                         0,
                                         vertex_touched.size()));
                    if (this->vertex_touched[cell_mm->vertex_index(v)] == false)
                      {
                        Assert(cell_soln->vertex_index(v) <
                                 vertex_touched.size(),
                               ExcIndexRange(cell_soln->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        Assert(
                          this->vertex_touched[cell_soln->vertex_index(v)] ==
                            false,
                          ExcMessage(
                            "Vertex should not have already been touched."));
                        for (unsigned int d = 0; d < dim; ++d)
                          {
                            const unsigned int idx_vtx_map = GridTools::
                              vertex_map_index_from_global_vertex_index<dim>(
                                cell_soln->vertex_index(v), d);
                            const unsigned int idx_soln =
                              cell_soln->vertex_dof_index(
                                v,
                                d + fe_soln_component_offset,
                                cell_soln_active_fe_index);
                            vertex_displacement_map[idx_vtx_map] =
                              local_solution_with_mesh_motion_field(idx_soln);
                          }

                        Assert(cell_mm->vertex_index(v) < vertex_touched.size(),
                               ExcIndexRange(cell_mm->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        this->vertex_touched[cell_mm->vertex_index(v)] = true;
                        Assert(cell_soln->vertex_index(v) <
                                 vertex_touched.size(),
                               ExcIndexRange(cell_soln->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        Assert(
                          this->vertex_touched[cell_soln->vertex_index(v)] ==
                            true,
                          ExcMessage(
                            "Vertex should have already been touched."));
                      }
                  }
              }
            else
              {
                // Get the local DoF numbering for each DoFHandler
                const unsigned int n_dofs_per_cell_soln =
                  cell_soln->get_fe().dofs_per_cell;
                const unsigned int n_dofs_per_cell_mm =
                  cell_mm->get_fe().dofs_per_cell;
                std::vector<types::global_dof_index> dof_indices_mm(
                  n_dofs_per_cell_mm);
                cell_mm->get_dof_indices(dof_indices_mm);

                // Get the local solution from the input DoFHandler
                Vector<typename VectorType::value_type> local_dof_values_soln(
                  n_dofs_per_cell_soln);
                cell_soln->get_dof_values(local_solution_with_mesh_motion_field,
                                          local_dof_values_soln);

                // For safely, first filter out the values that have nothing to
                // do with the mesh motion field. This just ensures that we
                // don't make any mistakes when extracting the mesh motion
                // solution for the sub-fe in the next step.
                for (unsigned int I = 0; I < n_dofs_per_cell_soln; ++I)
                  if (!is_mesh_motion_dof_index(cell_soln->get_fe(),
                                                I,
                                                first_mesh_motion_component))
                    local_dof_values_soln[I] =
                      numbers::signaling_nan<typename VectorType::value_type>();

                // Get the finite element related only to the mesh-motion field
                // for the input DoFHandler
                const FESystem<dim, spacedim> &fe_u_soln =
                  get_local_view_of_displacement_fe(
                    cell_soln, first_mesh_motion_component);
                const unsigned int n_dofs_per_cell_u_soln =
                  fe_u_soln.dofs_per_cell;

                // Next, extract the local sub-solution from the complete local
                // solution
                Vector<typename VectorType::value_type> local_dof_values_u_soln(
                  n_dofs_per_cell_u_soln);
                for (unsigned int I = 0, sub_I = 0; I < n_dofs_per_cell_soln;
                     ++I)
                  {
                    if (is_mesh_motion_dof_index(cell_soln->get_fe(),
                                                 I,
                                                 first_mesh_motion_component))
                      {
                        Assert(sub_I < local_dof_values_u_soln.size(),
                               ExcIndexRange(sub_I,
                                             0,
                                             local_dof_values_u_soln.size()));
                        local_dof_values_u_soln[sub_I++] =
                          local_dof_values_soln[I];
                      }

                    if (I == n_dofs_per_cell_soln - 1)
                      Assert(sub_I == local_dof_values_u_soln.size(),
                             ExcInternalError());
                  }

                // Transfer the solution to the mesh motion DoFHandler
                Assert(
                  cell_mm->get_fe().n_components() == fe_u_soln.n_components(),
                  ExcMessage(
                    "Finite elements for the mesh motion field are expected "
                    "to thave the same number of components."));
                FullMatrix<typename VectorType::value_type>
                  interpolation_matrix(n_dofs_per_cell_mm,
                                       n_dofs_per_cell_u_soln);
                FETools::get_interpolation_matrix(fe_u_soln,
                                                  cell_mm->get_fe(),
                                                  interpolation_matrix);
                Vector<typename VectorType::value_type> local_dof_values_mm(
                  n_dofs_per_cell_mm);
                interpolation_matrix.vmult(local_dof_values_mm,
                                           local_dof_values_u_soln);

                // Now to transfer the localised solution to the output
                // vector... First, extract the Eulerian vertex map from the
                // current cell
                const std::map<unsigned int, unsigned int>
                  map_dof_index_mm_to_vtx_comp_index =
                    get_mesh_motion_dof_index_to_vertex_component_index(
                      cell_mm);
                Assert(map_dof_index_mm_to_vtx_comp_index.size() ==
                         dof_indices_mm.size(),
                       ExcDimensionMismatch(
                         map_dof_index_mm_to_vtx_comp_index.size(),
                         dof_indices_mm.size()));

                // Now copy the data
                std::map<unsigned int, std::set<unsigned int>>
                  vertex_index_components_touched;
                for (unsigned int I = 0; I < n_dofs_per_cell_mm; ++I)
                  {
                    const unsigned int idx_mm = dof_indices_mm[I];
                    Assert(map_dof_index_mm_to_vtx_comp_index.find(idx_mm) !=
                             map_dof_index_mm_to_vtx_comp_index.end(),
                           ExcMessage("Index not in map."));
                    const auto it =
                      map_dof_index_mm_to_vtx_comp_index.find(idx_mm);
                    Assert(it != map_dof_index_mm_to_vtx_comp_index.end(),
                           ExcInternalError());
                    const unsigned int idx_vtx_map = it->second;

                    // Take note of which components of each index have been
                    // touched
                    const std::pair<unsigned int, unsigned int> idx_vtx_comp =
                      GridTools::
                        global_vertex_index_component_from_vertex_map_index<
                          dim>(idx_vtx_map);
                    Assert(idx_vtx_comp.first < n_vertices,
                           ExcIndexRange(idx_vtx_comp.first, 0, n_vertices));
                    Assert(idx_vtx_comp.second < dim,
                           ExcIndexRange(idx_vtx_comp.second, 0, dim));
                    Assert(vertex_index_components_touched[idx_vtx_comp.first]
                               .find(idx_vtx_comp.second) ==
                             vertex_index_components_touched[idx_vtx_comp.first]
                               .end(),
                           ExcInternalError());
                    vertex_index_components_touched[idx_vtx_comp.first].insert(
                      idx_vtx_comp.second);
                  }

                // Mark the touched components in the global vector
                Assert(
                  vertex_index_components_touched.size() ==
                    GeometryInfo<dim>::vertices_per_cell,
                  ExcDimensionMismatch(vertex_index_components_touched.size(),
                                       GeometryInfo<dim>::vertices_per_cell));
                for (const auto entry : vertex_index_components_touched)
                  {
                    // Check that there are the right number of components
                    // registered...
                    Assert(entry.second.size() == dim,
                           ExcDimensionMismatch(entry.second.size(), dim));
                    // ... and that they are within the right range.
                    for (const auto comp : entry.second)
                      {
                        (void)comp;
                        Assert(comp < dim, ExcIndexRange(comp, 0, dim));
                      }

                    const unsigned int idx_vtx = entry.first;
                    Assert(idx_vtx < vertex_touched.size(),
                           ExcIndexRange(idx_vtx, 0, vertex_touched.size()));
                    this->vertex_touched[idx_vtx] = true;
                  }
              }
          }

        Assert(internal::all_vertex_indices_touched(this->vertex_touched,
                                                    n_vertices),
               ExcMessage("Not all vertices have been touched."));

        return vertex_displacement_map;
      }


      template <template <int, int> class DoFHandlerType,
                typename VectorType,
                typename Number>
      void
      insert_vertex_displacement_map(
        VectorType &                         solution_with_mesh_motion_field,
        const DoFHandlerType<dim, spacedim> &dof_handler_soln,
        const Vector<Number> &               vertex_displacement_map,
        const unsigned int                   first_mesh_motion_component) const
      {
        insert_vertex_displacement_map(solution_with_mesh_motion_field,
                                       StaticMappingQ1<spacedim>::mapping,
                                       dof_handler_soln,
                                       vertex_displacement_map,
                                       first_mesh_motion_component);
      }


      template <template <int, int> class MappingType,
                template <int, int> class DoFHandlerType,
                typename VectorType,
                typename Number>
      void
      insert_vertex_displacement_map(
        VectorType &                         solution_with_mesh_motion_field,
        const MappingType<dim, spacedim> &   mapping,
        const DoFHandlerType<dim, spacedim> &dof_handler_soln,
        const Vector<Number> &               vertex_displacement_map,
        const unsigned int                   first_mesh_motion_component) const
      {
        reset();
        const unsigned int n_vertices =
          dof_handler_soln.get_triangulation().n_vertices();
        (void)n_vertices;
        Assert(n_vertices == this->dof_handler.get_triangulation().n_vertices(),
               ExcDimensionMismatch(
                 n_vertices,
                 this->dof_handler.get_triangulation().n_vertices()));
        Assert(vertex_displacement_map.size() == dim * n_vertices,
               ExcDimensionMismatch(vertex_displacement_map.size(),
                                    dim * n_vertices));

        // In the case of a distributed vector, we need to make a local copy
        // of the vector
        //        typedef typename internal::LocalVector<VectorType>::type
        //        LocalVectorType; LocalVectorType
        //        local_solution_with_mesh_motion_field
        //        (solution_with_mesh_motion_field);
        //        local_solution_with_mesh_motion_field = 0.0;

        const Triangulation<dim, spacedim> &tria =
          this->dof_handler.get_triangulation();
        Assert(&dof_handler_soln.get_triangulation() == &tria,
               ExcTriangulationMismatch());
        for (typename Triangulation<dim>::active_cell_iterator cell_tria =
               tria.begin_active();
             cell_tria != tria.end();
             ++cell_tria)
          {
            // When using an hp::DoFHanderl, the active cell iterator
            // type will be incompatible with that of the mesh motion
            // DoFHandler. So we use this conversion to get around that.
            typename DoFHandlerType<dim, spacedim>::active_cell_iterator
                               cell_soln(&tria,
                        cell_tria->level(),
                        cell_tria->index(),
                        &dof_handler_soln);
            const unsigned int cell_soln_active_fe_index =
              cell_soln->active_fe_index();

            // Get the equivalent cell for the mesh motion
            // DoF handler
            typename DoFHandler<dim, spacedim>::active_cell_iterator cell_mm(
              &tria,
              cell_tria->level(),
              cell_tria->index(),
              &this->dof_handler);

            // Check that the solution FE has the minimal number of required
            // components
            Assert(cell_soln->get_fe().n_components() >=
                     first_mesh_motion_component + dim,
                   ExcFENotEnoughComponents(cell_soln->get_fe().n_components(),
                                            first_mesh_motion_component + dim));

            // If the FE discretisation for the input solution exactly matches
            // that of the mesh motion field, then we can exploit it (namely the
            // continuous nature of the discretisation and the location of the
            // support points) bypass a lot of additional computations.
            if (fe_compatible_with_vertex_value_extraction(
                  cell_soln->get_fe(), first_mesh_motion_component) &&
                cell_soln->get_fe().tensor_degree() == 1)
              {
                // This offset helps deal with two thing:
                // 1. The input FE is a (complex) FESystem with some additional
                // components
                //    not related to the mesh motion field.
                // 2. The case where not all DoFs have support points at the
                // vertices.
                const unsigned int fe_soln_component_offset =
                  get_mesh_motion_vertex_dof_offset(
                    cell_soln->get_fe(), first_mesh_motion_component);
                Assert(fe_soln_component_offset <=
                         (cell_soln->get_fe().dofs_per_vertex - dim),
                       ExcIndexRange(fe_soln_component_offset,
                                     0,
                                     (cell_soln->get_fe().dofs_per_vertex -
                                      dim)));

                // Extract the Eulerian vertex map from the current cell
                for (unsigned int v = 0;
                     v < GeometryInfo<dim>::vertices_per_cell;
                     ++v)
                  {
                    Assert(cell_mm->vertex_index(v) < vertex_touched.size(),
                           ExcIndexRange(cell_mm->vertex_index(v),
                                         0,
                                         vertex_touched.size()));
                    if (this->vertex_touched[cell_mm->vertex_index(v)] == false)
                      {
                        Assert(cell_soln->vertex_index(v) <
                                 vertex_touched.size(),
                               ExcIndexRange(cell_soln->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        Assert(
                          this->vertex_touched[cell_soln->vertex_index(v)] ==
                            false,
                          ExcMessage(
                            "Vertex should not have already been touched."));
                        for (unsigned int d = 0; d < dim; ++d)
                          {
                            const unsigned int idx_vtx_map = GridTools::
                              vertex_map_index_from_global_vertex_index<dim>(
                                cell_soln->vertex_index(v), d);
                            const unsigned int idx_soln =
                              cell_soln->vertex_dof_index(
                                v,
                                d + fe_soln_component_offset,
                                cell_soln_active_fe_index);
                            Assert(
                              idx_vtx_map < vertex_displacement_map.size(),
                              ExcIndexRange(idx_vtx_map,
                                            0,
                                            vertex_displacement_map.size()));
                            if (solution_with_mesh_motion_field
                                  .locally_owned_elements()
                                  .is_element(idx_soln))
                              solution_with_mesh_motion_field(idx_soln) =
                                vertex_displacement_map(idx_vtx_map);
                            //                            local_solution_with_mesh_motion_field(idx_soln)
                            //                            =
                            //                            vertex_displacement_map(idx_vtx_map);
                          }

                        Assert(cell_mm->vertex_index(v) < vertex_touched.size(),
                               ExcIndexRange(cell_mm->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        this->vertex_touched[cell_mm->vertex_index(v)] = true;
                        Assert(cell_soln->vertex_index(v) <
                                 vertex_touched.size(),
                               ExcIndexRange(cell_soln->vertex_index(v),
                                             0,
                                             vertex_touched.size()));
                        Assert(
                          this->vertex_touched[cell_soln->vertex_index(v)] ==
                            true,
                          ExcMessage(
                            "Vertex should have already been touched."));
                      }
                  }
              }
            else
              {
                // Get the local DoF numbering for each DoFHandler
                const unsigned int n_dofs_per_cell_soln =
                  cell_soln->get_fe().dofs_per_cell;
                const unsigned int n_dofs_per_cell_mm =
                  cell_mm->get_fe().dofs_per_cell;
                std::vector<types::global_dof_index> dof_indices_soln(
                  n_dofs_per_cell_soln);
                std::vector<types::global_dof_index> dof_indices_mm(
                  n_dofs_per_cell_mm);
                cell_soln->get_dof_indices(dof_indices_soln);
                cell_mm->get_dof_indices(dof_indices_mm);

                // First, extract the Eulerian vertex map from the current
                // cell...
                const std::map<unsigned int, unsigned int>
                  map_dof_index_mm_to_vtx_comp_index =
                    get_mesh_motion_dof_index_to_vertex_component_index(
                      cell_mm);
                Assert(map_dof_index_mm_to_vtx_comp_index.size() ==
                         dof_indices_mm.size(),
                       ExcDimensionMismatch(
                         map_dof_index_mm_to_vtx_comp_index.size(),
                         dof_indices_mm.size()));

                // ... and then get the local solution from the mesh motion
                // DoFHandler
                Vector<Number> local_dof_values_mm(n_dofs_per_cell_mm);
                std::map<unsigned int, std::set<unsigned int>>
                  vertex_index_components_touched;
                for (unsigned int I = 0; I < n_dofs_per_cell_mm; ++I)
                  {
                    const unsigned int idx_mm = dof_indices_mm[I];
                    const auto         it =
                      map_dof_index_mm_to_vtx_comp_index.find(idx_mm);
                    Assert(it != map_dof_index_mm_to_vtx_comp_index.end(),
                           ExcInternalError());
                    const unsigned int idx_vtx_map = it->second;
                    Assert(idx_vtx_map < vertex_displacement_map.size(),
                           ExcIndexRange(idx_vtx_map,
                                         0,
                                         vertex_displacement_map.size()));
                    local_dof_values_mm[I] =
                      vertex_displacement_map(idx_vtx_map);

                    // Take note of which components of each index have been
                    // touched
                    const std::pair<unsigned int, unsigned int> idx_vtx_comp =
                      GridTools::
                        global_vertex_index_component_from_vertex_map_index<
                          dim>(idx_vtx_map);
                    Assert(idx_vtx_comp.first < n_vertices,
                           ExcIndexRange(idx_vtx_comp.first, 0, n_vertices));
                    Assert(idx_vtx_comp.second < dim,
                           ExcIndexRange(idx_vtx_comp.second, 0, dim));
                    Assert(vertex_index_components_touched[idx_vtx_comp.first]
                               .find(idx_vtx_comp.second) ==
                             vertex_index_components_touched[idx_vtx_comp.first]
                               .end(),
                           ExcInternalError());
                    vertex_index_components_touched[idx_vtx_comp.first].insert(
                      idx_vtx_comp.second);
                  }

                // Mark the touched components in the global vector
                Assert(
                  vertex_index_components_touched.size() ==
                    GeometryInfo<dim>::vertices_per_cell,
                  ExcDimensionMismatch(vertex_index_components_touched.size(),
                                       GeometryInfo<dim>::vertices_per_cell));
                for (const auto entry : vertex_index_components_touched)
                  {
                    // Check that there are the right number of components
                    // registered...
                    Assert(entry.second.size() == dim,
                           ExcDimensionMismatch(entry.second.size(), dim));
                    // ... and that they are within the right range.
                    for (const auto comp : entry.second)
                      {
                        (void)comp;
                        Assert(comp < dim, ExcIndexRange(comp, 0, dim));
                      }

                    const unsigned int idx_vtx = entry.first;
                    Assert(idx_vtx < vertex_touched.size(),
                           ExcIndexRange(idx_vtx, 0, vertex_touched.size()));
                    this->vertex_touched[idx_vtx] = true;
                  }

                // Get the finite element related only to the mesh-motion field
                // for the input DoFHandler
                const FESystem<dim, spacedim> &fe_u_soln =
                  get_local_view_of_displacement_fe(
                    cell_soln, first_mesh_motion_component);
                const unsigned int n_dofs_per_cell_u_soln =
                  fe_u_soln.dofs_per_cell;

                // Transfer the solution from the mesh motion DoFHandler to
                // the local displacement component of the solution vector
                Assert(
                  cell_mm->get_fe().n_components() == fe_u_soln.n_components(),
                  ExcMessage(
                    "Finite elements for the mesh motion field are expected "
                    "to thave the same number of components."));
                FullMatrix<Number> interpolation_matrix(n_dofs_per_cell_u_soln,
                                                        n_dofs_per_cell_mm);
                FETools::get_interpolation_matrix(cell_mm->get_fe(),
                                                  fe_u_soln,
                                                  interpolation_matrix);
                Vector<Number> local_dof_values_u_soln(n_dofs_per_cell_u_soln);
                interpolation_matrix.vmult(local_dof_values_u_soln,
                                           local_dof_values_mm);

                // Next, transfer the local sub-solution to the final solution
                for (unsigned int I = 0, sub_I = 0; I < n_dofs_per_cell_soln;
                     ++I)
                  {
                    if (is_mesh_motion_dof_index(cell_soln->get_fe(),
                                                 I,
                                                 first_mesh_motion_component))
                      {
                        Assert(sub_I < local_dof_values_u_soln.size(),
                               ExcIndexRange(sub_I,
                                             0,
                                             local_dof_values_u_soln.size()));
                        // TODO: Optimize this, as we end up writing into the
                        // same
                        //       entry multiple times (from each cell sharing
                        //       DoFs with support points on a cell face).

                        //                        local_solution_with_mesh_motion_field[dof_indices_soln[I]]
                        //                        =
                        //                        local_dof_values_u_soln[sub_I++];
                        const unsigned int idx_soln = dof_indices_soln[I];
                        if (solution_with_mesh_motion_field
                              .locally_owned_elements()
                              .is_element(idx_soln))
                          solution_with_mesh_motion_field(idx_soln) =
                            local_dof_values_u_soln[sub_I];

                        // Increment the local component counter, irrespective
                        // of whether it was used or not.
                        sub_I++;
                      }

                    if (I == n_dofs_per_cell_soln - 1)
                      Assert(sub_I == local_dof_values_u_soln.size(),
                             ExcInternalError());
                  }
              }
          }

        Assert(internal::all_vertex_indices_touched(this->vertex_touched,
                                                    n_vertices),
               ExcMessage("Not all vertices have been touched."));
        //        solution_with_mesh_motion_field =
        //        local_solution_with_mesh_motion_field;
        solution_with_mesh_motion_field.compress(VectorOperation::insert);
      }



      template <template <int, int> class DoFHandlerType, typename VectorType>
      Vector<typename VectorType::value_type>
      extract_eulerian_vertex_positions(
        const DoFHandlerType<dim, spacedim> &dof_handler_soln,
        const VectorType &                   solution_with_mesh_motion_field,
        const unsigned int                   first_mesh_motion_component) const
      {
        Vector<typename VectorType::value_type> eulerian_vertex_positions(
          extract_lagrangian_vertex_positions(dof_handler_soln,
                                              first_mesh_motion_component));

        eulerian_vertex_positions +=
          extract_vertex_displacement_map(dof_handler_soln,
                                          solution_with_mesh_motion_field,
                                          first_mesh_motion_component);
        return eulerian_vertex_positions;
      }

      //@}

      /**
       * @name Creating vertex maps with an input Mapping
       */
      //@{

      template <typename Number = double, template <int, int> class MappingType>
      Vector<Number>
      extract_lagrangian_vertex_positions(
        const MappingType<dim, spacedim> &mapping) const
      {
        return GridTools::get_lagrangian_vertex_positions(this->dof_handler);
      }



      template <typename Number = double, template <int, int> class MappingType>
      Vector<Number>
      extract_vertex_displacement_map(
        const MappingType<dim, spacedim> &mapping) const
      {
        Vector<Number> eulerian_vertex_displacement(
          extract_eulerian_vertex_positions(mapping));
        eulerian_vertex_displacement -=
          extract_lagrangian_vertex_positions(mapping);
        return eulerian_vertex_displacement;
      }



      template <typename Number = double, template <int, int> class MappingType>
      Vector<Number>
      extract_eulerian_vertex_positions(
        const MappingType<dim, spacedim> &mapping) const
      {
        Vector<Number> eulerian_vertex_displacement(
          extract_eulerian_vertex_positions(mapping));
        eulerian_vertex_displacement -=
          extract_lagrangian_vertex_positions(mapping);
        return eulerian_vertex_displacement;
      }

      //@}

      /**
       * @name Building constraints
       */
      //@{

      void
      make_hanging_vertex_constraints(
        ConstraintMatrix &hanging_vertex_constraints)
      {
        // Build the standard set of hanging node constraints
        ConstraintMatrix hanging_node_constraints;
        DoFTools::make_hanging_node_constraints(dof_handler,
                                                hanging_node_constraints);
        hanging_node_constraints.close();

        // And them translate them into hanging vertex constraints
        transfer_dof_constraints_to_vertex_constraints(
          hanging_vertex_constraints, dof_handler, hanging_node_constraints);
      }

      void
      make_periodic_vertex_constraints(
        ConstraintMatrix &periodic_vertex_constraints)
      {
        // Build the standard set of periodic constraints
        // This corresponds to the matching of displacement of vertices
        ConstraintMatrix periodic_constraints;
        AssertThrow(false, ExcNotImplemented());
        // DoFTools::make_periodicity_constraints(...);
        periodic_constraints.close();

        // And them translate them into hanging vertex constraints
        transfer_dof_constraints_to_vertex_constraints(
          periodic_vertex_constraints, dof_handler, periodic_constraints);
      }

      //@}

    private:
      // A FE that describes the mesh motion field
      // Since we're dealing with the vertices of a triangulation,
      // connected by straight edges, using a first-order FE_Q
      // it fine for this purpose.
      const unsigned int            degree;
      const FESystem<dim, spacedim> fe;
      DoFHandler<dim, spacedim>     dof_handler;
      const unsigned int            first_mm_dof;
      mutable std::vector<bool>     vertex_touched;

      /**
       * A cache for the views of the inputted FE's to mesh motion fields
       */
      mutable std::map<std::string, FESystem<dim, spacedim>>
        cache_cell_to_local_u_fe;


      /**
       * Check if the FE describing the various components of the solution field
       * all store values (not coefficients) at their support points, and that a
       * subset of those support points correspond to the vertices.
       *
       * @note It is assumed that the same FiniteElement (namely, FE_Q) is used
       * for all components of the displacement field.
       */
      bool
      fe_compatible_with_vertex_value_extraction(
        const FiniteElement<dim, spacedim> &fe_soln,
        const unsigned int                  first_mesh_motion_component) const
      {
        //        return false; // temporarily disable...
        //        try
        //        {
        //          const FEValuesExtractors::Vector extractor_u
        //          (first_mesh_motion_component); const ComponentMask mask_u =
        //          fe_soln.component_mask(extractor_u); if(dynamic_cast<const
        //          FiniteElement<dim,spacedim>*>(&fe_soln.get_sub_fe(mask_u)))
        //          // This function is desperate to throw and error
        //              return true;
        //        }
        //        catch (...)
        //        {
        //          // User didn't wrap their displacement FEs as an FESystem
        //        }
        //
        //        return false;

        // TODO: We have no idea whether the sub-FE's are collected into
        // FESystems or whether they're independently arranged. So we have to
        // increment over all of them until be find the ones that we want. See
        // get_mesh_motion_vertex_dof_offset()

        // Check that all FEs for the mesh motion DoFs are FE_Q's.
        for (unsigned int d = first_mesh_motion_component;
             d < dim + first_mesh_motion_component;
             ++d)
          if (!dynamic_cast<const FE_Q<dim, spacedim> *const>(
                &(fe_soln.get_sub_fe(d, 1))))
            return false;

        return true;

        //        return (fe_soln.n_base_elements() ==
        //        this->fe.n_base_elements() &&
        //                 dynamic_cast<const FE_Q<dim,spacedim>*
        //                 const>(&(fe_soln.base_element(0))));
      }



      unsigned int
      get_mesh_motion_vertex_dof_offset(
        const FiniteElement<dim, spacedim> &fe_soln,
        const unsigned int                  first_mesh_motion_component) const
      {
        // We have no idea whether the sub-FE's are collected into FESystems
        // or whether they're independently arranged. So we have to increment
        // over all of them until be find the ones that we want..
        unsigned int offset = 0;
        for (unsigned int f = 0; f < first_mesh_motion_component; ++f)
          {
            const FiniteElement<dim, spacedim> &sub_fe =
              fe_soln.get_sub_fe(f, 1);
            offset += sub_fe.dofs_per_vertex;
          }

        return offset;
      }



      /**
       * Returns a map that describes how DoFs related to the mesh motion
       * DoFHandler are mapped to the vertex DoFs on a given cell.
       */
      template <typename ActiveCellIterator>
      std::map<unsigned int, unsigned int>
      get_mesh_motion_dof_index_to_vertex_component_index(
        const ActiveCellIterator &cell_mm) const
      {
        std::map<unsigned int, unsigned int> map_dof_index_mm_to_vtx_comp_index;
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            for (unsigned int d = 0; d < dim; ++d)
              {
                const unsigned int idx_vtx_map =
                  GridTools::vertex_map_index_from_global_vertex_index<dim>(
                    cell_mm->vertex_index(v), d);
                Assert(idx_vtx_map < n_vertex_map_components(),
                       ExcIndexRange(idx_vtx_map, 0, n_vertices()));
                const unsigned int idx_mm =
                  cell_mm->vertex_dof_index(v, d + this->first_mm_dof);
                Assert(map_dof_index_mm_to_vtx_comp_index.find(idx_mm) ==
                         map_dof_index_mm_to_vtx_comp_index.end(),
                       ExcMessage("Index already in map."));
                map_dof_index_mm_to_vtx_comp_index[idx_mm] = idx_vtx_map;
              }
          }
        return map_dof_index_mm_to_vtx_comp_index;
      }



      /**
       * Get the finite element related only to the mesh-motion field for the
       * input DoFHandler.
       *
       * It would be nice to do the following:
       *     const FEValuesExtractors::Vector
       * extractor_u_soln(first_mesh_motion_component); const ComponentMask
       * mask_u_soln = cell_soln->get_fe().component_mask(extractor_u_soln);
       *     const FiniteElement<dim,spacedim>& fe_u_soln =
       * cell_soln->get_fe().get_sub_fe(mask_u_soln); Unfortunately, this does
       * not work for arbitrary finite element collections (e.g. if the user
       * does not wrap up the FE for the displacement field within its own
       * FESystem). So, we create our own sub-FESystem that mimics the
       * discretisation selected by the user.
       */
      template <typename ActiveCellIterator>
      const FESystem<dim, spacedim> &
      get_local_view_of_displacement_fe(
        const ActiveCellIterator &cell,
        const unsigned int        first_mesh_motion_component) const
      {
        const FiniteElement<dim, spacedim> &fe = cell->get_fe();
        const std::string                   key_str =
          fe.get_name() + "_fmm_dof_" +
          Utilities::int_to_string(first_mesh_motion_component);

        const auto it = cache_cell_to_local_u_fe.find(key_str);
        if (it != cache_cell_to_local_u_fe.end())
          {
            return it->second;
          }
        else
          {
            const FEValuesExtractors::Vector extractor_u_soln(
              first_mesh_motion_component);
            const ComponentMask mask_u_soln =
              fe.component_mask(extractor_u_soln);

            std::map<unsigned int, std::set<unsigned int>> map_base_components;
            std::set<unsigned int>                         found_u_components;
            for (unsigned int I = 0; I < fe.n_dofs_per_cell(); ++I)
              {
                if (is_mesh_motion_dof_index(fe,
                                             I,
                                             first_mesh_motion_component))
                  {
                    //              const unsigned int component          =
                    //              fe.system_to_component_index(I).first;
                    const unsigned int component =
                      get_mesh_motion_dof_index(fe,
                                                I,
                                                first_mesh_motion_component);
                    const unsigned int base_element_index =
                      fe.system_to_base_index(I).first.first;
                    Assert(mask_u_soln[component],
                           ExcMessage("Incorrect base element selected"));
                    map_base_components[base_element_index].insert(component);

                    found_u_components.insert(component);
                  }

                // Exit early if we've queried all of the necessary components
                if (found_u_components.size() == dim)
                  break;
              }

#  ifdef DEBUG
            unsigned int check_n_components = 0;
            for (const auto entry : map_base_components)
              check_n_components += entry.second.size();
            Assert(check_n_components == dim,
                   ExcMessage("Incorrect number of elements extracted."));
#  endif

            const unsigned int n_elements = map_base_components.size();
            std::vector<const FiniteElement<dim, spacedim> *> fes;
            std::vector<unsigned int>                         multiplicities;
            fes.reserve(n_elements);
            multiplicities.reserve(n_elements);
            for (auto entry : map_base_components)
              {
                const unsigned int base_element_index = entry.first;
                const unsigned int multiplicity       = entry.second.size();
                fes.push_back(&(fe.base_element(base_element_index)));
                multiplicities.push_back(multiplicity);
              }
            Assert(fes.size() == multiplicities.size(), ExcInternalError());

            cache_cell_to_local_u_fe.insert(
              std::make_pair(key_str,
                             FESystem<dim, spacedim>(fes, multiplicities)));
            return get_local_view_of_displacement_fe(
              cell, first_mesh_motion_component);
          }
      }


      /**
       * @note This function assumes that the same constraints will be enforced
       * on each component of the @p dof_handler.
       */
      void
      transfer_dof_constraints_to_vertex_constraints(
        ConstraintMatrix &               vertex_constraints,
        const DoFHandler<dim, spacedim> &dof_handler,
        const ConstraintMatrix &         dof_constraints)
      {
        reset();
        const unsigned int n_vertices =
          dof_handler.get_triangulation().n_vertices();
        (void)n_vertices;

        // We can't query the equivalent vertex index directly from a DoF index
        // extracted from the constraint matrix, so we get this using the
        // triangulation instead.
        std::map<types::global_dof_index, types::global_vertex_index>
          map_dof_to_vertex_index;
        for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
               dof_handler.begin_active();
             cell != dof_handler.end();
             ++cell)
          {
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 ++v)
              {
                for (unsigned int d = 0; d < dim; ++d)
                  map_dof_to_vertex_index[cell->vertex_dof_index(v, d)] =
                    cell->vertex_index(v);
              }
          }

        // Transfer the DoF constraints to the vertex index based equivalent
        // Only get the constraints for one component, since we don't
        // differentiate between vertex components
        const unsigned int component = 0;
        for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
               dof_handler.begin_active();
             cell != dof_handler.end();
             ++cell)
          {
            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 ++v)
              {
                const unsigned int idx_vtx_slave = cell->vertex_index(v);
                if (vertex_touched[idx_vtx_slave] == false)
                  {
                    const unsigned int idx_dof_slave =
                      cell->vertex_dof_index(v, component);
                    if (dof_constraints.is_constrained(idx_dof_slave))
                      {
                        vertex_constraints.add_line(idx_vtx_slave);

                        typedef const typename ConstraintMatrix::
                          ConstraintLine::Entries * ConstraintEntriesType;
                        const ConstraintEntriesType constraint_entries =
                          dof_constraints.get_constraint_entries(idx_dof_slave);
                        for (auto constraint : *constraint_entries)
                          {
                            const types::global_dof_index idx_dof_master =
                              constraint.first;
                            const auto it_m =
                              map_dof_to_vertex_index.find(idx_dof_master);
                            Assert(
                              it_m != map_dof_to_vertex_index.end(),
                              ExcMessage(
                                "Global index for master vertex not found"));
                            const double       v_coeff_m = constraint.second;
                            const unsigned int idx_vtx_master = it_m->second;

                            vertex_constraints.add_entry(idx_vtx_slave,
                                                         idx_vtx_master,
                                                         v_coeff_m);
                          }
                      }

                    vertex_touched[idx_vtx_slave] = true;
                  }
              }
          }

        Assert(
          vertex_constraints.n_constraints() * dim ==
            dof_constraints.n_constraints(),
          ExcMessage(
            "Not all of the constraints have been found and transferred."));
        Assert(vertex_constraints.n_constraints() < n_vertices,
               ExcMessage("Too many constraints have been generated"));
      }
    };

  } // namespace internal



  /**
   * An enumeration that defines the available selection of mesh
   * improvement wrappers that are provided by Mesquite.
   */
  enum class MesquiteWrapperTypes
  {
    /**
     * A standard patch-based Laplace smoother that aims to produce
     * a smooth mesh.
     *
     * Further implementation details are as follows:
     * - Notes: This is a local patch relaxation-solver. A 'smart'
     *   Laplacian solver is also available in Mesquite, but it is
     *   not used in this wrapper.
     * - Limitations/assumptions: No invertibility guarantee.
     * - Input Termination Criterion: Stop after 10 global iterations.
     * - Hardwired Parameters: None
     * - Global/Local: Local Patch with Culling
     */
    laplace = 0,
    /**
     * Makes the shape of an element as close as possible to that of
     * the ideal/regular element shape.
     *
     * Further implementation details are as follows:
     * - Notes: The wrapper will use a non-barrier metric on meshes that
     *   contain inverted elements and will use a barrier metric if the
     *   mesh does not contain inverted elements.
     * - Limitations/assumptions: There is no guarantee that the wrapper
     *   will be able to successfully untangle a mesh that contains
     *   inverted elements.
     * - Input Termination Criterion: CPU time limit of 300 seconds or
     *   maximum absolute vertex movement of 10 percent of the minimum
     *   edge length.
     * - Metric: TMPQualityMetric (Shape/ShapeBarrier)
     * - Objective Function: Algebraic mean of quality metric values
     * - Solver: Conjugate Gradient
     * - Global/Local: Global
     */
    shape_improvement,
    /**
     * Untangle elements. Prioritizes untangling over element shape or
     * other mesh quality measures.
     *
     * Further implementation details are as follows:
     * - Notes: A second optimization to improve element quality after
     *   untangling is often necessary.
     * - Limitations/assumptions: There is no guarantee that the optimal
     *   mesh computed using this wrapper will, in fact, be untangled.
     * - Input Termination Criterion: CPU time limit (not used if input
     *   value is non-positive) or fraction of mean edge length (default
     *   is 0.005). It also terminates if all elements are untangled, such
     *   that it should not modify an input mesh with no inverted elements.
     * - Metric: TUntangleBeta or TUntangleMu(TSizeNB1) or
     * TUntangleMu(TShapeSizeNB1)
     * - Objective Function: Algebraic mean of quality metric values
     * - Solver: Steepest Descent
     * - Global/Local: Local with culling, optionally Jacobi
     */
    untangler,
    /**
     * Make all the edges in the mesh of equal length while moving toward
     * the ideal shape. Intended for explicit PDE codes whose time-step
     * limination is governed by the minimum edge-length in the mesh.
     *
     * Further implementation details are as follows:
     * - Notes: Based on Target-matrix paradigm.
     * - Limitations/assumptions: Initial mesh must be non-inverted.
     *   User does not want to preserve or create anisotropic elements.
     * - Input Termination Criterion: maximum iterations (default=50),
     *   maximum absolute vertex movement
     * - Hardwired Parameters: None
     * - Metric: Target2DShapeSizeBarrier or Target3DShapeSizeBarrier
     * - Tradeoff Coefficient: 1.0
     * - Objective Function: Linear Average over the Sample Points
     * - Solver: Trust Region
     * - Global/Local: Global
     */
    minimum_edge_length_improvement,
    /**
     * Make the shape of an element as close as possible to that of
     * the ideal element shape, while preserving, as much as possible,
     * the size of each element in the mesh. To be used on isotropic
     * initial meshes that are already size-adapted.
     *
     * Further implementation details are as follows:
     * - Notes: Based on Target-matrix Paradigm.
     * - Limitations/assumptions: Initial mesh must be non-inverted.
     *   User wants to preserve sizes of elements in initial mesh and
     *   does not want to preserve or create anisotropic elements.
     * - Input Termination Criterion: maximum iterations (default=50),
     *   maximum absolute vertex movement
     * - Hardwired Parameters: None
     * - Metric: Target2DShapeSizeBarrier or Target3DShapeSizeBarrier
     * - Tradeoff Coefficient: 1.0
     * - Objective Function: Linear Average over the Sample Points
     * - Solver: Trust Region
     * - Global/Local: Global
     */
    size_adapted_shape_improvement,
    /**
     * Uses the initial mesh on undeformed geometric domain to guide
     * optimization of mesh moved to deformed geometric domain.
     *
     * Further implementation details are as follows:
     * - Notes: Uses a non-barrier metric which means that the wrapper
     *   could potentially invert/tangle elements.
     * - Limitations/assumptions: Application responsible for explicit
     *   handling of mesh on geometric curves and points. Initial mesh
     *   before moving to deformed domain must be available.
     * - Input Termination Criterion: CPU time limit (not used if input
     *   value is non-positive) or fraction of mean edge length
     *   (default is 0.005).
     * - Metric: TMPQualityMetric (TShapeNB1 or TShapeSizeNB1 or
     *   TShapeSizeOrientNB1)
     * - Objective Function: Algebraic mean of quality metric values
     * - Solver: Steepest Descent
     * - Global/Local: Local with culling
     */
    deforming_domain
  };



  /**
   * A smoother that can be used to enforce a set of constraints on
   * the vertices of a given mesh. It can be used either as a preconditioner
   * or a master smoother in any Mesquite2::InstructionQueue.
   *
   * @author Jean-Paul Pelteret, 2018
   */
  class VertexConstraintSlaver : public Mesquite2::VertexSlaver
  {
  public:
    /**
     * Class constructor
     */
    template <int dim, int spacedim, template <int, int> class MeshType>
    VertexConstraintSlaver(const MeshType<dim, spacedim> &mesh,
                           const ConstraintMatrix &       vertex_constraints)
      : vertex_constraints(vertex_constraints)
      , n_vertices(mesh.get_triangulation().n_vertices())
    {
      // We cannot handle periodic constraints just yet...
      Assert(
        mesh.get_triangulation().get_periodic_face_map().size() == 0,
        ExcMessage(
          "Enforcement of periodic constraints have not been implemented yet."))
    }


    /**
     * Destructor
     */
    virtual ~VertexConstraintSlaver()
    {}


    virtual double
    loop_over_mesh(Mesquite2::MeshDomainAssoc *mesh_and_domain,
                   const Mesquite2::Settings * settings,
                   Mesquite2::MsqError &       err)
    {
      Mesquite2::Mesh *mesh = mesh_and_domain->get_mesh();

      if (settings->get_slaved_ho_node_mode() !=
          Mesquite2::Settings::SLAVE_CALCULATED)
        {
          MSQ_SETERR(err)
          ("Request to calculate higher-order node slaved status "
           "when Settings::get_get_slaved_ho_node_mode() "
           "!= SLAVE_CALCUALTED",
           Mesquite2::MsqError::INVALID_STATE);
          return 0.0;
        }

      // Extract all of the vertices
      std::vector<Mesquite2::Mesh::VertexHandle> vertices;
      mesh->get_all_vertices(vertices, err);
      Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      Assert(vertices.size() == n_vertices,
             ExcDimensionMismatch(vertices.size(), n_vertices));

      // Extract the data flag (?) for each vertex
      std::vector<unsigned char> bytes(n_vertices);
      mesh->vertices_get_byte(vertices.data(), bytes.data(), n_vertices, err);
      Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

      // Set the slave flag for constrained vertices
      for (unsigned int v = 0; v < n_vertices; ++v)
        {
          if (vertex_constraints.is_constrained(v) == false)
            continue;
          bytes[v] |= Mesquite2::MsqVertex::MSQ_DEPENDENT;
        }

      // Push the flag changes
      mesh->vertices_set_byte(vertices.data(), bytes.data(), n_vertices, err);
      Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

      return 0.0;
    }

    virtual std::string
    get_name() const
    {
      return "VertexConstraintSlaver";
    }

    virtual void
    initialize_queue(Mesquite2::MeshDomainAssoc *mesh_and_domain,
                     const Mesquite2::Settings * settings,
                     Mesquite2::MsqError &       err)
    {}


  private:
    /**
     * A set of constraints to be applied
     */
    const ConstraintMatrix vertex_constraints;

    /**
     * The total number of vertices in the mesh
     */
    const unsigned int n_vertices;
  };



  /**
   * A smoother that can be used to enforce a set of constraints on
   * the vertices of a given mesh. It can be used either as a preconditioner
   * or a master smoother in any Mesquite2::InstructionQueue.
   *
   * @author Jean-Paul Pelteret, 2018
   */
  class VertexConstraintInstruction : public Mesquite2::VertexMover
  {
    // See the various classes in
    // https://github.com/trilinos/mesquite/tree/master/src/QualityImprover/Relaxation
    // https://github.com/trilinos/mesquite/blob/master/src/Mesh/Mesquite_MeshDecorator.cpp
    // for salvation, I mean, inspiration... Once you've studied those,
    // then dig deep and muster up the courage to continue.
    // Make sure that your will to live is strong because the documentation on
    // how to do mesh adjustment is just miserable.
    // For example, if you think that following SlaveBoundaryVertices
    // (which is derived from Mesquite2::VertexSlaver) was the right
    // way to go to enforce constraints, then you're wrong...
    // Its simply too limited in is interface.

  public:
    /**
     * Class constructor
     */
    VertexConstraintInstruction(
      Mesquite2::Mesh &             mesquite_mesh,
      const ConstraintMatrix &      vertex_constraints,
      Mesquite2::ObjectiveFunction *objective_function = nullptr);

    /**
     * Class destructor
     */
    virtual ~VertexConstraintInstruction();

    /**
     * An name identifying this class
     */
    virtual std::string
    get_name() const;

    /**
     * Return a pointer to a store PatchSet
     */
    virtual Mesquite2::PatchSet *
    get_patch_set();

  protected:
    /**
     * Operations performed upon initialization of a class instance.
     */
    virtual void
    initialize(Mesquite2::PatchData &, Mesquite2::MsqError &);

    /**
     * Operations performed upon destruction of a class instance.
     */
    virtual void
    cleanup();

    /**
     * Operations performed before mesh optimization commenses.
     */
    virtual void
    initialize_mesh_iteration(Mesquite2::PatchData &, Mesquite2::MsqError &);

    /**
     * Operations performed top optmize the coordinates of a patch of vertices.
     */
    virtual void
    optimize_vertex_positions(Mesquite2::PatchData &, Mesquite2::MsqError &);

    /**
     * Operations performed after mesh optimization is completed.
     */
    virtual void
    terminate_mesh_iteration(Mesquite2::PatchData &patch_data,
                             Mesquite2::MsqError & error);

  private:
    /**
     * The mesh to be adjusted
     */
    Mesquite2::Mesh *const mesquite_mesh;

    /**
     * A set of constraints to be applied
     */
    const ConstraintMatrix vertex_constraints;

    /**
     * We will work on the global view of the mesh,
     * not a patch. So this does nothing, but we
     * need to pass this back in one of the pure
     * virtual functions.
     */
    Mesquite2::VertexPatches patchSet;
  };



  /**
   * Mesh adaption using the Mesquite mesh adaptation library.
   *
   * This class provides an interface between deal.II's Triangulation class and
   * DoFHandler classes, and the Mesquite2::Mesh class implemented by Mesquite.
   * It provides the an interface to the Mesquite solvers and optimizers, and is
   * suitable for smoothing on meshes with hanging nodes.
   *
   * More information on the Mesquite suite can be found in
   * @code{.bib}
   * @Manual{Knupp2013a,
   * title       = {Mesquite Mesh Quality Improvement Toolkit User's Guide},
   * author      = {Knupp, P. and Freitag-Diachin, L. and Tidwell, B.},
   * year        = {2013},
   * institution = {Sandia National Laboratories},
   * keywords    = {Trilinos, Mesh Quality},
   * url         =
   * {https://software.sandia.gov/mesquite/doc-2.99/users-guide.pdf},
   * }
   * @endcode
   *
   * @warning This class is not yet compatible with parallel::distributed::Triangulation.
   *
   * @warning This class is not yet compatible with codimension 1 triangulations.
   *
   * @author Jean-Paul Pelteret, 2018
   */
  template <int dim, int spacedim = dim>
  class MesquiteMeshInterface
  {
    static_assert(
      dim == spacedim,
      "This class has not been implemented for the codimension 1 case.");

  public:
    // Forward declaration
    struct Settings;

    /**
     * Type of Mesquite mesh used
     */
    typedef Mesquite2::Mesh BaseMeshType;

    /**
     * Type of Mesquite mesh used
     */
    typedef Mesquite2::MeshImpl SerialMeshType;

    /**
     * Type of Mesquite mesh used
     */
    typedef Mesquite2::ParallelMeshImpl ParallelMeshType;

    /**
     * Types used to represent coordinate entries.
     */
    typedef double CoordinateType;

    /**
     * Types used to represent element indices.
     */
    typedef int ElementIndexType;

    /**
     * Types used to mark constraints
     */
    typedef bool ConstraintType;

    /**
     * Type used for vertex displacement maps
     */
    typedef Vector<double> VectorDisplacementMapType;

    /**
     * Type used for vertex position maps
     */
    typedef Vector<double> VectorPositionMapType;

    /**
     * Spatial dimension used to store the underlying mesh. For the given
     * @p BaseMeshType, it is necessary that this value is fixed and always
     * equal to 3.
     */
    static constexpr unsigned int meshdim = 3;

    /*
     * Default constructor.
     *
     * The initialize() function needs to be called before the class instance
     * is usable.
     */
    MesquiteMeshInterface();

    /**
     * Class constructor.
     *
     * This constructor extracts the necessary information and, most importantly
     * the vertex positions, directly from the triangulation and stores it in
     * its internal data structures.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the exact (reference) state of the
     * triangulation.
     *
     * @note The instance of the class is immediately usable after this constructor
     * is called.
     *
     * @tparam mesh A variable of a type that satisfies the requirements of the
     *         @ref ConceptMeshType "MeshType concept".
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     * @param additional_fixed_vertices An IndexSet that represents vertices
     *        that are to be considered fixed. If the @p fix_all_boundary_vertices
     *        flag is set then the boundary vertices will already be taken care
     * of.
     */
    template <template <int, int> class MeshType>
    MesquiteMeshInterface(
      const MeshType<dim, spacedim> &mesh,
      const bool                     fix_all_boundary_vertices = true,
      const IndexSet                 additional_fixed_vertices = IndexSet());

    /**
     * Class constructor.
     *
     * This constructor extracts the necessary information and, most importantly
     * the vertex positions, from the underlying triangulation as well as the
     * input displacement map, and stores it in its internal data structures
     * for later manipulation.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the Eulerian (spatial) view of the
     * triangulation, as determined by vector of
     * @p eulerian_vertex_displacements . If the Mapping does not modify the
     * position of the mesh vertices (like, for example, MappingQEulerian does),
     * then this function is equivalent to the one with the same name, and
     * without the mapping argument.
     *
     * @note The instance of the class is immediately usable after this constructor
     * is called.
     *
     * @tparam mesh A variable of a type that satisfies the requirements of the
     *         @ref ConceptMeshType "MeshType concept".
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     * @param additional_fixed_vertices An IndexSet that represents vertices
     *        that are to be considered fixed. If the @p fix_all_boundary_vertices
     *        flag is set then the boundary vertices will already be taken care
     * of.
     */
    template <template <int, int> class MappingType,
              template <int, int> class MeshType>
    MesquiteMeshInterface(
      const MappingType<dim, spacedim> &mapping,
      const MeshType<dim, spacedim> &   mesh,
      const bool                        fix_all_boundary_vertices = true,
      const IndexSet                    additional_fixed_vertices = IndexSet());

    /**
     * Class constructor.
     *
     * This constructor extracts the necessary information and, most importantly
     * the vertex positions, from the underlying triangulation as well as the
     * input vertex displacement vector, and stores it in its internal data
     * structures for later manipulation.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the Eulerian (spatial) view of the
     * triangulation, as determined by vector of
     * @p eulerian_vertex_displacements .
     *
     * @note The instance of the class is immediately usable after this constructor
     * is called.
     *
     * @note It is assumed that the entries for the Eulerian  vertex displacement
     *       map @p eulerian_vertex_positions are sequential, i.e. with the layout (for 3d)
     *       $\left[ v^{0}_{x} v^{0}_{y} v^{0}_{z},
     *               v^{1}_{x} v^{1}_{y} v^{1}_{z},
     *               ... ,
     *               v^{n}_{x} v^{n}_{y} v^{n}_{z} \right]$
     *       where $n$ is the number of vertices in the @p mesh.
     *
     * @tparam mesh A variable of a type that satisfies the requirements of the
     *         @ref ConceptMeshType "MeshType concept".
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     * @param additional_fixed_vertices An IndexSet that represents vertices
     *        that are to be considered fixed. If the @p fix_all_boundary_vertices
     *        flag is set then the boundary vertices will already be taken care
     * of.
     */
    template <template <int, int> class MeshType, typename Number>
    MesquiteMeshInterface(
      const MeshType<dim, spacedim> &mesh,
      const Vector<Number> &         eulerian_vertex_positions,
      const bool                     fix_all_boundary_vertices = true,
      const IndexSet                 additional_fixed_vertices = IndexSet());


    /**
     *
     * For some common cases, the "mesh motion field" would typically correspond
     * to the following:
     * - Finite elasticity in a Lagrangian framework: The total displacement
     * field
     * - Finite elasticity in a Eulerian framework: An incremental update of the
     *   displacement field
     * - ALE formulations: The mesh motion field (with a Lagrangian or Eularian
     *   flavor based on the implementation)
     */
    template <template <int, int> class DoFHandlerType, typename VectorType>
    MesquiteMeshInterface(
      const DoFHandlerType<dim, spacedim> &dof_handler,
      const VectorType &                   solution_with_mesh_motion_field,
      const unsigned int                   first_mesh_motion_component,
      const bool                           fix_all_boundary_vertices = true,
      const IndexSet additional_fixed_vertices = IndexSet());

    /**
     * Class destructor
     */
    virtual ~MesquiteMeshInterface();


    /**
     * Clear all data stored in this class.
     */
    void
    clear();


    /**
     * Returns the generated Mesquite view of a mesh.
     */
    BaseMeshType &
    get_mesquite_mesh();

    /**
     * Returns the generated Mesquite view of a mesh.
     */
    const BaseMeshType &
    get_mesquite_mesh() const;

    /**
     * Outputs a VTK file of the Mesquite mesh to a file with the
     * name @p filename.
     */
    void
    write_vtk(const std::string &filename) const;

    /**
     * Add to fixed vertex constraints to be applied to the mesh.
     *
     * @note This function can only be called after this class has been
     *       initialized, either by choosing the appropriate constructor or
     *       through a call to initialize() .
     */
    void
    add_mesh_fixed_vertices(const IndexSet &fixed_vertices);

    /**
     * Provide the set of fixed vertex constraints to be applied to the mesh.
     *
     * @note This function can only be called after this class has been
     *       initialized, either by choosing the appropriate constructor or
     *       through a call to initialize() .
     */
    void
    set_mesh_fixed_vertices(const IndexSet &fixed_vertices);

    /**
     * Use the built-in shape improvement algorithms to increase the quality of
     * the mesh. The result is stored within this class instance and is not
     * applied elsewhere until requested. This allows one to make multiple calls
     * to
     * @p execute and adjust the mesh several times before the result is used.
     *
     * @param max_vertex_movement Termination optimization if no vertex is moved
     *        by more than this distance in the previous solver step. This
     * parameter is only used by some of the algorithms.
     */
    void
    execute(const enum MesquiteWrapperTypes type,
            const Settings &                settings = Settings());


    /**
     * Use the built-in shape improvement algorithms to increase the quality of
     * the mesh. The result is stored within this class instance and is not
     * applied elsewhere until requested. This allows one to make multiple calls
     * to
     * @p execute and adjust the mesh several times before the result is used.
     *
     * @param max_vertex_movement Termination optimization if no vertex is moved
     *        by more than this distance in the previous solver step. This
     * parameter is only used by some of the algorithms.
     */
    void
    execute(const enum MesquiteWrapperTypes type,
            const MPI_Comm                  mpi_communicator,
            const Settings &                settings = Settings());

    /**
     * Use a custom shape improvement algorithm to increase the quality of the
     * mesh. The result is stored within this class instance and is not applied
     * elsewhere until requested. This allows one to make multiple calls to
     * @p execute and adjust the mesh several times before the result is used.
     *
     * @note The typical input type for this generic function is a
     *       Mesquite2::InstructionQueue, although Mesquite::Wrapper
     *       is also supported.
     *
     * @note By default, Mesquite will print the results of quality assessment
     *       to screen. If you would like to prevent this, add this line of code
     *       when you initialize the quality assessor:
     *       @code
     *       Mesquite2::QualityAssessor qa(&qa_metric);
     *       qa.disable_printing_results();
     *       @endcode
     */
    void
    execute(Mesquite2::IQInterface &queue,
            const Settings &        settings = Settings());

    /**
     * Use a custom shape improvement algorithm to increase the quality of the
     * mesh. The result is stored within this class instance and is not applied
     * elsewhere until requested. This allows one to make multiple calls to
     * @p execute and adjust the mesh several times before the result is used.
     *
     * @note The typical input type for this generic function is a
     *       Mesquite2::InstructionQueue, although Mesquite::Wrapper
     *       is also supported.
     *
     * @note By default, Mesquite will print the results of quality assessment
     *       to screen. If you would like to prevent this, add this line of code
     *       when you initialize the quality assessor:
     *       @code
     *       Mesquite2::QualityAssessor qa(&qa_metric);
     *       qa.disable_printing_results();
     *       @endcode
     */
    void
    execute(Mesquite2::IQInterface &queue,
            const MPI_Comm          mpi_communicator,
            const Settings &        settings);

    /**
     * @name Directly manipulating Triangulations
     */
    //@{

    /**
     * A function that extracts the necessary information and, most importantly
     * the vertex positions, directly from the triangulation and stores it in
     * its internal data structures for later manipulation.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the exact (reference) state of the
     * triangulation.
     *
     * @tparam mesh A variable of a type that satisfies the requirements of the
     *         @ref ConceptMeshType "MeshType concept".
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     * @param additional_fixed_vertices An IndexSet that represents vertices
     *        that are to be considered fixed. If the @p fix_all_boundary_vertices
     *        flag is set then the boundary vertices will already be taken care
     * of.
     */
    template <template <int, int> class MeshType>
    void
    initialize(const MeshType<dim, spacedim> &mesh,
               const bool                     fix_all_boundary_vertices = true,
               const IndexSet additional_fixed_vertices = IndexSet());

    /**
     * A function that extracts the necessary information and, most importantly
     * the vertex positions, from the triangulation as well as the
     * input vertex displacement vector, and stores it in its internal data
     * structures for later manipulation.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the Eulerian (spatial) view of the
     * triangulation, as determined by vector of
     * @p eulerian_vertex_displacements . If we denote the original @p mesh
     * vertex coordiates at $\mathbf{X}$ and the @p eulerian_vertex_displacements
     * as $\mathbf{u}$, then the new Eulerian vertex coordinates are
     * $\mathbf{x} = \mathbf{X} + \mathbf{u}$.
     *
     * @note It is assumed that the entries for the Eulerian vertex displacement
     *       map @p eulerian_vertex_positions are sequential, i.e. with the layout (for 3d)
     *       $\left[ u^{0}_{x} u^{0}_{y} u^{0}_{z},
     *               u^{1}_{x} u^{1}_{y} u^{1}_{z},
     *               ... ,
     *               u^{n}_{x} u^{n}_{y} u^{n}_{z} \right]$
     *       where $n$ is the number of vertices in the @p mesh.
     *
     * @tparam mesh A variable of a type that satisfies the requirements of the
     *         @ref ConceptMeshType "MeshType concept".
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     * @param additional_fixed_vertices An IndexSet that represents vertices
     *        that are to be considered fixed. If the @p fix_all_boundary_vertices
     *        flag is set then the boundary vertices will already be taken care
     * of.
     */
    template <template <int, int> class MeshType, typename Number>
    void
    initialize(const MeshType<dim, spacedim> &mesh,
               const Vector<Number> &         eulerian_vertex_displacements,
               const bool                     fix_all_boundary_vertices = true,
               const IndexSet additional_fixed_vertices = IndexSet());

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class MeshType, typename Number>
    void
    update_vertex_displacement_map(Vector<Number> &vertex_displacement_map,
                                   const MeshType<dim, spacedim> &mesh) const;

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class MeshType, typename Number>
    void
    update_eulerian_vertex_positions(Vector<Number> &eulerian_vertex_positions,
                                     const MeshType<dim, spacedim> &mesh) const;

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class MeshType>
    void
    move_triangulation_vertices(MeshType<dim, spacedim> &mesh) const;

    //@}

    /**
     * @name Producing Eulerian vectors of vertex displacements
     */
    //@{

    /**
     *
     */
    template <template <int, int> class DoFHandlerType, typename VectorType>
    void
    initialize(const DoFHandlerType<dim, spacedim> &dof_handler,
               const VectorType &                   solution,
               const unsigned int first_displacement_component,
               const bool         fix_all_boundary_vertices = true,
               const IndexSet     additional_fixed_vertices = IndexSet());

    /**
     * A function that extracts the necessary information and, most importantly
     * the vertex positions, from the underlying triangulation as well as the
     * input displacement map, and stores it in its internal data structures
     * for later manipulation.
     *
     * Calling this function implies that the view of the mesh to be built and
     * optimized is that given by the Eulerian (spatial) view of the
     * triangulation, as determined by vector of
     * @p eulerian_vertex_displacements . If the Mapping does not modify the
     * position of the mesh vertices (like, for example, MappingQEulerian does),
     * then this function is equivalent to the one with the same name, and
     * without the mapping argument.
     *
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     */
    template <template <int, int> class MappingType,
              template <int, int> class MeshType>
    void
    initialize(const MappingType<dim, spacedim> &mapping,
               const MeshType<dim, spacedim> &   mesh,
               const bool     fix_all_boundary_vertices = true,
               const IndexSet additional_fixed_vertices = IndexSet());

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class DoFHandlerType, typename VectorType>
    void
    update_mesh_displacement_solution(
      VectorType &                         solution,
      const DoFHandlerType<dim, spacedim> &dof_handler,
      const unsigned int                   first_displacement_component) const;

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class DoFHandlerType, typename VectorType>
    void
    update_eulerian_mesh_solution(
      VectorType &                         solution,
      const DoFHandlerType<dim, spacedim> &dof_handler,
      const unsigned int                   first_displacement_component) const;

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class MappingType,
              template <int, int> class DoFHandlerType,
              typename VectorType>
    void
    update_mesh_displacement_solution(
      VectorType &                         solution,
      const MappingType<dim, spacedim> &   mapping,
      const DoFHandlerType<dim, spacedim> &dof_handler,
      const unsigned int                   first_displacement_component) const;

    /**
     * Adjust the vertices of the triangulation according to the result of the
     * smoothing operations performed by this object.
     */
    template <template <int, int> class MappingType,
              template <int, int> class DoFHandlerType,
              typename VectorType>
    void
    update_eulerian_mesh_solution(
      VectorType &                         solution,
      const MappingType<dim, spacedim> &   mapping,
      const DoFHandlerType<dim, spacedim> &dof_handler,
      const unsigned int                   first_displacement_component) const;

    //@}

#  ifdef DEAL_II_WITH_VERDICT

    /**
     * Interrogating the current state of the mesh
     */
    //@{

    /**
     * Returns the scaled jacobian for all cells.
     */
    Vector<double>
    get_cell_quality() const;

    /**
     * Returns the minumum scaled jacobian for the entire mesh.
     */
    double
    get_minimum_cell_quality() const;

    //@}

#  endif

  private:
    /**
     * A translation between deal.II's vertex ordering and that used by
     * Mesquite.
     */
    static const std::array<unsigned int, GeometryInfo<dim>::vertices_per_cell>
      vertex_ordering;

    /**
     * An object that stores the description of the mesh in a form that
     * Mesquite understands.
     *
     * The @p MeshImpl class, as opposed to the @p ArrayMesh class, provides a
     * mechanism to implement constraints for hanging nodes. In theory
     * it would be better to use an @p  ArrayMesh here, since it can do in-place
     * mesh adaption (that is, moving vertices without copying data).
     * However, the Triangulation classes do not satisfy the data contiguity
     * requirements of the @p  ArrayMesh class, so we need to make a copy of the
     * salient data to perform the adaption. Also it is an advantage to be
     * able to adapt grids with hanging nodes and prescribe some more complex
     * constraints that can only be done with the @p MeshImpl class.
     */
    std::unique_ptr<BaseMeshType> mesquite_mesh;

    /**
     * A pointer to the triangulation that underlies the mesh to be manipulated.
     *
     * We store this pointer simply to check that the user doesn't initialize
     * the mesh with one triangulation and then apply the optimized mesh
     * to another.
     */
    const Triangulation<dim, spacedim> *ptr_triangulation;

    /**
     * A set of constraints that will be applied through a set of smoothing
     * instructions as a pre-processing step.
     */
    ConstraintMatrix pre_smoothing_vertex_constraints;

    /**
     * A set of constraints that will be applied through a set of smoothing
     * instructions as a post-processing step.
     */
    ConstraintMatrix post_smoothing_vertex_constraints;

    /**
     * An object to help translate between displacement field data as referenced
     * by DoFHandlers, and Triangulation vertex data.
     *
     * @note This object is only created when necessary, and is cached until the
     * class goes out of scope.
     */
    mutable std::unique_ptr<internal::MeshMotionData<dim>> ptr_mesh_motion_data;


    /**
     * Initialize the cached mesh motion data object, if necessary.
     */
    template <template <int, int> class MeshType>
    void
    initialize_mesh_motion_data(const MeshType<dim, spacedim> &mesh) const;

    /**
     * Returns the entry index for the cell connectivity vector
     * for a given cell and vertex number. This takes into account
     * the difference in vertex ordering between deal.II and Mesquite.
     */
    unsigned int
    get_local_cell_connectivity_ordering(const unsigned int &cell_index,
                                         const unsigned int &vertex) const;

    /**
     * A function that extracts the necessary information from the @p triangulation
     * and stores it in its internal data structures. The vertex (reference)
     * positions are taken directly from the triangulation.
     *
     * @note It is assumed that the entries for the Eulerian  vertex displacement
     *       map @p eulerian_vertex_positions are sequential, i.e. with the layout (for 3d)
     *       $\left[ v^{0}_{x} v^{0}_{y} v^{0}_{z},
     *               v^{1}_{x} v^{1}_{y} v^{1}_{z},
     *               ... ,
     *               v^{n}_{x} v^{n}_{y} v^{n}_{z} \right]$
     *       where $n$ is the number of vertices in the @p mesh.
     *
     * @param vertex_coords A list of mesh vertex coordinates with a stride of
     * dim between vertices.
     * @param cell_connectivity A vector describing the element connectivity
     * wtih a stride of n_vertices_per_cell between elements.
     */
    template <template <int, int> class MeshType, typename Number>
    void
    build_vertices_and_cell_connectivity(
      std::vector<CoordinateType> &  vertex_coords,
      std::vector<ElementIndexType> &cell_connectivity,
      const MeshType<dim, spacedim> &mesh,
      const Vector<Number> &         eulerian_vertex_displacements) const;

    /**
     * A function that extracts the necessary information from the @p triangulation
     * and stores it in its internal data structures.
     *
     * @param fixed_vertices A list of fixed vertices with a stride of 1 between
     *        vertices.
     * @param fix_all_boundary_vertices A flag to indicate whether or not the
     *        vertices on the triangulation boundary should be considered fixed
     *        in space or moved during mesh adaption.
     */
    IndexSet
    make_boundary_vertex_constraints(
      const Triangulation<dim, spacedim> &triangulation,
      const bool                          fix_all_boundary_vertices) const;

    /**
     * Returns an instruction that can enforce vertex constraints that are
     * required by a Triangulation. Specifically, these are the hanging
     * and periodic vertex constraints.
     *
     * @note Unfortunately this object cannot be copied, and is therefore
     * returned as a pointer.
     */
    template <template <int, int> class MeshType>
    ConstraintMatrix
    make_vertex_constraints(const MeshType<dim, spacedim> &mesh) const;


    /**
     * Returns the constraints that enforce hanging vertex positions
     */
    template <template <int, int> class MeshType>
    ConstraintMatrix
    make_hanging_vertex_constraints(const MeshType<dim, spacedim> &mesh) const;

    /**
     * Returns the constraints that enforce periodic vertex positions
     */
    template <template <int, int> class MeshType>
    ConstraintMatrix
    make_periodic_vertex_constraints(const MeshType<dim, spacedim> &mesh) const;

    /**
     * Create the Mesquite mesh with the given vertex and connectivity data.
     * The fixed vertex constraints are initialized seperately using the
     * set_mesh_fixed_vertices() function.
     *
     * @note We purposely ensure that the vertex indexing given to Mesquite
     * aligns with the indexing used by the input triangulation. This will be
     * exploited in several places.
     */
    void
    build_mesquite_mesh(std::unique_ptr<BaseMeshType> &      mesquite_mesh,
                        const std::vector<CoordinateType> &  vertex_coords,
                        const std::vector<ElementIndexType> &cell_connectivity,
                        const IndexSet &fixed_vertices) const;


    /**
     * Instructions to be run when the setting Settings::vertex_level_algorithm
     * is set to VertexConsideration::all_at_once .
     */
    void
    execute_optimization_multigrid(Mesquite2::IQInterface &queue,
                                   const MPI_Comm          mpi_communicator,
                                   const Settings &        settings);


    /**
     * Instructions to be run when the setting Settings::vertex_level_algorithm
     * is set to VertexConsideration::all_at_once .
     */
    void
    execute_optimization_global(Mesquite2::IQInterface &queue,
                                const MPI_Comm          mpi_communicator,
                                const Settings &        settings);

    /**
     * Perform some adjustments to the mesh before the primary smoothing
     * operation is performed.
     */
    template <typename MesquiteMeshType>
    void
    execute_pre_smoothing_instructions(MesquiteMeshType &     mesh,
                                       Mesquite2::MeshDomain *domain = nullptr);


    /**
     * Perform some final adjustments to ensure that hanging and
     * periodic vertices are correctly positioned.
     */
    template <typename MesquiteMeshType>
    void
    execute_post_smoothing_instructions(
      MesquiteMeshType &     mesh,
      Mesquite2::MeshDomain *domain = nullptr);
  };



  /**
   * Some settings to be applied when using the built in mesh optmization
   * wrappers.
   *
   * @author Jean-Paul Pelteret, 2018
   */
  template <int dim, int spacedim>
  struct MesquiteMeshInterface<dim, spacedim>::Settings
  {
    /**
     * An enumeration that dictates which vertices are to be considered during
     * smoothing.
     */
    enum class VertexConsideration
    {
      /**
       * With this flag chosen, the smoother will be applied multiple times,
       * once at each level of the trianguation. This means that the coarse grid
       * vertices will be optimized first, and then the those belonging to the
       * first refinement level (with the coarsest fixed in space) and so-on.
       *
       * This option is recommended for grids that are complex in shape and
       * discretization because Mesquite cannot take hanging vertices into
       * consideration during smoothing. Hanging vertces are therefore moved
       * into a sub-optimal position during smoothing, and we subsequently
       * correct their position as a post-processing step. This may lead to
       * skewed cells being generated on the refined levels.
       */
      multigrid,
      /**
       * This performs the same operation as @p multigrid on the coarse mesh,
       * with the optmizer chosen at that given by the user. At the finer
       * levels, standard Laplace smoothing is employed. This can help prevent
       * some oscillations in the vertex positions that may occur in complex
       * meshes due to the enforcement of constraints on the coarse mesh
       * positions during the fine mesh sweeps.
       */
      multigrid_laplace,
      /**
       * The vertices belonging to cells on all refinement levels are considered
       * during smoothing. Since the mesh generation and smoothing operations
       * happen only once, this should be the quickest option. It should also
       * work adequately for predominantly Cartesian grids.
       */
      global
    };


    /**
     * @param log_statistics Log the quality asessment statistics to deallog
     * @param max_vertex_movement Termination optimization if no vertex is moved
     *        by more than this distance in the previous solver step. This
     * parameter is only used by some of the algorithms.
     */
    Settings(const enum VertexConsideration vertex_level_algorithm =
               VertexConsideration::multigrid,
             const bool   log_statistics             = false,
             const bool   print_statistics_to_screen = false,
             const double max_vertex_movement        = 1e-9)
      : vertex_level_algorithm(vertex_level_algorithm)
      , log_statistics(log_statistics)
      , print_statistics_to_screen(print_statistics_to_screen)
      , max_vertex_movement(max_vertex_movement)
    {}

    /**
     * The algorithm applied to perform smoothing on the various refinement
     * levels of the triangulation.
     */
    enum VertexConsideration vertex_level_algorithm;

    /**
     * Log any generated statistics to deallog
     */
    bool log_statistics;

    /**
     * Print any generated statistics to screen
     */
    bool print_statistics_to_screen;

    /**
     * Maximum allowed vertex movement at the
     */
    double max_vertex_movement;
  };

} // namespace GridTools



/* -------------- inline functions and template definitions -------------- */



namespace GridTools
{
  // As per usual, Trilinos is a little light on the documentation. You can find
  // the sparse doxygen documentation for Mesquite here:
  // https://software.sandia.gov/mesquite/doc-2.99/html/


  template <int dim, int spacedim>
  const std::array<unsigned int, GeometryInfo<dim>::vertices_per_cell>
    MesquiteMeshInterface<dim, spacedim>::vertex_ordering =
      internal::MesquiteGeometryInfoHelper<dim>::vertex_ordering;



  VertexConstraintInstruction::VertexConstraintInstruction(
    Mesquite2::Mesh &             mesquite_mesh,
    const ConstraintMatrix &      vertex_constraints,
    Mesquite2::ObjectiveFunction *objective_function)
    : Mesquite2::VertexMover(objective_function)
    , mesquite_mesh(&mesquite_mesh)
    , vertex_constraints(vertex_constraints)
  {}



  VertexConstraintInstruction::~VertexConstraintInstruction()
  {}



  std::string
  VertexConstraintInstruction::get_name() const
  {
    return "VertexConstraintInstruction";
  }



  Mesquite2::PatchSet *
  VertexConstraintInstruction::get_patch_set()
  {
    return &patchSet;
  }



  void
  VertexConstraintInstruction::VertexConstraintInstruction::initialize(
    Mesquite2::PatchData &,
    Mesquite2::MsqError &)
  {}



  void
  VertexConstraintInstruction::cleanup()
  {}



  void
  VertexConstraintInstruction::optimize_vertex_positions(Mesquite2::PatchData &,
                                                         Mesquite2::MsqError &)
  {
    // This give a local patch view, which isn't really useful
    // for what we want to do. We want to move any vertices on
    // the global patch, which is not necessarily guarenteed to
    // be what's given as an input to this function
  }



  void
  VertexConstraintInstruction::initialize_mesh_iteration(Mesquite2::PatchData &,
                                                         Mesquite2::MsqError &)
  {}



  void
  VertexConstraintInstruction::terminate_mesh_iteration(
    Mesquite2::PatchData &pd,
    Mesquite2::MsqError & err)
  {
    // Extract all of the vertices
    std::vector<Mesquite2::Mesh::VertexHandle> vertices;
    mesquite_mesh->get_all_vertices(vertices, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
    const unsigned int n_vertices = vertices.size();

    // Retrieve the coordinates of the mesh points
    std::vector<Mesquite2::MsqVertex> vertex_coords(n_vertices);
    mesquite_mesh->vertices_get_coordinates(vertices.data(),
                                            vertex_coords.data(),
                                            n_vertices,
                                            err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    for (unsigned int v = 0; v < n_vertices; ++v)
      {
        if (vertex_constraints.is_constrained(v) == false)
          continue;

        // Determine the new coordiates of the constrained vertex
        // based on those its master vertices
        Mesquite2::Vector3D new_pos;
        typedef const typename ConstraintMatrix::ConstraintLine::Entries
          *                         ConstraintEntriesType;
        const ConstraintEntriesType constraint_entries =
          vertex_constraints.get_constraint_entries(v);
        for (auto constraint : *constraint_entries)
          {
            const types::global_dof_index v_idx_m   = constraint.first;
            const double                  v_coeff_m = constraint.second;
            new_pos += v_coeff_m * vertex_coords[v_idx_m];
          }

        mesquite_mesh->vertex_set_coordinates(vertices[v], new_pos, err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }
  }



  namespace internal
  {
    template <int dim, typename RangeNumberType>
    class ReferencePosition : public Function<dim, RangeNumberType>
    {
    public:
      ReferencePosition(unsigned int n_components      = dim,
                        unsigned int first_u_component = 0)
        : Function<dim, RangeNumberType>(n_components)
        , first_u_component(first_u_component)
      {}

      virtual ~ReferencePosition()
      {}

      virtual void
      vector_value(const Point<dim> &       p,
                   Vector<RangeNumberType> &return_value) const
      {
        Assert(return_value.size() >= (first_u_component + dim),
               ExcInternalError());
        for (unsigned int d = 0; d < dim; ++d)
          return_value[d + first_u_component] = p[d];
      }

      virtual void
      vector_value_list(const std::vector<Point<dim>> &       points,
                        std::vector<Vector<RangeNumberType>> &values) const
      {
        Assert(points.size() == values.size(),
               ExcDimensionMismatch(points.size(), values.size()));
        for (unsigned int p = 0; p < points.size(); ++p)
          vector_value(points[p], values[p]);
      }

    private:
      const unsigned int first_u_component;
    };
  } // namespace internal



  template <int dim, int spacedim>
  MesquiteMeshInterface<dim, spacedim>::MesquiteMeshInterface()
    : ptr_triangulation(nullptr)
  {
    AssertThrow(
      dim > 1, ExcMessage("This class has not been implemented for dim == 1."));
  }



  template <int dim, int spacedim>
  MesquiteMeshInterface<dim, spacedim>::~MesquiteMeshInterface()
  {
    clear();
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::clear()
  {
    ptr_mesh_motion_data.reset();
    ptr_triangulation = nullptr;
    mesquite_mesh.reset();
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  MesquiteMeshInterface<dim, spacedim>::MesquiteMeshInterface(
    const MeshType<dim, spacedim> &mesh,
    const bool                     fix_all_boundary_vertices,
    const IndexSet                 additional_fixed_vertices)
    : ptr_triangulation(nullptr)
  {
    AssertThrow(
      dim > 1, ExcMessage("This class has not been implemented for dim == 1."));

    clear();
    initialize(mesh, fix_all_boundary_vertices, additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType, typename Number>
  MesquiteMeshInterface<dim, spacedim>::MesquiteMeshInterface(
    const MeshType<dim, spacedim> &mesh,
    const Vector<Number> &         eulerian_vertex_positions,
    const bool                     fix_all_boundary_vertices,
    const IndexSet                 additional_fixed_vertices)
    : ptr_triangulation(nullptr)
  {
    AssertThrow(
      dim > 1, ExcMessage("This class has not been implemented for dim == 1."));

    clear();
    initialize(mesh,
               eulerian_vertex_positions,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class DoFHandlerType, typename VectorType>
  MesquiteMeshInterface<dim, spacedim>::MesquiteMeshInterface(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const VectorType &                   solution,
    const unsigned int                   first_displacement_component,
    const bool                           fix_all_boundary_vertices,
    const IndexSet                       additional_fixed_vertices)
    : ptr_triangulation(nullptr)
  {
    AssertThrow(
      dim > 1, ExcMessage("This class has not been implemented for dim == 1."));

    clear();
    initialize(dof_handler,
               solution,
               first_displacement_component,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MappingType,
            template <int, int> class MeshType>
  MesquiteMeshInterface<dim, spacedim>::MesquiteMeshInterface(
    const MappingType<dim, spacedim> &mapping,
    const MeshType<dim, spacedim> &   mesh,
    const bool                        fix_all_boundary_vertices,
    const IndexSet                    additional_fixed_vertices)
    : ptr_triangulation(nullptr)
  {
    AssertThrow(
      dim > 1, ExcMessage("This class has not been implemented for dim == 1."));

    clear();
    initialize(mapping,
               mesh,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::initialize(
    const MeshType<dim, spacedim> &mesh,
    const bool                     fix_all_boundary_vertices,
    const IndexSet                 additional_fixed_vertices)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();
    // Initialize an empty Eulerian vertex displacement map
    const Vector<double> eulerian_vertex_displacements(dim * tria.n_vertices());
    initialize(mesh,
               eulerian_vertex_displacements,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MappingType,
            template <int, int> class MeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::initialize(
    const MappingType<dim, spacedim> &mapping,
    const MeshType<dim, spacedim> &   mesh,
    const bool                        fix_all_boundary_vertices,
    const IndexSet                    additional_fixed_vertices)
  {
    Vector<double> eulerian_vertex_displacements =
      GridTools::get_vertex_positions(mapping, mesh);
    eulerian_vertex_displacements -=
      GridTools::get_lagrangian_vertex_positions(mesh);
    initialize(mesh,
               eulerian_vertex_displacements,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType, typename Number>
  void
  MesquiteMeshInterface<dim, spacedim>::initialize(
    const MeshType<dim, spacedim> &mesh,
    const Vector<Number> &         eulerian_vertex_displacements,
    const bool                     fix_all_boundary_vertices,
    const IndexSet                 additional_fixed_vertices)
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();
    ptr_triangulation                        = &tria;

    Assert(
      !(dynamic_cast<
        const parallel::distributed::Triangulation<dim, spacedim> *>(&tria)),
      ExcMessage("Parallel distributed triangulations are not yet supported."));

    std::vector<CoordinateType>   vertex_coords;
    std::vector<ElementIndexType> cell_connectivity;
    build_vertices_and_cell_connectivity(vertex_coords,
                                         cell_connectivity,
                                         tria,
                                         eulerian_vertex_displacements);

    // Check number of vertices and cells
    Assert(vertex_coords.size() / meshdim == tria.n_vertices(),
           ExcDimensionMismatch(vertex_coords.size() / meshdim,
                                tria.n_vertices()));
    Assert(cell_connectivity.size() / GeometryInfo<dim>::vertices_per_cell ==
             tria.n_active_cells(),
           ExcDimensionMismatch(cell_connectivity.size() /
                                  GeometryInfo<dim>::vertices_per_cell,
                                tria.n_active_cells()));

    // We initially create a mesh that already contains the fixed boundary
    // constraints if they're requested. Otherwise it is necessary to
    // provide some other constraints via a call to set_mesh_fixed_vertices();
    IndexSet fixed_vertices =
      make_boundary_vertex_constraints(tria, fix_all_boundary_vertices);
    fixed_vertices.add_indices(additional_fixed_vertices);

    // Finally create the mesh to pass to Mesquite
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert((vertex_coords.size() / meshdim) == ptr_triangulation->n_vertices(),
           ExcDimensionMismatch((vertex_coords.size() / meshdim),
                                ptr_triangulation->n_vertices()));
    build_mesquite_mesh(mesquite_mesh,
                        vertex_coords,
                        cell_connectivity,
                        fixed_vertices);

    // Finalise the other vertex constraints to be enforced after smoothing.
    // For now, this is only hanging node and periodic constraints
    post_smoothing_vertex_constraints = make_vertex_constraints(mesh);
    if (post_smoothing_vertex_constraints.n_constraints() > 0)
      {
        Assert(ptr_triangulation->has_hanging_nodes() ||
                 ptr_triangulation->get_periodic_face_map().size() > 0,
               ExcMessage("Unknown set of constraints added"));
      }
  }



  template <int dim, int spacedim>
  template <template <int, int> class DoFHandlerType, typename VectorType>
  void
  MesquiteMeshInterface<dim, spacedim>::initialize(
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const VectorType &                   solution_with_mesh_motion_field,
    const unsigned int                   first_mesh_motion_component,
    const bool                           fix_all_boundary_vertices,
    const IndexSet                       additional_fixed_vertices)
  {
    initialize_mesh_motion_data(dof_handler);
    Assert(ptr_mesh_motion_data, ExcNotInitialized());

    const Vector<typename VectorType::value_type> vertex_displacement_map =
      ptr_mesh_motion_data->extract_vertex_displacement_map(
        dof_handler,
        solution_with_mesh_motion_field,
        first_mesh_motion_component);

    initialize(dof_handler,
               vertex_displacement_map,
               fix_all_boundary_vertices,
               additional_fixed_vertices);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::initialize_mesh_motion_data(
    const MeshType<dim, spacedim> &mesh) const
  {
    if (!ptr_mesh_motion_data)
      ptr_mesh_motion_data.reset(new internal::MeshMotionData<dim>(mesh));
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType, typename Number>
  void
  MesquiteMeshInterface<dim, spacedim>::build_vertices_and_cell_connectivity(
    std::vector<CoordinateType> &  vertex_coords,
    std::vector<ElementIndexType> &cell_connectivity,
    const MeshType<dim, spacedim> &mesh,
    const Vector<Number> &         eulerian_vertex_displacements) const
  {
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();

    // Mesh vertex coordinates
    vertex_coords.resize(meshdim * tria.n_vertices());
    const std::vector<Point<dim>> &vertices = tria.get_vertices();
    for (unsigned int v = 0; v < vertices.size(); ++v)
      {
        for (unsigned int d = 0; d < dim; ++d)
          {
            const unsigned int c =
              vertex_map_index_from_global_vertex_index<meshdim>(
                v, d); // Flattened coordinate index (Mesquite)
            const unsigned int e =
              vertex_map_index_from_global_vertex_index<dim>(
                v, d); // Flattened coordinate index (deal.II)
            Assert(c < vertex_coords.size(),
                   ExcIndexRange(c, 0, vertex_coords.size()));
            Assert(e < eulerian_vertex_displacements.size(),
                   ExcIndexRange(e, 0, eulerian_vertex_displacements.size()));
            vertex_coords[c] =
              vertices[v][d] + eulerian_vertex_displacements[e];
          }
      }

    // Element connectivity
    cell_connectivity.resize(GeometryInfo<dim>::vertices_per_cell *
                             tria.n_active_cells());
    unsigned int cell_count = 0;
    for (auto cell : tria.active_cell_iterators())
      {
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            const unsigned int c = get_local_cell_connectivity_ordering(
              cell_count, v); // Flattened global vertex index
            Assert(c < cell_connectivity.size(),
                   ExcIndexRange(c, 0, cell_connectivity.size()));
            cell_connectivity[c] = cell->vertex_index(v);
          }
        ++cell_count;
      }
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::add_mesh_fixed_vertices(
    const IndexSet &fixed_vertices)
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(!fixed_vertices.is_empty(), ExcMessage("Input indexset is empty."));

    const unsigned int last_element = fixed_vertices.n_elements() - 1;
    const unsigned int max_idx      = *(fixed_vertices.at(last_element));
    (void)max_idx;
    Assert(max_idx < ptr_triangulation->n_vertices(),
           ExcIndexRange(max_idx, 0, ptr_triangulation->n_vertices()));

    // Error tracker
    Mesquite2::MsqError err;

    // Extract all of the vertices
    const unsigned int n_vertices = ptr_triangulation->n_vertices();
    std::vector<Mesquite2::Mesh::VertexHandle> vertices;
    mesquite_mesh->get_all_vertices(vertices, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
    Assert(vertices.size() == n_vertices,
           ExcDimensionMismatch(vertices.size(), n_vertices));

    // Get the existing flags
    std::vector<bool> fixed_flag_array;
    mesquite_mesh->vertices_get_fixed_flag(vertices.data(),
                                           fixed_flag_array,
                                           n_vertices,
                                           err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    // Add to the existing flags
    for (auto it = fixed_vertices.begin(); it != fixed_vertices.end(); ++it)
      fixed_flag_array[*it] = true;

    // Retrieve the coordinates of the mesh points
    // For some reason vertices_set_fixed_flag is not declared virtual inside
    // the base class :-/
    // Also, each call takes a different data type for the fixed vertex flag!
    if (SerialMeshType *const mesquite_mesh_serial =
          dynamic_cast<SerialMeshType *>(mesquite_mesh.get()))
      {
        mesquite_mesh_serial->vertices_set_fixed_flag(vertices.data(),
                                                      fixed_flag_array,
                                                      n_vertices,
                                                      err);
      }
    else if (ParallelMeshType *const mesquite_mesh_parallel =
               dynamic_cast<ParallelMeshType *>(mesquite_mesh.get()))
      {
        // Currently, there's an issue with this:
        // Undefined symbols for architecture x86_64:
        // "Mesquite2::ParallelMeshImpl::vertices_set_fixed_flag(void* const*,
        // bool const*, unsigned long, Mesquite2::MsqError&)"
        (void)mesquite_mesh_parallel;
        AssertThrow(
          false,
          ExcMessage(
            "Up till Mesquite version 2.3, the vertices_set_fixed_flag for the "
            "Mesquite2::ParallelMeshImpl class is declared but not defined. You should "
            "call the class constructor that takes in an IndexSet that represents the "
            "fixed vertices."));

        // A working implementation would look something like this:
        // std::unique_ptr<ConstraintType> fixed_flag_array(new ConstraintType
        // [n_vertices]); for (unsigned int c=0; c<n_vertices; ++c)
        //   if (fixed_vertices.is_element(c))
        //     fixed_flag_array.get()[c] = true;
        //   else
        //     fixed_flag_array.get()[c] = false;
        //
        // mesquite_mesh_parallel->vertices_set_fixed_flag(
        //     vertices.data(),
        //     fixed_flag_array.get(),
        //     n_vertices,
        //     err);
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::set_mesh_fixed_vertices(
    const IndexSet &fixed_vertices)
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(!fixed_vertices.is_empty(), ExcMessage("Input indexset is empty."));

    const unsigned int last_element = fixed_vertices.n_elements() - 1;
    const unsigned int max_idx      = *(fixed_vertices.at(last_element));
    (void)max_idx;
    Assert(max_idx < ptr_triangulation->n_vertices(),
           ExcIndexRange(max_idx, 0, ptr_triangulation->n_vertices()));

    // Error tracker
    Mesquite2::MsqError err;

    // Extract all of the vertices
    const unsigned int n_vertices = ptr_triangulation->n_vertices();
    std::vector<Mesquite2::Mesh::VertexHandle> vertices;
    mesquite_mesh->get_all_vertices(vertices, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
    Assert(vertices.size() == n_vertices,
           ExcDimensionMismatch(vertices.size(), n_vertices));

    // Retrieve the coordinates of the mesh points
    // For some reason vertices_set_fixed_flag is not declared virtual inside
    // the base class :-/
    // Also, each call takes a different data type for the fixed vertex flag!
    if (SerialMeshType *const mesquite_mesh_serial =
          dynamic_cast<SerialMeshType *>(mesquite_mesh.get()))
      {
        std::vector<ConstraintType> fixed_flag_array(n_vertices, false);
        for (auto it = fixed_vertices.begin(); it != fixed_vertices.end(); ++it)
          fixed_flag_array[*it] = true;

        mesquite_mesh_serial->vertices_set_fixed_flag(vertices.data(),
                                                      fixed_flag_array,
                                                      n_vertices,
                                                      err);
      }
    else if (ParallelMeshType *const mesquite_mesh_parallel =
               dynamic_cast<ParallelMeshType *>(mesquite_mesh.get()))
      {
        // Currently, there's an issue with this:
        // Undefined symbols for architecture x86_64:
        // "Mesquite2::ParallelMeshImpl::vertices_set_fixed_flag(void* const*,
        // bool const*, unsigned long, Mesquite2::MsqError&)"
        (void)mesquite_mesh_parallel;
        AssertThrow(
          false,
          ExcMessage(
            "Up till Mesquite version 2.3, the vertices_set_fixed_flag for the "
            "Mesquite2::ParallelMeshImpl class is declared but not defined. You should "
            "call the class constructor that takes in an IndexSet that represents the "
            "fixed vertices."));

        // A working implementation would look something like this:
        // std::unique_ptr<ConstraintType> fixed_flag_array(new ConstraintType
        // [n_vertices]); for (unsigned int c=0; c<n_vertices; ++c)
        //   if (fixed_vertices.is_element(c))
        //     fixed_flag_array.get()[c] = true;
        //   else
        //     fixed_flag_array.get()[c] = false;
        //
        // mesquite_mesh_parallel->vertices_set_fixed_flag(
        //     vertices.data(),
        //     fixed_flag_array.get(),
        //     n_vertices,
        //     err);
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  }



  template <int dim, int spacedim>
  IndexSet
  MesquiteMeshInterface<dim, spacedim>::make_boundary_vertex_constraints(
    const Triangulation<dim, spacedim> &tria,
    const bool                          fix_all_boundary_vertices) const
  {
    IndexSet fixed_vertices(tria.n_vertices());

    if (fix_all_boundary_vertices)
      fixed_vertices.add_indices(GridTools::get_boundary_vertex_indices(tria));

    return fixed_vertices;
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  ConstraintMatrix
  MesquiteMeshInterface<dim, spacedim>::make_vertex_constraints(
    const MeshType<dim, spacedim> &mesh) const
  {
    const ConstraintMatrix hanging_vertex_constraints(
      make_hanging_vertex_constraints(mesh));
    const ConstraintMatrix periodic_vertex_constraints(
      make_periodic_vertex_constraints(mesh));
    ConstraintMatrix all_vertex_constraints;
    all_vertex_constraints.merge(periodic_vertex_constraints,
                                 ConstraintMatrix::right_object_wins);
    all_vertex_constraints.merge(hanging_vertex_constraints,
                                 ConstraintMatrix::right_object_wins);
    all_vertex_constraints.close();

    return all_vertex_constraints;
  }



  namespace internal
  {
    template <int dim>
    void
    make_hanging_vertex_constraint_line(
      ConstraintMatrix &            hanging_vertex_constraints,
      const unsigned int &          slave,
      const std::set<unsigned int> &masters)
    {
      Assert(hanging_vertex_constraints.is_constrained(slave) == false,
             ExcMessage("Vertex constraint is already set"));
      hanging_vertex_constraints.add_line(slave);

      const unsigned int recip_coeff = masters.size();
      if (dim == 2)
        {
          Assert(recip_coeff == GeometryInfo<dim>::vertices_per_face,
                 ExcMessage(
                   "Incorrect number of constraints for this slaved vertex."));
        }
      else if (dim == 3)
        {
          Assert(recip_coeff == GeometryInfo<dim>::vertices_per_face ||
                   recip_coeff == GeometryInfo<dim - 1>::vertices_per_face,
                 ExcMessage(
                   "Incorrect number of constraints for this slaved vertex."));
        }
      const double coefficient = 1.0 / recip_coeff;

      for (auto master : masters)
        hanging_vertex_constraints.add_entry(slave, master, coefficient);
    }
  } // namespace internal



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  ConstraintMatrix
  MesquiteMeshInterface<dim, spacedim>::make_hanging_vertex_constraints(
    const MeshType<dim, spacedim> &mesh) const
  {
    // Lets use our helper data structure to do this for us.
    // This minimises the effort taken to implment this, as we could easily
    // mess this up, e.g. introduce constraint cycles or miss constraints.
    ConstraintMatrix hanging_vertex_constraints;

    if (mesh.get_triangulation().has_hanging_nodes())
      {
        initialize_mesh_motion_data(mesh);
        Assert(ptr_mesh_motion_data, ExcNotInitialized());
        ptr_mesh_motion_data->make_hanging_vertex_constraints(
          hanging_vertex_constraints);
      }

    return hanging_vertex_constraints;
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  ConstraintMatrix
  MesquiteMeshInterface<dim, spacedim>::make_periodic_vertex_constraints(
    const MeshType<dim, spacedim> &mesh) const
  {
    ConstraintMatrix periodic_vertex_constraints;

    if (mesh.get_triangulation().get_periodic_face_map().size() > 0)
      {
        initialize_mesh_motion_data(mesh);
        Assert(ptr_mesh_motion_data, ExcNotInitialized());
        ptr_mesh_motion_data->make_periodic_vertex_constraints(
          periodic_vertex_constraints);
      }

    periodic_vertex_constraints.close();
    return periodic_vertex_constraints;
  }


  template <int dim, int spacedim>
  unsigned int
  MesquiteMeshInterface<dim, spacedim>::get_local_cell_connectivity_ordering(
    const unsigned int &cell_index,
    const unsigned int &vertex) const
  {
    return GeometryInfo<dim>::vertices_per_cell * cell_index +
           vertex_ordering[vertex];
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::build_mesquite_mesh(
    std::unique_ptr<BaseMeshType> &      mesquite_mesh,
    const std::vector<CoordinateType> &  vertex_coords,
    const std::vector<ElementIndexType> &cell_connectivity,
    const IndexSet &                     fixed_vertices) const
  {
    // Note: We keep all checks against the state of ptr_triangulation ojut
    // of this function so that it can be used to generate a generic mesh

    const Mesquite2::EntityTopology topology =
      (dim == 2 ? Mesquite2::QUADRILATERAL : Mesquite2::HEXAHEDRON);

    const unsigned int n_vertices = vertex_coords.size() / meshdim;
    const unsigned int n_cells =
      cell_connectivity.size() / GeometryInfo<dim>::vertices_per_cell;

    // A pointer to a contiguous array of booleans that duplicates the
    // set of input constraints. We need this because @p MeshType expects
    // a pointer to an array of booleans, and std::vector<bool> does
    // not provide this.
    std::unique_ptr<ConstraintType> fixed_vertex_flags(
      new ConstraintType[n_vertices]);
    for (unsigned int c = 0; c < n_vertices; ++c)
      {
        if (fixed_vertices.is_element(c))
          fixed_vertex_flags.get()[c] = true;
        else
          fixed_vertex_flags.get()[c] = false;
      }

    // Always create a serial version of the mesh.
    // If the triangulation is a standard Triangulation or a
    // parallel::shared::Triangulation, then this holds all of the vertices and
    // cells in the triangulation. It the triangulation is a
    //  by a parallel::distributed::Triangulation, then this holds a local view
    // to the vertices and cells owned by the triangulation. The union object
    // will be created and utilized later.
    mesquite_mesh.reset(
      new SerialMeshType(n_vertices,               // number of vertices
                         n_cells,                  // number of elements
                         topology,                 // element type
                         fixed_vertex_flags.get(), // constraint flags
                         vertex_coords.data(),     // vertex coordinates
                         cell_connectivity.data()  // element connectivity
                         ));

    Mesquite2::MsqError err;
    (void)err;
    Assert(mesquite_mesh->get_geometric_dimension(err) == meshdim,
           ExcDimensionMismatch(mesquite_mesh->get_geometric_dimension(err),
                                meshdim));
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
  }



  template <int dim, int spacedim>
  typename MesquiteMeshInterface<dim, spacedim>::BaseMeshType &
  MesquiteMeshInterface<dim, spacedim>::get_mesquite_mesh()
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    return *mesquite_mesh;
  }


  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::write_vtk(
    const std::string &filename) const
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    if (SerialMeshType *const mesquite_mesh_serial =
          dynamic_cast<SerialMeshType *>(mesquite_mesh.get()))
      {
        Mesquite2::MsqError err;
        (void)err;
        mesquite_mesh_serial->write_vtk(filename.c_str(), err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }
    else if (ParallelMeshType *const mesquite_mesh_parallel =
               dynamic_cast<ParallelMeshType *>(mesquite_mesh.get()))
      {
        (void)mesquite_mesh_parallel;
        AssertThrow(false,
                    ExcMessage("Cannot write a mesh for a parallel object."));
      }
  }



  template <int dim, int spacedim>
  const typename MesquiteMeshInterface<dim, spacedim>::BaseMeshType &
  MesquiteMeshInterface<dim, spacedim>::get_mesquite_mesh() const
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    return *mesquite_mesh;
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute(
    const enum MesquiteWrapperTypes type,
    const Settings &                settings)
  {
    execute(type, MPI_COMM_SELF, settings);
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute(Mesquite2::IQInterface &queue,
                                                const Settings &settings)
  {
    execute(queue, MPI_COMM_SELF, settings);
  }



  namespace internal
  {
    namespace
    {
      template <typename MesquiteMeshType>
      void
      execute_instructions(MesquiteMeshType &      mesh,
                           Mesquite2::IQInterface &queue)
      {
        // Error tracker
        Mesquite2::MsqError err;

        // Apply the algorithm
        queue.run_instructions(&mesh, err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }

      template <typename MesquiteMeshType>
      void
      execute_instructions(MesquiteMeshType &      mesh,
                           Mesquite2::MeshDomain * domain,
                           Mesquite2::IQInterface &queue)
      {
        Assert(domain == nullptr, ExcInternalError());
        execute_instructions(mesh, queue);
      }

      void
      execute_instructions(Mesquite2::ParallelMesh &mesh,
                           Mesquite2::MeshDomain *  domain,
                           Mesquite2::IQInterface & queue)
      {
        Assert(domain != nullptr, ExcNotInitialized());

        // Error tracker
        Mesquite2::MsqError err;

        // Apply the algorithm
        queue.run_instructions(&mesh, domain, err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }
    } // namespace
  }   // namespace internal


  //  // ==== TESTING ===
  //  template<int dim, int spacedim>
  //  template<template<int, int> class MeshType>
  //  void
  //  MesquiteMeshInterface<dim,spacedim>::allow_floatng_hanging_vertices (
  //      Mesquite2::IQInterface       &queue,
  //      const MeshType<dim,spacedim> &mesh,
  //      const bool                    flag)
  //  {
  //    const ConstraintMatrix hanging_vertex_constraints
  //    (make_hanging_vertex_constraints(*ptr_triangulation)); const
  //    ConstraintMatrix periodic_vertex_constraints; ConstraintMatrix
  //    hnp_constraints; hnp_constraints.merge(periodic_vertex_constraints,
  //    ConstraintMatrix::right_object_wins);
  //    hnp_constraints.merge(hanging_vertex_constraints,
  //    ConstraintMatrix::right_object_wins); hnp_constraints.close();
  //
  //    if (hnp_constraints.n_constraints() > 0)
  //    {
  //      const unsigned int n_vertices = ptr_triangulation->n_vertices();
  //
  //      // Error tracker
  //      Mesquite2::MsqError err;
  //
  //      // Extract all of the vertices
  //      std::vector<Mesquite2::Mesh::VertexHandle> vertices;
  //      mesquite_mesh->get_all_vertices(vertices, err);
  //      Assert(err==Mesquite2::MsqError::NO_ERROR,
  //             ExcMesquiteError(err));
  //      Assert(vertices.size() == n_vertices,
  //             ExcDimensionMismatch(vertices.size(),n_vertices));
  //
  //
  //      // Retrieve the coordinates of the mesh points
  //      std::vector<Mesquite2::MsqVertex> vertex_coords (n_vertices);
  //      mesquite_mesh->vertices_get_coordinates(
  //        vertices.data(),
  //        vertex_coords.data(),
  //        n_vertices,
  //        err);
  //      Assert(err==Mesquite2::MsqError::NO_ERROR,
  //             ExcMesquiteError(err));
  //
  //      //
  //      https://software.sandia.gov/mesquite/doc-2.99/html/classMesquite_1_1MsqVertex.html#c861d45f03030d8b1d9e2c4f65f168fe
  //      for (unsigned int v=0; v<n_vertices; ++v)
  //      {
  //        if (hnp_constraints.is_constrained(v))
  //        {
  //          if (flag == true)
  //            vertex_coords[v].remove_soft_fixed_flag();
  //          else
  //            vertex_coords[v].set_soft_fixed_flag();
  //        }
  //      }
  //
  ////      // Extract the data flag (?) for each vertex
  ////      std::vector<unsigned char> bytes(n_vertices);
  ////      mesquite_mesh->vertices_get_byte(
  ////          vertices.data(),
  ////          bytes.data(),
  ////          n_vertices,
  ////          err );
  ////      Assert(err==Mesquite2::MsqError::NO_ERROR,
  ////             ExcMesquiteError(err));
  ////
  ////      // Set the slave flag for constrained vertices
  ////      for (unsigned int v=0; v<n_vertices; ++v)
  ////      {
  ////        if (hnp_constraints.is_constrained(v))
  ////        {
  ////          if (flag == true)
  ////            bytes[v] &= ~Mesquite2::MsqVertex::MSQ_CULLED;
  ////          else
  ////            bytes[v] |= Mesquite2::MsqVertex::MSQ_CULLED;
  ////        }
  ////      }
  ////
  ////      // Push the flag changes
  ////      mesquite_mesh->vertices_set_byte(
  ////          vertices.data(),
  ////          bytes.data(),
  ////          n_vertices,
  ////          err );
  ////      Assert(err==Mesquite2::MsqError::NO_ERROR,
  ////             ExcMesquiteError(err));
  //
  //      queue.set_fixed_vertex_mode (Mesquite2::Settings::FIXED_FLAG);
  //      queue.set_slaved_ho_node_mode(Mesquite2::Settings::SLAVE_ALL);
  //    }
  //  }
  //  //=================



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute(
    const enum MesquiteWrapperTypes type,
    const MPI_Comm                  mpi_communicator,
    const Settings &                settings)
  {
    // Error tracker
    Mesquite2::MsqError err;

    // Generate the specified algorithm type
    std::unique_ptr<Mesquite2::Wrapper> mesh_quality_algorithm;
    switch (type)
      {
        case (MesquiteWrapperTypes::laplace):
          mesh_quality_algorithm.reset(new Mesquite2::LaplaceWrapper());
          break;
        case (MesquiteWrapperTypes::shape_improvement):
          mesh_quality_algorithm.reset(new Mesquite2::ShapeImprover());
          break;
        case (MesquiteWrapperTypes::untangler):
          mesh_quality_algorithm.reset(new Mesquite2::UntangleWrapper());
          break;
        case (MesquiteWrapperTypes::minimum_edge_length_improvement):
          mesh_quality_algorithm.reset(new Mesquite2::PaverMinEdgeLengthWrapper(
            settings.max_vertex_movement));
          break;
        case (MesquiteWrapperTypes::size_adapted_shape_improvement):
          mesh_quality_algorithm.reset(
            new Mesquite2::SizeAdaptShapeWrapper(settings.max_vertex_movement));
          break;
        case (MesquiteWrapperTypes::deforming_domain):
          mesh_quality_algorithm.reset(new Mesquite2::DeformingDomainWrapper());
          break;
        default:
          AssertThrow(false,
                      ExcMessage("Mesquite wrapper type is unsupported."));
          break;
      }

    Assert(mesh_quality_algorithm, ExcNotInitialized());

    // Execuate all instructions
    execute(*mesh_quality_algorithm, mpi_communicator, settings);

    //    // ==== TESTING ===
    //    allow_floatng_hanging_vertices(*mesh_quality_algorithm,
    //                                   *ptr_triangulation,
    //                                   true);
    //    // Execuate all instructions
    //    execute(*mesh_quality_algorithm, mpi_communicator);
    //
    ////    allow_floatng_hanging_vertices(*mesh_quality_algorithm,
    ////                                   *ptr_triangulation,
    ////                                   false);
    ////    // Execuate all instructions
    ////    execute(*mesh_quality_algorithm, mpi_communicator);
    //    //=================
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute_optimization_multigrid(
    Mesquite2::IQInterface &queue,
    const MPI_Comm          mpi_communicator,
    const Settings &        settings)
  {
    Assert(settings.vertex_level_algorithm ==
               Settings::VertexConsideration::multigrid ||
             settings.vertex_level_algorithm ==
               Settings::VertexConsideration::multigrid_laplace,
           ExcMessage("Incorrect algorithm setting chosen."));
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(mesquite_mesh, ExcNotInitialized());
    const Triangulation<dim, spacedim> &tria =
      ptr_triangulation->get_triangulation();

    // Early return if we're performing operations on a triangulation with
    // no refined cells
    const unsigned int n_levels = ptr_triangulation->n_levels();
    if (n_levels == 0)
      {
        execute(queue, mpi_communicator, settings);
        return;
      }

    // We will use the same base settings as sent into this function,
    // with the exception that we will always perform evaluations
    // with all vertices considered (albeit with many fixed in space),
    Settings settings_level(settings);
    settings_level.vertex_level_algorithm =
      Settings::VertexConsideration::global;

    // To build our reduced mesh, we will use the mesh originally constructed
    // from the call to initialize(). We can then use the information, such
    // as the Eulerian vertex positions and status of fixed vertices. We also
    // know that their indexing is guarenteed to be aligned, which is a useful
    // property to exploit. Furthermore, by building a map between the indexing
    // used in the reduced and global mesh, we will be able to directly adjust
    // the vertex coordinates in the global mesh one we've performed optmization
    // of a view of one of its coarser levels.
    const unsigned int n_vertices = ptr_triangulation->n_vertices();

    // Error tracker
    Mesquite2::MsqError err;

    // A map of the global vertex indices encountered on a level and those
    // below it.
    IndexSet vertices_up_to_level(n_vertices);
    for (unsigned int l = 0; l < n_levels; ++l)
      {
        deallog.push("Level " + Utilities::int_to_string(l));

        // Extract all of the vertices
        std::vector<Mesquite2::Mesh::VertexHandle> vertices;
        mesquite_mesh->get_all_vertices(vertices, err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
        Assert(vertices.size() == n_vertices,
               ExcDimensionMismatch(vertices.size(), n_vertices));

        // Retrieve the coordinates of the mesh points
        std::vector<Mesquite2::MsqVertex> vertex_coords(n_vertices);
        mesquite_mesh->vertices_get_coordinates(vertices.data(),
                                                vertex_coords.data(),
                                                n_vertices,
                                                err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

        // Get some flags
        std::vector<bool> fixed_flag_array;
        mesquite_mesh->vertices_get_fixed_flag(vertices.data(),
                                               fixed_flag_array,
                                               n_vertices,
                                               err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

        // Freeze the position of all of the vertices used to construct the
        // problem at a lower level. This is done because, for each sweep
        // of the smoother, we consider these vertices to have been optimally
        // placed.
        const IndexSet vertices_before_level = vertices_up_to_level;

        // A map of the global vertex indices used to build up the reduced
        // problem on a given level
        IndexSet vertices_on_level(n_vertices);
        for (auto cell : ptr_triangulation->cell_iterators_on_level(l))
          {
            for (unsigned int v = 0;
                 v < GeometryInfo<spacedim>::vertices_per_cell;
                 ++v)
              {
                const unsigned int idx_global = cell->vertex_index(v);

                // Add all of the vertices on this level
                vertices_on_level.add_index(idx_global);
                // This stored indexset will be incrementally modified so that
                // it containts all of the global indices of the vertices on
                // this level as well as all of those of the levels below it. It
                // will also be used to make up a translation map for this new
                // view of the triangulation.
                vertices_up_to_level.add_index(idx_global);
              }
          }

        // Indexing maps
        std::map<unsigned int, unsigned int> vertex_index_global_to_reduced;
        // Mesh vertex coordinates
        std::vector<CoordinateType> vertex_coords_reduced(
          meshdim * vertices_on_level.n_elements(), 0.0);
        // Fixed vertices
        IndexSet fixed_vertices_reduced(n_vertices);
        // Now perform three operations with the reduced set of vertices:
        for (unsigned int idx_reduced = 0;
             idx_reduced < vertices_on_level.n_elements();
             ++idx_reduced)
          {
            // 1. Compact them into a map that will link the indexing used for
            // the level-reduced problem to the main one and vice versa
            const unsigned int idx_global =
              vertices_on_level.nth_index_in_set(idx_reduced);
            vertex_index_global_to_reduced[idx_global] = idx_reduced;

            // 2. Make the mesh coordinates
            for (unsigned int d = 0; d < dim; ++d)
              {
                const unsigned int c =
                  vertex_map_index_from_global_vertex_index<meshdim>(
                    idx_reduced, d); // Flattened coordinate index (Mesquite)
                Assert(idx_global < vertex_coords.size(),
                       ExcIndexRange(idx_global, 0, vertex_coords.size()));
                vertex_coords_reduced[c] = vertex_coords[idx_global][d];
              }

            // 3. Fixed vertices
            // Note: Additional fine level vertices may be constrained by the
            // user,
            //       or may be on the boundary and thus subject to constraint at
            //       initialization time by the fix_boundary_vertices flag.
            Assert(idx_global < fixed_flag_array.size(),
                   ExcIndexRange(idx_global, 0, fixed_flag_array.size()));
            if (fixed_flag_array[idx_global] == true ||
                vertices_before_level.is_element(idx_global))
              fixed_vertices_reduced.add_index(idx_reduced);
          }

        // Element connectivity
        std::vector<ElementIndexType> cell_connectivity_reduced;
        cell_connectivity_reduced.reserve(GeometryInfo<dim>::vertices_per_cell *
                                          tria.n_cells(l));
        unsigned int cell_count = 0;
        for (auto cell : tria.cell_iterators_on_level(l))
          {
            // Due to the fact that we're now working with a different looking
            // grid to the previous level iteration, is hard how many cells we
            // expect to create without iterating through the triangulation
            // first. So we preallocatd the maximum size for the
            // cell_connectivity vector and can resize this vector without much
            // penalty.
            cell_connectivity_reduced.resize(
              cell_connectivity_reduced.size() +
              GeometryInfo<dim>::vertices_per_cell);

            for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell;
                 ++v)
              {
                const unsigned int idx_global = cell->vertex_index(v);
                const auto it = vertex_index_global_to_reduced.find(idx_global);
                Assert(it != vertex_index_global_to_reduced.end(),
                       ExcMessage("Global index " +
                                  Utilities::int_to_string(idx_global) +
                                  " "
                                  "not found when setting connectivity for "
                                  "active cell on level " +
                                  Utilities::int_to_string(cell->level())));
                const unsigned int idx_reduced = it->second;

                const unsigned int c =
                  this->get_local_cell_connectivity_ordering(
                    cell_count, v); // Flattened global vertex index
                Assert(c < cell_connectivity_reduced.size(),
                       ExcIndexRange(c, 0, cell_connectivity_reduced.size()));

                cell_connectivity_reduced[c] = idx_reduced;
              }
            ++cell_count;
          }

        // Create the post-smoothing constraints
        // These can simply be copied from the existing set of constraints, with
        // the appropriate index transformation
        ConstraintMatrix pre_smoothing_vertex_constraints_reduced;
        for (auto idx_slave_global : vertices_on_level)
          {
            if (this->post_smoothing_vertex_constraints.is_constrained(
                  idx_slave_global))
              {
                const auto it_s =
                  vertex_index_global_to_reduced.find(idx_slave_global);
                Assert(it_s != vertex_index_global_to_reduced.end(),
                       ExcMessage("Global index for slave vertex not found"));
                const unsigned int idx_slave_reduced = it_s->second;
                pre_smoothing_vertex_constraints_reduced.add_line(
                  idx_slave_reduced);

                typedef const typename ConstraintMatrix::ConstraintLine::Entries
                  *                         ConstraintEntriesType;
                const ConstraintEntriesType constraint_entries =
                  this->post_smoothing_vertex_constraints
                    .get_constraint_entries(idx_slave_global);
                for (auto constraint : *constraint_entries)
                  {
                    const types::global_dof_index idx_master_global =
                      constraint.first;
                    const auto it_m =
                      vertex_index_global_to_reduced.find(idx_master_global);
                    Assert(it_m != vertex_index_global_to_reduced.end(),
                           ExcMessage(
                             "Global index for master vertex not found"));
                    const unsigned int idx_master_reduced = it_m->second;
                    const double       v_coeff_m          = constraint.second;
                    pre_smoothing_vertex_constraints_reduced.add_entry(
                      idx_slave_reduced, idx_master_reduced, v_coeff_m);
                  }

                // At the same time we should mark this hanging vertex as fixed
                // since we'll move it manually and then expect its position not
                // to be adjusted by the optmizer.
                fixed_vertices_reduced.add_index(idx_slave_reduced);
              }
          }


        // Create our reduced mesh
        std::unique_ptr<BaseMeshType> mesquite_mesh_reduced;
        build_mesquite_mesh(mesquite_mesh_reduced,
                            vertex_coords_reduced,
                            cell_connectivity_reduced,
                            fixed_vertices_reduced);

        // Check that all of the vertices are in the correct position
        {
          const unsigned int n_vertices_reduced =
            vertices_on_level.n_elements();
          std::vector<Mesquite2::Mesh::VertexHandle> vertices_reduced;
          mesquite_mesh_reduced->get_all_vertices(vertices_reduced, err);
          Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
          Assert(vertices_reduced.size() == n_vertices_reduced,
                 ExcDimensionMismatch(vertices_reduced.size(),
                                      n_vertices_reduced));

          // Retrieve the coordinates of the mesh points
          std::vector<Mesquite2::MsqVertex> vertex_coords_reduced(
            n_vertices_reduced);
          mesquite_mesh_reduced->vertices_get_coordinates(
            vertices_reduced.data(),
            vertex_coords_reduced.data(),
            n_vertices_reduced,
            err);
          Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

          for (auto idx_global : vertices_on_level)
            {
              const auto it = vertex_index_global_to_reduced.find(idx_global);
              Assert(
                it != vertex_index_global_to_reduced.end(),
                ExcMessage(
                  "Global index not found when checking vertex positions"));
              const unsigned int idx_slave = it->second;
              (void)idx_slave;
              Assert(idx_slave < n_vertices_reduced,
                     ExcIndexRange(idx_slave, 0, n_vertices_reduced));

              for (unsigned int d = 0; d < meshdim; ++d)
                {
                  Assert(std::abs(vertex_coords_reduced[idx_slave][d] -
                                  vertex_coords[idx_global][d]) < 1e-12,
                         ExcMessage("Vertex positions do not match"));
                }
            }
        }

        // Build another interface object to perform mesh smoothing on the
        // reduced problem
        {
          // Instead of getting it to configure itself, we prescribe the mesh
          // and pre-smoothing constraints
          GridTools::MesquiteMeshInterface<dim> mesquite_interface_reduced;
          mesquite_interface_reduced.mesquite_mesh.swap(mesquite_mesh_reduced);
          mesquite_mesh_reduced.reset();

          // Here we use pre-smoothing constraints to correct any hanging
          // vertices that form in the islands on a level. We can do this
          // because the sweep on the previous level would have determined the
          // (fixed) position of the master vertices.
          mesquite_interface_reduced.pre_smoothing_vertex_constraints.copy_from(
            pre_smoothing_vertex_constraints_reduced);
          mesquite_interface_reduced.pre_smoothing_vertex_constraints.close();

          bool write_vtk_debug = false;
          if (write_vtk_debug)
            {
              static unsigned int it = -1;
              if (l == 0)
                ++it;
              const std::string filename =
                "mesquite_mesh_multigrid.it_" + Utilities::int_to_string(it) +
                ".level_" + Utilities::int_to_string(l) + ".0.vtk";

              if (SerialMeshType *const mesquite_mesh_serial =
                    dynamic_cast<SerialMeshType *>(
                      mesquite_interface_reduced.mesquite_mesh.get()))
                {
                  Mesquite2::MsqError err;
                  (void)err;
                  mesquite_mesh_serial->write_vtk(filename.c_str(), err);
                  Assert(err == Mesquite2::MsqError::NO_ERROR,
                         ExcMesquiteError(err));
                }
              else if (ParallelMeshType *const mesquite_mesh_parallel =
                         dynamic_cast<ParallelMeshType *>(
                           mesquite_interface_reduced.mesquite_mesh.get()))
                {
                  (void)mesquite_mesh_parallel;
                  AssertThrow(false,
                              ExcMessage(
                                "Cannot write a mesh for a parallel object."));
                }
            }

          // Execute the same function using the other vertex consideration
          // setting This will actually result in some mesh adjustment being
          // performed.
          if (settings.vertex_level_algorithm ==
              Settings::VertexConsideration::multigrid_laplace)
            {
              Assert(settings_level.vertex_level_algorithm ==
                       Settings::VertexConsideration::global,
                     ExcInternalError());
              mesquite_interface_reduced.execute(MesquiteWrapperTypes::laplace,
                                                 mpi_communicator,
                                                 settings_level);
            }
          else
            mesquite_interface_reduced.execute_optimization_global(
              queue, mpi_communicator, settings_level);

          if (write_vtk_debug)
            {
              static unsigned int it = -1;
              if (l == 0)
                ++it;
              const std::string filename =
                "mesquite_mesh_multigrid.it_" + Utilities::int_to_string(it) +
                ".level_" + Utilities::int_to_string(l) + ".1.vtk";

              if (SerialMeshType *const mesquite_mesh_serial =
                    dynamic_cast<SerialMeshType *>(
                      mesquite_interface_reduced.mesquite_mesh.get()))
                {
                  Mesquite2::MsqError err;
                  (void)err;
                  mesquite_mesh_serial->write_vtk(filename.c_str(), err);
                  Assert(err == Mesquite2::MsqError::NO_ERROR,
                         ExcMesquiteError(err));
                }
              else if (ParallelMeshType *const mesquite_mesh_parallel =
                         dynamic_cast<ParallelMeshType *>(
                           mesquite_interface_reduced.mesquite_mesh.get()))
                {
                  (void)mesquite_mesh_parallel;
                  AssertThrow(false,
                              ExcMessage(
                                "Cannot write a mesh for a parallel object."));
                }
            }

          // Copy the updated mesh coordinates back to the master Mesquite mesh
          // Extract all of the vertices
          const unsigned int n_vertices_reduced =
            vertices_on_level.n_elements();
          std::vector<Mesquite2::Mesh::VertexHandle> vertices_reduced_optimized;
          mesquite_interface_reduced.mesquite_mesh->get_all_vertices(
            vertices_reduced_optimized, err);
          Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
          Assert(vertices_reduced_optimized.size() == n_vertices_reduced,
                 ExcDimensionMismatch(vertices_reduced_optimized.size(),
                                      n_vertices_reduced));

          // Retrieve the coordinates of the mesh points
          std::vector<Mesquite2::MsqVertex> vertex_coords_reduced_optimized(
            n_vertices_reduced);
          mesquite_interface_reduced.mesquite_mesh->vertices_get_coordinates(
            vertices_reduced_optimized.data(),
            vertex_coords_reduced_optimized.data(),
            n_vertices_reduced,
            err);
          Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

          // Make the update
          for (auto idx_global : vertices_on_level)
            {
              const auto it = vertex_index_global_to_reduced.find(idx_global);
              Assert(
                it != vertex_index_global_to_reduced.end(),
                ExcMessage(
                  "Global index not found when copying vertex positions to global mesh"));
              const unsigned int idx_slave = it->second;
              Assert(idx_slave < n_vertices_reduced,
                     ExcIndexRange(idx_slave, 0, n_vertices_reduced));

              Mesquite2::Vector3D new_pos;
              for (unsigned int d = 0; d < meshdim; ++d)
                new_pos[d] = vertex_coords_reduced_optimized[idx_slave][d];
              mesquite_mesh->vertex_set_coordinates(vertices[idx_global],
                                                    new_pos,
                                                    err);
              Assert(err == Mesquite2::MsqError::NO_ERROR,
                     ExcMesquiteError(err));
            }
        }

        deallog.pop();
      }
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute_optimization_global(
    Mesquite2::IQInterface &queue,
    const MPI_Comm          mpi_communicator,
    const Settings &        settings)
  {
    Assert(settings.vertex_level_algorithm ==
             Settings::VertexConsideration::global,
           ExcMessage("Incorrect algorithm setting chosen."));
    //    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(mesquite_mesh, ExcNotInitialized());

    // Disable obnoxious status printed to screen
    if (Mesquite2::Wrapper *mesh_quality_algorithm =
          dynamic_cast<Mesquite2::Wrapper *>(&queue))
      {
        if (!settings.print_statistics_to_screen)
          mesh_quality_algorithm->quality_assessor().disable_printing_results();
      }

    // Error tracker
    Mesquite2::MsqError err;

    const bool is_parallel_distributed =
      (dynamic_cast<const parallel::distributed::Triangulation<dim, spacedim>
                      *>(ptr_triangulation));
    AssertThrow(!is_parallel_distributed,
                ExcMessage(
                  "Parallel distributed triangulation are not supported yet."));

    // Improve the mesh
    // Mesquite has a number of interfaces to the instruction queue, and
    // we must consider all of them based on whether we have a planar (2d)
    // domain, and are running in serial or parallel. Yay!
    // Note: There should be some constraints added to the constraints
    // vector, otherwise this will likely not converge on a solution.
    if (dim == 2)
      {
        // Surface to constrain the 2d elements to
        // TODO: This would not necessarily be correct for the codimension 1
        // case.
        Assert(dim == spacedim,
               ExcMessage("Codimension 1 case is not yet supported."));
        Mesquite2::PlanarDomain planar_domain(Mesquite2::PlanarDomain::XY);

        if (!is_parallel_distributed)
          {
            // Build a view of the domain
            Mesquite2::MeshDomainAssoc mesh_and_domain(mesquite_mesh.get(),
                                                       &planar_domain);

            // Perform whatever pre-smoothing operations we're asked to do.
            execute_pre_smoothing_instructions(mesh_and_domain);

            // Improve the mesh
            internal::execute_instructions(mesh_and_domain, queue);

            // We ensure that all required mesh constraints are enforced, even
            // if the user does not add them to the optimization instruction
            // queue.
            execute_post_smoothing_instructions(mesh_and_domain);
          }
        else
          {
            // Build a parallel view to the grid
            const unsigned int n_mpi_processes =
              Utilities::MPI::n_mpi_processes(mpi_communicator);
            const unsigned int this_mpi_process =
              Utilities::MPI::this_mpi_process(mpi_communicator);
            const std::string gid_name =
              "Mesquite_ID_" + Utilities::int_to_string(n_mpi_processes);
            const std::string pid_name =
              "Mesquite_ID_" + Utilities::int_to_string(this_mpi_process);
            ParallelMeshType parallel_mesh(mesquite_mesh.get(),
                                           gid_name.c_str(),
                                           pid_name.c_str());

            // Perform whatever pre-smoothing operations we're asked to do.
            execute_pre_smoothing_instructions(parallel_mesh, &planar_domain);

            // Improve the mesh
            internal::execute_instructions(parallel_mesh,
                                           &planar_domain,
                                           queue);

            // We ensure that all required mesh constraints are enforced, even
            // if the user does not add them to the optimization instruction
            // queue.
            execute_post_smoothing_instructions(parallel_mesh, &planar_domain);
            AssertThrow(false, ExcNotImplemented());
          }
      }
    else
      {
        if (!is_parallel_distributed)
          {
            // Perform whatever pre-smoothing operations we're asked to do.
            execute_pre_smoothing_instructions(*mesquite_mesh);

            // Improve the mesh
            internal::execute_instructions(*mesquite_mesh, queue);

            // We ensure that all required mesh constraints are enforced, even
            // if the user does not add them to the optimization instruction
            // queue.
            execute_post_smoothing_instructions(*mesquite_mesh);
          }
        else
          {
            // Build a parallel view to the grid
            const unsigned int n_mpi_processes =
              Utilities::MPI::n_mpi_processes(mpi_communicator);
            const unsigned int this_mpi_process =
              Utilities::MPI::this_mpi_process(mpi_communicator);
            const std::string gid_name =
              "Mesquite_ID_" + Utilities::int_to_string(n_mpi_processes);
            const std::string pid_name =
              "Mesquite_ID_" + Utilities::int_to_string(this_mpi_process);
            ParallelMeshType parallel_mesh(mesquite_mesh.get(),
                                           gid_name.c_str(),
                                           pid_name.c_str());

            // Perform whatever pre-smoothing operations we're asked to do.
            execute_pre_smoothing_instructions(parallel_mesh);

            // Improve the mesh
            internal::execute_instructions(parallel_mesh, queue);

            // We ensure that all required mesh constraints are enforced, even
            // if the user does not add them to the optimization instruction
            // queue.
            execute_post_smoothing_instructions(parallel_mesh);
          }
      }

    // Print the output to logstream
    if (Mesquite2::Wrapper *mesh_quality_algorithm =
          dynamic_cast<Mesquite2::Wrapper *>(&queue))
      {
        if (settings.log_statistics)
          {
            // Make a magic black-hole-stream
            // https://stackoverflow.com/a/14232652
            class Buffer : public std::stringbuf
            {
            public:
              virtual ~Buffer()
              {}
              virtual int
              sync()
              {
                // add this->str() to database here
                // (optionally clear buffer afterwards)
                return 0;
              }
            } buffer;

            std::ostream black_hole_stream(&buffer); // Feed me garbage!
            mesh_quality_algorithm->quality_assessor().print_summary(
              black_hole_stream);
            deallog << black_hole_stream.rdbuf() << std::flush;
          }
      }
  }



  template <int dim, int spacedim>
  void
  MesquiteMeshInterface<dim, spacedim>::execute(Mesquite2::IQInterface &queue,
                                                const MPI_Comm mpi_communicator,
                                                const Settings &settings)
  {
    if (settings.vertex_level_algorithm ==
          Settings::VertexConsideration::multigrid ||
        settings.vertex_level_algorithm ==
          Settings::VertexConsideration::multigrid_laplace)
      {
        deallog.push("Mesquite multigrid smoother");
        execute_optimization_multigrid(queue, mpi_communicator, settings);
        deallog.pop();
      }
    else if (settings.vertex_level_algorithm ==
             Settings::VertexConsideration::global)
      {
        deallog.push("Mesquite global smoother");
        execute_optimization_global(queue, mpi_communicator, settings);
        deallog.pop();
      }
    else
      {
        Assert(
          false,
          ExcMessage(
            "Algorithm dictating the treatment of vertices is not implemented."));
      }
  }



  namespace internal
  {
    /**
     * It is required that all Mesquite2::InstructionQueue have a prescribed
     * master smoother, but in some cases we might just want to perform a
     * whole queue of "preconditioning" operations. This dummy smoother can
     * therefore be passed as the master smoother to any InstructionQueue that
     * doesn't actually need a smoother.
     */
    class DummySmoother : public Mesquite2::VertexMover
    {
    public:
      DummySmoother()
        : Mesquite2::VertexMover(nullptr)
      {}

      virtual ~DummySmoother()
      {}

      virtual std::string
      get_name() const
      {
        return "DummySmoother";
      }

      virtual Mesquite2::PatchSet *
      get_patch_set()
      {
        return &patchSet;
      }

    protected:
      virtual void
      initialize(Mesquite2::PatchData &, Mesquite2::MsqError &)
      {}

      virtual void
      cleanup()
      {}

      virtual void
      optimize_vertex_positions(Mesquite2::PatchData &, Mesquite2::MsqError &)
      {}

      virtual void
      initialize_mesh_iteration(Mesquite2::PatchData &, Mesquite2::MsqError &)
      {}

      virtual void
      terminate_mesh_iteration(Mesquite2::PatchData &, Mesquite2::MsqError &)
      {}

    private:
      /**
       * We will work on the global view of the mesh,
       * not a patch. So this does nothing, but we
       * need to pass this back in one of the pure
       * virtual functions.
       */
      Mesquite2::VertexPatches patchSet;
    };
  } // namespace internal



  template <int dim, int spacedim>
  template <typename MesquiteMeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::execute_pre_smoothing_instructions(
    MesquiteMeshType &     mesh,
    Mesquite2::MeshDomain *domain)
  {
    //    Assert(ptr_triangulation, ExcNotInitialized());

    // Error tracker
    Mesquite2::MsqError err;

    // A set of instructions to be run
    Mesquite2::InstructionQueue pre_smoothing_instructions;

    // There's no simple way to ensure that hanging and periodic vertices are
    // always correctly condsidered, other than to insert some sort of
    // VertexMover into the optimizer's instruction list. Its a pain to get a
    // user to do this, so we instead create this class to hijack the set of
    // solver instructions and inject in the necessary VertexMovers that will
    // do the job for us.
    std::unique_ptr<GridTools::VertexConstraintInstruction> ptr_hnsp;
    if (pre_smoothing_vertex_constraints.n_constraints() > 0)
      {
        ptr_hnsp.reset(
          new VertexConstraintInstruction(*mesquite_mesh,
                                          pre_smoothing_vertex_constraints));
        pre_smoothing_instructions.add_preconditioner(ptr_hnsp.get(), err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }

    // It is a strict requirement that one set a master smoother, so we just
    // use a dummy one. We could in fact have set the
    // VertexConstraintInstruction as the master smoother, since its also a
    // Mesquite2::VertexMover, but we will use this dummy just so retain the
    // flexibility of adding more instructions to the list in future.
    internal::DummySmoother smoother;
    pre_smoothing_instructions.set_master_quality_improver(&smoother, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    // Apply customised post-smoothing instructions,
    // e.g. enforce hanging and periodic vertices
    internal::execute_instructions(mesh, domain, pre_smoothing_instructions);
  }



  template <int dim, int spacedim>
  template <typename MesquiteMeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::execute_post_smoothing_instructions(
    MesquiteMeshType &     mesh,
    Mesquite2::MeshDomain *domain)
  {
    //    Assert(ptr_triangulation, ExcNotInitialized());

    // Error tracker
    Mesquite2::MsqError err;

    // A set of instructions to be run
    Mesquite2::InstructionQueue post_smoothing_instructions;

    // There's no simple way to ensure that hanging and periodic vertices are
    // always correctly condsidered, other than to insert some sort of
    // VertexMover into the optimizer's instruction list. Its a pain to get a
    // user to do this, so we instead create this class to hijack the set of
    // solver instructions and inject in the necessary VertexMovers that will
    // do the job for us.
    std::unique_ptr<GridTools::VertexConstraintInstruction> ptr_hnsp;
    if (post_smoothing_vertex_constraints.n_constraints() > 0)
      {
        ptr_hnsp.reset(
          new VertexConstraintInstruction(*mesquite_mesh,
                                          post_smoothing_vertex_constraints));
        post_smoothing_instructions.add_preconditioner(ptr_hnsp.get(), err);
        Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
      }

    // It is a strict requirement that one set a master smoother, so we just
    // use a dummy one. We could in fact have set the
    // VertexConstraintInstruction as the master smoother, since its also a
    // Mesquite2::VertexMover, but we will use this dummy just so retain the
    // flexibility of adding more instructions to the list in future.
    internal::DummySmoother smoother;
    post_smoothing_instructions.set_master_quality_improver(&smoother, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    // Apply customised post-smoothing instructions,
    // e.g. enforce hanging and periodic vertices
    internal::execute_instructions(mesh, domain, post_smoothing_instructions);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType, typename Number>
  void
  MesquiteMeshInterface<dim, spacedim>::update_vertex_displacement_map(
    Vector<Number> &               vertex_displacement_map,
    const MeshType<dim, spacedim> &mesh) const
  {
    update_eulerian_vertex_positions(vertex_displacement_map, mesh);
    vertex_displacement_map -= GridTools::get_lagrangian_vertex_positions(mesh);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType, typename Number>
  void
  MesquiteMeshInterface<dim, spacedim>::update_eulerian_vertex_positions(
    Vector<Number> &               eulerian_vertex_positions,
    const MeshType<dim, spacedim> &mesh) const
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    Assert(ptr_triangulation, ExcNotInitialized());
    const Triangulation<dim, spacedim> &tria = mesh.get_triangulation();
    Assert(&tria == ptr_triangulation, ExcTriangulationMismatch());

    // Error tracker
    Mesquite2::MsqError err;

    // Extract all of the vertices
    std::vector<Mesquite2::Mesh::VertexHandle> vertices;
    mesquite_mesh->get_all_vertices(vertices, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
    Assert(vertices.size() == tria.n_vertices(),
           ExcDimensionMismatch(vertices.size(), tria.n_vertices()));

    // Retrieve the coordinates of the mesh points
    std::vector<Mesquite2::MsqVertex> vertex_coords(tria.n_vertices());
    mesquite_mesh->vertices_get_coordinates(vertices.data(),
                                            vertex_coords.data(),
                                            tria.n_vertices(),
                                            err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    // Reinitialise the Eulerian map
    if (eulerian_vertex_positions.size() != dim * tria.n_vertices())
      eulerian_vertex_positions.reinit(dim * tria.n_vertices());
    eulerian_vertex_positions = 0.0;

    // Build the Eulerian map
    std::vector<bool> vertex_touched(tria.n_vertices(), false);
    for (auto cell : tria.active_cell_iterators())
      {
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            const unsigned int idx_vtx = cell->vertex_index(v);
            Assert(idx_vtx < vertex_touched.size(),
                   ExcIndexRange(idx_vtx, 0, vertex_touched.size()));
            if (!vertex_touched[idx_vtx])
              {
                for (unsigned int d = 0; d < dim; ++d)
                  {
                    const unsigned int idx_vtx_map =
                      vertex_map_index_from_global_vertex_index<dim>(idx_vtx,
                                                                     d);
                    Assert(idx_vtx_map < eulerian_vertex_positions.size(),
                           ExcIndexRange(idx_vtx_map,
                                         0,
                                         eulerian_vertex_positions.size()));
                    eulerian_vertex_positions[idx_vtx_map] =
                      vertex_coords[idx_vtx][d];
                  }

                vertex_touched[idx_vtx] = true;
              }
          }
      }
  }



  template <int dim, int spacedim>
  template <template <int, int> class MeshType>
  void
  MesquiteMeshInterface<dim, spacedim>::move_triangulation_vertices(
    MeshType<dim, spacedim> &mesh) const
  {
    Triangulation<dim, spacedim> &tria = mesh.get_triangulation();
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(&tria == ptr_triangulation, ExcTriangulationMismatch());

    Vector<double> eulerian_vertex_positions(dim * tria.n_vertices());
    update_eulerian_vertex_positions(eulerian_vertex_positions, tria);

    GridTools::move_triangulation_vertices(tria, eulerian_vertex_positions);
  }



  namespace internal
  {
    template <int dim, int spacedim, typename VectorType>
    VectorType
    get_lagrangian_mesh_solution(
      const DoFHandler<dim, spacedim> &dof_handler,
      const VectorType &               prototype_vector,
      const unsigned int               first_displacement_component)
    {
      const FEValuesExtractors::Vector displacement(
        first_displacement_component);
      VectorType solution_X;
      solution_X.reinit(prototype_vector);
      VectorTools::interpolate(
        dof_handler,
        GridTools::internal::ReferencePosition<dim>(
          first_displacement_component, dof_handler.get_fe().n_components()),
        solution_X,
        dof_handler.get_fe().component_mask(displacement));
      return solution_X;
    }


    template <int dim,
              int spacedim,
              template <int, int> class DoFHandlerType,
              typename VectorType>
    VectorType
    get_mesh_solution(const DoFHandlerType<dim, spacedim> &dof_handler,
                      const VectorType &                   solution,
                      const unsigned int first_mesh_motion_component_index)
    {
      // Duplicate out solution vector, for which we will
      // remove all entries not related to the mesh motion
      // field
      VectorType solution_component_u = solution;

      std::vector<types::global_dof_index> dof_indices;
      for (typename DoFHandler<dim, spacedim>::active_cell_iterator cell =
             dof_handler.begin_active();
           cell != dof_handler.end();
           ++cell)
        {
          const unsigned int n_dofs_per_cell = cell->get_fe().dofs_per_cell;
          if (dof_indices.size() != n_dofs_per_cell)
            dof_indices.resize(n_dofs_per_cell);
          cell->get_dof_indices(dof_indices);

          // Clear entries not related to mesh motion field.
          for (unsigned int I = 0; I < n_dofs_per_cell; ++I)
            {
              if (!is_mesh_motion_dof_index(cell->get_fe(),
                                            I,
                                            first_mesh_motion_component_index))
                solution_component_u[dof_indices[I]] = 0.0;
            }
        }

      solution_component_u.compress(VectorOperation::insert);
      return solution_component_u;
    }

  } // namespace internal



  template <int dim, int spacedim>
  template <template <int, int> class DoFHandlerType, typename VectorType>
  void
  MesquiteMeshInterface<dim, spacedim>::update_mesh_displacement_solution(
    VectorType &                         solution,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_displacement_component) const
  {
    update_mesh_displacement_solution(solution,
                                      StaticMappingQ1<dim, spacedim>::mapping,
                                      dof_handler,
                                      first_displacement_component);
  }



  template <int dim, int spacedim>
  template <template <int, int> class DoFHandlerType, typename VectorType>
  void
  MesquiteMeshInterface<dim, spacedim>::update_eulerian_mesh_solution(
    VectorType &                         solution,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_displacement_component) const
  {
    update_eulerian_mesh_solution(solution,
                                  StaticMappingQ1<dim, spacedim>::mapping,
                                  dof_handler,
                                  first_displacement_component);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MappingType,
            template <int, int> class DoFHandlerType,
            typename VectorType>
  void
  MesquiteMeshInterface<dim, spacedim>::update_mesh_displacement_solution(
    VectorType &                         solution,
    const MappingType<dim, spacedim> &   mapping,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_displacement_component) const
  {
    const Triangulation<dim, spacedim> &tria = dof_handler.get_triangulation();
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(&tria == ptr_triangulation, ExcTriangulationMismatch());
    Assert(ptr_mesh_motion_data, ExcNotInitialized());

    VectorDisplacementMapType vertex_displacement_map(dim * tria.n_vertices());
    update_vertex_displacement_map(vertex_displacement_map, tria);

    // Copy the mesh motion field into the appropriate entries in
    // the output vector
    ptr_mesh_motion_data->insert_vertex_displacement_map(
      solution,
      mapping,
      dof_handler,
      vertex_displacement_map,
      first_displacement_component);
  }



  template <int dim, int spacedim>
  template <template <int, int> class MappingType,
            template <int, int> class DoFHandlerType,
            typename VectorType>
  void
  MesquiteMeshInterface<dim, spacedim>::update_eulerian_mesh_solution(
    VectorType &                         solution,
    const MappingType<dim, spacedim> &   mapping,
    const DoFHandlerType<dim, spacedim> &dof_handler,
    const unsigned int                   first_displacement_component) const
  {
    const Triangulation<dim, spacedim> &tria = dof_handler.get_triangulation();
    Assert(ptr_triangulation, ExcNotInitialized());
    Assert(&tria == ptr_triangulation, ExcTriangulationMismatch());
    Assert(ptr_mesh_motion_data, ExcNotInitialized());

    VectorPositionMapType eulerian_vertex_positions(dim * tria.n_vertices());
    update_eulerian_vertex_positions(eulerian_vertex_positions, tria);

    // Copy the mesh motion field into the appropriate entries in
    // the output vector
    ptr_mesh_motion_data->insert_vertex_displacement_map(
      solution,
      mapping,
      dof_handler,
      eulerian_vertex_positions,
      first_displacement_component);
  }


#  ifdef DEAL_II_WITH_VERDICT


  template <int dim, int spacedim>
  Vector<double>
  MesquiteMeshInterface<dim, spacedim>::get_cell_quality() const
  {
    Assert(mesquite_mesh, ExcNotInitialized());
    Assert(ptr_triangulation, ExcNotInitialized());

    // Error tracker
    Mesquite2::MsqError err;

    // Extract all of the vertices
    const unsigned int n_vertices = ptr_triangulation->n_vertices();
    std::vector<Mesquite2::Mesh::VertexHandle> vertices;
    mesquite_mesh->get_all_vertices(vertices, err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));
    Assert(vertices.size() == n_vertices,
           ExcDimensionMismatch(vertices.size(), n_vertices));

    // Retrieve the (current) coordinates of the mesh points
    std::vector<Mesquite2::MsqVertex> vertex_coords(n_vertices);
    mesquite_mesh->vertices_get_coordinates(vertices.data(),
                                            vertex_coords.data(),
                                            n_vertices,
                                            err);
    Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

    // See
    // https://github.com/Kitware/VTK/blob/master/ThirdParty/verdict/vtkverdict/verdict.h.in#L188
    unsigned long metrics_flag = 0;
    if (dim == 2)
      {
        //    metrics_flag += V_QUAD_SHAPE;
        //    metrics_flag += V_QUAD_DISTORTION;
        //    metrics_flag += V_QUAD_AREA;
        metrics_flag += V_QUAD_SCALED_JACOBIAN;
      }
    else if (dim == 3)
      {
        metrics_flag += V_HEX_SCALED_JACOBIAN;
      }
    else
      {
        AssertThrow(false, ExcNotImplemented());
      }

    // Create the metrics for the 2d case, and initialise them
    // since Verdict can't be bothered to do it itself.
    QuadMetricVals quad_metrics;
    quad_metrics.edge_ratio            = std::numeric_limits<double>::min();
    quad_metrics.max_edge_ratio        = std::numeric_limits<double>::min();
    quad_metrics.aspect_ratio          = std::numeric_limits<double>::min();
    quad_metrics.radius_ratio          = std::numeric_limits<double>::min();
    quad_metrics.med_aspect_frobenius  = std::numeric_limits<double>::min();
    quad_metrics.max_aspect_frobenius  = std::numeric_limits<double>::min();
    quad_metrics.skew                  = std::numeric_limits<double>::min();
    quad_metrics.taper                 = std::numeric_limits<double>::min();
    quad_metrics.warpage               = std::numeric_limits<double>::min();
    quad_metrics.area                  = std::numeric_limits<double>::min();
    quad_metrics.stretch               = std::numeric_limits<double>::min();
    quad_metrics.minimum_angle         = std::numeric_limits<double>::min();
    quad_metrics.maximum_angle         = std::numeric_limits<double>::min();
    quad_metrics.oddy                  = std::numeric_limits<double>::min();
    quad_metrics.condition             = std::numeric_limits<double>::min();
    quad_metrics.jacobian              = std::numeric_limits<double>::min();
    quad_metrics.scaled_jacobian       = std::numeric_limits<double>::min();
    quad_metrics.shear                 = std::numeric_limits<double>::min();
    quad_metrics.shape                 = std::numeric_limits<double>::min();
    quad_metrics.relative_size_squared = std::numeric_limits<double>::min();
    quad_metrics.shape_and_size        = std::numeric_limits<double>::min();
    quad_metrics.shear_and_size        = std::numeric_limits<double>::min();
    quad_metrics.distortion            = std::numeric_limits<double>::min();

    // Create the metrics for the 3d case, and initialise them
    // since Verdict can't be bothered to do it itself.
    HexMetricVals hex_metrics;
    hex_metrics.edge_ratio            = std::numeric_limits<double>::min();
    hex_metrics.max_edge_ratio        = std::numeric_limits<double>::min();
    hex_metrics.skew                  = std::numeric_limits<double>::min();
    hex_metrics.taper                 = std::numeric_limits<double>::min();
    hex_metrics.volume                = std::numeric_limits<double>::min();
    hex_metrics.stretch               = std::numeric_limits<double>::min();
    hex_metrics.diagonal              = std::numeric_limits<double>::min();
    hex_metrics.dimension             = std::numeric_limits<double>::min();
    hex_metrics.oddy                  = std::numeric_limits<double>::min();
    hex_metrics.med_aspect_frobenius  = std::numeric_limits<double>::min();
    hex_metrics.max_aspect_frobenius  = std::numeric_limits<double>::min();
    hex_metrics.condition             = std::numeric_limits<double>::min();
    hex_metrics.jacobian              = std::numeric_limits<double>::min();
    hex_metrics.scaled_jacobian       = std::numeric_limits<double>::min();
    hex_metrics.shear                 = std::numeric_limits<double>::min();
    hex_metrics.shape                 = std::numeric_limits<double>::min();
    hex_metrics.relative_size_squared = std::numeric_limits<double>::min();
    hex_metrics.shape_and_size        = std::numeric_limits<double>::min();
    hex_metrics.shear_and_size        = std::numeric_limits<double>::min();
    hex_metrics.distortion            = std::numeric_limits<double>::min();

    Vector<double> quality_metric_per_cell(ptr_triangulation->n_active_cells());

    // We build the element connectivity and vertex positions as stored in the
    // Mesquite mesh. Note that the mesh may or may not have been smoothed at
    // this point; this allows one to perform smoothing, check the quality and
    // potentially smooth again based on the result prior to updating the
    // triangulation vertex positions or Eulerian displacement map.
    std::vector<ElementIndexType> cell_connectivity(
      GeometryInfo<dim>::vertices_per_cell *
      ptr_triangulation->n_active_cells());
    unsigned int cell_count = 0;
    double cell_vertex_coords[GeometryInfo<dim>::vertices_per_cell][meshdim];
    for (auto cell : ptr_triangulation->active_cell_iterators())
      {
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
          {
            // It appears that Mesquite and Verdict use the same vertex ordering
            // scheme, which is different to the lexicographical scheme used
            // by deal.II. We need to account for this.
            const unsigned int idx_vtx = cell->vertex_index(vertex_ordering[v]);
            for (unsigned int d = 0; d < dim; ++d)
              cell_vertex_coords[v][d] = vertex_coords[idx_vtx][d];

            // Check cell quality according to the given metrics
            if (dim == 2)
              {
                v_quad_quality(GeometryInfo<dim>::vertices_per_cell,
                               cell_vertex_coords,
                               metrics_flag,
                               &quad_metrics);
                quality_metric_per_cell[cell_count] =
                  quad_metrics.scaled_jacobian;
              }
            else if (dim == 3)
              {
                v_hex_quality(GeometryInfo<dim>::vertices_per_cell,
                              cell_vertex_coords,
                              metrics_flag,
                              &hex_metrics);
                quality_metric_per_cell[cell_count] =
                  hex_metrics.scaled_jacobian;
              }
            else
              {
                AssertThrow(false, ExcNotImplemented());
              }
          }
        ++cell_count;
      }

    return quality_metric_per_cell;
  }


  template <int dim, int spacedim>
  double
  MesquiteMeshInterface<dim, spacedim>::get_minimum_cell_quality() const
  {
    const Vector<double> cell_quality = get_cell_quality();
    return *std::min_element(cell_quality.begin(), cell_quality.end());
  }


#  endif

} // namespace GridTools


DEAL_II_NAMESPACE_CLOSE

#endif // DEAL_II_WITH_MESQUITE

#endif // dealii_grid_tools_mesquite_h
