// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2023 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

#ifndef dealii_vtk_utilities_h
#define dealii_vtk_utilities_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/vector.h>

#ifdef DEAL_II_WITH_VTK

#  include <vtkDataObject.h>
#  include <vtkDoubleArray.h>
#  include <vtkSmartPointer.h>
#  include <vtkUnstructuredGrid.h>

#  include <string>

#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_VTK
/**
 * Interface to the Visualization Toolkit (VTK).
 *
 * VTK is an open-source, freely available software system for 3D computer
 * graphics, modeling, image processing, volume rendering, scientific
 * visualization, and 2D plotting.

 * It supports a wide variety of visualization algorithms and advanced
 * modeling techniques, and it takes advantage of both threaded and distributed
 * memory parallel processing for speed and scalability, respectively.
 *
 * You can learn more about the VTK library at https://vtk.org/
 */
namespace VTKWrappers
{
  namespace internal
  {
    /**
     * Load a VTK file containing an unstructured grid and return the
     * unstructured grid object.
     *
     * If cleanup is true, overlapping points in the VTK file are merged using
     * VTK's cleaning utilities (only available in VTK >= 9.3).
     *
     * If the read operation fails, an exception is thrown.
     *
     * @param vtk_filename
     * @param cleanup
     * @param relative_tolerance
     * @return vtkUnstructuredGrid*
     */
    vtkSmartPointer<vtkUnstructuredGrid>
    load_vtk_file(const std::string &vtk_filename,
                  const bool         cleanup            = true,
                  const double       relative_tolerance = 0.0);

    /**
     * Convert any VTK dataset that can be represented as an unstructured grid
     * into a @p vtkUnstructuredGrid.
     *
     * Supported input types:
     * - @p vtkUnstructuredGrid → returned as-is (shallow copy)
     * - @p vtkStructuredGrid → hexahedral cells extracted from the structured
     *   connectivity
     * - @p vtkRectilinearGrid / @p vtkStructuredPoints → hexahedral cells
     *   reconstructed from grid dimensions
     * - @p vtkPolyData → triangles, lines, and polygons converted to their
     *   unstructured grid equivalents
     *
     * If the input is not a supported dataset type, an exception is thrown.
     * Point and cell data are copied to the output.
     *
     * This function is used internally by @p load_vtk_file and @p write_vtk
     * to normalize all VTK data to a single unstructured grid representation.
     *
     * @param data_object The VTK data object to convert.
     * @return A smart pointer to a newly created @p vtkUnstructuredGrid.
     */
    vtkSmartPointer<vtkUnstructuredGrid>
    convert_to_unstructured_grid(
      const vtkSmartPointer<vtkDataObject> &data_object);

    /**
     * Write a VTK file to disk from any VTK dataset that can be converted to
     * an unstructured grid.
     *
     * The output format is determined by the file extension:
     * - `.vtu` → VTK XML UnstructuredGrid format
     * - `.vtk` → legacy VTK UnstructuredGrid format
     *
     * Supported input types are the same as those handled by
     * @p convert_to_unstructured_grid: unstructured grids, structured grids,
     * rectilinear grids, structured points, and polydata. Each type is
     * converted to an unstructured grid internally before writing.
     *
     * If the write operation fails, an exception is thrown.
     *
     * @param[in] vtk_filename Path to the output file. The extension
     *         determines the format: `.vtu` for XML format, `.vtk` for legacy
     *         format.
     * @param[in] data_object  The VTK data object to write. Must
     *         point to a valid dataset that can be converted to an
     *         unstructured grid.
     */
    void
    write_vtk(const std::string                    &vtk_filename,
              const vtkSmartPointer<vtkDataObject> &data_object);
  } // namespace internal

  /**
   * Convert from a deal.II Point to a VTK double array.
   *
   * @tparam dim Dimension of the point
   * @param [in] p An input deal.II Point<dim>
   * @return A VTK smart pointer to the data array.
   */
  template <int dim>
  inline vtkSmartPointer<vtkDoubleArray>
  dealii_point_to_vtk_array(const dealii::Point<dim> &p);

  /**
   * Convert a VTK unstructured grid to a deal.II triangulation.
   *
   * This function translates the points and cells in @p unstructured_grid into
   * a serial deal.II Triangulation.
   *
   * The optional string parameters specify names of VTK cell-data arrays
   * to read and how they are interpreted:
   * - @p material_id_field: when non-empty, the named VTK cell-data scalar
   *   array is read and used to set per-codim-0-cell `material_id` values in
   *   the resulting Triangulation.
   * - @p boundary_id_field: when non-empty, the named VTK cell-data scalar
   *   array is read and used to set `boundary_id`s for codim-1 subcells (faces)
   *   reconstructed from VTK cells whose dimension is smaller than @p dim.
   * - @p manifold_id_field: when non-empty, the named VTK cell-data scalar
   *   array is read and used to set `manifold_id` for both codim-0 cells and
   *   appended subcells.
   *
   * Negative values read from VTK are interpreted as deal.II defaults
   * (e.g. `numbers::flat_manifold_id`) to avoid unsigned wrap when mapping
   * into deal.II types.
   */
  template <int dim, int spacedim>
  void
  unstructured_grid_to_dealii_triangulation(
    const vtkUnstructuredGrid    &unstructured_grid,
    Triangulation<dim, spacedim> &tria,
    const std::string            &material_id_field = "",
    const std::string            &boundary_id_field = "",
    const std::string            &manifold_id_field = "");

  /**
   * Convert a serial deal.II triangulation to a VTK unstructured grid.
   *
   * This function exports the active cells of @p tria to a newly allocated
   * VTK unstructured grid and returns it.
   *
   * The optional string parameters specify names of VTK cell-data arrays to
   * create and the export behaviour:
   * - @p material_id_field: when non-empty, a VTK cell-data scalar array with
   *   this name is created and filled with per-codim-0-cell `material_id`
   *   values. Appended subcells (boundary faces) receive a value of 0.
   * - @p boundary_id_field: when non-empty, a VTK cell-data scalar array with
   *   this name is created. Only appended codim-1 subcells carry meaningful
   *   boundary id values; codim-0 cells receive 0.
   * - @p manifold_id_field: when non-empty, a VTK cell-data scalar array with
   *   this name is created. Codim-0 cells are written with
   *   `numbers::flat_manifold_id` (i.e. -1) when they have the default
   *   manifold; non-flat manifold ids for appended subcells are written
   *   alongside boundary ids.
   */
  template <int dim, int spacedim>
  vtkSmartPointer<vtkUnstructuredGrid>
  dealii_triangulation_to_unstructured_grid(
    const Triangulation<dim, spacedim> &tria,
    const std::string                  &material_id_field = "",
    const std::string                  &boundary_id_field = "",
    const std::string                  &manifold_id_field = "");

  /**
   * Write a serial deal.II triangulation to a VTK file.
   *
   * This function first converts @p tria to a VTK unstructured grid through
   * dealii_triangulation_to_unstructured_grid() and then writes it to disk.
   *
   * The optional string parameters specify names of VTK cell-data arrays to
   * create and fill with deal.II ids:
   * - @p material_id_field: when non-empty, write per-codim-0-cell
   *   `material_id` values;
   * - @p boundary_id_field: when non-empty, write `boundary_id` values on
   *   appended codim-1 subcells;
   * - @p manifold_id_field: when non-empty, write `manifold_id` values on
   *   appended codim-1 subcells and `numbers::flat_manifold_id` on codim-0
   *   cells.
   *
   * The output format is selected from the file extension: `.vtk` for legacy
   * VTK and `.vtu` for VTK XML.
   */
  template <int dim, int spacedim>
  void
  write_vtk(const std::string                  &vtk_filename,
            const Triangulation<dim, spacedim> &tria,
            const std::string                  &material_id_field = "",
            const std::string                  &boundary_id_field = "",
            const std::string                  &manifold_id_field = "");

  /**
   * Read a VTK mesh file and populate a deal.II Triangulation.
   *
   * This function reads the mesh from the specified VTK file and fills the
   * given Triangulation object. If cleanup is true, overlapping points in the
   * VTK file are merged using VTK's cleaning utilities.
   *
   * The optional string parameters are forwarded to
   * unstructured_grid_to_dealii_triangulation() and control whether VTK
   * cell-data arrays are interpreted as deal.II ids:
   * - @p material_id_field: per-cell material ids for codim-0 cells;
   * - @p boundary_id_field: per-subcell boundary ids for codim-1 cells;
   * - @p manifold_id_field: manifold ids for both codim-0 and codim-1 cells.
   *
   * @param vtk_filename The name of the input VTK file.
   * @param tria The Triangulation object to populate.
   * @param material_id_field Optional VTK cell-data scalar array name to use
   * for material ids (default: empty string, disabled).
   * @param boundary_id_field Optional VTK cell-data scalar array name to use
   * for boundary ids (default: empty string, disabled).
   * @param manifold_id_field Optional VTK cell-data scalar array name to use
   * for manifold ids (default: empty string, disabled).
   * @param cleanup If true, merge overlapping points in the VTK file (default: true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  template <int dim, int spacedim>
  void
  read_tria(const std::string            &vtk_filename,
            Triangulation<dim, spacedim> &tria,
            const bool                    cleanup            = true,
            const double                  relative_tolerance = 0.0,
            const std::string            &material_id_field  = "",
            const std::string            &boundary_id_field  = "",
            const std::string            &manifold_id_field  = "");

  /**
   * Read cell data (scalar or vector) from a VTK file and store it in the
   * output vector.
   *
   * This function reads the specified cell data array (scalar or vector) from
   * the given VTK file and stores it in the provided output vector. For vector
   * data, all components are stored in row-major order (cell0_comp0,
   * cell0_comp1, ..., cell1_comp0, ...).
   *
   * @param vtk_filename The name of the input VTK file.
   * @param cell_data_name The name of the cell data array to read.
   * @param output_vector The vector to store the cell data values.
   * @param cleanup If true, merge overlapping points in the VTK file (default:
   * true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  void
  read_cell_data(const std::string &vtk_filename,
                 const std::string &cell_data_name,
                 Vector<double>    &output_vector,
                 const bool         cleanup            = true,
                 const double       relative_tolerance = 0.0);

  /**
   * Read vertex data from a VTK file and store it in the output vector.
   *
   * This function reads the specified vertex data array (scalar or vector) from
   * the given VTK file and stores it in the provided output vector.
   *
   * For vector data, all components are stored in row-major order
   * (vertex0_comp0, vertex0_comp1, ..., vertex1_comp0, ...).
   *
   * @param vtk_filename The name of the input VTK file.
   * @param vertex_data_name The name of the vertex data array to read.
   * @param output_vector The vector to store the vertex data values.
   * @param cleanup If true, merge overlapping points in the VTK file (default:
   * true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  void
  read_vertex_data(const std::string &vtk_filename,
                   const std::string &vertex_data_name,
                   Vector<double>    &output_vector,
                   const bool         cleanup            = true,
                   const double       relative_tolerance = 0.0);

  /**
   * Read all field data from a VTK file and store it in the output vector.
   *
   * This function reads all field data arrays (scalar or vector, cell or point
   * data) from the given VTK file and stores it in the provided output vector.
   *
   * The data is output in the following way:
   * - first all vertex data (point data) in the order they are found in the
   *   VTK file, with all components stored in row-major order (vertex0_comp0,
   *   vertex0_comp1, ..., vertex1_comp0, ...)
   * - then all cell data (cell data) in the order they are found in the VTK
   * file, with all components stored in row-major order (cell0_comp0,
   * cell0_comp1, ...).
   *
   * This is equivalent to calling read_vertex_data() for each vertex data
   * field, and then read_cell_data() for each cell data field, and
   * concatenating the resulting vectors in a single long vector.
   *
   * @param vtk_filename The name of the input VTK file.
   * @param output_vector The vector to store the vertex data values.
   * @param cleanup If true, merge overlapping points in the VTK file (default:
   * true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  void
  read_all_data(const std::string &vtk_filename,
                Vector<double>    &output_vector,
                const bool         cleanup            = true,
                const double       relative_tolerance = 0.0);

  /**
   * Create a FiniteElement representation for data stored in a VTK file.
   *
   * This function inspects the specified VTK file and builds a FiniteElement
   * (typically an FESystem) that can store the data fields contained in the
   * file. It returns a pair consisting of a std::unique_ptr to the constructed
   * FiniteElement and a vector of field names corresponding to the VTK data
   * arrays.
   *
   * Only meshes with a single cell type are supported (either all quads/hexes
   * or all triangles/tetrahedra). Mixed meshes are not supported. The mesh type
   * is inferred from the first cell encountered in the VTK file.
   *
   * Mapping rules:
   * - Vertex data (VTK point data):
   *   - quad/hex meshes: represented by FE_Q(1) blocks;
   *   - simplex meshes: represented by FE_SimplexP(1) blocks.
   * - Cell data (VTK cell data):
   *   - quad/hex meshes: represented by FE_DGQ(0) blocks;
   *   - simplex meshes: represented by FE_SimplexDGP(0) blocks.
   *
   *   Vector-valued cell or vertex data are represented as an FESystem of the
   *   scalar cell element repeated for each component.
   *
   * @param vtk_filename The path to the input VTK file.
   */
  template <int dim, int spacedim>
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
            std::vector<std::string>>
  vtk_to_finite_element(const std::string &vtk_filename);

  /**
   * Translate a vtk data file (obtained through read_all_data() above) to a
   * dealii vector type, associated with the given DoFHandler object.
   *
   * Prerequisites:
   * - The input triangulation must be a serial triangulation, obtained through
   *   read_tria() above.
   * - The DoFHandler must have been initialized with the finite element
   *   obtained through the vtk_to_finite_element() method, and degrees of
   *   freedom must have been already distributed.
   * - The triangulation associated with the DoFHandler doest not need to be
   *   same as the input triangulation, but it must be obtained from it. It may
   *   be, for example, a parallel::fullydistributed::Triangulation obtained by
   *   distributing the serial triangulation.
   * - The input data must refer to the serial triangulation (obtained through
   *   the read_vtk() method), and must have been obtained by calling the
   *   read_all_data() method above.
   *
   * The DoFHandler object may be serial or parallel, and you may renumber it to
   * your liking before calling this method. An example usage is the following:
   * @code
   * // Read serial triangulation from VTK file Triangulation<dim, spacedim>
   * serial_tria; VTKWrappers::read_tria(vtk_filename, serial_tria);
   *
   * // Read all data from VTK file Vector<double> serial_data;
   * VTKWrappers::read_all_data(vtk_filename, serial_data);
   *
   * // Read finite element from VTK file auto [fe, data_names] =
   * VTKWrappers::vtk_to_finite_element<dim, spacedim>(vtk_filename);
   *
   * // Split serial triangulation into a parallel one
   * parallel::fullydistributed::Triangulation<dim, spacedim> parallel_tria;
   * parallel_tria.copy_triangulation(serial_tria); // use default partitioner
   *
   * // Setup DoFHandler on parallel triangulation DoFHandler<dim, spacedim>
   * dof_handler(parallel_tria); dof_handler.distribute_dofs(*fe); // Optionally
   * renumber dofs here
   * ...
   * // Map serial data to distributed vector
   * LinearAlgebra::distributed::Vector<double> distributed_data;
   * VTKWrappers::data_to_dealii_vector( serial_tria, serial_data, dof_handler,
   * distributed_data);
   * @endcode
   *
   * This function exists to allow manipulation of the DoFHandler (e.g., through
   * block wise renumbering) before mapping the (serial) data to (a possibly
   * distributed) vector.
   *
   * @tparam dim
   * @tparam spacedim
   * @tparam VectorType
   * @param serial_tria
   * @param data
   * @param dh
   * @param output_vector
   */
  template <int dim, int spacedim, typename VectorType>
  void
  data_to_dealii_vector(const Triangulation<dim, spacedim> &serial_tria,
                        const Vector<double>               &data,
                        const DoFHandler<dim, spacedim>    &dh,
                        VectorType                         &output_vector);

  /**
   * Read a VTK mesh and all data fields into a DoFHandler and output vector.
   *
   * This function reads the mesh from the specified VTK file, populates the
   * Triangulation associated to the given DoFHandler, and queries all cell and
   * vertex data fields. For each data field, a suitable FESystem is constructed
   * (using FE_DGQ for cell data and FE_Q for vertex data, with the correct
   * number of components). DoFs are distributed and renumbered block-wise. All
   * data is read into the output_vector, and the names of the fields are stored
   * in data_names.
   *
   * @param vtk_filename The name of the input VTK file.
   * @param dof_handler The DoFHandler to distribute DoFs on the mesh.
   * @param output_vector The vector to store all data field values.
   * @param data_names The vector to store the names of all data fields found in
   * the VTK file.
   * @param material_id_field Optional VTK cell-data scalar array name to use
   * for material ids when reconstructing the Triangulation (default: empty
   * string, disabled).
   * @param boundary_id_field Optional VTK cell-data scalar array name to use
   * for boundary ids when reconstructing the Triangulation (default: empty
   * string, disabled).
   * @param manifold_id_field Optional VTK cell-data scalar array name to use
   * for manifold ids when reconstructing the Triangulation (default: empty
   * string, disabled).
   * @param cleanup If true, merge overlapping points in the VTK file (default:
   * true).
   * @param relative_tolerance Relative tolerance used when merging points via
   * VTK's cleaning utilities (default: 0).
   */
  template <int dim, int spacedim>
  void
  read_vtk(const std::string         &vtk_filename,
           DoFHandler<dim, spacedim> &dof_handler,
           Vector<double>            &output_vector,
           std::vector<std::string>  &data_names,
           const bool                 cleanup            = true,
           const double               relative_tolerance = 0.0,
           const std::string         &material_id_field  = "",
           const std::string         &boundary_id_field  = "",
           const std::string         &manifold_id_field  = "");

#  ifndef DOXYGEN
  // Template implementations

  template <int dim>
  inline vtkSmartPointer<vtkDoubleArray>
  dealii_point_to_vtk_array(const dealii::Point<dim> &p)
  {
    vtkSmartPointer<vtkDoubleArray> p_vtk =
      vtkSmartPointer<vtkDoubleArray>::New();

    p_vtk->SetNumberOfComponents(dim);
    p_vtk->SetNumberOfTuples(1);

    for (int d = 0; d < dim; ++d)
      p_vtk->FillComponent(d, p[d]);

    return p_vtk;
  }



  template <int dim, int spacedim, typename VectorType>
  void
  data_to_dealii_vector(const Triangulation<dim, spacedim> &serial_tria,
                        const Vector<double>               &data,
                        const DoFHandler<dim, spacedim>    &dh,
                        VectorType                         &output_vector)
  {
    AssertDimension(dh.n_dofs(), output_vector.size());
    const auto &fe = dh.get_fe();

    const auto dist_to_serial_vertices =
      GridTools::parallel_to_serial_vertex_indices(serial_tria,
                                                   dh.get_triangulation());

    const auto &locally_owned_dofs = dh.locally_owned_dofs();

    types::global_dof_index dofs_offset        = 0;
    unsigned int            vertex_comp_offset = 0;
    unsigned int            cell_comp_offset   = 0;
    for (unsigned int field = 0; field < fe.n_blocks(); ++field)
      {
        const auto        &field_fe = fe.base_element(field);
        const unsigned int n_comps  = field_fe.n_components();
        if (field_fe.n_dofs_per_vertex() > 0)
          {
            // This is a vertex data field
            const types::global_dof_index n_local_dofs =
              n_comps * serial_tria.n_vertices();
            for (const auto &cell : dh.active_cell_iterators())
              if (cell->is_locally_owned())
                for (const auto v : cell->vertex_indices())
                  {
                    const types::global_dof_index serial_vertex_index =
                      dist_to_serial_vertices[cell->vertex_index(v)];
                    if (serial_vertex_index != numbers::invalid_unsigned_int)
                      for (unsigned int c = 0; c < n_comps; ++c)
                        {
                          const types::global_dof_index dof_index =
                            cell->vertex_dof_index(v, vertex_comp_offset + c);
                          Assert(locally_owned_dofs.is_element(dof_index),
                                 ExcInternalError());
                          output_vector[dof_index] =
                            data[dofs_offset + n_comps * serial_vertex_index +
                                 c];
                        }
                  }
            dofs_offset += n_local_dofs;
            vertex_comp_offset += n_comps;
          }
        else if (field_fe.template n_dofs_per_object<dim>() > 0)
          {
            // this is a cell data field
            const types::global_dof_index n_local_dofs =
              n_comps * serial_tria.n_global_active_cells();

            // Assumption: serial and parallel meshes have the same ordering of
            // cells.
            auto serial_cell   = serial_tria.begin_active();
            auto parallel_cell = dh.begin_active();
            for (; parallel_cell != dh.end(); ++parallel_cell)
              if (parallel_cell->is_locally_owned())
                {
                  // Advance serial cell until we reach the same cell index of
                  // the parallel cell
                  while (serial_cell->id() < parallel_cell->id())
                    ++serial_cell;
                  const auto serial_cell_index =
                    serial_cell->global_active_cell_index();
                  for (unsigned int c = 0; c < n_comps; ++c)
                    {
                      const types::global_dof_index dof_index =
                        parallel_cell->dof_index(cell_comp_offset + c);
                      Assert(locally_owned_dofs.is_element(dof_index),
                             ExcInternalError());
                      output_vector[dof_index] =
                        data[dofs_offset + n_comps * serial_cell_index + c];
                    }
                }
            dofs_offset += n_local_dofs;
            cell_comp_offset += n_comps;
          }
      }
  }
#  endif

} // namespace VTKWrappers

#endif

DEAL_II_NAMESPACE_CLOSE
#endif
