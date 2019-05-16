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

/**
 * A collection of utilities for the Mesquite tests
 */

#include <deal.II/base/function.h>
#include <deal.II/base/point.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/mapping_q_eulerian.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools_mesquite.h>
#include <deal.II/grid/tria.h>

#include <deal.II/hp/dof_handler.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/block_vector_base.h>
#include <deal.II/lac/vector.h>

#include <deal.II/numerics/data_out.h>

#include <Mesquite_all_headers.hpp>

#include <fstream>
#include <string>
#include <vector>

#include "../tests.h"


void smooth_mesquite_in_place_arraymesh(Triangulation<2> &tria)
{
  const unsigned int dim = 2;

  // ---------------------
  // Convert triangulation to a format that Mesquite understands

  // Mesh vertex coordinates
  // Note: This has a stride of dim between vertices
  std::vector<double>            coords(dim * tria.n_vertices());
  const std::vector<Point<dim>> &vertices = tria.get_vertices();
  for (unsigned int v = 0; v < vertices.size(); ++v)
    for (unsigned int d = 0; d < dim; ++d)
      {
        const unsigned int c = dim * v + d; // Flattened coordinate index
        Assert(c < coords.size(), ExcIndexRange(c, 0, coords.size()));
        coords[c] = vertices[v][d];
      }

  // Element connectivity
  // Note: This has a stride of n_vertices_per_cell between elements
  std::vector<unsigned long> cells(GeometryInfo<dim>::vertices_per_cell *
                                   tria.n_active_cells());
  unsigned int               cell_count = 0;
  for (auto cell : tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          const unsigned int c =
            GeometryInfo<dim>::vertices_per_cell * cell_count +
            v; // Flattened global vertex index
          Assert(c < cells.size(), ExcIndexRange(c, 0, cells.size()));
          cells[c] = cell->vertex_index(v);
        }
      ++cell_count;
    }

  // Constraints: All boundary vertices are fixed.
  // Note: This has a stride of 1 between vertices
  std::vector<int> constraints(tria.n_vertices(), 0);
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
              Assert(c < constraints.size(),
                     ExcIndexRange(c, 0, constraints.size()));
              constraints[c] = 1;
            }
        }
    }

  // Pass the above into Mesquite
  // Note: This Mesquite2::ArrayMesh object understands both
  //       2d and 3d entities, while a Mesquite2::MeshImpl can
  //       only work in 3d.
  Mesquite2::ArrayMesh mesh(
    dim,                   // specify a mesh (dim coordinate values per vertex)
    tria.n_vertices(),     // number of vertices
    coords.data(),         // vertex coordinates
    constraints.data(),    // constraint flags
    tria.n_active_cells(), // number of elements
    Mesquite2::QUADRILATERAL, // element type
    cells.data()              // element connectivity
  );

  // ---------------------

  // Error tracker
  Mesquite2::MsqError err;

  // Surface to constrain the 2d elements to
  Mesquite2::PlanarDomain domain(Mesquite2::PlanarDomain::XY);

  // Build a view of the domain
  Mesquite2::MeshDomainAssoc mesh_and_domain(&mesh, &domain);

  // Improve the mesh
  Mesquite2::ShapeImprover shape_wrapper;
  shape_wrapper.run_instructions(&mesh_and_domain, err);
  Assert(err == Mesquite2::MsqError::NO_ERROR, ExcMesquiteError(err));

  // ---------------------
  // Move vertices in the triangulation as required
  std::vector<bool> vertex_touched(tria.n_vertices(), false);
  cell_count = 0;
  for (auto cell : tria.active_cell_iterators())
    {
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        {
          const unsigned int idx = cell->vertex_index(v);
          Assert(idx < vertex_touched.size(),
                 ExcIndexRange(idx, 0, vertex_touched.size()));
          if (!vertex_touched[idx])
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const unsigned int c =
                    dim * idx + d; // Flattened coordinate index
                  Assert(c < coords.size(), ExcIndexRange(c, 0, coords.size()));
                  cell->vertex(v)[d] = coords[c];
                }

              vertex_touched[idx] = true;
            }
        }
      ++cell_count;
    }
}


template <int dim, typename RangeNumberType = double>
class TransformationSinusoid : public Function<dim, RangeNumberType>
{
public:
  TransformationSinusoid(bool         allow_perturbation,
                         double       perturbation_size = 0.75,
                         unsigned int n_components      = dim,
                         unsigned int first_u_component = 0)
    : Function<dim, RangeNumberType>(n_components)
    , allow_perturbation(allow_perturbation)
    , perturbation_size(perturbation_size)
    , first_u_component(first_u_component)
  {}

  virtual ~TransformationSinusoid()
  {}

  // Add this operator so that it can be used in GridTools::transform
  Point<dim>
  operator()(const Point<dim> &in) const
  {
    Point<dim> pt;
    pt[0] = in(0);
    pt[1] = in(1) + std::sin(in(0) / 5.0 * 3.14159);
    if (dim == 3)
      pt[2] = in(2);

    if (allow_perturbation)
      {
        for (unsigned int d = 0; d < dim; ++d)
          pt[d] += perturbation_size * random_value();
      }

    return pt;
  }

  virtual void
  vector_value(const Point<dim> &p, Vector<RangeNumberType> &return_value) const
  {
    AssertThrow(dim > 1, ExcNotImplemented());

    const Point<dim> pt = this->operator()(p);

    Assert(return_value.size() >= (first_u_component + dim),
           ExcInternalError());
    for (unsigned int d = 0; d < dim; ++d)
      return_value[d + first_u_component] = pt[d];
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
  const bool         allow_perturbation;
  const double       perturbation_size;
  const unsigned int first_u_component;
};



enum class OutputFormat
{
  vtk,
  vtu
};



template <int dim>
void
output_mesh(const Triangulation<dim> &triangulation,
            std::string               filename,
            const bool                write_to_deallog = false,
            const enum OutputFormat   format           = OutputFormat::vtk)
{
  GridOut grid_out;

  if (format == OutputFormat::vtk)
    {
      filename += ".vtk";
      std::ofstream out(filename);
      grid_out.write_vtk(triangulation, out);

      if (write_to_deallog)
        grid_out.write_vtk(triangulation, deallog.get_file_stream());
    }
  else if (format == OutputFormat::vtu)
    {
      filename += ".vtu";
      std::ofstream out(filename);
      grid_out.write_vtu(triangulation, out);

      if (write_to_deallog)
        grid_out.write_vtu(triangulation, deallog.get_file_stream());
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }
}


template <int dim, int spacedim, template <int, int> class MeshType>
void
output_mesh_vertex_numbering(const MeshType<dim, spacedim> &mesh,
                             std::string                    filename)
{
  if (spacedim == 3)
    return;

  const std::string filename_orig(filename);

  filename += ".gnuplot";
  std::ofstream out(filename);

  std::map<types::global_dof_index, Point<dim>> vertex_positions;
  for (typename MeshType<dim, spacedim>::active_cell_iterator cell =
         mesh.begin_active();
       cell != mesh.end();
       ++cell)
    {
      for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; ++v)
        vertex_positions[cell->vertex_index(v)] = cell->vertex(v);
    }

  out
    << "set terminal png size 400,410 enhanced font \"Helvetica,8\""
    << std::endl
    << "set output \"" << filename_orig << ".png\"" << std::endl
    << "set size square" << std::endl
    << "set view equal xy" << std::endl
    << "unset xtics" << std::endl
    << "unset ytics" << std::endl
    << "plot '-' using 1:2 with lines notitle, '-' with labels point pt 2 offset 1,1 notitle"
    << std::endl;
  GridOut().write_gnuplot(mesh.get_triangulation(), out);
  out << "e" << std::endl;
  DoFTools::write_gnuplot_dof_support_point_info(out, vertex_positions);
  out << "e" << std::endl;
}


// Undo the perturbations of the boundary vertices
template <int dim, typename Vector>
void
fix_eulerian_boundary_vertex_map(Vector &                  eulerian_vertex_map,
                                 const Triangulation<dim> &triangulation,
                                 const Vector &eulerian_vertex_map_exact,
                                 const Vector &eulerian_vertex_map_perturbed)
{
  Assert(eulerian_vertex_map_exact.size() == dim * triangulation.n_vertices(),
         ExcDimensionMismatch(eulerian_vertex_map_exact.size(),
                              dim * triangulation.n_vertices()));
  Assert(eulerian_vertex_map_perturbed.size() ==
           dim * triangulation.n_vertices(),
         ExcDimensionMismatch(eulerian_vertex_map_perturbed.size(),
                              dim * triangulation.n_vertices()));

  // First assume the perturbed state...
  typename GridTools::internal::LocalVector<Vector>::type local_vector(
    eulerian_vertex_map_perturbed);

  // ... then perform the corrections
  std::vector<bool> vertex_touched(triangulation.n_vertices(), false);
  for (typename Triangulation<dim>::active_cell_iterator cell =
         triangulation.begin_active();
       cell != triangulation.end();
       ++cell)
    for (unsigned int f = 0; f < GeometryInfo<dim>::faces_per_cell; ++f)
      if (cell->face(f)->at_boundary())
        for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face; ++v)
          if (vertex_touched[cell->face(f)->vertex_index(v)] == false)
            {
              for (unsigned int d = 0; d < dim; ++d)
                {
                  const unsigned int idx =
                    GridTools::vertex_map_index_from_global_vertex_index<dim>(
                      cell->face(f)->vertex_index(v), d);
                  local_vector[idx] = eulerian_vertex_map_exact[idx];
                }

              vertex_touched[cell->face(f)->vertex_index(v)] = true;
            }

  eulerian_vertex_map = local_vector;
}



template <int dim, int spacedim, typename Number>
void
output_mesh(const Triangulation<dim, spacedim> &triangulation,
            const Vector<Number> &              vertex_displacement_map,
            std::string                         filename,
            const bool                          write_to_deallog = false,
            const enum OutputFormat             format = OutputFormat::vtk)
{
  Triangulation<dim, spacedim> tria_eulerian;
  tria_eulerian.copy_triangulation(triangulation);
  GridTools::move_triangulation_vertices(tria_eulerian,
                                         vertex_displacement_map);
  output_mesh(tria_eulerian, filename, write_to_deallog, format);
}



template <int dim, int spacedim>
unsigned int
get_tensor_degree(const DoFHandler<dim, spacedim> &dof_handler)
{
  return dof_handler.get_fe().tensor_degree();
}



template <int dim, int spacedim>
unsigned int
get_tensor_degree(const hp::DoFHandler<dim, spacedim> &dof_handler)
{
  unsigned int tensor_degree = 1;
  for (unsigned int f = 0; f < dof_handler.get_fe().size(); ++f)
    {
      tensor_degree =
        std::max(tensor_degree, dof_handler.get_fe()[f].tensor_degree());
    }
  return tensor_degree;
}



template <int dim,
          int spacedim,
          template <int, int> class DoFHandlerType,
          typename Vector>
void
output_mesh(
  const DoFHandlerType<dim, spacedim> &dof_handler,
  const Vector &                       solution,
  const unsigned int                   first_u_dof,
  std::string                          filename,
  const bool                           write_to_deallog = false,
  const enum OutputFormat              format           = OutputFormat::vtk,
  typename std::enable_if<!IsBlockVector<Vector>::value>::type * = nullptr)
{
  const unsigned int degree = get_tensor_degree(dof_handler);
  DataOut<dim, DoFHandlerType<dim, spacedim>> data_out;
  data_out.attach_dof_handler(dof_handler);

  const unsigned int n_dof_components =
    dof_handler.get_fe_collection().n_components();
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      n_dof_components, DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> solution_name(n_dof_components, "displacement");
  for (unsigned int d = first_u_dof; d < dim + first_u_dof; ++d)
    data_component_interpretation[d] =
      DataComponentInterpretation::component_is_part_of_vector;
  for (unsigned int i = 0; i < n_dof_components; ++i)
    {
      if (i >= first_u_dof && i < dim + first_u_dof)
        continue;
      solution_name[i] = "dummy_" + Utilities::int_to_string(i);
    }

  data_out.add_data_vector(
    solution,
    solution_name,
    DataOut<dim, DoFHandlerType<dim, spacedim>>::type_dof_data,
    data_component_interpretation);
  if (n_dof_components == dim)
    {
      // MappingQEulerian always expects displacement
      // DoFs as the first 3 components, and in some of
      // the tests we will purposely ensure that this is
      // not the case.
      const MappingQEulerian<dim> mapping(degree, dof_handler, solution);
      data_out.build_patches(mapping, degree);
    }
  else
    {
      data_out.build_patches(degree);
    }

  if (format == OutputFormat::vtk)
    {
      filename += ".vtk";
      std::ofstream out(filename);
      data_out.write_vtk(out);

      if (write_to_deallog)
        data_out.write_vtk(deallog.get_file_stream());
    }
  else if (format == OutputFormat::vtu)
    {
      filename += ".vtu";
      std::ofstream out(filename);
      data_out.write_vtu(out);

      if (write_to_deallog)
        data_out.write_vtu(deallog.get_file_stream());
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }
}



template <int dim, int spacedim, template <int, int> class DoFHandlerType>
struct EulerianMapping;

template <int dim, int spacedim>
struct EulerianMapping<dim, spacedim, DoFHandler>
{
  static const bool supported = true;

  template <typename BlockVector>
  static void
  build_patches(DataOut<dim, DoFHandler<dim, spacedim>> &data_out,
                const unsigned int                       degree,
                const DoFHandler<dim, spacedim> &        dof_handler,
                const BlockVector &                      solution,
                const unsigned int                       u_block)
  {
    const Vector<double>        solution_u(solution.block(u_block));
    const MappingQEulerian<dim> mapping(degree, dof_handler, solution_u);
    data_out.build_patches(mapping, degree);
  }
};

template <int dim, int spacedim>
struct EulerianMapping<dim, spacedim, hp::DoFHandler>
{
  static const bool supported = false;

  template <typename BlockVector>
  static void
  build_patches(DataOut<dim, hp::DoFHandler<dim, spacedim>> &data_out,
                const unsigned int                           degree,
                const hp::DoFHandler<dim, spacedim> &        dof_handler,
                const BlockVector &                          solution,
                const unsigned int                           u_block)
  {
    AssertThrow(false, ExcNotImplemented());
  }
};



template <int dim,
          int spacedim,
          template <int, int> class DoFHandlerType,
          typename BlockVector>
void
output_mesh(
  const DoFHandlerType<dim, spacedim> &dof_handler,
  const BlockVector &                  solution,
  const unsigned int                   u_block,
  const unsigned int                   first_u_dof,
  std::string                          filename,
  const bool                           write_to_deallog = false,
  const enum OutputFormat              format           = OutputFormat::vtk,
  typename std::enable_if<IsBlockVector<BlockVector>::value>::type * = nullptr)
{
  const unsigned int degree = get_tensor_degree(dof_handler);
  DataOut<dim, DoFHandlerType<dim, spacedim>> data_out;
  data_out.attach_dof_handler(dof_handler);

  const unsigned int n_dof_components =
    dof_handler.get_fe_collection().n_components();
  std::vector<DataComponentInterpretation::DataComponentInterpretation>
    data_component_interpretation(
      n_dof_components, DataComponentInterpretation::component_is_scalar);
  std::vector<std::string> solution_name(n_dof_components, "displacement");
  for (unsigned int d = first_u_dof; d < dim + first_u_dof; ++d)
    data_component_interpretation[d] =
      DataComponentInterpretation::component_is_part_of_vector;
  for (unsigned int i = 0; i < n_dof_components; ++i)
    {
      if (i >= first_u_dof && i < dim + first_u_dof)
        continue;
      solution_name[i] = "dummy_" + Utilities::int_to_string(i);
    }

  data_out.add_data_vector(
    solution,
    solution_name,
    DataOut<dim, DoFHandlerType<dim, spacedim>>::type_dof_data,
    data_component_interpretation);

  // MappingQEulerian always expects displacement
  // DoFs as the first 3 components, and in some of
  // the tests we will purposely ensure that this is
  // not the case.
  typedef EulerianMapping<dim, spacedim, DoFHandlerType> EulerianMappingType;
  if (n_dof_components == dim && EulerianMappingType::supported)
    EulerianMappingType::build_patches(
      data_out, degree, dof_handler, solution, u_block);
  else
    data_out.build_patches(degree);

  if (format == OutputFormat::vtk)
    {
      filename += ".vtk";
      std::ofstream out(filename);
      data_out.write_vtk(out);

      if (write_to_deallog)
        data_out.write_vtk(deallog.get_file_stream());
    }
  else if (format == OutputFormat::vtu)
    {
      filename += ".vtu";
      std::ofstream out(filename);
      data_out.write_vtu(out);

      if (write_to_deallog)
        data_out.write_vtu(deallog.get_file_stream());
    }
  else
    {
      AssertThrow(false, ExcNotImplemented());
    }
}



template <int dim,
          int spacedim,
          template <int, int> class DoFHandlerType,
          typename BlockVector>
void
output_mesh(
  const DoFHandlerType<dim, spacedim> &dof_handler,
  const BlockVector &                  solution,
  const std::vector<IndexSet>          locally_relevant_partitioning,
  const unsigned int                   u_block,
  const unsigned int                   first_u_dof,
  std::string                          filename,
  const bool                           write_to_deallog = false,
  const enum OutputFormat              format           = OutputFormat::vtk,
  typename std::enable_if<IsBlockVector<BlockVector>::value>::type * = nullptr)
{
  BlockVector locally_relevant_solution(locally_relevant_partitioning);
  locally_relevant_solution = solution;
  output_mesh(dof_handler,
              locally_relevant_solution,
              u_block,
              first_u_dof,
              filename,
              write_to_deallog,
              format);
}

template <int dim,
          int spacedim,
          template <int, int> class DoFHandlerType,
          typename BlockVector>
void
output_mesh(
  const DoFHandlerType<dim, spacedim> &dof_handler,
  const BlockVector &                  solution,
  const std::vector<IndexSet>          locally_relevant_partitioning,
  const unsigned int                   u_block,
  const unsigned int                   first_u_dof,
  std::string                          filename,
  const MPI_Comm &                     mpi_communicator,
  const bool                           write_to_deallog = false,
  const enum OutputFormat              format           = OutputFormat::vtk,
  typename std::enable_if<IsBlockVector<BlockVector>::value>::type * = nullptr)
{
  const unsigned int this_mpi_process =
    Utilities::MPI::this_mpi_process(mpi_communicator);
  const unsigned int n_mpi_processes =
    Utilities::MPI::n_mpi_processes(mpi_communicator);

  auto get_filename_proc = [](const std::string       filename,
                              const unsigned int      proccess,
                              const enum OutputFormat format,
                              const bool              with_suffix = false) {
    std::string out = filename + "." + Utilities::int_to_string(proccess);
    if (with_suffix && format == OutputFormat::vtk)
      out += ".vtk";
    else if (with_suffix && format == OutputFormat::vtu)
      out += ".vtu";
    return out;
  };
  const std::string filename_proc =
    get_filename_proc(filename, this_mpi_process, format);

  output_mesh(dof_handler,
              solution,
              locally_relevant_partitioning,
              u_block,
              first_u_dof,
              filename_proc,
              write_to_deallog,
              format);

  if (format == OutputFormat::vtu && this_mpi_process == 0)
    {
      std::vector<std::string> piece_names;
      for (unsigned int p = 0; p < n_mpi_processes; ++p)
        piece_names.emplace_back(get_filename_proc(filename, p, format, true));
      std::ofstream pvtu_output(filename + ".pvtu");
      DataOut<dim>().write_pvtu_record(pvtu_output, piece_names);
    }
}
