// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------
#include <deal.II/vtk/utilities.h>

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

#ifdef DEAL_II_WITH_VTK
#  include <deal.II/distributed/fully_distributed_tria.h>

#  include <deal.II/dofs/dof_renumbering.h>

#  include <deal.II/fe/fe.h>
#  include <deal.II/fe/fe_dgq.h>
#  include <deal.II/fe/fe_nothing.h>
#  include <deal.II/fe/fe_q.h>
#  include <deal.II/fe/fe_simplex_p.h>
#  include <deal.II/fe/fe_system.h>

#  include <deal.II/grid/grid_in.h>
#  include <deal.II/grid/grid_out.h>
#  include <deal.II/grid/grid_tools.h>
#  include <deal.II/grid/tria_accessor.h>
#  include <deal.II/grid/tria_description.h>
#  include <deal.II/grid/tria_iterator.h>

#  include <deal.II/lac/vector.h>

// Make sure that the VTK version macros are available.
#  include <vtkVersion.h> // Do not convert for module purposes

// VTK_VERSION_CHECK is defined by the header above for 9.3.0 and above, but
// we provide a fallback older versions.
#  ifndef VTK_VERSION_CHECK
#    define VTK_VERSION_CHECK(major, minor, build) \
      (10000000000ULL * (major) + 100000000ULL * (minor) + (build))
#  endif

// Normalize the version macro name to cover both the quick (>=9.3) and
// legacy (<9.3) headers.
#  if defined(VTK_VERSION_NUMBER_QUICK)
#    define DEAL_II_VTK_VERSION_NUMBER VTK_VERSION_NUMBER_QUICK
#  elif defined(VTK_VERSION_NUMBER)
#    define DEAL_II_VTK_VERSION_NUMBER VTK_VERSION_NUMBER
#  else
#    define DEAL_II_VTK_VERSION_NUMBER 0
#  endif

#  include <vtkAppendDataSets.h>
#  include <vtkCellData.h>
#  if DEAL_II_VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 3, 0)
#    include <vtkCleanUnstructuredGrid.h>
#    include <vtkGenericDataObjectReader.h>
#    include <vtkXMLGenericDataObjectReader.h>
#  endif
#  include <vtkCellType.h>
#  include <vtkDataArray.h>
#  include <vtkDataSet.h>
#  include <vtkDataSetReader.h>
#  include <vtkFieldData.h>
#  include <vtkGenericCell.h>
#  include <vtkHexahedron.h>
#  include <vtkIntArray.h>
#  include <vtkLine.h>
#  include <vtkPointData.h>
#  include <vtkPoints.h>
#  include <vtkPolyData.h>
#  include <vtkPolygon.h>
#  include <vtkRectilinearGrid.h>
#  include <vtkSmartPointer.h>
#  include <vtkStructuredGrid.h>
#  include <vtkStructuredPoints.h>
#  include <vtkUnstructuredGrid.h>
#  include <vtkUnstructuredGridWriter.h>
#  include <vtkXMLUnstructuredGridWriter.h>

#  include <fstream>
#  include <stdexcept>

#endif

DEAL_II_NAMESPACE_OPEN

#ifdef DEAL_II_WITH_VTK
namespace VTKWrappers
{
  namespace internal
  {
    vtkSmartPointer<vtkUnstructuredGrid>
    convert_to_unstructured_grid(
      const vtkSmartPointer<vtkDataObject> &data_object)
    {
      auto *dataset = dynamic_cast<vtkDataSet *>(data_object.Get());
      AssertThrow(dataset,
                  ExcMessage(
                    std::string(
                      "VTK file does not contain a supported dataset type: ") +
                    data_object->GetClassName()));

      auto append_filter = vtkSmartPointer<vtkAppendDataSets>::New();
      append_filter->SetOutputDataSetType(VTK_UNSTRUCTURED_GRID);
      append_filter->AddInputData(dataset);
      append_filter->Update();

      auto out = vtkSmartPointer<vtkUnstructuredGrid>::New();
      out->ShallowCopy(append_filter->GetOutput());

      AssertThrow(out->GetNumberOfPoints() == dataset->GetNumberOfPoints() ||
                    dataset->GetNumberOfPoints() == 0,
                  ExcMessage("VTK dataset could not be converted to an "
                             "unstructured grid without losing points."));

      AssertThrow(out->GetNumberOfCells() == dataset->GetNumberOfCells() ||
                    dataset->GetNumberOfCells() == 0,
                  ExcMessage("VTK dataset could not be converted to an "
                             "unstructured grid without losing cells."));

      out->GetPointData()->PassData(dataset->GetPointData());
      out->GetCellData()->PassData(dataset->GetCellData());
      out->GetFieldData()->PassData(dataset->GetFieldData());

      return out;
    }



    vtkSmartPointer<vtkUnstructuredGrid>
    load_vtk_file(const std::string &vtk_filename,
                  const bool         cleanup,
                  const double       relative_tolerance)
    {
      // check that the file exists
      std::ifstream file(vtk_filename);
      AssertThrow(file.good(),
                  ExcMessage("VTK file not found: " + vtk_filename));

      vtkSmartPointer<vtkDataObject> data_object;
      const auto                     dot_pos = vtk_filename.find_last_of('.');
      const std::string              ext =
        (dot_pos == std::string::npos ? "" : vtk_filename.substr(dot_pos + 1));

#  if DEAL_II_VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 3, 0)
      if (ext == "vtu")
        {
          auto reader = vtkSmartPointer<vtkXMLGenericDataObjectReader>::New();
          reader->SetFileName(vtk_filename.c_str());
          reader->Update();
          data_object = reader->GetOutput();
        }
      else if (ext == "vtk")
        {
          // Read the legacy VTK format without assuming a specific dataset
          // type; the result is normalized to an unstructured grid below.
          auto reader = vtkSmartPointer<vtkGenericDataObjectReader>::New();
          reader->SetFileName(vtk_filename.c_str());
          reader->Update();
          data_object = reader->GetOutput();
        }
      else
        AssertThrow(false,
                    ExcMessage("Unsupported file extension '" + ext +
                               "'. Use '.vtu' for VTK XML format or '.vtk' "
                               "for legacy VTK format."));
#  else
      AssertThrow(ext == "vtk",
                  ExcMessage("Unsupported file extension '" + ext +
                             "'. With VTK < 9.3 only '.vtk' legacy files can "
                             "be read."));
      {
        auto reader = vtkSmartPointer<vtkDataSetReader>::New();
        reader->SetFileName(vtk_filename.c_str());
        reader->Update();
        data_object = reader->GetOutput();
      }
#  endif

      auto out = convert_to_unstructured_grid(data_object);

#  if DEAL_II_VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 3, 0)
      if (cleanup)
        {
          auto cleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();
          cleaner->SetToleranceIsAbsolute(false);
          cleaner->SetTolerance(relative_tolerance);
          cleaner->SetInputData(out);
          cleaner->Update();
          out->ShallowCopy(cleaner->GetOutput());
        }
#  else
      (void)cleanup;            // avoid unused variable warning
      (void)relative_tolerance; // avoid unused variable warning
      deallog << "VTK version < 9.3: skipping cleanup step." << std::endl;
#  endif

      AssertThrow(out, ExcMessage("Failed to read VTK file: " + vtk_filename));
      return out;
    }



    void
    write_vtk(const std::string                    &vtk_filename,
              const vtkSmartPointer<vtkDataObject> &data_object)
    {
      // Convert any supported dataset to unstructured grid
      auto grid = convert_to_unstructured_grid(data_object);

      // Determine output format from file extension
      const auto dot_pos = vtk_filename.find_last_of('.');
      Assert(dot_pos != std::string::npos,
             ExcMessage("Cannot determine output format: no file extension "
                        "found in filename '" +
                        vtk_filename + "'."));
      const std::string ext = vtk_filename.substr(dot_pos + 1);

      if (ext == "vtu")
        {
          auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
          writer->SetFileName(vtk_filename.c_str());
          writer->SetInputData(grid);
          writer->Write();
        }
      else if (ext == "vtk")
        {
          auto writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
          writer->SetFileName(vtk_filename.c_str());
          writer->SetInputData(grid);
          writer->Write();
        }
      else
        {
          AssertThrow(false,
                      ExcMessage("Unsupported file extension '" + ext +
                                 "'. Use '.vtu' for VTK XML format or '.vtk' "
                                 "for legacy VTK format."));
        }
    }

  } // namespace internal



  template <int dim, int spacedim>
  void
  unstructured_grid_to_dealii_triangulation(
    const vtkUnstructuredGrid    &unstructured_grid,
    Triangulation<dim, spacedim> &tria,
    const std::string            &material_id_field,
    const std::string            &boundary_id_field,
    const std::string            &manifold_id_field)
  {
    auto &grid = const_cast<vtkUnstructuredGrid &>(unstructured_grid);

    auto get_cell_data_array =
      [&grid](const std::string &name) -> vtkDataArray * {
      if (name.empty())
        return nullptr;

      vtkCellData *cell_data = grid.GetCellData();
      if (cell_data == nullptr)
        return nullptr;

      vtkDataArray *data_array = cell_data->GetArray(name.c_str());
      if (data_array != nullptr)
        AssertThrow(data_array->GetNumberOfComponents() == 1,
                    ExcMessage("The VTK cell data array '" + name +
                               "' must be scalar."));

      return data_array;
    };

    // Get observing pointers to the arrays that store
    // material/boundary/manifold ids.
    vtkDataArray *material_ids = get_cell_data_array(material_id_field);
    vtkDataArray *boundary_ids = get_cell_data_array(boundary_id_field);
    vtkDataArray *manifold_ids = get_cell_data_array(manifold_id_field);

    auto material_id_from_vtk = [](vtkDataArray   &material_ids,
                                   const vtkIdType vtk_id) {
      const double material_value = material_ids.GetComponent(vtk_id, 0);
      if (material_value < 0)
        return types::material_id(0);

      return static_cast<types::material_id>(material_value);
    };

    // Read points
    vtkPoints                   *vtk_points = grid.GetPoints();
    const vtkIdType              n_points   = vtk_points->GetNumberOfPoints();
    std::vector<Point<spacedim>> points(n_points);
    for (vtkIdType i = 0; i < n_points; ++i)
      {
        std::array<double, 3> coords = {{0, 0, 0}};
        vtk_points->GetPoint(i, coords.data());
        for (unsigned int d = 0; d < spacedim; ++d)
          points[i][d] = coords[d];
        for (unsigned int d = spacedim; d < 3; ++d)
          AssertThrow(
            coords[d] == 0.0,
            ExcMessage(
              "VTK grid has non-zero coordinate in unused dimension."));
      }

    // Read cells
    std::vector<CellData<dim>> cells;
    SubCellData                subcell_data;
    const vtkIdType            n_cells = grid.GetNumberOfCells();
    std::vector<bool>          boundary_vertex_ids_present(n_points, false);
    std::vector<bool>          manifold_vertex_ids_present(n_points, false);
    std::vector<types::boundary_id> boundary_vertex_ids(
      n_points, numbers::internal_face_boundary_id);
    std::vector<types::manifold_id> manifold_vertex_ids(
      n_points, numbers::flat_manifold_id);
    for (vtkIdType vtk_id = 0; vtk_id < n_cells; ++vtk_id)
      {
        vtkCell *cell = grid.GetCell(vtk_id);
        if (cell->GetCellDimension() < dim)
          {
            if (boundary_ids == nullptr && manifold_ids == nullptr)
              continue;

            auto set_subcell_ids = [vtk_id, boundary_ids, manifold_ids](
                                     auto &cell_data) {
              if (boundary_ids != nullptr)
                {
                  double bval = boundary_ids->GetComponent(vtk_id, 0);
                  if (bval < 0)
                    cell_data.boundary_id = numbers::internal_face_boundary_id;
                  else
                    cell_data.boundary_id =
                      static_cast<types::boundary_id>(bval);
                }
              else if (manifold_ids != nullptr)
                cell_data.boundary_id = numbers::internal_face_boundary_id;

              if (manifold_ids != nullptr)
                {
                  double mval = manifold_ids->GetComponent(vtk_id, 0);
                  if (mval < 0)
                    cell_data.manifold_id = numbers::flat_manifold_id;
                  else
                    cell_data.manifold_id =
                      static_cast<types::manifold_id>(mval);
                }
            };

            if constexpr (dim == 1)
              {
                if (cell->GetCellType() == VTK_VERTEX)
                  {
                    AssertThrow(
                      cell->GetNumberOfPoints() == 1,
                      ExcMessage(
                        "Only vertex subcells with 1 point are supported."));

                    const vtkIdType vertex_index = cell->GetPointId(0);
                    AssertIndexRange(vertex_index, n_points);

                    if (boundary_ids != nullptr)
                      {
                        const double boundary_value =
                          boundary_ids->GetComponent(vtk_id, 0);
                        if (boundary_value >= 0)
                          {
                            boundary_vertex_ids_present[vertex_index] = true;
                            boundary_vertex_ids[vertex_index] =
                              static_cast<types::boundary_id>(boundary_value);
                          }
                      }

                    if (manifold_ids != nullptr)
                      {
                        const double manifold_value =
                          manifold_ids->GetComponent(vtk_id, 0);
                        manifold_vertex_ids_present[vertex_index] = true;
                        manifold_vertex_ids[vertex_index] =
                          (manifold_value < 0 ?
                             numbers::flat_manifold_id :
                             static_cast<types::manifold_id>(manifold_value));
                      }
                  }
              }
            else if constexpr (dim == 2)
              {
                if (cell->GetCellType() == VTK_LINE)
                  {
                    AssertThrow(
                      cell->GetNumberOfPoints() == 2,
                      ExcMessage(
                        "Only line subcells with 2 points are supported."));

                    CellData<1> cell_data(2);
                    for (unsigned int j = 0; j < 2; ++j)
                      cell_data.vertices[j] = cell->GetPointId(j);

                    set_subcell_ids(cell_data);
                    subcell_data.boundary_lines.push_back(cell_data);
                  }
              }
            else if constexpr (dim == 3)
              {
                if (cell->GetCellType() == VTK_LINE)
                  {
                    AssertThrow(
                      cell->GetNumberOfPoints() == 2,
                      ExcMessage(
                        "Only line subcells with 2 points are supported."));

                    CellData<1> cell_data(2);
                    for (unsigned int j = 0; j < 2; ++j)
                      cell_data.vertices[j] = cell->GetPointId(j);

                    set_subcell_ids(cell_data);
                    subcell_data.boundary_lines.push_back(cell_data);
                  }
                else if (cell->GetCellType() == VTK_TRIANGLE ||
                         cell->GetCellType() == VTK_QUAD)
                  {
                    AssertThrow(
                      cell->GetNumberOfPoints() == 3 ||
                        cell->GetNumberOfPoints() == 4,
                      ExcMessage(
                        "Only triangle and quad face subcells are supported."));

                    CellData<2> cell_data(cell->GetNumberOfPoints());
                    for (unsigned int j = 0; j < cell->GetNumberOfPoints(); ++j)
                      cell_data.vertices[j] = cell->GetPointId(j);

                    // VTK and deal.II use different vertex ordering for
                    // quad faces; match deal.II ordering by swapping the
                    // last two vertices (consistent with quad cell handling
                    // above).
                    if (cell->GetCellType() == VTK_QUAD)
                      std::swap(cell_data.vertices[2], cell_data.vertices[3]);

                    set_subcell_ids(cell_data);
                    subcell_data.boundary_quads.push_back(cell_data);
                  }
              }

            continue;
          }

        AssertThrow(cell->GetCellDimension() == dim,
                    ExcMessage("Unsupported VTK cell dimension."));

        if constexpr (dim == 1)
          {
            if (cell->GetCellType() != VTKCellType::VTK_LINE)
              AssertThrow(false,
                          ExcMessage(
                            "Unsupported cell type in 1D VTK file: only "
                            "VTK_LINE is supported."));
            AssertThrow(cell->GetNumberOfPoints() == 2,
                        ExcMessage(
                          "Only line cells with 2 points are supported."));
            CellData<1> cell_data(2);
            for (unsigned int j = 0; j < 2; ++j)
              cell_data.vertices[j] = cell->GetPointId(j);
            cell_data.material_id = 0;
            if (material_ids != nullptr)
              cell_data.material_id =
                material_id_from_vtk(*material_ids, vtk_id);
            if (boundary_ids != nullptr)
              {
                const double boundary_value =
                  boundary_ids->GetComponent(vtk_id, 0);
                if (boundary_value >= 0)
                  cell_data.boundary_id =
                    static_cast<types::boundary_id>(boundary_value);
              }
            if (manifold_ids != nullptr)
              {
                double mval = manifold_ids->GetComponent(vtk_id, 0);
                if (mval < 0)
                  cell_data.manifold_id = numbers::flat_manifold_id;
                else
                  cell_data.manifold_id = static_cast<types::manifold_id>(mval);
              }
            cells.push_back(cell_data);
          }
        else if constexpr (dim == 2)
          {
            if (cell->GetCellType() == VTKCellType::VTK_QUAD)
              {
                AssertThrow(cell->GetNumberOfPoints() == 4,
                            ExcMessage(
                              "Only quad cells with 4 points are supported."));
                CellData<2> cell_data(4);
                for (unsigned int j = 0; j < 4; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                std::swap(cell_data.vertices[2], cell_data.vertices[3]);
                cell_data.material_id = 0;
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTKCellType::VTK_TRIANGLE)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 3,
                  ExcMessage(
                    "Only triangle cells with 3 points are supported."));
                CellData<2> cell_data(3);
                for (unsigned int j = 0; j < 3; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else
              AssertThrow(false,
                          ExcMessage(
                            "Unsupported cell type in 2D VTK file: only "
                            "VTK_QUAD and VTK_TRIANGLE are supported."));
          }
        else if constexpr (dim == 3)
          {
            if (cell->GetCellType() == VTKCellType::VTK_HEXAHEDRON)
              {
                AssertThrow(cell->GetNumberOfPoints() == 8,
                            ExcMessage(
                              "Only hex cells with 8 points are supported."));
                CellData<3> cell_data(8);
                for (unsigned int j = 0; j < 8; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                // Numbering of vertices in VTK files is different from
                // deal.II
                std::swap(cell_data.vertices[2], cell_data.vertices[3]);
                std::swap(cell_data.vertices[6], cell_data.vertices[7]);
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTKCellType::VTK_TETRA)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 4,
                  ExcMessage(
                    "Only tetrahedron cells with 4 points are supported."));
                CellData<3> cell_data(4);
                for (unsigned int j = 0; j < 4; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTKCellType::VTK_WEDGE)
              {
                AssertThrow(cell->GetNumberOfPoints() == 6,
                            ExcMessage(
                              "Only prism cells with 6 points are supported."));
                CellData<3> cell_data(6);
                for (unsigned int j = 0; j < 6; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTKCellType::VTK_PYRAMID)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 5,
                  ExcMessage(
                    "Only pyramid cells with 5 points are supported."));
                CellData<3> cell_data(5);
                for (unsigned int j = 0; j < 5; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                if (material_ids != nullptr)
                  cell_data.material_id =
                    material_id_from_vtk(*material_ids, vtk_id);
                if (manifold_ids != nullptr)
                  {
                    double mval = manifold_ids->GetComponent(vtk_id, 0);
                    if (mval < 0)
                      cell_data.manifold_id = numbers::flat_manifold_id;
                    else
                      cell_data.manifold_id =
                        static_cast<types::manifold_id>(mval);
                  }
                cells.push_back(cell_data);
              }
            else
              AssertThrow(
                false,
                ExcMessage(
                  "Unsupported cell type in 3D VTK file: only "
                  "VTK_HEXAHEDRON, VTK_TETRA, VTK_WEDGE, and VTK_PYRAMID are supported."));
          }
        else
          {
            AssertThrow(false, ExcMessage("Unsupported dimension."));
          }
      }

    // Create triangulation
    tria.create_triangulation(points, cells, subcell_data);

    if constexpr (dim == 1)
      for (const auto &cell : tria.active_cell_iterators())
        for (unsigned int f = 0; f < cell->n_faces(); ++f)
          if (cell->face(f)->at_boundary())
            {
              const unsigned int vertex_index = cell->face(f)->vertex_index(0);

              if (boundary_vertex_ids_present[vertex_index])
                cell->face(f)->set_boundary_id(
                  boundary_vertex_ids[vertex_index]);
              if (manifold_vertex_ids_present[vertex_index])
                cell->face(f)->set_manifold_id(
                  manifold_vertex_ids[vertex_index]);
            }
  }



  template <int dim, int spacedim>
  vtkSmartPointer<vtkUnstructuredGrid>
  dealii_triangulation_to_unstructured_grid(
    const Triangulation<dim, spacedim> &tria,
    const std::string                  &material_id_field,
    const std::string                  &boundary_id_field,
    const std::string                  &manifold_id_field)
  {
    auto grid   = vtkSmartPointer<vtkUnstructuredGrid>::New();
    auto points = vtkSmartPointer<vtkPoints>::New();

    points->SetNumberOfPoints(tria.n_vertices());
    for (unsigned int i = 0; i < tria.n_vertices(); ++i)
      {
        std::array<double, 3> coords = {{0.0, 0.0, 0.0}};
        for (unsigned int d = 0; d < spacedim; ++d)
          coords[d] = tria.get_vertices()[i][d];
        points->SetPoint(i, coords.data());
      }
    grid->SetPoints(points);

    // Prepare optional VTK arrays for cell data. We will append values in
    // the same order that we insert cells: first full (codim 0) cells, then
    // any subcells (boundary faces) we append afterwards.
    vtkSmartPointer<vtkIntArray> material_array;
    vtkSmartPointer<vtkIntArray> boundary_array;
    vtkSmartPointer<vtkIntArray> manifold_array;

    const bool output_material = !material_id_field.empty();
    const bool output_boundary = !boundary_id_field.empty();
    const bool output_manifold = !manifold_id_field.empty();

    if (output_material)
      {
        material_array = vtkSmartPointer<vtkIntArray>::New();
        material_array->SetName(material_id_field.c_str());
        material_array->SetNumberOfComponents(1);
      }
    if (output_boundary)
      {
        boundary_array = vtkSmartPointer<vtkIntArray>::New();
        boundary_array->SetName(boundary_id_field.c_str());
        boundary_array->SetNumberOfComponents(1);
      }
    if (output_manifold)
      {
        manifold_array = vtkSmartPointer<vtkIntArray>::New();
        manifold_array->SetName(manifold_id_field.c_str());
        manifold_array->SetNumberOfComponents(1);
      }

    for (const auto &cell : tria.active_cell_iterators())
      {
        const unsigned int     n_vertices = cell->n_vertices();
        std::vector<vtkIdType> point_ids(n_vertices);
        for (unsigned int i = 0; i < n_vertices; ++i)
          point_ids[i] = static_cast<vtkIdType>(cell->vertex_index(i));

        int vtk_cell_type = -1;
        if constexpr (dim == 1)
          {
            AssertThrow(n_vertices == 2,
                        ExcMessage("Unsupported 1D cell with != 2 vertices."));
            vtk_cell_type = VTKCellType::VTK_LINE;
          }
        else if constexpr (dim == 2)
          {
            if (n_vertices == 4)
              {
                vtk_cell_type = VTKCellType::VTK_QUAD;
                std::swap(point_ids[2], point_ids[3]);
              }
            else if (n_vertices == 3)
              vtk_cell_type = VTKCellType::VTK_TRIANGLE;
            else
              AssertThrow(false,
                          ExcMessage("Unsupported 2D cell type: only "
                                     "quads and triangles are supported."));
          }
        else if constexpr (dim == 3)
          {
            if (n_vertices == 8)
              {
                vtk_cell_type = VTKCellType::VTK_HEXAHEDRON;
                std::swap(point_ids[2], point_ids[3]);
                std::swap(point_ids[6], point_ids[7]);
              }
            else if (n_vertices == 4)
              vtk_cell_type = VTKCellType::VTK_TETRA;
            else if (n_vertices == 6)
              vtk_cell_type = VTKCellType::VTK_WEDGE;
            else if (n_vertices == 5)
              vtk_cell_type = VTKCellType::VTK_PYRAMID;
            else
              AssertThrow(false,
                          ExcMessage("Unsupported 3D cell type: only "
                                     "hexes, tets, wedges and pyramids are "
                                     "supported."));
          }
        else
          AssertThrow(false, ExcMessage("Unsupported dimension."));

        grid->InsertNextCell(vtk_cell_type, n_vertices, point_ids.data());

        // For codim-0 cells write material/manifold values and zero for
        // boundary id, since boundary ids live on codim-1 subcells (vertices
        // in 1d, edges in 2d, and faces in 3d).
        if (output_material)
          material_array->InsertNextValue(
            static_cast<int>(cell->material_id()));
        if (output_boundary)
          boundary_array->InsertNextValue(0);
        if (output_manifold)
          manifold_array->InsertNextValue(
            static_cast<int>(numbers::flat_manifold_id));
      }

    // Append boundary faces as separate VTK cells when requested. We only
    // append codim-1 subcells (faces / edges) for which either a non-zero
    // boundary id or a non-flat manifold id exists, depending on the fields
    // requested by the caller.
    for (const auto &cell : tria.active_cell_iterators())
      for (unsigned int f = 0; f < cell->n_faces(); ++f)
        if (cell->face(f)->at_boundary())
          {
            const types::boundary_id face_bid = cell->face(f)->boundary_id();
            const types::manifold_id face_mid = cell->face(f)->manifold_id();

            const bool include_by_boundary = output_boundary && (face_bid != 0);
            const bool include_by_manifold =
              output_manifold && (face_mid != numbers::flat_manifold_id);

            if (!(include_by_boundary || include_by_manifold))
              continue;

            const unsigned int     nfv = cell->face(f)->n_vertices();
            std::vector<vtkIdType> face_point_ids(nfv);
            for (unsigned int v = 0; v < nfv; ++v)
              face_point_ids[v] =
                static_cast<vtkIdType>(cell->face(f)->vertex_index(v));

            int vtk_face_type = -1;
            if constexpr (dim == 1)
              {
                Assert(nfv == 1, ExcInternalError());
                vtk_face_type = VTK_VERTEX;
              }
            else if constexpr (dim == 2)
              vtk_face_type = VTK_LINE;
            else if constexpr (dim == 3)
              {
                if (nfv == 3)
                  vtk_face_type = VTK_TRIANGLE;
                else if (nfv == 4)
                  {
                    vtk_face_type = VTK_QUAD;
                    std::swap(face_point_ids[2], face_point_ids[3]);
                  }
                else
                  DEAL_II_ASSERT_UNREACHABLE();
              }

            grid->InsertNextCell(vtk_face_type, nfv, face_point_ids.data());

            if (output_material)
              material_array->InsertNextValue(
                static_cast<int>(cell->material_id()));
            if (output_boundary)
              boundary_array->InsertNextValue(static_cast<int>(face_bid));
            if (output_manifold)
              manifold_array->InsertNextValue(static_cast<int>(face_mid));
          }

    // Attach arrays to cell data if created.
    if (output_material)
      grid->GetCellData()->AddArray(material_array);
    if (output_boundary)
      grid->GetCellData()->AddArray(boundary_array);
    if (output_manifold)
      grid->GetCellData()->AddArray(manifold_array);

    return grid;
  }



  template <int dim, int spacedim>
  void
  write_vtk(const std::string                  &vtk_filename,
            const Triangulation<dim, spacedim> &tria,
            const std::string                  &material_id_field,
            const std::string                  &boundary_id_field,
            const std::string                  &manifold_id_field)
  {
    const auto grid = dealii_triangulation_to_unstructured_grid(
      tria, material_id_field, boundary_id_field, manifold_id_field);
    internal::write_vtk(vtk_filename, grid);
  }



  template <int dim, int spacedim>
  void
  read_tria(const std::string            &vtk_filename,
            Triangulation<dim, spacedim> &tria,
            const bool                    cleanup,
            const double                  relative_tolerance,
            const std::string            &material_id_field,
            const std::string            &boundary_id_field,
            const std::string            &manifold_id_field)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
      internal::load_vtk_file(vtk_filename, cleanup, relative_tolerance);
    unstructured_grid_to_dealii_triangulation(
      *grid, tria, material_id_field, boundary_id_field, manifold_id_field);
  }



  void
  read_cell_data(const std::string &vtk_filename,
                 const std::string &cell_data_name,
                 Vector<double>    &output_vector,
                 const bool         cleanup,
                 const double       relative_tolerance)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
      internal::load_vtk_file(vtk_filename, cleanup, relative_tolerance);
    vtkDataArray *data_array =
      grid->GetCellData()->GetArray(cell_data_name.c_str());
    AssertThrow(data_array,
                ExcMessage("Cell data array '" + cell_data_name +
                           "' not found in VTK file: " + vtk_filename));
    vtkIdType n_tuples     = data_array->GetNumberOfTuples();
    int       n_components = data_array->GetNumberOfComponents();
    output_vector.reinit(n_tuples * n_components);
    for (vtkIdType i = 0; i < n_tuples; ++i)
      for (int j = 0; j < n_components; ++j)
        output_vector[i * n_components + j] = data_array->GetComponent(i, j);
  }



  void
  read_vertex_data(const std::string &vtk_filename,
                   const std::string &vertex_data_name,
                   Vector<double>    &output_vector,
                   const bool         cleanup,
                   const double       relative_tolerance)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
      internal::load_vtk_file(vtk_filename, cleanup, relative_tolerance);
    vtkDataArray *data_array =
      grid->GetPointData()->GetArray(vertex_data_name.c_str());
    AssertThrow(data_array,
                ExcMessage("Point data array '" + vertex_data_name +
                           "' not found in VTK file: " + vtk_filename));
    vtkIdType n_tuples     = data_array->GetNumberOfTuples();
    int       n_components = data_array->GetNumberOfComponents();
    output_vector.reinit(n_tuples * n_components);
    for (vtkIdType i = 0; i < n_tuples; ++i)
      for (int j = 0; j < n_components; ++j)
        output_vector[i * n_components + j] = data_array->GetComponent(i, j);
  }



  void
  read_all_data(const std::string &vtk_filename,
                Vector<double>    &output_vector,
                const bool         cleanup,
                const double       relative_tolerance)
  {
    vtkSmartPointer<vtkUnstructuredGrid> grid =
      internal::load_vtk_file(vtk_filename, cleanup, relative_tolerance);

    std::vector<double> data;

    vtkPointData *point_data = grid->GetPointData();
    if (point_data)
      {
        for (int i = 0; i < point_data->GetNumberOfArrays(); ++i)
          {
            vtkDataArray *data_array = point_data->GetArray(i);
            if (!data_array)
              continue;
            vtkIdType    n_tuples     = data_array->GetNumberOfTuples();
            int          n_components = data_array->GetNumberOfComponents();
            unsigned int current_size = data.size();
            data.resize(current_size + n_tuples * n_components, 0.0);
            for (vtkIdType tuple_idx = 0; tuple_idx < n_tuples; ++tuple_idx)
              for (int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                data[current_size + tuple_idx * n_components + comp_idx] =
                  data_array->GetComponent(tuple_idx, comp_idx);
          }
      }

    vtkCellData *cell_data = grid->GetCellData();
    if (cell_data)
      {
        for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
          {
            vtkDataArray *data_array = cell_data->GetArray(i);
            if (!data_array)
              continue;
            vtkIdType    n_tuples     = data_array->GetNumberOfTuples();
            int          n_components = data_array->GetNumberOfComponents();
            unsigned int current_size = data.size();
            data.resize(current_size + n_tuples * n_components, true);
            for (vtkIdType tuple_idx = 0; tuple_idx < n_tuples; ++tuple_idx)
              for (int comp_idx = 0; comp_idx < n_components; ++comp_idx)
                data[current_size + tuple_idx * n_components + comp_idx] =
                  data_array->GetComponent(tuple_idx, comp_idx);
          }
      }
    output_vector.reinit(data.size());
    std::copy(data.begin(), data.end(), output_vector.begin());
  }



  template <int dim, int spacedim>
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
            std::vector<std::string>>
  vtk_to_finite_element(const std::string &vtk_filename)
  {
    std::vector<std::string> data_names;

    vtkSmartPointer<vtkUnstructuredGrid> grid =
      internal::load_vtk_file(vtk_filename, false, 0.0);

    vtkCellData  *cell_data  = grid->GetCellData();
    vtkPointData *point_data = grid->GetPointData();

    std::vector<std::shared_ptr<FiniteElement<dim, spacedim>>> fe_collection;
    std::vector<unsigned int> n_components_collection;

    bool            is_simplex    = false;
    const vtkIdType n_cells_check = grid->GetNumberOfCells();
    for (vtkIdType i = 0; i < n_cells_check; ++i)
      {
        vtkCell *cell = grid->GetCell(i);
        if (!cell)
          continue;
        const int cell_type = cell->GetCellType();
        if constexpr (dim == 2)
          {
            if (cell_type == VTKCellType::VTK_TRIANGLE)
              {
                is_simplex = true;
                break;
              }
          }
        else if constexpr (dim == 3)
          {
            if (cell_type == VTKCellType::VTK_TETRA)
              {
                is_simplex = true;
                break;
              }
          }
        else
          {
            is_simplex = false;
            break;
          }
      }

    // Query point data fields
    for (int i = 0; i < point_data->GetNumberOfArrays(); ++i)
      {
        vtkDataArray *arr = point_data->GetArray(i);
        if (!arr)
          continue;
        std::string name   = arr->GetName();
        int         n_comp = arr->GetNumberOfComponents();

        if (is_simplex)
          {
            if (n_comp == 1)
              fe_collection.push_back(
                std::make_shared<FE_SimplexP<dim, spacedim>>(1));
            else
              // Use FESystem for vector fields
              fe_collection.push_back(std::make_shared<FESystem<dim, spacedim>>(
                FE_SimplexP<dim, spacedim>(1), n_comp));
          }
        else
          {
            if (n_comp == 1)
              fe_collection.push_back(std::make_shared<FE_Q<dim, spacedim>>(1));
            else
              // Use FESystem for vector fields
              fe_collection.push_back(std::make_shared<FESystem<dim, spacedim>>(
                FE_Q<dim, spacedim>(1), n_comp));
          }
        n_components_collection.push_back(n_comp);
        data_names.push_back(name);
      }

    // Query cell data fields
    for (int i = 0; i < cell_data->GetNumberOfArrays(); ++i)
      {
        vtkDataArray *arr = cell_data->GetArray(i);
        if (!arr)
          continue;
        std::string name   = arr->GetName();
        int         n_comp = arr->GetNumberOfComponents();
        if (is_simplex)
          {
            if (n_comp == 1)
              fe_collection.push_back(
                std::make_shared<FE_SimplexDGP<dim, spacedim>>(0));
            else
              // Use FESystem for vector fields
              fe_collection.push_back(std::make_shared<FESystem<dim, spacedim>>(
                FE_SimplexDGP<dim, spacedim>(0), n_comp));
          }
        else
          {
            if (n_comp == 1)
              fe_collection.push_back(
                std::make_shared<FE_DGQ<dim, spacedim>>(0));
            else
              // Use FESystem for vector fields
              fe_collection.push_back(std::make_shared<FESystem<dim, spacedim>>(
                FE_DGQ<dim, spacedim>(0), n_comp));
          }
        n_components_collection.push_back(n_comp);
        data_names.push_back(name);
      }


    // Build a FESystem with all fields
    std::vector<const FiniteElement<dim, spacedim> *> fe_ptrs;
    std::vector<unsigned int>                         multiplicities;
    for (const auto &fe : fe_collection)
      {
        fe_ptrs.push_back(fe.get());
        multiplicities.push_back(1);
      }
    if (fe_ptrs.empty())
      return std::make_pair(std::make_unique<FE_Nothing<dim, spacedim>>(),
                            std::vector<std::string>());
    else
      {
        return std::make_pair(
          std::make_unique<FESystem<dim, spacedim>>(fe_ptrs, multiplicities),
          data_names);
      }
  }


  template <int dim, int spacedim>
  void
  read_vtk(const std::string         &vtk_filename,
           DoFHandler<dim, spacedim> &dof_handler,
           Vector<double>            &output_vector,
           std::vector<std::string>  &data_names,
           const bool                 cleanup,
           const double               relative_tolerance,
           const std::string         &material_id_field,
           const std::string         &boundary_id_field,
           const std::string         &manifold_id_field)
  {
    // Get a non-const reference to the triangulation
    auto &tria = const_cast<Triangulation<dim, spacedim> &>(
      dof_handler.get_triangulation());

    // Make sure the triangulation is actually a serial triangulation
    auto parallel_tria =
      dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria);
    AssertThrow(parallel_tria == nullptr,
                ExcMessage(
                  "The input triangulation must be a serial triangulation."));

    // Clear the triangulation to ensure it is empty before reading
    tria.clear();
    // Read the mesh from the VTK file
    read_tria(vtk_filename,
              tria,
              cleanup,
              relative_tolerance,
              material_id_field,
              boundary_id_field,
              manifold_id_field);

    Vector<double> raw_data_vector;
    read_all_data(vtk_filename, raw_data_vector, cleanup, relative_tolerance);

    auto [fe, data_names_from_fe] =
      vtk_to_finite_element<dim, spacedim>(vtk_filename);

    dof_handler.distribute_dofs(*fe);
    output_vector.reinit(dof_handler.n_dofs());
    data_to_dealii_vector(tria, raw_data_vector, dof_handler, output_vector);

    AssertDimension(dof_handler.n_dofs(), output_vector.size());
    AssertDimension(dof_handler.get_fe().n_blocks(), data_names_from_fe.size());
    data_names = data_names_from_fe;
  }

#  include "vtk/utilities.inst"

} // namespace VTKWrappers

#endif

DEAL_II_NAMESPACE_CLOSE
