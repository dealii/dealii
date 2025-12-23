// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------
#include <deal.II/vtk/utilities.h>

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

// Other VTK headers
#  include <vtkCellData.h>
#  if DEAL_II_VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 3, 0)
#    include <vtkCleanUnstructuredGrid.h>
#  endif
#  include <vtkDataArray.h>
#  include <vtkPointData.h>
#  include <vtkSmartPointer.h>
#  include <vtkUnstructuredGrid.h>
#  include <vtkUnstructuredGridReader.h>

#  include <stdexcept>

DEAL_II_NAMESPACE_OPEN

namespace VTKWrappers
{
  template <int dim, int spacedim>
  void
  read_tria(const std::string            &vtk_filename,
            Triangulation<dim, spacedim> &tria,
            const bool                    cleanup,
            const double                  relative_tolerance)
  {
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    {
      // check that the file exists
      std::ifstream file(vtk_filename);
      AssertThrow(file.good(),
                  ExcMessage("VTK file not found: " + vtk_filename));
    }
    reader->SetFileName(vtk_filename.c_str());
    reader->Update();
    vtkUnstructuredGrid *grid = reader->GetOutput();
    AssertThrow(grid, ExcMessage("Failed to read VTK file: " + vtk_filename));

#  if DEAL_II_VTK_VERSION_NUMBER >= VTK_VERSION_CHECK(9, 3, 0)
    auto cleaner = vtkSmartPointer<vtkCleanUnstructuredGrid>::New();

    // Cleanup the triangulation if requested
    if (cleanup)
      {
        cleaner->SetToleranceIsAbsolute(false);
        cleaner->SetTolerance(relative_tolerance);
        cleaner->SetInputData(grid);
        cleaner->Update();
        grid = cleaner->GetOutput();
        AssertThrow(grid,
                    ExcMessage("Failed to clean VTK file: " + vtk_filename));
      }
#  else
    (void)cleanup;            // avoid unused variable warning
    (void)relative_tolerance; // avoid unused variable warning
    deallog << "VTK version < 9.3: skipping cleanup step." << std::endl;
#  endif

    // Read points
    vtkPoints                   *vtk_points = grid->GetPoints();
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
              "VTK file has non-zero coordinate in unused dimension."));
      }

    // Read cells
    std::vector<CellData<dim>> cells;
    const vtkIdType            n_cells = grid->GetNumberOfCells();
    for (vtkIdType i = 0; i < n_cells; ++i)
      {
        vtkCell *cell = grid->GetCell(i);
        if constexpr (dim == 1)
          {
            if (cell->GetCellType() != VTK_LINE)
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
            cells.push_back(cell_data);
          }
        else if constexpr (dim == 2)
          {
            if (cell->GetCellType() == VTK_QUAD)
              {
                AssertThrow(cell->GetNumberOfPoints() == 4,
                            ExcMessage(
                              "Only quad cells with 4 points are supported."));
                CellData<2> cell_data(4);
                for (unsigned int j = 0; j < 4; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTK_TRIANGLE)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 3,
                  ExcMessage(
                    "Only triangle cells with 3 points are supported."));
                CellData<2> cell_data(3);
                for (unsigned int j = 0; j < 3; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
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
            if (cell->GetCellType() == VTK_HEXAHEDRON)
              {
                AssertThrow(cell->GetNumberOfPoints() == 8,
                            ExcMessage(
                              "Only hex cells with 8 points are supported."));
                CellData<3> cell_data(8);
                for (unsigned int j = 0; j < 8; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                // Numbering of vertices in VTK files is different from deal.II
                std::swap(cell_data.vertices[2], cell_data.vertices[3]);
                std::swap(cell_data.vertices[6], cell_data.vertices[7]);
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTK_TETRA)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 4,
                  ExcMessage(
                    "Only tetrahedron cells with 4 points are supported."));
                CellData<3> cell_data(4);
                for (unsigned int j = 0; j < 4; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTK_WEDGE)
              {
                AssertThrow(cell->GetNumberOfPoints() == 6,
                            ExcMessage(
                              "Only prism cells with 6 points are supported."));
                CellData<3> cell_data(6);
                for (unsigned int j = 0; j < 6; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
                cells.push_back(cell_data);
              }
            else if (cell->GetCellType() == VTK_PYRAMID)
              {
                AssertThrow(
                  cell->GetNumberOfPoints() == 5,
                  ExcMessage(
                    "Only pyramid cells with 5 points are supported."));
                CellData<3> cell_data(5);
                for (unsigned int j = 0; j < 5; ++j)
                  cell_data.vertices[j] = cell->GetPointId(j);
                cell_data.material_id = 0;
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
    tria.create_triangulation(points, cells, SubCellData());
  }



  void
  read_cell_data(const std::string &vtk_filename,
                 const std::string &cell_data_name,
                 Vector<double>    &output_vector)
  {
    {
      // check that the file exists
      std::ifstream file(vtk_filename);
      AssertThrow(file.good(),
                  ExcMessage("VTK file not found: " + vtk_filename));
    }
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtk_filename.c_str());
    reader->Update();
    vtkUnstructuredGrid *grid = reader->GetOutput();
    AssertThrow(grid, ExcMessage("Failed to read VTK file: " + vtk_filename));
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
                   const std::string &point_data_name,
                   Vector<double>    &output_vector)
  {
    {
      // check that the file exists
      std::ifstream file(vtk_filename);
      AssertThrow(file.good(),
                  ExcMessage("VTK file not found: " + vtk_filename));
    }
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtk_filename.c_str());
    reader->Update();
    vtkUnstructuredGrid *grid = reader->GetOutput();
    AssertThrow(grid, ExcMessage("Failed to read VTK file: " + vtk_filename));
    vtkDataArray *data_array =
      grid->GetPointData()->GetArray(point_data_name.c_str());
    AssertThrow(data_array,
                ExcMessage("Point data array '" + point_data_name +
                           "' not found in VTK file: " + vtk_filename));
    vtkIdType n_tuples     = data_array->GetNumberOfTuples();
    int       n_components = data_array->GetNumberOfComponents();
    output_vector.reinit(n_tuples * n_components);
    for (vtkIdType i = 0; i < n_tuples; ++i)
      for (int j = 0; j < n_components; ++j)
        output_vector[i * n_components + j] = data_array->GetComponent(i, j);
  }



  template <int dim, int spacedim>
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
            std::vector<std::string>>
  vtk_to_finite_element(const std::string &vtk_filename)
  {
    std::vector<std::string> data_names;
    auto reader = vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(vtk_filename.c_str());
    reader->Update();
    vtkUnstructuredGrid *grid = reader->GetOutput();
    AssertThrow(grid, ExcMessage("Failed to read VTK file: " + vtk_filename));

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
            if (cell_type == VTK_TRIANGLE)
              {
                is_simplex = true;
                break;
              }
          }
        else if constexpr (dim == 3)
          {
            if (cell_type == VTK_TETRA)
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
#else
DEAL_II_NAMESPACE_OPEN

namespace VTKWrappers
{
  template <int dim, int spacedim>
  void
  read_tria(const std::string &,
            Triangulation<dim, spacedim> &,
            const bool,
            const double)
  {
    AssertThrow(false, ExcMessage("deal.II is not built with VTK support."));
  }

  void
  read_cell_data(const std::string &, const std::string &, Vector<double> &)
  {
    AssertThrow(false, ExcMessage("deal.II is not built with VTK support."));
  }

  void
  read_vertex_data(const std::string &, const std::string &, Vector<double> &)
  {
    AssertThrow(false, ExcMessage("deal.II is not built with VTK support."));
  }

  template <int dim, int spacedim>
  std::pair<std::unique_ptr<FiniteElement<dim, spacedim>>,
            std::vector<std::string>>
  vtk_to_finite_element(const std::string &)
  {
    AssertThrow(false, ExcMessage("deal.II is not built with VTK support."));
  }
#endif

#include "vtk/utilities.inst"

} // namespace VTKWrappers

DEAL_II_NAMESPACE_CLOSE
