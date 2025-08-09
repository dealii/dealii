// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/parameter_handler.h>

#include <deal.II/gmsh/utilities.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/opencascade/utilities.h>

#include <cstdio>

#ifdef DEAL_II_WITH_GMSH


#  ifdef DEAL_II_GMSH_WITH_API
#    include <deal.II/base/config.h>

#    include <deal.II/base/utilities.h>

#    include <deal.II/grid/cell_id.h>
#    include <deal.II/grid/tria_description.h>
#  endif


#  ifdef DEAL_II_WITH_OPENCASCADE
#    include <TopoDS_Edge.hxx>
#  endif

DEAL_II_NAMESPACE_OPEN

namespace Gmsh
{
  AdditionalParameters::AdditionalParameters(
    const double       characteristic_length,
    const std::string &output_base_name)
    : characteristic_length(characteristic_length)
    , output_base_name(output_base_name)
  {}



  void
  AdditionalParameters::add_parameters(ParameterHandler &prm)
  {
    prm.add_parameter("Characteristic length", characteristic_length);
    prm.add_parameter("Intermediate file name base",
                      output_base_name,
                      "Keep empty, if you want the program to generate "
                      "temporary files, and then remove them when they "
                      "are no longer used.");
  }



#  ifdef DEAL_II_WITH_OPENCASCADE
  template <int spacedim>
  void
  create_triangulation_from_boundary_curve(const TopoDS_Edge          &boundary,
                                           Triangulation<2, spacedim> &tria,
                                           const AdditionalParameters &prm)
  {
    std::string base_name = prm.output_base_name;

    // If necessary, create a temp directory to put files into. The
    // following variable will hold the name of the tmp dir; it is
    // initialized with a template, and 'mkdtemp' then overwrites it
    // with the name of the directory it creates.
    char tmp_dir_name[] = "ctfbc-XXXXXX";
    if (prm.output_base_name.empty())
      {
        const char *temp = mkdtemp(tmp_dir_name);
        AssertThrow(temp != nullptr,
                    ExcMessage("Creating temporary directory failed!"));
        base_name = temp;
        base_name += "tmp";
      }

    const std::string iges_file_name     = base_name + ".iges";
    const std::string geo_file_name      = base_name + ".geo";
    const std::string msh_file_name      = base_name + ".msh";
    const std::string log_file_name      = base_name + ".log";
    const std::string warnings_file_name = base_name + "_warn.log";

    dealii::OpenCASCADE::write_IGES(boundary, iges_file_name);

    std::ofstream geofile;
    geofile.open(geo_file_name);
    geofile << "Merge \"" << iges_file_name << "\";" << std::endl
            << "Line Loop (2) = {1};" << std::endl
            << "Plane Surface (3) = {2};" << std::endl
            << "Characteristic Length { 1 } = " << prm.characteristic_length
            << ";" << std::endl
            << "Mesh.RecombineAll = 1;" << std::endl
            << "Mesh.SubdivisionAlgorithm = 1;" << std::endl;
    geofile.close();

    std::stringstream command;
    command << DEAL_II_GMSH_EXECUTABLE_PATH << " -2 " << geo_file_name << " 1> "
            << log_file_name << " 2> " << warnings_file_name;

    const auto ret_value = std::system(command.str().c_str());
    AssertThrow(ret_value == 0,
                ExcMessage("Gmsh failed to run. Check the " + log_file_name +
                           " file."));

    std::ifstream grid_file(msh_file_name);
    Assert(grid_file.fail() == false, ExcIO());

    GridIn<2, spacedim> gridin;
    gridin.attach_triangulation(tria);
    gridin.read_msh(grid_file);

    // Clean up files if a tmp directory was used:
    if (prm.output_base_name.empty())
      {
        // declaring the list without a type, i.e.,
        //
        //     auto filenames = {{iges_file_name, geo_file_name, ...}})
        //
        // causes internal compiler errors with GCC's concepts implementation,
        // so give it an explicit type:
        const std::array<const std::string *, 5> filenames{
          {&iges_file_name,
           &geo_file_name,
           &msh_file_name,
           &log_file_name,
           &warnings_file_name}};
        for (const std::string *filename : filenames)
          {
            const auto ret_value = std::remove(filename->c_str());
            AssertThrow(ret_value == 0,
                        ExcMessage("Failed to remove " + *filename));
          }

        // Finally also remove the tmp directory:
        const auto ret_value = std::remove(tmp_dir_name);
        AssertThrow(ret_value == 0,
                    ExcMessage("Failed to remove " +
                               std::string(tmp_dir_name)));
      }
  }
#  endif

#  ifdef DEAL_II_GMSH_WITH_API
  template <int dim, int spacedim>
  void
  read_partitioned_msh(
    parallel::fullydistributed::Triangulation<dim, spacedim> &tria,
    const MPI_Comm                                            mpi_comm,
    const std::string                                        &file_prefix,
    const std::string                                        &file_suffix)
  {
    const unsigned int nprocs = Utilities::MPI::n_mpi_processes(mpi_comm);
    const unsigned int rank   = Utilities::MPI::this_mpi_process(mpi_comm);

    std::string fname =
      file_prefix + "_" + std::to_string(rank + 1) + "." + file_suffix;

    if (nprocs == 1)
      {
        fname = file_prefix + "." + file_suffix;
        std::ifstream f(fname);
        AssertThrow(f.good(), ExcMessage("Missing mesh file: " + fname));
      }
    else
      {
        for (unsigned int i = 1; i <= nprocs; ++i)
          {
            const std::string check_fname =
              file_prefix + "_" + std::to_string(i) + "." + file_suffix;
            std::ifstream f(check_fname);
            AssertThrow(f.good(),
                        ExcMessage("Missing mesh file: " + check_fname));
          }

        const std::string extra_fname =
          file_prefix + "_" + std::to_string(nprocs + 1) + "." + file_suffix;
        std::ifstream f(extra_fname);
        AssertThrow(!f.good(),
                    ExcMessage("Unexpected extra mesh file(s): " +
                               extra_fname));
      }

    const std::map<int, std::uint8_t> gmsh_to_dealii_type = {
      {15, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {7, 5}, {6, 6}, {5, 7}};

    const std::array<std::vector<unsigned int>, 8> gmsh_to_dealii = {
      {{0},
       {0, 1},
       {0, 1, 2},
       {0, 1, 3, 2},
       {0, 1, 2, 3},
       {0, 1, 3, 2, 4},
       {0, 1, 2, 3, 4, 5},
       {0, 1, 3, 2, 4, 5, 7, 6}}};

    gmsh::initialize();
    gmsh::option::setNumber("General.Verbosity", 0);
    gmsh::open(fname);

    std::map<unsigned long, unsigned int> ghost_map;

    std::vector<std::pair<int, int>> entities;
    gmsh::model::getEntities(entities);

    for (const auto &e : entities)
      {
        const int entity_dim = e.first;
        const int entity_tag = e.second;

        if (entity_dim == dim)
          {
            std::vector<std::size_t> element_tags;
            std::vector<int>         partitions;

            gmsh::model::mesh::getGhostElements(entity_dim,
                                                entity_tag,
                                                element_tags,
                                                partitions);

            for (std::size_t i = 0; i < element_tags.size(); ++i)
              ghost_map[element_tags[i]] =
                static_cast<unsigned int>(partitions[i] - 1);
          }
      }

    std::vector<std::size_t> node_tags;
    std::vector<double>      coords, parametric_coords;
    gmsh::model::mesh::getNodes(node_tags, coords, parametric_coords);

    TriangulationDescription::Description<dim, spacedim>
      triangulation_description;
    triangulation_description.comm = mpi_comm;

    auto &vertices = triangulation_description.coarse_cell_vertices;
    auto &cells    = triangulation_description.coarse_cells;
    auto &coarse_cell_ids =
      triangulation_description.coarse_cell_index_to_coarse_cell_id;
    auto &cell_infos = triangulation_description.cell_infos;

    cell_infos.resize(1);
    vertices.resize(node_tags.size(), Point<spacedim>());

    std::map<std::size_t, unsigned int> node_tag_to_index;
    for (unsigned int i = 0; i < node_tags.size(); ++i)
      {
        node_tag_to_index[node_tags[i]] = i;
        for (unsigned int d = 0; d < spacedim; ++d)
          vertices[i][d] = coords[3 * i + d];
      }

    for (const auto &e : entities)
      {
        const int entity_dim = e.first;
        const int entity_tag = e.second;

        if (entity_dim == dim)
          {
            std::vector<int>                      element_types;
            std::vector<std::vector<std::size_t>> element_ids, element_nodes;

            gmsh::model::mesh::getElements(element_types,
                                           element_ids,
                                           element_nodes,
                                           entity_dim,
                                           entity_tag);



            for (unsigned int i = 0; i < element_types.size(); ++i)
              {
                if (element_ids[i].empty())
                  continue;


                const unsigned int n_vertices =
                  element_nodes[i].size() / element_ids[i].size();

                for (unsigned int j = 0; j < element_ids[i].size(); ++j)
                  {
                    CellData<dim> cell(n_vertices);
                    cell.material_id = 0;

                    const auto &type = gmsh_to_dealii_type.at(element_types[i]);

                    for (unsigned int v = 0; v < n_vertices; ++v)
                      {
                        const std::size_t node_tag =
                          element_nodes[i][j * n_vertices +
                                           gmsh_to_dealii[type][v]];
                        AssertThrow(node_tag_to_index.find(node_tag) !=
                                      node_tag_to_index.end(),
                                    ExcMessage("Node tag " +
                                               std::to_string(node_tag) +
                                               " not found in node list!"));
                        cell.vertices[v] = node_tag_to_index[node_tag];
                      }

                    cells.push_back(cell);
                    coarse_cell_ids.push_back(element_ids[i][j]);

                    TriangulationDescription::CellData<dim> cell_info;
                    cell_info.id =
                      CellId(element_ids[i][j], {}).template to_binary<dim>();

                    auto it = ghost_map.find(element_ids[i][j]);
                    if (it != ghost_map.end())
                      cell_info.subdomain_id = it->second;
                    else
                      cell_info.subdomain_id = rank;

                    cell_info.level_subdomain_id = cell_info.subdomain_id;

                    cell_infos[0].push_back(cell_info);
                  }
              }
          }
      }

    triangulation_description.settings =
      TriangulationDescription::Settings::default_setting;

    tria.create_triangulation(triangulation_description);

    MPI_Barrier(mpi_comm);

    gmsh::clear();
    gmsh::finalize();
  }
#  endif // DEAL_II_GMSH_WITH_API

  // explicit instantiations
// We don't build the utilities.inst file if deal.II isn't configured
// with GMSH, but doxygen doesn't know that and tries to find that
// file anyway for parsing -- which then of course it fails on. So
// exclude the following from doxygen consideration.
#  ifndef DOXYGEN
#    include "gmsh/utilities.inst"
#  endif
} // namespace Gmsh

DEAL_II_NAMESPACE_CLOSE

#endif
