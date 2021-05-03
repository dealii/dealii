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


#include <deal.II/base/exceptions.h>
#include <deal.II/base/path_search.h>
#include <deal.II/base/patterns.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/io/ios_state.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/serialization/serialization.hpp>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <algorithm>
#include <cctype>
#include <fstream>
#include <functional>
#include <map>

#ifdef DEAL_II_WITH_ASSIMP
#  include <assimp/Importer.hpp>  // C++ importer interface
#  include <assimp/postprocess.h> // Post processing flags
#  include <assimp/scene.h>       // Output data structure
#endif

#ifdef DEAL_II_TRILINOS_WITH_SEACAS
#  include <exodusII.h>
#endif


DEAL_II_NAMESPACE_OPEN


namespace
{
  /**
   * In 1d, boundary indicators are associated with vertices, but this is not
   * currently passed through the SubcellData structure. This function sets
   * boundary indicators on vertices after the triangulation has already been
   * created.
   *
   * TODO: Fix this properly via SubcellData
   */
  template <int spacedim>
  void
  assign_1d_boundary_ids(
    const std::map<unsigned int, types::boundary_id> &boundary_ids,
    Triangulation<1, spacedim> &                      triangulation)
  {
    if (boundary_ids.size() > 0)
      for (const auto &cell : triangulation.active_cell_iterators())
        for (unsigned int f : GeometryInfo<1>::face_indices())
          if (boundary_ids.find(cell->vertex_index(f)) != boundary_ids.end())
            {
              AssertThrow(
                cell->at_boundary(f),
                ExcMessage(
                  "You are trying to prescribe boundary ids on the face "
                  "of a 1d cell (i.e., on a vertex), but this face is not actually at "
                  "the boundary of the mesh. This is not allowed."));
              cell->face(f)->set_boundary_id(
                boundary_ids.find(cell->vertex_index(f))->second);
            }
  }


  template <int dim, int spacedim>
  void
  assign_1d_boundary_ids(const std::map<unsigned int, types::boundary_id> &,
                         Triangulation<dim, spacedim> &)
  {
    // we shouldn't get here since boundary ids are not assigned to
    // vertices except in 1d
    Assert(dim != 1, ExcInternalError());
  }
} // namespace

template <int dim, int spacedim>
GridIn<dim, spacedim>::GridIn()
  : tria(nullptr, typeid(*this).name())
  , default_format(ucd)
{}



template <int dim, int spacedim>
GridIn<dim, spacedim>::GridIn(Triangulation<dim, spacedim> &t)
  : tria(&t, typeid(*this).name())
  , default_format(ucd)
{}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::attach_triangulation(Triangulation<dim, spacedim> &t)
{
  tria = &t;
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_vtk(std::istream &in)
{
  std::string line;

  // verify that the first, third and fourth lines match
  // expectations. the second line of the file may essentially be
  // anything the author of the file chose to identify what's in
  // there, so we just ensure that we can read it
  {
    std::string text[4];
    text[0] = "# vtk DataFile Version 3.0";
    text[1] = "****";
    text[2] = "ASCII";
    text[3] = "DATASET UNSTRUCTURED_GRID";

    for (unsigned int i = 0; i < 4; ++i)
      {
        getline(in, line);
        if (i != 1)
          AssertThrow(
            line.compare(text[i]) == 0,
            ExcMessage(
              std::string(
                "While reading VTK file, failed to find a header line with text <") +
              text[i] + ">"));
      }
  }

  ///////////////////Declaring storage and mappings//////////////////

  std::vector<Point<spacedim>> vertices;
  std::vector<CellData<dim>>   cells;
  SubCellData                  subcelldata;

  std::string keyword;

  in >> keyword;

  //////////////////Processing the POINTS section///////////////

  if (keyword == "POINTS")
    {
      unsigned int n_vertices;
      in >> n_vertices;

      in >> keyword; // float, double, int, char, etc.

      for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
        {
          // VTK format always specifies vertex coordinates with 3 components
          Point<3> x;
          in >> x(0) >> x(1) >> x(2);

          vertices.emplace_back();
          for (unsigned int d = 0; d < spacedim; ++d)
            vertices.back()(d) = x(d);
        }
    }

  else
    AssertThrow(false,
                ExcMessage(
                  "While reading VTK file, failed to find POINTS section"));

  in >> keyword;

  unsigned int n_geometric_objects = 0;
  unsigned int n_ints;

  bool is_quad_or_hex_mesh = false;
  bool is_tria_or_tet_mesh = false;

  if (keyword == "CELLS")
    {
      // jump to the `CELL_TYPES` section and read in cell types
      std::vector<unsigned int> cell_types;
      {
        std::streampos oldpos = in.tellg();


        while (in >> keyword)
          if (keyword == "CELL_TYPES")
            {
              in >> n_ints;

              cell_types.resize(n_ints);

              for (unsigned int i = 0; i < n_ints; ++i)
                in >> cell_types[i];

              break;
            }

        in.seekg(oldpos);
      }

      in >> n_geometric_objects;
      in >> n_ints; // Ignore this, since we don't need it.

      if (dim == 3)
        {
          for (unsigned int count = 0; count < n_geometric_objects; count++)
            {
              unsigned int n_vertices;
              in >> n_vertices;

              // VTK_TETRA is 10, VTK_HEXAHEDRON is 12
              if (cell_types[count] == 10 || cell_types[count] == 12)
                {
                  if (cell_types[count] == 10)
                    is_tria_or_tet_mesh = true;
                  if (cell_types[count] == 12)
                    is_quad_or_hex_mesh = true;

                  // we assume that the file contains first all cells,
                  // and only then any faces or lines
                  AssertThrow(subcelldata.boundary_quads.size() == 0 &&
                                subcelldata.boundary_lines.size() == 0,
                              ExcNotImplemented());

                  cells.emplace_back(n_vertices);

                  for (unsigned int j = 0; j < n_vertices;
                       j++) // loop to feed data
                    in >> cells.back().vertices[j];

                  // Hexahedra need a permutation to go from VTK numbering
                  // to deal numbering
                  if (cell_types[count] == 12)
                    {
                      std::swap(cells.back().vertices[2],
                                cells.back().vertices[3]);
                      std::swap(cells.back().vertices[6],
                                cells.back().vertices[7]);
                    }

                  cells.back().material_id = 0;
                }
              // VTK_TRIANGLE is 5, VTK_QUAD is 9
              else if (cell_types[count] == 5 || cell_types[count] == 9)
                {
                  if (cell_types[count] == 5)
                    is_tria_or_tet_mesh = true;
                  if (cell_types[count] == 9)
                    is_quad_or_hex_mesh = true;

                  // we assume that the file contains first all cells,
                  // then all faces, and finally all lines
                  AssertThrow(subcelldata.boundary_lines.size() == 0,
                              ExcNotImplemented());

                  subcelldata.boundary_quads.emplace_back(n_vertices);

                  for (unsigned int j = 0; j < n_vertices;
                       j++) // loop to feed the data to the boundary
                    in >> subcelldata.boundary_quads.back().vertices[j];

                  subcelldata.boundary_quads.back().material_id = 0;
                }
              // VTK_LINE is 3
              else if (cell_types[count] == 3)
                {
                  subcelldata.boundary_lines.emplace_back(n_vertices);

                  for (unsigned int j = 0; j < n_vertices;
                       j++) // loop to feed the data to the boundary
                    in >> subcelldata.boundary_lines.back().vertices[j];

                  subcelldata.boundary_lines.back().material_id = 0;
                }

              else
                AssertThrow(
                  false,
                  ExcMessage(
                    "While reading VTK file, unknown cell type encountered"));
            }
        }
      else if (dim == 2)
        {
          for (unsigned int count = 0; count < n_geometric_objects; count++)
            {
              unsigned int n_vertices;
              in >> n_vertices;

              // VTK_TRIANGLE is 5, VTK_QUAD is 9
              if (cell_types[count] == 5 || cell_types[count] == 9)
                {
                  // we assume that the file contains first all cells,
                  // and only then any faces
                  AssertThrow(subcelldata.boundary_lines.size() == 0,
                              ExcNotImplemented());

                  if (cell_types[count] == 5)
                    is_tria_or_tet_mesh = true;
                  if (cell_types[count] == 9)
                    is_quad_or_hex_mesh = true;

                  cells.emplace_back(n_vertices);

                  for (unsigned int j = 0; j < n_vertices;
                       j++) // loop to feed data
                    in >> cells.back().vertices[j];

                  // Quadrilaterals need a permutation to go from VTK numbering
                  // to deal numbering
                  if (cell_types[count] == 9)
                    {
                      // Like Hexahedra - the last two vertices need to be
                      // flipped
                      std::swap(cells.back().vertices[2],
                                cells.back().vertices[3]);
                    }

                  cells.back().material_id = 0;
                }
              // VTK_LINE is 3
              else if (cell_types[count] == 3)
                {
                  // If this is encountered, the pointer comes out of the loop
                  // and starts processing boundaries.
                  subcelldata.boundary_lines.emplace_back(n_vertices);

                  for (unsigned int j = 0; j < n_vertices;
                       j++) // loop to feed the data to the boundary
                    {
                      in >> subcelldata.boundary_lines.back().vertices[j];
                    }

                  subcelldata.boundary_lines.back().material_id = 0;
                }

              else
                AssertThrow(
                  false,
                  ExcMessage(
                    "While reading VTK file, unknown cell type encountered"));
            }
        }
      else if (dim == 1)
        {
          for (unsigned int count = 0; count < n_geometric_objects; count++)
            {
              unsigned int type;
              in >> type;

              AssertThrow(
                cell_types[count] == 3 && type == 2,
                ExcMessage(
                  "While reading VTK file, unknown cell type encountered"));
              cells.emplace_back(type);

              for (unsigned int j = 0; j < type; j++) // loop to feed data
                in >> cells.back().vertices[j];

              cells.back().material_id = 0;
            }
        }
      else
        AssertThrow(false,
                    ExcMessage(
                      "While reading VTK file, failed to find CELLS section"));

      /////////////////////Processing the CELL_TYPES
      /// section////////////////////////

      in >> keyword;

      AssertThrow(
        keyword == "CELL_TYPES",
        ExcMessage(std::string(
          "While reading VTK file, missing CELL_TYPES section. Found <" +
          keyword + "> instead.")));

      in >> n_ints;
      AssertThrow(
        n_ints == n_geometric_objects,
        ExcMessage("The VTK reader found a CELL_DATA statement "
                   "that lists a total of " +
                   Utilities::int_to_string(n_ints) +
                   " cell data objects, but this needs to "
                   "equal the number of cells (which is " +
                   Utilities::int_to_string(cells.size()) +
                   ") plus the number of quads (" +
                   Utilities::int_to_string(subcelldata.boundary_quads.size()) +
                   " in 3d or the number of lines (" +
                   Utilities::int_to_string(subcelldata.boundary_lines.size()) +
                   ") in 2d."));

      int tmp_int;
      for (unsigned int i = 0; i < n_ints; ++i)
        in >> tmp_int;

      // Ignore everything up to CELL_DATA
      while (in >> keyword)
        if (keyword == "CELL_DATA")
          {
            unsigned int n_ids;
            in >> n_ids;

            AssertThrow(n_ids == n_geometric_objects,
                        ExcMessage("The VTK reader found a CELL_DATA statement "
                                   "that lists a total of " +
                                   Utilities::int_to_string(n_ids) +
                                   " cell data objects, but this needs to "
                                   "equal the number of cells (which is " +
                                   Utilities::int_to_string(cells.size()) +
                                   ") plus the number of quads (" +
                                   Utilities::int_to_string(
                                     subcelldata.boundary_quads.size()) +
                                   " in 3d or the number of lines (" +
                                   Utilities::int_to_string(
                                     subcelldata.boundary_lines.size()) +
                                   ") in 2d."));

            const std::vector<std::string> data_sets{"MaterialID",
                                                     "ManifoldID"};

            for (unsigned int i = 0; i < data_sets.size(); ++i)
              {
                // Ignore everything until we get to a SCALARS data set
                while (in >> keyword)
                  if (keyword == "SCALARS")
                    {
                      // Now see if we know about this type of data set,
                      // if not, just ignore everything till the next SCALARS
                      // keyword
                      std::string field_name;
                      in >> field_name;
                      if (std::find(data_sets.begin(),
                                    data_sets.end(),
                                    field_name) == data_sets.end())
                        // The data set here is not one of the ones we know, so
                        // keep ignoring everything until the next SCALARS
                        // keyword.
                        continue;

                      // Now we got somewhere. Proceed from here, assert
                      // that the type of the table is int, and ignore the
                      // rest of the line.
                      // SCALARS MaterialID int 1
                      // (the last number is optional)
                      std::string line;
                      std::getline(in, line);
                      AssertThrow(
                        line.substr(1,
                                    std::min(static_cast<std::size_t>(3),
                                             line.size() - 1)) == "int",
                        ExcMessage(
                          "While reading VTK file, material- and manifold IDs can only have type 'int'."));

                      in >> keyword;
                      AssertThrow(
                        keyword == "LOOKUP_TABLE",
                        ExcMessage(
                          "While reading VTK file, missing keyword 'LOOKUP_TABLE'."));

                      in >> keyword;
                      AssertThrow(
                        keyword == "default",
                        ExcMessage(
                          "While reading VTK file, missing keyword 'default'."));

                      // read material or manifold ids first for all cells,
                      // then for all faces, and finally for all lines. the
                      // assumption that cells come before all faces and
                      // lines has been verified above via an assertion, so
                      // the order used in the following blocks makes sense
                      for (unsigned int i = 0; i < cells.size(); i++)
                        {
                          int id;
                          in >> id;
                          if (field_name == "MaterialID")
                            cells[i].material_id =
                              static_cast<types::material_id>(id);
                          else if (field_name == "ManifoldID")
                            cells[i].manifold_id =
                              static_cast<types::manifold_id>(id);
                          else
                            Assert(false, ExcInternalError());
                        }

                      if (dim == 3)
                        {
                          for (auto &boundary_quad : subcelldata.boundary_quads)
                            {
                              int id;
                              in >> id;
                              if (field_name == "MaterialID")
                                boundary_quad.material_id =
                                  static_cast<types::material_id>(id);
                              else if (field_name == "ManifoldID")
                                boundary_quad.manifold_id =
                                  static_cast<types::manifold_id>(id);
                              else
                                Assert(false, ExcInternalError());
                            }
                          for (auto &boundary_line : subcelldata.boundary_lines)
                            {
                              int id;
                              in >> id;
                              if (field_name == "MaterialID")
                                boundary_line.material_id =
                                  static_cast<types::material_id>(id);
                              else if (field_name == "ManifoldID")
                                boundary_line.manifold_id =
                                  static_cast<types::manifold_id>(id);
                              else
                                Assert(false, ExcInternalError());
                            }
                        }
                      else if (dim == 2)
                        {
                          for (auto &boundary_line : subcelldata.boundary_lines)
                            {
                              int id;
                              in >> id;
                              if (field_name == "MaterialID")
                                boundary_line.material_id =
                                  static_cast<types::material_id>(id);
                              else if (field_name == "ManifoldID")
                                boundary_line.manifold_id =
                                  static_cast<types::manifold_id>(id);
                              else
                                Assert(false, ExcInternalError());
                            }
                        }
                    }
              }
          }

      Assert(subcelldata.check_consistency(dim), ExcInternalError());


      // TODO: the functions below (GridTools::delete_unused_vertices(),
      // GridTools::invert_all_negative_measure_cells(),
      // GridTools::consistently_order_cells()) need to be
      // revisited for simplex/mixed meshes

      if (dim == 1 || (is_quad_or_hex_mesh && !is_tria_or_tet_mesh))
        {
          GridTools::delete_unused_vertices(vertices, cells, subcelldata);

          if (dim == spacedim)
            GridTools::invert_all_negative_measure_cells(vertices, cells);

          GridTools::consistently_order_cells(cells);
          tria->create_triangulation(vertices, cells, subcelldata);
        }
      else
        {
          // simplex or mixed mesh
          tria->create_triangulation(vertices, cells, subcelldata);
        }
    }
  else
    AssertThrow(false,
                ExcMessage(
                  "While reading VTK file, failed to find CELLS section"));
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_vtu(std::istream &in)
{
  namespace pt = boost::property_tree;
  pt::ptree tree;
  pt::read_xml(in, tree);
  auto section = tree.get_optional<std::string>("VTKFile.dealiiData");

  AssertThrow(section,
              ExcMessage(
                "While reading a VTU file, failed to find dealiiData section. "
                "Notice that we can only read grid files in .vtu format that "
                "were created by the deal.II library, using a call to "
                "GridOut::write_vtu(), where the flag "
                "GridOutFlags::Vtu::serialize_triangulation is set to true."));

  const auto decoded =
    Utilities::decode_base64({section->begin(), section->end() - 1});
  const auto string_archive =
    Utilities::decompress({decoded.begin(), decoded.end()});
  std::istringstream              in_stream(string_archive);
  boost::archive::binary_iarchive ia(in_stream);
  tria->load(ia, 0);
}


template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_unv(std::istream &in)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  Assert((dim == 2) || (dim == 3), ExcNotImplemented());

  AssertThrow(in, ExcIO());
  skip_comment_lines(in, '#'); // skip comments (if any) at beginning of file

  int tmp;

  AssertThrow(in, ExcIO());
  in >> tmp;
  AssertThrow(in, ExcIO());
  in >> tmp;

  // section 2411 describes vertices: see
  // http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2411
  AssertThrow(tmp == 2411, ExcUnknownSectionType(tmp));

  std::vector<Point<spacedim>> vertices; // vector of vertex coordinates
  std::map<int, int>
    vertex_indices; // # vert in unv (key) ---> # vert in deal.II (value)

  int no_vertex = 0; // deal.II

  while (tmp != -1) // we do until reach end of 2411
    {
      int    no; // unv
      int    dummy;
      double x[3];

      AssertThrow(in, ExcIO());
      in >> no;

      tmp = no;
      if (tmp == -1)
        break;

      in >> dummy >> dummy >> dummy;

      AssertThrow(in, ExcIO());
      in >> x[0] >> x[1] >> x[2];

      vertices.emplace_back();

      for (unsigned int d = 0; d < spacedim; d++)
        vertices.back()(d) = x[d];

      vertex_indices[no] = no_vertex;

      no_vertex++;
    }

  AssertThrow(in, ExcIO());
  in >> tmp;
  AssertThrow(in, ExcIO());
  in >> tmp;

  // section 2412 describes elements: see
  // http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2412
  AssertThrow(tmp == 2412, ExcUnknownSectionType(tmp));

  std::vector<CellData<dim>> cells; // vector of cells
  SubCellData                subcelldata;

  std::map<int, int>
    cell_indices; // # cell in unv (key) ---> # cell in deal.II (value)
  std::map<int, int>
    line_indices; // # line in unv (key) ---> # line in deal.II (value)
  std::map<int, int>
    quad_indices; // # quad in unv (key) ---> # quad in deal.II (value)

  int no_cell = 0; // deal.II
  int no_line = 0; // deal.II
  int no_quad = 0; // deal.II

  while (tmp != -1) // we do until reach end of 2412
    {
      int no; // unv
      int type;
      int dummy;

      AssertThrow(in, ExcIO());
      in >> no;

      tmp = no;
      if (tmp == -1)
        break;

      in >> type >> dummy >> dummy >> dummy >> dummy;

      AssertThrow((type == 11) || (type == 44) || (type == 94) || (type == 115),
                  ExcUnknownElementType(type));

      if ((((type == 44) || (type == 94)) && (dim == 2)) ||
          ((type == 115) && (dim == 3))) // cell
        {
          const auto reference_cell = ReferenceCells::get_hypercube<dim>();
          cells.emplace_back();

          AssertThrow(in, ExcIO());
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            in >> cells.back()
                    .vertices[reference_cell.unv_vertex_to_deal_vertex(v)];

          cells.back().material_id = 0;

          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            cells.back().vertices[v] = vertex_indices[cells.back().vertices[v]];

          cell_indices[no] = no_cell;

          no_cell++;
        }
      else if (((type == 11) && (dim == 2)) ||
               ((type == 11) && (dim == 3))) // boundary line
        {
          AssertThrow(in, ExcIO());
          in >> dummy >> dummy >> dummy;

          subcelldata.boundary_lines.emplace_back();

          AssertThrow(in, ExcIO());
          for (unsigned int &vertex :
               subcelldata.boundary_lines.back().vertices)
            in >> vertex;

          subcelldata.boundary_lines.back().material_id = 0;

          for (unsigned int &vertex :
               subcelldata.boundary_lines.back().vertices)
            vertex = vertex_indices[vertex];

          line_indices[no] = no_line;

          no_line++;
        }
      else if (((type == 44) || (type == 94)) && (dim == 3)) // boundary quad
        {
          const auto reference_cell = ReferenceCells::Quadrilateral;
          subcelldata.boundary_quads.emplace_back();

          AssertThrow(in, ExcIO());
          Assert(subcelldata.boundary_quads.back().vertices.size() ==
                   GeometryInfo<2>::vertices_per_cell,
                 ExcInternalError());
          for (const unsigned int v : GeometryInfo<2>::vertex_indices())
            in >> subcelldata.boundary_quads.back()
                    .vertices[reference_cell.unv_vertex_to_deal_vertex(v)];

          subcelldata.boundary_quads.back().material_id = 0;

          for (unsigned int &vertex :
               subcelldata.boundary_quads.back().vertices)
            vertex = vertex_indices[vertex];

          quad_indices[no] = no_quad;

          no_quad++;
        }
      else
        AssertThrow(false,
                    ExcMessage("Unknown element label <" +
                               Utilities::int_to_string(type) +
                               "> when running in dim=" +
                               Utilities::int_to_string(dim)));
    }

  // note that so far all materials and bcs are explicitly set to 0
  // if we do not need more info on materials and bcs - this is end of file
  // if we do - section 2467 or 2477 comes

  in >> tmp; // tmp can be either -1 or end-of-file

  if (!in.eof())
    {
      AssertThrow(in, ExcIO());
      in >> tmp;

      // section 2467 (2477) describes (materials - first and bcs - second) or
      // (bcs - first and materials - second) - sequence depends on which
      // group is created first: see
      // http://www.sdrl.uc.edu/sdrl/referenceinfo/universalfileformats/file-format-storehouse/universal-dataset-number-2467
      AssertThrow((tmp == 2467) || (tmp == 2477), ExcUnknownSectionType(tmp));

      while (tmp != -1) // we do until reach end of 2467 or 2477
        {
          int n_entities; // number of entities in group
          int id;         // id is either material or bc
          int no;         // unv
          int dummy;

          AssertThrow(in, ExcIO());
          in >> dummy;

          tmp = dummy;
          if (tmp == -1)
            break;

          in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >>
            n_entities;

          AssertThrow(in, ExcIO());
          in >> id;

          const unsigned int n_lines =
            (n_entities % 2 == 0) ? (n_entities / 2) : ((n_entities + 1) / 2);

          for (unsigned int line = 0; line < n_lines; line++)
            {
              unsigned int n_fragments;

              if (line == n_lines - 1)
                n_fragments = (n_entities % 2 == 0) ? (2) : (1);
              else
                n_fragments = 2;

              for (unsigned int no_fragment = 0; no_fragment < n_fragments;
                   no_fragment++)
                {
                  AssertThrow(in, ExcIO());
                  in >> dummy >> no >> dummy >> dummy;

                  if (cell_indices.count(no) > 0) // cell - material
                    cells[cell_indices[no]].material_id = id;

                  if (line_indices.count(no) > 0) // boundary line - bc
                    subcelldata.boundary_lines[line_indices[no]].material_id =
                      id;

                  if (quad_indices.count(no) > 0) // boundary quad - bc
                    subcelldata.boundary_quads[quad_indices[no]].material_id =
                      id;
                }
            }
        }
    }

  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  GridTools::delete_unused_vertices(vertices, cells, subcelldata);

  if (dim == spacedim)
    GridTools::invert_all_negative_measure_cells(vertices, cells);

  GridTools::consistently_order_cells(cells);

  tria->create_triangulation(vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_ucd(std::istream &in,
                                const bool    apply_all_indicators_to_manifolds)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  AssertThrow(in, ExcIO());

  // skip comments at start of file
  skip_comment_lines(in, '#');


  unsigned int n_vertices;
  unsigned int n_cells;
  int          dummy;

  in >> n_vertices >> n_cells >> dummy // number of data vectors
    >> dummy                           // cell data
    >> dummy;                          // model data
  AssertThrow(in, ExcIO());

  // set up array of vertices
  std::vector<Point<spacedim>> vertices(n_vertices);
  // set up mapping between numbering
  // in ucd-file (key) and in the
  // vertices vector
  std::map<int, int> vertex_indices;

  for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
    {
      int    vertex_number;
      double x[3];

      // read vertex
      AssertThrow(in, ExcIO());
      in >> vertex_number >> x[0] >> x[1] >> x[2];

      // store vertex
      for (unsigned int d = 0; d < spacedim; ++d)
        vertices[vertex](d) = x[d];
      // store mapping; note that
      // vertices_indices[i] is automatically
      // created upon first usage
      vertex_indices[vertex_number] = vertex;
    }

  // set up array of cells
  std::vector<CellData<dim>> cells;
  SubCellData                subcelldata;

  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      // note that since in the input
      // file we found the number of
      // cells at the top, there
      // should still be input here,
      // so check this:
      AssertThrow(in, ExcIO());

      std::string cell_type;

      // we use an unsigned int because we
      // fill this variable through an read-in process
      unsigned int material_id;

      in >> dummy // cell number
        >> material_id;
      in >> cell_type;

      if (((cell_type == "line") && (dim == 1)) ||
          ((cell_type == "quad") && (dim == 2)) ||
          ((cell_type == "hex") && (dim == 3)))
        // found a cell
        {
          // allocate and read indices
          cells.emplace_back();
          for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
            in >> cells.back().vertices[GeometryInfo<dim>::ucd_to_deal[i]];

          // to make sure that the cast won't fail
          Assert(material_id <= std::numeric_limits<types::material_id>::max(),
                 ExcIndexRange(material_id,
                               0,
                               std::numeric_limits<types::material_id>::max()));
          // we use only material_ids in the range from 0 to
          // numbers::invalid_material_id-1
          AssertIndexRange(material_id, numbers::invalid_material_id);

          if (apply_all_indicators_to_manifolds)
            cells.back().manifold_id =
              static_cast<types::manifold_id>(material_id);
          cells.back().material_id = material_id;

          // transform from ucd to
          // consecutive numbering
          for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
            if (vertex_indices.find(cells.back().vertices[i]) !=
                vertex_indices.end())
              // vertex with this index exists
              cells.back().vertices[i] =
                vertex_indices[cells.back().vertices[i]];
            else
              {
                // no such vertex index
                AssertThrow(false,
                            ExcInvalidVertexIndex(cell,
                                                  cells.back().vertices[i]));

                cells.back().vertices[i] = numbers::invalid_unsigned_int;
              }
        }
      else if ((cell_type == "line") && ((dim == 2) || (dim == 3)))
        // boundary info
        {
          subcelldata.boundary_lines.emplace_back();
          in >> subcelldata.boundary_lines.back().vertices[0] >>
            subcelldata.boundary_lines.back().vertices[1];

          // to make sure that the cast won't fail
          Assert(material_id <= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,
                               0,
                               std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to
          // numbers::internal_face_boundary_id-1
          AssertIndexRange(material_id, numbers::internal_face_boundary_id);

          // Make sure to set both manifold id and boundary id appropriately in
          // both cases:
          // numbers::internal_face_boundary_id and numbers::flat_manifold_id
          // are ignored in Triangulation::create_triangulation.
          if (apply_all_indicators_to_manifolds)
            {
              subcelldata.boundary_lines.back().boundary_id =
                numbers::internal_face_boundary_id;
              subcelldata.boundary_lines.back().manifold_id =
                static_cast<types::manifold_id>(material_id);
            }
          else
            {
              subcelldata.boundary_lines.back().boundary_id =
                static_cast<types::boundary_id>(material_id);
              subcelldata.boundary_lines.back().manifold_id =
                numbers::flat_manifold_id;
            }

          // transform from ucd to
          // consecutive numbering
          for (unsigned int &vertex :
               subcelldata.boundary_lines.back().vertices)
            if (vertex_indices.find(vertex) != vertex_indices.end())
              // vertex with this index exists
              vertex = vertex_indices[vertex];
            else
              {
                // no such vertex index
                AssertThrow(false, ExcInvalidVertexIndex(cell, vertex));
                vertex = numbers::invalid_unsigned_int;
              }
        }
      else if ((cell_type == "quad") && (dim == 3))
        // boundary info
        {
          subcelldata.boundary_quads.emplace_back();
          for (const unsigned int i : GeometryInfo<2>::vertex_indices())
            in >> subcelldata.boundary_quads.back()
                    .vertices[GeometryInfo<2>::ucd_to_deal[i]];

          // to make sure that the cast won't fail
          Assert(material_id <= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,
                               0,
                               std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to
          // numbers::internal_face_boundary_id-1
          AssertIndexRange(material_id, numbers::internal_face_boundary_id);

          // Make sure to set both manifold id and boundary id appropriately in
          // both cases:
          // numbers::internal_face_boundary_id and numbers::flat_manifold_id
          // are ignored in Triangulation::create_triangulation.
          if (apply_all_indicators_to_manifolds)
            {
              subcelldata.boundary_quads.back().boundary_id =
                numbers::internal_face_boundary_id;
              subcelldata.boundary_quads.back().manifold_id =
                static_cast<types::manifold_id>(material_id);
            }
          else
            {
              subcelldata.boundary_quads.back().boundary_id =
                static_cast<types::boundary_id>(material_id);
              subcelldata.boundary_quads.back().manifold_id =
                numbers::flat_manifold_id;
            }

          // transform from ucd to
          // consecutive numbering
          for (unsigned int &vertex :
               subcelldata.boundary_quads.back().vertices)
            if (vertex_indices.find(vertex) != vertex_indices.end())
              // vertex with this index exists
              vertex = vertex_indices[vertex];
            else
              {
                // no such vertex index
                Assert(false, ExcInvalidVertexIndex(cell, vertex));
                vertex = numbers::invalid_unsigned_int;
              }
        }
      else
        // cannot read this
        AssertThrow(false, ExcUnknownIdentifier(cell_type));
    }


  // check that no forbidden arrays are used
  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow(in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  // ... and cells
  if (dim == spacedim)
    GridTools::invert_all_negative_measure_cells(vertices, cells);
  GridTools::consistently_order_cells(cells);
  tria->create_triangulation(vertices, cells, subcelldata);
}

namespace
{
  template <int dim, int spacedim>
  class Abaqus_to_UCD
  {
  public:
    Abaqus_to_UCD();

    void
    read_in_abaqus(std::istream &in);
    void
    write_out_avs_ucd(std::ostream &out) const;

  private:
    const double tolerance;

    std::vector<double>
    get_global_node_numbers(const int face_cell_no,
                            const int face_cell_face_no) const;

    // NL: Stored as [ global node-id (int), x-coord, y-coord, z-coord ]
    std::vector<std::vector<double>> node_list;
    // CL: Stored as [ material-id (int), node1, node2, node3, node4, node5,
    // node6, node7, node8 ]
    std::vector<std::vector<double>> cell_list;
    // FL: Stored as [ sideset-id (int), node1, node2, node3, node4 ]
    std::vector<std::vector<double>> face_list;
    // ELSET: Stored as [ (std::string) elset_name = (std::vector) of cells
    // numbers]
    std::map<std::string, std::vector<int>> elsets_list;
  };
} // namespace

template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_abaqus(std::istream &in,
                                   const bool apply_all_indicators_to_manifolds)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  // This implementation has only been verified for:
  // - 2d grids with codimension 0
  // - 3d grids with codimension 0
  // - 3d grids with codimension 1
  Assert((spacedim == 2 && dim == spacedim) ||
           (spacedim == 3 && (dim == spacedim || dim == spacedim - 1)),
         ExcNotImplemented());
  AssertThrow(in, ExcIO());

  // Read in the Abaqus file into an intermediate object
  // that is to be passed along to the UCD reader
  Abaqus_to_UCD<dim, spacedim> abaqus_to_ucd;
  abaqus_to_ucd.read_in_abaqus(in);

  std::stringstream in_ucd;
  abaqus_to_ucd.write_out_avs_ucd(in_ucd);

  // This next call is wrapped in a try-catch for the following reason:
  // It ensures that if the Abaqus mesh is read in correctly but produces
  // an erroneous result then the user is alerted to the source of the problem
  // and doesn't think that they've somehow called the wrong function.
  try
    {
      read_ucd(in_ucd, apply_all_indicators_to_manifolds);
    }
  catch (std::exception &exc)
    {
      std::cerr << "Exception on processing internal UCD data: " << std::endl
                << exc.what() << std::endl;

      AssertThrow(
        false,
        ExcMessage(
          "Internal conversion from ABAQUS file to UCD format was unsuccessful. "
          "More information is provided in an error message printed above. "
          "Are you sure that your ABAQUS mesh file conforms with the requirements "
          "listed in the documentation?"));
    }
  catch (...)
    {
      AssertThrow(
        false,
        ExcMessage(
          "Internal conversion from ABAQUS file to UCD format was unsuccessful. "
          "Are you sure that your ABAQUS mesh file conforms with the requirements "
          "listed in the documentation?"));
    }
}


template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_dbmesh(std::istream &in)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  Assert(dim == 2, ExcNotImplemented());

  AssertThrow(in, ExcIO());

  // skip comments at start of file
  skip_comment_lines(in, '#');

  // first read in identifier string
  std::string line;
  getline(in, line);

  AssertThrow(line == "MeshVersionFormatted 0", ExcInvalidDBMESHInput(line));

  skip_empty_lines(in);

  // next read dimension
  getline(in, line);
  AssertThrow(line == "Dimension", ExcInvalidDBMESHInput(line));
  unsigned int dimension;
  in >> dimension;
  AssertThrow(dimension == dim, ExcDBMESHWrongDimension(dimension));
  skip_empty_lines(in);

  // now there are a lot of fields of
  // which we don't know the exact
  // meaning and which are far from
  // being properly documented in the
  // manual. we skip everything until
  // we find a comment line with the
  // string "# END". at some point in
  // the future, someone may have the
  // knowledge to parse and interpret
  // the other fields in between as
  // well...
  while (getline(in, line), line.find("# END") == std::string::npos)
    ;
  skip_empty_lines(in);


  // now read vertices
  getline(in, line);
  AssertThrow(line == "Vertices", ExcInvalidDBMESHInput(line));

  unsigned int n_vertices;
  double       dummy;

  in >> n_vertices;
  std::vector<Point<spacedim>> vertices(n_vertices);
  for (unsigned int vertex = 0; vertex < n_vertices; ++vertex)
    {
      // read vertex coordinates
      for (unsigned int d = 0; d < dim; ++d)
        in >> vertices[vertex][d];
      // read Ref phi_i, whatever that may be
      in >> dummy;
    }
  AssertThrow(in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);

  // read edges. we ignore them at
  // present, so just read them and
  // discard the input
  getline(in, line);
  AssertThrow(line == "Edges", ExcInvalidDBMESHInput(line));

  unsigned int n_edges;
  in >> n_edges;
  for (unsigned int edge = 0; edge < n_edges; ++edge)
    {
      // read vertex indices
      in >> dummy >> dummy;
      // read Ref phi_i, whatever that may be
      in >> dummy;
    }
  AssertThrow(in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);



  // read cracked edges (whatever
  // that may be). we ignore them at
  // present, so just read them and
  // discard the input
  getline(in, line);
  AssertThrow(line == "CrackedEdges", ExcInvalidDBMESHInput(line));

  in >> n_edges;
  for (unsigned int edge = 0; edge < n_edges; ++edge)
    {
      // read vertex indices
      in >> dummy >> dummy;
      // read Ref phi_i, whatever that may be
      in >> dummy;
    }
  AssertThrow(in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);


  // now read cells.
  // set up array of cells
  getline(in, line);
  AssertThrow(line == "Quadrilaterals", ExcInvalidDBMESHInput(line));

  std::vector<CellData<dim>> cells;
  SubCellData                subcelldata;
  unsigned int               n_cells;
  in >> n_cells;
  for (unsigned int cell = 0; cell < n_cells; ++cell)
    {
      // read in vertex numbers. they
      // are 1-based, so subtract one
      cells.emplace_back();
      for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
        {
          in >> cells.back().vertices[GeometryInfo<dim>::ucd_to_deal[i]];

          AssertThrow((cells.back().vertices[i] >= 1) &&
                        (static_cast<unsigned int>(cells.back().vertices[i]) <=
                         vertices.size()),
                      ExcInvalidVertexIndex(cell, cells.back().vertices[i]));

          --cells.back().vertices[i];
        }

      // read and discard Ref phi_i
      in >> dummy;
    }
  AssertThrow(in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);


  // then there are again a whole lot
  // of fields of which I have no
  // clue what they mean. skip them
  // all and leave the interpretation
  // to other implementors...
  while (getline(in, line), ((line.find("End") == std::string::npos) && (in)))
    ;
  // ok, so we are not at the end of
  // the file, that's it, mostly


  // check that no forbidden arrays are used
  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow(in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  // ...and cells
  GridTools::invert_all_negative_measure_cells(vertices, cells);
  GridTools::consistently_order_cells(cells);
  tria->create_triangulation(vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_xda(std::istream &in)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  AssertThrow(in, ExcIO());

  const auto reference_cell = ReferenceCells::get_hypercube<dim>();

  std::string line;
  // skip comments at start of file
  std::getline(in, line);

  unsigned int n_vertices;
  unsigned int n_cells;

  // read cells, throw away rest of line
  in >> n_cells;
  std::getline(in, line);

  in >> n_vertices;
  std::getline(in, line);

  // ignore following 8 lines
  for (unsigned int i = 0; i < 8; ++i)
    std::getline(in, line);

  // set up array of cells
  std::vector<CellData<dim>> cells(n_cells);
  SubCellData                subcelldata;

  for (CellData<dim> &cell : cells)
    {
      // note that since in the input file we found the number of cells at the
      // top, there should still be input here, so check this:
      AssertThrow(in, ExcIO());

      // XDA happens to use ExodusII's numbering because XDA/XDR is libMesh's
      // native format, and libMesh's node numberings come from ExodusII:
      for (unsigned int i = 0; i < GeometryInfo<dim>::vertices_per_cell; i++)
        in >> cell.vertices[reference_cell.exodusii_vertex_to_deal_vertex(i)];
    }

  // set up array of vertices
  std::vector<Point<spacedim>> vertices(n_vertices);
  for (Point<spacedim> &vertex : vertices)
    {
      for (unsigned int d = 0; d < spacedim; ++d)
        in >> vertex[d];
      for (unsigned int d = spacedim; d < 3; ++d)
        {
          // file is always in 3D
          double dummy;
          in >> dummy;
        }
    }
  AssertThrow(in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  // ... and cells
  GridTools::invert_all_negative_measure_cells(vertices, cells);
  GridTools::consistently_order_cells(cells);
  tria->create_triangulation(vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_msh(std::istream &in)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  AssertThrow(in, ExcIO());

  unsigned int n_vertices;
  unsigned int n_cells;
  unsigned int dummy;
  std::string  line;
  // This array stores maps from the 'entities' to the 'physical tags' for
  // points, curves, surfaces and volumes. We use this information later to
  // assign boundary ids.
  std::array<std::map<int, int>, 4> tag_maps;

  in >> line;

  // first determine file format
  unsigned int gmsh_file_format = 0;
  if (line == "$NOD")
    gmsh_file_format = 10;
  else if (line == "$MeshFormat")
    gmsh_file_format = 20;
  else
    AssertThrow(false, ExcInvalidGMSHInput(line));

  // if file format is 2.0 or greater then we also have to read the rest of the
  // header
  if (gmsh_file_format == 20)
    {
      double       version;
      unsigned int file_type, data_size;

      in >> version >> file_type >> data_size;

      Assert((version >= 2.0) && (version <= 4.1), ExcNotImplemented());
      gmsh_file_format = static_cast<unsigned int>(version * 10);

      Assert(file_type == 0, ExcNotImplemented());
      Assert(data_size == sizeof(double), ExcNotImplemented());

      // read the end of the header and the first line of the nodes description
      // to synch ourselves with the format 1 handling above
      in >> line;
      AssertThrow(line == "$EndMeshFormat", ExcInvalidGMSHInput(line));

      in >> line;
      // if the next block is of kind $PhysicalNames, ignore it
      if (line == "$PhysicalNames")
        {
          do
            {
              in >> line;
            }
          while (line != "$EndPhysicalNames");
          in >> line;
        }

      // if the next block is of kind $Entities, parse it
      if (line == "$Entities")
        {
          unsigned long n_points, n_curves, n_surfaces, n_volumes;

          in >> n_points >> n_curves >> n_surfaces >> n_volumes;
          for (unsigned int i = 0; i < n_points; ++i)
            {
              // parse point ids
              int          tag;
              unsigned int n_physicals;
              double box_min_x, box_min_y, box_min_z, box_max_x, box_max_y,
                box_max_z;

              // we only care for 'tag' as key for tag_maps[0]
              if (gmsh_file_format > 40)
                {
                  in >> tag >> box_min_x >> box_min_y >> box_min_z >>
                    n_physicals;
                  box_max_x = box_min_x;
                  box_max_y = box_min_y;
                  box_max_z = box_min_z;
                }
              else
                {
                  in >> tag >> box_min_x >> box_min_y >> box_min_z >>
                    box_max_x >> box_max_y >> box_max_z >> n_physicals;
                }
              // if there is a physical tag, we will use it as boundary id below
              AssertThrow(n_physicals < 2,
                          ExcMessage("More than one tag is not supported!"));
              // if there is no physical tag, use 0 as default
              int physical_tag = 0;
              for (unsigned int j = 0; j < n_physicals; ++j)
                in >> physical_tag;
              tag_maps[0][tag] = physical_tag;
            }
          for (unsigned int i = 0; i < n_curves; ++i)
            {
              // parse curve ids
              int          tag;
              unsigned int n_physicals;
              double box_min_x, box_min_y, box_min_z, box_max_x, box_max_y,
                box_max_z;

              // we only care for 'tag' as key for tag_maps[1]
              in >> tag >> box_min_x >> box_min_y >> box_min_z >> box_max_x >>
                box_max_y >> box_max_z >> n_physicals;
              // if there is a physical tag, we will use it as boundary id below
              AssertThrow(n_physicals < 2,
                          ExcMessage("More than one tag is not supported!"));
              // if there is no physical tag, use 0 as default
              int physical_tag = 0;
              for (unsigned int j = 0; j < n_physicals; ++j)
                in >> physical_tag;
              tag_maps[1][tag] = physical_tag;
              // we don't care about the points associated to a curve, but have
              // to parse them anyway because their format is unstructured
              in >> n_points;
              for (unsigned int j = 0; j < n_points; ++j)
                in >> tag;
            }

          for (unsigned int i = 0; i < n_surfaces; ++i)
            {
              // parse surface ids
              int          tag;
              unsigned int n_physicals;
              double box_min_x, box_min_y, box_min_z, box_max_x, box_max_y,
                box_max_z;

              // we only care for 'tag' as key for tag_maps[2]
              in >> tag >> box_min_x >> box_min_y >> box_min_z >> box_max_x >>
                box_max_y >> box_max_z >> n_physicals;
              // if there is a physical tag, we will use it as boundary id below
              AssertThrow(n_physicals < 2,
                          ExcMessage("More than one tag is not supported!"));
              // if there is no physical tag, use 0 as default
              int physical_tag = 0;
              for (unsigned int j = 0; j < n_physicals; ++j)
                in >> physical_tag;
              tag_maps[2][tag] = physical_tag;
              // we don't care about the curves associated to a surface, but
              // have to parse them anyway because their format is unstructured
              in >> n_curves;
              for (unsigned int j = 0; j < n_curves; ++j)
                in >> tag;
            }
          for (unsigned int i = 0; i < n_volumes; ++i)
            {
              // parse volume ids
              int          tag;
              unsigned int n_physicals;
              double box_min_x, box_min_y, box_min_z, box_max_x, box_max_y,
                box_max_z;

              // we only care for 'tag' as key for tag_maps[3]
              in >> tag >> box_min_x >> box_min_y >> box_min_z >> box_max_x >>
                box_max_y >> box_max_z >> n_physicals;
              // if there is a physical tag, we will use it as boundary id below
              AssertThrow(n_physicals < 2,
                          ExcMessage("More than one tag is not supported!"));
              // if there is no physical tag, use 0 as default
              int physical_tag = 0;
              for (unsigned int j = 0; j < n_physicals; ++j)
                in >> physical_tag;
              tag_maps[3][tag] = physical_tag;
              // we don't care about the surfaces associated to a volume, but
              // have to parse them anyway because their format is unstructured
              in >> n_surfaces;
              for (unsigned int j = 0; j < n_surfaces; ++j)
                in >> tag;
            }
          in >> line;
          AssertThrow(line == "$EndEntities", ExcInvalidGMSHInput(line));
          in >> line;
        }

      // if the next block is of kind $PartitionedEntities, ignore it
      if (line == "$PartitionedEntities")
        {
          do
            {
              in >> line;
            }
          while (line != "$EndPartitionedEntities");
          in >> line;
        }

      // but the next thing should,
      // in any case, be the list of
      // nodes:
      AssertThrow(line == "$Nodes", ExcInvalidGMSHInput(line));
    }

  // now read the nodes list
  int n_entity_blocks = 1;
  if (gmsh_file_format > 40)
    {
      int min_node_tag;
      int max_node_tag;
      in >> n_entity_blocks >> n_vertices >> min_node_tag >> max_node_tag;
    }
  else if (gmsh_file_format == 40)
    {
      in >> n_entity_blocks >> n_vertices;
    }
  else
    in >> n_vertices;
  std::vector<Point<spacedim>> vertices(n_vertices);
  // set up mapping between numbering
  // in msh-file (nod) and in the
  // vertices vector
  std::map<int, int> vertex_indices;

  {
    unsigned int global_vertex = 0;
    for (int entity_block = 0; entity_block < n_entity_blocks; ++entity_block)
      {
        int           parametric;
        unsigned long numNodes;

        if (gmsh_file_format < 40)
          {
            numNodes   = n_vertices;
            parametric = 0;
          }
        else
          {
            // for gmsh_file_format 4.1 the order of tag and dim is reversed,
            // but we are ignoring both anyway.
            int tagEntity, dimEntity;
            in >> tagEntity >> dimEntity >> parametric >> numNodes;
          }

        std::vector<int> vertex_numbers;
        int              vertex_number;
        if (gmsh_file_format > 40)
          for (unsigned long vertex_per_entity = 0;
               vertex_per_entity < numNodes;
               ++vertex_per_entity)
            {
              in >> vertex_number;
              vertex_numbers.push_back(vertex_number);
            }

        for (unsigned long vertex_per_entity = 0; vertex_per_entity < numNodes;
             ++vertex_per_entity, ++global_vertex)
          {
            int    vertex_number;
            double x[3];

            // read vertex
            if (gmsh_file_format > 40)
              {
                vertex_number = vertex_numbers[vertex_per_entity];
                in >> x[0] >> x[1] >> x[2];
              }
            else
              in >> vertex_number >> x[0] >> x[1] >> x[2];

            for (unsigned int d = 0; d < spacedim; ++d)
              vertices[global_vertex](d) = x[d];
            // store mapping
            vertex_indices[vertex_number] = global_vertex;

            // ignore parametric coordinates
            if (parametric != 0)
              {
                double u = 0.;
                double v = 0.;
                in >> u >> v;
                (void)u;
                (void)v;
              }
          }
      }
    AssertDimension(global_vertex, n_vertices);
  }

  // Assert we reached the end of the block
  in >> line;
  static const std::string end_nodes_marker[] = {"$ENDNOD", "$EndNodes"};
  AssertThrow(line == end_nodes_marker[gmsh_file_format == 10 ? 0 : 1],
              ExcInvalidGMSHInput(line));

  // Now read in next bit
  in >> line;
  static const std::string begin_elements_marker[] = {"$ELM", "$Elements"};
  AssertThrow(line == begin_elements_marker[gmsh_file_format == 10 ? 0 : 1],
              ExcInvalidGMSHInput(line));

  // now read the cell list
  if (gmsh_file_format > 40)
    {
      int min_node_tag;
      int max_node_tag;
      in >> n_entity_blocks >> n_cells >> min_node_tag >> max_node_tag;
    }
  else if (gmsh_file_format == 40)
    {
      in >> n_entity_blocks >> n_cells;
    }
  else
    {
      n_entity_blocks = 1;
      in >> n_cells;
    }

  // set up array of cells and subcells (faces). In 1d, there is currently no
  // standard way in deal.II to pass boundary indicators attached to individual
  // vertices, so do this by hand via the boundary_ids_1d array
  std::vector<CellData<dim>>                 cells;
  SubCellData                                subcelldata;
  std::map<unsigned int, types::boundary_id> boundary_ids_1d;
  bool                                       is_quad_or_hex_mesh = false;
  bool                                       is_tria_or_tet_mesh = false;

  {
    unsigned int global_cell = 0;
    for (int entity_block = 0; entity_block < n_entity_blocks; ++entity_block)
      {
        unsigned int  material_id;
        unsigned long numElements;
        int           cell_type;

        if (gmsh_file_format < 40)
          {
            material_id = 0;
            cell_type   = 0;
            numElements = n_cells;
          }
        else if (gmsh_file_format == 40)
          {
            int tagEntity, dimEntity;
            in >> tagEntity >> dimEntity >> cell_type >> numElements;
            material_id = tag_maps[dimEntity][tagEntity];
          }
        else
          {
            // for gmsh_file_format 4.1 the order of tag and dim is reversed,
            int tagEntity, dimEntity;
            in >> dimEntity >> tagEntity >> cell_type >> numElements;
            material_id = tag_maps[dimEntity][tagEntity];
          }

        for (unsigned int cell_per_entity = 0; cell_per_entity < numElements;
             ++cell_per_entity, ++global_cell)
          {
            // note that since in the input
            // file we found the number of
            // cells at the top, there
            // should still be input here,
            // so check this:
            AssertThrow(in, ExcIO());

            unsigned int nod_num;

            /*
              For file format version 1, the format of each cell is as follows:
                elm-number elm-type reg-phys reg-elem number-of-nodes
              node-number-list

              However, for version 2, the format reads like this:
                elm-number elm-type number-of-tags < tag > ... node-number-list

              For version 4, we have:
                tag(int) numVert(int) ...

              In the following, we will ignore the element number (we simply
              enumerate them in the order in which we read them, and we will
              take reg-phys (version 1) or the first tag (version 2, if any tag
              is given at all) as material id. For version 4, we already read
              the material and the cell type in above.
            */

            unsigned int elm_number = 0;
            if (gmsh_file_format < 40)
              {
                in >> elm_number // ELM-NUMBER
                  >> cell_type;  // ELM-TYPE
              }

            if (gmsh_file_format < 20)
              {
                in >> material_id // REG-PHYS
                  >> dummy        // reg_elm
                  >> nod_num;
              }
            else if (gmsh_file_format < 40)
              {
                // read the tags; ignore all but the first one which we will
                // interpret as the material_id (for cells) or boundary_id
                // (for faces)
                unsigned int n_tags;
                in >> n_tags;
                if (n_tags > 0)
                  in >> material_id;
                else
                  material_id = 0;

                for (unsigned int i = 1; i < n_tags; ++i)
                  in >> dummy;

                if (cell_type == 1) // line
                  nod_num = 2;
                else if (cell_type == 2) // tri
                  nod_num = 3;
                else if (cell_type == 3) // quad
                  nod_num = 4;
                else if (cell_type == 4) // tet
                  nod_num = 4;
                else if (cell_type == 5) // hex
                  nod_num = 8;
              }
            else
              {
                // ignore tag
                int tag;
                in >> tag;

                if (cell_type == 1) // line
                  nod_num = 2;
                else if (cell_type == 2) // tri
                  nod_num = 3;
                else if (cell_type == 3) // quad
                  nod_num = 4;
                else if (cell_type == 4) // tet
                  nod_num = 4;
                else if (cell_type == 5) // hex
                  nod_num = 8;
              }


            /*       `ELM-TYPE'
                     defines the geometrical type of the N-th element:
                     `1'
                     Line (2 nodes, 1 edge).

                     `2'
                     Triangle (3 nodes, 3 edges).

                     `3'
                     Quadrangle (4 nodes, 4 edges).

                     `4'
                     Tetrahedron (4 nodes, 6 edges, 6 faces).

                     `5'
                     Hexahedron (8 nodes, 12 edges, 6 faces).

                     `15'
                     Point (1 node).
            */

            if (((cell_type == 1) && (dim == 1)) ||
                ((cell_type == 2) && (dim == 2)) ||
                ((cell_type == 3) && (dim == 2)) ||
                ((cell_type == 4) && (dim == 3)) ||
                ((cell_type == 5) && (dim == 3)))
              // found a cell
              {
                unsigned int vertices_per_cell = 0;
                if (cell_type == 1) // line
                  vertices_per_cell = 2;
                else if (cell_type == 2) // tri
                  {
                    vertices_per_cell   = 3;
                    is_tria_or_tet_mesh = true;
                  }
                else if (cell_type == 3) // quad
                  {
                    vertices_per_cell   = 4;
                    is_quad_or_hex_mesh = true;
                  }
                else if (cell_type == 4) // tet
                  {
                    vertices_per_cell   = 4;
                    is_tria_or_tet_mesh = true;
                  }
                else if (cell_type == 5) // hex
                  {
                    vertices_per_cell   = 8;
                    is_quad_or_hex_mesh = true;
                  }

                AssertThrow(nod_num == vertices_per_cell,
                            ExcMessage(
                              "Number of nodes does not coincide with the "
                              "number required for this object"));

                // allocate and read indices
                cells.emplace_back();
                cells.back().vertices.resize(vertices_per_cell);
                for (unsigned int i = 0; i < vertices_per_cell; ++i)
                  {
                    // hypercube cells need to be reordered
                    if (vertices_per_cell ==
                        GeometryInfo<dim>::vertices_per_cell)
                      {
                        in >> cells.back()
                                .vertices[GeometryInfo<dim>::ucd_to_deal[i]];
                      }
                    else
                      {
                        in >> cells.back().vertices[i];
                      }
                  }

                // to make sure that the cast won't fail
                Assert(material_id <=
                         std::numeric_limits<types::material_id>::max(),
                       ExcIndexRange(
                         material_id,
                         0,
                         std::numeric_limits<types::material_id>::max()));
                // we use only material_ids in the range from 0 to
                // numbers::invalid_material_id-1
                AssertIndexRange(material_id, numbers::invalid_material_id);

                cells.back().material_id = material_id;

                // transform from gmsh to consecutive numbering
                for (unsigned int i = 0; i < vertices_per_cell; ++i)
                  {
                    AssertThrow(
                      vertex_indices.find(cells.back().vertices[i]) !=
                        vertex_indices.end(),
                      ExcInvalidVertexIndexGmsh(cell_per_entity,
                                                elm_number,
                                                cells.back().vertices[i]));

                    // vertex with this index exists
                    cells.back().vertices[i] =
                      vertex_indices[cells.back().vertices[i]];
                  }
              }
            else if ((cell_type == 1) && ((dim == 2) || (dim == 3)))
              // boundary info
              {
                subcelldata.boundary_lines.emplace_back();
                in >> subcelldata.boundary_lines.back().vertices[0] >>
                  subcelldata.boundary_lines.back().vertices[1];

                // to make sure that the cast won't fail
                Assert(material_id <=
                         std::numeric_limits<types::boundary_id>::max(),
                       ExcIndexRange(
                         material_id,
                         0,
                         std::numeric_limits<types::boundary_id>::max()));
                // we use only boundary_ids in the range from 0 to
                // numbers::internal_face_boundary_id-1
                AssertIndexRange(material_id,
                                 numbers::internal_face_boundary_id);

                subcelldata.boundary_lines.back().boundary_id =
                  static_cast<types::boundary_id>(material_id);

                // transform from ucd to
                // consecutive numbering
                for (unsigned int &vertex :
                     subcelldata.boundary_lines.back().vertices)
                  if (vertex_indices.find(vertex) != vertex_indices.end())
                    // vertex with this index exists
                    vertex = vertex_indices[vertex];
                  else
                    {
                      // no such vertex index
                      AssertThrow(false,
                                  ExcInvalidVertexIndex(cell_per_entity,
                                                        vertex));
                      vertex = numbers::invalid_unsigned_int;
                    }
              }
            else if ((cell_type == 2 || cell_type == 3) && (dim == 3))
              // boundary info
              {
                unsigned int vertices_per_cell = 0;
                // check cell type
                if (cell_type == 2) // tri
                  {
                    vertices_per_cell   = 3;
                    is_tria_or_tet_mesh = true;
                  }
                else if (cell_type == 3) // quad
                  {
                    vertices_per_cell   = 4;
                    is_quad_or_hex_mesh = true;
                  }

                subcelldata.boundary_quads.emplace_back();

                // resize vertices
                subcelldata.boundary_quads.back().vertices.resize(
                  vertices_per_cell);
                // for loop
                for (unsigned int i = 0; i < vertices_per_cell; ++i)
                  in >> subcelldata.boundary_quads.back().vertices[i];

                // to make sure that the cast won't fail
                Assert(material_id <=
                         std::numeric_limits<types::boundary_id>::max(),
                       ExcIndexRange(
                         material_id,
                         0,
                         std::numeric_limits<types::boundary_id>::max()));
                // we use only boundary_ids in the range from 0 to
                // numbers::internal_face_boundary_id-1
                AssertIndexRange(material_id,
                                 numbers::internal_face_boundary_id);

                subcelldata.boundary_quads.back().boundary_id =
                  static_cast<types::boundary_id>(material_id);

                // transform from gmsh to
                // consecutive numbering
                for (unsigned int &vertex :
                     subcelldata.boundary_quads.back().vertices)
                  if (vertex_indices.find(vertex) != vertex_indices.end())
                    // vertex with this index exists
                    vertex = vertex_indices[vertex];
                  else
                    {
                      // no such vertex index
                      Assert(false,
                             ExcInvalidVertexIndex(cell_per_entity, vertex));
                      vertex = numbers::invalid_unsigned_int;
                    }
              }
            else if (cell_type == 15)
              {
                // read the indices of nodes given
                unsigned int node_index = 0;
                if (gmsh_file_format < 20)
                  {
                    // For points (cell_type==15), we can only ever
                    // list one node index.
                    AssertThrow(nod_num == 1, ExcInternalError());
                    in >> node_index;
                  }
                else
                  {
                    in >> node_index;
                  }

                // we only care about boundary indicators assigned to individual
                // vertices in 1d (because otherwise the vertices are not faces)
                if (dim == 1)
                  boundary_ids_1d[vertex_indices[node_index]] = material_id;
              }
            else
              {
                AssertThrow(false, ExcGmshUnsupportedGeometry(cell_type));
              }
          }
      }
    AssertDimension(global_cell, n_cells);
  }
  // Assert we reached the end of the block
  in >> line;
  static const std::string end_elements_marker[] = {"$ENDELM", "$EndElements"};
  AssertThrow(line == end_elements_marker[gmsh_file_format == 10 ? 0 : 1],
              ExcInvalidGMSHInput(line));

  // check that no forbidden arrays are used
  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow(in, ExcIO());

  // check that we actually read some cells.
  AssertThrow(cells.size() > 0, ExcGmshNoCellInformation());

  // TODO: the functions below (GridTools::delete_unused_vertices(),
  // GridTools::invert_all_negative_measure_cells(),
  // GridTools::consistently_order_cells()) need to be revisited
  // for simplex/mixed meshes

  if (dim == 1 || (is_quad_or_hex_mesh && !is_tria_or_tet_mesh))
    {
      // do some clean-up on vertices...
      GridTools::delete_unused_vertices(vertices, cells, subcelldata);
      // ... and cells
      if (dim == spacedim)
        GridTools::invert_all_negative_measure_cells(vertices, cells);
      GridTools::consistently_order_cells(cells);
    }
  tria->create_triangulation(vertices, cells, subcelldata);

  // in 1d, we also have to attach boundary ids to vertices, which does not
  // currently work through the call above
  if (dim == 1)
    assign_1d_boundary_ids(boundary_ids_1d, *tria);
}



#ifdef DEAL_II_GMSH_WITH_API
template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_msh(const std::string &fname)
{
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  // gmsh -> deal.II types
  const std::map<int, std::uint8_t> gmsh_to_dealii_type = {
    {{15, 0}, {1, 1}, {2, 2}, {3, 3}, {4, 4}, {7, 5}, {6, 6}, {5, 7}}};

  // Vertex renumbering, by dealii type
  const std::array<std::vector<unsigned int>, 8> gmsh_to_dealii = {
    {{0},
     {{0, 1}},
     {{0, 1, 2}},
     {{0, 1, 3, 2}},
     {{0, 1, 2, 3}},
     {{0, 1, 3, 2, 4}},
     {{0, 1, 2, 3, 4, 5}},
     {{0, 1, 3, 2, 4, 5, 7, 6}}}};

  std::vector<Point<spacedim>>               vertices;
  std::vector<CellData<dim>>                 cells;
  SubCellData                                subcelldata;
  std::map<unsigned int, types::boundary_id> boundary_ids_1d;

  gmsh::initialize();
  gmsh::option::setNumber("General.Verbosity", 0);
  gmsh::open(fname);

  AssertThrow(gmsh::model::getDimension() == dim,
              ExcMessage("You are trying to read a gmsh file with dimension " +
                         std::to_string(gmsh::model::getDimension()) +
                         " into a grid of dimension " + std::to_string(dim)));

  // Read all nodes, and store them in our vector of vertices. Before we do
  // that, make sure all tags are consecutive
  {
    gmsh::model::mesh::removeDuplicateNodes();
    gmsh::model::mesh::renumberNodes();
    std::vector<std::size_t> node_tags;
    std::vector<double>      coord;
    std::vector<double>      parametricCoord;
    gmsh::model::mesh::getNodes(
      node_tags, coord, parametricCoord, -1, -1, false, false);
    vertices.resize(node_tags.size());
    for (unsigned int i = 0; i < node_tags.size(); ++i)
      {
        // Check that renumbering worked!
        AssertDimension(node_tags[i], i + 1);
        for (unsigned int d = 0; d < spacedim; ++d)
          vertices[i][d] = coord[i * 3 + d];
#  ifdef DEBUG
        // Make sure the embedded dimension is right
        for (unsigned int d = spacedim; d < 3; ++d)
          Assert(coord[i * 3 + d] == 0,
                 ExcMessage("The grid you are reading contains nodes that are "
                            "nonzero in the coordinate with index " +
                            std::to_string(d) +
                            ", but you are trying to save "
                            "it on a grid embedded in a " +
                            std::to_string(spacedim) + " dimensional space."));
#  endif
      }
  }

  // Get all the elementary entities in the model, as a vector of (dimension,
  // tag) pairs:
  std::vector<std::pair<int, int>> entities;
  gmsh::model::getEntities(entities);

  for (const auto &e : entities)
    {
      // Dimension and tag of the entity:
      const int &entity_dim = e.first;
      const int &entity_tag = e.second;

      types::manifold_id manifold_id = numbers::flat_manifold_id;
      types::boundary_id boundary_id = 0;

      // Get the physical tags, to deduce boundary, material, and manifold_id
      std::vector<int> physical_tags;
      gmsh::model::getPhysicalGroupsForEntity(entity_dim,
                                              entity_tag,
                                              physical_tags);

      // Now fill manifold id and boundary or material id
      if (physical_tags.size())
        for (auto physical_tag : physical_tags)
          {
            std::string name;
            gmsh::model::getPhysicalName(entity_dim, physical_tag, name);
            if (!name.empty())
              try
                {
                  std::map<std::string, int> id_names;
                  Patterns::Tools::to_value(name, id_names);
                  bool throw_anyway      = false;
                  bool found_boundary_id = false;
                  // If the above did not throw, we keep going, and retrieve
                  // all the information that we know how to translate.
                  for (const auto &it : id_names)
                    {
                      const auto &name = it.first;
                      const auto &id   = it.second;
                      if (entity_dim == dim && name == "MaterialID")
                        {
                          boundary_id = static_cast<types::boundary_id>(id);
                          found_boundary_id = true;
                        }
                      else if (entity_dim < dim && name == "BoundaryID")
                        {
                          boundary_id = static_cast<types::boundary_id>(id);
                          found_boundary_id = true;
                        }
                      else if (name == "ManifoldID")
                        manifold_id = static_cast<types::manifold_id>(id);
                      else
                        // We did not recognize one of the keys. We'll fall
                        // back to setting the boundary id to the physical tag
                        // after reading all strings.
                        throw_anyway = true;
                    }
                  // If we didn't find a BoundaryID:XX or MaterialID:XX, and
                  // something was found but not recognized, then we set the
                  // material id or boundary id in the catch block below, using
                  // directly the physical tag
                  if (throw_anyway && !found_boundary_id)
                    throw;
                }
              catch (...)
                {
                  // When the above didn't work, we revert to the old
                  // behaviour: the physical tag itself is interpreted either
                  // as a material_id or a boundary_id, and no manifold id is
                  // known
                  boundary_id = physical_tag;
                }
          }

      // Get the mesh elements for the entity (dim, tag):
      std::vector<int>                      element_types;
      std::vector<std::vector<std::size_t>> element_ids, element_nodes;
      gmsh::model::mesh::getElements(
        element_types, element_ids, element_nodes, entity_dim, entity_tag);

      for (unsigned int i = 0; i < element_types.size(); ++i)
        {
          const auto &type       = gmsh_to_dealii_type.at(element_types[i]);
          const auto  n_vertices = gmsh_to_dealii[type].size();
          const auto &elements   = element_ids[i];
          const auto &nodes      = element_nodes[i];
          for (unsigned int j = 0; j < elements.size(); ++j)
            {
              if (entity_dim == dim)
                {
                  cells.emplace_back(n_vertices);
                  auto &cell = cells.back();
                  for (unsigned int v = 0; v < n_vertices; ++v)
                    cell.vertices[v] =
                      nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;
                  cell.manifold_id = manifold_id;
                  cell.material_id = boundary_id;
                }
              else if (entity_dim == 2)
                {
                  subcelldata.boundary_quads.emplace_back(n_vertices);
                  auto &face = subcelldata.boundary_quads.back();
                  for (unsigned int v = 0; v < n_vertices; ++v)
                    face.vertices[v] =
                      nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;

                  face.manifold_id = manifold_id;
                  face.boundary_id = boundary_id;
                }
              else if (entity_dim == 1)
                {
                  subcelldata.boundary_lines.emplace_back(n_vertices);
                  auto &line = subcelldata.boundary_lines.back();
                  for (unsigned int v = 0; v < n_vertices; ++v)
                    line.vertices[v] =
                      nodes[n_vertices * j + gmsh_to_dealii[type][v]] - 1;

                  line.manifold_id = manifold_id;
                  line.boundary_id = boundary_id;
                }
              else if (entity_dim == 0)
                {
                  // This should only happen in one dimension.
                  AssertDimension(dim, 1);
                  for (unsigned int j = 0; j < elements.size(); ++j)
                    boundary_ids_1d[nodes[j] - 1] = boundary_id;
                }
            }
        }
    }

  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  tria->create_triangulation(vertices, cells, subcelldata);

  // in 1d, we also have to attach boundary ids to vertices, which does not
  // currently work through the call above
  if (dim == 1)
    assign_1d_boundary_ids(boundary_ids_1d, *tria);

  gmsh::clear();
  gmsh::finalize();
}
#endif



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::parse_tecplot_header(
  std::string &              header,
  std::vector<unsigned int> &tecplot2deal,
  unsigned int &             n_vars,
  unsigned int &             n_vertices,
  unsigned int &             n_cells,
  std::vector<unsigned int> &IJK,
  bool &                     structured,
  bool &                     blocked)
{
  Assert(tecplot2deal.size() == dim, ExcInternalError());
  Assert(IJK.size() == dim, ExcInternalError());
  // initialize the output variables
  n_vars     = 0;
  n_vertices = 0;
  n_cells    = 0;
  switch (dim)
    {
      case 3:
        IJK[2] = 0;
        DEAL_II_FALLTHROUGH;
      case 2:
        IJK[1] = 0;
        DEAL_II_FALLTHROUGH;
      case 1:
        IJK[0] = 0;
    }
  structured = true;
  blocked    = false;

  // convert the string to upper case
  std::transform(header.begin(), header.end(), header.begin(), ::toupper);

  // replace all tabs, commas, newlines by
  // whitespaces
  std::replace(header.begin(), header.end(), '\t', ' ');
  std::replace(header.begin(), header.end(), ',', ' ');
  std::replace(header.begin(), header.end(), '\n', ' ');

  // now remove whitespace in front of and
  // after '='
  std::string::size_type pos = header.find('=');

  while (pos != static_cast<std::string::size_type>(std::string::npos))
    if (header[pos + 1] == ' ')
      header.erase(pos + 1, 1);
    else if (header[pos - 1] == ' ')
      {
        header.erase(pos - 1, 1);
        --pos;
      }
    else
      pos = header.find('=', ++pos);

  // split the string into individual entries
  std::vector<std::string> entries =
    Utilities::break_text_into_lines(header, 1, ' ');

  // now go through the list and try to extract
  for (unsigned int i = 0; i < entries.size(); ++i)
    {
      if (Utilities::match_at_string_start(entries[i], "VARIABLES=\""))
        {
          ++n_vars;
          // we assume, that the first variable
          // is x or no coordinate at all (not y or z)
          if (Utilities::match_at_string_start(entries[i], "VARIABLES=\"X\""))
            {
              tecplot2deal[0] = 0;
            }
          ++i;
          while (entries[i][0] == '"')
            {
              if (entries[i] == "\"X\"")
                tecplot2deal[0] = n_vars;
              else if (entries[i] == "\"Y\"")
                {
                  // we assume, that y contains
                  // zero data in 1d, so do
                  // nothing
                  if (dim > 1)
                    tecplot2deal[1] = n_vars;
                }
              else if (entries[i] == "\"Z\"")
                {
                  // we assume, that z contains
                  // zero data in 1d and 2d, so
                  // do nothing
                  if (dim > 2)
                    tecplot2deal[2] = n_vars;
                }
              ++n_vars;
              ++i;
            }
          // set i back, so that the next
          // string is treated correctly
          --i;

          AssertThrow(
            n_vars >= dim,
            ExcMessage(
              "Tecplot file must contain at least one variable for each dimension"));
          for (unsigned int d = 1; d < dim; ++d)
            AssertThrow(
              tecplot2deal[d] > 0,
              ExcMessage(
                "Tecplot file must contain at least one variable for each dimension."));
        }
      else if (Utilities::match_at_string_start(entries[i], "ZONETYPE=ORDERED"))
        structured = true;
      else if (Utilities::match_at_string_start(entries[i],
                                                "ZONETYPE=FELINESEG") &&
               dim == 1)
        structured = false;
      else if (Utilities::match_at_string_start(entries[i],
                                                "ZONETYPE=FEQUADRILATERAL") &&
               dim == 2)
        structured = false;
      else if (Utilities::match_at_string_start(entries[i],
                                                "ZONETYPE=FEBRICK") &&
               dim == 3)
        structured = false;
      else if (Utilities::match_at_string_start(entries[i], "ZONETYPE="))
        // unsupported ZONETYPE
        {
          AssertThrow(false,
                      ExcMessage(
                        "The tecplot file contains an unsupported ZONETYPE."));
        }
      else if (Utilities::match_at_string_start(entries[i],
                                                "DATAPACKING=POINT"))
        blocked = false;
      else if (Utilities::match_at_string_start(entries[i],
                                                "DATAPACKING=BLOCK"))
        blocked = true;
      else if (Utilities::match_at_string_start(entries[i], "F=POINT"))
        {
          structured = true;
          blocked    = false;
        }
      else if (Utilities::match_at_string_start(entries[i], "F=BLOCK"))
        {
          structured = true;
          blocked    = true;
        }
      else if (Utilities::match_at_string_start(entries[i], "F=FEPOINT"))
        {
          structured = false;
          blocked    = false;
        }
      else if (Utilities::match_at_string_start(entries[i], "F=FEBLOCK"))
        {
          structured = false;
          blocked    = true;
        }
      else if (Utilities::match_at_string_start(entries[i],
                                                "ET=QUADRILATERAL") &&
               dim == 2)
        structured = false;
      else if (Utilities::match_at_string_start(entries[i], "ET=BRICK") &&
               dim == 3)
        structured = false;
      else if (Utilities::match_at_string_start(entries[i], "ET="))
        // unsupported ElementType
        {
          AssertThrow(
            false,
            ExcMessage(
              "The tecplot file contains an unsupported ElementType."));
        }
      else if (Utilities::match_at_string_start(entries[i], "I="))
        IJK[0] = Utilities::get_integer_at_position(entries[i], 2).first;
      else if (Utilities::match_at_string_start(entries[i], "J="))
        {
          IJK[1] = Utilities::get_integer_at_position(entries[i], 2).first;
          AssertThrow(
            dim > 1 || IJK[1] == 1,
            ExcMessage(
              "Parameter 'J=' found in tecplot, although this is only possible for dimensions greater than 1."));
        }
      else if (Utilities::match_at_string_start(entries[i], "K="))
        {
          IJK[2] = Utilities::get_integer_at_position(entries[i], 2).first;
          AssertThrow(
            dim > 2 || IJK[2] == 1,
            ExcMessage(
              "Parameter 'K=' found in tecplot, although this is only possible for dimensions greater than 2."));
        }
      else if (Utilities::match_at_string_start(entries[i], "N="))
        n_vertices = Utilities::get_integer_at_position(entries[i], 2).first;
      else if (Utilities::match_at_string_start(entries[i], "E="))
        n_cells = Utilities::get_integer_at_position(entries[i], 2).first;
    }

  // now we have read all the fields we are
  // interested in. do some checks and
  // calculate the variables
  if (structured)
    {
      n_vertices = 1;
      n_cells    = 1;
      for (unsigned int d = 0; d < dim; ++d)
        {
          AssertThrow(
            IJK[d] > 0,
            ExcMessage(
              "Tecplot file does not contain a complete and consistent set of parameters"));
          n_vertices *= IJK[d];
          n_cells *= (IJK[d] - 1);
        }
    }
  else
    {
      AssertThrow(
        n_vertices > 0,
        ExcMessage(
          "Tecplot file does not contain a complete and consistent set of parameters"));
      if (n_cells == 0)
        // this means an error, although
        // tecplot itself accepts entries like
        // 'J=20' instead of 'E=20'. therefore,
        // take the max of IJK
        n_cells = *std::max_element(IJK.begin(), IJK.end());
      AssertThrow(
        n_cells > 0,
        ExcMessage(
          "Tecplot file does not contain a complete and consistent set of parameters"));
    }
}



template <>
void
GridIn<2>::read_tecplot(std::istream &in)
{
  const unsigned int dim      = 2;
  const unsigned int spacedim = 2;
  Assert(tria != nullptr, ExcNoTriangulationSelected());
  AssertThrow(in, ExcIO());

  // skip comments at start of file
  skip_comment_lines(in, '#');

  // some strings for parsing the header
  std::string line, header;

  // first, concatenate all header lines
  // create a searchstring with almost all
  // letters. exclude e and E from the letters
  // to search, as they might appear in
  // exponential notation
  std::string letters = "abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ";

  getline(in, line);
  while (line.find_first_of(letters) != std::string::npos)
    {
      header += " " + line;
      getline(in, line);
    }

  // now create some variables holding
  // important information on the mesh, get
  // this information from the header string
  std::vector<unsigned int> tecplot2deal(dim);
  std::vector<unsigned int> IJK(dim);
  unsigned int              n_vars, n_vertices, n_cells;
  bool                      structured, blocked;

  parse_tecplot_header(header,
                       tecplot2deal,
                       n_vars,
                       n_vertices,
                       n_cells,
                       IJK,
                       structured,
                       blocked);

  // reserve space for vertices. note, that in
  // tecplot vertices are ordered beginning
  // with 1, whereas in deal all indices start
  // with 0. in order not to use -1 for all the
  // connectivity information, a 0th vertex
  // (unused) is inserted at the origin.
  std::vector<Point<spacedim>> vertices(n_vertices + 1);
  vertices[0] = Point<spacedim>();
  // reserve space for cells
  std::vector<CellData<dim>> cells(n_cells);
  SubCellData                subcelldata;

  if (blocked)
    {
      // blocked data format. first we get all
      // the values of the first variable for
      // all points, after that we get all
      // values for the second variable and so
      // on.

      // dummy variable to read in all the info
      // we do not want to use
      double dummy;
      // which is the first index to read in
      // the loop (see below)
      unsigned int next_index = 0;

      // note, that we have already read the
      // first line containing the first variable
      if (tecplot2deal[0] == 0)
        {
          // we need the information in this
          // line, so extract it
          std::vector<std::string> first_var =
            Utilities::break_text_into_lines(line, 1);
          char *endptr;
          for (unsigned int i = 1; i < first_var.size() + 1; ++i)
            vertices[i](0) = std::strtod(first_var[i - 1].c_str(), &endptr);

          // if there are many points, the data
          // for this var might continue in the
          // next line(s)
          for (unsigned int j = first_var.size() + 1; j < n_vertices + 1; ++j)
            in >> vertices[j](next_index);
          // now we got all values of the first
          // variable, so increase the counter
          next_index = 1;
        }

      // main loop over all variables
      for (unsigned int i = 1; i < n_vars; ++i)
        {
          // if we read all the important
          // variables and do not want to
          // read further, because we are
          // using a structured grid, we can
          // stop here (and skip, for
          // example, a whole lot of solution
          // variables)
          if (next_index == dim && structured)
            break;

          if ((next_index < dim) && (i == tecplot2deal[next_index]))
            {
              // we need this line, read it in
              for (unsigned int j = 1; j < n_vertices + 1; ++j)
                in >> vertices[j](next_index);
              ++next_index;
            }
          else
            {
              // we do not need this line, read
              // it in and discard it
              for (unsigned int j = 1; j < n_vertices + 1; ++j)
                in >> dummy;
            }
        }
      Assert(next_index == dim, ExcInternalError());
    }
  else
    {
      // the data is not blocked, so we get all
      // the variables for one point, then the
      // next and so on. create a vector to
      // hold these components
      std::vector<double> vars(n_vars);

      // now fill the first vertex. note, that we
      // have already read the first line
      // containing the first vertex
      std::vector<std::string> first_vertex =
        Utilities::break_text_into_lines(line, 1);
      char *endptr;
      for (unsigned int d = 0; d < dim; ++d)
        vertices[1](d) =
          std::strtod(first_vertex[tecplot2deal[d]].c_str(), &endptr);

      // read the remaining vertices from the
      // list
      for (unsigned int v = 2; v < n_vertices + 1; ++v)
        {
          for (unsigned int i = 0; i < n_vars; ++i)
            in >> vars[i];
          // fill the vertex
          // coordinates. respect the position
          // of coordinates in the list of
          // variables
          for (unsigned int i = 0; i < dim; ++i)
            vertices[v](i) = vars[tecplot2deal[i]];
        }
    }

  if (structured)
    {
      // this is the part of the code that only
      // works in 2d
      unsigned int I = IJK[0], J = IJK[1];

      unsigned int cell = 0;
      // set up array of cells
      for (unsigned int j = 0; j < J - 1; ++j)
        for (unsigned int i = 1; i < I; ++i)
          {
            cells[cell].vertices[0] = i + j * I;
            cells[cell].vertices[1] = i + 1 + j * I;
            cells[cell].vertices[2] = i + (j + 1) * I;
            cells[cell].vertices[3] = i + 1 + (j + 1) * I;
            ++cell;
          }
      Assert(cell == n_cells, ExcInternalError());
      std::vector<unsigned int> boundary_vertices(2 * I + 2 * J - 4);
      unsigned int              k = 0;
      for (unsigned int i = 1; i < I + 1; ++i)
        {
          boundary_vertices[k] = i;
          ++k;
          boundary_vertices[k] = i + (J - 1) * I;
          ++k;
        }
      for (unsigned int j = 1; j < J - 1; ++j)
        {
          boundary_vertices[k] = 1 + j * I;
          ++k;
          boundary_vertices[k] = I + j * I;
          ++k;
        }
      Assert(k == boundary_vertices.size(), ExcInternalError());
      // delete the duplicated vertices at the
      // boundary, which occur, e.g. in c-type
      // or o-type grids around a body
      // (airfoil). this automatically deletes
      // unused vertices as well.
      GridTools::delete_duplicated_vertices(vertices,
                                            cells,
                                            subcelldata,
                                            boundary_vertices);
    }
  else
    {
      // set up array of cells, unstructured
      // mode, so the connectivity is
      // explicitly given
      for (unsigned int i = 0; i < n_cells; ++i)
        {
          // note that since in the input file
          // we found the number of cells at
          // the top, there should still be
          // input here, so check this:
          AssertThrow(in, ExcIO());

          // get the connectivity from the
          // input file. the vertices are
          // ordered like in the ucd format
          for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
            in >> cells[i].vertices[GeometryInfo<dim>::ucd_to_deal[j]];
        }
      // do some clean-up on vertices
      GridTools::delete_unused_vertices(vertices, cells, subcelldata);
    }

  // check that no forbidden arrays are
  // used. as we do not read in any
  // subcelldata, nothing should happen here.
  Assert(subcelldata.check_consistency(dim), ExcInternalError());
  AssertThrow(in, ExcIO());

  // do some cleanup on cells
  GridTools::invert_all_negative_measure_cells(vertices, cells);
  GridTools::consistently_order_cells(cells);
  tria->create_triangulation(vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_tecplot(std::istream &)
{
  Assert(false, ExcNotImplemented());
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read_assimp(const std::string &filename,
                                   const unsigned int mesh_index,
                                   const bool         remove_duplicates,
                                   const double       tol,
                                   const bool         ignore_unsupported_types)
{
#ifdef DEAL_II_WITH_ASSIMP
  // Only good for surface grids.
  AssertThrow(dim < 3, ExcImpossibleInDim(dim));

  // Create an instance of the Importer class
  Assimp::Importer importer;

  // And have it read the given file with some  postprocessing
  const aiScene *scene =
    importer.ReadFile(filename.c_str(),
                      aiProcess_RemoveComponent |
                        aiProcess_JoinIdenticalVertices |
                        aiProcess_ImproveCacheLocality | aiProcess_SortByPType |
                        aiProcess_OptimizeGraph | aiProcess_OptimizeMeshes);

  // If the import failed, report it
  AssertThrow(scene != nullptr, ExcMessage(importer.GetErrorString()));

  AssertThrow(scene->mNumMeshes != 0,
              ExcMessage("Input file contains no meshes."));

  AssertThrow((mesh_index == numbers::invalid_unsigned_int) ||
                (mesh_index < scene->mNumMeshes),
              ExcMessage("Too few meshes in the file."));

  unsigned int start_mesh =
    (mesh_index == numbers::invalid_unsigned_int ? 0 : mesh_index);
  unsigned int end_mesh =
    (mesh_index == numbers::invalid_unsigned_int ? scene->mNumMeshes :
                                                   mesh_index + 1);

  // Deal.II objects are created empty, and then filled with imported file.
  std::vector<Point<spacedim>> vertices;
  std::vector<CellData<dim>>   cells;
  SubCellData                  subcelldata;

  // A series of counters to merge cells.
  unsigned int v_offset = 0;
  unsigned int c_offset = 0;

  // The index of the mesh will be used as a material index.
  for (unsigned int m = start_mesh; m < end_mesh; ++m)
    {
      const aiMesh *mesh = scene->mMeshes[m];

      // Check that we know what to do with this mesh, otherwise just
      // ignore it
      if ((dim == 2) && mesh->mPrimitiveTypes != aiPrimitiveType_POLYGON)
        {
          AssertThrow(ignore_unsupported_types,
                      ExcMessage("Incompatible mesh " + std::to_string(m) +
                                 "/" + std::to_string(scene->mNumMeshes)));
          continue;
        }
      else if ((dim == 1) && mesh->mPrimitiveTypes != aiPrimitiveType_LINE)
        {
          AssertThrow(ignore_unsupported_types,
                      ExcMessage("Incompatible mesh " + std::to_string(m) +
                                 "/" + std::to_string(scene->mNumMeshes)));
          continue;
        }
      // Vertices
      const unsigned int n_vertices = mesh->mNumVertices;
      const aiVector3D * mVertices  = mesh->mVertices;

      // Faces
      const unsigned int n_faces = mesh->mNumFaces;
      const aiFace *     mFaces  = mesh->mFaces;

      vertices.resize(v_offset + n_vertices);
      cells.resize(c_offset + n_faces);

      for (unsigned int i = 0; i < n_vertices; ++i)
        for (unsigned int d = 0; d < spacedim; ++d)
          vertices[i + v_offset][d] = mVertices[i][d];

      unsigned int valid_cell = c_offset;
      for (unsigned int i = 0; i < n_faces; ++i)
        {
          if (mFaces[i].mNumIndices == GeometryInfo<dim>::vertices_per_cell)
            {
              for (const unsigned int f : GeometryInfo<dim>::vertex_indices())
                {
                  cells[valid_cell]
                    .vertices[GeometryInfo<dim>::ucd_to_deal[f]] =
                    mFaces[i].mIndices[f] + v_offset;
                }
              cells[valid_cell].material_id = m;
              ++valid_cell;
            }
          else
            {
              AssertThrow(ignore_unsupported_types,
                          ExcMessage("Face " + std::to_string(i) + " of mesh " +
                                     std::to_string(m) + " has " +
                                     std::to_string(mFaces[i].mNumIndices) +
                                     " vertices. We expected only " +
                                     std::to_string(
                                       GeometryInfo<dim>::vertices_per_cell)));
            }
        }
      cells.resize(valid_cell);

      // The vertices are added all at once. Cells are checked for
      // validity, so only valid_cells are now present in the deal.II
      // list of cells.
      v_offset += n_vertices;
      c_offset = valid_cell;
    }

  // No cells were read
  if (cells.size() == 0)
    return;

  if (remove_duplicates)
    {
      // The function delete_duplicated_vertices() needs to be called more
      // than once if a vertex is duplicated more than once. So we keep
      // calling it until the number of vertices does not change any more.
      unsigned int n_verts = 0;
      while (n_verts != vertices.size())
        {
          n_verts = vertices.size();
          std::vector<unsigned int> considered_vertices;
          GridTools::delete_duplicated_vertices(
            vertices, cells, subcelldata, considered_vertices, tol);
        }
    }

  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  if (dim == spacedim)
    GridTools::invert_all_negative_measure_cells(vertices, cells);
  GridTools::consistently_order_cells(cells);
  tria->create_triangulation(vertices, cells, subcelldata);

#else
  (void)filename;
  (void)mesh_index;
  (void)remove_duplicates;
  (void)tol;
  (void)ignore_unsupported_types;
  AssertThrow(false, ExcNeedsAssimp());
#endif
}

#ifdef DEAL_II_TRILINOS_WITH_SEACAS
// Namespace containing some extra functions for reading ExodusII files
namespace
{
  // Convert ExodusII strings to cell types. Use the number of nodes per element
  // to disambiguate some cases.
  ReferenceCell
  exodusii_name_to_type(const std::string &type_name,
                        const int          n_nodes_per_element)
  {
    Assert(type_name.size() > 0, ExcInternalError());
    // Try to canonify the name by switching to upper case and removing trailing
    // numbers. This makes, e.g., pyramid, PYRAMID, PYRAMID5, and PYRAMID13 all
    // equal.
    std::string type_name_2 = type_name;
    std::transform(type_name_2.begin(),
                   type_name_2.end(),
                   type_name_2.begin(),
                   [](unsigned char c) { return std::toupper(c); });
    const std::string numbers = "0123456789";
    type_name_2.erase(std::find_first_of(type_name_2.begin(),
                                         type_name_2.end(),
                                         numbers.begin(),
                                         numbers.end()),
                      type_name_2.end());

    if (type_name_2 == "TRI" || type_name_2 == "TRIANGLE")
      return ReferenceCells::Triangle;
    else if (type_name_2 == "QUAD" || type_name_2 == "QUADRILATERAL")
      return ReferenceCells::Quadrilateral;
    else if (type_name_2 == "SHELL")
      {
        if (n_nodes_per_element == 3)
          return ReferenceCells::Triangle;
        else
          return ReferenceCells::Quadrilateral;
      }
    else if (type_name_2 == "TET" || type_name_2 == "TETRA" ||
             type_name_2 == "TETRAHEDRON")
      return ReferenceCells::Tetrahedron;
    else if (type_name_2 == "PYRA" || type_name_2 == "PYRAMID")
      return ReferenceCells::Pyramid;
    else if (type_name_2 == "WEDGE")
      return ReferenceCells::Wedge;
    else if (type_name_2 == "HEX" || type_name_2 == "HEXAHEDRON")
      return ReferenceCells::Hexahedron;

    Assert(false, ExcNotImplemented());
    return ReferenceCells::Invalid;
  }

  // Associate deal.II boundary ids with sidesets (a face can be in multiple
  // sidesets - to translate we assign each set of side set ids to a
  // boundary_id or manifold_id)
  template <int dim, int spacedim = dim>
  std::pair<SubCellData, std::vector<std::vector<int>>>
  read_exodusii_sidesets(const int                         ex_id,
                         const int                         n_side_sets,
                         const std::vector<CellData<dim>> &cells,
                         const bool apply_all_indicators_to_manifolds)
  {
    SubCellData                   subcelldata;
    std::vector<std::vector<int>> b_or_m_id_to_sideset_ids;
    // boundary id 0 is the default
    b_or_m_id_to_sideset_ids.emplace_back();
    // deal.II does not support assigning boundary ids with nonzero codimension
    // meshes so completely skip this information in that case.
    //
    // Exodus prints warnings if we try to get empty sets so always check first
    if (dim == spacedim && n_side_sets > 0)
      {
        std::vector<int> side_set_ids(n_side_sets);
        int ierr = ex_get_ids(ex_id, EX_SIDE_SET, side_set_ids.data());
        AssertThrowExodusII(ierr);

        // First collect all side sets on all boundary faces (indexed here as
        // max_faces_per_cell * cell_n + face_n). We then sort and uniquify the
        // side sets so that we can convert a set of side set indices into a
        // single deal.II boundary or manifold id (and save the correspondence).
        constexpr auto max_faces_per_cell = GeometryInfo<dim>::faces_per_cell;
        std::map<std::size_t, std::vector<int>> face_side_sets;
        for (const int side_set_id : side_set_ids)
          {
            int n_sides                = -1;
            int n_distribution_factors = -1;

            ierr = ex_get_set_param(ex_id,
                                    EX_SIDE_SET,
                                    side_set_id,
                                    &n_sides,
                                    &n_distribution_factors);
            AssertThrowExodusII(ierr);
            if (n_sides > 0)
              {
                std::vector<int> elements(n_sides);
                std::vector<int> faces(n_sides);
                ierr = ex_get_set(ex_id,
                                  EX_SIDE_SET,
                                  side_set_id,
                                  elements.data(),
                                  faces.data());
                AssertThrowExodusII(ierr);

                // According to the manual (subsection 4.8): "The internal
                // number of an element numbering is defined implicitly by the
                // order in which it appears in the file. Elements are numbered
                // internally (beginning with 1) consecutively across all
                // element blocks." Hence element i in Exodus numbering is entry
                // i - 1 in the cells array.
                for (int side_n = 0; side_n < n_sides; ++side_n)
                  {
                    const long        element_n = elements[side_n] - 1;
                    const long        face_n    = faces[side_n] - 1;
                    const std::size_t face_id =
                      element_n * max_faces_per_cell + face_n;
                    face_side_sets[face_id].push_back(side_set_id);
                  }
              }
          }

        // Collect into a sortable data structure:
        std::vector<std::pair<std::size_t, std::vector<int>>>
          face_id_to_side_sets;
        for (auto &pair : face_side_sets)
          {
            Assert(pair.second.size() > 0, ExcInternalError());
            face_id_to_side_sets.push_back(std::move(pair));
          }

        // sort by side sets:
        std::sort(face_id_to_side_sets.begin(),
                  face_id_to_side_sets.end(),
                  [](const auto &a, const auto &b) {
                    return std::lexicographical_compare(a.second.begin(),
                                                        a.second.end(),
                                                        b.second.begin(),
                                                        b.second.end());
                  });

        types::boundary_id current_b_or_m_id = 0;
        for (const auto &pair : face_id_to_side_sets)
          {
            const std::size_t       face_id          = pair.first;
            const std::vector<int> &face_sideset_ids = pair.second;
            if (face_sideset_ids != b_or_m_id_to_sideset_ids.back())
              {
                // Since we sorted by sideset ids we are guaranteed that if this
                // doesn't match the last set then it has not yet been seen
                ++current_b_or_m_id;
                b_or_m_id_to_sideset_ids.push_back(face_sideset_ids);
                Assert(current_b_or_m_id == b_or_m_id_to_sideset_ids.size() - 1,
                       ExcInternalError());
              }
            // Record the b_or_m_id of the current face.
            const unsigned int   local_face_n = face_id % max_faces_per_cell;
            const CellData<dim> &cell = cells[face_id / max_faces_per_cell];
            const ReferenceCell  cell_type =
              ReferenceCell::n_vertices_to_type(dim, cell.vertices.size());
            const unsigned int deal_face_n =
              cell_type.exodusii_face_to_deal_face(local_face_n);
            const ReferenceCell face_reference_cell =
              cell_type.face_reference_cell(deal_face_n);

            // The orientation we pick doesn't matter here since when we create
            // the Triangulation we will sort the vertices for each CellData
            // object created here.
            if (dim == 2)
              {
                CellData<1> boundary_line(face_reference_cell.n_vertices());
                if (apply_all_indicators_to_manifolds)
                  boundary_line.manifold_id = current_b_or_m_id;
                else
                  boundary_line.boundary_id = current_b_or_m_id;
                for (unsigned int j = 0; j < face_reference_cell.n_vertices();
                     ++j)
                  boundary_line.vertices[j] =
                    cell.vertices[cell_type.face_to_cell_vertices(
                      deal_face_n, j, 0)];

                subcelldata.boundary_lines.push_back(std::move(boundary_line));
              }
            else if (dim == 3)
              {
                CellData<2> boundary_quad(face_reference_cell.n_vertices());
                if (apply_all_indicators_to_manifolds)
                  boundary_quad.manifold_id = current_b_or_m_id;
                else
                  boundary_quad.boundary_id = current_b_or_m_id;
                for (unsigned int j = 0; j < face_reference_cell.n_vertices();
                     ++j)
                  boundary_quad.vertices[j] =
                    cell.vertices[cell_type.face_to_cell_vertices(
                      deal_face_n, j, 0)];

                subcelldata.boundary_quads.push_back(std::move(boundary_quad));
              }
          }
      }

    return std::make_pair(std::move(subcelldata),
                          std::move(b_or_m_id_to_sideset_ids));
  }
} // namespace
#endif

template <int dim, int spacedim>
typename GridIn<dim, spacedim>::ExodusIIData
GridIn<dim, spacedim>::read_exodusii(
  const std::string &filename,
  const bool         apply_all_indicators_to_manifolds)
{
#ifdef DEAL_II_TRILINOS_WITH_SEACAS
  // deal.II always uses double precision numbers for geometry
  int component_word_size = sizeof(double);
  // setting to zero uses the stored word size
  int   floating_point_word_size = 0;
  float ex_version               = 0.0;

  const int ex_id = ex_open(filename.c_str(),
                            EX_READ,
                            &component_word_size,
                            &floating_point_word_size,
                            &ex_version);
  AssertThrow(ex_id > 0,
              ExcMessage("ExodusII failed to open the specified input file."));

  // Read basic mesh information:
  std::vector<char> string_temp(MAX_LINE_LENGTH + 1, '\0');
  int               mesh_dimension   = 0;
  int               n_nodes          = 0;
  int               n_elements       = 0;
  int               n_element_blocks = 0;
  int               n_node_sets      = 0;
  int               n_side_sets      = 0;

  int ierr = ex_get_init(ex_id,
                         string_temp.data(),
                         &mesh_dimension,
                         &n_nodes,
                         &n_elements,
                         &n_element_blocks,
                         &n_node_sets,
                         &n_side_sets);
  AssertThrowExodusII(ierr);
  AssertDimension(mesh_dimension, spacedim);

  // Read nodes:
  std::vector<double> xs(n_nodes);
  std::vector<double> ys(n_nodes);
  std::vector<double> zs(n_nodes);

  ierr = ex_get_coord(ex_id, xs.data(), ys.data(), zs.data());
  AssertThrowExodusII(ierr);

  // Even if there is a node numbering array the values stored inside the
  // ExodusII file must use the contiguous, internal ordering (see Section 4.5
  // of the manual - "Internal (contiguously numbered) node and element IDs must
  // be used for all data structures that contain node or element numbers (IDs),
  // including node set node lists, side set element lists, and element
  // connectivity.")
  std::vector<Point<spacedim>> vertices;
  vertices.reserve(n_nodes);
  for (int vertex_n = 0; vertex_n < n_nodes; ++vertex_n)
    {
      switch (spacedim)
        {
          case 1:
            vertices.emplace_back(xs[vertex_n]);
            break;
          case 2:
            vertices.emplace_back(xs[vertex_n], ys[vertex_n]);
            break;
          case 3:
            vertices.emplace_back(xs[vertex_n], ys[vertex_n], zs[vertex_n]);
            break;
          default:
            Assert(spacedim <= 3, ExcNotImplemented());
        }
    }

  std::vector<int> element_block_ids(n_element_blocks);
  ierr = ex_get_ids(ex_id, EX_ELEM_BLOCK, element_block_ids.data());
  AssertThrowExodusII(ierr);

  std::vector<CellData<dim>> cells;
  // Elements are grouped together by same reference cell type in element
  // blocks. There may be multiple blocks for a single reference cell type,
  // but "each element block may contain only one element type".
  for (const int element_block_id : element_block_ids)
    {
      std::fill(string_temp.begin(), string_temp.end(), '\0');
      int n_block_elements         = 0;
      int n_nodes_per_element      = 0;
      int n_edges_per_element      = 0;
      int n_faces_per_element      = 0;
      int n_attributes_per_element = 0;

      // Extract element data.
      ierr = ex_get_block(ex_id,
                          EX_ELEM_BLOCK,
                          element_block_id,
                          string_temp.data(),
                          &n_block_elements,
                          &n_nodes_per_element,
                          &n_edges_per_element,
                          &n_faces_per_element,
                          &n_attributes_per_element);
      AssertThrowExodusII(ierr);
      const ReferenceCell type =
        exodusii_name_to_type(string_temp.data(), n_nodes_per_element);

      // The number of nodes per element may be larger than what we want to
      // read - for example, if the Exodus file contains a QUAD9 element, we
      // only want to read the first four values and ignore the rest.
      Assert(int(type.n_vertices()) <= n_nodes_per_element, ExcInternalError());

      std::vector<int> connection(n_nodes_per_element * n_block_elements);
      ierr = ex_get_conn(ex_id,
                         EX_ELEM_BLOCK,
                         element_block_id,
                         connection.data(),
                         nullptr,
                         nullptr);
      AssertThrowExodusII(ierr);

      for (unsigned int elem_n = 0; elem_n < connection.size();
           elem_n += n_nodes_per_element)
        {
          CellData<dim> cell(type.n_vertices());
          for (unsigned int i : type.vertex_indices())
            {
              cell.vertices[type.exodusii_vertex_to_deal_vertex(i)] =
                connection[elem_n + i] - 1;
            }
          cell.material_id = element_block_id;
          cells.push_back(cell);
        }
    }

  // Extract boundary data.
  auto pair = read_exodusii_sidesets<dim, spacedim>(
    ex_id, n_side_sets, cells, apply_all_indicators_to_manifolds);
  ierr = ex_close(ex_id);
  AssertThrowExodusII(ierr);

  tria->create_triangulation(vertices, cells, pair.first);
  ExodusIIData out;
  out.id_to_sideset_ids = std::move(pair.second);
  return out;
#else
  (void)filename;
  (void)apply_all_indicators_to_manifolds;
  AssertThrow(false, ExcNeedsExodusII());
  return {};
#endif
}


template <int dim, int spacedim>
void
GridIn<dim, spacedim>::skip_empty_lines(std::istream &in)
{
  std::string line;
  while (in)
    {
      // get line
      getline(in, line);

      // check if this is a line that
      // consists only of spaces, and
      // if not put the whole thing
      // back and return
      if (std::find_if(line.begin(), line.end(), [](const char c) {
            return c != ' ';
          }) != line.end())
        {
          in.putback('\n');
          for (int i = line.length() - 1; i >= 0; --i)
            in.putback(line[i]);
          return;
        }

      // else: go on with next line
    }
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::skip_comment_lines(std::istream &in,
                                          const char    comment_start)
{
  char c;
  // loop over the following comment
  // lines
  while (in.get(c) && c == comment_start)
    // loop over the characters after
    // the comment starter
    while (in.get() != '\n')
      ;


  // put back first character of
  // first non-comment line
  if (in)
    in.putback(c);

  // at last: skip additional empty lines, if present
  skip_empty_lines(in);
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::debug_output_grid(
  const std::vector<CellData<dim>> & /*cells*/,
  const std::vector<Point<spacedim>> & /*vertices*/,
  std::ostream & /*out*/)
{
  Assert(false, ExcNotImplemented());
}



template <>
void
GridIn<2>::debug_output_grid(const std::vector<CellData<2>> &cells,
                             const std::vector<Point<2>> &   vertices,
                             std::ostream &                  out)
{
  double min_x = vertices[cells[0].vertices[0]](0),
         max_x = vertices[cells[0].vertices[0]](0),
         min_y = vertices[cells[0].vertices[0]](1),
         max_y = vertices[cells[0].vertices[0]](1);

  for (unsigned int i = 0; i < cells.size(); ++i)
    {
      for (const auto vertex : cells[i].vertices)
        {
          const Point<2> &p = vertices[vertex];

          if (p(0) < min_x)
            min_x = p(0);
          if (p(0) > max_x)
            max_x = p(0);
          if (p(1) < min_y)
            min_y = p(1);
          if (p(1) > max_y)
            max_y = p(1);
        }

      out << "# cell " << i << std::endl;
      Point<2> center;
      for (const auto vertex : cells[i].vertices)
        center += vertices[vertex];
      center /= 4;

      out << "set label \"" << i << "\" at " << center(0) << ',' << center(1)
          << " center" << std::endl;

      // first two line right direction
      for (unsigned int f = 0; f < 2; ++f)
        out << "set arrow from " << vertices[cells[i].vertices[f]](0) << ','
            << vertices[cells[i].vertices[f]](1) << " to "
            << vertices[cells[i].vertices[(f + 1) % 4]](0) << ','
            << vertices[cells[i].vertices[(f + 1) % 4]](1) << std::endl;
      // other two lines reverse direction
      for (unsigned int f = 2; f < 4; ++f)
        out << "set arrow from " << vertices[cells[i].vertices[(f + 1) % 4]](0)
            << ',' << vertices[cells[i].vertices[(f + 1) % 4]](1) << " to "
            << vertices[cells[i].vertices[f]](0) << ','
            << vertices[cells[i].vertices[f]](1) << std::endl;
      out << std::endl;
    }


  out << std::endl
      << "set nokey" << std::endl
      << "pl [" << min_x << ':' << max_x << "][" << min_y << ':' << max_y
      << "] " << min_y << std::endl
      << "pause -1" << std::endl;
}



template <>
void
GridIn<3>::debug_output_grid(const std::vector<CellData<3>> &cells,
                             const std::vector<Point<3>> &   vertices,
                             std::ostream &                  out)
{
  for (const auto &cell : cells)
    {
      // line 0
      out << vertices[cell.vertices[0]] << std::endl
          << vertices[cell.vertices[1]] << std::endl
          << std::endl
          << std::endl;
      // line 1
      out << vertices[cell.vertices[1]] << std::endl
          << vertices[cell.vertices[2]] << std::endl
          << std::endl
          << std::endl;
      // line 2
      out << vertices[cell.vertices[3]] << std::endl
          << vertices[cell.vertices[2]] << std::endl
          << std::endl
          << std::endl;
      // line 3
      out << vertices[cell.vertices[0]] << std::endl
          << vertices[cell.vertices[3]] << std::endl
          << std::endl
          << std::endl;
      // line 4
      out << vertices[cell.vertices[4]] << std::endl
          << vertices[cell.vertices[5]] << std::endl
          << std::endl
          << std::endl;
      // line 5
      out << vertices[cell.vertices[5]] << std::endl
          << vertices[cell.vertices[6]] << std::endl
          << std::endl
          << std::endl;
      // line 6
      out << vertices[cell.vertices[7]] << std::endl
          << vertices[cell.vertices[6]] << std::endl
          << std::endl
          << std::endl;
      // line 7
      out << vertices[cell.vertices[4]] << std::endl
          << vertices[cell.vertices[7]] << std::endl
          << std::endl
          << std::endl;
      // line 8
      out << vertices[cell.vertices[0]] << std::endl
          << vertices[cell.vertices[4]] << std::endl
          << std::endl
          << std::endl;
      // line 9
      out << vertices[cell.vertices[1]] << std::endl
          << vertices[cell.vertices[5]] << std::endl
          << std::endl
          << std::endl;
      // line 10
      out << vertices[cell.vertices[2]] << std::endl
          << vertices[cell.vertices[6]] << std::endl
          << std::endl
          << std::endl;
      // line 11
      out << vertices[cell.vertices[3]] << std::endl
          << vertices[cell.vertices[7]] << std::endl
          << std::endl
          << std::endl;
    }
}



template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read(const std::string &filename, Format format)
{
  // Search file class for meshes
  PathSearch  search("MESH");
  std::string name;
  // Open the file and remember its name
  if (format == Default)
    name = search.find(filename);
  else
    name = search.find(filename, default_suffix(format));


  if (format == Default)
    {
      const std::string::size_type slashpos = name.find_last_of('/');
      const std::string::size_type dotpos   = name.find_last_of('.');
      if (dotpos < name.length() &&
          (dotpos > slashpos || slashpos == std::string::npos))
        {
          std::string ext = name.substr(dotpos + 1);
          format          = parse_format(ext);
        }
    }

  if (format == assimp)
    {
      read_assimp(name);
    }
  else if (format == exodusii)
    {
      read_exodusii(name);
    }
  else
    {
      std::ifstream in(name.c_str());
      read(in, format);
    }
}


template <int dim, int spacedim>
void
GridIn<dim, spacedim>::read(std::istream &in, Format format)
{
  if (format == Default)
    format = default_format;

  switch (format)
    {
      case dbmesh:
        read_dbmesh(in);
        return;

      case msh:
        read_msh(in);
        return;

      case vtk:
        read_vtk(in);
        return;

      case vtu:
        read_vtu(in);
        return;

      case unv:
        read_unv(in);
        return;

      case ucd:
        read_ucd(in);
        return;

      case abaqus:
        read_abaqus(in);
        return;

      case xda:
        read_xda(in);
        return;

      case tecplot:
        read_tecplot(in);
        return;

      case assimp:
        Assert(false,
               ExcMessage("There is no read_assimp(istream &) function. "
                          "Use the read_assimp(string &filename, ...) "
                          "functions, instead."));
        return;

      case exodusii:
        Assert(false,
               ExcMessage("There is no read_exodusii(istream &) function. "
                          "Use the read_exodusii(string &filename, ...) "
                          "function, instead."));
        return;

      case Default:
        break;
    }
  Assert(false, ExcInternalError());
}



template <int dim, int spacedim>
std::string
GridIn<dim, spacedim>::default_suffix(const Format format)
{
  switch (format)
    {
      case dbmesh:
        return ".dbmesh";
      case exodusii:
        return ".e";
      case msh:
        return ".msh";
      case vtk:
        return ".vtk";
      case vtu:
        return ".vtu";
      case unv:
        return ".unv";
      case ucd:
        return ".inp";
      case abaqus:
        return ".inp"; // Typical suffix for Abaqus mesh files conflicts with
                       // UCD.
      case xda:
        return ".xda";
      case tecplot:
        return ".dat";
      default:
        Assert(false, ExcNotImplemented());
        return ".unknown_format";
    }
}



template <int dim, int spacedim>
typename GridIn<dim, spacedim>::Format
GridIn<dim, spacedim>::parse_format(const std::string &format_name)
{
  if (format_name == "dbmesh")
    return dbmesh;

  if (format_name == "exodusii")
    return exodusii;

  if (format_name == "msh")
    return msh;

  if (format_name == "unv")
    return unv;

  if (format_name == "vtk")
    return vtk;

  if (format_name == "vtu")
    return vtu;

  // This is also the typical extension of Abaqus input files.
  if (format_name == "inp")
    return ucd;

  if (format_name == "ucd")
    return ucd;

  if (format_name == "xda")
    return xda;

  if (format_name == "tecplot")
    return tecplot;

  if (format_name == "dat")
    return tecplot;

  if (format_name == "plt")
    // Actually, this is the extension for the
    // tecplot binary format, which we do not
    // support right now. However, some people
    // tend to create tecplot ascii files with
    // the extension 'plt' instead of
    // 'dat'. Thus, include this extension
    // here. If it actually is a binary file,
    // the read_tecplot() function will fail
    // and throw an exception, anyway.
    return tecplot;

  AssertThrow(false, ExcInvalidState());
  // return something weird
  return Format(Default);
}



template <int dim, int spacedim>
std::string
GridIn<dim, spacedim>::get_format_names()
{
  return "dbmesh|exodusii|msh|unv|vtk|vtu|ucd|abaqus|xda|tecplot|assimp";
}



namespace
{
  template <int dim, int spacedim>
  Abaqus_to_UCD<dim, spacedim>::Abaqus_to_UCD()
    : tolerance(5e-16) // Used to offset Cubit tolerance error when outputting
                       // value close to zero
  {
    AssertThrow(spacedim == 2 || spacedim == 3, ExcNotImplemented());
  }



  // Convert from a string to some other data type
  // Reference: http://www.codeguru.com/forum/showthread.php?t=231054
  template <class T>
  bool
  from_string(T &t, const std::string &s, std::ios_base &(*f)(std::ios_base &))
  {
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
  }



  // Extract an integer from a string
  int
  extract_int(const std::string &s)
  {
    std::string tmp;
    for (const char c : s)
      {
        if (isdigit(c))
          {
            tmp += c;
          }
      }

    int number = 0;
    from_string(number, tmp, std::dec);
    return number;
  }



  template <int dim, int spacedim>
  void
  Abaqus_to_UCD<dim, spacedim>::read_in_abaqus(std::istream &input_stream)
  {
    // References:
    // http://www.egr.msu.edu/software/abaqus/Documentation/docs/v6.7/books/usb/default.htm?startat=pt01ch02.html
    // http://www.cprogramming.com/tutorial/string.html

    AssertThrow(input_stream, ExcIO());
    std::string line;

    while (std::getline(input_stream, line))
      {
      cont:
        std::transform(line.begin(), line.end(), line.begin(), ::toupper);

        if (line.compare("*HEADING") == 0 || line.compare(0, 2, "**") == 0 ||
            line.compare(0, 5, "*PART") == 0)
          {
            // Skip header and comments
            while (std::getline(input_stream, line))
              {
                if (line[0] == '*')
                  goto cont; // My eyes, they burn!
              }
          }
        else if (line.compare(0, 5, "*NODE") == 0)
          {
            // Extract list of vertices
            // Header line might be:
            // *NODE, NSET=ALLNODES
            // *NODE

            // Contains lines in the form:
            // Index, x, y, z
            while (std::getline(input_stream, line))
              {
                if (line[0] == '*')
                  goto cont;

                std::vector<double> node(spacedim + 1);

                std::istringstream iss(line);
                char               comma;
                for (unsigned int i = 0; i < spacedim + 1; ++i)
                  iss >> node[i] >> comma;

                node_list.push_back(node);
              }
          }
        else if (line.compare(0, 8, "*ELEMENT") == 0)
          {
            // Element construction.
            // There are different header formats, the details
            // of which we're not particularly interested in except
            // whether they represent quads or hexahedrals.
            // *ELEMENT, TYPE=S4R, ELSET=EB<material id>
            // *ELEMENT, TYPE=C3D8R, ELSET=EB<material id>
            // *ELEMENT, TYPE=C3D8
            // Elements itself (n=4 or n=8):
            // Index, i[0], ..., i[n]

            int material = 0;
            // Scan for material id
            {
              const std::string before_material = "ELSET=EB";
              const std::size_t idx             = line.find(before_material);
              if (idx != std::string::npos)
                {
                  from_string(material,
                              line.substr(idx + before_material.size()),
                              std::dec);
                }
            }

            // Read ELEMENT definition
            while (std::getline(input_stream, line))
              {
                if (line[0] == '*')
                  goto cont;

                std::istringstream iss(line);
                char               comma;

                // We will store the material id in the zeroth entry of the
                // vector and the rest of the elements represent the global
                // node numbers
                const unsigned int n_data_per_cell =
                  1 + GeometryInfo<dim>::vertices_per_cell;
                std::vector<double> cell(n_data_per_cell);
                for (unsigned int i = 0; i < n_data_per_cell; ++i)
                  iss >> cell[i] >> comma;

                // Overwrite cell index from file by material
                cell[0] = static_cast<double>(material);
                cell_list.push_back(cell);
              }
          }
        else if (line.compare(0, 8, "*SURFACE") == 0)
          {
            // Extract the definitions of boundary surfaces
            // Old format from Cubit:
            // *SURFACE, NAME=SS<boundary indicator>
            //    <element index>,     S<face number>
            // Abaqus default format:
            // *SURFACE, TYPE=ELEMENT, NAME=SURF-<indicator>

            // Get name of the surface and extract id from it;
            // this will be the boundary indicator
            const std::string name_key = "NAME=";
            const std::size_t name_idx_start =
              line.find(name_key) + name_key.size();
            std::size_t name_idx_end = line.find(',', name_idx_start);
            if (name_idx_end == std::string::npos)
              {
                name_idx_end = line.size();
              }
            const int b_indicator = extract_int(
              line.substr(name_idx_start, name_idx_end - name_idx_start));

            // Read SURFACE definition
            // Note that the orientation of the faces is embedded within the
            // definition of each "set" of faces that comprise the surface
            // These are either marked by an "S" or "E" in 3d or 2d
            // respectively.
            while (std::getline(input_stream, line))
              {
                if (line[0] == '*')
                  goto cont;

                // Change all characters to upper case
                std::transform(line.begin(),
                               line.end(),
                               line.begin(),
                               ::toupper);

                // Surface can be created from ELSET, or directly from cells
                // If elsets_list contains a key with specific name - refers to
                // that ELSET, otherwise refers to cell
                std::istringstream iss(line);
                int                el_idx;
                int                face_number;
                char               temp;

                // Get relevant faces, taking into account the element
                // orientation
                std::vector<double> quad_node_list;
                const std::string   elset_name = line.substr(0, line.find(','));
                if (elsets_list.count(elset_name) != 0)
                  {
                    // Surface refers to ELSET
                    std::string stmp;
                    iss >> stmp >> temp >> face_number;

                    const std::vector<int> cells = elsets_list[elset_name];
                    for (const int cell : cells)
                      {
                        el_idx = cell;
                        quad_node_list =
                          get_global_node_numbers(el_idx, face_number);
                        quad_node_list.insert(quad_node_list.begin(),
                                              b_indicator);

                        face_list.push_back(quad_node_list);
                      }
                  }
                else
                  {
                    // Surface refers directly to elements
                    char comma;
                    iss >> el_idx >> comma >> temp >> face_number;
                    quad_node_list =
                      get_global_node_numbers(el_idx, face_number);
                    quad_node_list.insert(quad_node_list.begin(), b_indicator);

                    face_list.push_back(quad_node_list);
                  }
              }
          }
        else if (line.compare(0, 6, "*ELSET") == 0)
          {
            // Get ELSET name.
            // Materials are attached to elsets with specific name
            std::string elset_name;
            {
              const std::string elset_key = "*ELSET, ELSET=";
              const std::size_t idx       = line.find(elset_key);
              if (idx != std::string::npos)
                {
                  const std::string comma       = ",";
                  const std::size_t first_comma = line.find(comma);
                  const std::size_t second_comma =
                    line.find(comma, first_comma + 1);
                  const std::size_t elset_name_start =
                    line.find(elset_key) + elset_key.size();
                  elset_name = line.substr(elset_name_start,
                                           second_comma - elset_name_start);
                }
            }

            // There are two possibilities of storing cells numbers in ELSET:
            // 1. If the header contains the 'GENERATE' keyword, then the next
            // line describes range of cells as:
            //    cell_id_start, cell_id_end, cell_step
            // 2. If the header does not contain the 'GENERATE' keyword, then
            // the next lines contain cells numbers
            std::vector<int>  elements;
            const std::size_t generate_idx = line.find("GENERATE");
            if (generate_idx != std::string::npos)
              {
                // Option (1)
                std::getline(input_stream, line);
                std::istringstream iss(line);
                char               comma;
                int                elid_start;
                int                elid_end;
                int elis_step = 1; // Default if case stride not provided

                // Some files don't have the stride size
                // Compare mesh test cases ./grids/abaqus/3d/other_simple.inp to
                // ./grids/abaqus/2d/2d_test_abaqus.inp
                iss >> elid_start >> comma >> elid_end;
                AssertThrow(comma == ',',
                            ExcMessage(
                              std::string(
                                "While reading an ABAQUS file, the reader "
                                "expected a comma but found a <") +
                              comma + "> in the line <" + line + ">."));
                AssertThrow(
                  elid_start <= elid_end,
                  ExcMessage(
                    std::string(
                      "While reading an ABAQUS file, the reader encountered "
                      "a GENERATE statement in which the upper bound <") +
                    Utilities::int_to_string(elid_end) +
                    "> for the element numbers is not larger or equal "
                    "than the lower bound <" +
                    Utilities::int_to_string(elid_start) + ">."));

                // https://stackoverflow.com/questions/8046357/how-do-i-check-if-a-stringstream-variable-is-empty-null
                if (iss.rdbuf()->in_avail() != 0)
                  iss >> comma >> elis_step;
                AssertThrow(comma == ',',
                            ExcMessage(
                              std::string(
                                "While reading an ABAQUS file, the reader "
                                "expected a comma but found a <") +
                              comma + "> in the line <" + line + ">."));

                for (int i = elid_start; i <= elid_end; i += elis_step)
                  elements.push_back(i);
                elsets_list[elset_name] = elements;

                std::getline(input_stream, line);
              }
            else
              {
                // Option (2)
                while (std::getline(input_stream, line))
                  {
                    if (line[0] == '*')
                      break;

                    std::istringstream iss(line);
                    char               comma;
                    int                elid;
                    while (!iss.eof())
                      {
                        iss >> elid >> comma;
                        AssertThrow(
                          comma == ',',
                          ExcMessage(
                            std::string(
                              "While reading an ABAQUS file, the reader "
                              "expected a comma but found a <") +
                            comma + "> in the line <" + line + ">."));

                        elements.push_back(elid);
                      }
                  }

                elsets_list[elset_name] = elements;
              }

            goto cont;
          }
        else if (line.compare(0, 5, "*NSET") == 0)
          {
            // Skip nodesets; we have no use for them
            while (std::getline(input_stream, line))
              {
                if (line[0] == '*')
                  goto cont;
              }
          }
        else if (line.compare(0, 14, "*SOLID SECTION") == 0)
          {
            // The ELSET name, which describes a section for particular material
            const std::string elset_key = "ELSET=";
            const std::size_t elset_start =
              line.find("ELSET=") + elset_key.size();
            const std::size_t elset_end = line.find(',', elset_start + 1);
            const std::string elset_name =
              line.substr(elset_start, elset_end - elset_start);

            // Solid material definition.
            // We assume that material id is taken from material name,
            // eg. "Material-1" -> ID=1
            const std::string material_key = "MATERIAL=";
            const std::size_t last_equal =
              line.find("MATERIAL=") + material_key.size();
            const std::size_t material_id_start = line.find('-', last_equal);
            int               material_id       = 0;
            from_string(material_id,
                        line.substr(material_id_start + 1),
                        std::dec);

            // Assign material id to cells
            const std::vector<int> &elset_cells = elsets_list[elset_name];
            for (const int elset_cell : elset_cells)
              {
                const int cell_id     = elset_cell - 1;
                cell_list[cell_id][0] = material_id;
              }
          }
        // Note: All other lines / entries are ignored
      }
  }

  template <int dim, int spacedim>
  std::vector<double>
  Abaqus_to_UCD<dim, spacedim>::get_global_node_numbers(
    const int face_cell_no,
    const int face_cell_face_no) const
  {
    std::vector<double> quad_node_list(GeometryInfo<dim>::vertices_per_face);

    // These orderings were reverse engineered by hand and may
    // conceivably be erroneous.
    // TODO: Currently one test (2d unstructured mesh) in the test
    // suite fails, presumably because of an ordering issue.
    if (dim == 2)
      {
        if (face_cell_face_no == 1)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][1];
            quad_node_list[1] = cell_list[face_cell_no - 1][2];
          }
        else if (face_cell_face_no == 2)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][2];
            quad_node_list[1] = cell_list[face_cell_no - 1][3];
          }
        else if (face_cell_face_no == 3)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][3];
            quad_node_list[1] = cell_list[face_cell_no - 1][4];
          }
        else if (face_cell_face_no == 4)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][4];
            quad_node_list[1] = cell_list[face_cell_no - 1][1];
          }
        else
          {
            AssertThrow(face_cell_face_no <= 4,
                        ExcMessage("Invalid face number in 2d"));
          }
      }
    else if (dim == 3)
      {
        if (face_cell_face_no == 1)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][1];
            quad_node_list[1] = cell_list[face_cell_no - 1][4];
            quad_node_list[2] = cell_list[face_cell_no - 1][3];
            quad_node_list[3] = cell_list[face_cell_no - 1][2];
          }
        else if (face_cell_face_no == 2)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][5];
            quad_node_list[1] = cell_list[face_cell_no - 1][8];
            quad_node_list[2] = cell_list[face_cell_no - 1][7];
            quad_node_list[3] = cell_list[face_cell_no - 1][6];
          }
        else if (face_cell_face_no == 3)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][1];
            quad_node_list[1] = cell_list[face_cell_no - 1][2];
            quad_node_list[2] = cell_list[face_cell_no - 1][6];
            quad_node_list[3] = cell_list[face_cell_no - 1][5];
          }
        else if (face_cell_face_no == 4)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][2];
            quad_node_list[1] = cell_list[face_cell_no - 1][3];
            quad_node_list[2] = cell_list[face_cell_no - 1][7];
            quad_node_list[3] = cell_list[face_cell_no - 1][6];
          }
        else if (face_cell_face_no == 5)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][3];
            quad_node_list[1] = cell_list[face_cell_no - 1][4];
            quad_node_list[2] = cell_list[face_cell_no - 1][8];
            quad_node_list[3] = cell_list[face_cell_no - 1][7];
          }
        else if (face_cell_face_no == 6)
          {
            quad_node_list[0] = cell_list[face_cell_no - 1][1];
            quad_node_list[1] = cell_list[face_cell_no - 1][5];
            quad_node_list[2] = cell_list[face_cell_no - 1][8];
            quad_node_list[3] = cell_list[face_cell_no - 1][4];
          }
        else
          {
            AssertThrow(face_cell_no <= 6,
                        ExcMessage("Invalid face number in 3d"));
          }
      }
    else
      {
        AssertThrow(dim == 2 || dim == 3, ExcNotImplemented());
      }

    return quad_node_list;
  }

  template <int dim, int spacedim>
  void
  Abaqus_to_UCD<dim, spacedim>::write_out_avs_ucd(std::ostream &output) const
  {
    // References:
    // http://www.dealii.org/developer/doxygen/deal.II/structGeometryInfo.html
    // http://people.scs.fsu.edu/~burkardt/data/ucd/ucd.html

    AssertThrow(output, ExcIO());

    // save old formatting options
    const boost::io::ios_base_all_saver formatting_saver(output);

    // Write out title - Note: No other commented text can be inserted below the
    // title in a UCD file
    output << "# Abaqus to UCD mesh conversion" << std::endl;
    output << "# Mesh type: AVS UCD" << std::endl;

    // ========================================================
    // ASCII UCD File Format
    // The input file cannot contain blank lines or lines with leading blanks.
    // Comments, if present, must precede all data in the file.
    // Comments within the data will cause read errors.
    // The general order of the data is as follows:
    // 1. Numbers defining the overall structure, including the number of nodes,
    //    the number of cells, and the length of the vector of data associated
    //    with the nodes, cells, and the model.
    //     e.g. 1:
    //        <num_nodes> <num_cells> <num_ndata> <num_cdata> <num_mdata>
    //     e.g. 2:
    //        n_elements = n_hex_cells + n_bc_quads + n_quad_cells + n_bc_edges
    //        outfile.write(str(n_nodes) + " " + str(n_elements) + " 0 0 0\n")
    // 2. For each node, its node id and the coordinates of that node in space.
    //    Node-ids must be integers, but any number including non sequential
    //    numbers can be used. Mid-edge nodes are treated like any other node.
    // 3. For each cell: its cell-id, material, cell type (hexahedral, pyramid,
    //    etc.), and the list of node-ids that correspond to each of the cell's
    //    vertices. The below table specifies the different cell types and the
    //    keyword used to represent them in the file.

    // Write out header
    output << node_list.size() << "\t" << (cell_list.size() + face_list.size())
           << "\t0\t0\t0" << std::endl;

    output.width(16);
    output.precision(8);

    // Write out node numbers
    // Loop over all nodes
    for (const auto &node : node_list)
      {
        // Node number
        output << node[0] << "\t";

        // Node coordinates
        output.setf(std::ios::scientific, std::ios::floatfield);
        for (unsigned int jj = 1; jj < spacedim + 1; ++jj)
          {
            // invoke tolerance -> set points close to zero equal to zero
            if (std::abs(node[jj]) > tolerance)
              output << static_cast<double>(node[jj]) << "\t";
            else
              output << 0.0 << "\t";
          }
        if (spacedim == 2)
          output << 0.0 << "\t";

        output << std::endl;
        output.unsetf(std::ios::floatfield);
      }

    // Write out cell node numbers
    for (unsigned int ii = 0; ii < cell_list.size(); ++ii)
      {
        output << ii + 1 << "\t" << cell_list[ii][0] << "\t"
               << (dim == 2 ? "quad" : "hex") << "\t";
        for (unsigned int jj = 1; jj < GeometryInfo<dim>::vertices_per_cell + 1;
             ++jj)
          output << cell_list[ii][jj] << "\t";

        output << std::endl;
      }

    // Write out quad node numbers
    for (unsigned int ii = 0; ii < face_list.size(); ++ii)
      {
        output << ii + 1 << "\t" << face_list[ii][0] << "\t"
               << (dim == 2 ? "line" : "quad") << "\t";
        for (unsigned int jj = 1; jj < GeometryInfo<dim>::vertices_per_face + 1;
             ++jj)
          output << face_list[ii][jj] << "\t";

        output << std::endl;
      }
  }
} // namespace


// explicit instantiations
#include "grid_in.inst"

DEAL_II_NAMESPACE_CLOSE
