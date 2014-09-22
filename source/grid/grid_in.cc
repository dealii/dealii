// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2014 by the deal.II authors
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


#include <deal.II/base/path_search.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/exceptions.h>

#include <deal.II/grid/grid_in.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_reordering.h>
#include <deal.II/grid/grid_tools.h>

#include <map>
#include <algorithm>
#include <fstream>
#include <functional>
#include <cctype>


#ifdef DEAL_II_WITH_NETCDF
#include <netcdfcpp.h>
#endif


DEAL_II_NAMESPACE_OPEN

template <int dim, int spacedim>
GridIn<dim, spacedim>::GridIn () :
  tria(0, typeid(*this).name()), default_format(ucd)
{}


template <int dim, int spacedim>
void GridIn<dim, spacedim>::attach_triangulation (Triangulation<dim, spacedim> &t)
{
  tria = &t;
}



template<int dim, int spacedim>
void GridIn<dim, spacedim>::read_vtk(std::istream &in)
{
  Assert((dim == 2)||(dim == 3), ExcNotImplemented());
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
        getline(in,line);
        if (i != 1)
          AssertThrow (line.compare(text[i]) == 0,
                       ExcMessage(std::string("While reading VTK file, failed to find a header line with text <") +
                                  text[i] + ">"));
      }
  }

  ///////////////////Declaring storage and mappings//////////////////

  std::vector< Point<spacedim> > vertices;//vector of vertices
  std::vector< CellData<dim> > cells;//vector of cells
  SubCellData subcelldata;//subcell data that includes bounds and material IDs.
  std::map<int, int> vertex_indices; // # vert in unv (key) ---> # vert in deal.II (value)
  std::map<int, int> cell_indices; // # cell in unv (key) ---> # cell in deal.II (value)
  std::map<int, int> quad_indices; // # quad in unv (key) ---> # quad in deal.II (value)
  std::map<int, int> line_indices; // # line in unv(key) ---> # line in deal.II (value)

  unsigned int no_vertices, no_quads=0, no_lines=0;

  std::string keyword;

  in >> keyword;

  //////////////////Processing the POINTS section///////////////

  if (keyword == "POINTS")
    {
      in>>no_vertices;// taking the no. of vertices
      in.ignore(256, '\n');//ignoring the number beside the total no. of points.

      for (unsigned int count = 0; count < no_vertices; count++) //loop to read three values till the no . vertices is satisfied
        {
          // VTK format always specifies vertex coordinates with 3 components
          Point<3> x;
          in >> x(0) >> x(1) >> x(2);

          vertices.push_back(Point<spacedim>());
          for (unsigned int d=0; d<spacedim; ++d)
            vertices.back()(d) = x(d);

          vertex_indices[count] = count;
        }
    }

  else
    AssertThrow (false,
                 ExcMessage ("While reading VTK file, failed to find POINTS section"));


  //////////////////ignoring space between points and cells sections////////////////////
  std::string checkline;
  int no;
  in.ignore(256, '\n');//this move pointer to the next line ignoring unwanted no.
  no = in.tellg();
  getline(in,checkline);
  if (checkline.compare("") != 0)
    {
      in.seekg(no);
    }

  in >> keyword;

  unsigned int total_cells, no_cells = 0, type;// declaring counters, refer to the order of declaring variables for an idea of what is what!

  ///////////////////Processing the CELLS section that contains cells(cells) and bound_quads(subcelldata)///////////////////////

  if (keyword == "CELLS")
    {
      in>>total_cells;
      in.ignore(256,'\n');

      if (dim == 3)
        {
          for (unsigned int count = 0; count < total_cells; count++)
            {
              in>>type;

              if (type == 8)
                {

                  cells.push_back(CellData<dim>());

                  for (unsigned int j = 0; j < type; j++) //loop to feed data
                    in >> cells.back().vertices[j];


                  cells.back().material_id = 0;

                  for (unsigned int j = 0; j < type; j++) //loop to feed the data of the vertices to the cell
                    {
                      cells.back().vertices[j] = vertex_indices[cells.back().vertices[j]];
                    }
                  cell_indices[count] = count;
                  no_cells++;
                }

              else if ( type == 4)
                {

                  subcelldata.boundary_quads.push_back(CellData<2>());

                  for (unsigned int j = 0; j < type; j++) //loop to feed the data to the boundary
                    {
                      in >> subcelldata.boundary_quads.back().vertices[j];
                    }
                  subcelldata.boundary_quads.back().material_id = 0;
                  for (unsigned int j = 0; j < type; j++)
                    {
                      subcelldata.boundary_quads.back().vertices[j] = vertex_indices[subcelldata.boundary_quads.back().vertices[j]];
                    }
                  quad_indices[no_quads] = no_quads + 1;
                  no_quads++;
                }

              else
                AssertThrow (false,
                             ExcMessage ("While reading VTK file, unknown file type encountered"));
            }
        }

      else if (dim == 2)
        {
          for (unsigned int count = 0; count < total_cells; count++)
            {
              in>>type;

              if (type == 4)
                {
                  cells.push_back(CellData<dim>());

                  for (unsigned int j = 0; j < type; j++) //loop to feed data
                    in >> cells.back().vertices[j];

                  cells.back().material_id = 0;

                  for (unsigned int j = 0; j < type; j++) //loop to feed the data of the vertices to the cell
                    {
                      cells.back().vertices[j] = vertex_indices[cells.back().vertices[j]];
                    }
                  cell_indices[count] = count;
                  no_cells++;
                }

              else if (type == 2)
                {
                  //If this is encountered, the pointer comes out of the loop
                  //and starts processing boundaries.
                  subcelldata.boundary_lines.push_back(CellData<1>());

                  for (unsigned int j = 0; j < type; j++) //loop to feed the data to the boundary
                    {
                      in >> subcelldata.boundary_lines.back().vertices[j];
                    }
                  subcelldata.boundary_lines.back().material_id = 0;
                  for (unsigned int j = 0; j < type; j++)
                    {
                      subcelldata.boundary_lines.back().vertices[j] = vertex_indices[subcelldata.boundary_lines.back().vertices[j]];
                    }
                  line_indices[no_lines] = no_lines + 1;
                  no_lines++;
                }

              else
                AssertThrow (false,
                             ExcMessage ("While reading VTK file, unknown file type encountered"));
            }
        }
      else
        AssertThrow (false,
                     ExcMessage ("While reading VTK file, failed to find CELLS section"));

      /////////////////////Processing the CELL_TYPES section////////////////////////

      in >> keyword;

      if (keyword == "CELL_TYPES")//Entering the cell_types section and ignoring data.
        {
          in.ignore(256, '\n');

          while (!in.eof())
            {
              in>>keyword;
              if (keyword != "12" && keyword != "9")
                {
                  break;
                }
            }
        }

      ////////////////////////Processing the CELL_DATA section/////////////////////////////

      if (keyword == "CELL_DATA")
        {
          int no_ids;
          in>>no_ids;

          std::string linenew;
          std::string textnew[2];
          textnew[0] = "SCALARS MaterialID double";
          textnew[1] = "LOOKUP_TABLE default";

          in.ignore(256, '\n');

          for (unsigned int i = 0; i < 2; i++)
            {
              getline(in, linenew);
              if (i == 0)
                if (linenew.size() > textnew[0].size())
                  linenew.resize(textnew[0].size());

              AssertThrow (linenew.compare(textnew[i]) == 0,
                           ExcMessage (std::string("While reading VTK file, failed to find <") +
                                       textnew[i] + "> section"));
            }

          for (unsigned int i = 0; i < no_cells; i++) //assigning IDs to cells.
            {
              int id;
              in>>id;
              cells[cell_indices[i]].material_id = id;
            }

          if (dim == 3)
            {
              for (unsigned int i = 0; i < no_quads; i++) //assigning IDs to bounds.
                {
                  int id;
                  in>>id;
                  subcelldata.boundary_quads[quad_indices[i]].material_id = id;
                }
            }
          else if (dim == 2)
            {
              for (unsigned int i = 0; i < no_lines; i++) //assigning IDs to bounds.
                {
                  int id;
                  in>>id;
                  subcelldata.boundary_lines[line_indices[i]].material_id = id;
                }
            }
        }

      Assert(subcelldata.check_consistency(dim), ExcInternalError());

      GridTools::delete_unused_vertices(vertices,
                                        cells,
                                        subcelldata);

      if (dim == spacedim)
        GridReordering<dim, spacedim>::invert_all_cells_of_negative_grid(vertices,
            cells);

      GridReordering<dim, spacedim>::reorder_cells(cells);
      tria->create_triangulation_compatibility(vertices,
                                               cells,
                                               subcelldata);

      return;
    }
  else
    AssertThrow (false,
                 ExcMessage ("While reading VTK file, failed to find CELLS section"));
}



template<int dim, int spacedim>
void GridIn<dim, spacedim>::read_unv(std::istream &in)
{
  Assert(tria != 0, ExcNoTriangulationSelected());
  Assert((dim == 2)||(dim == 3), ExcNotImplemented());

  AssertThrow(in, ExcIO());
  skip_comment_lines(in, '#'); // skip comments (if any) at beginning of file

  int tmp;

  AssertThrow(in, ExcIO());
  in >> tmp;
  AssertThrow(in, ExcIO());
  in >> tmp;

  AssertThrow(tmp == 2411, ExcUnknownSectionType(tmp)); // section 2411 describes vertices http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse/unv_2411.htm/

  std::vector< Point<spacedim> > vertices; // vector of vertex coordinates
  std::map<int, int> vertex_indices; // # vert in unv (key) ---> # vert in deal.II (value)

  int no_vertex = 0; // deal.II

  while (tmp != -1) // we do until reach end of 2411
    {
      int no; // unv
      int dummy;
      double x[3];

      AssertThrow(in, ExcIO());
      in >> no;

      tmp = no;
      if (tmp == -1)
        break;

      in >> dummy >> dummy >> dummy;

      AssertThrow(in, ExcIO());
      in >> x[0] >> x[1] >> x[2];

      vertices.push_back(Point<spacedim>());

      for (unsigned int d = 0; d < spacedim; d++)
        vertices.back()(d) = x[d];

      vertex_indices[no] = no_vertex;

      no_vertex++;
    }

  AssertThrow(in, ExcIO());
  in >> tmp;
  AssertThrow(in, ExcIO());
  in >> tmp;

  AssertThrow(tmp == 2412, ExcUnknownSectionType(tmp)); // section 2412 describes elements http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse/unv_2412.htm/

  std::vector< CellData<dim> > cells; // vector of cells
  SubCellData subcelldata;

  std::map<int, int> cell_indices; // # cell in unv (key) ---> # cell in deal.II (value)
  std::map<int, int> line_indices; // # line in unv (key) ---> # line in deal.II (value)
  std::map<int, int> quad_indices; // # quad in unv (key) ---> # quad in deal.II (value)

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

      AssertThrow((type == 11)||(type == 44)||(type == 115), ExcUnknownElementType(type));

      if ( ((type == 44)&&(dim == 2)) || ((type == 115)&&(dim == 3)) ) // cell
        {
          cells.push_back(CellData<dim>());

          AssertThrow(in, ExcIO());
          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; v++)
            in >> cells.back().vertices[v];

          cells.back().material_id = 0;

          for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_cell; v++)
            cells.back().vertices[v] = vertex_indices[cells.back().vertices[v]];

          cell_indices[no] = no_cell;

          no_cell++;
        }
      else if ( ((type == 11)&&(dim == 2)) || ((type == 11)&&(dim == 3)) ) // boundary line
        {
          AssertThrow(in, ExcIO());
          in >> dummy >> dummy >> dummy;

          subcelldata.boundary_lines.push_back(CellData<1>());

          AssertThrow(in, ExcIO());
          for (unsigned int v = 0; v < 2; v++)
            in >> subcelldata.boundary_lines.back().vertices[v];

          subcelldata.boundary_lines.back().material_id = 0;

          for (unsigned int v = 0; v < 2; v++)
            subcelldata.boundary_lines.back().vertices[v] = vertex_indices[subcelldata.boundary_lines.back().vertices[v]];

          line_indices[no] = no_line;

          no_line++;
        }
      else if ( (type == 44) && (dim == 3) ) // boundary quad
        {
          subcelldata.boundary_quads.push_back(CellData<2>());

          AssertThrow(in, ExcIO());
          for (unsigned int v = 0; v < 4; v++)
            in >> subcelldata.boundary_quads.back().vertices[v];

          subcelldata.boundary_quads.back().material_id = 0;

          for (unsigned int v = 0; v < 4; v++)
            subcelldata.boundary_quads.back().vertices[v] = vertex_indices[subcelldata.boundary_quads.back().vertices[v]];

          quad_indices[no] = no_quad;

          no_quad++;
        }
    }

// note that so far all materials and bcs are explicitly set to 0
// if we do not need more info on materials and bcs - this is end of file
// if we do - section 2467 comes

  in >> tmp; // tmp can be either -1 or end-of-file

  if ( !in.eof() )
    {
      AssertThrow(in, ExcIO());
      in >> tmp;

      AssertThrow(tmp == 2467, ExcUnknownSectionType(tmp)); // section 2467 describes (materials - first and bcs - second) or (bcs - first and materials - second) - sequence depends on which group is created first
      // http://www.sdrl.uc.edu/universal-file-formats-for-modal-analysis-testing-1/file-format-storehouse/unv_2467.htm/

      while (tmp != -1) // we do until reach end of 2467
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

          in >> dummy >> dummy >> dummy >> dummy >> dummy >> dummy >> n_entities;

          AssertThrow(in, ExcIO());
          in >> id;

          const unsigned int n_lines = (n_entities%2 == 0)?(n_entities/2):((n_entities+1)/2);

          for (unsigned int line = 0; line < n_lines; line++)
            {
              unsigned int n_fragments;

              if (line == n_lines-1)
                n_fragments = (n_entities%2 == 0)?(2):(1);
              else
                n_fragments = 2;

              for (unsigned int no_fragment = 0; no_fragment < n_fragments; no_fragment++)
                {
                  AssertThrow(in, ExcIO());
                  in >> dummy >> no >> dummy >> dummy;

                  if ( cell_indices.count(no) > 0 ) // cell - material
                    cells[cell_indices[no]].material_id = id;

                  if ( line_indices.count(no) > 0 ) // boundary line - bc
                    subcelldata.boundary_lines[line_indices[no]].material_id = id;

                  if ( quad_indices.count(no) > 0 ) // boundary quad - bc
                    subcelldata.boundary_quads[quad_indices[no]].material_id = id;
                }
            }
        }
    }

  Assert(subcelldata.check_consistency(dim), ExcInternalError());

  GridTools::delete_unused_vertices(vertices,
                                    cells,
                                    subcelldata);

  if (dim == spacedim)
    GridReordering<dim, spacedim>::invert_all_cells_of_negative_grid(vertices,
        cells);

  GridReordering<dim, spacedim>::reorder_cells(cells);

  tria->create_triangulation_compatibility(vertices,
                                           cells,
                                           subcelldata);
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::read_ucd (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  // skip comments at start of file
  skip_comment_lines (in, '#');


  unsigned int n_vertices;
  unsigned int n_cells;
  int dummy;

  in >> n_vertices
     >> n_cells
     >> dummy         // number of data vectors
     >> dummy         // cell data
     >> dummy;        // model data
  AssertThrow (in, ExcIO());

  // set up array of vertices
  std::vector<Point<spacedim> >     vertices (n_vertices);
  // set up mapping between numbering
  // in ucd-file (key) and in the
  // vertices vector
  std::map<int,int> vertex_indices;

  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
      int vertex_number;
      double x[3];

      // read vertex
      AssertThrow (in, ExcIO());
      in >> vertex_number
         >> x[0] >> x[1] >> x[2];

      // store vertex
      for (unsigned int d=0; d<spacedim; ++d)
        vertices[vertex](d) = x[d];
      // store mapping; note that
      // vertices_indices[i] is automatically
      // created upon first usage
      vertex_indices[vertex_number] = vertex;
    };

  // set up array of cells
  std::vector<CellData<dim> > cells;
  SubCellData                 subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      // note that since in the input
      // file we found the number of
      // cells at the top, there
      // should still be input here,
      // so check this:
      AssertThrow (in, ExcIO());

      std::string cell_type;

      // we use an unsigned int because we
      // fill this variable through an read-in process
      unsigned int material_id;

      in >> dummy          // cell number
         >> material_id;
      in >> cell_type;

      if (((cell_type == "line") && (dim == 1)) ||
          ((cell_type == "quad") && (dim == 2)) ||
          ((cell_type == "hex" ) && (dim == 3)))
        // found a cell
        {
          // allocate and read indices
          cells.push_back (CellData<dim>());
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
            in >> cells.back().vertices[i];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::material_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::material_id>::max()));
          // we use only material_ids in the range from 0 to numbers::invalid_material_id-1
          Assert(material_id < numbers::invalid_material_id,
                 ExcIndexRange(material_id,0,numbers::invalid_material_id));

          cells.back().material_id = static_cast<types::material_id>(material_id);

          // transform from ucd to
          // consecutive numbering
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
            if (vertex_indices.find (cells.back().vertices[i]) != vertex_indices.end())
              // vertex with this index exists
              cells.back().vertices[i] = vertex_indices[cells.back().vertices[i]];
            else
              {
                // no such vertex index
                AssertThrow (false,
                             ExcInvalidVertexIndex(cell, cells.back().vertices[i]));

                cells.back().vertices[i] = numbers::invalid_unsigned_int;
              }
        }
      else if ((cell_type == "line") && ((dim == 2) || (dim == 3)))
        // boundary info
        {
          subcelldata.boundary_lines.push_back (CellData<1>());
          in >> subcelldata.boundary_lines.back().vertices[0]
             >> subcelldata.boundary_lines.back().vertices[1];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
          Assert(material_id < numbers::internal_face_boundary_id,
                 ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));

          subcelldata.boundary_lines.back().boundary_id
            = static_cast<types::boundary_id>(material_id);

          // transform from ucd to
          // consecutive numbering
          for (unsigned int i=0; i<2; ++i)
            if (vertex_indices.find (subcelldata.boundary_lines.back().vertices[i]) !=
                vertex_indices.end())
              // vertex with this index exists
              subcelldata.boundary_lines.back().vertices[i]
                = vertex_indices[subcelldata.boundary_lines.back().vertices[i]];
            else
              {
                // no such vertex index
                AssertThrow (false,
                             ExcInvalidVertexIndex(cell,
                                                   subcelldata.boundary_lines.back().vertices[i]));
                subcelldata.boundary_lines.back().vertices[i]
                  = numbers::invalid_unsigned_int;
              };
        }
      else if ((cell_type == "quad") && (dim == 3))
        // boundary info
        {
          subcelldata.boundary_quads.push_back (CellData<2>());
          in >> subcelldata.boundary_quads.back().vertices[0]
             >> subcelldata.boundary_quads.back().vertices[1]
             >> subcelldata.boundary_quads.back().vertices[2]
             >> subcelldata.boundary_quads.back().vertices[3];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
          Assert(material_id < numbers::internal_face_boundary_id,
                 ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));

          subcelldata.boundary_quads.back().boundary_id
            = static_cast<types::boundary_id>(material_id);

          // transform from ucd to
          // consecutive numbering
          for (unsigned int i=0; i<4; ++i)
            if (vertex_indices.find (subcelldata.boundary_quads.back().vertices[i]) !=
                vertex_indices.end())
              // vertex with this index exists
              subcelldata.boundary_quads.back().vertices[i]
                = vertex_indices[subcelldata.boundary_quads.back().vertices[i]];
            else
              {
                // no such vertex index
                Assert (false,
                        ExcInvalidVertexIndex(cell,
                                              subcelldata.boundary_quads.back().vertices[i]));
                subcelldata.boundary_quads.back().vertices[i] =
                  numbers::invalid_unsigned_int;
              };

        }
      else
        // cannot read this
        AssertThrow (false, ExcUnknownIdentifier(cell_type));
    };


  // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  // ... and cells
  if (dim==spacedim)
    GridReordering<dim,spacedim>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::read_dbmesh (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  Assert (dim==2, ExcNotImplemented());

  AssertThrow (in, ExcIO());

  // skip comments at start of file
  skip_comment_lines (in, '#');


  // first read in identifier string
  std::string line;
  getline (in, line);

  AssertThrow (line=="MeshVersionFormatted 0",
               ExcInvalidDBMESHInput(line));

  skip_empty_lines (in);

  // next read dimension
  getline (in, line);
  AssertThrow (line=="Dimension", ExcInvalidDBMESHInput(line));
  unsigned int dimension;
  in >> dimension;
  AssertThrow (dimension == dim, ExcDBMESHWrongDimension(dimension));
  skip_empty_lines (in);

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
  while (getline(in,line), line.find("# END")==std::string::npos)
    ;
  skip_empty_lines (in);


  // now read vertices
  getline (in, line);
  AssertThrow (line=="Vertices", ExcInvalidDBMESHInput(line));

  unsigned int n_vertices;
  double dummy;

  in >> n_vertices;
  std::vector<Point<spacedim> >     vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
      // read vertex coordinates
      for (unsigned int d=0; d<dim; ++d)
        in >> vertices[vertex][d];
      // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);

  // read edges. we ignore them at
  // present, so just read them and
  // discard the input
  getline (in, line);
  AssertThrow (line=="Edges", ExcInvalidDBMESHInput(line));

  unsigned int n_edges;
  in >> n_edges;
  for (unsigned int edge=0; edge<n_edges; ++edge)
    {
      // read vertex indices
      in >> dummy >> dummy;
      // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);



  // read cracked edges (whatever
  // that may be). we ignore them at
  // present, so just read them and
  // discard the input
  getline (in, line);
  AssertThrow (line=="CrackedEdges", ExcInvalidDBMESHInput(line));

  in >> n_edges;
  for (unsigned int edge=0; edge<n_edges; ++edge)
    {
      // read vertex indices
      in >> dummy >> dummy;
      // read Ref phi_i, whatever that may be
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);


  // now read cells.
  // set up array of cells
  getline (in, line);
  AssertThrow (line=="Quadrilaterals", ExcInvalidDBMESHInput(line));

  std::vector<CellData<dim> > cells;
  SubCellData            subcelldata;
  unsigned int n_cells;
  in >> n_cells;
  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      // read in vertex numbers. they
      // are 1-based, so subtract one
      cells.push_back (CellData<dim>());
      for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
        {
          in >> cells.back().vertices[i];

          AssertThrow ((cells.back().vertices[i] >= 1)
                       &&
                       (static_cast<unsigned int>(cells.back().vertices[i]) <= vertices.size()),
                       ExcInvalidVertexIndex(cell, cells.back().vertices[i]));

          --cells.back().vertices[i];
        };

      // read and discard Ref phi_i
      in >> dummy;
    };
  AssertThrow (in, ExcInvalidDBMeshFormat());

  skip_empty_lines(in);


  // then there are again a whole lot
  // of fields of which I have no
  // clue what they mean. skip them
  // all and leave the interpretation
  // to other implementors...
  while (getline(in,line), ((line.find("End")==std::string::npos) && (in)))
    ;
  // ok, so we are not at the end of
  // the file, that's it, mostly


  // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  // ...and cells
  GridReordering<dim,spacedim>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::read_xda (std::istream &)
{
  Assert (false, ExcNotImplemented());
}



template <>
void GridIn<2>::read_xda (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  std::string line;
  // skip comments at start of file
  getline (in, line);


  unsigned int n_vertices;
  unsigned int n_cells;

  // read cells, throw away rest of line
  in >> n_cells;
  getline (in, line);

  in >> n_vertices;
  getline (in, line);

  // ignore following 8 lines
  for (unsigned int i=0; i<8; ++i)
    getline (in, line);

  // set up array of cells
  std::vector<CellData<2> > cells (n_cells);
  SubCellData subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      // note that since in the input
      // file we found the number of
      // cells at the top, there
      // should still be input here,
      // so check this:
      AssertThrow (in, ExcIO());
      Assert (GeometryInfo<2>::vertices_per_cell == 4,
              ExcInternalError());

      for (unsigned int i=0; i<4; ++i)
        in >> cells[cell].vertices[i];
    };



  // set up array of vertices
  std::vector<Point<2> > vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
      double x[3];

      // read vertex
      in >> x[0] >> x[1] >> x[2];

      // store vertex
      for (unsigned int d=0; d<2; ++d)
        vertices[vertex](d) = x[d];
    };
  AssertThrow (in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  // ... and cells
  GridReordering<2>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<2>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}



template <>
void GridIn<3>::read_xda (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  static const unsigned int xda_to_dealII_map[] = {0,1,5,4,3,2,6,7};

  std::string line;
  // skip comments at start of file
  getline (in, line);


  unsigned int n_vertices;
  unsigned int n_cells;

  // read cells, throw away rest of line
  in >> n_cells;
  getline (in, line);

  in >> n_vertices;
  getline (in, line);

  // ignore following 8 lines
  for (unsigned int i=0; i<8; ++i)
    getline (in, line);

  // set up array of cells
  std::vector<CellData<3> > cells (n_cells);
  SubCellData subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      // note that since in the input
      // file we found the number of
      // cells at the top, there
      // should still be input here,
      // so check this:
      AssertThrow (in, ExcIO());
      Assert(GeometryInfo<3>::vertices_per_cell == 8,
             ExcInternalError());

      unsigned int xda_ordered_nodes[8];

      for (unsigned int i=0; i<8; ++i)
        in >> xda_ordered_nodes[i];

      for (unsigned int i=0; i<8; i++)
        cells[cell].vertices[i] = xda_ordered_nodes[xda_to_dealII_map[i]];
    };



  // set up array of vertices
  std::vector<Point<3> > vertices (n_vertices);
  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
      double x[3];

      // read vertex
      in >> x[0] >> x[1] >> x[2];

      // store vertex
      for (unsigned int d=0; d<3; ++d)
        vertices[vertex](d) = x[d];
    };
  AssertThrow (in, ExcIO());

  // do some clean-up on vertices...
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  // ... and cells
  GridReordering<3>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<3>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}




template <int dim, int spacedim>
void GridIn<dim, spacedim>::read_msh (std::istream &in)
{
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  unsigned int n_vertices;
  unsigned int n_cells;
  unsigned int dummy;
  std::string line;

  in >> line;

  // first determine file format
  unsigned int gmsh_file_format = 0;
  if (line == "$NOD")
    gmsh_file_format = 1;
  else if (line == "$MeshFormat")
    gmsh_file_format = 2;
  else
    AssertThrow (false, ExcInvalidGMSHInput(line));

  // if file format is 2 or greater
  // then we also have to read the
  // rest of the header
  if (gmsh_file_format == 2)
    {
      double version;
      unsigned int file_type, data_size;

      in >> version >> file_type >> data_size;

      Assert ( (version >= 2.0) &&
               (version <= 2.2), ExcNotImplemented());
      Assert (file_type == 0, ExcNotImplemented());
      Assert (data_size == sizeof(double), ExcNotImplemented());

      // read the end of the header
      // and the first line of the
      // nodes description to synch
      // ourselves with the format 1
      // handling above
      in >> line;
      AssertThrow (line == "$EndMeshFormat",
                   ExcInvalidGMSHInput(line));

      in >> line;
      // if the next block is of kind
      // $PhysicalNames, ignore it
      if (line == "$PhysicalNames")
        {
          do
            {
              in >> line;
            }
          while (line != "$EndPhysicalNames");
          in >> line;
        }

      // but the next thing should,
      // in any case, be the list of
      // nodes:
      AssertThrow (line == "$Nodes",
                   ExcInvalidGMSHInput(line));
    }

  // now read the nodes list
  in >> n_vertices;
  std::vector<Point<spacedim> >     vertices (n_vertices);
  // set up mapping between numbering
  // in msh-file (nod) and in the
  // vertices vector
  std::map<int,int> vertex_indices;

  for (unsigned int vertex=0; vertex<n_vertices; ++vertex)
    {
      int vertex_number;
      double x[3];

      // read vertex
      in >> vertex_number
         >> x[0] >> x[1] >> x[2];

      for (unsigned int d=0; d<spacedim; ++d)
        vertices[vertex](d) = x[d];
      // store mapping
      vertex_indices[vertex_number] = vertex;
    }

  // Assert we reached the end of the block
  in >> line;
  static const std::string end_nodes_marker[] = {"$ENDNOD", "$EndNodes" };
  AssertThrow (line==end_nodes_marker[gmsh_file_format-1],
               ExcInvalidGMSHInput(line));

  // Now read in next bit
  in >> line;
  static const std::string begin_elements_marker[] = {"$ELM", "$Elements" };
  AssertThrow (line==begin_elements_marker[gmsh_file_format-1],
               ExcInvalidGMSHInput(line));

  in >> n_cells;

  // set up array of cells
  std::vector<CellData<dim> > cells;
  SubCellData                 subcelldata;

  for (unsigned int cell=0; cell<n_cells; ++cell)
    {
      // note that since in the input
      // file we found the number of
      // cells at the top, there
      // should still be input here,
      // so check this:
      AssertThrow (in, ExcIO());

      unsigned int cell_type;
      unsigned int material_id;
      unsigned int nod_num;

      /*
        For file format version 1, the format of each cell is as follows:
          elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list

        However, for version 2, the format reads like this:
          elm-number elm-type number-of-tags < tag > ... node-number-list

        In the following, we will ignore the element number (we simply enumerate
        them in the order in which we read them, and we will take reg-phys
        (version 1) or the first tag (version 2, if any tag is given at all) as
        material id.
      */

      in >> dummy          // ELM-NUMBER
         >> cell_type;     // ELM-TYPE

      switch (gmsh_file_format)
        {
        case 1:
        {
          in >> material_id  // REG-PHYS
             >> dummy        // reg_elm
             >> nod_num;
          break;
        }

        case 2:
        {
          // read the tags; ignore
          // all but the first one
          unsigned int n_tags;
          in >> n_tags;
          if (n_tags > 0)
            in >> material_id;
          else
            material_id = 0;

          for (unsigned int i=1; i<n_tags; ++i)
            in >> dummy;

          nod_num = GeometryInfo<dim>::vertices_per_cell;

          break;
        }

        default:
          AssertThrow (false, ExcNotImplemented());
        }


      /*       `ELM-TYPE'
               defines the geometrical type of the N-th element:
               `1'
               Line (2 nodes, 1 edge).

               `3'
               Quadrangle (4 nodes, 4 edges).

               `5'
               Hexahedron (8 nodes, 12 edges, 6 faces).

               `15'
               Point (1 node).
      */

      if (((cell_type == 1) && (dim == 1)) ||
          ((cell_type == 3) && (dim == 2)) ||
          ((cell_type == 5) && (dim == 3)))
        // found a cell
        {
          AssertThrow (nod_num == GeometryInfo<dim>::vertices_per_cell,
                       ExcMessage ("Number of nodes does not coincide with the "
                                   "number required for this object"));

          // allocate and read indices
          cells.push_back (CellData<dim>());
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
            in >> cells.back().vertices[i];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::material_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::material_id>::max()));
          // we use only material_ids in the range from 0 to numbers::invalid_material_id-1
          Assert(material_id < numbers::invalid_material_id,
                 ExcIndexRange(material_id,0,numbers::invalid_material_id));

          cells.back().material_id = static_cast<types::material_id>(material_id);

          // transform from ucd to
          // consecutive numbering
          for (unsigned int i=0; i<GeometryInfo<dim>::vertices_per_cell; ++i)
            {
              AssertThrow (vertex_indices.find (cells.back().vertices[i]) !=
                           vertex_indices.end(),
                           ExcInvalidVertexIndex(cell, cells.back().vertices[i]));

              // vertex with this index exists
              cells.back().vertices[i] = vertex_indices[cells.back().vertices[i]];
            }
        }
      else if ((cell_type == 1) && ((dim == 2) || (dim == 3)))
        // boundary info
        {
          subcelldata.boundary_lines.push_back (CellData<1>());
          in >> subcelldata.boundary_lines.back().vertices[0]
             >> subcelldata.boundary_lines.back().vertices[1];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
          Assert(material_id < numbers::internal_face_boundary_id,
                 ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));

          subcelldata.boundary_lines.back().boundary_id
            = static_cast<types::boundary_id>(material_id);

          // transform from ucd to
          // consecutive numbering
          for (unsigned int i=0; i<2; ++i)
            if (vertex_indices.find (subcelldata.boundary_lines.back().vertices[i]) !=
                vertex_indices.end())
              // vertex with this index exists
              subcelldata.boundary_lines.back().vertices[i]
                = vertex_indices[subcelldata.boundary_lines.back().vertices[i]];
            else
              {
                // no such vertex index
                AssertThrow (false,
                             ExcInvalidVertexIndex(cell,
                                                   subcelldata.boundary_lines.back().vertices[i]));
                subcelldata.boundary_lines.back().vertices[i] =
                  numbers::invalid_unsigned_int;
              };
        }
      else if ((cell_type == 3) && (dim == 3))
        // boundary info
        {
          subcelldata.boundary_quads.push_back (CellData<2>());
          in >> subcelldata.boundary_quads.back().vertices[0]
             >> subcelldata.boundary_quads.back().vertices[1]
             >> subcelldata.boundary_quads.back().vertices[2]
             >> subcelldata.boundary_quads.back().vertices[3];

          // to make sure that the cast wont fail
          Assert(material_id<= std::numeric_limits<types::boundary_id>::max(),
                 ExcIndexRange(material_id,0,std::numeric_limits<types::boundary_id>::max()));
          // we use only boundary_ids in the range from 0 to numbers::internal_face_boundary_id-1
          Assert(material_id < numbers::internal_face_boundary_id,
                 ExcIndexRange(material_id,0,numbers::internal_face_boundary_id));

          subcelldata.boundary_quads.back().boundary_id
            = static_cast<types::boundary_id>(material_id);

          // transform from gmsh to
          // consecutive numbering
          for (unsigned int i=0; i<4; ++i)
            if (vertex_indices.find (subcelldata.boundary_quads.back().vertices[i]) !=
                vertex_indices.end())
              // vertex with this index exists
              subcelldata.boundary_quads.back().vertices[i]
                = vertex_indices[subcelldata.boundary_quads.back().vertices[i]];
            else
              {
                // no such vertex index
                Assert (false,
                        ExcInvalidVertexIndex(cell,
                                              subcelldata.boundary_quads.back().vertices[i]));
                subcelldata.boundary_quads.back().vertices[i] =
                  numbers::invalid_unsigned_int;
              };

        }
      else if (cell_type == 15)
        {
          // Ignore vertices
          // but read the
          // number of nodes
          // given
          switch (gmsh_file_format)
            {
            case 1:
            {
              for (unsigned int i=0; i<nod_num; ++i)
                in >> dummy;
              break;
            }
            case 2:
            {
              in >> dummy;
              break;
            }
            }
        }
      else
        // cannot read this, so throw
        // an exception. treat
        // triangles and tetrahedra
        // specially since this
        // deserves a more explicit
        // error message
        {
          AssertThrow (cell_type != 2,
                       ExcMessage("Found triangles while reading a file "
                                  "in gmsh format. deal.II does not "
                                  "support triangles"));
          AssertThrow (cell_type != 11,
                       ExcMessage("Found tetrahedra while reading a file "
                                  "in gmsh format. deal.II does not "
                                  "support tetrahedra"));

          AssertThrow (false, ExcGmshUnsupportedGeometry(cell_type));
        }
    };

  // Assert we reached the end of the block
  in >> line;
  static const std::string end_elements_marker[] = {"$ENDELM", "$EndElements" };
  AssertThrow (line==end_elements_marker[gmsh_file_format-1],
               ExcInvalidGMSHInput(line));

  // check that no forbidden arrays are used
  Assert (subcelldata.check_consistency(dim), ExcInternalError());

  AssertThrow (in, ExcIO());

  // check that we actually read some
  // cells.
  AssertThrow(cells.size() > 0, ExcGmshNoCellInformation());

  // do some clean-up on
  // vertices...
  GridTools::delete_unused_vertices (vertices, cells, subcelldata);
  // ... and cells
  if (dim==spacedim)
    GridReordering<dim,spacedim>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}


template <>
void GridIn<1>::read_netcdf (const std::string &)
{
  AssertThrow(false, ExcImpossibleInDim(1));
}

template <>
void GridIn<1,2>::read_netcdf (const std::string &)
{
  AssertThrow(false, ExcImpossibleInDim(1));
}


template <>
void GridIn<2, 3>::read_netcdf (const std::string &)
{
  Assert(false, ExcNotImplemented());
}

template <>
void GridIn<2>::read_netcdf (const std::string &filename)
{
#ifndef DEAL_II_WITH_NETCDF
  // do something with unused
  // filename
  filename.c_str();
  AssertThrow(false, ExcNeedsNetCDF());
#else
  const unsigned int dim=2;
  const unsigned int spacedim=2;
  const bool output=false;
  Assert (tria != 0, ExcNoTriangulationSelected());
  // this function assumes the TAU
  // grid format.
  //
  // This format stores 2d grids as
  // 3d grids. In particular, a 2d
  // grid of n_cells quadrilaterals
  // in the y=0 plane is duplicated
  // to y=1 to build n_cells
  // hexaeders.  The surface
  // quadrilaterals of this 3d grid
  // are marked with boundary
  // marker. In the following we read
  // in all data required, find the
  // boundary marker associated with
  // the plane y=0, and extract the
  // corresponding 2d data to build a
  // Triangulation<2>.

  // In the following, we assume that
  // the 2d grid lies in the x-z
  // plane (y=0). I.e. we choose:
  // point[coord]=0, with coord=1
  const unsigned int coord=1;
  // Also x-y-z (0-1-2) point
  // coordinates will be transformed
  // to x-y (x2d-y2d) coordinates.
  // With coord=1 as above, we have
  // x-z (0-2) -> (x2d-y2d)
  const unsigned int x2d=0;
  const unsigned int y2d=2;
  // For the case, the 2d grid lies
  // in x-y or y-z plane instead, the
  // following code must be extended
  // to find the right value for
  // coord, and setting x2d and y2d
  // accordingly.

  // First, open the file
  NcFile nc (filename.c_str());
  AssertThrow(nc.is_valid(), ExcIO());

  // then read n_cells
  NcDim *elements_dim=nc.get_dim("no_of_elements");
  AssertThrow(elements_dim->is_valid(), ExcIO());
  const unsigned int n_cells=elements_dim->size();

  // then we read
  //   int marker(no_of_markers)
  NcDim *marker_dim=nc.get_dim("no_of_markers");
  AssertThrow(marker_dim->is_valid(), ExcIO());
  const unsigned int n_markers=marker_dim->size();

  NcVar *marker_var=nc.get_var("marker");
  AssertThrow(marker_var->is_valid(), ExcIO());
  AssertThrow(marker_var->num_dims()==1, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                marker_var->get_dim(0)->size())==n_markers, ExcIO());

  std::vector<int> marker(n_markers);
  // use &* to convert
  // vector<int>::iterator to int *
  marker_var->get(&*marker.begin(), n_markers);

  if (output)
    {
      std::cout << "n_cell=" << n_cells << std::endl;
      std::cout << "marker: ";
      for (unsigned int i=0; i<n_markers; ++i)
        std::cout << marker[i] << " ";
      std::cout << std::endl;
    }

  // next we read
  // int boundarymarker_of_surfaces(
  //   no_of_surfaceelements)
  NcDim *bquads_dim=nc.get_dim("no_of_surfacequadrilaterals");
  AssertThrow(bquads_dim->is_valid(), ExcIO());
  const unsigned int n_bquads=bquads_dim->size();

  NcVar *bmarker_var=nc.get_var("boundarymarker_of_surfaces");
  AssertThrow(bmarker_var->is_valid(), ExcIO());
  AssertThrow(bmarker_var->num_dims()==1, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                bmarker_var->get_dim(0)->size())==n_bquads, ExcIO());

  std::vector<int> bmarker(n_bquads);
  bmarker_var->get(&*bmarker.begin(), n_bquads);

  // for each marker count the
  // number of boundary quads
  // which carry this marker
  std::map<int, unsigned int> n_bquads_per_bmarker;
  for (unsigned int i=0; i<n_markers; ++i)
    {
      // the markers should all be
      // different
      AssertThrow(n_bquads_per_bmarker.find(marker[i])==
                  n_bquads_per_bmarker.end(), ExcIO());

      n_bquads_per_bmarker[marker[i]]=
        count(bmarker.begin(), bmarker.end(), marker[i]);
    }
  // Note: the n_bquads_per_bmarker
  // map could be used to find the
  // right coord by finding the
  // marker0 such that
  // a/ n_bquads_per_bmarker[marker0]==n_cells
  // b/ point[coord]==0,
  // Condition a/ would hold for at
  // least two markers, marker0 and
  // marker1, whereas b/ would hold
  // for marker0 only. For marker1 we
  // then had point[coord]=constant
  // with e.g. constant=1 or -1
  if (output)
    {
      std::cout << "n_bquads_per_bmarker: " << std::endl;
      std::map<int, unsigned int>::const_iterator
      iter=n_bquads_per_bmarker.begin();
      for (; iter!=n_bquads_per_bmarker.end(); ++iter)
        std::cout << "  n_bquads_per_bmarker[" << iter->first
                  << "]=" << iter->second << std::endl;
    }

  // next we read
  // int points_of_surfacequadrilaterals(
  //   no_of_surfacequadrilaterals,
  //   points_per_surfacequadrilateral)
  NcDim *quad_vertices_dim=nc.get_dim("points_per_surfacequadrilateral");
  AssertThrow(quad_vertices_dim->is_valid(), ExcIO());
  const unsigned int vertices_per_quad=quad_vertices_dim->size();
  AssertThrow(vertices_per_quad==GeometryInfo<dim>::vertices_per_cell, ExcIO());

  NcVar *vertex_indices_var=nc.get_var("points_of_surfacequadrilaterals");
  AssertThrow(vertex_indices_var->is_valid(), ExcIO());
  AssertThrow(vertex_indices_var->num_dims()==2, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                vertex_indices_var->get_dim(0)->size())==n_bquads, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                vertex_indices_var->get_dim(1)->size())==vertices_per_quad, ExcIO());

  std::vector<int> vertex_indices(n_bquads*vertices_per_quad);
  vertex_indices_var->get(&*vertex_indices.begin(), n_bquads, vertices_per_quad);

  for (unsigned int i=0; i<vertex_indices.size(); ++i)
    AssertThrow(vertex_indices[i]>=0, ExcInternalError());

  if (output)
    {
      std::cout << "vertex_indices:" << std::endl;
      for (unsigned int i=0, v=0; i<n_bquads; ++i)
        {
          for (unsigned int j=0; j<vertices_per_quad; ++j)
            std::cout << vertex_indices[v++] << " ";
          std::cout << std::endl;
        }
    }

  // next we read
  //   double points_xc(no_of_points)
  //   double points_yc(no_of_points)
  //   double points_zc(no_of_points)
  NcDim *vertices_dim=nc.get_dim("no_of_points");
  AssertThrow(vertices_dim->is_valid(), ExcIO());
  const unsigned int n_vertices=vertices_dim->size();
  if (output)
    std::cout << "n_vertices=" << n_vertices << std::endl;

  NcVar *points_xc=nc.get_var("points_xc");
  NcVar *points_yc=nc.get_var("points_yc");
  NcVar *points_zc=nc.get_var("points_zc");
  AssertThrow(points_xc->is_valid(), ExcIO());
  AssertThrow(points_yc->is_valid(), ExcIO());
  AssertThrow(points_zc->is_valid(), ExcIO());
  AssertThrow(points_xc->num_dims()==1, ExcIO());
  AssertThrow(points_yc->num_dims()==1, ExcIO());
  AssertThrow(points_zc->num_dims()==1, ExcIO());
  AssertThrow(points_yc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  AssertThrow(points_zc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  AssertThrow(points_xc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  std::vector<std::vector<double> > point_values(
    3, std::vector<double> (n_vertices));
  points_xc->get(&*point_values[0].begin(), n_vertices);
  points_yc->get(&*point_values[1].begin(), n_vertices);
  points_zc->get(&*point_values[2].begin(), n_vertices);

  // and fill the vertices
  std::vector<Point<spacedim> > vertices (n_vertices);
  for (unsigned int i=0; i<n_vertices; ++i)
    {
      vertices[i](0)=point_values[x2d][i];
      vertices[i](1)=point_values[y2d][i];
    }

  // For all boundary quads in the
  // point[coord]=0 plane add the
  // bmarker to zero_plane_markers
  std::map<int, bool> zero_plane_markers;
  for (unsigned int quad=0; quad<n_bquads; ++quad)
    {
      bool zero_plane=true;
      for (unsigned int i=0; i<vertices_per_quad; ++i)
        if (point_values[coord][vertex_indices[quad*vertices_per_quad+i]]!=0)
          {
            zero_plane=false;
            break;
          }

      if (zero_plane)
        zero_plane_markers[bmarker[quad]]=true;
    }
  unsigned int sum_of_zero_plane_cells=0;
  for (std::map<int, bool>::const_iterator iter=zero_plane_markers.begin();
       iter != zero_plane_markers.end(); ++iter)
    {
      sum_of_zero_plane_cells+=n_bquads_per_bmarker[iter->first];
      if (output)
        std::cout << "bmarker=" << iter->first << std::endl;
    }
  AssertThrow(sum_of_zero_plane_cells==n_cells, ExcIO());

  // fill cells with all quads
  // associated with
  // zero_plane_markers
  std::vector<CellData<dim> > cells(n_cells);
  for (unsigned int quad=0, cell=0; quad<n_bquads; ++quad)
    {
      bool zero_plane=false;
      for (std::map<int, bool>::const_iterator iter=zero_plane_markers.begin();
           iter != zero_plane_markers.end(); ++iter)
        if (bmarker[quad]==iter->first)
          {
            zero_plane=true;
            break;
          }

      if (zero_plane)
        {
          for (unsigned int i=0; i<vertices_per_quad; ++i)
            {
              Assert(point_values[coord][vertex_indices[
                                           quad*vertices_per_quad+i]]==0, ExcNotImplemented());
              cells[cell].vertices[i]=vertex_indices[quad*vertices_per_quad+i];
            }
          ++cell;
        }
    }

  SubCellData subcelldata;
  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
#endif
}


template <>
void GridIn<3>::read_netcdf (const std::string &filename)
{
#ifndef DEAL_II_WITH_NETCDF
  // do something with the function argument
  // to make sure it at least looks used,
  // even if it is not
  (void)filename;
  AssertThrow(false, ExcNeedsNetCDF());
#else
  const unsigned int dim=3;
  const unsigned int spacedim=3;
  const bool output=false;
  Assert (tria != 0, ExcNoTriangulationSelected());
  // this function assumes the TAU
  // grid format.

  // First, open the file
  NcFile nc (filename.c_str());
  AssertThrow(nc.is_valid(), ExcIO());

  // then read n_cells
  NcDim *elements_dim=nc.get_dim("no_of_elements");
  AssertThrow(elements_dim->is_valid(), ExcIO());
  const unsigned int n_cells=elements_dim->size();
  if (output)
    std::cout << "n_cell=" << n_cells << std::endl;
  // and n_hexes
  NcDim *hexes_dim=nc.get_dim("no_of_hexaeders");
  AssertThrow(hexes_dim->is_valid(), ExcIO());
  const unsigned int n_hexes=hexes_dim->size();
  AssertThrow(n_hexes==n_cells,
              ExcMessage("deal.II can handle purely hexaedral grids, only."));

  // next we read
  // int points_of_hexaeders(
  //   no_of_hexaeders,
  //   points_per_hexaeder)
  NcDim *hex_vertices_dim=nc.get_dim("points_per_hexaeder");
  AssertThrow(hex_vertices_dim->is_valid(), ExcIO());
  const unsigned int vertices_per_hex=hex_vertices_dim->size();
  AssertThrow(vertices_per_hex==GeometryInfo<dim>::vertices_per_cell, ExcIO());

  NcVar *vertex_indices_var=nc.get_var("points_of_hexaeders");
  AssertThrow(vertex_indices_var->is_valid(), ExcIO());
  AssertThrow(vertex_indices_var->num_dims()==2, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                vertex_indices_var->get_dim(0)->size())==n_cells, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                vertex_indices_var->get_dim(1)->size())==vertices_per_hex, ExcIO());

  std::vector<int> vertex_indices(n_cells*vertices_per_hex);
  // use &* to convert
  // vector<int>::iterator to int *
  vertex_indices_var->get(&*vertex_indices.begin(), n_cells, vertices_per_hex);

  for (unsigned int i=0; i<vertex_indices.size(); ++i)
    AssertThrow(vertex_indices[i]>=0, ExcInternalError());

  if (output)
    {
      std::cout << "vertex_indices:" << std::endl;
      for (unsigned int cell=0, v=0; cell<n_cells; ++cell)
        {
          for (unsigned int i=0; i<vertices_per_hex; ++i)
            std::cout << vertex_indices[v++] << " ";
          std::cout << std::endl;
        }
    }

  // next we read
  //   double points_xc(no_of_points)
  //   double points_yc(no_of_points)
  //   double points_zc(no_of_points)
  NcDim *vertices_dim=nc.get_dim("no_of_points");
  AssertThrow(vertices_dim->is_valid(), ExcIO());
  const unsigned int n_vertices=vertices_dim->size();
  if (output)
    std::cout << "n_vertices=" << n_vertices << std::endl;

  NcVar *points_xc=nc.get_var("points_xc");
  NcVar *points_yc=nc.get_var("points_yc");
  NcVar *points_zc=nc.get_var("points_zc");
  AssertThrow(points_xc->is_valid(), ExcIO());
  AssertThrow(points_yc->is_valid(), ExcIO());
  AssertThrow(points_zc->is_valid(), ExcIO());
  AssertThrow(points_xc->num_dims()==1, ExcIO());
  AssertThrow(points_yc->num_dims()==1, ExcIO());
  AssertThrow(points_zc->num_dims()==1, ExcIO());
  AssertThrow(points_yc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  AssertThrow(points_zc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  AssertThrow(points_xc->get_dim(0)->size()==
              static_cast<int>(n_vertices), ExcIO());
  std::vector<std::vector<double> > point_values(
    3, std::vector<double> (n_vertices));
  // we switch y and z
  const bool switch_y_z=false;
  points_xc->get(&*point_values[0].begin(), n_vertices);
  if (switch_y_z)
    {
      points_yc->get(&*point_values[2].begin(), n_vertices);
      points_zc->get(&*point_values[1].begin(), n_vertices);
    }
  else
    {
      points_yc->get(&*point_values[1].begin(), n_vertices);
      points_zc->get(&*point_values[2].begin(), n_vertices);
    }

  // and fill the vertices
  std::vector<Point<spacedim> > vertices (n_vertices);
  for (unsigned int i=0; i<n_vertices; ++i)
    {
      vertices[i](0)=point_values[0][i];
      vertices[i](1)=point_values[1][i];
      vertices[i](2)=point_values[2][i];
    }

  // and cells
  std::vector<CellData<dim> > cells(n_cells);
  for (unsigned int cell=0; cell<n_cells; ++cell)
    for (unsigned int i=0; i<vertices_per_hex; ++i)
      cells[cell].vertices[i]=vertex_indices[cell*vertices_per_hex+i];

  // for setting up the SubCellData
  // we read the vertex indices of
  // the boundary quadrilaterals and
  // their boundary markers

  // first we read
  // int points_of_surfacequadrilaterals(
  //   no_of_surfacequadrilaterals,
  //   points_per_surfacequadrilateral)
  NcDim *quad_vertices_dim=nc.get_dim("points_per_surfacequadrilateral");
  AssertThrow(quad_vertices_dim->is_valid(), ExcIO());
  const unsigned int vertices_per_quad=quad_vertices_dim->size();
  AssertThrow(vertices_per_quad==GeometryInfo<dim>::vertices_per_face, ExcIO());

  NcVar *bvertex_indices_var=nc.get_var("points_of_surfacequadrilaterals");
  AssertThrow(bvertex_indices_var->is_valid(), ExcIO());
  AssertThrow(bvertex_indices_var->num_dims()==2, ExcIO());
  const unsigned int n_bquads=bvertex_indices_var->get_dim(0)->size();
  AssertThrow(static_cast<unsigned int>(
                bvertex_indices_var->get_dim(1)->size())==
              GeometryInfo<dim>::vertices_per_face, ExcIO());

  std::vector<int> bvertex_indices(n_bquads*vertices_per_quad);
  bvertex_indices_var->get(&*bvertex_indices.begin(), n_bquads, vertices_per_quad);

  if (output)
    {
      std::cout << "bquads: ";
      for (unsigned int i=0; i<n_bquads; ++i)
        {
          for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
            std::cout << bvertex_indices[
                        i*GeometryInfo<dim>::vertices_per_face+v] << " ";
          std::cout << std::endl;
        }
    }

  // next we read
  // int boundarymarker_of_surfaces(
  //   no_of_surfaceelements)
  NcDim *bquads_dim=nc.get_dim("no_of_surfacequadrilaterals");
  AssertThrow(bquads_dim->is_valid(), ExcIO());
  AssertThrow(static_cast<unsigned int>(
                bquads_dim->size())==n_bquads, ExcIO());

  NcVar *bmarker_var=nc.get_var("boundarymarker_of_surfaces");
  AssertThrow(bmarker_var->is_valid(), ExcIO());
  AssertThrow(bmarker_var->num_dims()==1, ExcIO());
  AssertThrow(static_cast<unsigned int>(
                bmarker_var->get_dim(0)->size())==n_bquads, ExcIO());

  std::vector<int> bmarker(n_bquads);
  bmarker_var->get(&*bmarker.begin(), n_bquads);
  // we only handle boundary
  // indicators that fit into an
  // types::boundary_id. Also, we don't
  // take numbers::internal_face_boundary_id
  // as it denotes an internal face
  for (unsigned int i=0; i<bmarker.size(); ++i)
    Assert(0<=bmarker[i] && bmarker[i]<numbers::internal_face_boundary_id, ExcIO());

  // finally we setup the boundary
  // information
  SubCellData subcelldata;
  subcelldata.boundary_quads.resize(n_bquads);
  for (unsigned int i=0; i<n_bquads; ++i)
    {
      for (unsigned int v=0; v<GeometryInfo<dim>::vertices_per_face; ++v)
        subcelldata.boundary_quads[i].vertices[v]=bvertex_indices[
                                                    i*GeometryInfo<dim>::vertices_per_face+v];
      subcelldata.boundary_quads[i].boundary_id
        = static_cast<types::boundary_id>(bmarker[i]);
    }

  GridTools::delete_unused_vertices(vertices, cells, subcelldata);
  GridReordering<dim,spacedim>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
#endif
}


template <int dim, int spacedim>
void GridIn<dim, spacedim>::parse_tecplot_header(std::string &header,
                                                 std::vector<unsigned int> &tecplot2deal,
                                                 unsigned int &n_vars,
                                                 unsigned int &n_vertices,
                                                 unsigned int &n_cells,
                                                 std::vector<unsigned int> &IJK,
                                                 bool &structured,
                                                 bool &blocked)
{
  Assert(tecplot2deal.size()==dim, ExcInternalError());
  Assert(IJK.size()==dim, ExcInternalError());
  // initialize the output variables
  n_vars=0;
  n_vertices=0;
  n_cells=0;
  switch (dim)
    {
    case 3:
      IJK[2]=0;
    case 2:
      IJK[1]=0;
    case 1:
      IJK[0]=0;
    }
  structured=true;
  blocked=false;

  // convert the string to upper case
  std::transform(header.begin(),header.end(),header.begin(),::toupper);

  // replace all tabs, commas, newlines by
  // whitespaces
  std::replace(header.begin(),header.end(),'\t',' ');
  std::replace(header.begin(),header.end(),',',' ');
  std::replace(header.begin(),header.end(),'\n',' ');

  // now remove whitespace in front of and
  // after '='
  std::string::size_type pos=header.find("=");

  while (pos!=static_cast<std::string::size_type>(std::string::npos))
    if (header[pos+1]==' ')
      header.erase(pos+1,1);
    else if (header[pos-1]==' ')
      {
        header.erase(pos-1,1);
        --pos;
      }
    else
      pos=header.find("=",++pos);

  // split the string into individual entries
  std::vector<std::string> entries=Utilities::break_text_into_lines(header,1,' ');

  // now go through the list and try to extract
  for (unsigned int i=0; i<entries.size(); ++i)
    {
      if (Utilities::match_at_string_start(entries[i],"VARIABLES=\""))
        {
          ++n_vars;
          // we assume, that the first variable
          // is x or no coordinate at all (not y or z)
          if (Utilities::match_at_string_start(entries[i],"VARIABLES=\"X\""))
            {
              tecplot2deal[0]=0;
            }
          ++i;
          while (entries[i][0]=='"')
            {
              if (entries[i]=="\"X\"")
                tecplot2deal[0]=n_vars;
              else if (entries[i]=="\"Y\"")
                {
                  // we assume, that y contains
                  // zero data in 1d, so do
                  // nothing
                  if (dim>1)
                    tecplot2deal[1]=n_vars;
                }
              else if (entries[i]=="\"Z\"")
                {
                  // we assume, that z contains
                  // zero data in 1d and 2d, so
                  // do nothing
                  if (dim>2)
                    tecplot2deal[2]=n_vars;
                }
              ++n_vars;
              ++i;
            }
          // set i back, so that the next
          // string is treated correctly
          --i;

          AssertThrow(n_vars>=dim,
                      ExcMessage("Tecplot file must contain at least one variable for each dimension"));
          for (unsigned int d=1; d<dim; ++d)
            AssertThrow(tecplot2deal[d]>0,
                        ExcMessage("Tecplot file must contain at least one variable for each dimension."));
        }
      else if (Utilities::match_at_string_start(entries[i],"ZONETYPE=ORDERED"))
        structured=true;
      else if (Utilities::match_at_string_start(entries[i],"ZONETYPE=FELINESEG") && dim==1)
        structured=false;
      else if (Utilities::match_at_string_start(entries[i],"ZONETYPE=FEQUADRILATERAL") && dim==2)
        structured=false;
      else if (Utilities::match_at_string_start(entries[i],"ZONETYPE=FEBRICK") && dim==3)
        structured=false;
      else if (Utilities::match_at_string_start(entries[i],"ZONETYPE="))
        // unsupported ZONETYPE
        {
          AssertThrow(false,ExcMessage("The tecplot file contains an unsupported ZONETYPE."));
        }
      else if (Utilities::match_at_string_start(entries[i],"DATAPACKING=POINT"))
        blocked=false;
      else if (Utilities::match_at_string_start(entries[i],"DATAPACKING=BLOCK"))
        blocked=true;
      else if (Utilities::match_at_string_start(entries[i],"F=POINT"))
        {
          structured=true;
          blocked=false;
        }
      else if (Utilities::match_at_string_start(entries[i],"F=BLOCK"))
        {
          structured=true;
          blocked=true;
        }
      else if (Utilities::match_at_string_start(entries[i],"F=FEPOINT"))
        {
          structured=false;
          blocked=false;
        }
      else if (Utilities::match_at_string_start(entries[i],"F=FEBLOCK"))
        {
          structured=false;
          blocked=true;
        }
      else if (Utilities::match_at_string_start(entries[i],"ET=QUADRILATERAL") && dim==2)
        structured=false;
      else if (Utilities::match_at_string_start(entries[i],"ET=BRICK") && dim==3)
        structured=false;
      else if (Utilities::match_at_string_start(entries[i],"ET="))
        // unsupported ElementType
        {
          AssertThrow(false,ExcMessage("The tecplot file contains an unsupported ElementType."));
        }
      else if (Utilities::match_at_string_start(entries[i],"I="))
        IJK[0]=Utilities::get_integer_at_position(entries[i],2).first;
      else if (Utilities::match_at_string_start(entries[i],"J="))
        {
          IJK[1]=Utilities::get_integer_at_position(entries[i],2).first;
          AssertThrow(dim>1 || IJK[1]==1,
                      ExcMessage("Parameter 'J=' found in tecplot, although this is only possible for dimensions greater than 1."));
        }
      else if (Utilities::match_at_string_start(entries[i],"K="))
        {
          IJK[2]=Utilities::get_integer_at_position(entries[i],2).first;
          AssertThrow(dim>2 || IJK[2]==1,
                      ExcMessage("Parameter 'K=' found in tecplot, although this is only possible for dimensions greater than 2."));
        }
      else if (Utilities::match_at_string_start(entries[i],"N="))
        n_vertices=Utilities::get_integer_at_position(entries[i],2).first;
      else if (Utilities::match_at_string_start(entries[i],"E="))
        n_cells=Utilities::get_integer_at_position(entries[i],2).first;
    }

  // now we have read all the fields we are
  // interested in. do some checks and
  // calculate the variables
  if (structured)
    {
      n_vertices=1;
      n_cells=1;
      for (unsigned int d=0; d<dim; ++d)
        {
          AssertThrow(IJK[d]>0,
                      ExcMessage("Tecplot file does not contain a complete and consistent set of parameters"));
          n_vertices*=IJK[d];
          n_cells*=(IJK[d]-1);
        }
    }
  else
    {
      AssertThrow(n_vertices>0,
                  ExcMessage("Tecplot file does not contain a complete and consistent set of parameters"));
      if (n_cells==0)
        // this means an error, although
        // tecplot itself accepts entries like
        // 'J=20' instead of 'E=20'. therefore,
        // take the max of IJK
        n_cells=*std::max_element(IJK.begin(),IJK.end());
      AssertThrow(n_cells>0,
                  ExcMessage("Tecplot file does not contain a complete and consistent set of parameters"));
    }
}




template <>
void GridIn<2>::read_tecplot (std::istream &in)
{
  const unsigned int dim=2;
  const unsigned int spacedim=2;
  Assert (tria != 0, ExcNoTriangulationSelected());
  AssertThrow (in, ExcIO());

  // skip comments at start of file
  skip_comment_lines (in, '#');

  // some strings for parsing the header
  std::string line, header;

  // first, concatenate all header lines
  // create a searchstring with almost all
  // letters. exclude e and E from the letters
  // to search, as they might appear in
  // exponential notation
  std::string letters ="abcdfghijklmnopqrstuvwxyzABCDFGHIJKLMNOPQRSTUVWXYZ";

  getline(in,line);
  while (line.find_first_of(letters)!=std::string::npos)
    {
      header+=" "+line;
      getline(in,line);
    }

  // now create some variables holding
  // important information on the mesh, get
  // this information from the header string
  std::vector<unsigned int> tecplot2deal(dim);
  std::vector<unsigned int> IJK(dim);
  unsigned int n_vars,
           n_vertices,
           n_cells;
  bool structured,
       blocked;

  parse_tecplot_header(header,
                       tecplot2deal,n_vars,n_vertices,n_cells,IJK,
                       structured,blocked);

  // reserve space for vertices. note, that in
  // tecplot vertices are ordered beginning
  // with 1, whereas in deal all indices start
  // with 0. in order not to use -1 for all the
  // connectivity information, a 0th vertex
  // (unused) is inserted at the origin.
  std::vector<Point<spacedim> > vertices(n_vertices+1);
  vertices[0]=Point<spacedim>();
  // reserve space for cells
  std::vector<CellData<dim> > cells(n_cells);
  SubCellData                 subcelldata;

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
      unsigned int next_index=0;

      // note, that we have already read the
      // first line containing the first variable
      if (tecplot2deal[0]==0)
        {
          // we need the information in this
          // line, so extract it
          std::vector<std::string> first_var=Utilities::break_text_into_lines(line,1);
          char *endptr;
          for (unsigned int i=1; i<first_var.size()+1; ++i)
            vertices[i](0) = std::strtod (first_var[i-1].c_str(), &endptr);

          // if there are many points, the data
          // for this var might continue in the
          // next line(s)
          for (unsigned int j=first_var.size()+1; j<n_vertices+1; ++j)
            in>>vertices[j](next_index);
          // now we got all values of the first
          // variable, so increase the counter
          next_index=1;
        }

      // main loop over all variables
      for (unsigned int i=1; i<n_vars; ++i)
        {
          // if we read all the important
          // variables and do not want to
          // read further, because we are
          // using a structured grid, we can
          // stop here (and skip, for
          // example, a whole lot of solution
          // variables)
          if (next_index==dim && structured)
            break;

          if ((next_index<dim) && (i==tecplot2deal[next_index]))
            {
              // we need this line, read it in
              for (unsigned int j=1; j<n_vertices+1; ++j)
                in>>vertices[j](next_index);
              ++next_index;
            }
          else
            {
              // we do not need this line, read
              // it in and discard it
              for (unsigned int j=1; j<n_vertices+1; ++j)
                in>>dummy;
            }
        }
      Assert(next_index==dim, ExcInternalError());
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
      std::vector<std::string> first_vertex=Utilities::break_text_into_lines(line,1);
      char *endptr;
      for (unsigned int d=0; d<dim; ++d)
        vertices[1](d) = std::strtod (first_vertex[tecplot2deal[d]].c_str(), &endptr);

      // read the remaining vertices from the
      // list
      for (unsigned int v=2; v<n_vertices+1; ++v)
        {
          for (unsigned int i=0; i<n_vars; ++i)
            in>>vars[i];
          // fill the vertex
          // coordinates. respect the position
          // of coordinates in the list of
          // variables
          for (unsigned int i=0; i<dim; ++i)
            vertices[v](i)=vars[tecplot2deal[i]];
        }
    }

  if (structured)
    {
      // this is the part of the code that only
      // works in 2d
      unsigned int I=IJK[0],
                   J=IJK[1];

      unsigned int cell=0;
      // set up array of cells
      for (unsigned int j=0; j<J-1; ++j)
        for (unsigned int i=1; i<I; ++i)
          {
            cells[cell].vertices[0]=i+  j    *I;
            cells[cell].vertices[1]=i+1+j    *I;
            cells[cell].vertices[2]=i+1+(j+1)*I;
            cells[cell].vertices[3]=i  +(j+1)*I;
            ++cell;
          }
      Assert(cell=n_cells, ExcInternalError());
      std::vector<unsigned int> boundary_vertices(2*I+2*J-4);
      unsigned int k=0;
      for (unsigned int i=1; i<I+1; ++i)
        {
          boundary_vertices[k]=i;
          ++k;
          boundary_vertices[k]=i+(J-1)*I;
          ++k;
        }
      for (unsigned int j=1; j<J-1; ++j)
        {
          boundary_vertices[k]=1+j*I;
          ++k;
          boundary_vertices[k]=I+j*I;
          ++k;
        }
      Assert(k==boundary_vertices.size(), ExcInternalError());
      // delete the duplicated vertices at the
      // boundary, which occur, e.g. in c-type
      // or o-type grids around a body
      // (airfoil). this automatically deletes
      // unused vertices as well.
      GridTools::delete_duplicated_vertices(vertices,cells,subcelldata,boundary_vertices);
    }
  else
    {
      // set up array of cells, unstructured
      // mode, so the connectivity is
      // explicitly given
      for (unsigned int i=0; i<n_cells; ++i)
        {
          // note that since in the input file
          // we found the number of cells at
          // the top, there should still be
          // input here, so check this:
          AssertThrow (in, ExcIO());

          // get the connectivity from the
          // input file. the vertices are
          // ordered like in the ucd format
          for (unsigned int j=0; j<GeometryInfo<dim>::vertices_per_cell; ++j)
            in>>cells[i].vertices[j];
        }
      // do some clean-up on vertices
      GridTools::delete_unused_vertices (vertices, cells, subcelldata);
    }

  // check that no forbidden arrays are
  // used. as we do not read in any
  // subcelldata, nothing should happen here.
  Assert (subcelldata.check_consistency(dim), ExcInternalError());
  AssertThrow (in, ExcIO());

  // do some cleanup on cells
  GridReordering<dim,spacedim>::invert_all_cells_of_negative_grid (vertices, cells);
  GridReordering<dim,spacedim>::reorder_cells (cells);
  tria->create_triangulation_compatibility (vertices, cells, subcelldata);
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::read_tecplot(std::istream &)
{
  Assert(false, ExcNotImplemented());
}


template <int dim, int spacedim>
void GridIn<dim, spacedim>::skip_empty_lines (std::istream &in)
{
  std::string line;
  while (in)
    {
      // get line
      getline (in, line);

      // check if this is a line that
      // consists only of spaces, and
      // if not put the whole thing
      // back and return
      if (std::find_if (line.begin(), line.end(),
                        std::bind2nd (std::not_equal_to<char>(),' '))
          != line.end())
        {
          in.putback ('\n');
          for (int i=line.length()-1; i>=0; --i)
            in.putback (line[i]);
          return;
        }

      // else: go on with next line
    }
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::skip_comment_lines (std::istream &in,
                                                const char    comment_start)
{
  char c;
  // loop over the following comment
  // lines
  while ((c=in.get()) == comment_start)
    // loop over the characters after
    // the comment starter
    while (in.get() != '\n')
      ;


  // put back first character of
  // first non-comment line
  in.putback (c);

  // at last: skip additional empty lines, if present
  skip_empty_lines(in);
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::debug_output_grid (const std::vector<CellData<dim> > &/*cells*/,
                                               const std::vector<Point<spacedim> >    &/*vertices*/,
                                               std::ostream                      &/*out*/)
{
  Assert (false, ExcNotImplemented());
}



template <>
void
GridIn<2>::debug_output_grid (const std::vector<CellData<2> > &cells,
                              const std::vector<Point<2> >    &vertices,
                              std::ostream                    &out)
{
  double min_x = vertices[cells[0].vertices[0]](0),
         max_x = vertices[cells[0].vertices[0]](0),
         min_y = vertices[cells[0].vertices[0]](1),
         max_y = vertices[cells[0].vertices[0]](1);

  for (unsigned int i=0; i<cells.size(); ++i)
    {
      for (unsigned int v=0; v<4; ++v)
        {
          const Point<2> &p = vertices[cells[i].vertices[v]];

          if (p(0) < min_x)
            min_x = p(0);
          if (p(0) > max_x)
            max_x = p(0);
          if (p(1) < min_y)
            min_y = p(1);
          if (p(1) > max_y)
            max_y = p(1);
        };

      out << "# cell " << i << std::endl;
      Point<2> center;
      for (unsigned int f=0; f<4; ++f)
        center += vertices[cells[i].vertices[f]];
      center /= 4;

      out << "set label \"" << i << "\" at "
          << center(0) << ',' << center(1)
          << " center"
          << std::endl;

      // first two line right direction
      for (unsigned int f=0; f<2; ++f)
        out << "set arrow from "
            << vertices[cells[i].vertices[f]](0) << ','
            << vertices[cells[i].vertices[f]](1)
            << " to "
            << vertices[cells[i].vertices[(f+1)%4]](0) << ','
            << vertices[cells[i].vertices[(f+1)%4]](1)
            << std::endl;
      // other two lines reverse direction
      for (unsigned int f=2; f<4; ++f)
        out << "set arrow from "
            << vertices[cells[i].vertices[(f+1)%4]](0) << ','
            << vertices[cells[i].vertices[(f+1)%4]](1)
            << " to "
            << vertices[cells[i].vertices[f]](0) << ','
            << vertices[cells[i].vertices[f]](1)
            << std::endl;
      out << std::endl;
    };


  out << std::endl
      << "set nokey" << std::endl
      << "pl [" << min_x << ':' << max_x << "]["
      << min_y << ':' << max_y <<  "] "
      << min_y << std::endl
      << "pause -1" << std::endl;
}



template <>
void
GridIn<3>::debug_output_grid (const std::vector<CellData<3> > &cells,
                              const std::vector<Point<3> >    &vertices,
                              std::ostream                    &out)
{
  for (unsigned int cell=0; cell<cells.size(); ++cell)
    {
      // line 0
      out << vertices[cells[cell].vertices[0]]
          << std::endl
          << vertices[cells[cell].vertices[1]]
          << std::endl << std::endl << std::endl;
      // line 1
      out << vertices[cells[cell].vertices[1]]
          << std::endl
          << vertices[cells[cell].vertices[2]]
          << std::endl << std::endl << std::endl;
      // line 2
      out << vertices[cells[cell].vertices[3]]
          << std::endl
          << vertices[cells[cell].vertices[2]]
          << std::endl << std::endl << std::endl;
      // line 3
      out << vertices[cells[cell].vertices[0]]
          << std::endl
          << vertices[cells[cell].vertices[3]]
          << std::endl << std::endl << std::endl;
      // line 4
      out << vertices[cells[cell].vertices[4]]
          << std::endl
          << vertices[cells[cell].vertices[5]]
          << std::endl << std::endl << std::endl;
      // line 5
      out << vertices[cells[cell].vertices[5]]
          << std::endl
          << vertices[cells[cell].vertices[6]]
          << std::endl << std::endl << std::endl;
      // line 6
      out << vertices[cells[cell].vertices[7]]
          << std::endl
          << vertices[cells[cell].vertices[6]]
          << std::endl << std::endl << std::endl;
      // line 7
      out << vertices[cells[cell].vertices[4]]
          << std::endl
          << vertices[cells[cell].vertices[7]]
          << std::endl << std::endl << std::endl;
      // line 8
      out << vertices[cells[cell].vertices[0]]
          << std::endl
          << vertices[cells[cell].vertices[4]]
          << std::endl << std::endl << std::endl;
      // line 9
      out << vertices[cells[cell].vertices[1]]
          << std::endl
          << vertices[cells[cell].vertices[5]]
          << std::endl << std::endl << std::endl;
      // line 10
      out << vertices[cells[cell].vertices[2]]
          << std::endl
          << vertices[cells[cell].vertices[6]]
          << std::endl << std::endl << std::endl;
      // line 11
      out << vertices[cells[cell].vertices[3]]
          << std::endl
          << vertices[cells[cell].vertices[7]]
          << std::endl << std::endl << std::endl;
    };
}



template <int dim, int spacedim>
void GridIn<dim, spacedim>::read (const std::string &filename,
                                  Format format)
{
  // Search file class for meshes
  PathSearch search("MESH");
  std::string name;
  // Open the file and remember its name
  if (format == Default)
    name = search.find(filename);
  else
    name = search.find(filename, default_suffix(format));

  std::ifstream in(name.c_str());

  if (format == Default)
    {
      const std::string::size_type slashpos = name.find_last_of('/');
      const std::string::size_type dotpos = name.find_last_of('.');
      if (dotpos < name.length()
          && (dotpos > slashpos || slashpos == std::string::npos))
        {
          std::string ext = name.substr(dotpos+1);
          format = parse_format(ext);
        }
    }
  if (format == netcdf)
    read_netcdf(filename);
  else
    read(in, format);
}


template <int dim, int spacedim>
void GridIn<dim, spacedim>::read (std::istream &in,
                                  Format format)
{
  if (format == Default)
    format = default_format;

  switch (format)
    {
    case dbmesh:
      read_dbmesh (in);
      return;

    case msh:
      read_msh (in);
      return;

    case vtk:
      read_vtk (in);
      return;

    case unv:
      read_unv (in);
      return;

    case ucd:
      read_ucd (in);
      return;

    case xda:
      read_xda (in);
      return;

    case netcdf:
      Assert(false, ExcMessage("There is no read_netcdf(istream &) function. "
                               "Use the read(_netcdf)(string &filename) "
                               "functions, instead."));
      return;

    case tecplot:
      read_tecplot (in);
      return;

    case Default:
      break;
    }
  Assert (false, ExcInternalError());
}



template <int dim, int spacedim>
std::string
GridIn<dim, spacedim>::default_suffix (const Format format)
{
  switch (format)
    {
    case dbmesh:
      return ".dbmesh";
    case msh:
      return ".msh";
    case vtk:
      return ".vtk";
    case unv:
      return ".unv";
    case ucd:
      return ".inp";
    case xda:
      return ".xda";
    case netcdf:
      return ".nc";
    case tecplot:
      return ".dat";
    default:
      Assert (false, ExcNotImplemented());
      return ".unknown_format";
    }
}



template <int dim, int spacedim>
typename GridIn<dim, spacedim>::Format
GridIn<dim, spacedim>::parse_format (const std::string &format_name)
{
  if (format_name == "dbmesh")
    return dbmesh;

  if (format_name == "msh")
    return msh;

  if (format_name == "unv")
    return unv;

  if (format_name == "vtk")
    return vtk;

  if (format_name == "inp")
    return ucd;

  if (format_name == "ucd")
    return ucd;

  if (format_name == "xda")
    return xda;

  if (format_name == "netcdf")
    return netcdf;

  if (format_name == "nc")
    return netcdf;

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

  AssertThrow (false, ExcInvalidState ());
  // return something weird
  return Format(Default);
}



template <int dim, int spacedim>
std::string GridIn<dim, spacedim>::get_format_names ()
{
  return "dbmesh|msh|unv|vtk|ucd|xda|netcdf|tecplot";
}



//explicit instantiations
#include "grid_in.inst"

DEAL_II_NAMESPACE_CLOSE
