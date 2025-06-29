// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1999 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/base/exceptions.h>
#include <deal.II/base/geometry_info.h>
#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/point.h>
#include <deal.II/base/qprojector.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/distributed/tria_base.h>

#include <deal.II/fe/mapping.h>

#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/numerics/data_out.h>

#include <boost/archive/binary_oarchive.hpp>

#ifdef DEAL_II_GMSH_WITH_API
#  include <gmsh.h>
#endif

#include <algorithm>
#include <cmath>
#include <cstring>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <list>
#include <set>

DEAL_II_NAMESPACE_OPEN


namespace GridOutFlags
{
  DX::DX(const bool write_cells,
         const bool write_faces,
         const bool write_diameter,
         const bool write_measure,
         const bool write_all_faces)
    : write_cells(write_cells)
    , write_faces(write_faces)
    , write_diameter(write_diameter)
    , write_measure(write_measure)
    , write_all_faces(write_all_faces)
  {}

  void
  DX::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Write cells",
                        "true",
                        Patterns::Bool(),
                        "Write the mesh connectivity as DX grid cells");
    param.declare_entry("Write faces",
                        "false",
                        Patterns::Bool(),
                        "Write faces of cells. These may be boundary faces "
                        "or all faces between mesh cells, according to "
                        "\"Write all faces\"");
    param.declare_entry("Write diameter",
                        "false",
                        Patterns::Bool(),
                        "If cells are written, additionally write their"
                        " diameter as data for visualization");
    param.declare_entry("Write measure",
                        "false",
                        Patterns::Bool(),
                        "Write the volume of each cell as data");
    param.declare_entry("Write all faces",
                        "true",
                        Patterns::Bool(),
                        "Write all faces, not only boundary");
  }

  void
  DX::parse_parameters(ParameterHandler &param)
  {
    write_cells     = param.get_bool("Write cells");
    write_faces     = param.get_bool("Write faces");
    write_diameter  = param.get_bool("Write diameter");
    write_measure   = param.get_bool("Write measure");
    write_all_faces = param.get_bool("Write all faces");
  }


  Msh::Msh(const bool write_faces, const bool write_lines)
    : write_faces(write_faces)
    , write_lines(write_lines)
  {}

  void
  Msh::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Write faces", "false", Patterns::Bool());
    param.declare_entry("Write lines", "false", Patterns::Bool());
  }


  void
  Msh::parse_parameters(ParameterHandler &param)
  {
    write_faces = param.get_bool("Write faces");
    write_lines = param.get_bool("Write lines");
  }


  Ucd::Ucd(const bool write_preamble,
           const bool write_faces,
           const bool write_lines)
    : write_preamble(write_preamble)
    , write_faces(write_faces)
    , write_lines(write_lines)
  {}



  void
  Ucd::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Write preamble", "true", Patterns::Bool());
    param.declare_entry("Write faces", "false", Patterns::Bool());
    param.declare_entry("Write lines", "false", Patterns::Bool());
  }


  void
  Ucd::parse_parameters(ParameterHandler &param)
  {
    write_preamble = param.get_bool("Write preamble");
    write_faces    = param.get_bool("Write faces");
    write_lines    = param.get_bool("Write lines");
  }



  Gnuplot::Gnuplot(const bool         write_cell_numbers,
                   const unsigned int n_extra_curved_line_points,
                   const bool         curved_inner_cells,
                   const bool         write_additional_boundary_lines)
    : write_cell_numbers(write_cell_numbers)
    , n_extra_curved_line_points(n_extra_curved_line_points)
    , curved_inner_cells(curved_inner_cells)
    , write_additional_boundary_lines(write_additional_boundary_lines)
  {}



  void
  Gnuplot::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Cell number", "false", Patterns::Bool());
    param.declare_entry("Boundary points", "2", Patterns::Integer());
  }


  void
  Gnuplot::parse_parameters(ParameterHandler &param)
  {
    write_cell_numbers         = param.get_bool("Cell number");
    n_extra_curved_line_points = param.get_integer("Boundary points");
  }


  EpsFlagsBase::EpsFlagsBase(const SizeType     size_type,
                             const unsigned int size,
                             const double       line_width,
                             const bool         color_lines_on_user_flag,
                             const unsigned int n_boundary_face_points,
                             const bool         color_lines_level)
    : size_type(size_type)
    , size(size)
    , line_width(line_width)
    , color_lines_on_user_flag(color_lines_on_user_flag)
    , n_boundary_face_points(n_boundary_face_points)
    , color_lines_level(color_lines_level)
  {}


  void
  EpsFlagsBase::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Size by",
                        "width",
                        Patterns::Selection("width|height"),
                        "Depending on this parameter, either the "
                        "width or height "
                        "of the eps is scaled to \"Size\"");
    param.declare_entry("Size",
                        "300",
                        Patterns::Integer(),
                        "Size of the output in points");
    param.declare_entry("Line width",
                        "0.5",
                        Patterns::Double(),
                        "Width of the lines drawn in points");
    param.declare_entry("Color by flag",
                        "false",
                        Patterns::Bool(),
                        "Draw lines with user flag set in different color");
    param.declare_entry("Boundary points",
                        "2",
                        Patterns::Integer(),
                        "Number of points on boundary edges. "
                        "Increase this beyond 2 to see curved boundaries.");
    param.declare_entry("Color by level",
                        "false",
                        Patterns::Bool(),
                        "Draw different colors according to grid level.");
  }


  void
  EpsFlagsBase::parse_parameters(ParameterHandler &param)
  {
    if (param.get("Size by") == "width")
      size_type = width;
    else if (param.get("Size by") == "height")
      size_type = height;
    size                     = param.get_integer("Size");
    line_width               = param.get_double("Line width");
    color_lines_on_user_flag = param.get_bool("Color by flag");
    n_boundary_face_points   = param.get_integer("Boundary points");
    color_lines_level        = param.get_bool("Color by level");
  }



  Eps<1>::Eps(const SizeType     size_type,
              const unsigned int size,
              const double       line_width,
              const bool         color_lines_on_user_flag,
              const unsigned int n_boundary_face_points)
    : EpsFlagsBase(size_type,
                   size,
                   line_width,
                   color_lines_on_user_flag,
                   n_boundary_face_points)
  {}


  void
  Eps<1>::declare_parameters(ParameterHandler &)
  {}


  void
  Eps<1>::parse_parameters(ParameterHandler &param)
  {
    EpsFlagsBase::parse_parameters(param);
  }



  Eps<2>::Eps(const SizeType     size_type,
              const unsigned int size,
              const double       line_width,
              const bool         color_lines_on_user_flag,
              const unsigned int n_boundary_face_points,
              const bool         write_cell_numbers,
              const bool         write_cell_number_level,
              const bool         write_vertex_numbers,
              const bool         color_lines_level)
    : EpsFlagsBase(size_type,
                   size,
                   line_width,
                   color_lines_on_user_flag,
                   n_boundary_face_points,
                   color_lines_level)
    , write_cell_numbers(write_cell_numbers)
    , write_cell_number_level(write_cell_number_level)
    , write_vertex_numbers(write_vertex_numbers)
  {}


  void
  Eps<2>::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Cell number",
                        "false",
                        Patterns::Bool(),
                        "(2d only) Write cell numbers"
                        " into the centers of cells");
    param.declare_entry("Level number",
                        "false",
                        Patterns::Bool(),
                        "(2d only) if \"Cell number\" is true, write "
                        "numbers in the form level.number");
    param.declare_entry("Vertex number",
                        "false",
                        Patterns::Bool(),
                        "Write numbers for each vertex");
  }


  void
  Eps<2>::parse_parameters(ParameterHandler &param)
  {
    EpsFlagsBase::parse_parameters(param);
    write_cell_numbers      = param.get_bool("Cell number");
    write_cell_number_level = param.get_bool("Level number");
    write_vertex_numbers    = param.get_bool("Vertex number");
  }



  Eps<3>::Eps(const SizeType     size_type,
              const unsigned int size,
              const double       line_width,
              const bool         color_lines_on_user_flag,
              const unsigned int n_boundary_face_points,
              const double       azimut_angle,
              const double       turn_angle)
    : EpsFlagsBase(size_type,
                   size,
                   line_width,
                   color_lines_on_user_flag,
                   n_boundary_face_points)
    , azimut_angle(azimut_angle)
    , turn_angle(turn_angle)
  {}


  void
  Eps<3>::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Azimuth",
                        "30",
                        Patterns::Double(),
                        "Azimuth of the view point, that is, the angle "
                        "in the plane from the x-axis.");
    param.declare_entry("Elevation",
                        "30",
                        Patterns::Double(),
                        "Elevation of the view point above the xy-plane.");
  }


  void
  Eps<3>::parse_parameters(ParameterHandler &param)
  {
    EpsFlagsBase::parse_parameters(param);
    azimut_angle = 90 - param.get_double("Elevation");
    turn_angle   = param.get_double("Azimuth");
  }



  XFig::XFig()
    : draw_boundary(true)
    , color_by(material_id)
    , level_depth(true)
    , n_boundary_face_points(0)
    , scaling(1., 1.)
    , fill_style(20)
    , line_style(0)
    , line_thickness(1)
    , boundary_style(0)
    , boundary_thickness(3)
  {}


  void
  XFig::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Boundary", "true", Patterns::Bool());
    param.declare_entry("Level color", "false", Patterns::Bool());
    param.declare_entry("Level depth", "true", Patterns::Bool());
    // TODO: Unify this number with other output formats
    param.declare_entry("Boundary points", "0", Patterns::Integer());
    param.declare_entry("Fill style", "20", Patterns::Integer());
    param.declare_entry("Line style", "0", Patterns::Integer());
    param.declare_entry("Line width", "1", Patterns::Integer());
    param.declare_entry("Boundary style", "0", Patterns::Integer());
    param.declare_entry("Boundary width", "3", Patterns::Integer());
  }


  void
  XFig::parse_parameters(ParameterHandler &param)
  {
    draw_boundary          = param.get_bool("Boundary");
    level_depth            = param.get_bool("Level depth");
    n_boundary_face_points = param.get_integer("Boundary points");
    fill_style             = param.get_integer("Fill style");
    line_style             = param.get_integer("Line style");
    line_thickness         = param.get_integer("Line width");
    boundary_style         = param.get_integer("Boundary style");
    boundary_thickness     = param.get_integer("Boundary width");
  }

  Svg::Svg(const unsigned int line_thickness,
           const unsigned int boundary_line_thickness,
           bool               margin,
           const Background   background,
           const int          azimuth_angle,
           const int          polar_angle,
           const Coloring     coloring,
           const bool         convert_level_number_to_height,
           const bool         label_level_number,
           const bool         label_cell_index,
           const bool         label_material_id,
           const bool         label_subdomain_id,
           const bool         draw_colorbar,
           const bool         draw_legend,
           const bool         label_boundary_id)
    : height(1000)
    , width(0)
    , line_thickness(line_thickness)
    , boundary_line_thickness(boundary_line_thickness)
    , margin(margin)
    , background(background)
    , azimuth_angle(azimuth_angle)
    , polar_angle(polar_angle)
    , coloring(coloring)
    , convert_level_number_to_height(convert_level_number_to_height)
    , level_height_factor(0.3f)
    , cell_font_scaling(1.f)
    , label_level_number(label_level_number)
    , label_cell_index(label_cell_index)
    , label_material_id(label_material_id)
    , label_subdomain_id(label_subdomain_id)
    , label_level_subdomain_id(false)
    , label_boundary_id(label_boundary_id)
    , draw_colorbar(draw_colorbar)
    , draw_legend(draw_legend)
  {}

  MathGL::MathGL()
    : draw_bounding_box(false) // box
  {}

  void
  MathGL::declare_parameters(ParameterHandler &param)
  {
    param.declare_entry("Draw bounding box", "false", Patterns::Bool());
  }

  void
  MathGL::parse_parameters(ParameterHandler &param)
  {
    draw_bounding_box = param.get_bool("Draw bounding box");
  }
} // end namespace GridOutFlags



GridOut::GridOut()
  : default_format(none)
{}


void
GridOut::set_flags(const GridOutFlags::DX &flags)
{
  dx_flags = flags;
}



void
GridOut::set_flags(const GridOutFlags::Msh &flags)
{
  msh_flags = flags;
}


void
GridOut::set_flags(const GridOutFlags::Ucd &flags)
{
  ucd_flags = flags;
}



void
GridOut::set_flags(const GridOutFlags::Gnuplot &flags)
{
  gnuplot_flags = flags;
}



void
GridOut::set_flags(const GridOutFlags::Eps<1> &flags)
{
  eps_flags_1 = flags;
}



void
GridOut::set_flags(const GridOutFlags::Eps<2> &flags)
{
  eps_flags_2 = flags;
}



void
GridOut::set_flags(const GridOutFlags::Eps<3> &flags)
{
  eps_flags_3 = flags;
}



void
GridOut::set_flags(const GridOutFlags::XFig &flags)
{
  xfig_flags = flags;
}


void
GridOut::set_flags(const GridOutFlags::Svg &flags)
{
  svg_flags = flags;
}


void
GridOut::set_flags(const GridOutFlags::MathGL &flags)
{
  mathgl_flags = flags;
}

void
GridOut::set_flags(const GridOutFlags::Vtk &flags)
{
  vtk_flags = flags;
}

void
GridOut::set_flags(const GridOutFlags::Vtu &flags)
{
  vtu_flags = flags;
}

std::string
GridOut::default_suffix(const OutputFormat output_format)
{
  switch (output_format)
    {
      case none:
        return "";
      case dx:
        return ".dx";
      case gnuplot:
        return ".gnuplot";
      case ucd:
        return ".inp";
      case eps:
        return ".eps";
      case xfig:
        return ".fig";
      case msh:
        return ".msh";
      case svg:
        return ".svg";
      case mathgl:
        return ".mathgl";
      case vtk:
        return ".vtk";
      case vtu:
        return ".vtu";
      default:
        DEAL_II_NOT_IMPLEMENTED();
        return "";
    }
}



std::string
GridOut::default_suffix() const
{
  return default_suffix(default_format);
}



GridOut::OutputFormat
GridOut::parse_output_format(const std::string &format_name)
{
  if (format_name == "none" || format_name == "false")
    return none;

  if (format_name == "dx")
    return dx;

  if (format_name == "ucd")
    return ucd;

  if (format_name == "gnuplot")
    return gnuplot;

  if (format_name == "eps")
    return eps;

  if (format_name == "xfig")
    return xfig;

  if (format_name == "msh")
    return msh;

  if (format_name == "svg")
    return svg;

  if (format_name == "mathgl")
    return mathgl;

  if (format_name == "vtk")
    return vtk;

  if (format_name == "vtu")
    return vtu;

  AssertThrow(false, ExcInvalidState());
  // return something weird
  return OutputFormat(-1);
}



std::string
GridOut::get_output_format_names()
{
  return "none|dx|gnuplot|eps|ucd|xfig|msh|svg|mathgl|vtk|vtu";
}


void
GridOut::declare_parameters(ParameterHandler &param)
{
  param.declare_entry("Format",
                      "none",
                      Patterns::Selection(get_output_format_names()));

  param.enter_subsection("DX");
  GridOutFlags::DX::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Msh");
  GridOutFlags::Msh::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Ucd");
  GridOutFlags::Ucd::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Gnuplot");
  GridOutFlags::Gnuplot::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Eps");
  GridOutFlags::EpsFlagsBase::declare_parameters(param);
  GridOutFlags::Eps<1>::declare_parameters(param);
  GridOutFlags::Eps<2>::declare_parameters(param);
  GridOutFlags::Eps<3>::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("XFig");
  GridOutFlags::XFig::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("MathGL");
  GridOutFlags::MathGL::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Vtk");
  GridOutFlags::Vtk::declare_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Vtu");
  GridOutFlags::Vtu::declare_parameters(param);
  param.leave_subsection();
}



void
GridOut::parse_parameters(ParameterHandler &param)
{
  default_format = parse_output_format(param.get("Format"));

  param.enter_subsection("DX");
  dx_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Msh");
  msh_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Ucd");
  ucd_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Gnuplot");
  gnuplot_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Eps");
  eps_flags_1.parse_parameters(param);
  eps_flags_2.parse_parameters(param);
  eps_flags_3.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("XFig");
  xfig_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("MathGL");
  mathgl_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Vtk");
  vtk_flags.parse_parameters(param);
  param.leave_subsection();

  param.enter_subsection("Vtu");
  vtu_flags.parse_parameters(param);
  param.leave_subsection();
}



std::size_t
GridOut::memory_consumption() const
{
  return (sizeof(dx_flags) + sizeof(msh_flags) + sizeof(ucd_flags) +
          sizeof(gnuplot_flags) + sizeof(eps_flags_1) + sizeof(eps_flags_2) +
          sizeof(eps_flags_3) + sizeof(xfig_flags) + sizeof(svg_flags) +
          sizeof(mathgl_flags) + sizeof(vtk_flags) + sizeof(vtu_flags));
}



template <>
void
GridOut::write_dx(const Triangulation<1> &, std::ostream &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}

template <>
void
GridOut::write_dx(const Triangulation<1, 2> &, std::ostream &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}

template <>
void
GridOut::write_dx(const Triangulation<1, 3> &, std::ostream &) const
{
  DEAL_II_NOT_IMPLEMENTED();
}



template <int dim, int spacedim>
void
GridOut::write_dx(const Triangulation<dim, spacedim> &tria,
                  std::ostream                       &out) const
{
  // TODO:[GK] allow for boundary faces only
  Assert(dx_flags.write_all_faces, ExcNotImplemented());
  AssertThrow(out.fail() == false, ExcIO());
  // Copied and adapted from write_ucd
  const std::vector<Point<spacedim>> &vertices    = tria.get_vertices();
  const std::vector<bool>            &vertex_used = tria.get_used_vertices();

  const unsigned int n_vertices = tria.n_used_vertices();

  // vertices are implicitly numbered from 0 to
  // n_vertices-1. we have to renumber the
  // vertices, because otherwise we would end
  // up with wrong results, if there are unused
  // vertices
  std::vector<unsigned int> renumber(vertices.size());
  // fill this vector with new vertex numbers
  // ranging from 0 to n_vertices-1
  unsigned int new_number = 0;
  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertex_used[i])
      renumber[i] = new_number++;
  Assert(new_number == n_vertices, ExcInternalError());

  // write the vertices
  out << "object \"vertices\" class array type float rank 1 shape " << dim
      << " items " << n_vertices << " data follows" << '\n';

  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertex_used[i])
      out << '\t' << vertices[i] << '\n';

  // write cells or faces
  const bool write_cells = dx_flags.write_cells;
  const bool write_faces = (dim > 1) ? dx_flags.write_faces : false;

  const unsigned int n_cells = tria.n_active_cells();
  const unsigned int n_faces =
    tria.n_active_cells() * GeometryInfo<dim>::faces_per_cell;

  const unsigned int n_vertices_per_cell = GeometryInfo<dim>::vertices_per_cell;
  const unsigned int n_vertices_per_face = GeometryInfo<dim>::vertices_per_face;

  if (write_cells)
    {
      out << "object \"cells\" class array type int rank 1 shape "
          << n_vertices_per_cell << " items " << n_cells << " data follows"
          << '\n';

      for (const auto &cell : tria.active_cell_iterators())
        {
          for (const unsigned int v : GeometryInfo<dim>::vertex_indices())
            out
              << '\t'
              << renumber[cell->vertex_index(GeometryInfo<dim>::dx_to_deal[v])];
          out << '\n';
        }
      out << "attribute \"element type\" string \"";
      if (dim == 1)
        out << "lines";
      if (dim == 2)
        out << "quads";
      if (dim == 3)
        out << "cubes";
      out << "\"" << '\n'
          << "attribute \"ref\" string \"positions\"" << '\n'
          << '\n';

      // Additional cell information

      out << "object \"material\" class array type int rank 0 items " << n_cells
          << " data follows" << '\n';
      for (const auto &cell : tria.active_cell_iterators())
        out << ' ' << cell->material_id();
      out << '\n' << "attribute \"dep\" string \"connections\"" << '\n' << '\n';

      out << "object \"level\" class array type int rank 0 items " << n_cells
          << " data follows" << '\n';
      for (const auto &cell : tria.active_cell_iterators())
        out << ' ' << cell->level();
      out << '\n' << "attribute \"dep\" string \"connections\"" << '\n' << '\n';

      if (dx_flags.write_measure)
        {
          out << "object \"measure\" class array type float rank 0 items "
              << n_cells << " data follows" << '\n';
          for (const auto &cell : tria.active_cell_iterators())
            out << '\t' << cell->measure();
          out << '\n'
              << "attribute \"dep\" string \"connections\"" << '\n'
              << '\n';
        }

      if (dx_flags.write_diameter)
        {
          out << "object \"diameter\" class array type float rank 0 items "
              << n_cells << " data follows" << '\n';
          for (const auto &cell : tria.active_cell_iterators())
            out << '\t' << cell->diameter();
          out << '\n'
              << "attribute \"dep\" string \"connections\"" << '\n'
              << '\n';
        }
    }

  if (write_faces)
    {
      out << "object \"faces\" class array type int rank 1 shape "
          << n_vertices_per_face << " items " << n_faces << " data follows"
          << '\n';

      for (const auto &cell : tria.active_cell_iterators())
        {
          for (const unsigned int f : cell->face_indices())
            {
              typename Triangulation<dim, spacedim>::face_iterator face =
                cell->face(f);

              for (unsigned int v = 0; v < GeometryInfo<dim>::vertices_per_face;
                   ++v)
                out << '\t'
                    << renumber[face->vertex_index(
                         GeometryInfo<dim - 1>::dx_to_deal[v])];
              out << '\n';
            }
        }
      out << "attribute \"element type\" string \"";
      if (dim == 2)
        out << "lines";
      if (dim == 3)
        out << "quads";
      out << "\"" << '\n'
          << "attribute \"ref\" string \"positions\"" << '\n'
          << '\n';


      // Additional face information

      out << "object \"boundary\" class array type int rank 0 items " << n_faces
          << " data follows" << '\n';
      for (const auto &cell : tria.active_cell_iterators())
        {
          // Little trick to get -1 for the interior
          for (const unsigned int f : GeometryInfo<dim>::face_indices())
            {
              out << ' '
                  << static_cast<std::make_signed_t<types::boundary_id>>(
                       cell->face(f)->boundary_id());
            }
          out << '\n';
        }
      out << "attribute \"dep\" string \"connections\"" << '\n' << '\n';

      if (dx_flags.write_measure)
        {
          out << "object \"face measure\" class array type float rank 0 items "
              << n_faces << " data follows" << '\n';
          for (const auto &cell : tria.active_cell_iterators())
            {
              for (const unsigned int f : GeometryInfo<dim>::face_indices())
                out << ' ' << cell->face(f)->measure();
              out << '\n';
            }
          out << "attribute \"dep\" string \"connections\"" << '\n' << '\n';
        }

      if (dx_flags.write_diameter)
        {
          out << "object \"face diameter\" class array type float rank 0 items "
              << n_faces << " data follows" << '\n';
          for (const auto &cell : tria.active_cell_iterators())
            {
              for (const unsigned int f : GeometryInfo<dim>::face_indices())
                out << ' ' << cell->face(f)->diameter();
              out << '\n';
            }
          out << "attribute \"dep\" string \"connections\"" << '\n' << '\n';
        }
    }


  // Write additional face information

  if (write_faces)
    {
    }
  else
    {
    }

  // The wrapper
  out << "object \"deal data\" class field" << '\n'
      << "component \"positions\" value \"vertices\"" << '\n'
      << "component \"connections\" value \"cells\"" << '\n';

  if (write_cells)
    {
      out << "object \"cell data\" class field" << '\n'
          << "component \"positions\" value \"vertices\"" << '\n'
          << "component \"connections\" value \"cells\"" << '\n';
      out << "component \"material\" value \"material\"" << '\n';
      out << "component \"level\" value \"level\"" << '\n';
      if (dx_flags.write_measure)
        out << "component \"measure\" value \"measure\"" << '\n';
      if (dx_flags.write_diameter)
        out << "component \"diameter\" value \"diameter\"" << '\n';
    }

  if (write_faces)
    {
      out << "object \"face data\" class field" << '\n'
          << "component \"positions\" value \"vertices\"" << '\n'
          << "component \"connections\" value \"faces\"" << '\n';
      out << "component \"boundary\" value \"boundary\"" << '\n';
      if (dx_flags.write_measure)
        out << "component \"measure\" value \"face measure\"" << '\n';
      if (dx_flags.write_diameter)
        out << "component \"diameter\" value \"face diameter\"" << '\n';
    }

  out << '\n' << "object \"grid data\" class group" << '\n';
  if (write_cells)
    out << "member \"cells\" value \"cell data\"" << '\n';
  if (write_faces)
    out << "member \"faces\" value \"face data\"" << '\n';
  out << "end" << '\n';

  // make sure everything now gets to
  // disk
  out.flush();

  AssertThrow(out.fail() == false, ExcIO());
}



template <int dim, int spacedim>
void
GridOut::write_msh(const Triangulation<dim, spacedim> &tria,
                   std::ostream                       &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // get the positions of the
  // vertices and whether they are
  // used.
  const std::vector<Point<spacedim>> &vertices    = tria.get_vertices();
  const std::vector<bool>            &vertex_used = tria.get_used_vertices();

  const unsigned int n_vertices = tria.n_used_vertices();

  // Write Header
  // The file format is:
  /*


  $NOD
  number-of-nodes
  node-number x-coord y-coord z-coord
  ...
  $ENDNOD
  $ELM
  number-of-elements
  elm-number elm-type reg-phys reg-elem number-of-nodes node-number-list
  ...
  $ENDELM
  */
  out << "$NOD" << '\n' << n_vertices << '\n';

  // actually write the vertices.
  // note that we shall number them
  // with first index 1 instead of 0
  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertex_used[i])
      {
        out << i + 1 // vertex index
            << "  " << vertices[i];
        for (unsigned int d = spacedim + 1; d <= 3; ++d)
          out << " 0"; // fill with zeroes
        out << '\n';
      }

  // Write cells preamble
  out << "$ENDNOD" << '\n'
      << "$ELM" << '\n'
      << tria.n_active_cells() +
           ((msh_flags.write_faces ? n_boundary_faces(tria) : 0) +
            (msh_flags.write_lines ? n_boundary_lines(tria) : 0))
      << '\n';

  static constexpr std::array<unsigned int, 8> local_vertex_numbering = {
    {0, 1, 5, 4, 2, 3, 7, 6}};

  // write cells. Enumerate cells
  // consecutively, starting with 1
  for (const auto &cell : tria.active_cell_iterators())
    {
      out << cell->active_cell_index() + 1 << ' '
          << cell->reference_cell().gmsh_element_type() << ' '
          << cell->material_id() << ' ' << cell->subdomain_id() << ' '
          << cell->n_vertices() << ' ';

      // Vertex numbering follows UCD conventions.

      for (const unsigned int vertex : cell->vertex_indices())
        {
          if (cell->reference_cell() == ReferenceCells::get_hypercube<dim>())
            out << cell->vertex_index(
                     dim == 3 ? local_vertex_numbering[vertex] :
                                GeometryInfo<dim>::ucd_to_deal[vertex]) +
                     1
                << ' ';
          else if (cell->reference_cell() == ReferenceCells::get_simplex<dim>())
            out << cell->vertex_index(vertex) + 1 << ' ';
          else
            DEAL_II_NOT_IMPLEMENTED();
        }
      out << '\n';
    }

  // write faces and lines with non-zero boundary indicator
  unsigned int next_element_index = tria.n_active_cells() + 1;
  if (msh_flags.write_faces)
    {
      next_element_index = write_msh_faces(tria, next_element_index, out);
    }
  if (msh_flags.write_lines)
    {
      next_element_index = write_msh_lines(tria, next_element_index, out);
    }

  out << "$ENDELM\n";

  // make sure everything now gets to
  // disk
  out.flush();

  AssertThrow(out.fail() == false, ExcIO());
}


template <int dim, int spacedim>
void
GridOut::write_ucd(const Triangulation<dim, spacedim> &tria,
                   std::ostream                       &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // get the positions of the
  // vertices and whether they are
  // used.
  const std::vector<Point<spacedim>> &vertices    = tria.get_vertices();
  const std::vector<bool>            &vertex_used = tria.get_used_vertices();

  const unsigned int n_vertices = tria.n_used_vertices();

  // write preamble
  if (ucd_flags.write_preamble)
    {
      // block this to have local
      // variables destroyed after
      // use
      std::time_t time1 = std::time(nullptr);
      std::tm    *time  = std::localtime(&time1);
      out
        << "# This file was generated by the deal.II library." << '\n'
        << "# Date =  " << time->tm_year + 1900 << "/" << time->tm_mon + 1
        << "/" << time->tm_mday << '\n'
        << "# Time =  " << time->tm_hour << ":" << std::setw(2) << time->tm_min
        << ":" << std::setw(2) << time->tm_sec << '\n'
        << "#" << '\n'
        << "# For a description of the UCD format see the AVS Developer's guide."
        << '\n'
        << "#" << '\n';
    }

  // start with ucd data
  out << n_vertices << ' '
      << tria.n_active_cells() +
           ((ucd_flags.write_faces ? n_boundary_faces(tria) : 0) +
            (ucd_flags.write_lines ? n_boundary_lines(tria) : 0))
      << " 0 0 0" // no data
      << '\n';

  // actually write the vertices.
  // note that we shall number them
  // with first index 1 instead of 0
  for (unsigned int i = 0; i < vertices.size(); ++i)
    if (vertex_used[i])
      {
        out << i + 1 // vertex index
            << "  " << vertices[i];
        for (unsigned int d = spacedim + 1; d <= 3; ++d)
          out << " 0"; // fill with zeroes
        out << '\n';
      }

  // write cells. Enumerate cells
  // consecutively, starting with 1
  for (const auto &cell : tria.active_cell_iterators())
    {
      out << cell->active_cell_index() + 1 << ' ' << cell->material_id() << ' ';
      switch (dim)
        {
          case 1:
            out << "line    ";
            break;
          case 2:
            out << "quad    ";
            break;
          case 3:
            out << "hex     ";
            break;
          default:
            DEAL_II_NOT_IMPLEMENTED();
        }

      // it follows a list of the
      // vertices of each cell. in 1d
      // this is simply a list of the
      // two vertices, in 2d its counter
      // clockwise, as usual in this
      // library. in 3d, the same applies
      // (special thanks to AVS for
      // numbering their vertices in a
      // way compatible to deal.II!)
      //
      // technical reference:
      // AVS Developer's Guide, Release 4,
      // May, 1992, p. E6
      //
      // note: vertex numbers are 1-base
      for (const unsigned int vertex : GeometryInfo<dim>::vertex_indices())
        out << cell->vertex_index(GeometryInfo<dim>::ucd_to_deal[vertex]) + 1
            << ' ';
      out << '\n';
    }

  // write faces and lines with non-zero boundary indicator
  unsigned int next_element_index = tria.n_active_cells() + 1;
  if (ucd_flags.write_faces)
    {
      next_element_index = write_ucd_faces(tria, next_element_index, out);
    }
  if (ucd_flags.write_lines)
    {
      next_element_index = write_ucd_lines(tria, next_element_index, out);
    }

  // make sure everything now gets to
  // disk
  out.flush();

  AssertThrow(out.fail() == false, ExcIO());
}



template <int dim, int spacedim>
void
GridOut::write_xfig(const Triangulation<dim, spacedim> &,
                    std::ostream &,
                    const Mapping<dim, spacedim> *) const
{
  DEAL_II_NOT_IMPLEMENTED();
}


// TODO:[GK] Obey parameters
template <>
void
GridOut::write_xfig(const Triangulation<2> &tria,
                    std::ostream           &out,
                    const Mapping<2> * /*mapping*/) const
{
  const int dim      = 2;
  const int spacedim = 2;

  const unsigned int nv = GeometryInfo<dim>::vertices_per_cell;

  // The following text was copied
  // from an existing XFig file.
  out << "#FIG 3.2\nLandscape\nCenter\nInches" << std::endl
      << "A4\n100.00\nSingle"
      << std::endl
      // Background is transparent
      << "-3" << std::endl
      << "# generated by deal.II GridOut class" << std::endl
      << "# reduce first number to scale up image" << std::endl
      << "1200 2" << std::endl;
  // Write custom palette
  // grey
  unsigned int colno = 32;
  out << "0 " << colno++ << " #ff0000" << std::endl;
  out << "0 " << colno++ << " #ff8000" << std::endl;
  out << "0 " << colno++ << " #ffd000" << std::endl;
  out << "0 " << colno++ << " #ffff00" << std::endl;
  out << "0 " << colno++ << " #c0ff00" << std::endl;
  out << "0 " << colno++ << " #80ff00" << std::endl;
  out << "0 " << colno++ << " #00f000" << std::endl;
  out << "0 " << colno++ << " #00f0c0" << std::endl;
  out << "0 " << colno++ << " #00f0ff" << std::endl;
  out << "0 " << colno++ << " #00c0ff" << std::endl;
  out << "0 " << colno++ << " #0080ff" << std::endl;
  out << "0 " << colno++ << " #0040ff" << std::endl;
  out << "0 " << colno++ << " #0000c0" << std::endl;
  out << "0 " << colno++ << " #5000ff" << std::endl;
  out << "0 " << colno++ << " #8000ff" << std::endl;
  out << "0 " << colno++ << " #b000ff" << std::endl;
  out << "0 " << colno++ << " #ff00ff" << std::endl;
  out << "0 " << colno++ << " #ff80ff" << std::endl;
  // grey
  for (unsigned int i = 0; i < 8; ++i)
    out << "0 " << colno++ << " #" << std::hex << 32 * i + 31 << 32 * i + 31
        << 32 * i + 31 << std::dec << std::endl;
  // green
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #00" << std::hex << 16 * i + 15 << std::dec
        << "00" << std::endl;
  // yellow
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #" << std::hex << 16 * i + 15 << 16 * i + 15
        << std::dec << "00" << std::endl;
  // red
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #" << std::hex << 16 * i + 15 << std::dec
        << "0000" << std::endl;
  // purple
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #" << std::hex << 16 * i + 15 << "00"
        << 16 * i + 15 << std::dec << std::endl;
  // blue
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #0000" << std::hex << 16 * i + 15 << std::dec
        << std::endl;
  // cyan
  for (unsigned int i = 1; i < 16; ++i)
    out << "0 " << colno++ << " #00" << std::hex << 16 * i + 15 << 16 * i + 15
        << std::dec << std::endl;

  // We write all cells and cells on
  // coarser levels are behind cells
  // on finer levels. Level 0
  // corresponds to a depth of 900,
  // each level subtracting 1
  for (const auto &cell : tria.cell_iterators())
    {
      // If depth is not encoded, write finest level only
      if (!xfig_flags.level_depth && !cell->is_active())
        continue;
      // Code for polygon
      out << "2 3  " << xfig_flags.line_style << ' '
          << xfig_flags.line_thickness
          // with black line
          << " 0 ";
      // Fill color
      switch (xfig_flags.color_by)
        {
            // TODO[GK]: Simplify after deprecation period is over
          case GridOutFlags::XFig::material_id:
            out << cell->material_id() + 32;
            break;
          case GridOutFlags::XFig::level_number:
            out << cell->level() + 8;
            break;
          case GridOutFlags::XFig::subdomain_id:
            out << cell->subdomain_id() + 32;
            break;
          case GridOutFlags::XFig::level_subdomain_id:
            out << cell->level_subdomain_id() + 32;
            break;
          default:
            DEAL_II_ASSERT_UNREACHABLE();
        }

      // Depth, unused, fill
      out << ' '
          << (xfig_flags.level_depth ? (900 - cell->level()) :
                                       (900 + cell->material_id()))
          << " 0 " << xfig_flags.fill_style
          << " 0.0 "
          // some style parameters
          << " 0 0 -1 0 0 "
          // number of points
          << nv + 1 << std::endl;

      // For each point, write scaled
      // and shifted coordinates
      // multiplied by 1200
      // (dots/inch)
      for (unsigned int k = 0; k <= nv; ++k)
        {
          const Point<dim> &p =
            cell->vertex(GeometryInfo<dim>::ucd_to_deal[k % nv]);
          for (unsigned int d = 0; d < static_cast<unsigned int>(dim); ++d)
            {
              int val = static_cast<int>(1200 * xfig_flags.scaling[d] *
                                         (p[d] - xfig_flags.offset[d]));
              out << '\t' << ((d == 0) ? val : -val);
            }
          out << std::endl;
        }
      // Now write boundary edges
      static const unsigned int face_reorder[4] = {2, 1, 3, 0};
      if (xfig_flags.draw_boundary)
        for (const unsigned int f : face_reorder)
          {
            Triangulation<dim, spacedim>::face_iterator face = cell->face(f);
            const types::boundary_id bi = face->boundary_id();
            if (bi != numbers::internal_face_boundary_id)
              {
                // Code for polyline
                out << "2 1 "
                    // with line style and thickness
                    << xfig_flags.boundary_style << ' '
                    << xfig_flags.boundary_thickness << ' ' << 1 + bi;
                // Fill color
                out << " -1 ";
                // Depth 100 less than cells
                out << (xfig_flags.level_depth ? (800 - cell->level()) :
                                                 800 + bi)
                    // unused, no fill
                    << " 0 -1 0.0 "
                    // some style parameters
                    << " 0 0 -1 0 0 "
                    // number of points
                    << GeometryInfo<dim>::vertices_per_face << std::endl;

                // For each point, write scaled
                // and shifted coordinates
                // multiplied by 1200
                // (dots/inch)

                for (unsigned int k = 0;
                     k < GeometryInfo<dim>::vertices_per_face;
                     ++k)
                  {
                    const Point<dim> &p = face->vertex(k % nv);
                    for (unsigned int d = 0; d < static_cast<unsigned int>(dim);
                         ++d)
                      {
                        int val =
                          static_cast<int>(1200 * xfig_flags.scaling[d] *
                                           (p[d] - xfig_flags.offset[d]));
                        out << '\t' << ((d == 0) ? val : -val);
                      }
                    out << std::endl;
                  }
              }
          }
    }

  // make sure everything now gets to
  // disk
  out.flush();

  AssertThrow(out.fail() == false, ExcIO());
}



#ifdef DEAL_II_GMSH_WITH_API
template <int dim, int spacedim>
void
GridOut::write_msh(const Triangulation<dim, spacedim> &tria,
                   const std::string                  &filename) const
{
  // mesh Type renumbering
  const std::array<int, 8> dealii_to_gmsh_type = {{15, 1, 2, 3, 4, 7, 6, 5}};

  // Vertex renumbering, by dealii type
  const std::array<std::vector<unsigned int>, 8> dealii_to_gmsh = {
    {{0},
     {{0, 1}},
     {{0, 1, 2}},
     {{0, 1, 3, 2}},
     {{0, 1, 2, 3}},
     {{0, 1, 3, 2, 4}},
     {{0, 1, 2, 3, 4, 5}},
     {{0, 1, 3, 2, 4, 5, 7, 6}}}};

  // Extract all vertices (nodes in gmsh terminology), and store their three
  // dimensional coordinates (regardless of dim).
  const auto              &vertices = tria.get_vertices();
  std::vector<double>      coords(3 * vertices.size());
  std::vector<std::size_t> nodes(vertices.size());

  // Each node has a strictly positive tag. We assign simply its index+1.
  std::size_t i = 0;
  for (const auto &p : vertices)
    {
      for (unsigned int d = 0; d < spacedim; ++d)
        coords[i * 3 + d] = p[d];
      nodes[i] = i + 1;
      ++i;
    }

  // Construct one entity tag per boundary and manifold id pair.
  // We need to be smart here, in order to save some disk space. All cells need
  // to be written, but only faces and lines that have non default boundary ids
  // and/or manifold ids. We collect them into pairs, and for each unique pair,
  // we create a gmsh entity where we store the elements. Pre-count all the
  // entities, and make sure we know which pair refers to what entity and
  // vice-versa.
  using IdPair = std::pair<types::material_id, types::manifold_id>;
  std::map<IdPair, int> id_pair_to_entity_tag;
  std::vector<IdPair>   all_pairs;
  {
    std::set<IdPair> set_of_pairs;
    for (const auto &cell : tria.active_cell_iterators())
      {
        set_of_pairs.insert({cell->material_id(), cell->manifold_id()});
        for (const auto &f : cell->face_iterators())
          if (f->manifold_id() != numbers::flat_manifold_id ||
              (f->boundary_id() != 0 &&
               f->boundary_id() != numbers::internal_face_boundary_id))
            set_of_pairs.insert({f->boundary_id(), f->manifold_id()});
        if (dim > 2)
          for (const auto l : cell->line_indices())
            {
              const auto &f = cell->line(l);
              if (f->manifold_id() != numbers::flat_manifold_id ||
                  (f->boundary_id() != 0 &&
                   f->boundary_id() != numbers::internal_face_boundary_id))
                set_of_pairs.insert({f->boundary_id(), f->manifold_id()});
            }
      }
    all_pairs = {set_of_pairs.begin(), set_of_pairs.end()};

    int entity = 1;
    for (const auto &p : set_of_pairs)
      id_pair_to_entity_tag[p] = entity++;
  }

  const auto n_entity_tags = id_pair_to_entity_tag.size();

  // All elements in the mesh, by entity tag, and by dealii type.
  std::vector<std::vector<std::vector<std::size_t>>> element_ids(
    n_entity_tags, std::vector<std::vector<std::size_t>>(8));
  std::vector<std::vector<std::vector<std::size_t>>> element_nodes(
    n_entity_tags, std::vector<std::vector<std::size_t>>(8));

  // One element id counter for all dimensions.
  std::size_t element_id = 1;

  const auto add_element = [&](const auto &element, const int &entity_tag) {
    const auto type = element->reference_cell();

    Assert(entity_tag > 0, ExcInternalError());
    // Add all vertex ids. Make sure we renumber to gmsh, and we add 1 to the
    // global index.
    for (const auto v : element->vertex_indices())
      element_nodes[entity_tag - 1][type].emplace_back(
        element->vertex_index(dealii_to_gmsh[type][v]) + 1);

    // Save the element id.
    element_ids[entity_tag - 1][type].emplace_back(element_id);
    ++element_id;
  };

  // Will create a separate gmsh entity, only  if it's a cell, or if the
  // boundary and/or the manifold ids are not the default ones.
  // In the meanwhile, also store each pair of dimension and entity tag that was
  // requested.
  std::set<std::pair<int, int>> dim_entity_tag;

  auto maybe_add_element =
    [&](const auto               &element,
        const types::boundary_id &boundary_or_material_id) {
      const auto struct_dim  = element->structure_dimension;
      const auto manifold_id = element->manifold_id();

      // Exclude default boundary/manifold id or invalid/flag
      const bool non_default_boundary_or_material_id =
        (boundary_or_material_id != 0 &&
         boundary_or_material_id != numbers::internal_face_boundary_id);
      const bool non_default_manifold =
        manifold_id != numbers::flat_manifold_id;
      if (struct_dim == dim || non_default_boundary_or_material_id ||
          non_default_manifold)
        {
          const auto entity_tag =
            id_pair_to_entity_tag[{boundary_or_material_id, manifold_id}];
          add_element(element, entity_tag);
          dim_entity_tag.insert({struct_dim, entity_tag});
        }
    };

  // Loop recursively over all cells, faces, and possibly lines.
  for (const auto &cell : tria.active_cell_iterators())
    {
      maybe_add_element(cell, cell->material_id());
      for (const auto &face : cell->face_iterators())
        maybe_add_element(face, face->boundary_id());
      if (dim > 2)
        for (const auto l : cell->line_indices())
          maybe_add_element(cell->line(l), cell->line(l)->boundary_id());
    }

  // Now that we collected everything, plug them into gmsh
  gmsh::initialize();
  gmsh::option::setNumber("General.Verbosity", 0);
  gmsh::model::add("Grid generated in deal.II");
  for (const auto &p : dim_entity_tag)
    {
      gmsh::model::addDiscreteEntity(p.first, p.second);
      gmsh::model::mesh::addNodes(p.first, p.second, nodes, coords);
    }

  for (unsigned int entity_tag = 0; entity_tag < n_entity_tags; ++entity_tag)
    for (unsigned int t = 1; t < 8; ++t)
      {
        const auto all_element_ids   = element_ids[entity_tag][t];
        const auto all_element_nodes = element_nodes[entity_tag][t];
        const auto gmsh_t            = dealii_to_gmsh_type[t];
        if (all_element_ids.size() > 0)
          gmsh::model::mesh::addElementsByType(entity_tag + 1,
                                               gmsh_t,
                                               all_element_ids,
                                               all_element_nodes);
      }

  // Now for each individual pair of dim and entry, add a physical group, if
  // necessary
  for (const auto &it : dim_entity_tag)
    {
      const auto &d           = it.first;
      const auto &entity_tag  = it.second;
      const auto &boundary_id = all_pairs[entity_tag - 1].first;
      const auto &manifold_id = all_pairs[entity_tag - 1].second;

      std::string physical_name;
      if (d == dim && boundary_id != 0)
        physical_name += "MaterialID:" + Utilities::int_to_string(
                                           static_cast<int>(boundary_id));
      else if (d < dim && boundary_id != 0)
        physical_name +=
          "BoundaryID:" +
          (boundary_id == numbers::internal_face_boundary_id ?
             "-1" :
             Utilities::int_to_string(static_cast<int>(boundary_id)));

      std::string sep = physical_name != "" ? ", " : "";
      if (manifold_id != numbers::flat_manifold_id)
        physical_name +=
          sep + "ManifoldID:" +
          Utilities::int_to_string(static_cast<int>(manifold_id));
      const auto physical_tag =
        gmsh::model::addPhysicalGroup(d, {entity_tag}, -1);
      if (physical_name != "")
        gmsh::model::setPhysicalName(d, physical_tag, physical_name);
    }


  gmsh::write(filename);
  gmsh::clear();
  gmsh::finalize();
}
#endif



namespace
{
  /**
   * This function projects a three-dimensional point (Point<3> point) onto a
   * two-dimensional image plane, specified by the position of the camera
   * viewing system (Point<3> camera_position), camera direction (Point<3>
   * camera_position), camera horizontal (Point<3> camera_horizontal,
   * necessary for the correct alignment of the later images), and the focus
   * of the camera (float camera_focus).
   *
   * For SVG output of grids.
   */
  Point<2>
  svg_project_point(const Point<3>     &point,
                    const Point<3>     &camera_position,
                    const Tensor<1, 3> &camera_direction,
                    const Tensor<1, 3> &camera_horizontal,
                    const float         camera_focus)
  {
    const Tensor<1, 3> camera_vertical =
      cross_product_3d(camera_horizontal, camera_direction);

    const float phi =
      camera_focus / ((point - camera_position) * camera_direction);

    const Point<3> projection =
      camera_position + phi * (point - camera_position);

    return {(projection - camera_position - camera_focus * camera_direction) *
              camera_horizontal,
            (projection - camera_position - camera_focus * camera_direction) *
              camera_vertical};
  }
} // namespace



template <int dim, int spacedim>
void
GridOut::write_svg(const Triangulation<dim, spacedim> &,
                   std::ostream & /*out*/) const
{
  Assert(false,
         ExcMessage("Mesh output in SVG format is not implemented for anything "
                    "other than two-dimensional meshes in two-dimensional "
                    "space. That's because three-dimensional meshes are best "
                    "viewed in programs that allow changing the viewpoint, "
                    "but SVG format does not allow this: It is an inherently "
                    "2d format, and for three-dimensional meshes would "
                    "require choosing one, fixed viewpoint."
                    "\n\n"
                    "You probably want to output your mesh in a format such "
                    "as VTK, VTU, or gnuplot."));
}


void
GridOut::write_svg(const Triangulation<2, 2> &tria, std::ostream &out) const
{
  unsigned int n = 0;

  unsigned int min_level, max_level;

  // Svg files require an underlying drawing grid. The size of this
  // grid is provided in the parameters height and width. Each of them
  // may be zero, such that it is computed from the other. Obviously,
  // both of them zero does not produce reasonable output.
  unsigned int height = svg_flags.height;
  unsigned int width  = svg_flags.width;
  Assert(height != 0 || width != 0,
         ExcMessage("You have to set at least one of width and height"));

  unsigned int margin_in_percent = 0;
  if (svg_flags.margin || svg_flags.background == GridOutFlags::Svg::dealii)
    margin_in_percent = 8;

  // initial font size for cell labels
  unsigned int cell_label_font_size;

  // get date and time
  // time_t time_stamp;
  // tm *now;
  // time_stamp = time(0);
  // now = localtime(&time_stamp);

  float camera_focus;

  Point<3> point;
  Point<2> projection_decomposition;

  float x_max_perspective, x_min_perspective;
  float y_max_perspective, y_min_perspective;

  float x_dimension_perspective, y_dimension_perspective;


  // auxiliary variables for the bounding box and the range of cell levels
  double x_min = tria.begin()->vertex(0)[0];
  double x_max = x_min;
  double y_min = tria.begin()->vertex(0)[1];
  double y_max = y_min;

  double x_dimension, y_dimension;

  min_level = max_level = tria.begin()->level();

  // auxiliary set for the materials being used
  std::set<unsigned int> materials;

  // auxiliary set for the levels being used
  std::set<unsigned int> levels;

  // auxiliary set for the subdomains being used
  std::set<unsigned int> subdomains;

  // auxiliary set for the level subdomains being used
  std::set<int> level_subdomains;

  // We use an active cell iterator to determine the
  // bounding box of the given triangulation and check
  // the cells for material id, level number, subdomain id
  // (, and level subdomain id).
  for (const auto &cell : tria.cell_iterators())
    {
      for (unsigned int vertex_index = 0; vertex_index < cell->n_vertices();
           ++vertex_index)
        {
          if (cell->vertex(vertex_index)[0] < x_min)
            x_min = cell->vertex(vertex_index)[0];
          if (cell->vertex(vertex_index)[0] > x_max)
            x_max = cell->vertex(vertex_index)[0];

          if (cell->vertex(vertex_index)[1] < y_min)
            y_min = cell->vertex(vertex_index)[1];
          if (cell->vertex(vertex_index)[1] > y_max)
            y_max = cell->vertex(vertex_index)[1];
        }

      if (static_cast<unsigned int>(cell->level()) < min_level)
        min_level = cell->level();
      if (static_cast<unsigned int>(cell->level()) > max_level)
        max_level = cell->level();

      materials.insert(cell->material_id());
      levels.insert(cell->level());
      if (cell->is_active())
        subdomains.insert(cell->subdomain_id() + 2);
      level_subdomains.insert(cell->level_subdomain_id() + 2);
    }

  x_dimension = x_max - x_min;
  y_dimension = y_max - y_min;

  // count the materials being used
  const unsigned int n_materials = materials.size();

  // count the levels being used
  const unsigned int n_levels = levels.size();

  // count the subdomains being used
  const unsigned int n_subdomains = subdomains.size();

  // count the level subdomains being used
  const unsigned int n_level_subdomains = level_subdomains.size();

  switch (svg_flags.coloring)
    {
      case GridOutFlags::Svg::material_id:
        n = n_materials;
        break;
      case GridOutFlags::Svg::level_number:
        n = n_levels;
        break;
      case GridOutFlags::Svg::subdomain_id:
        n = n_subdomains;
        break;
      case GridOutFlags::Svg::level_subdomain_id:
        n = n_level_subdomains;
        break;
      default:
        break;
    }

  // set the camera position to top view, targeting at the origin
  // vectors and variables for the perspective view
  Point<3> camera_position;
  camera_position[0] = 0;
  camera_position[1] = 0;
  camera_position[2] = 2. * std::max(x_dimension, y_dimension);

  Tensor<1, 3> camera_direction;
  camera_direction[0] = 0;
  camera_direction[1] = 0;
  camera_direction[2] = -1;

  Tensor<1, 3> camera_horizontal;
  camera_horizontal[0] = 1;
  camera_horizontal[1] = 0;
  camera_horizontal[2] = 0;

  camera_focus = .5 * std::max(x_dimension, y_dimension);

  Point<3> camera_position_temp;
  Point<3> camera_direction_temp;
  Point<3> camera_horizontal_temp;

  const double angle_factor = 3.14159265 / 180.;

  // (I) rotate the camera to the chosen polar angle
  camera_position_temp[1] =
    std::cos(angle_factor * svg_flags.polar_angle) * camera_position[1] -
    std::sin(angle_factor * svg_flags.polar_angle) * camera_position[2];
  camera_position_temp[2] =
    std::sin(angle_factor * svg_flags.polar_angle) * camera_position[1] +
    std::cos(angle_factor * svg_flags.polar_angle) * camera_position[2];

  camera_direction_temp[1] =
    std::cos(angle_factor * svg_flags.polar_angle) * camera_direction[1] -
    std::sin(angle_factor * svg_flags.polar_angle) * camera_direction[2];
  camera_direction_temp[2] =
    std::sin(angle_factor * svg_flags.polar_angle) * camera_direction[1] +
    std::cos(angle_factor * svg_flags.polar_angle) * camera_direction[2];

  camera_horizontal_temp[1] =
    std::cos(angle_factor * svg_flags.polar_angle) * camera_horizontal[1] -
    std::sin(angle_factor * svg_flags.polar_angle) * camera_horizontal[2];
  camera_horizontal_temp[2] =
    std::sin(angle_factor * svg_flags.polar_angle) * camera_horizontal[1] +
    std::cos(angle_factor * svg_flags.polar_angle) * camera_horizontal[2];

  camera_position[1] = camera_position_temp[1];
  camera_position[2] = camera_position_temp[2];

  camera_direction[1] = camera_direction_temp[1];
  camera_direction[2] = camera_direction_temp[2];

  camera_horizontal[1] = camera_horizontal_temp[1];
  camera_horizontal[2] = camera_horizontal_temp[2];

  // (II) rotate the camera to the chosen azimuth angle
  camera_position_temp[0] =
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_position[0] -
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_position[1];
  camera_position_temp[1] =
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_position[0] +
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_position[1];

  camera_direction_temp[0] =
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_direction[0] -
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_direction[1];
  camera_direction_temp[1] =
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_direction[0] +
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_direction[1];

  camera_horizontal_temp[0] =
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_horizontal[0] -
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_horizontal[1];
  camera_horizontal_temp[1] =
    std::sin(angle_factor * svg_flags.azimuth_angle) * camera_horizontal[0] +
    std::cos(angle_factor * svg_flags.azimuth_angle) * camera_horizontal[1];

  camera_position[0] = camera_position_temp[0];
  camera_position[1] = camera_position_temp[1];

  camera_direction[0] = camera_direction_temp[0];
  camera_direction[1] = camera_direction_temp[1];

  camera_horizontal[0] = camera_horizontal_temp[0];
  camera_horizontal[1] = camera_horizontal_temp[1];

  // translate the camera to the given triangulation
  camera_position[0] = x_min + .5 * x_dimension;
  camera_position[1] = y_min + .5 * y_dimension;

  camera_position[0] += 2. * std::max(x_dimension, y_dimension) *
                        std::sin(angle_factor * svg_flags.polar_angle) *
                        std::sin(angle_factor * svg_flags.azimuth_angle);
  camera_position[1] -= 2. * std::max(x_dimension, y_dimension) *
                        std::sin(angle_factor * svg_flags.polar_angle) *
                        std::cos(angle_factor * svg_flags.azimuth_angle);


  // determine the bounding box of the given triangulation on the projection
  // plane of the camera viewing system
  point[0] = tria.begin()->vertex(0)[0];
  point[1] = tria.begin()->vertex(0)[1];
  point[2] = 0;

  float min_level_min_vertex_distance = 0;

  if (svg_flags.convert_level_number_to_height)
    {
      point[2] = svg_flags.level_height_factor *
                 (static_cast<float>(tria.begin()->level()) /
                  static_cast<float>(n_levels)) *
                 std::max(x_dimension, y_dimension);
    }

  projection_decomposition = svg_project_point(
    point, camera_position, camera_direction, camera_horizontal, camera_focus);

  x_max_perspective = projection_decomposition[0];
  x_min_perspective = projection_decomposition[0];

  y_max_perspective = projection_decomposition[1];
  y_min_perspective = projection_decomposition[1];

  for (const auto &cell : tria.cell_iterators())
    {
      point[0] = cell->vertex(0)[0];
      point[1] = cell->vertex(0)[1];
      point[2] = 0;

      if (svg_flags.convert_level_number_to_height)
        {
          point[2] =
            svg_flags.level_height_factor *
            (static_cast<float>(cell->level()) / static_cast<float>(n_levels)) *
            std::max(x_dimension, y_dimension);
        }

      projection_decomposition = svg_project_point(point,
                                                   camera_position,
                                                   camera_direction,
                                                   camera_horizontal,
                                                   camera_focus);

      if (x_max_perspective < projection_decomposition[0])
        x_max_perspective = projection_decomposition[0];
      if (x_min_perspective > projection_decomposition[0])
        x_min_perspective = projection_decomposition[0];

      if (y_max_perspective < projection_decomposition[1])
        y_max_perspective = projection_decomposition[1];
      if (y_min_perspective > projection_decomposition[1])
        y_min_perspective = projection_decomposition[1];

      point[0] = cell->vertex(1)[0];
      point[1] = cell->vertex(1)[1];

      projection_decomposition = svg_project_point(point,
                                                   camera_position,
                                                   camera_direction,
                                                   camera_horizontal,
                                                   camera_focus);

      if (x_max_perspective < projection_decomposition[0])
        x_max_perspective = projection_decomposition[0];
      if (x_min_perspective > projection_decomposition[0])
        x_min_perspective = projection_decomposition[0];

      if (y_max_perspective < projection_decomposition[1])
        y_max_perspective = projection_decomposition[1];
      if (y_min_perspective > projection_decomposition[1])
        y_min_perspective = projection_decomposition[1];

      point[0] = cell->vertex(2)[0];
      point[1] = cell->vertex(2)[1];

      projection_decomposition = svg_project_point(point,
                                                   camera_position,
                                                   camera_direction,
                                                   camera_horizontal,
                                                   camera_focus);

      if (x_max_perspective < projection_decomposition[0])
        x_max_perspective = projection_decomposition[0];
      if (x_min_perspective > projection_decomposition[0])
        x_min_perspective = projection_decomposition[0];

      if (y_max_perspective < projection_decomposition[1])
        y_max_perspective = projection_decomposition[1];
      if (y_min_perspective > projection_decomposition[1])
        y_min_perspective = projection_decomposition[1];

      if (cell->n_vertices() == 4) // in case of quadrilateral
        {
          point[0] = cell->vertex(3)[0];
          point[1] = cell->vertex(3)[1];

          projection_decomposition = svg_project_point(point,
                                                       camera_position,
                                                       camera_direction,
                                                       camera_horizontal,
                                                       camera_focus);

          if (x_max_perspective < projection_decomposition[0])
            x_max_perspective = projection_decomposition[0];
          if (x_min_perspective > projection_decomposition[0])
            x_min_perspective = projection_decomposition[0];

          if (y_max_perspective < projection_decomposition[1])
            y_max_perspective = projection_decomposition[1];
          if (y_min_perspective > projection_decomposition[1])
            y_min_perspective = projection_decomposition[1];
        }

      if (static_cast<unsigned int>(cell->level()) == min_level)
        min_level_min_vertex_distance = cell->minimum_vertex_distance();
    }

  x_dimension_perspective = x_max_perspective - x_min_perspective;
  y_dimension_perspective = y_max_perspective - y_min_perspective;

  // create the svg file with an internal style sheet
  if (width == 0)
    width = static_cast<unsigned int>(
      .5 + height * (x_dimension_perspective / y_dimension_perspective));
  else if (height == 0)
    height = static_cast<unsigned int>(
      .5 + width * (y_dimension_perspective / x_dimension_perspective));
  unsigned int additional_width = 0;
  // font size for date, time, legend, and colorbar
  unsigned int font_size =
    static_cast<unsigned int>(.5 + (height / 100.) * 1.75);
  cell_label_font_size = static_cast<unsigned int>(
    .5 + (height * .15 * svg_flags.cell_font_scaling *
          min_level_min_vertex_distance / std::min(x_dimension, y_dimension)));

  if (svg_flags.draw_legend &&
      (svg_flags.label_level_number || svg_flags.label_cell_index ||
       svg_flags.label_material_id || svg_flags.label_subdomain_id ||
       svg_flags.label_level_subdomain_id))
    {
      additional_width = static_cast<unsigned int>(
        .5 + height * .4); // additional width for legend
    }
  else if (svg_flags.draw_colorbar && (svg_flags.coloring != 0u))
    {
      additional_width = static_cast<unsigned int>(
        .5 + height * .175); // additional width for colorbar
    }

  // out << "<!-- deal.ii GridOut " << now->tm_mday << '/' << now->tm_mon + 1 <<
  // '/' << now->tm_year + 1900
  //    << ' ' << now->tm_hour << ':';
  //
  // if (now->tm_min < 10) out << '0';
  //
  // out << now->tm_min << " -->" << '\n';

  // basic svg header
  out << "<svg width=\"" << width + additional_width << "\" height=\"" << height
      << "\" xmlns=\"http://www.w3.org/2000/svg\" version=\"1.1\">" << '\n'
      << '\n';


  if (svg_flags.background == GridOutFlags::Svg::dealii)
    {
      out
        << " <linearGradient id=\"background_gradient\" gradientUnits=\"userSpaceOnUse\" x1=\"0\" y1=\"0\" x2=\"0\" y2=\""
        << height << "\">" << '\n'
        << "  <stop offset=\"0\" style=\"stop-color:white\"/>" << '\n'
        << "  <stop offset=\"1\" style=\"stop-color:lightsteelblue\"/>" << '\n'
        << " </linearGradient>" << '\n';
    }

  out << '\n';

  // header for the internal style sheet
  out << "<!-- internal style sheet -->" << '\n'
      << "<style type=\"text/css\"><![CDATA[" << '\n';

  // set the background of the output graphic
  if (svg_flags.background == GridOutFlags::Svg::dealii)
    out << " rect.background{fill:url(#background_gradient)}" << '\n';
  else if (svg_flags.background == GridOutFlags::Svg::white)
    out << " rect.background{fill:white}" << '\n';
  else
    out << " rect.background{fill:none}" << '\n';

  // basic svg graphic element styles
  out << " rect{fill:none; stroke:rgb(25,25,25); stroke-width:"
      << svg_flags.line_thickness << '}' << '\n'
      << " text{font-family:Helvetica; text-anchor:middle; fill:rgb(25,25,25)}"
      << '\n'
      << " line{stroke:rgb(25,25,25); stroke-width:"
      << svg_flags.boundary_line_thickness << '}' << '\n'
      << " path{fill:none; stroke:rgb(25,25,25); stroke-width:"
      << svg_flags.line_thickness << '}' << '\n'
      << " circle{fill:white; stroke:black; stroke-width:2}" << '\n'
      << '\n';

  // polygon styles with respect to the chosen cell coloring
  if (svg_flags.coloring != 0u)
    {
      unsigned int labeling_index      = 0;
      auto         materials_it        = materials.begin();
      auto         levels_it           = levels.begin();
      auto         subdomains_it       = subdomains.begin();
      auto         level_subdomains_it = level_subdomains.begin();

      for (unsigned int index = 0; index < n; ++index)
        {
          double h;

          if (n != 1)
            {
              // The assert is a workaround for a compiler bug in ROCm 5.7 which
              // evaluated index/(n-1) when n == 1 in debug mode. When adding
              // the assert the ratio is not evaluated.
              Assert((n - 1.) != 0., ExcInvalidState());
              h = .6 - (index / (n - 1.)) * .6;
            }
          else
            h = .6;

          unsigned int r = 0;
          unsigned int g = 0;
          unsigned int b = 0;

          unsigned int i = static_cast<unsigned int>(h * 6);

          double f = h * 6 - i;
          double q = 1 - f;
          double t = f;

          switch (i % 6)
            {
              case 0:
                r = 255, g = static_cast<unsigned int>(.5 + 255 * t);
                break;
              case 1:
                r = static_cast<unsigned int>(.5 + 255 * q), g = 255;
                break;
              case 2:
                g = 255, b = static_cast<unsigned int>(.5 + 255 * t);
                break;
              case 3:
                g = static_cast<unsigned int>(.5 + 255 * q), b = 255;
                break;
              case 4:
                r = static_cast<unsigned int>(.5 + 255 * t), b = 255;
                break;
              case 5:
                r = 255, b = static_cast<unsigned int>(.5 + 255 * q);
                break;
              default:
                break;
            }

          switch (svg_flags.coloring)
            {
              case GridOutFlags::Svg::material_id:
                labeling_index = *materials_it++;
                break;
              case GridOutFlags::Svg::level_number:
                labeling_index = *levels_it++;
                break;
              case GridOutFlags::Svg::subdomain_id:
                labeling_index = *subdomains_it++;
                break;
              case GridOutFlags::Svg::level_subdomain_id:
                labeling_index = *level_subdomains_it++;
                break;
              default:
                break;
            }

          out << " path.p" << labeling_index << "{fill:rgb(" << r << ',' << g
              << ',' << b << "); "
              << "stroke:rgb(25,25,25); stroke-width:"
              << svg_flags.line_thickness << '}' << '\n';

          out << " path.ps" << labeling_index << "{fill:rgb("
              << static_cast<unsigned int>(.5 + .75 * r) << ','
              << static_cast<unsigned int>(.5 + .75 * g) << ','
              << static_cast<unsigned int>(.5 + .75 * b) << "); "
              << "stroke:rgb(20,20,20); stroke-width:"
              << svg_flags.line_thickness << '}' << '\n';

          out << " rect.r" << labeling_index << "{fill:rgb(" << r << ',' << g
              << ',' << b << "); "
              << "stroke:rgb(25,25,25); stroke-width:"
              << svg_flags.line_thickness << '}' << '\n';

          ++labeling_index;
        }
    }

  out << "]]></style>" << '\n' << '\n';

  // background rectangle
  out << " <rect class=\"background\" width=\"" << width << "\" height=\""
      << height << "\"/>" << '\n';

  if (svg_flags.background == GridOutFlags::Svg::dealii)
    {
      unsigned int x_offset = 0;

      if (svg_flags.margin)
        x_offset = static_cast<unsigned int>(.5 + (height / 100.) *
                                                    (margin_in_percent / 2.));
      else
        x_offset = static_cast<unsigned int>(.5 + height * .025);

      out
        << " <text x=\"" << x_offset << "\" y=\""
        << static_cast<unsigned int>(.5 + height * .0525) << '\"'
        << " style=\"font-weight:100; fill:lightsteelblue; text-anchor:start; font-family:Courier; font-size:"
        << static_cast<unsigned int>(.5 + height * .045) << "px\">"
        << "deal.II"
        << "</text>" << '\n';

      // out << " <text x=\"" << x_offset + static_cast<unsigned int>(.5 +
      // height * .045 * 4.75) << "\" y=\"" << static_cast<unsigned int>(.5 +
      // height * .0525) << '\"'
      //     << " style=\"fill:lightsteelblue; text-anchor:start; font-size:" <<
      //     font_size << "\">"
      //     << now->tm_mday << '/' << now->tm_mon + 1 << '/' << now->tm_year +
      //     1900
      //     << " - " << now->tm_hour << ':';
      //
      // if(now->tm_min < 10) out << '0';
      //
      // out << now->tm_min
      //     << "</text>"<< '\n' << '\n';
    }

  // draw the cells, starting out from the minimal level (in order to guaranty a
  // correct perspective view)
  out << "  <!-- cells -->" << '\n';

  for (unsigned int level_index = min_level; level_index <= max_level;
       level_index++)
    {
      for (const auto &cell : tria.cell_iterators_on_level(level_index))
        {
          if (!svg_flags.convert_level_number_to_height && !cell->is_active())
            continue;

          // draw the current cell
          out << "  <path";

          if (svg_flags.coloring != 0u)
            {
              out << " class=\"p";

              if (!cell->is_active() &&
                  svg_flags.convert_level_number_to_height)
                out << 's';

              switch (svg_flags.coloring)
                {
                  case GridOutFlags::Svg::material_id:
                    out << cell->material_id();
                    break;
                  case GridOutFlags::Svg::level_number:
                    out << static_cast<unsigned int>(cell->level());
                    break;
                  case GridOutFlags::Svg::subdomain_id:
                    if (cell->is_active())
                      out << cell->subdomain_id() + 2;
                    else
                      out << 'X';
                    break;
                  case GridOutFlags::Svg::level_subdomain_id:
                    out << cell->level_subdomain_id() + 2;
                    break;
                  default:
                    break;
                }

              out << '\"';
            }

          out << " d=\"M ";

          point[0] = cell->vertex(0)[0];
          point[1] = cell->vertex(0)[1];
          point[2] = 0;

          if (svg_flags.convert_level_number_to_height)
            {
              point[2] = svg_flags.level_height_factor *
                         (static_cast<float>(cell->level()) /
                          static_cast<float>(n_levels)) *
                         std::max(x_dimension, y_dimension);
            }

          projection_decomposition = svg_project_point(point,
                                                       camera_position,
                                                       camera_direction,
                                                       camera_horizontal,
                                                       camera_focus);

          out << static_cast<unsigned int>(
                   .5 +
                   ((projection_decomposition[0] - x_min_perspective) /
                    x_dimension_perspective) *
                     (width - (width / 100.) * 2. * margin_in_percent) +
                   ((width / 100.) * margin_in_percent))
              << ' '
              << static_cast<unsigned int>(
                   .5 + height - (height / 100.) * margin_in_percent -
                   ((projection_decomposition[1] - y_min_perspective) /
                    y_dimension_perspective) *
                     (height - (height / 100.) * 2. * margin_in_percent));

          out << " L ";

          point[0] = cell->vertex(1)[0];
          point[1] = cell->vertex(1)[1];

          projection_decomposition = svg_project_point(point,
                                                       camera_position,
                                                       camera_direction,
                                                       camera_horizontal,
                                                       camera_focus);

          out << static_cast<unsigned int>(
                   .5 +
                   ((projection_decomposition[0] - x_min_perspective) /
                    x_dimension_perspective) *
                     (width - (width / 100.) * 2. * margin_in_percent) +
                   ((width / 100.) * margin_in_percent))
              << ' '
              << static_cast<unsigned int>(
                   .5 + height - (height / 100.) * margin_in_percent -
                   ((projection_decomposition[1] - y_min_perspective) /
                    y_dimension_perspective) *
                     (height - (height / 100.) * 2. * margin_in_percent));

          out << " L ";

          if (cell->n_vertices() == 4) // in case of quadrilateral
            {
              point[0] = cell->vertex(3)[0];
              point[1] = cell->vertex(3)[1];

              projection_decomposition = svg_project_point(point,
                                                           camera_position,
                                                           camera_direction,
                                                           camera_horizontal,
                                                           camera_focus);

              out << static_cast<unsigned int>(
                       .5 +
                       ((projection_decomposition[0] - x_min_perspective) /
                        x_dimension_perspective) *
                         (width - (width / 100.) * 2. * margin_in_percent) +
                       ((width / 100.) * margin_in_percent))
                  << ' '
                  << static_cast<unsigned int>(
                       .5 + height - (height / 100.) * margin_in_percent -
                       ((projection_decomposition[1] - y_min_perspective) /
                        y_dimension_perspective) *
                         (height - (height / 100.) * 2. * margin_in_percent));

              out << " L ";
            }

          point[0] = cell->vertex(2)[0];
          point[1] = cell->vertex(2)[1];

          projection_decomposition = svg_project_point(point,
                                                       camera_position,
                                                       camera_direction,
                                                       camera_horizontal,
                                                       camera_focus);

          out << static_cast<unsigned int>(
                   .5 +
                   ((projection_decomposition[0] - x_min_perspective) /
                    x_dimension_perspective) *
                     (width - (width / 100.) * 2. * margin_in_percent) +
                   ((width / 100.) * margin_in_percent))
              << ' '
              << static_cast<unsigned int>(
                   .5 + height - (height / 100.) * margin_in_percent -
                   ((projection_decomposition[1] - y_min_perspective) /
                    y_dimension_perspective) *
                     (height - (height / 100.) * 2. * margin_in_percent));

          out << " L ";

          point[0] = cell->vertex(0)[0];
          point[1] = cell->vertex(0)[1];

          projection_decomposition = svg_project_point(point,
                                                       camera_position,
                                                       camera_direction,
                                                       camera_horizontal,
                                                       camera_focus);

          out << static_cast<unsigned int>(
                   .5 +
                   ((projection_decomposition[0] - x_min_perspective) /
                    x_dimension_perspective) *
                     (width - (width / 100.) * 2. * margin_in_percent) +
                   ((width / 100.) * margin_in_percent))
              << ' '
              << static_cast<unsigned int>(
                   .5 + height - (height / 100.) * margin_in_percent -
                   ((projection_decomposition[1] - y_min_perspective) /
                    y_dimension_perspective) *
                     (height - (height / 100.) * 2. * margin_in_percent));

          out << "\"/>" << '\n';

          // label the current cell
          if (svg_flags.label_level_number || svg_flags.label_cell_index ||
              svg_flags.label_material_id || svg_flags.label_subdomain_id ||
              svg_flags.label_level_subdomain_id)
            {
              point[0] = cell->center()[0];
              point[1] = cell->center()[1];
              point[2] = 0;

              if (svg_flags.convert_level_number_to_height)
                {
                  point[2] = svg_flags.level_height_factor *
                             (static_cast<float>(cell->level()) /
                              static_cast<float>(n_levels)) *
                             std::max(x_dimension, y_dimension);
                }

              const double distance_to_camera =
                std::hypot(point[0] - camera_position[0],
                           point[1] - camera_position[1],
                           point[2] - camera_position[2]);
              const double distance_factor =
                distance_to_camera / (2. * std::max(x_dimension, y_dimension));

              projection_decomposition = svg_project_point(point,
                                                           camera_position,
                                                           camera_direction,
                                                           camera_horizontal,
                                                           camera_focus);

              const unsigned int font_size_this_cell =
                static_cast<unsigned int>(
                  .5 +
                  cell_label_font_size *
                    std::pow(.5, cell->level() - 4. + 3.5 * distance_factor));

              out << "  <text"
                  << " x=\""
                  << static_cast<unsigned int>(
                       .5 +
                       ((projection_decomposition[0] - x_min_perspective) /
                        x_dimension_perspective) *
                         (width - (width / 100.) * 2. * margin_in_percent) +
                       ((width / 100.) * margin_in_percent))
                  << "\" y=\""
                  << static_cast<unsigned int>(
                       .5 + height - (height / 100.) * margin_in_percent -
                       ((projection_decomposition[1] - y_min_perspective) /
                        y_dimension_perspective) *
                         (height - (height / 100.) * 2. * margin_in_percent) +
                       0.5 * font_size_this_cell)
                  << "\" style=\"font-size:" << font_size_this_cell << "px\">";

              if (svg_flags.label_level_number)
                {
                  out << cell->level();
                }

              if (svg_flags.label_cell_index)
                {
                  if (svg_flags.label_level_number)
                    out << '.';
                  out << cell->index();
                }

              if (svg_flags.label_material_id)
                {
                  if (svg_flags.label_level_number ||
                      svg_flags.label_cell_index)
                    out << ',';
                  out << static_cast<std::make_signed_t<types::material_id>>(
                    cell->material_id());
                }

              if (svg_flags.label_subdomain_id)
                {
                  if (svg_flags.label_level_number ||
                      svg_flags.label_cell_index || svg_flags.label_material_id)
                    out << ',';
                  if (cell->is_active())
                    out << static_cast<std::make_signed_t<types::subdomain_id>>(
                      cell->subdomain_id());
                  else
                    out << 'X';
                }

              if (svg_flags.label_level_subdomain_id)
                {
                  if (svg_flags.label_level_number ||
                      svg_flags.label_cell_index ||
                      svg_flags.label_material_id ||
                      svg_flags.label_subdomain_id)
                    out << ',';
                  out << static_cast<std::make_signed_t<types::subdomain_id>>(
                    cell->level_subdomain_id());
                }

              out << "</text>" << '\n';
            }

          // if the current cell lies at the boundary of the triangulation, draw
          // the additional boundary line
          if (svg_flags.boundary_line_thickness != 0u)
            {
              for (auto faceIndex : cell->face_indices())
                {
                  if (cell->at_boundary(faceIndex))
                    {
                      point[0] = cell->face(faceIndex)->vertex(0)[0];
                      point[1] = cell->face(faceIndex)->vertex(0)[1];
                      point[2] = 0;

                      if (svg_flags.convert_level_number_to_height)
                        {
                          point[2] = svg_flags.level_height_factor *
                                     (static_cast<float>(cell->level()) /
                                      static_cast<float>(n_levels)) *
                                     std::max(x_dimension, y_dimension);
                        }

                      projection_decomposition =
                        svg_project_point(point,
                                          camera_position,
                                          camera_direction,
                                          camera_horizontal,
                                          camera_focus);

                      out << "  <line x1=\""
                          << static_cast<unsigned int>(
                               .5 +
                               ((projection_decomposition[0] -
                                 x_min_perspective) /
                                x_dimension_perspective) *
                                 (width -
                                  (width / 100.) * 2. * margin_in_percent) +
                               ((width / 100.) * margin_in_percent))
                          << "\" y1=\""
                          << static_cast<unsigned int>(
                               .5 + height -
                               (height / 100.) * margin_in_percent -
                               ((projection_decomposition[1] -
                                 y_min_perspective) /
                                y_dimension_perspective) *
                                 (height -
                                  (height / 100.) * 2. * margin_in_percent));

                      point[0] = cell->face(faceIndex)->vertex(1)[0];
                      point[1] = cell->face(faceIndex)->vertex(1)[1];
                      point[2] = 0;

                      if (svg_flags.convert_level_number_to_height)
                        {
                          point[2] = svg_flags.level_height_factor *
                                     (static_cast<float>(cell->level()) /
                                      static_cast<float>(n_levels)) *
                                     std::max(x_dimension, y_dimension);
                        }

                      projection_decomposition =
                        svg_project_point(point,
                                          camera_position,
                                          camera_direction,
                                          camera_horizontal,
                                          camera_focus);

                      out << "\" x2=\""
                          << static_cast<unsigned int>(
                               .5 +
                               ((projection_decomposition[0] -
                                 x_min_perspective) /
                                x_dimension_perspective) *
                                 (width -
                                  (width / 100.) * 2. * margin_in_percent) +
                               ((width / 100.) * margin_in_percent))
                          << "\" y2=\""
                          << static_cast<unsigned int>(
                               .5 + height -
                               (height / 100.) * margin_in_percent -
                               ((projection_decomposition[1] -
                                 y_min_perspective) /
                                y_dimension_perspective) *
                                 (height -
                                  (height / 100.) * 2. * margin_in_percent))
                          << "\"/>" << '\n';


                      if (svg_flags.label_boundary_id)
                        {
                          const double distance_to_camera =
                            std::hypot(point[0] - camera_position[0],
                                       point[1] - camera_position[1],
                                       point[2] - camera_position[2]);
                          const double distance_factor =
                            distance_to_camera /
                            (2. * std::max(x_dimension, y_dimension));

                          const unsigned int font_size_this_edge =
                            static_cast<unsigned int>(
                              .5 + .5 * cell_label_font_size *
                                     std::pow(.5,
                                              cell->level() - 4. +
                                                3.5 * distance_factor));

                          point[0] = cell->face(faceIndex)->center()[0];
                          point[1] = cell->face(faceIndex)->center()[1];
                          point[2] = 0;

                          if (svg_flags.convert_level_number_to_height)
                            {
                              point[2] = svg_flags.level_height_factor *
                                         (static_cast<float>(cell->level()) /
                                          static_cast<float>(n_levels)) *
                                         std::max(x_dimension, y_dimension);
                            }

                          projection_decomposition =
                            svg_project_point(point,
                                              camera_position,
                                              camera_direction,
                                              camera_horizontal,
                                              camera_focus);

                          const unsigned int xc = static_cast<unsigned int>(
                            .5 +
                            ((projection_decomposition[0] - x_min_perspective) /
                             x_dimension_perspective) *
                              (width -
                               (width / 100.) * 2. * margin_in_percent) +
                            ((width / 100.) * margin_in_percent));
                          const unsigned int yc = static_cast<unsigned int>(
                            .5 + height - (height / 100.) * margin_in_percent -
                            ((projection_decomposition[1] - y_min_perspective) /
                             y_dimension_perspective) *
                              (height -
                               (height / 100.) * 2. * margin_in_percent));

                          out << "    <circle cx=\"" << xc << "\" cy=\"" << yc
                              << "\" r=\"" << font_size_this_edge << "\" />"
                              << '\n';

                          out << "    <text x=\"" << xc << "\" y=\"" << yc
                              << "\" style=\"font-size:" << font_size_this_edge
                              << "px\" dominant-baseline=\"middle\">"
                              << static_cast<int>(
                                   cell->face(faceIndex)->boundary_id())
                              << "</text>" << '\n';
                        }
                    }
                }
            }
        }
    }



  // draw the legend
  if (svg_flags.draw_legend)
    out << '\n' << " <!-- legend -->" << '\n';

  additional_width = 0;
  if (!svg_flags.margin)
    additional_width = static_cast<unsigned int>(.5 + (height / 100.) * 2.5);

  // explanation of the cell labeling
  if (svg_flags.draw_legend &&
      (svg_flags.label_level_number || svg_flags.label_cell_index ||
       svg_flags.label_material_id || svg_flags.label_subdomain_id ||
       svg_flags.label_level_subdomain_id || svg_flags.label_boundary_id))
    {
      unsigned int line_offset = 0;
      out << " <rect x=\"" << width + additional_width << "\" y=\""
          << static_cast<unsigned int>(.5 + (height / 100.) * margin_in_percent)
          << "\" width=\""
          << static_cast<unsigned int>(.5 + (height / 100.) *
                                              (40. - margin_in_percent))
          << "\" height=\"" << static_cast<unsigned int>(.5 + height * .215)
          << "\"/>" << '\n';

      out << " <text x=\""
          << width + additional_width +
               static_cast<unsigned int>(.5 + (height / 100.) * 1.25)
          << "\" y=\""
          << static_cast<unsigned int>(.5 +
                                       (height / 100.) * margin_in_percent +
                                       (++line_offset) * 1.5 * font_size)
          << "\" style=\"text-anchor:start; font-weight:bold; font-size:"
          << font_size << "px\">"
          << "cell label"
          << "</text>" << '\n';

      if (svg_flags.label_level_number)
        {
          out << "  <text x=\""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "cell_level";

          if (svg_flags.label_cell_index || svg_flags.label_material_id ||
              svg_flags.label_subdomain_id ||
              svg_flags.label_level_subdomain_id)
            out << '.';

          out << "</text>" << '\n';
        }

      if (svg_flags.label_cell_index)
        {
          out << "  <text x=\""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "cell_index";

          if (svg_flags.label_material_id || svg_flags.label_subdomain_id ||
              svg_flags.label_level_subdomain_id)
            out << ',';

          out << "</text>" << '\n';
        }

      if (svg_flags.label_material_id)
        {
          out << "  <text x=\""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "material_id";

          if (svg_flags.label_subdomain_id ||
              svg_flags.label_level_subdomain_id)
            out << ',';

          out << "</text>" << '\n';
        }

      if (svg_flags.label_subdomain_id)
        {
          out << "  <text x= \""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "subdomain_id";

          if (svg_flags.label_level_subdomain_id)
            out << ',';

          out << "</text>" << '\n';
        }

      if (svg_flags.label_level_subdomain_id)
        {
          out << "  <text x= \""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "level_subdomain_id"
              << "</text>" << '\n';
        }

      if (svg_flags.label_boundary_id)
        {
          out << " <text x=\""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 1.25)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-weight:bold; font-size:"
              << font_size << "px\">"
              << "edge label"
              << "</text>" << '\n';

          out << "  <text x= \""
              << width + additional_width +
                   static_cast<unsigned int>(.5 + (height / 100.) * 2.)
              << "\" y=\""
              << static_cast<unsigned int>(.5 +
                                           (height / 100.) * margin_in_percent +
                                           (++line_offset) * 1.5 * font_size)
              << "\" style=\"text-anchor:start; font-style:oblique; font-size:"
              << font_size << "px\">"
              << "boundary_id"
              << "</text>" << '\n';
        }
    }

  // show azimuth angle and polar angle as text below the explanation of the
  // cell labeling
  if (svg_flags.draw_legend)
    {
      out << "  <text x=\"" << width + additional_width << "\" y=\""
          << static_cast<unsigned int>(
               .5 + (height / 100.) * margin_in_percent + 13.75 * font_size)
          << "\" style=\"text-anchor:start; font-size:" << font_size << "px\">"
          << "azimuth: " << svg_flags.azimuth_angle
          << "°, polar: " << svg_flags.polar_angle << "°</text>" << '\n';
    }


  // draw the colorbar
  if (svg_flags.draw_colorbar && (svg_flags.coloring != 0u))
    {
      out << '\n' << " <!-- colorbar -->" << '\n';

      out << " <text x=\"" << width + additional_width << "\" y=\""
          << static_cast<unsigned int>(
               .5 + (height / 100.) * (margin_in_percent + 29.) -
               (font_size / 1.25))
          << "\" style=\"text-anchor:start; font-weight:bold; font-size:"
          << font_size << "px\">";

      switch (svg_flags.coloring)
        {
          case 1:
            out << "material_id";
            break;
          case 2:
            out << "level_number";
            break;
          case 3:
            out << "subdomain_id";
            break;
          case 4:
            out << "level_subdomain_id";
            break;
          default:
            break;
        }

      out << "</text>" << '\n';

      unsigned int element_height = static_cast<unsigned int>(
        ((height / 100.) * (71. - 2. * margin_in_percent)) / n);
      unsigned int element_width =
        static_cast<unsigned int>(.5 + (height / 100.) * 2.5);

      int  labeling_index      = 0;
      auto materials_it        = materials.begin();
      auto levels_it           = levels.begin();
      auto subdomains_it       = subdomains.begin();
      auto level_subdomains_it = level_subdomains.begin();

      for (unsigned int index = 0; index < n; ++index)
        {
          switch (svg_flags.coloring)
            {
              case GridOutFlags::Svg::material_id:
                labeling_index = *materials_it++;
                break;
              case GridOutFlags::Svg::level_number:
                labeling_index = *levels_it++;
                break;
              case GridOutFlags::Svg::subdomain_id:
                labeling_index = *subdomains_it++;
                break;
              case GridOutFlags::Svg::level_subdomain_id:
                labeling_index = *level_subdomains_it++;
                break;
              default:
                break;
            }

          out << "  <rect class=\"r" << labeling_index << "\" x=\""
              << width + additional_width << "\" y=\""
              << static_cast<unsigned int>(.5 + (height / 100.) *
                                                  (margin_in_percent + 29)) +
                   (n - index - 1) * element_height
              << "\" width=\"" << element_width << "\" height=\""
              << element_height << "\"/>" << '\n';

          out << "  <text x=\""
              << width + additional_width + 1.5 * element_width << "\" y=\""
              << static_cast<unsigned int>(.5 + (height / 100.) *
                                                  (margin_in_percent + 29)) +
                   (n - index - 1 + .5) * element_height +
                   static_cast<unsigned int>(.5 + font_size * .35)
              << "\""
              << " style=\"text-anchor:start; font-size:"
              << static_cast<unsigned int>(.5 + font_size) << "px";

          if (index == 0 || index == n - 1)
            out << "; font-weight:bold";

          out << "\">" << labeling_index;

          if (index == n - 1)
            out << " max";
          if (index == 0)
            out << " min";

          out << "</text>" << '\n';

          ++labeling_index;
        }
    }


  // finalize the svg file
  out << '\n' << "</svg>";
  out.flush();
}



template <>
void
GridOut::write_mathgl(const Triangulation<1> &, std::ostream &) const
{
  // 1d specialization not done yet
  DEAL_II_NOT_IMPLEMENTED();
}



template <int dim, int spacedim>
void
GridOut::write_mathgl(const Triangulation<dim, spacedim> &tria,
                      std::ostream                       &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // (i) write header
  {
    // block this to have local variables destroyed after use
    const std::time_t time1 = std::time(nullptr);
    const std::tm    *time  = std::localtime(&time1);

    out
      << "\n#"
      << "\n# This file was generated by the deal.II library."
      << "\n#   Date =  " << time->tm_year + 1900 << "/" << std::setfill('0')
      << std::setw(2) << time->tm_mon + 1 << "/" << std::setfill('0')
      << std::setw(2) << time->tm_mday << "\n#   Time =  " << std::setfill('0')
      << std::setw(2) << time->tm_hour << ":" << std::setfill('0')
      << std::setw(2) << time->tm_min << ":" << std::setfill('0')
      << std::setw(2) << time->tm_sec << "\n#"
      << "\n# For a description of the MathGL script format see the MathGL manual.  "
      << "\n#"
      << "\n# Note: This file is understood by MathGL v2.1 and higher only, and can "
      << "\n#       be quickly viewed in a graphical environment using \'mglview\'. "
      << "\n#" << '\n';
  }

  // define a helper to keep loops approximately dim-independent
  // since MathGL labels axes as x, y, z
  const std::string axes = "xyz";

  // (ii) write preamble and graphing tweaks
  out << "\n#"
      << "\n#   Preamble."
      << "\n#" << '\n';

  if (mathgl_flags.draw_bounding_box)
    out << "\nbox";

  // deal with dimension dependent preamble; eg. default sizes and
  // views for MathGL (cf. gnuplot).
  switch (dim)
    {
      case 2:
        out << "\nsetsize 800 800";
        out << "\nrotate 0 0";
        break;
      case 3:
        out << "\nsetsize 800 800";
        out << "\nrotate 60 40";
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }
  out << '\n';

  // (iii) write vertex ordering
  out << "\n#"
      << "\n#   Vertex ordering."
      << "\n#   list <vertex order> <vertex indices>"
      << "\n#" << '\n';

  // todo: This denotes the natural ordering of vertices, but it needs
  // to check this is really always true for a given grid (it's not
  // true in step-1 grid-2 for instance).
  switch (dim)
    {
      case 2:
        out << "\nlist f 0 1 2 3" << '\n';
        break;
      case 3:
        out
          << "\nlist f 0 2 4 6 | 1 3 5 7 | 0 4 1 5 | 2 6 3 7 | 0 1 2 3 | 4 5 6 7"
          << '\n';
        break;
      default:
        DEAL_II_NOT_IMPLEMENTED();
    }

  // (iv) write a list of vertices of cells
  out << "\n#"
      << "\n#   List of vertices."
      << "\n#   list <id> <vertices>"
      << "\n#" << '\n';

  // run over all active cells and write out a list of
  // xyz-coordinates that correspond to vertices
  // No global indices in deal.II, so we make one up here.
  for (const auto &cell : tria.active_cell_iterators())
    {
      for (unsigned int i = 0; i < dim; ++i)
        {
          // if (cell->direction_flag ()==true)
          //   out << "\ntrue";
          // else
          //   out << "\nfalse";

          out << "\nlist " << axes[i] << cell->active_cell_index() << " ";
          for (const unsigned int j : GeometryInfo<dim>::vertex_indices())
            out << cell->vertex(j)[i] << " ";
        }
      out << '\n';
    }

  // (v) write out cells to plot as quadplot objects
  out << "\n#"
      << "\n#   List of cells to quadplot."
      << "\n#   quadplot <vertex order> <id> <style>"
      << "\n#" << '\n';
  for (unsigned int i = 0; i < tria.n_active_cells(); ++i)
    {
      out << "\nquadplot f ";
      for (unsigned int j = 0; j < dim; ++j)
        out << axes[j] << i << " ";
      out << "\'k#\'";
    }
  out << '\n';

  // (vi) write footer
  out << "\n#"
      << "\n#"
      << "\n#" << '\n';

  // make sure everything now gets to the output stream
  out.flush();
  AssertThrow(out.fail() == false, ExcIO());
}



namespace
{
  /**
   * A function that is able to convert each cell of a triangulation into
   * a patch that can then be output by the functions in DataOutBase.
   * This is made particularly simple because the patch only needs to
   * contain geometry info and additional properties of cells
   */
  template <int dim, int spacedim, typename IteratorType>
  void
  generate_triangulation_patches(
    std::vector<DataOutBase::Patch<dim, spacedim>> &patches,
    const IteratorType                             &begin,
    const std_cxx20::type_identity_t<IteratorType> &end)
  {
    // convert each of the active cells into a patch
    for (auto cell = begin; cell != end; ++cell)
      {
        DataOutBase::Patch<dim, spacedim> patch;
        patch.reference_cell = cell->reference_cell();
        patch.n_subdivisions = 1;
        patch.data.reinit(5, cell->n_vertices());

        for (const unsigned int v : cell->vertex_indices())
          {
            patch.vertices[v] = cell->vertex(v);
            patch.data(0, v)  = cell->level();
            patch.data(1, v) =
              static_cast<std::make_signed_t<types::manifold_id>>(
                cell->manifold_id());
            patch.data(2, v) =
              static_cast<std::make_signed_t<types::material_id>>(
                cell->material_id());
            if (cell->is_active())
              patch.data(3, v) =
                static_cast<std::make_signed_t<types::subdomain_id>>(
                  cell->subdomain_id());
            else
              patch.data(3, v) = -1;
            patch.data(4, v) =
              static_cast<std::make_signed_t<types::subdomain_id>>(
                cell->level_subdomain_id());
          }
        patches.push_back(patch);
      }
  }



  std::vector<std::string>
  triangulation_patch_data_names()
  {
    std::vector<std::string> v(5);
    v[0] = "level";
    v[1] = "manifold";
    v[2] = "material";
    v[3] = "subdomain";
    v[4] = "level_subdomain";
    return v;
  }

  /**
   * Return all boundary lines of non-internal faces in three dimension.
   */
  std::vector<typename Triangulation<3, 3>::active_line_iterator>
  get_boundary_edge_iterators(const Triangulation<3, 3> &tria)
  {
    std::vector<typename Triangulation<3, 3>::active_line_iterator> res;

    std::vector<bool> flags;
    tria.save_user_flags_line(flags);
    const_cast<Triangulation<3, 3> &>(tria).clear_user_flags_line();

    for (auto face : tria.active_face_iterators())
      for (const auto l : face->line_indices())
        {
          const auto line = face->line(l);
          if (line->user_flag_set() || line->has_children())
            continue;
          else
            line->set_user_flag();
          if (line->at_boundary())
            res.emplace_back(line);
        }
    const_cast<Triangulation<3, 3> &>(tria).load_user_flags_line(flags);
    return res;
  }



  /**
   * Same as above, for 1 and 2 dimensional grids. Does nothing.
   */
  template <int dim, int spacedim>
  std::vector<typename Triangulation<dim, spacedim>::active_line_iterator>
  get_boundary_edge_iterators(const Triangulation<dim, spacedim> &)
  {
    return {};
  }



  /**
   * Return all lines of a face in three dimension that have a non-standard
   * boundary indicator (!=0), or a non-flat manifold indicator.
   */
  std::vector<typename Triangulation<3, 3>::active_line_iterator>
  get_relevant_edge_iterators(const Triangulation<3, 3> &tria)
  {
    std::vector<typename Triangulation<3, 3>::active_line_iterator> res;

    std::vector<bool> flags;
    tria.save_user_flags_line(flags);
    const_cast<Triangulation<3, 3> &>(tria).clear_user_flags_line();

    for (auto face : tria.active_face_iterators())
      for (const auto l : face->line_indices())
        {
          const auto line = face->line(l);
          if (line->user_flag_set() || line->has_children())
            continue;
          else
            line->set_user_flag();
          if (line->manifold_id() != numbers::flat_manifold_id ||
              (line->boundary_id() != 0 &&
               line->boundary_id() != numbers::invalid_boundary_id))
            res.emplace_back(line);
        }
    const_cast<Triangulation<3, 3> &>(tria).load_user_flags_line(flags);
    return res;
  }


  /**
   * Same as above, for 1 and 2 dimensional grids. Does nothing.
   */
  template <int dim, int spacedim>
  std::vector<typename Triangulation<dim, spacedim>::active_line_iterator>
  get_relevant_edge_iterators(const Triangulation<dim, spacedim> &)
  {
    return {};
  }



  /**
   * Return all boundary faces of a triangulation.
   */
  template <int dim, int spacedim>
  std::vector<typename Triangulation<dim, spacedim>::active_face_iterator>
  get_boundary_face_iterators(const Triangulation<dim, spacedim> &tria)
  {
    std::vector<typename Triangulation<dim, spacedim>::active_face_iterator>
      res;
    if (dim == 1)
      return res;
    for (auto face : tria.active_face_iterators())
      {
        if (face->boundary_id() != numbers::invalid_boundary_id)
          res.push_back(face);
      }
    return res;
  }



  /**
   * Return all faces of a triangulation that have a non-standard
   * boundary indicator (!=0), or a non-flat manifold indicator.
   */
  template <int dim, int spacedim>
  std::vector<typename Triangulation<dim, spacedim>::active_face_iterator>
  get_relevant_face_iterators(const Triangulation<dim, spacedim> &tria)
  {
    std::vector<typename Triangulation<dim, spacedim>::active_face_iterator>
      res;
    if (dim == 1)
      return res;
    for (auto face : tria.active_face_iterators())
      {
        if (face->manifold_id() != numbers::flat_manifold_id ||
            (face->boundary_id() != 0 &&
             face->boundary_id() != numbers::invalid_boundary_id))
          res.push_back(face);
      }
    return res;
  }
} // namespace



template <int dim, int spacedim>
void
GridOut::write_vtk(const Triangulation<dim, spacedim> &tria,
                   std::ostream                       &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // get the positions of the vertices
  const std::vector<Point<spacedim>> &vertices = tria.get_vertices();

  const auto n_vertices = vertices.size();

  out << "# vtk DataFile Version 3.0\n"
      << "Triangulation generated with deal.II\n"
      << "ASCII\n"
      << "DATASET UNSTRUCTURED_GRID\n"
      << "POINTS " << n_vertices << " double\n";

  // actually write the vertices.
  for (const auto &v : vertices)
    {
      out << v;
      for (unsigned int d = spacedim + 1; d <= 3; ++d)
        out << " 0"; // fill with zeroes
      out << '\n';
    }

  const auto faces = vtk_flags.output_only_relevant ?
                       get_relevant_face_iterators(tria) :
                       get_boundary_face_iterators(tria);
  const auto edges = vtk_flags.output_only_relevant ?
                       get_relevant_edge_iterators(tria) :
                       get_boundary_edge_iterators(tria);

  AssertThrow(
    vtk_flags.output_cells || (dim >= 2 && vtk_flags.output_faces) ||
      (dim >= 3 && vtk_flags.output_edges),
    ExcMessage(
      "At least one of the flags (output_cells, output_faces, output_edges) has to be enabled!"));

  // Write cells preamble
  const int n_cells = (vtk_flags.output_cells ? tria.n_active_cells() : 0) +
                      (vtk_flags.output_faces ? faces.size() : 0) +
                      (vtk_flags.output_edges ? edges.size() : 0);

  // VTK now expects a number telling the total storage requirement to read all
  // cell connectivity information. The connectivity information is read cell by
  // cell, first specifying how many vertices are required to describe the cell,
  // and then specifying the index of every vertex. This means that for every
  // deal.II object type, we always need n_vertices + 1 integer per cell.
  // Compute the total number here.
  int cells_size = 0;

  if (vtk_flags.output_cells)
    for (const auto &cell : tria.active_cell_iterators())
      cells_size += cell->n_vertices() + 1;

  if (vtk_flags.output_faces)
    for (const auto &face : faces)
      cells_size += face->n_vertices() + 1;

  if (vtk_flags.output_edges)
    for (const auto &edge : edges)
      cells_size += edge->n_vertices() + 1;

  AssertThrow(cells_size > 0, ExcMessage("No cells given to be output!"));

  out << "\nCELLS " << n_cells << ' ' << cells_size << '\n';
  /*
   * VTK cells:
   *
   *  1 VTK_VERTEX
   *  3 VTK_LINE
   *  5 VTK_TRIANGLE
   *  9 VTK_QUAD
   * 10 VTK_TETRA
   * 14 VTK_PYRAMID
   * 13 VTK_WEDGE
   * 12 VTK_HEXAHEDRON
   *
   * see also: https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
   */
  static const std::array<int, 8> deal_to_vtk_cell_type = {
    {1, 3, 5, 9, 10, 14, 13, 12}};
  static const std::array<unsigned int, 8> vtk_to_deal_hypercube = {
    {0, 1, 3, 2, 4, 5, 7, 6}};

  // write cells.
  if (vtk_flags.output_cells)
    for (const auto &cell : tria.active_cell_iterators())
      {
        out << cell->n_vertices();
        for (const unsigned int i : cell->vertex_indices())
          {
            out << ' ';
            const auto reference_cell = cell->reference_cell();

            if ((reference_cell == ReferenceCells::Vertex) ||
                (reference_cell == ReferenceCells::Line) ||
                (reference_cell == ReferenceCells::Quadrilateral) ||
                (reference_cell == ReferenceCells::Hexahedron))
              out << cell->vertex_index(vtk_to_deal_hypercube[i]);
            else if ((reference_cell == ReferenceCells::Triangle) ||
                     (reference_cell == ReferenceCells::Tetrahedron) ||
                     (reference_cell == ReferenceCells::Wedge))
              out << cell->vertex_index(i);
            else if (reference_cell == ReferenceCells::Pyramid)
              {
                static const std::array<unsigned int, 5> permutation_table{
                  {0, 1, 3, 2, 4}};
                out << cell->vertex_index(permutation_table[i]);
              }
            else
              DEAL_II_NOT_IMPLEMENTED();
          }
        out << '\n';
      }
  if (vtk_flags.output_faces)
    for (const auto &face : faces)
      {
        out << face->n_vertices();
        for (const unsigned int i : face->vertex_indices())
          {
            out << ' '
                << face->vertex_index(GeometryInfo<dim>::vertices_per_face ==
                                          face->n_vertices() ?
                                        vtk_to_deal_hypercube[i] :
                                        i);
          }
        out << '\n';
      }
  if (vtk_flags.output_edges)
    for (const auto &edge : edges)
      {
        out << 2;
        for (const unsigned int i : edge->vertex_indices())
          out << ' ' << edge->vertex_index(i);
        out << '\n';
      }

  // write cell types
  out << "\nCELL_TYPES " << n_cells << '\n';
  if (vtk_flags.output_cells)
    {
      for (const auto &cell : tria.active_cell_iterators())
        out << deal_to_vtk_cell_type[static_cast<int>(cell->reference_cell())]
            << ' ';
      out << '\n';
    }
  if (vtk_flags.output_faces)
    {
      for (const auto &face : faces)
        out << deal_to_vtk_cell_type[static_cast<int>(face->reference_cell())]
            << ' ';
      out << '\n';
    }
  if (vtk_flags.output_edges)
    {
      for (const auto &edge : edges)
        out << deal_to_vtk_cell_type[static_cast<int>(edge->reference_cell())]
            << ' ';
    }
  out << "\n\nCELL_DATA " << n_cells << '\n'
      << "SCALARS MaterialID int 1\n"
      << "LOOKUP_TABLE default\n";

  // Now material id and boundary id
  if (vtk_flags.output_cells)
    {
      for (const auto &cell : tria.active_cell_iterators())
        {
          out << static_cast<std::make_signed_t<types::material_id>>(
                   cell->material_id())
              << ' ';
        }
      out << '\n';
    }
  if (vtk_flags.output_faces)
    {
      for (const auto &face : faces)
        {
          out << static_cast<std::make_signed_t<types::boundary_id>>(
                   face->boundary_id())
              << ' ';
        }
      out << '\n';
    }
  if (vtk_flags.output_edges)
    {
      for (const auto &edge : edges)
        {
          out << static_cast<std::make_signed_t<types::boundary_id>>(
                   edge->boundary_id())
              << ' ';
        }
    }

  out << "\n\nSCALARS ManifoldID int 1\n"
      << "LOOKUP_TABLE default\n";

  // Now manifold id
  if (vtk_flags.output_cells)
    {
      for (const auto &cell : tria.active_cell_iterators())
        {
          out << static_cast<std::make_signed_t<types::manifold_id>>(
                   cell->manifold_id())
              << ' ';
        }
      out << '\n';
    }
  if (vtk_flags.output_faces)
    {
      for (const auto &face : faces)
        {
          out << static_cast<std::make_signed_t<types::manifold_id>>(
                   face->manifold_id())
              << ' ';
        }
      out << '\n';
    }
  if (vtk_flags.output_edges)
    {
      for (const auto &edge : edges)
        {
          out << static_cast<std::make_signed_t<types::manifold_id>>(
                   edge->manifold_id())
              << ' ';
        }
      out << '\n';
    }

  out.flush();

  AssertThrow(out.fail() == false, ExcIO());
}



template <int dim, int spacedim>
void
GridOut::write_vtu(const Triangulation<dim, spacedim> &tria,
                   std::ostream                       &out) const
{
  AssertThrow(out.fail() == false, ExcIO());

  // convert the cells of the triangulation into a set of patches
  // and then have them output. since there is no data attached to
  // the geometry, we also do not have to provide any names, identifying
  // information, etc.
  std::vector<DataOutBase::Patch<dim, spacedim>> patches;
  patches.reserve(tria.n_active_cells());
  generate_triangulation_patches(patches, tria.begin_active(), tria.end());

  DataOutBase::write_vtu_header(out, vtu_flags);
  DataOutBase::write_vtu_main(
    patches,
    triangulation_patch_data_names(),
    std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>(),
    vtu_flags,
    out);
  if (vtu_flags.serialize_triangulation)
    {
      out << " </UnstructuredGrid>\n";
      out << "<dealiiData  encoding=\"base64\">";
      std::stringstream               outstring;
      boost::archive::binary_oarchive ia(outstring);
      tria.save(ia, 0);
      const auto compressed = Utilities::compress(outstring.str());
      out << Utilities::encode_base64({compressed.begin(), compressed.end()});
      out << "\n</dealiiData>\n";
      out << "</VTKFile>\n";
    }
  else
    DataOutBase::write_vtu_footer(out);

  out << std::flush;
  AssertThrow(out.fail() == false, ExcIO());
}



template <int dim, int spacedim>
void
GridOut::write_mesh_per_processor_as_vtu(
  const Triangulation<dim, spacedim> &tria,
  const std::string                  &filename_without_extension,
  const bool                          view_levels,
  const bool                          include_artificial) const
{
  std::vector<DataOutBase::Patch<dim, spacedim>> patches;
  const unsigned int                             n_datasets = 4;
  std::vector<std::string>                       data_names;
  data_names.emplace_back("level");
  data_names.emplace_back("subdomain");
  data_names.emplace_back("level_subdomain");
  data_names.emplace_back("proc_writing");

  const auto &reference_cells = tria.get_reference_cells();

  AssertDimension(reference_cells.size(), 1);

  const auto &reference_cell = reference_cells[0];

  const unsigned int n_q_points = reference_cell.n_vertices();

  for (const auto &cell : tria.cell_iterators())
    {
      if (!view_levels)
        {
          if (cell->has_children())
            continue;
          if (!include_artificial &&
              cell->subdomain_id() == numbers::artificial_subdomain_id)
            continue;
        }
      else if (!include_artificial)
        {
          if (cell->has_children() &&
              cell->level_subdomain_id() == numbers::artificial_subdomain_id)
            continue;
          else if (cell->is_active() &&
                   cell->level_subdomain_id() ==
                     numbers::artificial_subdomain_id &&
                   cell->subdomain_id() == numbers::artificial_subdomain_id)
            continue;
        }

      DataOutBase::Patch<dim, spacedim> patch;
      patch.data.reinit(n_datasets, n_q_points);
      patch.points_are_available = false;
      patch.reference_cell       = reference_cell;

      for (unsigned int vertex = 0; vertex < n_q_points; ++vertex)
        {
          patch.vertices[vertex] = cell->vertex(vertex);
          patch.data(0, vertex)  = cell->level();
          if (cell->is_active())
            patch.data(1, vertex) = static_cast<double>(
              static_cast<std::make_signed_t<types::subdomain_id>>(
                cell->subdomain_id()));
          else
            patch.data(1, vertex) = -1.0;
          patch.data(2, vertex) = static_cast<double>(
            static_cast<std::make_signed_t<types::subdomain_id>>(
              cell->level_subdomain_id()));
          patch.data(3, vertex) = tria.locally_owned_subdomain();
        }

      for (auto f : reference_cell.face_indices())
        patch.neighbors[f] = numbers::invalid_unsigned_int;
      patches.push_back(patch);
    }

  // only create .pvtu file if running in parallel
  // if not, just create a .vtu file with no reference
  // to the processor number
  std::string new_file = filename_without_extension + ".vtu";
  if (const parallel::TriangulationBase<dim, spacedim> *tr =
        dynamic_cast<const parallel::TriangulationBase<dim, spacedim> *>(&tria))
    {
      new_file = filename_without_extension + ".proc" +
                 Utilities::int_to_string(tr->locally_owned_subdomain(), 4) +
                 ".vtu";

      // create .pvtu record
      if (tr->locally_owned_subdomain() == 0)
        {
          std::vector<std::string> filenames;

          // .pvtu needs to reference the files without a relative path because
          // it will be written in the same directory. For this, remove any
          // paths from filename.
          std::size_t pos = filename_without_extension.find_last_of('/');
          if (pos == std::string::npos)
            pos = 0;
          else
            pos += 1;
          const unsigned int n_procs =
            Utilities::MPI::n_mpi_processes(tr->get_mpi_communicator());
          filenames.reserve(n_procs);
          for (unsigned int i = 0; i < n_procs; ++i)
            filenames.push_back(filename_without_extension.substr(pos) +
                                ".proc" + Utilities::int_to_string(i, 4) +
                                ".vtu");

          const std::string pvtu_filename =
            (filename_without_extension + ".pvtu");
          std::ofstream pvtu_output(pvtu_filename);

          DataOut<dim, spacedim> data_out;
          data_out.attach_triangulation(*tr);

          // We need a dummy vector with the names of the data values in the
          // .vtu files in order that the .pvtu contains reference these values
          const Vector<float> dummy_vector(tr->n_active_cells());
          data_out.add_data_vector(dummy_vector, "level");
          data_out.add_data_vector(dummy_vector, "subdomain");
          data_out.add_data_vector(dummy_vector, "level_subdomain");
          data_out.add_data_vector(dummy_vector, "proc_writing");

          data_out.build_patches();

          data_out.write_pvtu_record(pvtu_output, filenames);
        }
    }

  std::ofstream out(new_file);
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    vector_data_ranges;
  DataOutBase::write_vtu(
    patches, data_names, vector_data_ranges, vtu_flags, out);
}



unsigned int
GridOut::n_boundary_faces(const Triangulation<1, 1> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_lines(const Triangulation<1, 1> &) const
{
  return 0;
}


unsigned int
GridOut::n_boundary_faces(const Triangulation<1, 2> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_lines(const Triangulation<1, 2> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_faces(const Triangulation<1, 3> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_lines(const Triangulation<1, 3> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_lines(const Triangulation<2, 2> &) const
{
  return 0;
}

unsigned int
GridOut::n_boundary_lines(const Triangulation<2, 3> &) const
{
  return 0;
}



template <int dim, int spacedim>
unsigned int
GridOut::n_boundary_faces(const Triangulation<dim, spacedim> &tria) const
{
  typename Triangulation<dim, spacedim>::active_face_iterator face, endf;
  unsigned int                                                n_faces = 0;

  for (const auto &face : tria.active_face_iterators())
    if ((face->at_boundary()) && (face->boundary_id() != 0))
      ++n_faces;

  return n_faces;
}



template <int dim, int spacedim>
unsigned int
GridOut::n_boundary_lines(const Triangulation<dim, spacedim> &tria) const
{
  // save the user flags for lines so
  // we can use these flags to track
  // which ones we've already counted
  std::vector<bool> line_flags;
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).save_user_flags_line(
    line_flags);
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
    .clear_user_flags_line();

  unsigned int n_lines = 0;

  for (const auto &cell : tria.active_cell_iterators())
    for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
      if (cell->line(l)->at_boundary() && (cell->line(l)->boundary_id() != 0) &&
          (cell->line(l)->user_flag_set() == false))
        {
          ++n_lines;
          cell->line(l)->set_user_flag();
        }

  // at the end, restore the user
  // flags for the lines
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).load_user_flags_line(
    line_flags);

  return n_lines;
}



unsigned int
GridOut::write_msh_faces(const Triangulation<1, 1> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_msh_faces(const Triangulation<1, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_msh_faces(const Triangulation<1, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_msh_lines(const Triangulation<1, 1> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_msh_lines(const Triangulation<1, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_msh_lines(const Triangulation<1, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_msh_lines(const Triangulation<2, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_msh_lines(const Triangulation<2, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}



template <int dim, int spacedim>
unsigned int
GridOut::write_msh_faces(const Triangulation<dim, spacedim> &tria,
                         const unsigned int                  next_element_index,
                         std::ostream                       &out) const
{
  unsigned int current_element_index = next_element_index;

  for (const auto &face : tria.active_face_iterators())
    if (face->at_boundary() && (face->boundary_id() != 0))
      {
        out << current_element_index << ' '
            << face->reference_cell().gmsh_element_type() << ' ';
        out << static_cast<unsigned int>(face->boundary_id()) << ' '
            << static_cast<unsigned int>(face->boundary_id()) << ' '
            << face->n_vertices();
        // note: vertex numbers are 1-base
        for (const unsigned int vertex : face->vertex_indices())
          {
            if (face->reference_cell() == ReferenceCells::Quadrilateral)
              out << ' '
                  << face->vertex_index(
                       GeometryInfo<dim - 1>::ucd_to_deal[vertex]) +
                       1;
            else if ((face->reference_cell() == ReferenceCells::Triangle) ||
                     (face->reference_cell() == ReferenceCells::Line))
              out << ' ' << face->vertex_index(vertex) + 1;
            else
              DEAL_II_ASSERT_UNREACHABLE();
          }
        out << '\n';

        ++current_element_index;
      }
  return current_element_index;
}



template <int dim, int spacedim>
unsigned int
GridOut::write_msh_lines(const Triangulation<dim, spacedim> &tria,
                         const unsigned int                  next_element_index,
                         std::ostream                       &out) const
{
  unsigned int current_element_index = next_element_index;
  // save the user flags for lines so
  // we can use these flags to track
  // which ones we've already taken
  // care of
  std::vector<bool> line_flags;
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).save_user_flags_line(
    line_flags);
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
    .clear_user_flags_line();

  for (const auto &cell : tria.active_cell_iterators())
    for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
      if (cell->line(l)->at_boundary() && (cell->line(l)->boundary_id() != 0) &&
          (cell->line(l)->user_flag_set() == false))
        {
          out << next_element_index << ' '
              << ReferenceCells::Line.gmsh_element_type() << ' ';
          out << static_cast<unsigned int>(cell->line(l)->boundary_id()) << ' '
              << static_cast<unsigned int>(cell->line(l)->boundary_id())
              << " 2 "; // two vertex indices to follow
          // note: vertex numbers are 1-base
          for (unsigned int vertex = 0; vertex < 2; ++vertex)
            out << ' '
                << cell->line(l)->vertex_index(
                     GeometryInfo<dim - 2>::ucd_to_deal[vertex]) +
                     1;
          out << '\n';

          // move on to the next line
          // but mark the current one
          // as taken care of
          ++current_element_index;
          cell->line(l)->set_user_flag();
        }

  // at the end, restore the user
  // flags for the lines
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).load_user_flags_line(
    line_flags);

  return current_element_index;
}



unsigned int
GridOut::write_ucd_faces(const Triangulation<1, 1> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_ucd_faces(const Triangulation<1, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_ucd_faces(const Triangulation<1, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_ucd_lines(const Triangulation<1, 1> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_ucd_lines(const Triangulation<1, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_ucd_lines(const Triangulation<1, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}


unsigned int
GridOut::write_ucd_lines(const Triangulation<2, 2> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}

unsigned int
GridOut::write_ucd_lines(const Triangulation<2, 3> &,
                         const unsigned int next_element_index,
                         std::ostream &) const
{
  return next_element_index;
}



template <int dim, int spacedim>
unsigned int
GridOut::write_ucd_faces(const Triangulation<dim, spacedim> &tria,
                         const unsigned int                  next_element_index,
                         std::ostream                       &out) const
{
  unsigned int current_element_index = next_element_index;
  typename Triangulation<dim, spacedim>::active_face_iterator face, endf;

  for (const auto &face : tria.active_face_iterators())
    if (face->at_boundary() && (face->boundary_id() != 0))
      {
        out << current_element_index << "  "
            << static_cast<unsigned int>(face->boundary_id()) << "  ";
        switch (dim)
          {
            case 2:
              out << "line    ";
              break;
            case 3:
              out << "quad    ";
              break;
            default:
              DEAL_II_NOT_IMPLEMENTED();
          }
        // note: vertex numbers are 1-base
        for (unsigned int vertex = 0;
             vertex < GeometryInfo<dim>::vertices_per_face;
             ++vertex)
          out << face->vertex_index(
                   GeometryInfo<dim - 1>::ucd_to_deal[vertex]) +
                   1
              << ' ';
        out << '\n';

        ++current_element_index;
      }
  return current_element_index;
}



template <int dim, int spacedim>
unsigned int
GridOut::write_ucd_lines(const Triangulation<dim, spacedim> &tria,
                         const unsigned int                  next_element_index,
                         std::ostream                       &out) const
{
  unsigned int current_element_index = next_element_index;
  // save the user flags for lines so
  // we can use these flags to track
  // which ones we've already taken
  // care of
  std::vector<bool> line_flags;
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).save_user_flags_line(
    line_flags);
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria)
    .clear_user_flags_line();

  for (const auto &cell : tria.active_cell_iterators())
    for (unsigned int l = 0; l < GeometryInfo<dim>::lines_per_cell; ++l)
      if (cell->line(l)->at_boundary() && (cell->line(l)->boundary_id() != 0) &&
          (cell->line(l)->user_flag_set() == false))
        {
          out << current_element_index << "  "
              << static_cast<unsigned int>(cell->line(l)->boundary_id())
              << "  line    ";
          // note: vertex numbers in ucd format are 1-base
          for (unsigned int vertex = 0; vertex < 2; ++vertex)
            out << cell->line(l)->vertex_index(
                     GeometryInfo<dim - 2>::ucd_to_deal[vertex]) +
                     1
                << ' ';
          out << '\n';

          // move on to the next line
          // but mark the current one
          // as taken care of
          ++current_element_index;
          cell->line(l)->set_user_flag();
        }

  // at the end, restore the user
  // flags for the lines
  const_cast<dealii::Triangulation<dim, spacedim> &>(tria).load_user_flags_line(
    line_flags);
  return current_element_index;
}



namespace internal
{
  namespace
  {
    /**
     * GNUPlot output can, optionally, output multiple line segments on each
     * grid line to make curved lines look curved. However, this is very
     * wasteful when the line itself is straight since we do not need multiple
     * line segments to draw a straight line. This function tries to identify
     * whether or not the collection of points corresponds to a straight line:
     * if it does then all but the first and last points can be removed.
     */
    template <int spacedim>
    void
    remove_collinear_points(std::vector<Point<spacedim>> &points)
    {
      while (points.size() > 2)
        {
          Tensor<1, spacedim> first_difference = points[1] - points[0];
          first_difference /= first_difference.norm();
          Tensor<1, spacedim> second_difference = points[2] - points[1];
          second_difference /= second_difference.norm();
          // If the three points are collinear then remove the middle one.
          if ((first_difference - second_difference).norm() < 1e-10)
            points.erase(points.begin() + 1);
          else
            break;
        }
    }



    template <int spacedim>
    void
    write_gnuplot(const dealii::Triangulation<1, spacedim> &tria,
                  std::ostream                             &out,
                  const Mapping<1, spacedim> *,
                  const GridOutFlags::Gnuplot &gnuplot_flags)
    {
      AssertThrow(out.fail() == false, ExcIO());

      for (const auto &cell : tria.active_cell_iterators())
        {
          if (gnuplot_flags.write_cell_numbers)
            out << "# cell " << cell << '\n';

          out << cell->vertex(0) << ' ' << cell->level() << ' '
              << cell->material_id() << '\n'
              << cell->vertex(1) << ' ' << cell->level() << ' '
              << cell->material_id() << '\n'
              << "\n\n";
        }

      // make sure everything now gets to
      // disk
      out.flush();

      AssertThrow(out.fail() == false, ExcIO());
    }



    template <int spacedim>
    void
    write_gnuplot(const dealii::Triangulation<2, spacedim> &tria,
                  std::ostream                             &out,
                  const Mapping<2, spacedim>               *mapping,
                  const GridOutFlags::Gnuplot              &gnuplot_flags)
    {
      AssertThrow(out.fail() == false, ExcIO());

      const int dim = 2;

      const unsigned int n_additional_points =
        gnuplot_flags.n_extra_curved_line_points;
      const unsigned int n_points = 2 + n_additional_points;

      // If we need to plot curved lines then generate a quadrature formula to
      // place points via the mapping
      Quadrature<dim>             q_projector;
      std::vector<Point<dim - 1>> boundary_points;
      if (mapping != nullptr)
        {
          boundary_points.resize(n_points);
          boundary_points[0][0]            = 0;
          boundary_points[n_points - 1][0] = 1;
          for (unsigned int i = 1; i < n_points - 1; ++i)
            boundary_points[i][0] = 1. * i / (n_points - 1);

          const std::vector<double> dummy_weights(n_points, 1. / n_points);
          const Quadrature<dim - 1> quadrature(boundary_points, dummy_weights);

          q_projector = QProjector<dim>::project_to_all_faces(
            ReferenceCells::get_hypercube<dim>(), quadrature);
        }

      static constexpr std::array<unsigned int, 8> local_vertex_numbering = {
        {0, 1, 5, 4, 2, 3, 7, 6}};
      for (const auto &cell : tria.active_cell_iterators())
        {
          if (gnuplot_flags.write_cell_numbers)
            out << "# cell " << cell << '\n';

          if (mapping == nullptr ||
              (dim == spacedim ?
                 (!cell->at_boundary() && !gnuplot_flags.curved_inner_cells) :
                 // ignore checking for boundary or interior cells in the codim
                 // 1 case: 'or false' is a no-op
                 false))
            {
              // write out the four sides of this cell by putting the four
              // points (+ the initial point again) in a row and lifting the
              // drawing pencil at the end
              for (const unsigned int i : GeometryInfo<dim>::vertex_indices())
                out << cell->vertex(dim == 3 ?
                                      local_vertex_numbering[i] :
                                      GeometryInfo<dim>::ucd_to_deal[i])
                    << ' ' << cell->level() << ' ' << cell->material_id()
                    << '\n';
              out << cell->vertex(0) << ' ' << cell->level() << ' '
                  << cell->material_id() << '\n'
                  << '\n' // double new line for gnuplot 3d plots
                  << '\n';
            }
          else
            // cell is at boundary and we are to treat curved boundaries. so
            // loop over all faces and draw them as small pieces of lines
            {
              for (const unsigned int face_no :
                   GeometryInfo<dim>::face_indices())
                {
                  const typename dealii::Triangulation<dim,
                                                       spacedim>::face_iterator
                    face = cell->face(face_no);
                  if (dim != spacedim || face->at_boundary() ||
                      gnuplot_flags.curved_inner_cells)
                    {
                      // Save the points on each face to a vector and then try
                      // to remove collinear points that won't show up in the
                      // generated plot.
                      std::vector<Point<spacedim>> line_points;
                      // compute offset of quadrature points within set of
                      // projected points
                      const auto offset =
                        QProjector<dim>::DataSetDescriptor::face(
                          cell->reference_cell(),
                          face_no,
                          cell->combined_face_orientation(face_no),
                          n_points);
                      line_points.reserve(n_points);
                      for (unsigned int i = 0; i < n_points; ++i)
                        line_points.push_back(
                          mapping->transform_unit_to_real_cell(
                            cell, q_projector.point(offset + i)));
                      internal::remove_collinear_points(line_points);

                      for (const Point<spacedim> &point : line_points)
                        out << point << ' ' << cell->level() << ' '
                            << cell->material_id() << '\n';

                      out << '\n' << '\n';
                    }
                  else
                    {
                      // if, however, the face is not at the boundary and we
                      // don't want to curve anything, then draw it as usual
                      out << face->vertex(0) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n'
                          << face->vertex(1) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n'
                          << '\n'
                          << '\n';
                    }
                }
            }
        }

      // make sure everything now gets to disk
      out.flush();

      AssertThrow(out.fail() == false, ExcIO());
    }



    template <int spacedim>
    void
    write_gnuplot(const dealii::Triangulation<3, spacedim> &tria,
                  std::ostream                             &out,
                  const Mapping<3, spacedim>               *mapping,
                  const GridOutFlags::Gnuplot              &gnuplot_flags)
    {
      AssertThrow(out.fail() == false, ExcIO());

      const int dim = 3;

      const unsigned int n_additional_points =
        gnuplot_flags.n_extra_curved_line_points;
      const unsigned int n_points = 2 + n_additional_points;

      // If we need to plot curved lines then generate a quadrature formula to
      // place points via the mapping
      std::unique_ptr<Quadrature<dim>> q_projector;
      std::vector<Point<1>>            boundary_points;
      if (mapping != nullptr)
        {
          boundary_points.resize(n_points);
          boundary_points[0][0]            = 0;
          boundary_points[n_points - 1][0] = 1;
          for (unsigned int i = 1; i < n_points - 1; ++i)
            boundary_points[i][0] = 1. * i / (n_points - 1);

          const std::vector<double> dummy_weights(n_points, 1. / n_points);
          const Quadrature<1> quadrature1d(boundary_points, dummy_weights);

          // tensor product of points, only one copy
          const QIterated<dim - 1> quadrature(quadrature1d, 1);
          q_projector = std::make_unique<Quadrature<dim>>(
            QProjector<dim>::project_to_all_faces(
              ReferenceCells::get_hypercube<dim>(), quadrature));
        }

      for (const auto &cell : tria.active_cell_iterators())
        {
          if (gnuplot_flags.write_cell_numbers)
            out << "# cell " << cell << '\n';

          if (mapping == nullptr || n_points == 2 ||
              (!cell->has_boundary_lines() &&
               !gnuplot_flags.curved_inner_cells))
            {
              if (cell->reference_cell() == ReferenceCells::Hexahedron)
                {
                  // front face
                  out << cell->vertex(0) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(1) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(5) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(4) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(0) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';
                  // back face
                  out << cell->vertex(2) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(3) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(7) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(6) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(2) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';

                  // now for the four connecting lines
                  out << cell->vertex(0) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(2) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';
                  out << cell->vertex(1) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(3) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';
                  out << cell->vertex(5) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(7) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';
                  out << cell->vertex(4) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << cell->vertex(6) << ' ' << cell->level() << ' '
                      << cell->material_id() << '\n'
                      << '\n';
                }
              else if (cell->reference_cell() == ReferenceCells::Tetrahedron)
                {
                  // Draw the tetrahedron as a two collections of lines.
                  for (const unsigned int v : {0, 1, 2, 0, 3, 2})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of first line

                  for (const unsigned int v : {3, 1})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of second line
                }
              else if (cell->reference_cell() == ReferenceCells::Wedge)
                {
                  // Draw the wedge as a collection of three
                  // lines. The first one wraps around the base,
                  // goes up to the top, and wraps around that. The
                  // second and third are just individual lines
                  // going from base to top.
                  for (const unsigned int v : {0, 1, 2, 0, 3, 4, 5, 3})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of first line

                  for (const unsigned int v : {1, 4})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of second line

                  for (const unsigned int v : {2, 5})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of third line
                }
              else if (cell->reference_cell() == ReferenceCells::Pyramid)
                {
                  // Draw the pyramid as a collections of two lines.
                  for (const unsigned int v : {0, 1, 3, 2, 0, 4, 1})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of first line

                  for (const unsigned int v : {2, 4, 3})
                    {
                      out << cell->vertex(v) << ' ' << cell->level() << ' '
                          << cell->material_id() << '\n';
                    }
                  out << '\n'; // end of second line
                }
              else
                DEAL_II_NOT_IMPLEMENTED();
            }
          else // need to handle curved boundaries
            {
              Assert(cell->reference_cell() == ReferenceCells::Hexahedron,
                     ExcNotImplemented());
              for (const unsigned int face_no :
                   GeometryInfo<dim>::face_indices())
                {
                  const typename dealii::Triangulation<dim,
                                                       spacedim>::face_iterator
                    face = cell->face(face_no);

                  if (face->at_boundary() &&
                      gnuplot_flags.write_additional_boundary_lines)
                    {
                      const auto offset =
                        QProjector<dim>::DataSetDescriptor::face(
                          cell->reference_cell(),
                          face_no,
                          cell->combined_face_orientation(face_no),
                          n_points * n_points);
                      for (unsigned int i = 0; i < n_points - 1; ++i)
                        for (unsigned int j = 0; j < n_points - 1; ++j)
                          {
                            const Point<spacedim> p0 =
                              mapping->transform_unit_to_real_cell(
                                cell,
                                q_projector->point(offset + i * n_points + j));
                            out << p0 << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';
                            out << (mapping->transform_unit_to_real_cell(
                                     cell,
                                     q_projector->point(
                                       offset + (i + 1) * n_points + j)))
                                << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';
                            out << (mapping->transform_unit_to_real_cell(
                                     cell,
                                     q_projector->point(
                                       offset + (i + 1) * n_points + j + 1)))
                                << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';
                            out << (mapping->transform_unit_to_real_cell(
                                     cell,
                                     q_projector->point(offset + i * n_points +
                                                        j + 1)))
                                << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';
                            // and the first point again
                            out << p0 << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';
                            out << '\n' << '\n';
                          }
                    }
                  else
                    {
                      for (unsigned int l = 0;
                           l < GeometryInfo<dim>::lines_per_face;
                           ++l)
                        {
                          const typename dealii::Triangulation<dim, spacedim>::
                            line_iterator line = face->line(l);

                          const Point<spacedim> &v0 = line->vertex(0),
                                                &v1 = line->vertex(1);
                          if (line->at_boundary() ||
                              gnuplot_flags.curved_inner_cells)
                            {
                              // Save the points on each face to a vector and
                              // then try to remove collinear points that won't
                              // show up in the generated plot.
                              std::vector<Point<spacedim>> line_points;
                              // transform_real_to_unit_cell could be replaced
                              // by using QProjector<dim>::project_to_line
                              // which is not yet implemented
                              const Point<spacedim>
                                u0 = mapping->transform_real_to_unit_cell(cell,
                                                                          v0),
                                u1 = mapping->transform_real_to_unit_cell(cell,
                                                                          v1);
                              line_points.reserve(n_points);
                              for (unsigned int i = 0; i < n_points; ++i)
                                line_points.push_back(
                                  mapping->transform_unit_to_real_cell(
                                    cell,
                                    (1 - boundary_points[i][0]) * u0 +
                                      boundary_points[i][0] * u1));
                              internal::remove_collinear_points(line_points);
                              for (const Point<spacedim> &point : line_points)
                                out << point << ' ' << cell->level() << ' '
                                    << static_cast<unsigned int>(
                                         cell->material_id())
                                    << '\n';
                            }
                          else
                            out << v0 << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n'
                                << v1 << ' ' << cell->level() << ' '
                                << cell->material_id() << '\n';

                          out << '\n' << '\n';
                        }
                    }
                }
            }
        }

      // make sure everything now gets to disk
      out.flush();

      AssertThrow(out.fail() == false, ExcIO());
    }
  } // namespace
} // namespace internal



template <int dim, int spacedim>
void
GridOut::write_gnuplot(const Triangulation<dim, spacedim> &tria,
                       std::ostream                       &out,
                       const Mapping<dim, spacedim>       *mapping) const
{
  internal::write_gnuplot(tria, out, mapping, gnuplot_flags);
}



namespace internal
{
  namespace
  {
    struct LineEntry
    {
      Point<2>     first;
      Point<2>     second;
      bool         colorize;
      unsigned int level;
      LineEntry(const Point<2>    &f,
                const Point<2>    &s,
                const bool         c,
                const unsigned int l)
        : first(f)
        , second(s)
        , colorize(c)
        , level(l)
      {}
    };


    void
    write_eps(const dealii::Triangulation<1> &,
              std::ostream &,
              const Mapping<1> *,
              const GridOutFlags::Eps<2> &,
              const GridOutFlags::Eps<3> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

    void
    write_eps(const dealii::Triangulation<1, 2> &,
              std::ostream &,
              const Mapping<1, 2> *,
              const GridOutFlags::Eps<2> &,
              const GridOutFlags::Eps<3> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

    void
    write_eps(const dealii::Triangulation<1, 3> &,
              std::ostream &,
              const Mapping<1, 3> *,
              const GridOutFlags::Eps<2> &,
              const GridOutFlags::Eps<3> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }

    void
    write_eps(const dealii::Triangulation<2, 3> &,
              std::ostream &,
              const Mapping<2, 3> *,
              const GridOutFlags::Eps<2> &,
              const GridOutFlags::Eps<3> &)
    {
      DEAL_II_NOT_IMPLEMENTED();
    }



    template <int dim, int spacedim>
    void
    write_eps(const dealii::Triangulation<dim, spacedim> &tria,
              std::ostream                               &out,
              const Mapping<dim, spacedim>               *mapping,
              const GridOutFlags::Eps<2>                 &eps_flags_2,
              const GridOutFlags::Eps<3>                 &eps_flags_3)
    {
      using LineList = std::list<LineEntry>;

      // We should never get here in 1d since this function is overloaded for
      // all dim == 1 cases.
      Assert(dim == 2 || dim == 3, ExcInternalError());

      // Copy, with an object slice, something containing the flags common to
      // all dimensions in order to avoid the recurring distinctions between
      // the different eps_flags present.
      const GridOutFlags::EpsFlagsBase eps_flags_base =
        dim == 2 ?
          static_cast<const GridOutFlags::EpsFlagsBase &>(eps_flags_2) :
          static_cast<const GridOutFlags::EpsFlagsBase &>(eps_flags_3);

      AssertThrow(out.fail() == false, ExcIO());
      const unsigned int n_points = eps_flags_base.n_boundary_face_points;

      // make up a list of lines by which
      // we will construct the triangulation
      //
      // this part unfortunately is a bit
      // dimension dependent, so we have to
      // treat every dimension different.
      // however, by directly producing
      // the lines to be printed, i.e. their
      // 2d images, we can later do the
      // actual output dimension independent
      // again
      LineList line_list;

      switch (dim)
        {
          case 1:
            {
              DEAL_II_ASSERT_UNREACHABLE();
              break;
            }

          case 2:
            {
              for (const auto &cell : tria.active_cell_iterators())
                for (const unsigned int line_no : cell->line_indices())
                  {
                    typename dealii::Triangulation<dim, spacedim>::line_iterator
                      line = cell->line(line_no);

                    // first treat all
                    // interior lines and
                    // make up a list of
                    // them. if curved
                    // lines shall not be
                    // supported (i.e. no
                    // mapping is
                    // provided), then also
                    // treat all other
                    // lines
                    if (!line->has_children() &&
                        (mapping == nullptr || !line->at_boundary()))
                      // one would expect
                      // make_pair(line->vertex(0),
                      //           line->vertex(1))
                      // here, but that is
                      // not dimension
                      // independent, since
                      // vertex(i) is
                      // Point<dim>, but we
                      // want a Point<2>.
                      // in fact, whenever
                      // we're here, the
                      // vertex is a
                      // Point<dim>, but
                      // the compiler does
                      // not know
                      // this. hopefully,
                      // the compiler will
                      // optimize away this
                      // little kludge
                      line_list.emplace_back(
                        Point<2>(line->vertex(0)[0], line->vertex(0)[1]),
                        Point<2>(line->vertex(1)[0], line->vertex(1)[1]),
                        line->user_flag_set(),
                        cell->level());
                  }

              // next if we are to treat
              // curved boundaries
              // specially, then add lines
              // to the list consisting of
              // pieces of the boundary
              // lines
              if (mapping != nullptr)
                {
                  // to do so, first
                  // generate a sequence of
                  // points on a face and
                  // project them onto the
                  // faces of a unit cell
                  std::vector<Point<dim - 1>> boundary_points(n_points);

                  for (unsigned int i = 0; i < n_points; ++i)
                    boundary_points[i][0] = 1. * (i + 1) / (n_points + 1);

                  const Quadrature<dim - 1> quadrature(boundary_points);
                  const Quadrature<dim>     q_projector(
                    QProjector<dim>::project_to_all_faces(
                      ReferenceCells::get_hypercube<dim>(), quadrature));

                  // next loop over all
                  // boundary faces and
                  // generate the info from
                  // them
                  for (const auto &cell : tria.active_cell_iterators())
                    for (const unsigned int face_no :
                         GeometryInfo<dim>::face_indices())
                      {
                        const typename dealii::Triangulation<dim, spacedim>::
                          face_iterator face = cell->face(face_no);

                        if (face->at_boundary())
                          {
                            Point<dim> p0_dim(face->vertex(0));
                            Point<2>   p0(p0_dim[0], p0_dim[1]);

                            // loop over
                            // all pieces
                            // of the line
                            // and generate
                            // line-lets
                            const auto offset =
                              QProjector<dim>::DataSetDescriptor::face(
                                cell->reference_cell(),
                                face_no,
                                cell->combined_face_orientation(face_no),
                                n_points);
                            for (unsigned int i = 0; i < n_points; ++i)
                              {
                                const Point<dim> p1_dim(
                                  mapping->transform_unit_to_real_cell(
                                    cell, q_projector.point(offset + i)));
                                const Point<2> p1(p1_dim[0], p1_dim[1]);

                                line_list.emplace_back(p0,
                                                       p1,
                                                       face->user_flag_set(),
                                                       cell->level());
                                p0 = p1;
                              }

                            // generate last piece
                            const Point<dim> p1_dim(face->vertex(1));
                            const Point<2>   p1(p1_dim[0], p1_dim[1]);
                            line_list.emplace_back(p0,
                                                   p1,
                                                   face->user_flag_set(),
                                                   cell->level());
                          }
                      }
                }

              break;
            }

          case 3:
            {
              // curved boundary output
              // presently not supported
              Assert(mapping == nullptr, ExcNotImplemented());

              // loop over all lines and compute their
              // projection on the plane perpendicular
              // to the direction of sight

              // direction of view equals the unit
              // vector of the position of the
              // spectator to the origin.
              //
              // we chose here the viewpoint as in
              // gnuplot as default.
              //
              // TODO:[WB] Fix a potential problem with viewing angles in 3d Eps
              // GridOut
              // note: the following might be wrong
              // if one of the base vectors below
              // is in direction of the viewer, but
              // I am too tired at present to fix
              // this
              const double     pi         = numbers::PI;
              const double     z_angle    = eps_flags_3.azimut_angle;
              const double     turn_angle = eps_flags_3.turn_angle;
              const Point<dim> view_direction(
                -std::sin(z_angle * 2. * pi / 360.) *
                  std::sin(turn_angle * 2. * pi / 360.),
                +std::sin(z_angle * 2. * pi / 360.) *
                  std::cos(turn_angle * 2. * pi / 360.),
                -std::cos(z_angle * 2. * pi / 360.));

              // decide about the two unit vectors
              // in this plane. we chose the first one
              // to be the projection of the z-axis
              // to this plane
              const Tensor<1, dim> vector1 =
                Point<dim>(0, 0, 1) -
                ((Point<dim>(0, 0, 1) * view_direction) * view_direction);
              const Tensor<1, dim> unit_vector1 = vector1 / vector1.norm();

              // now the third vector is fixed. we
              // chose the projection of a more or
              // less arbitrary vector to the plane
              // perpendicular to the first one
              const Tensor<1, dim> vector2 =
                (Point<dim>(1, 0, 0) -
                 ((Point<dim>(1, 0, 0) * view_direction) * view_direction) -
                 ((Point<dim>(1, 0, 0) * unit_vector1) * unit_vector1));
              const Tensor<1, dim> unit_vector2 = vector2 / vector2.norm();


              for (const auto &cell : tria.active_cell_iterators())
                for (const unsigned int line_no : cell->line_indices())
                  {
                    typename dealii::Triangulation<dim, spacedim>::line_iterator
                      line = cell->line(line_no);
                    line_list.emplace_back(
                      Point<2>(line->vertex(0) * unit_vector2,
                               line->vertex(0) * unit_vector1),
                      Point<2>(line->vertex(1) * unit_vector2,
                               line->vertex(1) * unit_vector1),
                      line->user_flag_set(),
                      cell->level());
                  }

              break;
            }

          default:
            DEAL_II_NOT_IMPLEMENTED();
        }



      // find out minimum and maximum x and
      // y coordinates to compute offsets
      // and scaling factors
      double       x_min     = tria.begin_active()->vertex(0)[0];
      double       x_max     = x_min;
      double       y_min     = tria.begin_active()->vertex(0)[1];
      double       y_max     = y_min;
      unsigned int max_level = line_list.begin()->level;

      for (LineList::const_iterator line = line_list.begin();
           line != line_list.end();
           ++line)
        {
          x_min = std::min(x_min, line->first[0]);
          x_min = std::min(x_min, line->second[0]);

          x_max = std::max(x_max, line->first[0]);
          x_max = std::max(x_max, line->second[0]);

          y_min = std::min(y_min, line->first[1]);
          y_min = std::min(y_min, line->second[1]);

          y_max = std::max(y_max, line->first[1]);
          y_max = std::max(y_max, line->second[1]);

          max_level = std::max(max_level, line->level);
        }

      // scale in x-direction such that
      // in the output 0 <= x <= 300.
      // don't scale in y-direction to
      // preserve the shape of the
      // triangulation
      const double scale =
        (eps_flags_base.size /
         (eps_flags_base.size_type == GridOutFlags::EpsFlagsBase::width ?
            x_max - x_min :
            y_min - y_max));


      // now write preamble
      {
        // block this to have local
        // variables destroyed after
        // use
        std::time_t time1 = std::time(nullptr);
        std::tm    *time  = std::localtime(&time1);
        out << "%!PS-Adobe-2.0 EPSF-1.2" << '\n'
            << "%%Title: deal.II Output" << '\n'
            << "%%Creator: the deal.II library" << '\n'
            << "%%Creation Date: " << time->tm_year + 1900 << "/"
            << time->tm_mon + 1 << "/" << time->tm_mday << " - "
            << time->tm_hour << ":" << std::setw(2) << time->tm_min << ":"
            << std::setw(2) << time->tm_sec << '\n'
            << "%%BoundingBox: "
            // lower left corner
            << "0 0 "
            // upper right corner
            << static_cast<unsigned int>(
                 std::floor(((x_max - x_min) * scale) + 1))
            << ' '
            << static_cast<unsigned int>(
                 std::floor(((y_max - y_min) * scale) + 1))
            << '\n';

        // define some abbreviations to keep
        // the output small:
        // m=move turtle to
        // x=execute line stroke
        // b=black pen
        // r=red pen
        out << "/m {moveto} bind def" << '\n'
            << "/x {lineto stroke} bind def" << '\n'
            << "/b {0 0 0 setrgbcolor} def" << '\n'
            << "/r {1 0 0 setrgbcolor} def" << '\n';

        // calculate colors for level
        // coloring; level 0 is black,
        // other levels are blue
        // ... red
        if (eps_flags_base.color_lines_level)
          out << "/l  { neg " << (max_level) << " add "
              << (0.66666 / std::max(1U, (max_level - 1)))
              << " mul 1 0.8 sethsbcolor} def" << '\n';

        // in 2d, we can also plot cell
        // and vertex numbers, but this
        // requires a somewhat more
        // lengthy preamble. please
        // don't ask me what most of
        // this means, it is reverse
        // engineered from what GNUPLOT
        // uses in its output
        if ((dim == 2) && (eps_flags_2.write_cell_numbers ||
                           eps_flags_2.write_vertex_numbers))
          {
            out
              << ("/R {rmoveto} bind def\n"
                  "/Symbol-Oblique /Symbol findfont [1 0 .167 1 0 0] makefont\n"
                  "dup length dict begin {1 index /FID eq {pop pop} {def} ifelse} forall\n"
                  "currentdict end definefont\n"
                  "/MFshow {{dup dup 0 get findfont exch 1 get scalefont setfont\n"
                  "[ currentpoint ] exch dup 2 get 0 exch rmoveto dup dup 5 get exch 4 get\n"
                  "{show} {stringwidth pop 0 rmoveto}ifelse dup 3 get\n"
                  "{2 get neg 0 exch rmoveto pop} {pop aload pop moveto}ifelse} forall} bind def\n"
                  "/MFwidth {0 exch {dup 3 get{dup dup 0 get findfont exch 1 get scalefont setfont\n"
                  "5 get stringwidth pop add}\n"
                  "{pop} ifelse} forall} bind def\n"
                  "/MCshow { currentpoint stroke m\n"
                  "exch dup MFwidth -2 div 3 -1 roll R MFshow } def\n")
              << '\n';
          }

        out << "%%EndProlog" << '\n' << '\n';

        // set fine lines
        out << eps_flags_base.line_width << " setlinewidth" << '\n';
      }

      // now write the lines
      const Point<2> offset(x_min, y_min);

      for (LineList::const_iterator line = line_list.begin();
           line != line_list.end();
           ++line)
        if (eps_flags_base.color_lines_level && (line->level > 0))
          out << line->level << " l " << (line->first - offset) * scale << " m "
              << (line->second - offset) * scale << " x" << '\n';
        else
          out << ((line->colorize && eps_flags_base.color_lines_on_user_flag) ?
                    "r " :
                    "b ")
              << (line->first - offset) * scale << " m "
              << (line->second - offset) * scale << " x" << '\n';

      // finally write the cell numbers
      // in 2d, if that is desired
      if ((dim == 2) && (eps_flags_2.write_cell_numbers == true))
        {
          out << "(Helvetica) findfont 140 scalefont setfont" << '\n';

          for (const auto &cell : tria.active_cell_iterators())
            {
              out << (cell->center()[0] - offset[0]) * scale << ' '
                  << (cell->center()[1] - offset[1]) * scale << " m" << '\n'
                  << "[ [(Helvetica) 12.0 0.0 true true (";
              if (eps_flags_2.write_cell_number_level)
                out << cell;
              else
                out << cell->index();

              out << ")] "
                  << "] -6 MCshow" << '\n';
            }
        }

      // and the vertex numbers
      if ((dim == 2) && (eps_flags_2.write_vertex_numbers == true))
        {
          out << "(Helvetica) findfont 140 scalefont setfont" << '\n';

          // have a list of those
          // vertices which we have
          // already tracked, to avoid
          // doing this multiply
          std::set<unsigned int> treated_vertices;
          for (const auto &cell : tria.active_cell_iterators())
            for (const unsigned int vertex_no : cell->vertex_indices())
              if (treated_vertices.find(cell->vertex_index(vertex_no)) ==
                  treated_vertices.end())
                {
                  treated_vertices.insert(cell->vertex_index(vertex_no));

                  out << (cell->vertex(vertex_no)[0] - offset[0]) * scale << ' '
                      << (cell->vertex(vertex_no)[1] - offset[1]) * scale
                      << " m" << '\n'
                      << "[ [(Helvetica) 10.0 0.0 true true ("
                      << cell->vertex_index(vertex_no) << ")] "
                      << "] -6 MCshow" << '\n';
                }
        }

      out << "showpage" << '\n';

      // make sure everything now gets to
      // disk
      out.flush();

      AssertThrow(out.fail() == false, ExcIO());
    }
  } // namespace
} // namespace internal


template <int dim, int spacedim>
void
GridOut::write_eps(const Triangulation<dim, spacedim> &tria,
                   std::ostream                       &out,
                   const Mapping<dim, spacedim>       *mapping) const
{
  internal::write_eps(tria, out, mapping, eps_flags_2, eps_flags_3);
}


template <int dim, int spacedim>
void
GridOut::write(const Triangulation<dim, spacedim> &tria,
               std::ostream                       &out,
               const OutputFormat                  output_format,
               const Mapping<dim, spacedim>       *mapping) const
{
  switch (output_format)
    {
      case none:
        return;

      case dx:
        write_dx(tria, out);
        return;

      case ucd:
        write_ucd(tria, out);
        return;

      case gnuplot:
        write_gnuplot(tria, out, mapping);
        return;

      case eps:
        write_eps(tria, out, mapping);
        return;

      case xfig:
        write_xfig(tria, out, mapping);
        return;

      case msh:
        write_msh(tria, out);
        return;

      case svg:
        write_svg(tria, out);
        return;

      case mathgl:
        write_mathgl(tria, out);
        return;

      case vtk:
        write_vtk(tria, out);
        return;

      case vtu:
        write_vtu(tria, out);
        return;
    }

  DEAL_II_ASSERT_UNREACHABLE();
}


template <int dim, int spacedim>
void
GridOut::write(const Triangulation<dim, spacedim> &tria,
               std::ostream                       &out,
               const Mapping<dim, spacedim>       *mapping) const
{
  write(tria, out, default_format, mapping);
}


// explicit instantiations
#include "grid/grid_out.inst"


DEAL_II_NAMESPACE_CLOSE
