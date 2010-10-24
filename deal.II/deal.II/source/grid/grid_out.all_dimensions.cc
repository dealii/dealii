//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2005, 2006, 2008, 2010 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------


#include <grid/grid_out.h>
#include <base/parameter_handler.h>

DEAL_II_NAMESPACE_OPEN


namespace GridOutFlags
{
  DX::DX (const bool write_cells,
	  const bool write_faces,
	  const bool write_diameter,
	  const bool write_measure,
	  const bool write_all_faces) :
    write_cells (write_cells),
    write_faces (write_faces),
    write_diameter (write_diameter),
    write_measure (write_measure),
    write_all_faces (write_all_faces)
  {}

  void DX::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Write cells", "true", Patterns::Bool(),
			"Write the mesh connectivity as DX grid cells");
    param.declare_entry("Write faces", "false", Patterns::Bool(),
			"Write faces of cells. These may be boundary faces "
			"or all faces between mesh cells, according to "
			"\"Write all faces\"");
    param.declare_entry("Write diameter", "false", Patterns::Bool(),
			"If cells are written, additionally write their"
			" diameter as data for visualization");
    param.declare_entry("Write measure", "false", Patterns::Bool(),
			"Write the volume of each cell as data");
    param.declare_entry("Write all faces", "true", Patterns::Bool(),
			"Write all faces, not only boundary");
  }

  void DX::parse_parameters (ParameterHandler& param)
  {
    write_cells = param.get_bool("Write cells");
    write_faces = param.get_bool("Write faces");
    write_diameter = param.get_bool("Write diameter");
    write_measure = param.get_bool("Write measure");
    write_all_faces = param.get_bool("Write all faces");
  }


  Msh::Msh (const bool write_faces,
	    const bool write_lines) :
    write_faces (write_faces),
    write_lines (write_lines)
  {}

  void Msh::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Write faces", "false", Patterns::Bool());
    param.declare_entry("Write lines", "false", Patterns::Bool());
  }


  void Msh::parse_parameters (ParameterHandler& param)
  {
    write_faces = param.get_bool("Write faces");
    write_lines = param.get_bool("Write lines");
  }


  Ucd::Ucd (const bool write_preamble,
	    const bool write_faces,
	    const bool write_lines) :
		  write_preamble (write_preamble),
		  write_faces (write_faces),
		  write_lines (write_lines)
  {}



  void Ucd::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Write preamble", "true", Patterns::Bool());
    param.declare_entry("Write faces", "false", Patterns::Bool());
    param.declare_entry("Write lines", "false", Patterns::Bool());
  }


  void Ucd::parse_parameters (ParameterHandler& param)
  {
    write_preamble = param.get_bool("Write preamble");
    write_faces = param.get_bool("Write faces");
    write_lines = param.get_bool("Write lines");
  }


  Gnuplot::Gnuplot (const bool write_cell_numbers,
		    const unsigned int n_boundary_face_points,
		    const bool         curved_inner_cells) :
		  write_cell_numbers (write_cell_numbers),
		  n_boundary_face_points(n_boundary_face_points),
		  curved_inner_cells(curved_inner_cells)
  {}



  void Gnuplot::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Cell number", "false", Patterns::Bool());
    param.declare_entry("Boundary points", "2", Patterns::Integer());
  }


  void Gnuplot::parse_parameters (ParameterHandler& param)
  {
    write_cell_numbers = param.get_bool("Cell number");
    n_boundary_face_points = param.get_integer("Boundary points");
  }


  EpsFlagsBase::EpsFlagsBase (const SizeType     size_type,
			      const unsigned int size,
			      const double       line_width,
			      const bool color_lines_on_user_flag,
			      const unsigned int n_boundary_face_points,
			      const bool color_lines_level) :
		  size_type (size_type),
		  size (size),
		  line_width (line_width),
		  color_lines_on_user_flag(color_lines_on_user_flag),
		  n_boundary_face_points(n_boundary_face_points),
		  color_lines_level(color_lines_level)
  {}


  void EpsFlagsBase::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Size by", "width",
			    Patterns::Selection("width|height"),
			"Depending on this parameter, either the"
			"width or height "
			"of the eps is scaled to \"Size\"");
    param.declare_entry("Size", "300", Patterns::Integer(),
			"Size of the output in points");
    param.declare_entry("Line width", "0.5", Patterns::Double(),
			"Width of the lines drawn in points");
    param.declare_entry("Color by flag", "false", Patterns::Bool(),
			"Draw lines with user flag set in different color");
    param.declare_entry("Boundary points", "2", Patterns::Integer(),
			"Number of points on boundary edges. "
			"Increase this beyond 2 to see curved boundaries.");
    param.declare_entry("Color by level", "false", Patterns::Bool(),
			"Draw different colors according to grid level.");
  }


  void EpsFlagsBase::parse_parameters (ParameterHandler& param)
  {
    if (param.get("Size by") == std::string("width"))
      size_type = width;
    else if (param.get("Size by") == std::string("height"))
      size_type = height;
    size = param.get_integer("Size");
    line_width = param.get_double("Line width");
    color_lines_on_user_flag = param.get_bool("Color by flag");
    n_boundary_face_points = param.get_integer("Boundary points");
    color_lines_level = param.get_bool("Color by level");
  }



  Eps<1>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points)
  {}


  void Eps<1>::declare_parameters (ParameterHandler&)
  {}


  void Eps<1>::parse_parameters (ParameterHandler& param)
  {
    EpsFlagsBase::parse_parameters(param);
  }



  Eps<2>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points,
	       const bool         write_cell_numbers,
	       const bool         write_cell_number_level,
	       const bool         write_vertex_numbers,
	       const bool         color_lines_level
      )
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points,
			       color_lines_level),
		  write_cell_numbers (write_cell_numbers),
		  write_cell_number_level (write_cell_number_level),
		  write_vertex_numbers (write_vertex_numbers)
  {}


  void Eps<2>::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Cell number", "false", Patterns::Bool(),
			"(2D only) Write cell numbers"
			" into the centers of cells");
    param.declare_entry("Level number", "false", Patterns::Bool(),
			"(2D only) if \"Cell number\" is true, write"
			"numbers in the form level.number");
    param.declare_entry("Vertex number", "false", Patterns::Bool(),
			"Write numbers for each vertex");
  }


  void Eps<2>::parse_parameters (ParameterHandler& param)
  {
    EpsFlagsBase::parse_parameters(param);
    write_cell_numbers = param.get_bool("Cell number");
    write_cell_number_level = param.get_bool("Level number");
    write_vertex_numbers = param.get_bool("Vertex number");
  }



  Eps<3>::Eps (const SizeType     size_type,
	       const unsigned int size,
	       const double       line_width,
	       const bool color_lines_on_user_flag,
	       const unsigned int n_boundary_face_points,
	       const double        azimut_angle,
	       const double        turn_angle)
		  :
		  EpsFlagsBase(size_type, size, line_width,
			       color_lines_on_user_flag,
			       n_boundary_face_points),
		  azimut_angle (azimut_angle),
		  turn_angle (turn_angle)
  {}


  void Eps<3>::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Azimuth", "30", Patterns::Double(),
			"Azimuth of the viw point, that is, the angle "
			"in the plane from the x-axis.");
    param.declare_entry("Elevation", "30", Patterns::Double(),
			"Elevation of the view point above the xy-plane.");
  }


  void Eps<3>::parse_parameters (ParameterHandler& param)
  {
    EpsFlagsBase::parse_parameters(param);
    azimut_angle = 90- param.get_double("Elevation");
    turn_angle = param.get_double("Azimuth");
  }



  XFig::XFig ()
		  :
    draw_boundary(true),
    level_color(false),
    level_depth(true),
    n_boundary_face_points(0),
    scaling(1.,1.),
    fill_style (20),
    line_style(0),
    line_thickness(1),
    boundary_style(0),
    boundary_thickness(3)
  {}


  void XFig::declare_parameters (ParameterHandler& param)
  {
    param.declare_entry("Boundary", "true", Patterns::Bool());
    param.declare_entry("Level color", "false", Patterns::Bool());
    param.declare_entry("Level depth", "true", Patterns::Bool());
//TODO: Unify this number with other output formats
    param.declare_entry("Boundary points", "0", Patterns::Integer());
    param.declare_entry("Fill style", "20", Patterns::Integer());
    param.declare_entry("Line style", "0", Patterns::Integer());
    param.declare_entry("Line width", "1", Patterns::Integer());
    param.declare_entry("Boundary style", "0", Patterns::Integer());
    param.declare_entry("Boundary width", "3", Patterns::Integer());
  }


  void XFig::parse_parameters (ParameterHandler& param)
  {
    draw_boundary = param.get_bool("Boundary");
    level_color = param.get_bool("Level color");
    level_depth = param.get_bool("Level depth");
    n_boundary_face_points = param.get_integer("Boundary points");
    fill_style = param.get_integer("Fill style");
    line_style = param.get_integer("Line style");
    line_thickness = param.get_integer("Line width");
    boundary_style = param.get_integer("Boundary style");
    boundary_thickness = param.get_integer("Boundary width");
  }


}  // end namespace GridOutFlags



GridOut::GridOut ()
		:
		default_format (none)
{}


void GridOut::set_flags (const GridOutFlags::DX &flags)
{
  dx_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Msh &flags)
{
  msh_flags = flags;
}


void GridOut::set_flags (const GridOutFlags::Ucd &flags)
{
  ucd_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Gnuplot &flags)
{
  gnuplot_flags = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<1> &flags)
{
  eps_flags_1 = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<2> &flags)
{
  eps_flags_2 = flags;
}



void GridOut::set_flags (const GridOutFlags::Eps<3> &flags)
{
  eps_flags_3 = flags;
}



void GridOut::set_flags (const GridOutFlags::XFig &flags)
{
  xfig_flags = flags;
}



std::string
GridOut::default_suffix (const OutputFormat output_format)
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
      default:
	    Assert (false, ExcNotImplemented());
	    return "";
    };
}



std::string
GridOut::default_suffix () const
{
  return default_suffix(default_format);
}



GridOut::OutputFormat
GridOut::parse_output_format (const std::string &format_name)
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

  AssertThrow (false, ExcInvalidState ());
				   // return something weird
  return OutputFormat(-1);
}



std::string GridOut::get_output_format_names ()
{
  return "none|dx|gnuplot|eps|ucd|xfig|msh";
}


void
GridOut::declare_parameters(ParameterHandler& param)
{
  param.declare_entry("Format", "none",
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
}



void
GridOut::parse_parameters(ParameterHandler& param)
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
}



unsigned int
GridOut::memory_consumption () const
{
  return (sizeof(dx_flags) +
	  sizeof(msh_flags) +
	  sizeof(ucd_flags) +
	  sizeof(gnuplot_flags) +
	  sizeof(eps_flags_1) +
	  sizeof(eps_flags_2) +
	  sizeof(eps_flags_3) +
	  sizeof(xfig_flags));
}

DEAL_II_NAMESPACE_CLOSE
