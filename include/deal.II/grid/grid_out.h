// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2013 by the deal.II authors
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

#ifndef __deal2__grid_out_h
#define __deal2__grid_out_h



#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/point.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/mapping.h>

#include <string>

DEAL_II_NAMESPACE_OPEN

class ParameterHandler;
template <int dim, int spacedim> class Triangulation;
template <int dim, int spacedim> class Mapping;


/**
 * Within this namespace, we define several structures that are used
 * to describe flags that can be given to grid output routines to
 * modify the default outfit of the grids written into a file. See the
 * different subclasses and the documentation of the GridOut
 * class for more details.
 *
 * @ingroup output
 */
namespace GridOutFlags
{
  /**
   * Flags for grid output in OpenDX format.
   *
   * @ingroup output
   */
  struct DX
  {
    /**
     * Write cells.
     */
    bool write_cells;

    /**
     * Write faces.
     */
    bool write_faces;

    /**
     * Write field with diameters.
     */
    bool write_diameter;

    /**
     * Write field with area/volume.
     */
    bool write_measure;

    /**
     * Write all faces, including interior faces. If <tt>false</tt>,
     * only boundary faces are written.
     */
    bool write_all_faces;

    /**
     * Constructor.
     */
    DX (const bool write_cells = true,
        const bool write_faces = false,
        const bool write_diameter = false,
        const bool write_measure = false,
        const bool write_all_faces = true);

    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };

  /**
   * Flags describing the details of output in MSH format.
   *
   * @ingroup output
   */
  struct Msh
  {
    /**
     * When writing a mesh, write boundary faces explicitly if their
     * boundary indicator is not the default boundary indicator, which
     * is zero.  This is necessary if you later want to re-read the
     * grid and want to get the same boundary indicators for the
     * different parts of the boundary of the triangulation.
     *
     * It is not necessary if you only want to write the triangulation
     * to view or print it.
     *
     * Default: @p false.
     */
    bool write_faces;
    /**
      * When writing a mesh, write boundary lines explicitly if their
      * boundary indicator is not the default boundary indicator,
      * which is zero.  This is necessary if you later want to re-read
      * the grid and want to get the same boundary indicators for the
      * different parts of the boundary of the triangulation.
      *
      * It is not necessary if you only want to write the
      * triangulation to view or print it.
      *
      * This is used only if <tt>dim==3</tt>, and ignored in all other
      * cases.
      *
      * Default: @p false.
      */
    bool write_lines;

    /**
     * Constructor.
     */
    Msh (const bool write_faces    = false,
         const bool write_lines    = false);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags describing the details of output in UCD format.
   *
   * @ingroup output
   */
  struct Ucd
  {
    /**
     * Write a comment at the beginning of the file stating the date
     * of creation and some other data.  While this is supported by
     * the UCD format (and the AVS program), some other programs get
     * confused by this, so the default is to not write a
     * preamble. However, a preamble can be written using this flag.
     *
     * Default: <code>false</code>.
     */
    bool write_preamble;

    /**
     * When writing a mesh, write boundary faces explicitly if their
     * boundary indicator is not the default boundary indicator, which
     * is zero.  This is necessary if you later want to re-read the
     * grid and want to get the same boundary indicators for the
     * different parts of the boundary of the triangulation.
     *
     * It is not necessary if you only want to write the triangulation
     * to view or print it.
     *
     * Default: @p false.
     */
    bool write_faces;

    /**
     * When writing a mesh, write boundary lines explicitly if their
     * boundary indicator is not the default boundary indicator, which
     * is zero.  This is necessary if you later want to re-read the
     * grid and want to get the same boundary indicators for the
     * different parts of the boundary of the triangulation.
     *
     * It is not necessary if you only want to write the triangulation
     * to view or print it.
     *
     * This directive is ignored if
     * <tt>dim!=3</tt>.
     *
     *  Default: @p false.
     */
    bool write_lines;

    /**
     * Constructor.
     */
    Ucd (const bool write_preamble = false,
         const bool write_faces    = false,
         const bool write_lines    = false);

    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags describing the details of output in GNUPLOT format.
   *
   * @ingroup output
   */
  struct Gnuplot
  {
    /**
     * Write the number of each cell into the output file before
     * starting with the lines it is composed of, as a comment. This
     * might be useful if you want to find out details about the grid,
     * for example the position of cells of which you know the
     * number. It enlarges the size of the output significantly,
     * however.
     *
     * Default: @p false.
     */
    bool write_cell_numbers;

    /**
     * Based on the vertices of the face and #n_boundary_face_points
     * additional points a tensor product mesh (tranformed to the real
     * space) of (#n_boundary_face_points+2)<sup>dim-1</sup> points is
     * plotted on each boundary face.
     */
    unsigned int n_boundary_face_points;

    /**
     * Flag. If true also inner cells are plotted with curved
     * boundaries. This is useful when for e.g.  MappingQEulerian with
     * #n_boundary_face_points&gt;.
     */
    bool curved_inner_cells;

    /**
     * Constructor.
     */
    Gnuplot (const bool         write_cell_number = false,
             const unsigned int n_boundary_face_points = 2,
             const bool         curved_inner_cells = false);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };

  /**
   * Flags describing the details of output for encapsulated
   * postscript.  In this structure, the flags common to all
   * dimensions are listed. Flags which are specific to one space
   * dimension only are listed in derived classes.
   *
   * By default, the size of the picture is scaled such that the width
   * equals 300 units.
   *
   * @ingroup output
   */
  struct EpsFlagsBase
  {
    /**
     * Enum denoting the possibilities whether the scaling should be
     * done such that the given @p size equals the width or the height
     * of the resulting picture.
     */
    enum SizeType
    {
      width, height
    };

    /**
     * See above. Default is @p width.
     */
    SizeType size_type;

    /**
     * Width or height of the output as given in postscript units This
     * usually is given by the strange unit 1/72 inch. Whether this is
     * height or width is specified by the flag @p size_type.
     *
     * Default is 300.
     */
    unsigned int size;

    /**
     * Width of a line in postscript units. Default is 0.5.
     */
    double line_width;

    /**
     * Should lines with a set @p user_flag be drawn in a different
     * color (red)?  See @ref GlossUserFlags for information about
     * user flags.
    */
    bool color_lines_on_user_flag;

    /**
     * This is the number of points on a boundary face, that are
     * ploted additionally to the vertices of the face.
     *
     * This is used if the mapping used is not the standard @p
     * MappingQ1 mapping.
     */
    unsigned int n_boundary_face_points;

    /**
     * Should lines be colored according to their refinement level?
     * This overrides color_lines_on_user_flag for all levels except
     * level 0.  Colors are: level 0: black, other levels: rainbow
     * scale from blue to red.
     */
    bool color_lines_level;

    /**
     * Constructor.
     */
    EpsFlagsBase (const SizeType     size_type  = width,
                  const unsigned int size       = 300,
                  const double       line_width = 0.5,
                  const bool color_lines_on_user_flag = false,
                  const unsigned int n_boundary_face_points = 2,
                  const bool color_lines_level = false);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of
     * ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags describing the details of output for encapsulated
   * postscript for all dimensions not explicitly specialized
   * below. Some flags that are common to all dimensions are listed in
   * the base class.
   *
   * This class does not actually exist, we only here declare the
   * general template and declare explicit specializations below.
   *
   * @ingroup output
  */
  template <int dim>
  struct Eps
  {};

  /**
   * Flags specific to the output of grids in one space dimensions.
   *
   * @ingroup output
   */
  template <>
  struct Eps<1> : public EpsFlagsBase
  {
    /**
     * Constructor.
     */
    Eps (const SizeType     size_type  = width,
         const unsigned int size       = 300,
         const double       line_width = 0.5,
         const bool         color_lines_on_user_flag = false,
         const unsigned int n_boundary_face_points = 2);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags specific to the output of grids in two space dimensions.
   *
   * @ingroup output
   */
  template <>
  struct Eps<2> : public EpsFlagsBase
  {
    /**
     * If this flag is set, then we place the number of the cell into
     * the middle of each cell. The default value is to not do this.
     *
     * The format of the cell number written is <tt>level.index</tt>,
     * or simply @p index, depending on the value of the following
     * flag.
     */
    bool write_cell_numbers;
    /**
     * If the cell numbers shall be written, using the above flag,
     * then the value of this flag determines whether the format shall
     * be <tt>level.index</tt>, or simply @p index. If @p true, the
     * first format is taken. Default is @p true.
     *
     * The flag has obviously no effect if @p write_cell_numbers is @p
     * false.
     */
    bool write_cell_number_level;

    /**
     * Vertex numbers can be written onto the vertices. This is
     * controlled by the following flag. Default is @p false.
     */
    bool write_vertex_numbers;

    /**
     * Constructor.
     */
    Eps (const SizeType     size_type  = width,
         const unsigned int size       = 300,
         const double       line_width = 0.5,
         const bool         color_lines_on_user_flag = false,
         const unsigned int n_boundary_face_points = 2,
         const bool         write_cell_numbers = false,
         const bool         write_cell_number_level = true,
         const bool         write_vertex_numbers = false,
         const bool         color_lines_level = false);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };

  /**
   * Flags specific to the output of grids in three space dimensions.
   *
   * @ingroup output
   */
  template <>
  struct Eps<3> : public EpsFlagsBase
  {
    /**
     * Angle of the line origin-viewer against the z-axis in degrees.
     *
     * Default is the Gnuplot-default of 60.
     */
    double azimut_angle;

    /**
     * Angle by which the viewers position projected onto the
     * x-y-plane is rotated around the z-axis, in positive sense when
     * viewed from above. The unit are degrees, and zero equals a
     * position above or below the negative y-axis.
     *
     * Default is the Gnuplot-default of 30.
     */
    double turn_angle;

    /**
     * Constructor.
     */
    Eps (const SizeType     size_type  = width,
         const unsigned int size       = 300,
         const double       line_width = 0.5,
         const bool         color_lines_on_user_flag = false,
         const unsigned int n_boundary_face_points = 2,
         const double       azimut_angle    = 60,
         const double       turn_angle      = 30);
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };

  /**
   * Flags for XFig output.
   *
   * @ingroup output
   */
  struct XFig
  {
    /**
     * Draw boundary lines. Default is true.
     */
    bool draw_boundary;

    /**
     * An enum used for deciding which field is used for coloring the cells.
     */
    enum Coloring
    {
      /// Convert the material id into the cell color
      material_id,
      /// Convert the level into the cell color
      level_number,
      /// Convert the global subdomain id into the cell color
      subdomain_id,
      /// Convert the level subdomain id into the cell color
      level_subdomain_id
    } color_by;

    /**
     * @deprecated Use the color_by enum instead.
     *
     * Change color depending on level. Default is false, therefore,
     * color is coded by material or boundary id.
     */
    bool level_color DEAL_II_DEPRECATED;

    /**
     * Code level to depth. Default is true. If false, color depends
     * on material or boundary id.
     *
     * Depth of the object is 900-level, if this value is true.
     */
    bool level_depth;

    /**
     * Additional points for curved boundaries. Default is none.
     */
    unsigned int n_boundary_face_points;

    /**
     * Scaling of graph. The default is a unit length of one inch.
     */
    Point<2> scaling;

    /**
     * Offset of the graph. Before scaling, the coordinates are
     * shifted by this value. Default is zero in each direction.
     */
    Point<2> offset;

    /**
     * Style for filling cells. Default is solid fill (20). This value
     * is forwarded unchanged into the corresponding field
     * <tt>fill_style</tt> of the polyline object of XFig.
     */
    int fill_style;

    /**
     * Style for drawing border lines of polygons. Defaults to solid
     * (0) and is forwarded to XFig.
     */
    int line_style;

    /**
     * Thickness of border lines of polygons. Default is 1.
     *
     * Set this to zero to avoid border lines for very fine meshes.
     */
    int line_thickness;

    /**
     * Style for drawing lines at the boundary. Defaults to solid (0).
     */
    int boundary_style;

    /**
     * Thickness of boundary lines. Default is 3.
     */
    int boundary_thickness;

    /**
     * Constructor.
     */
    XFig();
    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags controlling SVG output.
   *
   * @ingroup output
   */
  struct Svg
  {
    /// Height of the plot in SVG units, computed from width if zero. Defaults to 1000
    unsigned int height;
    /// The width of the plot. Computed automatically from height if zero (default)
    unsigned int width;
    /// Thickness of the lines between cells
    unsigned int line_thickness;
    /// Thickness of lines at the boundary
    unsigned int boundary_line_thickness;

    /// Margin around the plotted area
    bool margin;

    /**
     * Background style.
     */
    enum Background
    {
      /// Use transparent value of SVG
      transparent,
      /// Use white background
      white,
      /// Use a gradient from white (top) to steelblue (bottom), and add date and time plus a deal.II logo. Automatically draws a margin.
      dealii
    };

    Background background;

    // View angles for the perspective view of the grid; Default is 0, 0 (top view).
    /**
     * The azimuth angle measured from ??? in degrees. Default is 0.
     */
    int azimuth_angle;
    /**
     * The angle from vertically above the xy-plane. Default is 0.
     */
    int polar_angle;

    /**
     * Cell coloring.
     */
    enum Coloring
    {
      /// No cell coloring
      none,
      /// Convert the material id into the cell color (default)
      material_id,
      /// Convert the level number into the cell color
      level_number,
      /// Convert the subdomain id into the cell color
      subdomain_id,
      /// Convert the level subdomain id into the cell color
      level_subdomain_id
    };

    Coloring coloring;

    /// Interpret the level number of the cells as altitude over the x-y-plane (useful in the perpspective view).
    bool convert_level_number_to_height;

    /// The factor determining the vertical distance between levels (default = 0.3)
    float level_height_factor;

    /**
     * Cell labeling (fixed order).
     *
     * The following booleans determine which properties of the cell
     * shall be displayed as text in the middle of each cell.
     */
    bool label_level_number;    // default: true
    bool label_cell_index;      // default: true
    bool label_material_id;     // default: false
    bool label_subdomain_id;    // default: false
    bool label_level_subdomain_id;  // default: false

    /// Draw a colorbar next to the plotted grid with respect to the chosen coloring of the cells
    bool draw_colorbar;
    /// Draw a legend next to the plotted grid, explaining the label of the cells
    bool draw_legend;

    /**
     * Constructor.
     */
    Svg(const unsigned int line_thickness = 2,
        const unsigned int boundary_line_thickness = 4,
        bool margin = true,
        const Background background = white,
        const int azimuth_angle = 0,
        const int polar_angle = 0,
        const Coloring coloring = level_number,
        const bool convert_level_number_to_height = false,
        const bool label_level_number = true,
        const bool label_cell_index = true,
        const bool label_material_id = false,
        const bool label_subdomain_id = false,
        const bool draw_colorbar = true,
        const bool draw_legend = true);
  };

  /**
   * Flags for grid output in MathGL format.
   *
   * @ingroup output
   */
  struct MathGL
  {
    /**
     * Constructor.
     */
    MathGL ();

    /**
     * Draw a bounding box around the graph.
     */
    bool draw_bounding_box;

    /**
     * Declare parameters in ParameterHandler.
     */
    static void declare_parameters (ParameterHandler &param);

    /**
     * Parse parameters of ParameterHandler.
     */
    void parse_parameters (ParameterHandler &param);
  };


  /**
   * Flags for grid output in Vtk format. These flags are the same as
   * those declared in DataOutBase::VtkFlags.
   *
   * @ingroup output
   */
  struct Vtk : public DataOutBase::VtkFlags
  {};


  /**
   * Flags for grid output in Vtu format. These flags are the same as
   * those declared in DataOutBase::VtuFlags.
   *
   * @ingroup output
   */
  struct Vtu : public DataOutBase::VtkFlags
  {};
}



/**
 * This class provides a means to output a triangulation to a file in
 * different formats. See the enum GridOut::OutputFormat for a list of formats
 * and the corresponding output function names.
 *
 * Usage is simple: either you use the direct form
 * @code
 *   ofstream output_file("some_filename");
 *   GridOut().write_gnuplot (tria, output_file);
 * @endcode
 * if you know which format you want to have, or if you want the format to be
 * a runtime parameter, you can write
 * @code
 *   GridOut::OutputFormat grid_format =
 *                   GridOut::parse_output_format(get_format_name_from_somewhere());
 *   ofstream output_file("some_filename" + GridOut::default_suffix(output_format));
 *   GridOut().write (tria, output_file, output_format);
 * @endcode
 * The function <tt>get_output_format_names()</tt> provides a list of possible names of
 * output formats in a string that is understandable by the ParameterHandler class.
 *
 * Note that here, we have created an unnamed object of type GridOut and called
 * one of its <tt>write_*</tt> functions. This looks like as if the respective function
 * could really be made @p static. This was not done in order to allow for
 * parameters to be passed to the different output functions in a way compatible
 * with the scheme of allowing the right output format to be selected at run-time
 * through the generic @p write function.
 *
 * In order to explain this, consider each function had one or more additional
 * parameters giving the details of output, for example position of the spectator
 * for 3d meshed, line thicknesses, etc. While this would allow each output
 * function any flexibility it needs, it would not allow us to use the generic
 * function @p write which is given a parameter determining the output format,
 * since it is impractical to give it a list of parameters for each and every
 * output format supported which it may then pass on to the respective output
 * function.
 *
 * Rather, we have chosen to let each object of this class
 * GridOut have a set of parameters for each supported output
 * format. These are collected in structures GridOutFlags::Eps(),
 * GridOutFlags::Gnuplot(), etc declared in the GridOutFlags
 * namespace, and you can set your preferred flags like this:
 * @code
 *   GridOut grid_out;
 *   GridOutFlags::Ucd ucd_flags;
 *   ...    // set some fields in ucd_flags
 *   grid_out.set_flags (ucd_flags);
 *   ...
 *   ...    // write some file with data_out
 * @endcode
 * The respective output function then use the so-set flags. By default, they
 * are set to reasonable values as described above and in the documentation
 * of the different flags structures. Resetting the flags can
 * be done by calling <tt>grid_out.set_flags (GridOutFlags::Ucd());</tt>, since the
 * default constructor of each of the flags structures sets the parameters
 * to their initial values.
 *
 * The advantage of this approach is that it is possible to change the flags
 * of one or more output formats according to your needs and later use the
 * generic @p write function; the actual output function then called will
 * use the flags as set before.
 *
 * Note that some of the structures describing the flags of the different
 * output formats are empty since the respective format does not support
 * any flags. The structure and the @p set_flags function are provided
 * anyway. Note also that some of the structures may differ between the
 * dimensions supported by this class; they then have a template parameter,
 * as usual.
 *
 * @ingroup grid
 * @ingroup output
 * @author Wolfgang Bangerth, Guido Kanschat, Luca Heltai, Stefan Nauber, Christian Wülker
 * @date 1999 - 2013
 */
class GridOut
{
public:
  /**
   * Declaration of a name for each
   * of the different output
   * formats. These are used by the
   * generic output function
   * write() to determine the
   * actual output format.
   */
  enum OutputFormat
  {
    /// Do nothing in write()
    none,
    /// write() calls write_dx()
    dx,
    /// write() calls write_gnuplot()
    gnuplot,
    /// write() calls write_eps()
    eps,
    /// write() calls write_ucd()
    ucd,
    /// write() calls write_xfig()
    xfig,
    /// write() calls write_msh()
    msh,
    /// write() calls write_svg()
    svg,
    /// write() calls write_mathgl()
    mathgl,
    /// write() calls write_vtk()
    vtk,
    /// write() calls write_vtu()
    vtu
  };

  /**
   * Constructor.
   */
  GridOut ();

  /**
   * Write triangulation in OpenDX
   * format.
   *
   * Cells or faces are written
   * together with their level and
   * their material id or boundary
   * indicator, resp.
   *
   * Not implemented for the
   * codimension one case.
   */
  template <int dim, int spacedim>
  void write_dx (const Triangulation<dim,spacedim> &tria,
                 std::ostream             &out) const;

  /**
   * Write the triangulation in the
   * gnuplot format.
   *
   * In GNUPLOT format, each cell
   * is written as a sequence of
   * its confining lines. Apart
   * from the coordinates of the
   * line's end points, the level
   * and the material of the cell
   * are appended to each line of
   * output. Therefore, if you let
   * GNUPLOT draw a 2d grid as a 3d
   * plot, you will see more
   * refined cells being raised
   * against cells with less
   * refinement.  Also, if you draw
   * a cut through a 3d grid, you
   * can extrude the refinement
   * level in the direction
   * orthogonal to the cut
   * plane. The same can be done
   * with the material id, which is
   * plotted after the level.
   *
   * A more useful application of
   * this feature is the following:
   * if you use the GNUPLOT
   * command (for a 2d grid here)
   * @verbatim
   * splot [:][:][2.5:3.5] "grid_file.gnuplot" *
   * @endverbatim
   * then the
   * whole x- and y-range will be
   * plotted, i.e. the whole grid,
   * but only those lines with a
   * z-value between 2.5 and
   * 3.5. Since the z-values were
   * chosen to be the level to
   * which a cell belongs, this
   * results in a plot of those
   * cells only that belong to
   * level 3 in this example. This
   * way, it is easy to produce
   * plots of the different levels
   * of grid.
   *
   * @p mapping is a pointer to a
   * mapping used for the
   * transformation of cells at the
   * boundary. If zero, then use
   * standard Q1 mapping.
   *
   * Names and values of additional
   * flags controlling the output
   * can be found in the
   * documentation of the
   * GridOutFlags::Gnuplot() class.
   *
   * Not implemented for the
   * codimension one case.
   */
  template <int dim, int spacedim>
  void write_gnuplot (const Triangulation<dim,spacedim> &tria,
                      std::ostream           &out,
                      const Mapping<dim,spacedim> *mapping=0) const;

  /**
   * Write the triangulation in the msh format.
   *
   * Msh is the format used by Gmsh and it is described in the gmsh
   * user's guide. Besides the usual output of the grid only, you can
   * decide through additional flags (see below, and the documentation
   * of the GridOutFlags::Msh() class) whether boundary faces with
   * non-zero boundary indicator shall be written to the file
   * explicitly. This is useful, if you want to re-read the grid later
   * on, since <tt>deal.II</tt> sets the boundary indicator to zero by
   * default; therefore, to obtain the same triangulation as before,
   * you have to specify faces with differing boundary indicators
   * explicitly, which is done by this flag.
   *
   * Names and values of further flags controlling the output can be
   * found in the documentation of the GridOut::Msh() class.
   *
   * Works also in the codimension one case.
   */
  template <int dim, int spacedim>
  void write_msh (const Triangulation<dim,spacedim> &tria,
                  std::ostream             &out) const;

  /**
   * Write the triangulation in the ucd format.
   *
   * UCD (unstructured cell data) is the format used by AVS and some
   * other programs. It is described in the AVS developer's
   * guide. Besides the usual output of the grid only, you can decide
   * through additional flags (see below, and the documentation of the
   * GridOutFlags::Ucd() class) whether boundary faces with non-zero
   * boundary indicator shall be written to the file explicitly. This
   * is useful, if you want to re-read the grid later on, since
   * <tt>deal.II</tt> sets the boundary indicator to zero by default;
   * therefore, to obtain the same triangulation as before, you have
   * to specify faces with differing boundary indicators explicitly,
   * which is done by this flag.
   *
   * Names and values of further flags controlling the output can be
   * found in the documentation of the GridOut::Ucd() class.
   *
   * Works also for the codimension one case.
   */
  template <int dim, int spacedim>
  void write_ucd (const Triangulation<dim,spacedim> &tria,
                  std::ostream             &out) const;

  /**
   * Write the triangulation in the encapsulated postscript format.
   *
   * In this format, each line of the triangulation is written
   * separately. We scale the picture such that either x-values or
   * y-values range between zero and a fixed size. The other axis is
   * scaled by the same factor. Which axis is taken to compute the
   * scale and the size of the box it shall fit into is determined by
   * the output flags (see below, and the documentation of the
   * GridOutFlags::Eps() class).
   *
   * The bounding box is close to the triangulation on all four sides,
   * without an extra frame. The line width is chosen to be 0.5 by
   * default, but can be changed. The line width is to be compared
   * with the extension of the picture, of which the default is 300.
   *
   * The flag @p color_lines_on_user_flag allows to draw lines with
   * the @p user_flag set to be drawn in red. The colors black and red
   * are defined as @p b and @p r in the preamble of the output file
   * and can be changed there according to need.
   *
   * @p mapping is a pointer to a mapping used for the transformation
   * of cells at the boundary. If zero, then use standard Q1 mapping.
   *
   * Names and values of additional flags controlling the output can
   * be found in the documentation of the GridOutFlags::Eps()
   * class. Especially the viewpoint for three dimensional grids is of
   * importance here.
   *
   * Not implemented for the codimension one case.
   */
  template <int dim, int spacedim>
  void write_eps (const Triangulation<dim, spacedim> &tria,
                  std::ostream             &out,
                  const Mapping<dim, spacedim>       *mapping=0) const;

  /**
   * Write two-dimensional XFig-file.
   *
   * This function writes all grid cells as polygons and optionally
   * boundary lines. Several parameters can be adjusted by the
   * XFigFlags control object.
   *
   * If levels are coded to depth, the complete grid hierarchy is
   * plotted with fine cells before their parents. This way, levels
   * can be switched on and off in xfig by selecting levels.
   *
   * Polygons are either at depth 900-level or at 900+@p material_id,
   * depending on the flag @p level_depth. Accordingly, boundary edges
   * are at depth 800-level or at 800+@p boundary_id. Therefore,
   * boundary edges are alway in front of cells.
   *
   * Not implemented for the codimension one case.
   */
  template <int dim, int spacedim>
  void write_xfig (const Triangulation<dim, spacedim> &tria,
                   std::ostream              &out,
                   const Mapping<dim, spacedim>        *mapping=0) const;

  /**
   * Write the triangulation in the SVG format.
   *
   * SVG (Scalable Vector Graphics) is
   * an XML-based vector image format
   * developed and maintained by the
   * World Wide Web Consortium (W3C).
   * This function conforms to the
   * latest specification SVG 1.1,
   * released on August 16, 2011.
   *
   * The cells of the triangulation are written as polygons with
   * additional lines at the boundary of the triangulation. A coloring
   * of the cells is further possible in order to visualize a certain
   * property of the cells such as their level or material id. A
   * colorbar can be drawn to encode the chosen coloring.  Moreover, a
   * cell label can be added, showing level index, etc.
   *
   * @note This function is currently only implemented for
   * two-dimensional grids in two
   * space dimensions.
   */
  void write_svg (const Triangulation<2,2> &tria,
                  std::ostream             &out) const;

  /**
   * Declaration of the same function as above for all other dimensions and
   * space dimensions. This function is not currently implemented and is only
   * declared to exist to support dimension independent programming.
   */
  template <int dim, int spacedim>
  void write_svg (const Triangulation<dim,spacedim> &tria,
                  std::ostream                      &out) const;


  /**
   * Write triangulation in MathGL script format. To interpret this
   * file a version of MathGL>=2.0.0 is required.
   *
   * To get a handle on the resultant MathGL script within a graphical
   * environment an interpreter is needed. A suggestion to start with
   * is <code>mglview</code>, which is bundled with MathGL.
   * <code>mglview</code> can interpret and display small-to-medium
   * MathGL scripts in a graphical window and enables conversion to
   * other formats such as EPS, PNG, JPG, SVG, as well as view/display
   * animations. Some minor editing, such as modifying the lighting or
   * alpha channels, can also be done.
   *
   * @note Not implemented for the codimension one case.
   */
  template <int dim, int spacedim>
  void write_mathgl (const Triangulation<dim, spacedim> &tria,
                     std::ostream             &out) const;

  /**
   * Write triangulation in VTK format.
   */
  template <int dim, int spacedim>
  void write_vtk (const Triangulation<dim,spacedim> &tria,
                  std::ostream                      &out) const;

  /**
   * Write triangulation in VTU format.
   */
  template <int dim, int spacedim>
  void write_vtu (const Triangulation<dim,spacedim> &tria,
                  std::ostream                      &out) const;

  /**
   * Write grid to @p out according to the given data format. This
   * function simply calls the appropriate <tt>write_*</tt> function.
   */
  template <int dim, int spacedim>
  void write (const Triangulation<dim,spacedim> &tria,
              std::ostream                      &out,
              const OutputFormat                 output_format,
              const Mapping<dim,spacedim>       *mapping=0) const;

  /**
   * Write mesh in default format set by ParameterHandler.
   */
  template <int dim, int spacedim>
  void write (const Triangulation<dim,spacedim> &tria,
              std::ostream                      &out,
              const Mapping<dim,spacedim>       *mapping=0) const;

  /**
   * Set flags for DX output
   */
  void set_flags (const GridOutFlags::DX &flags);

  /**
   * Set flags for GMSH output
   */
  void set_flags (const GridOutFlags::Msh &flags);

  /**
   * Set flags for UCD output
   */
  void set_flags (const GridOutFlags::Ucd &flags);

  /**
   * Set flags for GNUPLOT output
   */
  void set_flags (const GridOutFlags::Gnuplot &flags);

  /**
   * Set flags for EPS output of a
   * one-dimensional triangulation
   */
  void set_flags (const GridOutFlags::Eps<1> &flags);

  /**
   * Set flags for EPS output of a
   * two-dimensional triangulation
   */
  void set_flags (const GridOutFlags::Eps<2> &flags);

  /**
   * Set flags for EPS output of a three-dimensional triangulation
   */
  void set_flags (const GridOutFlags::Eps<3> &flags);

  /**
   * Set flags for EPS output of a three-dimensional triangulation
   */
  void set_flags (const GridOutFlags::XFig &flags);

  /**
   * Set flags for SVG output
   */
  void set_flags (const GridOutFlags::Svg &flags);

  /**
   * Set flags for MathGL output
   */
  void set_flags (const GridOutFlags::MathGL &flags);

  /**
   * Set flags for VTK output
   */
  void set_flags (const GridOutFlags::Vtk &flags);

  /**
   * Set flags for VTU output
   */
  void set_flags (const GridOutFlags::Vtu &flags);

  /**
   * Provide a function that can tell us which
   * suffix a given output format
   * usually has. For example, it defines the following mappings:
   * <ul>
   * <li> @p OpenDX: <tt>.dx</tt>
   * <li> @p gnuplot: <tt>.gnuplot</tt>
   * <li> @p ucd: <tt>.inp</tt>
   * <li> @p eps: <tt>.eps</tt>.
   * </ul>
   * Similar mappings are provided for all implemented formats.
   *
   * Since this function does not need data from this object, it is
   * static and can thus be called without creating an object of this
   * class.
   */
  static std::string default_suffix (const OutputFormat output_format);

  /**
   * Default suffix for the default output format selected throuw
   * ParameterHandler.
   */
  std::string default_suffix () const;

  /**
   * Return the @p OutputFormat value corresponding to the given
   * string. If the string does not match any known format, an
   * exception is thrown.
   *
   * Since this function does not need data from this object, it is
   * static and can thus be called without creating an object of this
   * class. Its main purpose is to allow a program to use any
   * implemented output format without the need to extend the
   * program's parser each time a new format is implemented.
   *
   * To get a list of presently available format names, e.g. to give
   * it to the ParameterHandler class, use the function
   * get_output_format_names().
   */
  static OutputFormat parse_output_format (const std::string &format_name);

  /**
   * Return a list of implemented output formats. The different names
   * are separated by vertical bar signs (<tt>`|'</tt>) as used by the
   * ParameterHandler classes.
   */
  static std::string get_output_format_names ();

  /**
   * Declare parameters in ParameterHandler.
   */
  static void declare_parameters (ParameterHandler &param);

  /**
   * Parse parameters of ParameterHandler.
   */
  void parse_parameters (ParameterHandler &param);

  /**
   * Determine an estimate for the memory consumption (in bytes) of
   * this object.
   */
  std::size_t memory_consumption () const;

  /**
   * Exception
   */
  DeclException0 (ExcInvalidState);

private:
  /**
   * The default output format, set by a ParameterHandler.
   */
  OutputFormat default_format;

  /**
   * Flags for OpenDX output.
   */
  GridOutFlags::DX dx_flags;

  /**
   * Flags for GMSH output. Can be changed by using the set_flags(const GridOutFlags::Msh&)
   * function.
   */
  GridOutFlags::Msh     msh_flags;

  /**
   * Flags for UCD output. Can be changed by using the set_flags(const GridOutFlags::Ucd&)
   * function.
   */
  GridOutFlags::Ucd     ucd_flags;

  /**
   * Flags to be used upon output of GNUPLOT data. Can be changed by
   * using the set_flags(const GridOutFlags::Gnuplot&) function.
   */
  GridOutFlags::Gnuplot gnuplot_flags;

  /**
   * Flags to be used upon output of EPS
   * data in one space dimension. Can be
   * changed by using the set_flags(const GridOutFlags::Eps<1>&)
   * function.
   */
  GridOutFlags::Eps<1>  eps_flags_1;

  /**
   * Flags to be used upon output of EPS
   * data in two space dimensions. Can be
   * changed by using the @p set_flags
   * function.
   */
  GridOutFlags::Eps<2>  eps_flags_2;

  /**
   * Flags to be used upon output of EPS
   * data in three space dimensions. Can be
   * changed by using the @p set_flags
   * function.
   */
  GridOutFlags::Eps<3>  eps_flags_3;

  /**
   * Flags used for XFig output.
   */
  GridOutFlags::XFig xfig_flags;

  /**
   * Flags used for Svg output.
   */
  GridOutFlags::Svg svg_flags;

  /**
   * Flags for MathGL output.
   */
  GridOutFlags::MathGL mathgl_flags;

  /**
   * Flags for VTK output.
   */
  GridOutFlags::Vtk vtk_flags;

  /**
   * Flags for VTU output.
   */
  GridOutFlags::Vtu vtu_flags;

  /**
   * Write the grid information about faces to @p out. Only those
   * faces are printed which are on the boundary and which have a
   * boundary indicator not equal to zero, since the latter is the
   * default for boundary faces.
   *
   * Since cells and faces are continuously numbered, the @p
   * starting_index for the numbering of the faces is passed also.
   *
   * This function unfortunately can not be included in the regular @p
   * write_msh function, since it needs special treatment for the case
   * <tt>dim==1</tt>, in which case the face iterators are
   * <tt>void*</tt>'s and lack the member functions which are
   * called. We would not actually call these functions, but the
   * compiler would complain anyway when compiling the function for
   * <tt>dim==1</tt>. Bad luck.
   */
  template <int dim, int spacedim>
  void write_msh_faces (const Triangulation<dim,spacedim> &tria,
                        const unsigned int        starting_index,
                        std::ostream             &out) const;

  /**
   * Declaration of the specialization
   * of above function for 1d. Does
   * nothing.
   */
  void write_msh_faces (const Triangulation<1,1> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  /**
   * Declaration of the specialization
   * of above function for 1d, 2sd. Does
   * nothing.
   */
  void write_msh_faces (const Triangulation<1,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  void write_msh_faces (const Triangulation<1,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;



  /**
   * Write the grid information about lines to @p out. Only those
   * lines are printed which are on the boundary and which have a
   * boundary indicator not equal to zero, since the latter is the
   * default for boundary faces.
   *
   * Since cells and faces are continuously numbered, the @p
   * starting_index for the numbering of the lines is passed also.
   *
   * This function unfortunately can not be included in the regular @p
   * write_msh function, since it needs special treatment for the case
   * <tt>dim==1</tt> and <tt>dim==2</tt>, in which case the edge
   * iterators are <tt>void*</tt>'s and lack the member functions
   * which are called. We would not actually call these functions, but
   * the compiler would complain anyway when compiling the function
   * for <tt>dim==1/2</tt>. Bad luck.
   */
  template <int dim, int spacedim>
  void write_msh_lines (const Triangulation<dim,spacedim> &tria,
                        const unsigned int        starting_index,
                        std::ostream             &out) const;

  /**
   * Declaration of the specialization of above function for 1d. Does
   * nothing.
   */
  void write_msh_lines (const Triangulation<1,1> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;

  /**
   * Declaration of the specialization of above function for 1d,
   * 2sd. Does nothing.
   */
  void write_msh_lines (const Triangulation<1,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  void write_msh_lines (const Triangulation<1,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  /**
   * Declaration of the specialization of above function for 2d. Does
   * nothing.
   */
  void write_msh_lines (const Triangulation<2,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  /**
   * Declaration of the specialization of above function for 2d,
   * 3sd. Does nothing.
   */
  void write_msh_lines (const Triangulation<2,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;

  /**
   * Write the grid information about faces to @p out. Only those
   * faces are printed which are on the boundary and which have a
   * boundary indicator not equal to zero, since the latter is the
   * default for boundary faces.
   *
   * Since cells and faces are continuously numbered, the @p
   * starting_index for the numbering of the faces is passed also.
   *
   * This function unfortunately can not be included in the regular @p
   * write_ucd function, since it needs special treatment for the case
   * <tt>dim==1</tt>, in which case the face iterators are
   * <tt>void*</tt>'s and lack the member functions which are
   * called. We would not actually call these functions, but the
   * compiler would complain anyway when compiling the function for
   * <tt>dim==1</tt>. Bad luck.
   */


  template <int dim, int spacedim>
  void write_ucd_faces (const Triangulation<dim,spacedim> &tria,
                        const unsigned int        starting_index,
                        std::ostream             &out) const;

  /**
   * Declaration of the specialization of above function for 1d. Does
   * nothing.
   */
  void write_ucd_faces (const Triangulation<1,1> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;

  /**
   * Declaration of the specialization of above function for 1d,
   * 2sd. Does nothing.
   */
  void write_ucd_faces (const Triangulation<1,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  void write_ucd_faces (const Triangulation<1,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;


  /**
   * Write the grid information about lines to @p out. Only those
   * lines are printed which are on the boundary and which have a
   * boundary indicator not equal to zero, since the latter is the
   * default for boundary lines.
   *
   * Since cells, faces and lines are continuously numbered, the @p
   * starting_index for the numbering of the faces is passed also.
   *
   * This function unfortunately can not be included in the regular @p
   * write_ucd function, since it needs special treatment for the case
   * <tt>dim==1/2</tt>, in which case the edge iterators are
   * <tt>void*</tt>'s and lack the member functions which are
   * called. We would not actually call these functions, but the
   * compiler would complain anyway when compiling the function for
   * <tt>dim==1/2</tt>. Bad luck.
   */


  template <int dim, int spacedim>
  void write_ucd_lines (const Triangulation<dim,spacedim> &tria,
                        const unsigned int        starting_index,
                        std::ostream             &out) const;

  /**
   * Declaration of the specialization of above function for 1d. Does
   * nothing.
   */
  void write_ucd_lines (const Triangulation<1,1> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  /**
   * Declaration of the specialization of above function for 1d,
   * 2sd. Does nothing.
   */
  void write_ucd_lines (const Triangulation<1,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  void write_ucd_lines (const Triangulation<1,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;


  /**
   * Declaration of the specialization of above function for 2d. Does
   * nothing.
   */
  void write_ucd_lines (const Triangulation<2,2> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;
  /**
   * Declaration of the specialization of above function for 2d,
   * 3sd. Does nothing.
   */
  void write_ucd_lines (const Triangulation<2,3> &tria,
                        const unsigned int      starting_index,
                        std::ostream           &out) const;


  /**
   * This function projects a three-dimensional point (Point<3> point)
   * onto a two-dimensional image plane, specified by the position of
   * the camera viewing system (Point<3> camera_position), camera
   * direction (Point<3> camera_position), camera horizontal (Point<3>
   * camera_horizontal, necessary for the correct alignment of the
   * later images), and the focus of the camera (float camera_focus).
   *
   * For SVG output of grids.
   */
  static Point<2> svg_project_point(Point<3> point,
                                    Point<3> camera_position,
                                    Point<3> camera_direction,
                                    Point<3> camera_horizontal,
                                    float camera_focus);

  /**
   * Return the number of faces in the triangulation which have a
   * boundary indicator not equal to zero. Only these faces are
   * explicitly printed in the <tt>write_*</tt> functions; all faces
   * with indicator numbers::internal_face_boundary_id are interior
   * ones and an indicator with value zero for faces at the boundary

   * are considered default.
   *
   * This function always returns an empty list in one dimension.
   *
   * The reason for this function is the same as for
   * write_ucd_faces(). See there for more information.
   */
  template <int dim, int spacedim>
  unsigned int n_boundary_faces (const Triangulation<dim,spacedim> &tria) const;

  /**
   * Declaration of the specialization of above function for
   * 1d. Simply returns zero.
   */
  unsigned int n_boundary_faces (const Triangulation<1,1> &tria) const;

  /**
   * Declaration of the specialization of above function for 1d,
   * 2sd. Simply returns zero.
   */
  unsigned int n_boundary_faces (const Triangulation<1,2> &tria) const;
  unsigned int n_boundary_faces (const Triangulation<1,3> &tria) const;

  /**
   * Return the number of lines in the triangulation which have a
   * boundary indicator not equal to zero. Only these lines are
   * explicitly printed in the <tt>write_*</tt> functions; all lines
   * with indicator numbers::internal_face_boundary_id are interior
   * ones and an indicator with value zero for faces at the boundary
   * are considered default.
   *
   * This function always returns an empty list in one and two
   * dimensions.
   *
   * The reason for this function is the same as for
   * write_ucd_faces(). See there for more information.
   */
  template <int dim, int spacedim>
  unsigned int n_boundary_lines (const Triangulation<dim,spacedim> &tria) const;

  /**
   * Declaration of the specialization of above function for
   * 1d. Simply returns zero.
   */
  unsigned int n_boundary_lines (const Triangulation<1,1> &tria) const;

  /**
   * Declaration of the specialization of above function for 1d,
   * 2sd. Simply returns zero.
   */
  unsigned int n_boundary_lines (const Triangulation<1,2> &tria) const;
  unsigned int n_boundary_lines (const Triangulation<1,3> &tria) const;

  /**
   * Declaration of the specialization of above function for
   * 2d. Simply returns zero.
   */
  unsigned int n_boundary_lines (const Triangulation<2,2> &tria) const;
  /**
   * Declaration of the specialization of above function for 2d,
   * 3sd. Simply returns zero.
   */
  unsigned int n_boundary_lines (const Triangulation<2,3> &tria) const;
};



DEAL_II_NAMESPACE_CLOSE

#endif
