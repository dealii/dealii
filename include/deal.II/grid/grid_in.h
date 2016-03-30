// ---------------------------------------------------------------------
//
// Copyright (C) 1999 - 2016 by the deal.II authors
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

#ifndef dealii__grid_in_h
#define dealii__grid_in_h


#include <deal.II/base/config.h>
#include <deal.II/base/exceptions.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/point.h>
#include <iostream>
#include <vector>
#include <string>

DEAL_II_NAMESPACE_OPEN

template <int dim, int space_dim> class Triangulation;
template <int dim> struct CellData;

/**
 * This class implements an input mechanism for grid data. It allows to read a
 * grid structure into a triangulation object. At present, UCD (unstructured
 * cell data), DB Mesh, XDA, Gmsh, Tecplot, NetCDF, UNV, VTK, and Cubit are
 * supported as input format for grid data. Any numerical data other than
 * geometric (vertex locations) and topological (how vertices form cells,
 * faces, and edges) information is ignored, but the readers for the various
 * formats generally do read information that associates material ids or
 * boundary ids to cells or faces (see @ref GlossMaterialId "this" and
 * @ref GlossBoundaryIndicator "this" glossary entry for more information).
 *
 * @note Since deal.II only supports line, quadrilateral and hexahedral
 * meshes, the functions in this class can only read meshes that consist
 * exclusively of such cells. If you absolutely need to work with a mesh that
 * uses triangles or tetrahedra, then your only option is to convert the mesh
 * to quadrilaterals and hexahedra. A tool that can do this is tethex, see
 * http://code.google.com/p/tethex/wiki/Tethex .
 *
 * The mesh you read will form the coarsest level of a @p Triangulation
 * object. As such, it must not contain hanging nodes or other forms of
 * adaptive refinement, or strange things will happen if the mesh represented
 * by the input file does in fact have them. This is due to the fact that most
 * mesh description formats do not store neighborship information between
 * cells, so the grid reading functions have to regenerate it. They do so by
 * checking whether two cells have a common face. If there are hanging nodes
 * in a triangulation, adjacent cells have no common (complete) face, so the
 * grid reader concludes that the adjacent cells have no neighbors along these
 * faces and must therefore be at the boundary. In effect, an internal crack
 * of the domain is introduced this way. Since such cases are very hard to
 * detect (how is GridIn supposed to decide whether a place where the faces of
 * two small cells coincide with the face or a larger cell is in fact a
 * hanging node associated with local refinement, or is indeed meant to be a
 * crack in the domain?), the library does not make any attempt to catch such
 * situations, and you will get a triangulation that probably does not do what
 * you want. If your goal is to save and later read again a triangulation that
 * has been adaptively refined, then this class is not your solution; rather
 * take a look at the PersistentTriangulation class.
 *
 * To read grid data, the triangulation to be filled has to be empty.
 * Upon calling the functions of this class, the input file may
 * contain only lines in one dimension; lines and quads in two
 * dimensions; and lines, quads, and hexes in three dimensions. All
 * other cell types (e.g. triangles in two dimensions, triangles or
 * tetrahedra in 3d) are rejected.  (Here, the "dimension" refers to
 * the dimensionality of the mesh; it may be embedded in a higher
 * dimensional space, such as a mesh on the two-dimensional surface of
 * the sphere embedded in 3d, or a 1d mesh that discretizes a line in
 * 3d.) The result will be a triangulation that consists of the cells
 * described in the input file, and to the degree possible with
 * material indicators and boundary indicators correctly set as
 * described in the input file.
 *
 * @note You can not expect vertex and cell numbers in the triangulation
 * to match those in the input file. (This is already clear based on the
 * fact that we number cells and vertices separately, whereas this is not
 * the case for some input file formats; some formats also do not require
 * consecutive numbering, or start numbering at indices other than zero.)
 *
 *
 * <h3>Supported input formats</h3>
 *
 * At present, the following input formats are supported:
 * <ul>
 * <li> @p UCD (unstructured cell data) format: this format is used for grid
 * input as well as data output. If there are data vectors in the input file,
 * they are ignored, as we are only interested in the grid in this class. The
 * UCD format requires the vertices to be in following ordering: in 2d
 * @verbatim
 *      3-----2
 *      |     |
 *      |     |
 *      |     |
 *      0-----1
 * @endverbatim
 * and in 3d
 * @verbatim
 *         7-------6        7-------6
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     3   |       |    3-------2   |
 *     |   4-------5    |       |   5
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     0-------1        0-------1
 * @endverbatim
 * Note, that this ordering is different from the deal.II numbering scheme,
 * see the Triangulation class.  The exact description of the UCD format can
 * be found in the AVS Explorer manual (see http://www.avs.com).  The @p UCD
 * format can be read by the read_ucd() function.
 *
 * <li> <tt>DB mesh</tt> format: this format is used by the @p BAMG mesh
 * generator (see http://www-rocq.inria.fr/gamma/cdrom/www/bamg/eng.htm. The
 * documentation of the format in the @p BAMG manual is very incomplete, so we
 * don't actually parse many of the fields of the output since we don't know
 * their meaning, but the data that is read is enough to build up the mesh as
 * intended by the mesh generator. This format can be read by the
 * read_dbmesh() function.
 *
 * <li> @p XDA format: this is a rather simple format used by the MGF code. We
 * don't have an exact specification of the format, but the reader can read in
 * several example files. If the reader does not grok your files, it should be
 * fairly simple to extend it.
 *
 * <li> <tt>Gmsh 1.0 mesh</tt> format: this format is used by the @p GMSH mesh
 * generator (see http://www.geuz.org/gmsh/). The documentation in the @p GMSH
 * manual explains how to generate meshes compatible with the deal.II library
 * (i.e. quads rather than triangles). In order to use this format, Gmsh has
 * to output the file in the old format 1.0. This is done adding the line
 * "Mesh.MshFileVersion = 1" to the input file.
 *
 * <li> <tt>Gmsh 2.0 mesh</tt> format: this is a variant of the above format.
 * The read_msh() function automatically determines whether an input file is
 * version 1 or version 2.
 *
 * <li> <tt>Tecplot</tt> format: this format is used by @p TECPLOT and often
 * serves as a basis for data exchange between different applications. Note,
 * that currently only the ASCII format is supported, binary data cannot be
 * read.
 *
 * <li> <tt>UNV</tt> format: this format is generated by the Salome mesh
 * generator, see http://www.salome-platform.org/ . The sections of the format
 * that the GridIn::read_unv function supports are documented here:
 * <ul>
 * <li> section 2411: http://www.sdrl.uc.edu/universal-file-formats-for-modal-
 * analysis-testing-1/file-format-storehouse/unv_2411.htm
 * <li> section 2412: http://www.sdrl.uc.edu/universal-file-formats-for-modal-
 * analysis-testing-1/file-format-storehouse/unv_2412.htm
 * <li> section 2467: http://www.sdrl.uc.edu/universal-file-formats-for-modal-
 * analysis-testing-1/file-format-storehouse/unv_2467.htm
 * <li> all sections of this format, even if they may not be supported in our
 * reader, can be found here: http://www.sdrl.uc.edu/universal-file-formats-
 * for-modal-analysis-testing-1/file-format-storehouse/file-formats
 * </ul>
 * Note that Salome, let's say in 2D, can only make a quad mesh on an object
 * that has exactly 4 edges (or 4 pieces of the boundary). That means, that if
 * you have a more complicated object and would like to mesh it with quads,
 * you will need to decompose the object into >= 2 separate objects. Then 1)
 * each of these separate objects is meshed, 2) the appropriate groups of
 * cells and/or faces associated with each of these separate objects are
 * created, 3) a compound mesh is built up, and 4) all numbers that might be
 * associated with some of the internal faces of this compound mesh are
 * removed.
 *
 * <li> <tt>VTK</tt> format: VTK Unstructured Grid Legacy file reader
 * generator. The reader can handle only Unstructured Grid format of data at
 * present for 2D & 3D geometries. The documentation for the general legacy
 * vtk file, including Unstructured Grid format can be found here:
 * http://www.cacr.caltech.edu/~slombey/asci/vtk/vtk_formats.simple.html
 *
 * The VTK format requires the vertices to be in following ordering: in 2d
 * @verbatim
 *      3-----2
 *      |     |
 *      |     |
 *      |     |
 *      0-----1
 * @endverbatim
 * and in 3d
 * @verbatim
 *         7-------6        7-------6
 *        /|       |       /       /|
 *       / |       |      /       / |
 *      /  |       |     /       /  |
 *     4   |       |    4-------5   |
 *     |   3-------2    |       |   2
 *     |  /       /     |       |  /
 *     | /       /      |       | /
 *     |/       /       |       |/
 *     0-------1        0-------1
 * @endverbatim
 *
 *
 * <li> <tt>Cubit</tt> format: deal.II doesn't directly support importing from
 * Cubit at this time. However, Cubit can export in UCD format using a simple
 * plug-in, and the resulting UCD file can then be read by this class. The
 * plug-in script can be found on the deal.II wiki page,
 * http://code.google.com/p/dealii/wiki/MeshInputAndOutput .
 *
 * Alternatively, Cubit can generate ABAQUS files that can be read in via the
 * read_abaqus() function. This may be a better option for geometries with
 * complex boundary condition surfaces and multiple materials
 *  - information which is currently not easily obtained through
 * Cubit's python interface.
 *
 * </ul>
 *
 * <h3>Structure of input grid data. The GridReordering class</h3>
 *
 * It is your duty to use a correct numbering of vertices in the cell list,
 * i.e. for lines in 1d, you have to first give the vertex with the lower
 * coordinate value, then that with the higher coordinate value. For
 * quadrilaterals in two dimensions, the vertex indices in the @p quad list
 * have to be such that the vertices are numbered in counter-clockwise sense.
 *
 * In two dimensions, another difficulty occurs, which has to do with the
 * sense of a quadrilateral. A quad consists of four lines which have a
 * direction, which is by definition as follows:
 * @verbatim
 *   3-->--2
 *   |     |
 *   ^     ^
 *   |     |
 *   0-->--1
 * @endverbatim
 * Now, two adjacent cells must have a vertex numbering such that the
 * direction of the common side is the same. For example, the following two
 * quads
 * @verbatim
 *   3---4---5
 *   |   |   |
 *   0---1---2
 * @endverbatim
 * may be characterised by the vertex numbers <tt>(0 1 4 3)</tt> and <tt>(1 2
 * 5 4)</tt>, since the middle line would get the direction <tt>1->4</tt> when
 * viewed from both cells.  The numbering <tt>(0 1 4 3)</tt> and <tt>(5 4 1
 * 2)</tt> would not be allowed, since the left quad would give the common
 * line the direction <tt>1->4</tt>, while the right one would want to use
 * <tt>4->1</tt>, leading to an ambiguity. The Triangulation object is capable
 * of detecting this special case, which can be eliminated by rotating the
 * indices of the right quad by two. However, it would not know what to do if
 * you gave the vertex indices <tt>(4 1 2 5)</tt>, since then it would have to
 * rotate by one element or three, the decision which to take is not yet
 * implemented.
 *
 * There are more ambiguous cases, where the triangulation may not know what
 * to do at all without the use of sophisticated algorithms. Furthermore,
 * similar problems exist in three space dimensions, where faces and lines
 * have orientations that need to be taken care of.
 *
 * For this reason, the <tt>read_*</tt> functions of this class that read in
 * grids in various input formats call the GridReordering class to bring the
 * order of vertices that define the cells into an ordering that satisfies the
 * requirements of the Triangulation class. Be sure to read the documentation
 * of that class if you experience unexpected problems when reading grids
 * through this class.
 *
 *
 * <h3>Dealing with distorted mesh cells</h3>
 *
 * For each of the mesh reading functions, the last call is always to
 * Triangulation::create_triangulation(). That function checks whether all the
 * cells it creates as part of the coarse mesh are distorted or not (where
 * distortion here means that the Jacobian of the mapping from the reference
 * cell to the real cell has a non-positive determinant, i.e. the cell is
 * pinched or twisted; see the entry on
 * @ref GlossDistorted "distorted cells"
 * in the glossary). If it finds any such cells, it throws an exception. This
 * exception is not caught in the grid reader functions of the current class,
 * and so will propagate through to the function that called it. There, you
 * can catch and ignore the exception if you are certain that there is no harm
 * in dealing with such cells. If you were not aware that your mesh had such
 * cells, your results will likely be of dubious quality at best if you ignore
 * the exception.
 *
 *
 * @ingroup grid
 * @ingroup input
 * @author Wolfgang Bangerth, 1998, 2000, Luca Heltai, 2004, 2007, Jean-Paul
 * Pelteret 2015, Timo Heister 2015,  Krzysztof Bzowski, 2015
 */

template <int dim, int spacedim=dim>
class GridIn
{
public:
  /**
   * List of possible mesh input formats. These values are used when calling
   * the function read() in order to determine the actual reader to be called.
   */
  enum Format
  {
    /// Use GridIn::default_format stored in this object
    Default,
    /// Use read_unv()
    unv,
    /// Use read_ucd()
    ucd,
    /// Use read_abaqus()
    abaqus,
    /// Use read_dbmesh()
    dbmesh,
    /// Use read_xda()
    xda,
    /// Use read_msh()
    msh,
    /// Use read_netcdf()
    netcdf,
    /// Use read_tecplot()
    tecplot,
    /// Use read_vtk()
    vtk
  };

  /**
   * Constructor.
   */
  GridIn ();

  /**
   * Attach this triangulation to be fed with the grid data.
   */
  void attach_triangulation (Triangulation<dim,spacedim> &tria);

  /**
   * Read from the given stream. If no format is given,
   * GridIn::Format::Default is used.
   */
  void read (std::istream &in, Format format=Default);

  /**
   * Open the file given by the string and call the previous function read().
   * This function uses the PathSearch mechanism to find files. The file class
   * used is <code>MESH</code>.
   */
  void read (const std::string &in, Format format=Default);

  /**
   * Read grid data from an vtk file. Numerical data is ignored.
   *
   * @author Mayank Sabharwal, Andreas Putz, 2013
   */
  void read_vtk(std::istream &in);

  /**
   * Read grid data from an unv file as generated by the Salome mesh
   * generator. Numerical data is ignored.
   *
   * Note the comments on generating this file format in the general
   * documentation of this class.
   */
  void read_unv(std::istream &in);

  /**
   * Read grid data from an ucd file. Numerical data is ignored.
   */
  void read_ucd (std::istream &in);

  /**
   * Read grid data from an Abaqus file. Numerical and constitutive data is
   * ignored.
   *
   * @note The current implementation of this mesh reader is suboptimal, and
   * may therefore be slow for large meshes.
   *
   * @note Usage tips for Cubit:
   * - Multiple material-id's can be defined in the mesh.
   * This is done by specifying blocksets in the pre-processor.
   * - Arbitrary surface boundaries can be defined in the mesh.
   * This is done by specifying sidesets in the pre-processor. In particular,
   * boundaries are not confined to just surfaces (in 3d) individual element
   * faces can be added to the sideset as well. This is useful when a boundary
   * condition is to be applied on a complex shape boundary that is difficult
   * to define using "surfaces" alone. Similar can be done in 2d.
   *
   * @note Compatibility information for this file format is listed below.
   * - Files generated in Abaqus CAE 6.12 have been verified to be
   * correctly imported, but older (or newer) versions of Abaqus may also
   * generate valid input decks.
   * - Files generated using Cubit 11.x, 12.x and 13.x are valid, but only
   * when using a specific set of export steps. These are as follows:
   *     - Go to "Analysis setup mode" by clicking on the disc icon in the
   * toolbar on the right.
   *     - Select "Export Mesh" under "Operation" by clicking on the
   * necessary icon in the toolbar on the right.
   *     - Select an output file. In Cubit version 11.0 and 12.0 it might be
   * necessary to click on the browse button and type it in the dialogue that
   * pops up.
   *     - Select the dimension to output in.
   *     - Tick the overwrite box.
   *     - If using Cubit v12.0 onwards, uncheck the box "Export using Cubit
   * ID's". An invalid file will encounter errors if this box is left checked.
   *     - Click apply.
   */
  void read_abaqus (std::istream &in);

  /**
   * Read grid data from a file containing data in the DB mesh format.
   */
  void read_dbmesh (std::istream &in);

  /**
   * Read grid data from a file containing data in the XDA format.
   */
  void read_xda (std::istream &in);

  /**
   * Read grid data from an msh file, either version 1 or version 2 of that
   * file format. The GMSH formats are documented at
   * http://www.geuz.org/gmsh/.
   *
   * @note The input function of deal.II does not distinguish between newline
   * and other whitespace. Therefore, deal.II will be able to read files in a
   * slightly more general format than Gmsh.
   */
  void read_msh (std::istream &in);

  /**
   * Read grid data from a NetCDF file. The only data format currently
   * supported is the <tt>TAU grid format</tt>.
   *
   * This function requires the library to be linked with the NetCDF library.
   */
  void read_netcdf (const std::string &filename);

  /**
   * Read grid data from a file containing tecplot ASCII data. This also works
   * in the absence of any tecplot installation.
   */
  void read_tecplot (std::istream &in);

  /**
   * Returns the standard suffix for a file in this format.
   */
  static std::string default_suffix (const Format format);

  /**
   * Return the enum Format for the format name.
   */
  static Format parse_format (const std::string &format_name);

  /**
   * Return a list of implemented input formats. The different names are
   * separated by vertical bar signs (<tt>`|'</tt>) as used by the
   * ParameterHandler classes.
   */
  static std::string get_format_names ();

  /**
   * Exception
   */
  DeclException1(ExcUnknownSectionType,
                 int,
                 << "The section type <" << arg1 << "> in an UNV "
                 << "input file is not implemented.");

  /**
   * Exception
   */
  DeclException1(ExcUnknownElementType,
                 int,
                 << "The element type <" << arg1 << "> in an UNV "
                 << "input file is not implemented.");

  /**
   * Exception
   */
  DeclException1 (ExcUnknownIdentifier,
                  std::string,
                  << "The identifier <" << arg1 << "> as name of a "
                  << "part in an UCD input file is unknown or the "
                  << "respective input routine is not implemented."
                  << "(Maybe the space dimension of triangulation and "
                  << "input file do not match?");
  /**
   * Exception
   */
  DeclException0 (ExcNoTriangulationSelected);
  /**
   * Exception
   */
  DeclException2 (ExcInvalidVertexIndex,
                  int, int,
                  << "While creating cell " << arg1
                  << ", you are referencing a vertex with index " << arg2
                  << " but no vertex with this index has been described in the input file.");
  /**
   * Exception
   */
  DeclException0 (ExcInvalidDBMeshFormat);
  /**
   * Exception
   */
  DeclException1 (ExcInvalidDBMESHInput,
                  std::string,
                  << "The string <" << arg1 << "> is not recognized at the present"
                  << " position of a DB Mesh file.");

  /**
   * Exception
   */
  DeclException1 (ExcDBMESHWrongDimension,
                  int,
                  << "The specified dimension " << arg1
                  << " is not the same as that of the triangulation to be created.");

  DeclException1 (ExcInvalidGMSHInput,
                  std::string,
                  << "The string <" << arg1 << "> is not recognized at the present"
                  << " position of a Gmsh Mesh file.");

  DeclException1 (ExcGmshUnsupportedGeometry,
                  int,
                  << "The Element Identifier <" << arg1 << "> is not "
                  << "supported in the deal.II library when "
                  << "reading meshes in " << dim << " dimensions.\n"
                  << "Supported elements are: \n"
                  << "ELM-TYPE\n"
                  << "1 Line (2 nodes, 1 edge).\n"
                  << "3 Quadrilateral (4 nodes, 4 edges).\n"
                  << "5 Hexahedron (8 nodes, 12 edges, 6 faces) when in 3d.\n"
                  << "15 Point (1 node, ignored when read)");


  DeclException0 (ExcGmshNoCellInformation);
protected:
  /**
   * Store address of the triangulation to be fed with the data read in.
   */
  SmartPointer<Triangulation<dim,spacedim>,GridIn<dim,spacedim> > tria;

  /**
   * This function can write the raw cell data objects created by the
   * <tt>read_*</tt> functions in Gnuplot format to a stream. This is
   * sometimes handy if one would like to see what actually was created, if it
   * is known that the data is not correct in some way, but the Triangulation
   * class refuses to generate a triangulation because of these errors. In
   * particular, the output of this class writes out the cell numbers along
   * with the direction of the faces of each cell. In particular the latter
   * information is needed to verify whether the cell data objects follow the
   * requirements of the ordering of cells and their faces, i.e. that all
   * faces need to have unique directions and specified orientations with
   * respect to neighboring cells (see the documentations to this class and
   * the GridReordering class).
   *
   * The output of this function consists of vectors for each line bounding
   * the cells indicating the direction it has with respect to the orientation
   * of this cell, and the cell number. The whole output is in a form such
   * that it can be read in by Gnuplot and generate the full plot without
   * further ado by the user.
   */
  static void debug_output_grid (const std::vector<CellData<dim> > &cells,
                                 const std::vector<Point<spacedim> > &vertices,
                                 std::ostream &out);

private:

  /**
   * Skip empty lines in the input stream, i.e. lines that contain either
   * nothing or only whitespace.
   */
  static void skip_empty_lines (std::istream &in);

  /**
   * Skip lines of comment that start with the indicated character (e.g.
   * <tt>#</tt>) following the point where the given input stream presently
   * is. After the call to this function, the stream is at the start of the
   * first line after the comment lines, or at the same position as before if
   * there were no lines of comments.
   */
  static void skip_comment_lines (std::istream    &in,
                                  const char  comment_start);

  /**
   * This function does the nasty work (due to very lax conventions and
   * different versions of the tecplot format) of extracting the important
   * parameters from a tecplot header, contained in the string @p header. The
   * other variables are output variables, their value has no influence on the
   * function execution..
   */
  static void parse_tecplot_header(std::string   &header,
                                   std::vector<unsigned int> &tecplot2deal,
                                   unsigned int  &n_vars,
                                   unsigned int  &n_vertices,
                                   unsigned int  &n_cells,
                                   std::vector<unsigned int> &IJK,
                                   bool          &structured,
                                   bool          &blocked);

  /**
   * Input format used by read() if no format is given.
   */
  Format default_format;
};

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
GridIn<2>::debug_output_grid (const std::vector<CellData<2> > &cells,
                              const std::vector<Point<2> >    &vertices,
                              std::ostream                    &out);


template <>
void
GridIn<2,3>::debug_output_grid (const std::vector<CellData<2> > &cells,
                                const std::vector<Point<3> >    &vertices,
                                std::ostream                    &out);
template <>
void
GridIn<3>::debug_output_grid (const std::vector<CellData<3> > &cells,
                              const std::vector<Point<3> >    &vertices,
                              std::ostream                    &out);

#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
