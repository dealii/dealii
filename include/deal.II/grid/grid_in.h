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

#ifndef dealii_grid_in_h
#define dealii_grid_in_h


#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/observer_pointer.h>
#include <deal.II/base/point.h>

#include <deal.II/lac/vector.h>

#include <iostream>
#include <map>
#include <string>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
template <int dim, int spacedim>
DEAL_II_CXX20_REQUIRES((concepts::is_valid_dim_spacedim<dim, spacedim>))
class Triangulation;
template <int dim>
struct CellData;
#endif

/**
 * This class implements an input mechanism for mesh data. It allows to read a
 * mesh structure into a Triangulation object. At present, UCD (unstructured
 * cell data), DB Mesh, XDA, %Gmsh, Tecplot, UNV, VTK, ASSIMP, and Cubit
 * are supported as input format for grid data. Any numerical data other than
 * geometric (vertex locations) and topological (how vertices form cells,
 * faces, and edges) information is ignored, but the readers for the various
 * formats generally do read information that associates material ids or
 * boundary ids to cells or faces (see
 * @ref GlossMaterialId "this"
 * and
 * @ref GlossBoundaryIndicator "this"
 * glossary entry for more information).
 *
 * In practice, the list of formats this class supports is of course limited.
 * However, other people have written tools to convert mesh files from one
 * format to another, and if you have a mesh you cannot read with the functions
 * of this class, you may want to convert it into a format that is in fact
 * supported. One of these tools you may want to look up is
 * [meshio](https://github.com/nschloe/meshio).
 * A separate tool that can convert tetrahedral to hexahedral mesges is tethex,
 * available [here](https://github.com/martemyev/tethex).
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
 * <li> <tt>%Gmsh 1.0 mesh</tt> format: this format is used by the @p %Gmsh mesh
 * generator (see http://gmsh.info/). The documentation in the @p %Gmsh
 * manual explains how to generate meshes compatible with the deal.II library
 * (i.e. quads rather than triangles). In order to use this format, %Gmsh has
 * to output the file in the old format 1.0. This is done adding the line
 * "Mesh.MshFileVersion = 1" to the input file.
 *
 * <li> <tt>%Gmsh 2.0 mesh</tt> format: this is a variant of the above format.
 * The read_msh() function automatically determines whether an input file is
 * version 1 or version 2.
 *
 * <li> <tt>Tecplot</tt> format: this format is used by @p TECPLOT and often
 * serves as a basis for data exchange between different applications. Note,
 * that currently only the ASCII format is supported, binary data cannot be
 * read.
 *
 * <li> <tt>UNV</tt> format: this format is generated by the Salome mesh
 * generator, see http://www.salome-platform.org/ . The sections of this format
 * that the GridIn::read_unv function supports are (only) 2411, 2412 and 2467.
 * A detailed description of all the sections of the UNV file can be found in
 * the following links:
 * <ul>
 * <li> https://www.ceas3.uc.edu/sdrluff/
 * <li> https://docs.plm.automation.siemens.com/tdoc/nx/12/nx_help#uid:xid1128419:index_advanced:xid1404601:xid1404604
 * </ul>
 *
 * A Wikipedia page dedicated to Universal File Format is available here:
 * https://en.wikipedia.org/wiki/Universal_File_Format
 *
 * Note that Salome, let's say in 2d, can only make a quad mesh on an object
 * that has exactly 4 edges (or 4 pieces of the boundary). That means, that if
 * you have a more complicated object and would like to mesh it with quads,
 * you will need to decompose the object into >= 2 separate objects. Then 1)
 * each of these separate objects is meshed, 2) the appropriate groups of
 * cells and/or faces associated with each of these separate objects are
 * created, 3) a compound mesh is built up, and 4) all numbers that might be
 * associated with some of the internal faces of this compound mesh are
 * removed.
 *
 * @note It is essential that __all__ the mesh groups in Salome have purely
 * integer names for this function to work as expected. For reference, the
 * user guide of [OpenFCST](http://www.openfcst.mece.ualberta.ca/index.html)
 * (a project that uses deal.II) has a section on creating deal.II-compatible
 * meshes on Salome.
 *
 * <li> <tt>VTK</tt> format: VTK Unstructured Grid Legacy file reader
 * generator. The reader can handle only Unstructured Grid format of data at
 * present for 2d & 3d geometries. The documentation for the general legacy
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
 * plug-in script can be found on the deal.II wiki page under
 * <a href="https://github.com/dealii/dealii/wiki/Mesh-Input-And-Output">Mesh
 * Input and Output</a>.
 *
 * Alternatively, Cubit can generate ABAQUS files that can be read in via the
 * read_abaqus() function. This may be a better option for geometries with
 * complex boundary condition surfaces and multiple materials - information
 * which is currently not easily obtained through Cubit's python interface.
 *
 * </ul>
 *
 * <h3>Structure of input grid data.</h3>
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
 * grids in various input formats call the GridTools::consistently_order_cells()
 * function to bring the order of vertices that define the cells into an
 * ordering that satisfies the requirements of the Triangulation class. Be sure
 * to read the documentation of that class if you experience unexpected problems
 * when reading grids through this class.
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
 */

template <int dim, int spacedim = dim>
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
    /// Use read_tecplot()
    tecplot,
    /// Use read_vtk()
    vtk,
    /// Use read_vtu()
    vtu,
    /// Use read_assimp()
    assimp,
    /// Use read_exodusii()
    exodusii,
  };

  /**
   * Constructor.
   */
  GridIn();

  /**
   * Constructor. Attach this triangulation to be fed with the grid data.
   */
  GridIn(Triangulation<dim, spacedim> &tria);

  /**
   * Attach this triangulation to be fed with the grid data.
   */
  void
  attach_triangulation(Triangulation<dim, spacedim> &tria);

  /**
   * Read from the given stream. If no format is given,
   * GridIn::Format::Default is used.
   */
  void
  read(std::istream &in, Format format = Default);

  /**
   * Open the file given by the string and call the previous function
   * read() taking a std::istream argument.
   */
  void
  read(const std::string &in, Format format = Default);

  /**
   * Read grid data from a unstructured vtk file. The vtk file may contain
   * the following VTK cell types: VTK_HEXAHEDRON (12), VTK_TETRA (10),
   * VTK_QUAD (9), VTK_TRIANGLE (5), and VTK_LINE (3).
   *
   * Depending on the template dimension, only some of the above are accepted.
   *
   * In particular, in three dimensions, this function expects the file to
   * contain
   *
   * - VTK_HEXAHEDRON/VTK_TETRA cell types
   * - VTK_QUAD/VTK_TRIANGLE cell types, to specify optional boundary or
   *   interior quad faces
   * - VTK_LINE cell types, to specify optional boundary or interior edges
   *
   * In two dimensions:
   *
   * - VTK_QUAD/VTK_TRIANGLE cell types
   * - VTK_LINE cell types, to specify optional boundary or interior edges
   *
   * In one dimension
   *
   * - VTK_LINE cell types
   *
   * The input file may specify boundary ids, material ids, and manifold ids
   * using the CELL_DATA section of the
   * [VTK file format](http://www.vtk.org/VTK/img/file-formats.pdf).
   *
   * This function interprets two types of CELL_DATA contained in the input
   * file: `SCALARS MaterialID`, used to specify the material_id of the cells,
   * or the boundary_id of the faces and edges, and `SCALARS ManifoldID`, that
   * can be used to specify the manifold id of any Triangulation object (cell,
   * face, or edge).
   *
   * The companion GridOut::write_vtk function can be used to write VTK files
   * compatible with this method.
   *
   * Also see
   * @ref simplex "Simplex support".
   */
  void
  read_vtk(std::istream &in);

  /**
   * Read grid data from a unstructured vtu file, saved by deal.II using
   * GridOut::write_vtu(), with the flag
   * GridOutFlags::Vtu::serialize_triangulation set to true.
   *
   * Notice that this function does not support reading in arbitrary vtu files,
   * but only files that were written by deal.II itself, using the function
   * GridOut::write_vtu and setting GridOutFlags::Vtu::serialize_triangulation
   * to true.
   *
   * When this flag is set to true, the generated vtu file contains the
   * triangulation in a xml section which is ignored by general vtu readers.
   * If this section is absent, an exception is thrown.
   */
  void
  read_vtu(std::istream &in);


  /**
   * Read grid data from an unv file as generated by the Salome mesh
   * generator. Numerical data is ignored.
   *
   * Note the comments on generating this file format in the general
   * documentation of this class.
   */
  void
  read_unv(std::istream &in);

  /**
   * Read grid data from an ucd file. Numerical data is ignored.
   * It is not possible to use a ucd file to set both boundary_id and
   * manifold_id for the same cell. Yet it is possible to use
   * the flag apply_all_indicators_to_manifolds to decide if
   * the indicators in the file refer to manifolds (flag set to true)
   * or boundaries (flag set to false). If the flag is set, the
   * indicators are used for cells as manifold id, too.
   */
  void
  read_ucd(std::istream &in,
           const bool    apply_all_indicators_to_manifolds = false);

  /**
   * Read grid data from an Abaqus file. Numerical and constitutive data is
   * ignored. As in the case of the ucd file format, it is possible to use
   * the flag apply_all_indicators_to_manifolds to decide if
   * the indicators in the file refer to manifolds (flag set to true)
   * or boundaries (flag set to false).
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
   * - Files generated using Cubit 11.x, 12.x, 13.x, 14.x and 15.x are valid,
   * but only when using a specific set of export steps. These are as follows:
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
  void
  read_abaqus(std::istream &in,
              const bool    apply_all_indicators_to_manifolds = false);

  /**
   * Read grid data from a file containing data in the DB mesh format.
   */
  void
  read_dbmesh(std::istream &in);

  /**
   * Read grid data from a file containing data in the XDA format.
   */
  void
  read_xda(std::istream &in);

  /**
   * Read grid data from an msh file. The %Gmsh formats are documented at
   * http://www.gmsh.info/.
   *
   * Also see
   * @ref simplex "Simplex support".
   */
  void
  read_msh(std::istream &in);

#ifdef DEAL_II_GMSH_WITH_API
  /**
   * Read grid data using Gmsh API. Any file supported by Gmsh can be passed as
   * argument. The format is deduced from the filename extension.
   *
   * This function interprets non-named physical ids (gmsh format < 4.0) as
   * material or boundary ids (similarly to what happens with the other
   * read_msh() function). If you want to specify non default manifold or
   * boundary ids, you must group all entities that require a non default
   * boundary or manifold id into named physical groups, where the name is
   * interpreted using the function Patterns::Tools::to_value() applied to a
   * `std::map<std::string, int>`. The keys can be either `MaterialID` (if the
   * physical group refers to object of dimension `dim`), `BoundaryID` (if the
   * group refers to objects of dimension < `dim`), or `ManifoldID`.
   *
   * From the Gmsh documentation, the formats of the physical tags follows the
   * following conventions:
   * @code
   * \$PhysicalNames // same as MSH version 2
   *   numPhysicalNames(ASCII int)
   *   dimension(ASCII int) physicalTag(ASCII int) "name"(127 characters max)
   *   ...
   * \$EndPhysicalNames
   * @endcode
   *
   * For example, the following snippet of mesh file
   * @code
   * MeshFormat
   * 4.1 0 8
   * \$EndMeshFormat
   * \$PhysicalNames
   * 4
   * 1 1 "ManifoldID:0"
   * 1 2 "BoundaryID: -1, ManifoldID: 1"
   * 2 3 "ManifoldID: 1"
   * 2 4 "MaterialID: 2, ManifoldID: 1"
   * \$EndPhysicalNames
   * \$Entities
   * ...
   * @endcode
   *
   * refers to a two dimensional grid where:
   * - a portion of the boundary of dimension 1 has physical tag 1, and manifold
   *   id 0
   * - some internal faces (lines of dimension 1) have manifold id 1
   * - some elements have manifold id 1 (and material id equal to the default
   *   value, i.e., zero)
   * - some elements have manifold id 1 and material id equal to 2
   *
   * If the physical groups are not named, then the behavior is the same as
   * the other read_msh() function, i.e., the physical tag itself is interpreted
   * as a boundary or material id.  Physical surface numbers created in Gmsh,
   * which can be seen in the .geo file, become material IDs.
   *
   *
   * Also see
   * @ref simplex "Simplex support".
   */
  void
  read_msh(const std::string &filename);
#endif

  /**
   * Read grid data from a `.mphtxt` file. `.mphtxt` is one of the file formats
   * typically generated by COMSOL. The file format is described at
   * http://victorsndvg.github.io/FEconv/formats/mphtxt.xhtml .
   *
   * The reader interprets the "geometric entity indicators" that COMSOL
   * writes into these files as either boundary indicators (for edges and faces
   * of cells) or as material ids (for cells). See the glossary for a
   * description of
   * @ref GlossBoundaryIndicator "boundary indicators".
   * and
   * @ref GlossMaterialId "material indicators"
   *
   * COMSOL has a habit of assigning "geometric entity indicators" not only
   * to edges and faces on the actual boundary, but also to interior faces
   * and edges. For example, for the following volume mesh generated by
   * COMSOL,
   * @image html "comsol-mesh-boundary-volume-mesh.png"
   * the marked edges and faces are as follows:
   * @image html "comsol-mesh-marked-lines.png"
   * @image html "comsol-mesh-marked-triangles.png"
   * Here, some of the marked lines and faces with explicitly given
   * geometric entity indicators are in the *interior* of the domain -- an
   * artifact of the geometry description that was used to describe
   * the mesh. However, we can of course not assign boundary indicators to
   * interior edges and faces. As a consequence, this reader function simply
   * ignores the geometric entity indicator for edges and faces that
   * are not in fact on the boundary of the domain. The result is then a mesh
   * in which only the following edges and faces are explicitly assigned
   * boundary indicators:
   * @image html "comsol-mesh-boundary-lines.png"
   * @image html "comsol-mesh-boundary-triangles.png"
   *
   * Also see
   * @ref simplex "Simplex support".
   */
  void
  read_comsol_mphtxt(std::istream &in);

  /**
   * Read grid data from a file containing tecplot ASCII data. This also works
   * in the absence of any tecplot installation.
   */
  void
  read_tecplot(std::istream &in);

  /**
   * Read in a file supported by Assimp, and generate a Triangulation
   * out of it.  If you specify a @p mesh_index, only the mesh with
   * the given index will be extracted, otherwise all meshes which are
   * present in the file will be used to generate the Triangulation.
   *
   * This function can only be used to read two-dimensional meshes (possibly
   * embedded in three dimensions). This is the standard for graphical software
   * such as blender, or 3d studio max, and that is what the original Assimp
   * library was built for. We "bend" it to deal.II to support complex
   * co-dimension one meshes and complex two-dimensional meshes.
   *
   * If @p remove_duplicates is set to true (the default), then
   * duplicated vertices will be removed if their distance is lower
   * than @p tol.
   *
   * Only the elements compatible with the given dimension and space dimension
   * will be extracted from the mesh, and only those elements that are
   * compatible with deal.II are supported. If you set
   * `ignore_unsupported_element_types`, all the other element types are simply
   * ignored by this algorithm. If your mesh contains a mixture of triangles
   * and quadrilaterals, for example, only the quadrilaterals will be
   * extracted. The resulting mesh (as represented in the Triangulation object)
   * may not make any sense if you are mixing compatible and incompatible
   * element types. If `ignore_unsupported_element_types` is set to `false`,
   * then an exception is thrown when an unsupported type is encountered.
   *
   * @param filename The file to read from
   * @param mesh_index Index of the mesh within the file
   * @param remove_duplicates Remove duplicated vertices
   * @param tol Tolerance to use when removing vertices
   * @param ignore_unsupported_element_types Don't throw exceptions if we
   *        encounter unsupported types during parsing
   */
  void
  read_assimp(const std::string &filename,
              const unsigned int mesh_index = numbers::invalid_unsigned_int,
              const bool         remove_duplicates                = true,
              const double       tol                              = 1e-12,
              const bool         ignore_unsupported_element_types = true);

  /**
   * A structure containing some of the information provided by ExodusII that
   * doesn't have a direct representation in the Triangulation object.
   *
   * @note This struct exists to enable forward compatibility with future
   * versions of read_exodusii that may provide additional output data, but for
   * now it has a single field.
   */
  struct ExodusIIData
  {
    /**
     * A vector containing a mapping from deal.II boundary ids (or manifold ids)
     * to the provided ExodusII sideset ids.
     */
    std::vector<std::vector<int>> id_to_sideset_ids;
  };

  /**
   * Read in a mesh stored in the ExodusII file format.
   *
   * ExodusII is a feature-rich file format that supports many more features
   * (like node sets, finite element fields, quality assurance data, and more)
   * than most other grid formats supported by this class. Many of these
   * features do not have equivalent representations in deal.II and are
   * therefore not supported (for example, deal.II does not assign degrees of
   * freedom directly to nodes, so data stored in a nodal format is not loaded
   * by this function). At the current time only the following information is
   * extracted from the input file:
   *
   * <ol>
   *   <li>Block ids: the block id of an element is loaded as its material
   *     id.</li>
   *   <li>Elements and vertices: the core geometric information stored in the
   *     ExodusII file populates the attached Triangulation object. Higher-order
   *     elements are automatically truncated to lower-order elements since
   *     deal.II does not support this feature (e.g., there is no equivalent to
   *     the <code>QUAD9</code> element in deal.II since all quadrilaterals have
   *     four vertices and additional geometric information is either stored in
   *     a Manifold or something like MappingQEulerian).</li>
   *   <li>Sideset ids: these are interpreted as boundary ids or manifold ids
   *     (see the note on the output value below). An error will occur if you
   *     attempt to read an ExodusII file that assigns a sideset id to an
   *     internal face boundary id.</li>
   * </ol>
   *
   * Sideset ids are not translated for Triangulations with nonzero codimension
   * since those Triangulations do not support the setting of boundary ids.
   *
   * @param filename The name of the file to read from.
   *
   * @param apply_all_indicators_to_manifolds Boolean determining if the sideset
   * ids should be interpreted as manifold ids or boundary ids. The default
   * value is <tt>false</tt>, i.e., treat all sideset ids as boundary ids. If
   * your mesh sets sideset ids on internal faces then it will be necessary to
   * set this argument to <code>true</code> and then do some postprocessing to
   * set the boundary ids correctly.
   *
   * @return This function returns a struct containing some extra data stored by
   * the ExodusII file that cannot be loaded into a Triangulation - see
   * ExodusIIData for more information.
   *
   * A cell face in ExodusII can be in an arbitrary number of sidesets (i.e., it
   * can have an arbitrary number of sideset ids) - however, a boundary cell
   * face in deal.II has exactly one boundary id. All boundary faces that are
   * not in a sideset are given the (default) boundary id of $0$. This function
   * then groups sidesets together into unique sets and gives each one a
   * boundary id. For example: Consider a single-quadrilateral mesh whose left
   * side has no sideset id, right side has sideset ids $0$ and $1$, and whose
   * bottom and top sides have sideset ids of $0$. The left face will have a
   * boundary id of $0$, the top and bottom faces boundary ids of $1$, and the
   * right face a boundary id of $2$. Hence the vector returned by this function
   * in that case will be $\{\{\}, \{0\}, \{0, 1\}\}$.
   */
  ExodusIIData
  read_exodusii(const std::string &filename,
                const bool         apply_all_indicators_to_manifolds = false);

  /**
   * Return the standard suffix for a file in this format.
   */
  static std::string
  default_suffix(const Format format);

  /**
   * Return the enum Format for the format name.
   */
  static Format
  parse_format(const std::string &format_name);

  /**
   * Return a list of implemented input formats. The different names are
   * separated by vertical bar signs (<tt>`|'</tt>) as used by the
   * ParameterHandler classes.
   */
  static std::string
  get_format_names();

  /**
   * Return a map containing cell data associated with the elements of an
   * external vtk format mesh imported using read_vtk().
   * The format of the returned map is as
   * follows:
   * - std::string stores the name of the field data (identifier) as specified
   * in the external mesh
   * - Vector<double> stores value for the given identifier in each cell.
   * To access the value, use cell_data[name_field][cell->active_cell_index()].
   *
   * For example, if the vtk mesh contains field data "Density" defined on
   * cells, then `cell_data["Density"][0]` provides the density defined at cell
   * ID '0', which corresponds to index 0 of the vector. The length of the
   * vector in `cell_data["Density"]` equals the number of elements in the
   * coarse mesh.
   */
  const std::map<std::string, Vector<double>> &
  get_cell_data() const;

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
  DeclException1(ExcUnknownIdentifier,
                 std::string,
                 << "The identifier <" << arg1 << "> as name of a "
                 << "part in an UCD input file is unknown or the "
                 << "respective input routine is not implemented."
                 << "(Maybe the space dimension of triangulation and "
                 << "input file do not match?");
  /**
   * Exception
   */
  DeclExceptionMsg(ExcNoTriangulationSelected,
                   "No Triangulation has been attached to this GridIn object "
                   "so that nothing can be filled during any read function "
                   "calls.  Please pass a reference to the Triangulation tria "
                   "to be  filled in the constructor GridIn(tria) or attach "
                   "it with the function call GridIn::attach_triangulation().");
  /**
   * Exception
   */
  DeclException2(
    ExcInvalidVertexIndex,
    int,
    int,
    << "While creating cell " << arg1
    << ", you are referencing a vertex with index " << arg2
    << " but no vertex with this index has been described in the input file.");
  /**
   * Exception
   */
  DeclException3(
    ExcInvalidVertexIndexGmsh,
    int,
    int,
    int,
    << "While creating cell " << arg1 << " (which is numbered as " << arg2
    << " in the input file), you are referencing a vertex with index " << arg3
    << " but no vertex with this index has been described in the input file.");
  /**
   * Exception
   */
  DeclException0(ExcInvalidDBMeshFormat);
  /**
   * Exception
   */
  DeclException1(ExcInvalidDBMESHInput,
                 std::string,
                 << "The string <" << arg1
                 << "> is not recognized at the present"
                 << " position of a DB Mesh file.");

  /**
   * Exception
   */
  DeclException1(
    ExcDBMESHWrongDimension,
    int,
    << "The specified dimension " << arg1
    << " is not the same as that of the triangulation to be created.");

  DeclException1(ExcInvalidGMSHInput,
                 std::string,
                 << "The string <" << arg1
                 << "> is not recognized at the present"
                 << " position of a Gmsh Mesh file.");

  DeclException1(ExcGmshUnsupportedGeometry,
                 int,
                 << "The Element Identifier <" << arg1 << "> is not "
                 << "supported in the deal.II library when "
                 << "reading meshes in " << dim << " dimensions.\n"
                 << "Supported elements are: \n"
                 << "ELM-TYPE\n"
                 << "1 Line (2 nodes, 1 edge).\n"
                 << "2 Triangle (3 nodes, 3 edges).\n"
                 << "3 Quadrilateral (4 nodes, 4 edges).\n"
                 << "4 Tetrahedron (4 nodes, 6 edges, 4 faces) when in 3d.\n"
                 << "5 Hexahedron (8 nodes, 12 edges, 6 faces) when in 3d.\n"
                 << "15 Point (1 node, ignored when read).");


  DeclException2(
    ExcGmshNoCellInformation,
    unsigned int,
    unsigned int,
    "While reading a gmsh file, the reader function did not find "
    "any cells. This sometimes happens if the file only contains a "
    "surface mesh, but not a volume mesh."
    "\n\n"
    "The reader function did find " +
      std::to_string(arg1) + " lines and " + std::to_string(arg2) +
      " facets (surface triangles or quadrilaterals).");

protected:
  /**
   * Store address of the triangulation to be fed with the data read in.
   */
  ObserverPointer<Triangulation<dim, spacedim>, GridIn<dim, spacedim>> tria;

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
   * the GridTools::consistently_order_cells() function).
   *
   * The output of this function consists of vectors for each line bounding
   * the cells indicating the direction it has with respect to the orientation
   * of this cell, and the cell number. The whole output is in a form such
   * that it can be read in by Gnuplot and generate the full plot without
   * further ado by the user.
   */
  static void
  debug_output_grid(const std::vector<CellData<dim>>   &cells,
                    const std::vector<Point<spacedim>> &vertices,
                    std::ostream                       &out);

private:
  /**
   * Skip empty lines in the input stream, i.e. lines that contain either
   * nothing or only whitespace.
   */
  static void
  skip_empty_lines(std::istream &in);

  /**
   * Skip lines of comment that start with the indicated character (e.g.
   * <tt>#</tt>) following the point where the given input stream presently
   * is. After the call to this function, the stream is at the start of the
   * first line after the comment lines, or at the same position as before if
   * there were no lines of comments.
   */
  static void
  skip_comment_lines(std::istream &in, const char comment_start);

  /**
   * This function does the nasty work (due to very lax conventions and
   * different versions of the tecplot format) of extracting the important
   * parameters from a tecplot header, contained in the string @p header. The
   * other variables are output variables, their value has no influence on the
   * function execution..
   */
  static void
  parse_tecplot_header(std::string               &header,
                       std::vector<unsigned int> &tecplot2deal,
                       unsigned int              &n_vars,
                       unsigned int              &n_vertices,
                       unsigned int              &n_cells,
                       std::vector<unsigned int> &IJK,
                       bool                      &structured,
                       bool                      &blocked);

  /**
   * Input format used by read() if no format is given.
   */
  Format default_format;

  /**
   * Data member that stores field data defined at the cells of the mesh.
   * The format is as follows:
   * - std::string stores the name of the field data (identifier) as specified
   * in the external mesh
   * - Vector<double> stores value for the given identifier in each cell id.
   *
   * To access the value use cell_data[name_field][cell->active_cell_index()].
   */
  std::map<std::string, Vector<double>> cell_data;
};

/* -------------- declaration of explicit specializations ------------- */

#ifndef DOXYGEN

template <>
void
GridIn<2>::debug_output_grid(const std::vector<CellData<2>> &cells,
                             const std::vector<Point<2>>    &vertices,
                             std::ostream                   &out);


template <>
void
GridIn<2, 3>::debug_output_grid(const std::vector<CellData<2>> &cells,
                                const std::vector<Point<3>>    &vertices,
                                std::ostream                   &out);
template <>
void
GridIn<3>::debug_output_grid(const std::vector<CellData<3>> &cells,
                             const std::vector<Point<3>>    &vertices,
                             std::ostream                   &out);
#endif // DOXYGEN

DEAL_II_NAMESPACE_CLOSE

#endif
