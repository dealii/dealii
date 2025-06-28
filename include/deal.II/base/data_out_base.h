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

#ifndef dealii_data_out_base_h
#define dealii_data_out_base_h


#include <deal.II/base/config.h>

#include <deal.II/base/geometry_info.h>
#include <deal.II/base/mpi_stub.h>
#include <deal.II/base/point.h>
#include <deal.II/base/table.h>

#include <deal.II/grid/reference_cell.h>

#include <deal.II/numerics/data_component_interpretation.h>

// To be able to serialize XDMFEntry
#include <boost/serialization/map.hpp>

#include <limits>
#include <ostream>
#include <string>
#include <tuple>
#include <typeinfo>
#include <vector>

DEAL_II_NAMESPACE_OPEN

// Forward declarations
#ifndef DOXYGEN
class ParameterHandler;
class XDMFEntry;
#endif

/**
 * This is a base class for output of data on meshes of very general form.
 * Output data is expected as a set of <tt>patches</tt> and written to the
 * output stream in the format expected by the visualization tool. For a list
 * of output formats, check the enumeration #OutputFormat. For each format
 * listed there, this class contains a function <tt>write_format</tt>, writing
 * the output. Refer to the documentation of those functions for details on a
 * certain format.
 *
 * <h3>Structure of the output data</h3>
 *
 * Data is not written with the deal.II mesh structure. Instead, it relies on
 * a set of <tt>patches</tt> created by a derived class (for example the
 * DataOut, DataOutStack, DataOutFaces, DataOutRotation, or MatrixOut
 * classes).  Each Patch describes a single logical cell of a mesh, possibly
 * subdivided a number of times to represent higher order polynomials defined
 * on this cell. To this end, a patch consists of a <tt>dim</tt>-dimensional
 * regular grid with the same number of grid points in each direction. In the
 * simplest case it may consist of the corner points of a single mesh cell.
 * For each point of this local grid, the Patch contains an arbitrary number
 * of data values, though the number of data sets must be the same for each
 * point on each patch.
 *
 * By offering this interface to the different output formats, it is simple to
 * extend this class to new formats without depending on such things as actual
 * triangulations and handling of data vectors. These things shall be provided
 * by derived class which have a user callable interface then.
 *
 * Inside each patch, the data is organized in the usual lexicographical
 * order, <i>x</i> running fastest, then <i>y</i> and <i>z</i>. Nodes are
 * stored in this order and cells as well. Each cell in 3d is stored such that
 * the front face is in the <i>xz</i>-plane. In order to enhance
 * intelligibility of this concept, the following two sections are kept from a
 * previous version of this documentation.
 *
 *
 * <h4>Patches</h4>
 *
 * Grids can be thought of as a collection of cells; if you want to write out
 * data on such a grid, you can do so by writing them one cell at a time. The
 * functions in this class therefore take a list of objects describing the
 * data on one cell each. This data for each cell usually consists of a list
 * of vertices for this cell, and a list of data values (for example solution
 * data, error information, etc) at each of these vertices.
 *
 * In some cases, this interface to a cell is too restricted, however. For
 * example, you may have higher order elements and printing the values at the
 * vertices only is not enough. For this reason, we not only provide writing
 * the data on the vertices only, but the data is organizes as a tensor
 * product grid on each cell. The parameter <tt>n_subdivisions</tt>, which is
 * given for each patch separately, denotes how often the cell is to be
 * divided for output; for example, <tt>n_subdivisions==1</tt> yields no
 * subdivision of the cell, <tt>n_subdivisions==2</tt> will produce a grid of
 * 3 times 3 points in two spatial dimensions and 3 times 3 times 3 points in
 * three dimensions, <tt>n_subdivisions==3</tt> will yield 4 times 4 (times 4)
 * points, etc. The actual location of these points on the patch will be
 * computed by a multilinear transformation from the vertices given for this
 * patch.  For cells at the boundary, a mapping might be used to calculate the
 * position of the inner points. In that case the coordinates are stored
 * inside the Patch, as they cannot be easily recovered otherwise.
 *
 * Given these comments, the actual data to be printed on this patch of points
 * consists of several data sets each of which has a value at each of the
 * patch points. For example with <tt>n_subdivisions==2</tt> in two space
 * dimensions, each data set has to provide nine values, and since the patch
 * is to be printed as a tensor product (or its transformation to the real
 * space cell), its values are to be ordered like <i>(x0,y0) (x0,y1) (x0,y2)
 * (x1,y0) (x1,y1) (x1,y2) (x2,y0) (x2,y1) (x2,y2)</i>, i.e. the z-coordinate
 * runs fastest, then the y-coordinate, then x (if there are that many space
 * directions).
 *
 *
 * <h4>Generalized patches</h4>
 *
 * In general, the patches as explained above might be too restricted. For
 * example, one might want to draw only the outer faces of a domain in a
 * three-dimensional computation, if one is not interested in what happens
 * inside. Then, the objects that should be drawn are two-dimensional in a
 * three-dimensional world. The Patch class and associated output functions
 * handle these cases. The Patch class therefore takes two template
 * parameters, the first, named <tt>dim</tt> denoting the dimension of the
 * object (in the above example, this would be two), while the second, named
 * <tt>spacedim</tt>, denotes the dimension of the embedding space (this would
 * be three). The corner points of a patch have the dimension of the space,
 * while their number is determined by the dimension of the patch. By default,
 * the second template parameter has the same value as the first, which would
 * correspond to outputting a cell, rather than a face or something else.
 *
 * <h3>DataOutBaseInterface</h3>
 *
 * The members of this namespace are not usually called from user code
 * directly. Rather, classes that use the functions declared here are
 * typically derived from DataOutInterface.
 *
 * The interface of this class basically consists of the declaration of a data
 * type describing a patch and a bunch of functions taking a list of patches
 * and writing them in one format or other to the stream. It is in the
 * responsibility of the derived classes to provide this list of patches. In
 * addition to the list of patches, a name for each data set may be given.
 *
 *
 * <h3>Querying interface</h3>
 *
 * This class also provides a few functions (parse_output_format(),
 * get_output_format_names(), default_suffix()) that can be used to query
 * which output formats this class supports. The provide a list of names for
 * all the formats we can output, parse a string and return an enum indicating
 * each format, and provide a way to convert a value of this enum into the
 * usual suffix used for files of that name. Using these functions, one can
 * entirely free applications from knowledge which formats the library
 * presently allows to output; several of the example programs show how to do
 * this.
 *
 * <h3>Output parameters</h3>
 *
 * All functions take a parameter which is a structure of type
 * <tt>XFlags</tt>, where <tt>X</tt> is the name of the output format. To find
 * out what flags are presently supported, read the documentation of the
 * different structures.
 *
 * Note that usually the output formats used for scientific visualization
 * programs have no or very few parameters (apart from some compatibility
 * flags) because there the actual appearance of output is determined using
 * the visualization program and the files produced by this class store more
 * or less only raw data.
 *
 * The direct output formats, like Postscript or Povray need to be given a lot
 * more parameters, though, since there the output file has to contain all
 * details of the viewpoint, light source, etc.
 *
 * <h3>Writing backends</h3>
 *
 * An abstraction layer has been introduced to facilitate coding backends for
 * additional visualization tools. It is applicable for data formats
 * separating the information into a field of vertices, a field of connection
 * information for the grid cells and data fields.
 *
 * For each of these fields, output functions are implemented, namely
 * write_nodes(), write_cells() and write_data(). In order to use these
 * functions, a format specific output stream must be written, following the
 * examples of DXStream, GmvStream, VtkStream and so on, implemented in the
 * .cc file.
 *
 * In this framework, the implementation of a new output format is reduced to
 * writing the section headers and the new output stream class for writing a
 * single mesh object.
 *
 * @ingroup output
 */
namespace DataOutBase
{
  /**
   * An enum for different levels of compression used in several places
   * to determine zlib compression levels for binary output. At some
   * places, it is possible to output the data also as plain text as
   * an alternative, which is convenient for debugging. We use
   * this flag to indicate such an output as well.
   */
  enum class CompressionLevel
  {
    /**
     * Do not use any compression.
     */
    no_compression,
    /**
     * Use the fastest available compression algorithm.
     */
    best_speed,
    /**
     * Use the algorithm which results in the smallest compressed files.
     */
    best_compression,
    /**
     * Use the default compression algorithm. This is a compromise between
     * speed and file size.
     */
    default_compression,
    /**
     * Output as plain text (ASCII) if available.
     */
    plain_text
  };


  /**
   * Data structure describing a patch of data in <tt>dim</tt> space
   * dimensions.
   *
   * A patch consists of the following data:
   * <ul>
   * <li> the corner #vertices,
   * <li> the number #n_subdivisions of the number of cells the Patch has in
   * each space direction,
   * <li> the #data attached to each vertex, in the usual lexicographic
   * ordering,
   * <li> information on #neighbors.
   * </ul>
   *
   * See the general documentation of the DataOutBase class for more
   * information on its contents and purposes.  In the case of two dimensions,
   * the next picture is an example of <tt>n_subdivisions</tt> = 4 because the
   * number of (sub)cells within each patch is equal to
   * <tt>2<sup>dim</sup></tt>.
   *
   * @ingroup output
   */
  template <int dim, int spacedim = dim>
  struct Patch
  {
    /**
     * Make the <tt>spacedim</tt> template parameter available.
     */
    static const unsigned int space_dim = spacedim;

    /**
     * Corner points of a patch.  Interior points are computed by a multilinear
     * transformation of the unit cell to the cell specified by these corner
     * points, if <code>points_are_available==false</code>.
     *
     * On the other hand, if <code>points_are_available==true</code>, then
     * the coordinates of the points at which output is to be generated
     * is attached in additional rows to the <code>data</code> table.
     *
     * The order of points is the same as for cells in the
     * triangulation.
     *
     * @note This array is sized to accommodate the maximal number of vertices
     *   one might encounter in a cell, namely those for a hypercube cell.
     *   For other kinds of cells (triangles, tetrahedra, etc.), only the
     *   first few elements of this array will be used.
     */
    std::array<Point<spacedim>, GeometryInfo<dim>::vertices_per_cell> vertices;

    /**
     * Patch indices of neighbors of the current patch. This is made available
     * for the OpenDX format that requires neighbor
     * information for advanced output.
     *
     * @note This array is sized to accommodate the maximal number of faces
     *   one might encounter in a cell, namely those for a hypercube cell.
     *   For other kinds of cells (triangles, tetrahedra, etc.), only the
     *   first few elements of this array will be used.
     */
    std::array<unsigned int, GeometryInfo<dim>::faces_per_cell> neighbors;

    /**
     * Number of this patch. Since we are not sure patches are always
     * handled in the same order, we better store this.
     */
    unsigned int patch_index;

    /**
     * Number of subdivisions with which this patch is to be written.
     * <tt>1</tt> means no subdivision, <tt>2</tt> means bisection, <tt>3</tt>
     * trisection, etc.
     */
    unsigned int n_subdivisions;

    /**
     * Data vectors. The format is as follows: <tt>data(i,.)</tt> denotes the
     * data belonging to the <tt>i</tt>th data vector. <tt>data.n_cols()</tt>
     * therefore equals the number of output points; this number is
     * <tt>(subdivisions+1)^{dim}</tt>. <tt>data.n_rows()</tt> equals the number
     * of data vectors. For the current purpose, a data vector equals one
     * scalar, even if multiple scalars may later be interpreted as vectors.
     *
     * Within each column, <tt>data(.,j)</tt> are the data values at the
     * output point <tt>j</tt>, where <tt>j</tt> denotes the usual
     * lexicographic ordering in deal.II. This is also the order of points as
     * provided by the <tt>QIterated</tt> class when used with the
     * <tt>QTrapezoid</tt> class as subquadrature.
     *
     * Since the number of data vectors is usually the same for all patches to
     * be printed, <tt>data.size()</tt> should yield the same value for all
     * patches provided. The exception are patches for which
     * points_are_available are set, where the actual coordinates of the point
     * are appended to the 'data' field, see the documentation of the
     * points_are_available flag.
     */
    Table<2, float> data;

    /**
     * A flag indicating whether the coordinates of the interior patch points
     * (assuming that the patch is supposed to be subdivided further) are
     * appended to the @p data table (@p true) or not (@p false). The latter
     * is the default and in this case the locations of the points interior to
     * this patch are computed by (bi-, tri-)linear interpolation from the
     * vertices of the patch.
     *
     * This option exists since patch points may be evaluated using a Mapping
     * (rather than by a linear interpolation) and therefore have to be stored
     * in the Patch structure.
     */
    bool points_are_available;

    /**
     * Reference-cell type of the underlying cell of this patch.
     */
    ReferenceCell reference_cell;

    /**
     * Default constructor. Sets #n_subdivisions to one, #points_are_available
     * to false, and #patch_index to #no_neighbor.
     */
    Patch();

    /**
     * Compare the present patch for equality with another one. This is used
     * in a few of the automated tests in our testsuite.
     */
    bool
    operator==(const Patch &patch) const;

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object. This is not exact (but will usually be close) because
     * calculating the memory usage of trees (e.g., <tt>std::map</tt>) is
     * difficult.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Swap the current object's contents with those of the given argument.
     */
    void
    swap(Patch<dim, spacedim> &other_patch) noexcept;

    /**
     * Value to be used if this patch has no neighbor on one side.
     */
    static const unsigned int no_neighbor = numbers::invalid_unsigned_int;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception
     */
    DeclException2(
      ExcInvalidCombinationOfDimensions,
      int,
      int,
      << "It is not possible to have a structural dimension of " << arg1
      << " to be larger than the space dimension of the surrounding"
      << " space " << arg2);
    /** @} */
  };



  /**
   * A specialization of the general Patch<dim,spacedim> template that is
   * tailored to the case of points, i.e., zero-dimensional objects embedded
   * in @p spacedim dimensional space.
   *
   * The current class is compatible with the general template to allow for
   * using the same functions accessing patches of arbitrary dimensionality
   * in a generic way. However, it makes some variables that are nonsensical
   * for zero-dimensional patches into @p static variables that exist only
   * once in the entire program, as opposed to once per patch. Specifically,
   * this is the case for the @p neighbors array and the @p n_subdivisions
   * member variable that make no sense for zero-dimensional patches because
   * points have no natural neighbors across their non-existent faces, nor
   * can they reasonably be subdivided.
   */
  template <int spacedim>
  struct Patch<0, spacedim>
  {
    /**
     * Make the <tt>spacedim</tt> template parameter available.
     */
    static const unsigned int space_dim = spacedim;

    /**
     * Corner points of a patch.  For the current class of zero-dimensional
     * patches, there is of course only a single vertex.
     *
     * If <code>points_are_available==true</code>, then
     * the coordinates of the point at which output is to be generated
     * is attached as an additional row to the <code>data</code> table.
     */
    Point<spacedim> vertices[1];

    /**
     * An unused, @p static variable that exists only to allow access
     * from general code in a generic fashion.
     */
    static unsigned int neighbors[1];

    /**
     * Number of this patch. Since we are not sure patches are always
     * handled in the same order, we better store this.
     */
    unsigned int patch_index;

    /**
     * Number of subdivisions with which this patch is to be written.
     * <tt>1</tt> means no subdivision, <tt>2</tt> means bisection, <tt>3</tt>
     * trisection, etc.
     *
     * Since subdivision makes no sense for zero-dimensional patches,
     * this variable is not used but exists only to allow access
     * from general code in a generic fashion.
     */
    static const unsigned int n_subdivisions;

    /**
     * Data vectors. The format is as follows: <tt>data(i,.)</tt> denotes the
     * data belonging to the <tt>i</tt>th data vector. <tt>data.n_cols()</tt>
     * therefore equals the number of output points; this number is
     * of course one for the current class, given that we produce output on
     * points. <tt>data.n_rows()</tt> equals the number of
     * data vectors. For the current purpose, a data vector equals one scalar,
     * even if multiple scalars may later be interpreted as vectors.
     *
     * Within each column, <tt>data(.,j)</tt> are the data values at the
     * output point <tt>j</tt>; for the current class, @p j can only
     * be zero.
     *
     * Since the number of data vectors is usually the same for all patches to
     * be printed, <tt>data.size()</tt> should yield the same value for all
     * patches provided. The exception are patches for which
     * points_are_available are set, where the actual coordinates of the point
     * are appended to the 'data' field, see the documentation of the
     * points_are_available flag.
     */
    Table<2, float> data;

    /**
     * A flag indicating whether the coordinates of the interior patch points
     * (assuming that the patch is supposed to be subdivided further) are
     * appended to the @p data table (@p true) or not (@p false). The latter
     * is the default and in this case the locations of the points interior to
     * this patch are computed by (bi-, tri-)linear interpolation from the
     * vertices of the patch.
     *
     * This option exists since patch points may be evaluated using a Mapping
     * (rather than by a linear interpolation) and therefore have to be stored
     * in the Patch structure.
     */
    bool points_are_available;

    /**
     * Reference-cell type of the underlying cell of this patch. Since for
     * zero-dimensional objects, a patch can only refer to a vertex, this
     * field is always equal to ReferenceCells::Vertex and can not be changed.
     */
    static const ReferenceCell reference_cell;

    /**
     * Default constructor. Sets #points_are_available
     * to false, and #patch_index to #no_neighbor.
     */
    Patch();

    /**
     * Compare the present patch for equality with another one. This is used
     * in a few of the automated tests in our testsuite.
     */
    bool
    operator==(const Patch &patch) const;

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object. This is not exact (but will usually be close) because
     * calculating the memory usage of trees (e.g., <tt>std::map</tt>) is
     * difficult.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Swap the current object's contents with those of the given argument.
     */
    void
    swap(Patch<0, spacedim> &other_patch) noexcept;

    /**
     * Value to be used if this patch has no neighbor on one side.
     */
    static const unsigned int no_neighbor = numbers::invalid_unsigned_int;

    /**
     * @addtogroup Exceptions
     * @{
     */

    /**
     * Exception
     */
    DeclException2(
      ExcInvalidCombinationOfDimensions,
      int,
      int,
      << "It is not possible to have a structural dimension of " << arg1
      << " to be larger than the space dimension of the surrounding"
      << " space " << arg2);
    /** @} */
  };


  /**
   * Base class describing common functionality between different output
   * flags.
   *
   * This is implemented with the "Curiously Recurring Template Pattern";
   * derived classes use their own type to fill in the typename so that
   * <tt>memory_consumption</tt> works correctly. See the Wikipedia page on
   * the pattern for more information.
   *
   * @ingroup output
   */
  template <typename FlagsType>
  struct OutputFlagsBase
  {
    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     *
     * This method does nothing, but child classes may override this method to
     * add fields to <tt>prm</tt>.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in declare_parameters() and set the flags
     * for this output format accordingly.
     *
     * This method does nothing, but child classes may override this method to
     * add fields to <tt>prm</tt>.
     */
    void
    parse_parameters(const ParameterHandler &prm);

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object. This is not exact (but will usually be close) because
     * calculating the memory usage of trees (e.g., <tt>std::map</tt>) is
     * difficult.
     */
    std::size_t
    memory_consumption() const;
  };


  template <typename FlagsType>
  void
  OutputFlagsBase<FlagsType>::declare_parameters(ParameterHandler &)
  {}


  template <typename FlagsType>
  void
  OutputFlagsBase<FlagsType>::parse_parameters(const ParameterHandler &)
  {}


  template <typename FlagsType>
  std::size_t
  OutputFlagsBase<FlagsType>::memory_consumption() const
  {
    return sizeof(FlagsType);
  }


  /**
   * Flags controlling the details of output in OpenDX format.
   *
   * @ingroup output
   */
  struct DXFlags : public OutputFlagsBase<DXFlags>
  {
    /**
     * Write neighbor information. This information is necessary for instance,
     * if OpenDX is supposed to compute integral curves (streamlines). If it
     * is not present, streamlines end at cell boundaries.
     */
    bool write_neighbors;
    /**
     * Write integer values of the Triangulation in binary format.
     */
    bool int_binary;
    /**
     * Write coordinate vectors in binary format.
     */
    bool coordinates_binary;

    /**
     * Write data vectors in binary format.
     */
    bool data_binary;

    /**
     * Write binary coordinate vectors as double (64 bit) numbers instead of
     * float (32 bit).
     */
    bool data_double;

    /**
     * Constructor.
     */
    DXFlags(const bool write_neighbors    = false,
            const bool int_binary         = false,
            const bool coordinates_binary = false,
            const bool data_binary        = false);

    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in declare_parameters() and set the flags
     * for this output format accordingly.
     *
     * The flags thus obtained overwrite all previous contents of this object.
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * Flags controlling the details of output in UCD format for AVS.
   *
   * @ingroup output
   */
  struct UcdFlags : public OutputFlagsBase<UcdFlags>
  {
    /**
     * Write a comment at the beginning of the file stating the date of
     * creation and some other data.  While this is supported by the UCD
     * format and AVS, some other programs get confused by this, so the
     * default is to not write a preamble. However, a preamble can be written
     * using this flag.
     *
     * Default: <code>false</code>.
     */
    bool write_preamble;

    /**
     * Constructor.
     */
    UcdFlags(const bool write_preamble = false);

    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in declare_parameters() and set the flags
     * for this output format accordingly.
     *
     * The flags thus obtained overwrite all previous contents of this object.
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * Flags controlling the details of output in Gnuplot format.
   *
   * @ingroup output
   */
  struct GnuplotFlags : public OutputFlagsBase<GnuplotFlags>
  {
    /**
     * Default constructor. Sets up the dimension labels with the default values
     * of <tt>"x"</tt>, <tt>"y"</tt>, and <tt>"z"</tt>.
     */
    GnuplotFlags();

    /**
     * Constructor which sets up non-default values for the dimension labels.
     */
    GnuplotFlags(const std::vector<std::string> &space_dimension_labels);

    /**
     * Labels to use in each spatial dimension. These default to <tt>"x"</tt>,
     * <tt>"y"</tt>, and <tt>"z"</tt>. Labels are printed to the Gnuplot file
     * surrounded by angle brackets: For example, if the space dimension is 2
     * and the labels are <tt>"x"</tt> and <tt>"t"</tt>, then the relevant
     * line will start with
     * @verbatim
     * # <x> <t>
     * @endverbatim
     * Any extra labels will be ignored.
     *
     * If you specify these labels yourself then there should be at least
     * <tt>spacedim</tt> labels, where <tt>spacedim</tt> is the spatial
     * dimension of the output data.
     */
    std::vector<std::string> space_dimension_labels;

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object.
     */
    std::size_t
    memory_consumption() const;

    /**
     * Exception to raise when there are not enough specified dimension
     * labels.
     */
    DeclExceptionMsg(ExcNotEnoughSpaceDimensionLabels,
                     "There should be at least one space dimension per spatial "
                     "dimension (extras are ignored).");
  };

  /**
   * Flags controlling the details of output in Povray format. Several flags
   * are implemented, see their respective documentation.
   *
   * @ingroup output
   */
  struct PovrayFlags : public OutputFlagsBase<PovrayFlags>
  {
    /**
     * Normal vector interpolation, if set to true
     *
     * default = false
     */
    bool smooth;

    /**
     * Use bicubic patches (b-splines) instead of triangles.
     *
     * default = false
     */
    bool bicubic_patch;

    /**
     * include external "data.inc" with camera, light and texture definition
     * for the scene.
     *
     * default = false
     */
    bool external_data;

    /**
     * Constructor.
     */
    PovrayFlags(const bool smooth        = false,
                const bool bicubic_patch = false,
                const bool external_data = false);

    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in declare_parameters() and set the flags
     * for this output format accordingly.
     *
     * The flags thus obtained overwrite all previous contents of this object.
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };


  /**
   * Flags controlling the details of output in encapsulated postscript
   * format.
   *
   * @ingroup output
   */
  struct EpsFlags : public OutputFlagsBase<EpsFlags>
  {
    /**
     * This denotes the number of the data vector which shall be used for
     * generating the height information. By default, the first data vector is
     * taken, i.e. <tt>height_vector==0</tt>, if there is any data vector. If
     * there is no data vector, no height information is generated.
     */
    unsigned int height_vector;

    /**
     * Number of the vector which is to be taken to colorize cells. The same
     * applies as for #height_vector.
     */
    unsigned int color_vector;

    /**
     * Enum denoting the possibilities whether the scaling should be done such
     * that the given <tt>size</tt> equals the width or the height of the
     * resulting picture.
     */
    enum SizeType
    {
      /// Scale to given width
      width,
      /// Scale to given height
      height
    };

    /**
     * See above. Default is <tt>width</tt>.
     */
    SizeType size_type;

    /**
     * Width or height of the output as given in postscript units This usually
     * is given by the strange unit 1/72 inch. Whether this is height or width
     * is specified by the flag <tt>size_type</tt>.
     *
     * Default is 300, which represents a size of roughly 10 cm.
     */
    unsigned int size;

    /**
     * Width of a line in postscript units. Default is 0.5.
     */
    double line_width;

    /**
     * Angle of the line origin-viewer against the z-axis in degrees.
     *
     * Default is the Gnuplot-default of 60.
     */
    double azimut_angle;

    /**
     * Angle by which the viewers position projected onto the x-y-plane is
     * rotated around the z-axis, in positive sense when viewed from above.
     * The unit are degrees, and zero equals a position above or below the
     * negative y-axis.
     *
     * Default is the Gnuplot-default of 30.  An example of a Gnuplot-default
     * of 0 is the following:
     *
     * @verbatim
     *
     *          3________7
     *          /       /|
     *         /       / |
     *       2/______6/  |
     *       |   |   |   |
     * O-->  |   0___|___4
     *       |  /    |  /
     *       | /     | /
     *      1|/______5/
     *
     * @endverbatim
     */
    double turn_angle;

    /**
     * Factor by which the z-axis is to be stretched as compared to the x- and
     * y-axes. This is to compensate for the different sizes that coordinate
     * and solution values may have and to prevent that the plot looks to much
     * out-of-place (no elevation at all if solution values are much smaller
     * than coordinate values, or the common "extremely mountainous area" in
     * the opposite case.
     *
     * Default is <tt>1.0</tt>.
     */
    double z_scaling;

    /**
     * Flag the determines whether the lines bounding the cells (or the parts
     * of each patch) are to be plotted.
     *
     * Default: <tt>true</tt>.
     */
    bool draw_mesh;

    /**
     * Flag whether to fill the regions between the lines bounding the cells
     * or not. If not, no hidden line removal is performed, which in this
     * crude implementation is done through writing the cells in a
     * back-to-front order, thereby hiding the cells in the background by
     * cells in the foreground.
     *
     * If this flag is <tt>false</tt> and #draw_mesh is <tt>false</tt> as
     * well, nothing will be printed.
     *
     * If this flag is <tt>true</tt>, then the cells will be drawn either
     * colored by one of the data sets (if #shade_cells is <tt>true</tt>), or
     * pure white (if #shade_cells is false or if there are no data sets).
     *
     * Default is <tt>true</tt>.
     */
    bool draw_cells;

    /**
     * Flag to determine whether the cells shall be colorized by the data set
     * denoted by #color_vector, or simply be painted in white. This flag only
     * makes sense if <tt>#draw_cells==true</tt>. Colorization is done through
     * #color_function.
     *
     * Default is <tt>true</tt>.
     */
    bool shade_cells;

    /**
     * Structure keeping the three color values in the RGB system.
     */
    struct RgbValues
    {
      float red;
      float green;
      float blue;

      /**
       * Return <tt>true</tt> if the color represented by the three color
       * values is a grey scale, i.e. all components are equal.
       */
      bool
      is_grey() const;
    };

    /**
     * Definition of a function pointer type taking a value and returning a
     * triple of color values in RGB values.
     *
     * Besides the actual value by which the color is to be computed, min and
     * max values of the data to be colorized are given as well.
     */
    using ColorFunction = RgbValues (*)(const double value,
                                        const double min_value,
                                        const double max_value);

    /**
     * This is a pointer to the function which is used to colorize the cells.
     * By default, it points to the static function default_color_function()
     * which is a member of this class.
     */
    ColorFunction color_function;


    /**
     * Default colorization function. This one does what one usually wants: It
     * shifts colors from black (lowest value) through blue, green and red to
     * white (highest value). For the exact definition of the color scale
     * refer to the implementation.
     *
     * This function was originally written by Stefan Nauber.
     */
    static RgbValues
    default_color_function(const double value,
                           const double min_value,
                           const double max_value);

    /**
     * This is an alternative color function producing a grey scale between
     * black (lowest values) and white (highest values). You may use it by
     * setting the #color_function variable to the address of this function.
     */
    static RgbValues
    grey_scale_color_function(const double value,
                              const double min_value,
                              const double max_value);

    /**
     * This is one more alternative color function producing a grey scale
     * between white (lowest values) and black (highest values), i.e. the
     * scale is reversed to the previous one. You may use it by setting the
     * #color_function variable to the address of this function.
     */
    static RgbValues
    reverse_grey_scale_color_function(const double value,
                                      const double min_value,
                                      const double max_value);

    /**
     * Constructor.
     */
    EpsFlags(const unsigned int  height_vector  = 0,
             const unsigned int  color_vector   = 0,
             const SizeType      size_type      = width,
             const unsigned int  size           = 300,
             const double        line_width     = 0.5,
             const double        azimut_angle   = 60,
             const double        turn_angle     = 30,
             const double        z_scaling      = 1.0,
             const bool          draw_mesh      = true,
             const bool          draw_cells     = true,
             const bool          shade_cells    = true,
             const ColorFunction color_function = &default_color_function);

    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     *
     * For coloring, only the color functions declared in this class are
     * offered.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in declare_parameters() and set the flags
     * for this output format accordingly.
     *
     * The flags thus obtained overwrite all previous contents of this object.
     */
    void
    parse_parameters(const ParameterHandler &prm);
  };

  /**
   * Flags controlling the details of output in GMV format. At present no
   * flags are implemented.
   *
   * @ingroup output
   */
  struct GmvFlags : public OutputFlagsBase<GmvFlags>
  {};

  /**
   * Flags controlling the details of output in HDF5 format.
   *
   * @ingroup output
   */
  struct Hdf5Flags : public OutputFlagsBase<Hdf5Flags>
  {
    /**
     * Flag determining the compression level at which zlib, if available, is
     * run. The default is <tt>best_speed</tt>.
     */
    DataOutBase::CompressionLevel compression_level;

    explicit Hdf5Flags(
      const CompressionLevel compression_level = CompressionLevel::best_speed);
  };

  /**
   * Flags controlling the details of output in Tecplot format.
   *
   * @ingroup output
   */
  struct TecplotFlags : public OutputFlagsBase<TecplotFlags>
  {
    /**
     * Tecplot allows to assign names to zones. This variable stores this
     * name.
     */
    const char *zone_name;

    /**
     * Solution time for each zone in a strand. This value must be
     * non-negative, otherwise it will not be written to file. Do not assign any
     * value for this in case of a static zone.
     */
    double solution_time;

    /**
     * Constructor.
     */
    TecplotFlags(const char  *zone_name     = nullptr,
                 const double solution_time = -1.0);

    /**
     * Return an estimate for the memory consumption, in bytes, of this
     * object.
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * Flags controlling the details of output in VTK format.
   *
   * @ingroup output
   */
  struct VtkFlags : public OutputFlagsBase<VtkFlags>
  {
    /**
     * The time of the time step if this file is part of a time dependent
     * simulation.
     *
     * The value of this variable is written into the output file according to
     * the instructions provided in
     * http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
     * unless it is at its default value of
     * @verbatim std::numeric_limits<unsigned int>::min() @endverbatim.
     */
    double time;

    /**
     * The number of the time step if this file is part of a time dependent
     * simulation, or the cycle within a nonlinear or other iteration.
     *
     * The value of this variable is written into the output file according to
     * the instructions provided in
     * http://www.visitusers.org/index.php?title=Time_and_Cycle_in_VTK_files
     * unless it is at its default value of
     * @verbatim std::numeric_limits<unsigned int>::min() @endverbatim.
     */
    unsigned int cycle;

    /**
     * Flag to determine whether the current date and time shall be printed as
     * a comment in the file's second line.
     *
     * Default is <tt>true</tt>.
     */
    bool print_date_and_time;

    /**
     * Flag determining the compression level at which zlib, if available, is
     * run. The default is <tt>best_speed</tt>.
     */
    DataOutBase::CompressionLevel compression_level;

    /**
     * Flag determining whether to write patches as linear cells
     * or as a high-order Lagrange cell.
     *
     * Default is <tt>false</tt>.
     *
     * @note The ability to write data that corresponds to higher order
     * polynomials rather than simply linear or bilinear is a feature that was
     * only introduced in VTK 8.1.0 in December 2017. You will need at least
     * Paraview version 5.5.0 released in April 2018 or a similarly recent
     * version of VisIt for this feature to work (for example, VisIt 3.1.1,
     * released in February 2020, does not yet support this feature). Older
     * versions of these programs are likely going to result in errors when
     * trying to read files generated with this flag set to true. Experience
     * with these programs shows that these error messages are likely going to
     * be rather less descriptive and more obscure.
     */
    bool write_higher_order_cells;

    /**
     * A map that describes for (some or all) of the output quantities what
     * the physical units are. This field is ignored for VTK file format, but
     * used for VTU format where it is attached to the individual scalar,
     * vector, or tensor fields for use by visualization or other postprocessing
     * tools. The default is to not attach any physical units to fields at all,
     * i.e., an empty map.
     *
     * If the map does not contain an entry for a specific output variable, then
     * no unit will be written into the output file. In other words, it is not
     * an error to provide units for only some variables.
     *
     * step-19, step-44 and step-69 all demonstrate how to use this variable.
     *
     * @note While the functions that make use of this information do not care
     *   about how physical units are specified, downstream postprocessing tools
     *   should and do. As a consequence, these units should be specified in a
     *   format that is understandable to these postprocessing tools. As an
     *   example, the [unyt project](https://unyt.readthedocs.io/en/stable/)
     *   describes a standard for describing and converting units.
     */
    std::map<std::string, std::string> physical_units;

    /**
     * Constructor. Initializes the member variables with names corresponding
     * to the argument names of this function.
     */
    explicit VtkFlags(
      const double           time  = std::numeric_limits<double>::min(),
      const unsigned int     cycle = std::numeric_limits<unsigned int>::min(),
      const bool             print_date_and_time = true,
      const CompressionLevel compression_level   = CompressionLevel::best_speed,
      const bool             write_higher_order_cells          = false,
      const std::map<std::string, std::string> &physical_units = {});
  };


  /**
   * Flags for SVG output.
   *
   * @ingroup output
   */
  struct SvgFlags : public OutputFlagsBase<SvgFlags>
  {
    /**
     * Height of the image in SVG units. Default value is 4000.
     */
    unsigned int height;

    /**
     * Width of the image in SVG units. If left zero, the width is computed
     * from the height.
     */
    unsigned int width;

    /**
     * This denotes the number of the data vector which shall be used for
     * generating the height information. By default, the first data vector is
     * taken, i.e. <tt>#height_vector==0</tt>, if there is any data vector. If
     * there is no data vector, no height information is generated.
     */
    unsigned int height_vector;

    /**
     * Angles for the perspective view
     */
    int azimuth_angle, polar_angle;

    unsigned int line_thickness;

    /**
     * Draw a margin of 5% around the plotted area
     */
    bool margin;

    /**
     * Draw a colorbar encoding the cell coloring
     */
    bool draw_colorbar;

    /**
     * Constructor.
     */
    SvgFlags(const unsigned int height_vector  = 0,
             const int          azimuth_angle  = 37,
             const int          polar_angle    = 45,
             const unsigned int line_thickness = 1,
             const bool         margin         = true,
             const bool         draw_colorbar  = true);
  };


  /**
   * Flags controlling the details of output in deal.II intermediate format.
   * At present no flags are implemented.
   *
   * @ingroup output
   */
  struct Deal_II_IntermediateFlags
    : public OutputFlagsBase<Deal_II_IntermediateFlags>
  {
    /**
     * An indicator of the current file format version used to write
     * intermediate format. We do not attempt to be backward compatible, so
     * this number is used only to verify that the format we are writing is
     * what the current readers and writers understand.
     */
    static const unsigned int format_version;
  };

  /**
   * Flags controlling the behavior of the DataOutFilter class.
   *
   * @ingroup output
   */
  struct DataOutFilterFlags
  {
    /**
     * Whether or not to filter out duplicate vertices and associated values.
     * Setting this value to `true` will drastically
     * reduce the output data size but will result in an output file that
     * does not faithfully represent the actual data if the data corresponds
     * to discontinuous fields. In particular, along subdomain boundaries
     * the data will still be discontinuous, while it will look like a
     * continuous field inside of the subdomain.
     */
    bool filter_duplicate_vertices;

    /**
     * Whether the XDMF output refers to HDF5 files. This affects how output
     * is structured.
     */
    bool xdmf_hdf5_output;

    /**
     * Constructor.
     */
    DataOutFilterFlags(const bool filter_duplicate_vertices = false,
                       const bool xdmf_hdf5_output          = false);

    /**
     * Declare all flags with name and type as offered by this class, for use
     * in input files.
     */
    static void
    declare_parameters(ParameterHandler &prm);

    /**
     * Read the parameters declared in <tt>declare_parameters</tt> and set the
     * flags for this output format accordingly.
     *
     * The flags thus obtained overwrite all previous contents of this object.
     */
    void
    parse_parameters(const ParameterHandler &prm);

    /**
     * Determine an estimate for the memory consumption (in bytes) of this
     * object.
     */
    std::size_t
    memory_consumption() const;
  };

  /**
   * DataOutFilter provides a way to remove redundant vertices and values
   * generated by the deal.II output. By default, DataOutBase and the classes
   * that build on it output data at each vertex of each cell. This means that
   * data is output multiple times for each vertex of the mesh, once for each
   * cell adjacent to the vertex. The purpose of
   * this scheme is to support output of discontinuous quantities, either
   * because the finite element space is discontinuous or because the quantity
   * that is output is computed from a solution field and is discontinuous
   * across faces (for example for quantities computed via DataPostprocessor;
   * typical cases where output quantities are discontinuous are when a
   * postprocessor computes a quantity using the *gradient* of the solution,
   * which is generally discontinuous even if the element itself is
   * continuous). Other cases where the output is discontinuous are if the
   * data to be output is not a finite element field but, for example,
   * results from a class such as MatrixOut.
   *
   * This class is an attempt to rein in the amount of data that is written.
   * If the fields that are written to files are indeed discontinuous, the
   * only way to faithfully represent them is indeed to write multiple values
   * for each vertex (this is typically done by creating multiple logical nodes
   * in the output file, all of which have the same physical location; data is
   * then associated to nodes, allowing to have multiple values associated
   * with the same location). However,
   * for fine meshes, one may not necessarily be interested in an exact
   * representation of output fields that will likely only have small
   * discontinuities. Rather, it may be sufficient to just output one value
   * per vertex, which may be chosen arbitrarily from among those that are
   * defined at this vertex, i.e., chosen arbitrarily from any of the
   * adjacent cells.
   */
  class DataOutFilter
  {
  public:
    /**
     * Default constructor.
     */
    DataOutFilter();

    /**
     * Constructor with a given set of flags. See DataOutFilterFlags for
     * possible flags.
     */
    DataOutFilter(const DataOutBase::DataOutFilterFlags &flags);

    /**
     * Write a point with the specified index into the filtered data set. If
     * the point already exists and we are filtering redundant values, the
     * provided index will internally refer to another recorded point.
     */
    template <int dim>
    void
    write_point(const unsigned int index, const Point<dim> &p);

    /**
     * Record a deal.II cell in the internal reordered format.
     */
    template <int dim>
    void
    write_cell(const unsigned int                   index,
               const unsigned int                   start,
               const std::array<unsigned int, dim> &offsets);

    /**
     * Record a single deal.II cell without subdivisions (e.g. simplex) in the
     * internal reordered format.
     */
    void
    write_cell_single(const unsigned int   index,
                      const unsigned int   start,
                      const unsigned int   n_points,
                      const ReferenceCell &reference_cell);

    /**
     * Filter and record a data set. If there are multiple values at a given
     * vertex and redundant values are being removed, one is arbitrarily
     * chosen as the recorded value. In the future this can be expanded to
     * average/min/max multiple values at a given vertex.
     */
    void
    write_data_set(const std::string      &name,
                   const unsigned int      dimension,
                   const unsigned int      set_num,
                   const Table<2, double> &data_vectors);

    /**
     * Resize and fill a vector with all the filtered node vertex points, for
     * output to a file.
     */
    void
    fill_node_data(std::vector<double> &node_data) const;

    /**
     * Resize and fill a vector with all the filtered cell vertex indices, for
     * output to a file.
     */
    void
    fill_cell_data(const unsigned int         local_node_offset,
                   std::vector<unsigned int> &cell_data) const;

    /**
     * Get the name of the data set indicated by the set number.
     */
    std::string
    get_data_set_name(const unsigned int set_num) const;

    /**
     * Get the dimensionality of the data set indicated by the set number.
     */
    unsigned int
    get_data_set_dim(const unsigned int set_num) const;

    /**
     * Get the raw double valued data of the data set indicated by the set
     * number.
     */
    const double *
    get_data_set(const unsigned int set_num) const;

    /**
     * Return the number of nodes in this DataOutFilter. This may be smaller
     * than the original number of nodes if filtering is enabled.
     */
    unsigned int
    n_nodes() const;

    /**
     * Return the number of filtered cells in this DataOutFilter. Cells are
     * not filtered so this will be the original number of cells.
     */
    unsigned int
    n_cells() const;

    /**
     * Return the number of filtered data sets in this DataOutFilter. Data
     * sets are not filtered so this will be the original number of data sets.
     */
    unsigned int
    n_data_sets() const;

    /**
     * Empty functions to do base class inheritance.
     */
    void
    flush_points();

    /**
     * Empty functions to do base class inheritance.
     */
    void
    flush_cells();


  private:
    /**
     * Empty class to provide comparison function for Map3dPoint.
     */
    struct Point3Comp
    {
      bool
      operator()(const Point<3> &one, const Point<3> &two) const
      {
        /*
         * The return statement below is an optimized version of the following
         * code:
         *
         * for (unsigned int d=0; d<3; ++d)
         * {
         *   if (one(d) < two(d))
         *     return true;
         *   else if (one(d) > two(d))
         *     return false;
         * }
         * return false;
         */

        return (one[0] < two[0] ||
                (!(two[0] < one[0]) &&
                 (one[1] < two[1] || (!(two[1] < one[1]) && one[2] < two[2]))));
      }
    };

    using Map3DPoint = std::multimap<Point<3>, unsigned int, Point3Comp>;

    /**
     * Flags used to specify filtering behavior.
     */
    DataOutBase::DataOutFilterFlags flags;

    /**
     * The number of space dimensions in which the vertices represented
     * by the current object live. This corresponds to the usual
     * @p dim argument, but since this class is not templated on the
     * dimension, we need to store it here.
     */
    unsigned int node_dim;

    /**
     * The number of cells stored in
     * @ref filtered_cells.
     */
    unsigned int num_cells;

    /**
     * Map of points to an internal index.
     */
    Map3DPoint existing_points;

    /**
     * Map of actual point index to internal point index.
     */
    std::map<unsigned int, unsigned int> filtered_points;

    /**
     * Map of cells to the filtered points.
     */
    std::map<unsigned int, unsigned int> filtered_cells;

    /**
     * Data set names.
     */
    std::vector<std::string> data_set_names;

    /**
     * Data set dimensions.
     */
    std::vector<unsigned int> data_set_dims;

    /**
     * Data set data.
     */
    std::vector<std::vector<double>> data_sets;

    /**
     * Record a cell vertex index based on the internal reordering.
     */
    void
    internal_add_cell(const unsigned int cell_index,
                      const unsigned int pt_index);
  };


  /**
   * Provide a data type specifying the presently supported output formats.
   */
  enum OutputFormat
  {
    /**
     * Use the format already stored in the object.
     */
    default_format,

    /**
     * Do not write any output.
     */
    none,

    /**
     * Output for OpenDX.
     */
    dx,

    /**
     * Output in the UCD format for AVS.
     */
    ucd,

    /**
     * Output for the Gnuplot tool.
     */
    gnuplot,

    /**
     * Output for the Povray raytracer.
     */
    povray,

    /**
     * Output in encapsulated PostScript.
     */
    eps,

    /**
     * Output for GMV.
     */
    gmv,

    /**
     * Output for Tecplot in text format.
     */
    tecplot,

    /**
     * Output in VTK format.
     */
    vtk,

    /**
     * Output in VTK format.
     */
    vtu,

    /**
     * Output in SVG format.
     */
    svg,

    /**
     * Output in deal.II intermediate format.
     */
    deal_II_intermediate,

    /**
     * Output in HDF5 format.
     */
    hdf5
  };


  /**
   * Write the given list of patches to the output stream in OpenDX format.
   */
  template <int dim, int spacedim>
  void
  write_dx(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                  &nonscalar_data_ranges,
    const DXFlags &flags,
    std::ostream  &out);

  /**
   * Write the given list of patches to the output stream in eps format.
   *
   * Output in this format circumvents the use of auxiliary graphic programs
   * converting some output format into a graphics format. This has the
   * advantage that output is easy and fast, and the disadvantage that you
   * have to give a whole bunch of parameters which determine the direction of
   * sight, the mode of colorization, the scaling of the height axis, etc. (Of
   * course, all these parameters have reasonable default values, which you
   * may want to change.)
   *
   * This function only supports output for two-dimensional domains (i.e.,
   * with dim==2), with values in the vertical direction taken from a data
   * vector.
   *
   * Basically, output consists of the mesh and the cells in between them. You
   * can draw either of these, or both, or none if you are really interested
   * in an empty picture. If written, the mesh uses black lines. The cells in
   * between the mesh are either not printed (this will result in a loss of
   * hidden line removal, i.e.  you can "see through" the cells to lines
   * behind), printed in white (which does nothing apart from the hidden line
   * removal), or colorized using one of the data vectors (which need not be
   * the same as the one used for computing the height information) and a
   * customizable color function. The default color functions chooses the
   * color between black, blue, green, red and white, with growing values of
   * the data field chosen for colorization. At present, cells are displayed
   * with one color per cell only, which is taken from the value of the data
   * field at the center of the cell; bilinear interpolation of the color on a
   * cell is not used.
   *
   * By default, the viewpoint is chosen like the default viewpoint in
   * GNUPLOT, i.e.  with an angle of 60 degrees with respect to the positive
   * z-axis and rotated 30 degrees in positive sense (as seen from above) away
   * from the negative y-axis.  Of course you can change these settings.
   *
   * EPS output is written without a border around the picture, i.e. the
   * bounding box is close to the output on all four sides. Coordinates are
   * written using at most five digits, to keep picture size at a reasonable
   * size.
   *
   * All parameters along with their default values are listed in the
   * documentation of the <tt>EpsFlags</tt> member class of this class. See
   * there for more and detailed information.
   */
  template <int spacedim>
  void
  write_eps(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string>        &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const EpsFlags &flags,
    std::ostream   &out);

  /**
   * This is the same function as above except for domains that are not
   * two-dimensional. This function is not implemented (and will throw an error
   * if called) but is declared to allow for dimension-independent programs.
   */
  template <int dim, int spacedim>
  void
  write_eps(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const EpsFlags &flags,
    std::ostream   &out);


  /**
   * Write the given list of patches to the output stream in GMV format.
   *
   * Data is written in the following format: nodes are considered the points
   * of the patches. In spatial dimensions less than three, zeroes are
   * inserted for the missing coordinates. The data vectors are written as
   * node or cell data, where for the first the data space is interpolated to
   * (bi-,tri-)linear elements.
   */
  template <int dim, int spacedim>
  void
  write_gmv(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const GmvFlags &flags,
    std::ostream   &out);

  /**
   * Write the given list of patches to the output stream in gnuplot format.
   * Visualization of two-dimensional data can then be achieved by starting
   * <tt>gnuplot</tt> and entering the commands
   *
   * @verbatim
   * set data style lines
   * splot "filename" using 1:2:n
   * @endverbatim
   * This example assumes that the number of the data vector displayed is
   * <b>n-2</b>.
   *
   * The GNUPLOT format is not able to handle data on unstructured grids
   * directly. Directly would mean that you only give the vertices and the
   * solution values thereon and the program constructs its own grid to
   * represent the data. This is only possible for a structured tensor product
   * grid in two dimensions. However, it is possible to give several such
   * patches within one file, which is exactly what the respective function of
   * this class does: writing each cell's data as a patch of data, at least if
   * the patches as passed from derived classes represent cells. Note that the
   * functions on patches need not be continuous at interfaces between
   * patches, so this method also works for discontinuous elements. Note also,
   * that GNUPLOT can do hidden line removal for patched data.
   *
   * While this discussion applies to two spatial dimensions, it is more
   * complicated in 3d. The reason is that we could still use patches, but it
   * is difficult when trying to visualize them, since if we use a cut through
   * the data (by, for example, using x- and z-coordinates, a fixed y-value
   * and plot function values in z-direction, then the patched data is not a
   * patch in the sense GNUPLOT wants it any more. Therefore, we use another
   * approach, namely writing the data on the 3d grid as a sequence of lines,
   * i.e. two points each associated with one or more data sets.  There are
   * therefore 12 lines for each subcells of a patch.
   *
   * Given the lines as described above, a cut through this data in Gnuplot
   * can then be achieved like this:
   * @verbatim
   *   set data style lines
   *   splot [:][:][0:] "T" using 1:2:(\$3==.5 ? \$4 : -1)
   * @endverbatim
   *
   * This command plots data in $x$- and $y$-direction unbounded, but in
   * $z$-direction only those data points which are above the $x$-$y$-plane (we
   * assume here a positive solution, if it has negative values, you might
   * want to decrease the lower bound). Furthermore, it only takes the data
   * points with z-values (<tt>&3</tt>) equal to 0.5, i.e. a cut through the
   * domain at <tt>z=0.5</tt>. For the data points on this plane, the data
   * values of the first data set (<tt>&4</tt>) are raised in z-direction
   * above the x-y-plane; all other points are denoted the value <tt>-1</tt>
   * instead of the value of the data vector and are not plotted due to the
   * lower bound in z plotting direction, given in the third pair of brackets.
   *
   * More complex cuts are possible, including nonlinear ones. Note however,
   * that only those points which are actually on the cut-surface are plotted.
   */
  template <int dim, int spacedim>
  void
  write_gnuplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                       &nonscalar_data_ranges,
    const GnuplotFlags &flags,
    std::ostream       &out);

  /**
   * Write the given list of patches to the output stream for the Povray
   * raytracer.
   *
   * Output in this format creates a povray source file, include standard
   * camera and light source definition for rendering with povray 3.1 At
   * present, this format only supports output for two-dimensional data, with
   * values in the third direction taken from a data vector.
   *
   * The output uses two different povray-objects:
   *
   * <ul>
   * <li> <tt>BICUBIC_PATCH</tt> A <tt>bicubic_patch</tt> is a 3-dimensional
   * Bezier patch. It consists of 16 Points describing the surface. The 4
   * corner points are touched by the object, while the other 12 points pull
   * and stretch the patch into shape. One <tt>bicubic_patch</tt> is generated
   * on each patch. Therefore the number of subdivisions has to be 3 to provide
   * the patch with 16 points. A bicubic patch is not exact but generates very
   * smooth images.
   *
   * <li> <tt>MESH</tt> The mesh object is used to store large number of
   * triangles. Every square of the patch data is split into one upper-left
   * and one lower-right triangle. If the number of subdivisions is three, 32
   * triangle are generated for every patch.
   *
   * Using the smooth flag povray interpolates the normals on the triangles,
   * imitating a curved surface
   * </ul>
   *
   * All objects get one texture definition called Tex. This texture has to be
   * declared somewhere before the object data. This may be in an external
   * data file or at the beginning of the output file. Setting the
   * <tt>external_data</tt> flag to false, an standard camera, light and
   * texture (scaled to fit the scene) is added to the output file. Set to
   * true an include file "data.inc" is included. This file is not generated
   * by deal and has to include camera, light and the texture definition Tex.
   *
   * You need povray (>=3.0) to render the scene. The minimum options for
   * povray are:
   * @verbatim
   *   povray +I<inputfile> +W<horiz. size> +H<ver. size> +L<include path>
   * @endverbatim
   * If the external file "data.inc" is used, the path to this file has to be
   * included in the povray options.
   */
  template <int dim, int spacedim>
  void
  write_povray(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                      &nonscalar_data_ranges,
    const PovrayFlags &flags,
    std::ostream      &out);

  /**
   * Write the given list of patches to the output stream in Tecplot ASCII
   * format (FEBLOCK).
   *
   * For more information consult the Tecplot Users and Reference manuals.
   */
  template <int dim, int spacedim>
  void
  write_tecplot(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                       &nonscalar_data_ranges,
    const TecplotFlags &flags,
    std::ostream       &out);

  /**
   * Write the given list of patches to the output stream in UCD format
   * described in the AVS developer's guide (now AVS). Due to limitations in
   * the present format, only node based data can be output, which in one
   * reason why we invented the patch concept. In order to write higher order
   * elements, you may split them up into several subdivisions of each cell.
   * These subcells will then, however, also appear as different cells by
   * programs which understand the UCD format.
   *
   * No use is made of the possibility to give model data since these are not
   * supported by all UCD aware programs. You may give cell data in derived
   * classes by setting all values of a given data set on a patch to the same
   * value.
   */
  template <int dim, int spacedim>
  void
  write_ucd(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const UcdFlags &flags,
    std::ostream   &out);

  /**
   * Write the given list of patches to the output stream in VTK format. The
   * data is written in the traditional VTK format as opposed to the XML-based
   * format that write_vtu() produces.
   *
   * The nonscalar_data_ranges argument denotes ranges of components in the
   * output that are considered a vector, rather than simply a collection of
   * scalar fields. The VTK output format has special provisions that allow
   * these components to be output by a single name rather than having to
   * group several scalar fields into a vector later on in the visualization
   * program.
   *
   * @note VTK is a legacy format and has largely been supplanted by the VTU
   * format (an XML-structured version of VTK). In particular, VTU allows for
   * the compression of data and consequently leads to much smaller file sizes
   * that equivalent VTK files for large files. Since all visualization
   * programs that support VTK also support VTU, you should consider using the
   * latter file format instead, by using the write_vtu() function.
   */
  template <int dim, int spacedim>
  void
  write_vtk(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out);


  /**
   * Write the given list of patches to the output stream in VTU format. The
   * data is written in the XML-based VTK format as opposed to the traditional
   * format that write_vtk() produces.
   *
   * The nonscalar_data_ranges argument denotes ranges of components in the
   * output that are considered a vector, rather than simply a collection of
   * scalar fields. The VTK output format has special provisions that allow
   * these components to be output by a single name rather than having to
   * group several scalar fields into a vector later on in the visualization
   * program.
   *
   * Some visualization programs, such as ParaView, can read several separate
   * VTU files to parallelize visualization. In that case, you need a
   * <code>.pvtu</code> file that describes which VTU files form a group. The
   * DataOutInterface::write_pvtu_record() function can generate such a
   * centralized record. Likewise, DataOutInterface::write_visit_record() does
   * the same for VisIt (although VisIt can also read <code>pvtu</code> records
   * since version 2.5.1). Finally, for time dependent problems, you may also
   * want to look at DataOutInterface::write_pvd_record()
   *
   * The use of this function is explained in step-40.
   */
  template <int dim, int spacedim>
  void
  write_vtu(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out);

  /**
   * This writes the header for the xml based vtu file format. This routine is
   * used internally together with DataOutInterface::write_vtu_footer() and
   * DataOutInterface::write_vtu_main() by DataOutBase::write_vtu().
   */
  void
  write_vtu_header(std::ostream &out, const VtkFlags &flags);

  /**
   * This function writes the footer for the xml based vtu file format. This
   * routine is used internally together with
   * DataOutInterface::write_vtu_header() and DataOutInterface::write_vtu_main()
   * by DataOutBase::write_vtu().
   */
  void
  write_vtu_footer(std::ostream &out);

  /**
   * This function writes the main part for the xml based vtu file format. This
   * routine is used internally together with
   * DataOutInterface::write_vtu_header() and
   * DataOutInterface::write_vtu_footer() by DataOutBase::write_vtu().
   */
  template <int dim, int spacedim>
  void
  write_vtu_main(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags,
    std::ostream   &out);

  /**
   * Some visualization programs, such as ParaView, can read several separate
   * VTU files that all form part of the same simulation, in order to
   * parallelize visualization. In that case, you need a
   * <code>.pvtu</code> file that describes which VTU files (written, for
   * example, through the DataOutInterface::write_vtu() function) form a group.
   * The current function can generate such a centralized record.
   *
   * This function is typically not called by itself from user space, but
   * you may want to call it through DataOutInterface::write_pvtu_record()
   * since the DataOutInterface class has access to information that you
   * would have to provide to the current function by hand.
   *
   * In any case, whether this function is called directly or via
   * DataOutInterface::write_pvtu_record(), the central record file so
   * written contains a list of (scalar or vector) fields that describes which
   * fields can actually be found in the individual files that comprise the set
   * of parallel VTU files along with the names of these files. This function
   * gets the names and types of fields through the third and fourth
   * argument; you can determine these by hand, but in practice, this function
   * is most easily called by calling DataOutInterfaces::write_pvtu_record(),
   * which determines the last two arguments by calling
   * DataOutInterface::get_dataset_names() and
   * DataOutInterface::get_nonscalar_data_ranges() functions. The second
   * argument to this function specifies the names of the files that form the
   * parallel set.
   *
   * @note Use DataOutBase::write_vtu() and DataOutInterface::write_vtu()
   * for writing each piece. Also note that
   * only one parallel process needs to call the current function, listing the
   * names of the files written by all parallel processes.
   *
   * @note In order to tell Paraview to group together multiple
   * <code>pvtu</code> files that each describe one time step of a time
   * dependent simulation, see the DataOutBase::write_pvd_record()
   * function.
   *
   * @note Older versions of VisIt (before 2.5.1), can not read
   * <code>pvtu</code> records. However, it can read visit records as written
   * by the write_visit_record() function.
   */
  void
  write_pvtu_record(
    std::ostream                   &out,
    const std::vector<std::string> &piece_names,
    const std::vector<std::string> &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const VtkFlags &flags);

  /**
   * In ParaView it is possible to visualize time-dependent data tagged with
   * the current integration time of a time dependent simulation. To use this
   * feature you need a <code>.pvd</code> file that describes which VTU or
   * PVTU file belongs to which timestep. This function writes a file that
   * provides this mapping, i.e., it takes a list of pairs each of which
   * indicates a particular time instant and the corresponding file that
   * contains the graphical data for this time instant.
   *
   * A typical use case, in program that computes a time dependent solution,
   * would be the following (<code>time</code> and <code>time_step</code> are
   * member variables of the class with types <code>double</code> and
   * <code>unsigned int</code>, respectively; the variable
   * <code>times_and_names</code> is of type
   * <code>std::vector@<std::pair@<double,std::string@> @></code>):
   *
   * @code
   * template <int dim>
   * void MyEquation<dim>::output_results () const
   * {
   *   DataOut<dim> data_out;
   *
   *   data_out.attach_dof_handler(dof_handler);
   *   data_out.add_data_vector(solution, "U");
   *   data_out.build_patches();
   *
   *   const std::string filename = "solution-" +
   *                                Utilities::int_to_string (timestep_n, 3) +
   *                                ".vtu";
   *   std::ofstream output(filename);
   *   data_out.write_vtu(output);
   *
   *   times_and_names.emplace_back (time, filename);
   *   std::ofstream pvd_output ("solution.pvd");
   *   DataOutBase::write_pvd_record (pvd_output, times_and_names);
   * }
   * @endcode
   *
   * @note See DataOutInterface::write_vtu, DataOutInterface::write_pvtu_record,
   * and DataOutInterface::write_vtu_in_parallel
   * for writing solutions at each timestep.
   *
   * @note The second element of each pair, i.e., the file in which the
   * graphical data for each time is stored, may itself be again a file that
   * references other files. For example, it could be the name for a
   * <code>.pvtu</code> file that references multiple parts of a parallel
   * computation.
   */
  void
  write_pvd_record(
    std::ostream                                      &out,
    const std::vector<std::pair<double, std::string>> &times_and_names);

  /**
   * This function is the exact equivalent of the write_pvtu_record() function
   * but for older versions of the VisIt visualization program and for one
   * visualization graph (or one time step only). See there for the purpose of
   * this function.
   *
   * This function is documented in the "Creating a master file for parallel"
   * section (section 5.7) of the "Getting data into VisIt" report that can be
   * found here:
   * https://visit-dav.github.io/visit-website/pdfs/GettingDataIntoVisIt2.0.0.pdf
   */
  void
  write_visit_record(std::ostream                   &out,
                     const std::vector<std::string> &piece_names);

  /**
   * This function is equivalent to the write_visit_record() above but for
   * multiple time steps. Here is an example of how the function would be
   * used:
   * @code
   * const unsigned int number_of_time_steps = 3;
   * std::vector<std::vector<std::string > > piece_names(number_of_time_steps);
   *
   * piece_names[0].emplace_back("subdomain_01.time_step_0.vtk");
   * piece_names[0].emplace_back("subdomain_02.time_step_0.vtk");
   *
   * piece_names[1].emplace_back("subdomain_01.time_step_1.vtk");
   * piece_names[1].emplace_back("subdomain_02.time_step_1.vtk");
   *
   * piece_names[2].emplace_back("subdomain_01.time_step_2.vtk");
   * piece_names[2].emplace_back("subdomain_02.time_step_2.vtk");
   *
   * std::ofstream visit_output ("solution.visit");
   *
   * DataOutBase::write_visit_record(visit_output, piece_names);
   * @endcode
   *
   * This function is documented in the "Creating a master file for parallel"
   * section (section 5.7) of the "Getting data into VisIt" report that can be
   * found here:
   * https://visit-dav.github.io/visit-website/pdfs/GettingDataIntoVisIt2.0.0.pdf
   */
  void
  write_visit_record(std::ostream                                &out,
                     const std::vector<std::vector<std::string>> &piece_names);

  /**
   * This function is equivalent to the write_visit_record() above but for
   * multiple time steps and with additional information about the time for
   * each timestep. Here is an example of how the function would be
   * used:
   * @code
   * const unsigned int number_of_time_steps = 3;
   * std::vector<std::pair<double,std::vector<std::string > > >
   * times_and_piece_names(number_of_time_steps);
   *
   * times_and_piece_names[0].first = 0.0;
   * times_and_piece_names[0].second.emplace_back("subdomain_01.time_step_0.vtk");
   * times_and_piece_names[0].second.emplace_back("subdomain_02.time_step_0.vtk");
   *
   * times_and_piece_names[1].first = 0.5;
   * times_and_piece_names[1].second.emplace_back("subdomain_01.time_step_1.vtk");
   * times_and_piece_names[1].second.emplace_back("subdomain_02.time_step_1.vtk");
   *
   * times_and_piece_names[2].first = 1.0;
   * times_and_piece_names[2].second.emplace_back("subdomain_01.time_step_2.vtk");
   * times_and_piece_names[2].second.emplace_back("subdomain_02.time_step_2.vtk");
   *
   * std::ofstream visit_output ("solution.visit");
   *
   * DataOutBase::write_visit_record(visit_output, times_and_piece_names);
   * @endcode
   *
   * This function is documented in the "Creating a master file for parallel"
   * section (section 5.7) of the "Getting data into VisIt" report that can be
   * found here:
   * https://visit-dav.github.io/visit-website/pdfs/GettingDataIntoVisIt2.0.0.pdf
   */
  void
  write_visit_record(
    std::ostream &out,
    const std::vector<std::pair<double, std::vector<std::string>>>
      &times_and_piece_names);

  /**
   * Write the given list of patches to the output stream in SVG format.
   *
   * SVG (Scalable Vector Graphics) is an XML-based vector image format
   * developed and maintained by the World Wide Web Consortium (W3C). This
   * function conforms to the latest specification SVG 1.1, released on August
   * 16, 2011. Controlling the graphic output is possible by setting or
   * clearing the respective flags (see the SvgFlags struct). At present, this
   * format only supports output for two-dimensional data, with values in the
   * third direction taken from a data vector.
   *
   * For the output, each patch is subdivided into four triangles which are
   * then written as polygons and filled with a linear color gradient. The
   * arising coloring of the patches visualizes the data values at the
   * vertices taken from the specified data vector. A colorbar can be drawn to
   * encode the coloring.
   *
   * @note This function is so far only implemented for two dimensions with an
   * additional dimension reserved for data information.
   */
  template <int spacedim>
  void
  write_svg(
    const std::vector<Patch<2, spacedim>> &patches,
    const std::vector<std::string>        &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                   &nonscalar_data_ranges,
    const SvgFlags &flags,
    std::ostream   &out);

  /**
   * Write the given list of patches to the output stream in deal.II
   * intermediate format. This is not a format understood by any other
   * graphics program, but is rather a direct dump of the intermediate
   * internal format used by deal.II. This internal format is generated by the
   * various classes that can generate output using the DataOutBase class, for
   * example from a finite element solution, and is then converted in the
   * present class to the final graphics format.
   *
   * Note that the intermediate format is what its name suggests: a direct
   * representation of internal data. It isn't standardized and will change
   * whenever we change our internal representation. You can only expect to
   * process files written in this format using the same version of deal.II
   * that was used for writing.
   *
   * The reason why we offer to write out this intermediate format is that it
   * can be read back into a deal.II program using the DataOutReader class,
   * which is helpful in at least two contexts: First, this can be used to
   * later generate graphical output in any other graphics format presently
   * understood; this way, it is not necessary to know at run-time which
   * output format is requested, or if multiple output files in different
   * formats are needed. Secondly, in contrast to almost all other graphics
   * formats, it is possible to merge several files that contain intermediate
   * format data, and generate a single output file from it, which may be
   * again in intermediate format or any of the final formats. This latter
   * option is most helpful for parallel programs: as demonstrated in the
   * step-17 example program, it is possible to let only one processor
   * generate the graphical output for the entire parallel program, but this
   * can become vastly inefficient if many processors are involved, because
   * the load is no longer balanced. The way out is to let each processor
   * generate intermediate graphical output for its chunk of the domain, and
   * the later merge the different files into one, which is an operation that
   * is much cheaper than the generation of the intermediate data.
   *
   * Intermediate format deal.II data is usually stored in files with the
   * ending <tt>.d2</tt>.
   */
  template <int dim, int spacedim>
  void
  write_deal_II_intermediate(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                                    &nonscalar_data_ranges,
    const Deal_II_IntermediateFlags &flags,
    std::ostream                    &out);

  /**
   * Like write_deal_II_intermediate() but write all patches from all ranks
   * using MPI I/O
   * into a single file with name @p name. Compression using zlib is optional and controlled
   * by the @p compression argument.
   *
   * The files typically have the extension <tt>.pd2</tt>.
   */
  template <int dim, int spacedim>
  void
  write_deal_II_intermediate_in_parallel(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                                    &nonscalar_data_ranges,
    const Deal_II_IntermediateFlags &flags,
    const std::string               &filename,
    const MPI_Comm                   comm,
    const CompressionLevel           compression);

  /**
   * Write the data in @p data_filter to a single HDF5 file containing both the
   * mesh and solution values.
   */
  template <int dim, int spacedim>
  void
  write_hdf5_parallel(const std::vector<Patch<dim, spacedim>> &patches,
                      const DataOutFilter                     &data_filter,
                      const DataOutBase::Hdf5Flags            &flags,
                      const std::string                       &filename,
                      const MPI_Comm                           comm);

  /**
   * Write the data in @p data_filter to HDF5 file(s). If @p write_mesh_file is
   * false, the mesh data will not be written and the solution file will
   * contain only the solution values. If @p write_mesh_file is true and the
   * filenames are the same, the resulting file will contain both mesh data
   * and solution values.
   */
  template <int dim, int spacedim>
  void
  write_hdf5_parallel(const std::vector<Patch<dim, spacedim>> &patches,
                      const DataOutFilter                     &data_filter,
                      const DataOutBase::Hdf5Flags            &flags,
                      const bool                               write_mesh_file,
                      const std::string                       &mesh_filename,
                      const std::string &solution_filename,
                      const MPI_Comm     comm);

  /**
   * DataOutFilter is an intermediate data format that reduces the amount of
   * data that will be written to files. The object filled by this function
   * can then later be used again to write data in a concrete file format;
   * see, for example, DataOutBase::write_hdf5_parallel().
   */
  template <int dim, int spacedim>
  void
  write_filtered_data(
    const std::vector<Patch<dim, spacedim>> &patches,
    const std::vector<std::string>          &data_names,
    const std::vector<
      std::tuple<unsigned int,
                 unsigned int,
                 std::string,
                 DataComponentInterpretation::DataComponentInterpretation>>
                  &nonscalar_data_ranges,
    DataOutFilter &filtered_data);

  /**
   * Given an input stream that contains data written by
   * write_deal_II_intermediate(), determine the <tt>dim</tt> and
   * <tt>spacedim</tt> template parameters with which that function was
   * called, and return them as a pair of values.
   *
   * Note that this function eats a number of elements at the present position
   * of the stream, and therefore alters it. In order to read from it using,
   * for example, the DataOutReader class, you may wish to either reset the
   * stream to its previous position, or close and reopen it.
   */
  std::pair<unsigned int, unsigned int>
  determine_intermediate_format_dimensions(std::istream &input);

  /**
   * Return the OutputFormat value corresponding to the given string. If the
   * string does not match any known format, an exception is thrown.
   *
   * The main purpose of this function is to allow a program to use any
   * implemented output format without the need to extend the program's parser
   * each time a new format is implemented.
   *
   * To get a list of presently available format names, e.g. to give it to the
   * ParameterHandler class, use the function get_output_format_names().
   */
  OutputFormat
  parse_output_format(const std::string &format_name);

  /**
   * Return a list of implemented output formats. The different names are
   * separated by vertical bar signs (<tt>`|'</tt>) as used by the
   * ParameterHandler classes.
   */
  std::string
  get_output_format_names();

  /**
   * Provide a function which tells us which suffix a file with a given output
   * format usually has. At present the following formats are defined:
   * <ul>
   * <li> <tt>dx</tt>: <tt>.dx</tt>
   * <li> <tt>ucd</tt>: <tt>.inp</tt>
   * <li> <tt>gnuplot</tt>: <tt>.gnuplot</tt>
   * <li> <tt>povray</tt>: <tt>.pov</tt>
   * <li> <tt>eps</tt>: <tt>.eps</tt>
   * <li> <tt>gmv</tt>: <tt>.gmv</tt>
   * <li> <tt>tecplot</tt>: <tt>.dat</tt>
   * <li> <tt>vtk</tt>: <tt>.vtk</tt>
   * <li> <tt>vtu</tt>: <tt>.vtu</tt>
   * <li> <tt>svg</tt>: <tt>.svg</tt>
   * <li> <tt>deal_II_intermediate</tt>: <tt>.d2</tt>.
   * </ul>
   *
   * @deprecated Using Tecplot binary output is deprecated.
   */
  std::string
  default_suffix(const OutputFormat output_format);

  /**
   * @addtogroup Exceptions
   * @{
   */

  /**
   * Exception
   */
  DeclException2(ExcInvalidDatasetSize,
                 int,
                 int,
                 << "The number of points in this data set is " << arg1
                 << ", but we expected " << arg2
                 << " in each space direction.");
  /**
   * An output function did not receive any patches for writing.
   */
  DeclExceptionMsg(ExcNoPatches,
                   "You are trying to write graphical data into a file, but "
                   "no data is available in the intermediate format that "
                   "the DataOutBase functions require. Did you forget to "
                   "call a function such as DataOut::build_patches()?");
  /**
   * Exception
   */
  DeclExceptionMsg(ExcTecplotAPIError,
                   "The error code of one of the Tecplot functions was "
                   "not zero as expected.");
  /**
   * Exception
   */
  DeclException1(ExcErrorOpeningTecplotFile,
                 char *,
                 << "There was an error opening Tecplot file " << arg1
                 << " for output.");

  /** @} */
} // namespace DataOutBase



/**
 * This class is the interface to the functions in the DataOutBase namespace,
 * as already its name might suggest. It does not offer much functionality
 * apart from a way to access the implemented formats and a way to dynamically
 * dispatch what output format to chose.
 *
 * This class is thought as a base class to classes actually generating data
 * for output. It has two abstract virtual functions, get_patches() and
 * get_dataset_names() produce the data which is actually needed. These are
 * the only functions that need to be overloaded by a derived class. In
 * addition to that, it has a function for each output format supported by
 * the underlying base class which gets the output data using these two
 * virtual functions and passes them to the raw output functions.
 *
 * The purpose of this class is mainly two-fold: to support storing flags by
 * which the output in the different output formats are controlled, and means
 * to work with output in a way where output format, flags and other things
 * are determined at run time. In addition to that it offers the abstract
 * interface to derived classes briefly discussed above.
 *
 *
 * <h3>Output flags</h3>
 *
 * The way we treat flags in this class is very similar to that used in the
 * <tt>GridOut</tt> class. For detailed information on the why's and how's, as
 * well as an example of programming, we refer to the documentation of that
 * class.
 *
 * Basically, this class stores a set of flags for each output format
 * supported by the underlying <tt>DataOutBase</tt> class. These are used
 * whenever one of the <tt>write_*</tt> functions is used. By default, the
 * values of these flags are set to reasonable start-ups, but in case you want
 * to change them, you can create a structure holding the flags for one of the
 * output formats and set it using the <tt>set_flags</tt> functions of this
 * class to determine all future output the object might produce by that
 * output format.
 *
 * For information on what parameters are supported by different output
 * functions, please see the documentation of the <tt>DataOutBase</tt> class
 * and its member classes.
 *
 *
 * <h3>Run time selection of output parameters</h3>
 *
 * In the output flags classes, described above, many flags are defined for
 * output in the different formats. In order to make them available to the
 * input file handler class <tt>ParameterHandler</tt>, each of these has a
 * function declaring these flags to the parameter handler and to read them
 * back from an actual input file. In order to avoid that in user programs
 * these functions have to be called for each available output format and the
 * respective flag class, the present <tt>DataOutInterface</tt> class offers a
 * function <tt>declare_parameters</tt> which calls the respective function of
 * all known output format flags classes. The flags of each such format are
 * packed together in a subsection in the input file. Likewise, there is a
 * function <tt>parse_parameters</tt> which reads these parameters and stores
 * them in the flags associated with this object (see above).
 *
 * Using these functions, you do not have to track which formats are presently
 * implemented.
 *
 * Usage is as follows:
 * @code
 *   // within function declaring parameters:
 *   prm.enter_subsection("Output format options");
 *   DataOutInterface<dim>::declare_parameters(prm);
 *   prm.leave_subsection();
 *
 *   ...
 *   // within function doing the output:
 *   DataOut<dim> out;
 *   prm.enter_subsection("Output format options");
 *   out.parse_parameters(prm);
 *   prm.leave_subsection();
 * @endcode
 * Note that in the present example, the class <tt>DataOut</tt> was used.
 * However, any other class derived from <tt>DataOutInterface</tt> would work
 * alike.
 *
 *
 * <h3>Run time selection of formats</h3>
 *
 * This class, much like the <tt>GridOut</tt> class, has a set of functions
 * providing a list of supported output formats, an <tt>enum</tt> denoting all
 * these and a function to parse a string and return the respective
 * <tt>enum</tt> value if it is a valid output format's name (actually, these
 * functions are inherited from the base class). Finally, there is a function
 * <tt>write</tt>, which takes a value of this <tt>enum</tt> and dispatches to
 * one of the actual <tt>write_*</tt> functions depending on the output format
 * selected by this value.
 *
 * The functions offering the different output format names are, respectively,
 * <tt>default_suffix</tt>, <tt>parse_output_format</tt>, and
 * <tt>get_output_format_names</tt>. They make the selection of output formats
 * in parameter files much easier, and especially independent of the formats
 * presently implemented. User programs need therefore not be changed whenever
 * a new format is implemented.
 *
 * Additionally, objects of this class have a default format, which can be set
 * by the parameter "Output format" of the parameter file. Within a program,
 * this can be changed by the member function <tt>set_default_format</tt>.
 * Using this default format, it is possible to leave the format selection
 * completely to the parameter file. A suitable suffix for the output file
 * name can be obtained by <tt>default_suffix</tt> without arguments.
 *
 * @ingroup output
 */
template <int dim, int spacedim = dim>
class DataOutInterface
{
public:
  /**
   * Constructor.
   */
  DataOutInterface();

  /**
   * Destructor. Does nothing, but is declared virtual since this class has
   * virtual functions.
   */
  virtual ~DataOutInterface() = default;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in OpenDX
   * format. See DataOutBase::write_dx.
   */
  void
  write_dx(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in EPS
   * format. See DataOutBase::write_eps.
   */
  void
  write_eps(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in GMV
   * format. See DataOutBase::write_gmv.
   */
  void
  write_gmv(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in GNUPLOT
   * format. See DataOutBase::write_gnuplot.
   */
  void
  write_gnuplot(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in POVRAY
   * format. See DataOutBase::write_povray.
   */
  void
  write_povray(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in Tecplot
   * format. See DataOutBase::write_tecplot.
   */
  void
  write_tecplot(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in UCD
   * format for AVS. See DataOutBase::write_ucd.
   */
  void
  write_ucd(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in Vtk
   * format. See DataOutBase::write_vtk.
   *
   * @note VTK is a legacy format and has largely been supplanted by the VTU
   * format (an XML-structured version of VTK). In particular, VTU allows for
   * the compression of data and consequently leads to much smaller file sizes
   * that equivalent VTK files for large files. Since all visualization
   * programs that support VTK also support VTU, you should consider using the
   * latter file format instead, by using the write_vtu() function.
   */
  void
  write_vtk(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in Vtu
   * (VTK's XML) format. See DataOutBase::write_vtu.
   *
   * Some visualization programs, such as ParaView, can read several separate
   * VTU files to parallelize visualization. In that case, you need a
   * <code>.pvtu</code> file that describes which VTU files form a group. The
   * DataOutInterface::write_pvtu_record() function can generate such a
   * centralized record. Likewise, DataOutInterface::write_visit_record() does
   * the same for older versions of VisIt (although VisIt can also read
   * <code>pvtu</code> records since version 2.5.1). Finally,
   * DataOutInterface::write_pvd_record() can be used to group together the
   * files that jointly make up a time dependent simulation.
   */
  void
  write_vtu(std::ostream &out) const;

  /**
   * Collective MPI call to write the solution from all participating nodes
   * (those in the given communicator) to a single compressed .vtu file on a
   * shared file system.  The communicator can be a sub communicator of the
   * one used by the computation. This routine uses MPI I/O to achieve high
   * performance on parallel filesystems. In order to use this function,
   * you need to be using a file system that supports parallel MPI I/O,
   * and you will get error messages about failed MPI calls if you do not.
   * Also see DataOutInterface::write_vtu().
   */
  void
  write_vtu_in_parallel(const std::string &filename, const MPI_Comm comm) const;

  /**
   * Some visualization programs, such as ParaView and VisIt, can read several
   * separate VTU files that all form part of the same simulation, in order to
   * parallelize visualization. In that case, you need a
   * <code>.pvtu</code> file that describes which VTU files (written, for
   * example, through the DataOutInterface::write_vtu() function) form a group.
   * The current function can generate such a centralized record.
   *
   * The central record file generated by this function
   * contains a list of (scalar or vector) fields that describes which
   * fields can actually be found in the individual files that comprise the set
   * of parallel VTU files along with the names of these files. This function
   * gets the names and types of fields through the get_dataset_names() and
   * get_nonscalar_data_ranges() functions of this class. The second argument
   * to this function specifies the names of the files that form the parallel
   * set.
   *
   * @note Use DataOutBase::write_vtu() and DataOutInterface::write_vtu()
   * for writing each piece. Also note that
   * only one parallel process needs to call the current function, listing the
   * names of the files written by all parallel processes.
   *
   * @note The use of this function is explained in step-40.
   *
   * @note In order to tell Paraview to group together multiple
   * <code>pvtu</code> files that each describe one time step of a time
   * dependent simulation, see the DataOutBase::write_pvd_record()
   * function.
   *
   * @note Older versions of VisIt (before 2.5.1), can not read
   * <code>pvtu</code> records. However, it can read visit records as written
   * by the write_visit_record() function.
   */
  void
  write_pvtu_record(std::ostream                   &out,
                    const std::vector<std::string> &piece_names) const;

  /**
   * This function writes several .vtu files and a .pvtu record in parallel
   * and constructs the filenames automatically. It is a combination of
   * DataOutInterface::write_vtu() or
   * DataOutInterface::write_vtu_in_parallel(), and
   * DataOutInterface::write_pvtu_record().
   *
   * For example, running
   * <code> write_vtu_with_pvtu_record("output/", "solution", 3, comm, 4, 2)
   * </code> on 10 processes generates the files
   * @code
   * output/solution_0003.0.vtu
   * output/solution_0003.1.vtu
   * output/solution_0003.pvtu
   * @endcode
   * where the `.0.vtu` file contains the output of the first half of the
   * processes grouped together, and the `.1.vtu` the data from the remaining
   * half.
   *
   * A specified @p directory and a @p filename_without_extension
   * form the first part of the filename. The filename is then extended with
   * a @p counter labeling the current timestep/iteration/etc., the processor ID,
   * and finally the .vtu/.pvtu ending. Since the number of timesteps to be
   * written depends on the application, the number of digits to be reserved in
   * the filename can be specified as parameter @p n_digits_for_counter, and the number
   * is not padded with leading zeros if this parameter is left at its default
   * value numbers::invalid_unsigned_int. If more than one file identifier
   * is needed (e.g. time step number and iteration counter of solver), the
   * last identifier is used as @p counter, while all other identifiers have to be
   * added to @p filename_without_extension when calling this function.
   *
   * In a
   * parallel setting, several files are typically written per time step. The
   * number of files written in parallel depends on the number of MPI processes
   * (see parameter @p mpi_communicator), and a
   * specified number of @p n_groups with default value 0. The background is that
   * VTU file output supports grouping files from several CPUs into a given
   * number of files using MPI I/O when writing on a parallel filesystem. The
   * default value of @p n_groups is 0, meaning that every MPI rank will write one
   * file. A value of 1 will generate one big file containing the solution over
   * the whole domain, while a larger value will create @p n_groups files (but not
   * more than there are MPI ranks). For all values other than `n_groups==0`,
   * this function calls write_vtu_in_parallel(); for this function to work you
   * need to be using a file system that supports parallel MPI I/O, and you will
   * get error messages about failed MPI calls if you do not.
   *
   * Note that only one processor needs to
   * generate the .pvtu file, where processor zero is chosen to take over this
   * job.
   *
   * The return value is the filename of the centralized file for the pvtu
   * record.
   *
   * @note The code simply combines the strings @p directory and
   * @p filename_without_extension, i.e., the user has to make sure that
   * @p directory contains a trailing character, e.g. "/", that separates the
   * directory from the filename.
   *
   * @note Use an empty string "" for the first argument if output is to be
   * written in the current working directory.
   */
  std::string
  write_vtu_with_pvtu_record(
    const std::string &directory,
    const std::string &filename_without_extension,
    const unsigned int counter,
    const MPI_Comm     mpi_communicator,
    const unsigned int n_digits_for_counter = numbers::invalid_unsigned_int,
    const unsigned int n_groups             = 0) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in SVG
   * format. See DataOutBase::write_svg.
   */
  void
  write_svg(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it to <tt>out</tt> in deal.II
   * intermediate format. See DataOutBase::write_deal_II_intermediate.
   *
   * Note that the intermediate format is what its name suggests: a direct
   * representation of internal data. It isn't standardized and will change
   * whenever we change our internal representation. You can only expect to
   * process files written in this format using the same version of deal.II
   * that was used for writing.
   */
  void
  write_deal_II_intermediate(std::ostream &out) const;

  /**
   * Obtain data through get_patches() and write it using MPI I/O in parallel
   * to the file @p filename in the parallel
   * deal.II intermediate format. See
   * DataOutBase::write_deal_II_intermediate_in_parallel().
   */
  void
  write_deal_II_intermediate_in_parallel(
    const std::string                  &filename,
    const MPI_Comm                      comm,
    const DataOutBase::CompressionLevel compression) const;

  /**
   * Create an XDMFEntry based on the data in the data_filter. This assumes
   * the mesh and solution data were written to a single file. See
   * write_xdmf_file() for an example of usage.
   */
  XDMFEntry
  create_xdmf_entry(const DataOutBase::DataOutFilter &data_filter,
                    const std::string                &h5_filename,
                    const double                      cur_time,
                    const MPI_Comm                    comm) const;

  /**
   * Create an XDMFEntry based on the data in the data_filter. This assumes
   * the mesh and solution data were written to separate files. See
   * write_xdmf_file() for an example of usage.
   */
  XDMFEntry
  create_xdmf_entry(const DataOutBase::DataOutFilter &data_filter,
                    const std::string                &h5_mesh_filename,
                    const std::string                &h5_solution_filename,
                    const double                      cur_time,
                    const MPI_Comm                    comm) const;

  /**
   * Write an XDMF file based on the provided vector of XDMFEntry objects.
   * Below is an example of how to use this function with HDF5 and the
   * DataOutFilter:
   *
   * @code
   * DataOutBase::DataOutFilterFlags flags(true, true);
   * DataOutBase::DataOutFilter data_filter(flags);
   * std::vector<XDMFEntry> xdmf_entries;
   * // Filter the data and store it in data_filter
   * data_out.write_filtered_data(data_filter);
   * // Write the filtered data to HDF5
   * data_out.write_hdf5_parallel(data_filter, "solution.h5", MPI_COMM_WORLD);
   * // Create an XDMF entry detailing the HDF5 file
   * auto new_xdmf_entry = data_out.create_xdmf_entry(data_filter,
   *                                                  "solution.h5",
   *                                                  simulation_time,
   *                                                  MPI_COMM_WORLD);
   * // Add the XDMF entry to the list
   * xdmf_entries.push_back(new_xdmf_entry);
   * // Create an XDMF file from all stored entries
   * data_out.write_xdmf_file(xdmf_entries, "solution.xdmf", MPI_COMM_WORLD);
   * @endcode
   */
  void
  write_xdmf_file(const std::vector<XDMFEntry> &entries,
                  const std::string            &filename,
                  const MPI_Comm                comm) const;

  /**
   * Write the data in @p data_filter to a single HDF5 file containing both the
   * mesh and solution values. Below is an example of how to use this function
   * with the DataOutFilter:
   *
   * @code
   * DataOutBase::DataOutFilterFlags flags(true, true);
   * DataOutBase::DataOutFilter data_filter(flags);
   * // Filter the data and store it in data_filter
   * data_out.write_filtered_data(data_filter);
   * // Write the filtered data to HDF5
   * data_out.write_hdf5_parallel(data_filter, "solution.h5", MPI_COMM_WORLD);
   * @endcode
   */
  void
  write_hdf5_parallel(const DataOutBase::DataOutFilter &data_filter,
                      const std::string                &filename,
                      const MPI_Comm                    comm) const;

  /**
   * Write the data in data_filter to HDF5 file(s). If write_mesh_file is
   * false, the mesh data will not be written and the solution file will
   * contain only the solution values. If write_mesh_file is true and the
   * filenames are the same, the resulting file will contain both mesh data
   * and solution values.
   */
  void
  write_hdf5_parallel(const DataOutBase::DataOutFilter &data_filter,
                      const bool                        write_mesh_file,
                      const std::string                &mesh_filename,
                      const std::string                &solution_filename,
                      const MPI_Comm                    comm) const;

  /**
   * DataOutFilter is an intermediate data format that reduces the amount of
   * data that will be written to files. The object filled by this function
   * can then later be used again to write data in a concrete file format;
   * see, for example, DataOutBase::write_hdf5_parallel().
   */
  void
  write_filtered_data(DataOutBase::DataOutFilter &filtered_data) const;


  /**
   * Write data and grid to <tt>out</tt> according to the given data format.
   * This function simply calls the appropriate <tt>write_*</tt> function. If
   * no output format is requested, the <tt>default_format</tt> is written.
   *
   * An error occurs if no format is provided and the default format is
   * <tt>default_format</tt>.
   */
  void
  write(std::ostream                   &out,
        const DataOutBase::OutputFormat output_format =
          DataOutBase::default_format) const;

  /**
   * Set the default format. The value set here is used anytime, output for
   * format <tt>default_format</tt> is requested.
   */
  void
  set_default_format(const DataOutBase::OutputFormat default_format);


  /**
   * Set the flags to be used for output. This method expects <tt>flags</tt>
   * to be a member of one of the child classes of <tt>OutputFlagsBase</tt>.
   */
  template <typename FlagType>
  void
  set_flags(const FlagType &flags);


  /**
   * A function that returns the same string as the respective function in the
   * base class does; the only exception being that if the parameter is
   * omitted, then the value for the present default format is returned, i.e.
   * the correct suffix for the format that was set through
   * set_default_format() or parse_parameters() before calling this function.
   */
  std::string
  default_suffix(const DataOutBase::OutputFormat output_format =
                   DataOutBase::default_format) const;

  /**
   * Declare parameters for all output formats by declaring subsections within
   * the parameter file for each output format and call the respective
   * <tt>declare_parameters</tt> functions of the flag classes for each output
   * format.
   *
   * Some of the declared subsections may not contain entries, if the
   * respective format does not export any flags.
   *
   * Note that the top-level parameters denoting the number of subdivisions
   * per patch and the output format are not declared, since they are only
   * passed to virtual functions and are not stored inside objects of this
   * type. You have to declare them yourself.
   */
  static void
  declare_parameters(ParameterHandler &prm);

  /**
   * Read the parameters declared in declare_parameters() and set the flags
   * for the output formats accordingly.
   *
   * The flags thus obtained overwrite all previous contents of the flag
   * objects as default-constructed or set by the set_flags() function.
   */
  void
  parse_parameters(ParameterHandler &prm);

  /**
   * Return an estimate for the memory consumption, in bytes, of this object.
   * This is not exact (but will usually be close) because calculating the
   * memory usage of trees (e.g., <tt>std::map</tt>) is difficult.
   */
  std::size_t
  memory_consumption() const;

protected:
  /**
   * This is the abstract function through which derived classes propagate
   * preprocessed data in the form of Patch structures (declared in the base
   * class DataOutBase) to the actual output function. You need to overload
   * this function to allow the output functions to know what they shall
   * print.
   */
  virtual const std::vector<DataOutBase::Patch<dim, spacedim>> &
  get_patches() const = 0;

  /**
   * Abstract virtual function through which the names of data sets are
   * obtained by the output functions of the base class.
   */
  virtual std::vector<std::string>
  get_dataset_names() const = 0;

  /**
   * This functions returns information about how the individual components of
   * output files that consist of more than one data set are to be
   * interpreted.
   *
   * It returns a list of index pairs and corresponding name and type indicating
   * which components of the output are to be considered vector- or
   * tensor-valued rather than just a collection of scalar data. The index pairs
   * are inclusive; for example, if we have a Stokes problem in 2d with
   * components (u,v,p), then the corresponding vector data range should be
   * (0,1), and the returned list would consist of only a single element with a
   * tuple such as (0,1,"velocity",component_is_part_of_vector).
   *
   * Since some of the derived classes do not know about non-scalar data, this
   * function has a default implementation that simply returns an empty
   * string, meaning that all data is to be considered a collection of scalar
   * fields.
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const;

  /**
   * Validate that the names of the datasets returned by get_dataset_names() and
   * get_nonscalar_data_ranges() are valid. This currently consists of checking
   * that names are not used more than once. If an invalid state is encountered,
   * an Assert() will be triggered in debug mode.
   */
  void
  validate_dataset_names() const;


  /**
   * The default number of subdivisions for patches. This is filled by
   * parse_parameters() and should be obeyed by build_patches() in derived
   * classes.
   */
  unsigned int default_subdivisions;

private:
  /**
   * Standard output format.  Use this format, if output format default_format
   * is requested. It can be changed by the <tt>set_format</tt> function or in
   * a parameter file.
   */
  DataOutBase::OutputFormat default_fmt;

  /**
   * Flags to be used upon output of OpenDX data. Can be changed by using the
   * <tt>set_flags</tt> function.
   */
  DataOutBase::DXFlags dx_flags;

  /**
   * Flags to be used upon output of UCD data. Can be changed by using the
   * <tt>set_flags</tt> function.
   */
  DataOutBase::UcdFlags ucd_flags;

  /**
   * Flags to be used upon output of GNUPLOT data. Can be changed by using the
   * <tt>set_flags</tt> function.
   */
  DataOutBase::GnuplotFlags gnuplot_flags;

  /**
   * Flags to be used upon output of POVRAY data. Can be changed by using the
   * <tt>set_flags</tt> function.
   */
  DataOutBase::PovrayFlags povray_flags;

  /**
   * Flags to be used upon output of EPS data in one space dimension. Can be
   * changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::EpsFlags eps_flags;

  /**
   * Flags to be used upon output of gmv data in one space dimension. Can be
   * changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::GmvFlags gmv_flags;

  /**
   * Flags to be used upon output of hdf5 data in one space dimension. Can be
   * changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::Hdf5Flags hdf5_flags;

  /**
   * Flags to be used upon output of Tecplot data in one space dimension. Can
   * be changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::TecplotFlags tecplot_flags;

  /**
   * Flags to be used upon output of vtk data in one space dimension. Can be
   * changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::VtkFlags vtk_flags;

  /**
   * Flags to be used upon output of svg data in one space dimension. Can be
   * changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::SvgFlags svg_flags;

  /**
   * Flags to be used upon output of deal.II intermediate data in one space
   * dimension. Can be changed by using the <tt>set_flags</tt> function.
   */
  DataOutBase::Deal_II_IntermediateFlags deal_II_intermediate_flags;
};



/**
 * A class that is used to read data written in deal.II intermediate format
 * back in, so that it can be written out in any of the other supported
 * graphics formats. This class has two main purposes:
 *
 * The first use of this class is so that application programs can defer the
 * decision of which graphics format to use until after the program has been
 * run. The data is written in intermediate format into a file, and later on
 * it can then be converted into any graphics format you wish. This may be
 * useful, for example, if you want to convert it to gnuplot format to get a
 * quick glimpse and later on want to convert it to OpenDX format as well to
 * get a high quality version of the data. The present class allows to read
 * this intermediate format back into the program, and allows it to be written
 * in any other supported format using the relevant functions of the base
 * class.
 *
 * The second use is mostly useful in parallel programs: rather than having
 * one central process generate the graphical output for the entire program,
 * one can let each process generate the graphical data for the cells it owns,
 * and write it into a separate file in intermediate format. Later on, all
 * these intermediate files can then be read back in and merged together, a
 * process that is fast compared to generating the data in the first place.
 * The use of the intermediate format is mostly because it allows separate
 * files to be merged, while this is almost impossible once the data has been
 * written out in any of the supported established graphics formats.
 *
 * This second use scenario is explained in some detail in the step-18 example
 * program.
 *
 * In order to read data back into this object, you have to know the template
 * parameters for the space dimension which were used when writing the
 * data. If this knowledge is available at compile time, then this is no
 * problem. However, if it is not (such as in a simple format converter), then
 * it needs to be figured out at run time, even though the compiler already
 * needs it at compile time. A way around using the
 * DataOutBase::determine_intermediate_format_dimensions() function.
 *
 * Note that the intermediate format is what its name suggests: a direct
 * representation of internal data. It isn't standardized and will change
 * whenever we change our internal representation. You can only expect to
 * process files written in this format using the same version of deal.II that
 * was used for writing.
 *
 * @ingroup input output
 */
template <int dim, int spacedim = dim>
class DataOutReader : public DataOutInterface<dim, spacedim>
{
public:
  /**
   * Read a sequence of patches as written previously by
   * <tt>DataOutBase::write_deal_II_intermediate</tt> and store them in the
   * present object. This overwrites any previous content.
   */
  void
  read(std::istream &in);

  /**
   * Read all data previously written using
   * DataOutBase::write_deal_II_intermediate_in_parallel() from all
   * MPI ranks into this data structure.
   */
  void
  read_whole_parallel_file(std::istream &in);

  /**
   * This function can be used to merge the patches read by the other object
   * into the patches that this present object stores. This is sometimes handy
   * if one has, for example, a domain decomposition algorithm where each
   * block is represented by a DoFHandler of its own, but one wants to output
   * the solution on all the blocks at the same time. Alternatively, it may
   * also be used for parallel programs, where each process only generates
   * output for its share of the cells, even if all processes can see all
   * cells.
   *
   * For this to work, the input files for the present object and the given
   * argument need to have the same number of output vectors, and they need to
   * use the same number of subdivisions per patch. The output will probably
   * look rather funny if patches in both objects overlap in space.
   *
   * If you call read() for this object after merging in patches, the previous
   * state is overwritten, and the merged-in patches are lost.
   *
   * This function will fail if either this or the other object did not yet
   * set up any patches.
   */
  void
  merge(const DataOutReader<dim, spacedim> &other);

  /**
   * Exception
   */
  DeclExceptionMsg(ExcIncompatibleDatasetNames,
                   "You are trying to merge two sets of patches for which the "
                   "declared names of the variables do not match.");
  /**
   * Exception
   */
  DeclExceptionMsg(ExcIncompatiblePatchLists,
                   "You are trying to merge two sets of patches for which the "
                   "number of subdivisions or the number of vector components "
                   "do not match.");
  /**
   * Exception
   */
  DeclException4(ExcIncompatibleDimensions,
                 int,
                 int,
                 int,
                 int,
                 << "Either the dimensions <" << arg1 << "> and <" << arg2
                 << "> or the space dimensions <" << arg3 << "> and <" << arg4
                 << "> do not match!");

protected:
  /**
   * This is the function through which this class propagates preprocessed
   * data in the form of Patch structures (declared in the base class
   * DataOutBase) to the actual output function.
   *
   * It returns the patches as read the last time a stream was given to the
   * read() function.
   */
  virtual const std::vector<dealii::DataOutBase::Patch<dim, spacedim>> &
  get_patches() const override;

  /**
   * Abstract virtual function through which the names of data sets are
   * obtained by the output functions of the base class.
   *
   * Return the names of the variables as read the last time we read a file.
   */
  virtual std::vector<std::string>
  get_dataset_names() const override;

  /**
   * This functions returns information about how the individual components of
   * output files that consist of more than one data set are to be
   * interpreted.
   *
   * It returns a list of index pairs and corresponding name indicating which
   * components of the output are to be considered vector-valued rather than
   * just a collection of scalar data. The index pairs are inclusive; for
   * example, if we have a Stokes problem in 2d with components (u,v,p), then
   * the corresponding vector data range should be (0,1), and the returned
   * list would consist of only a single element with a tuple such as
   * (0,1,"velocity").
   *
   * Since some of the derived classes do not know about vector data, this
   * function has a default implementation that simply returns an empty
   * string, meaning that all data is to be considered a collection of scalar
   * fields.
   */
  virtual std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
  get_nonscalar_data_ranges() const override;

private:
  /**
   * Arrays holding the set of patches as well as the names of output
   * variables, all of which we read from an input stream.
   */
  std::vector<dealii::DataOutBase::Patch<dim, spacedim>> patches;
  std::vector<std::string>                               dataset_names;

  /**
   * Information about whether certain components of the output field are to
   * be considered vectors.
   */
  std::vector<
    std::tuple<unsigned int,
               unsigned int,
               std::string,
               DataComponentInterpretation::DataComponentInterpretation>>
    nonscalar_data_ranges;
};



/**
 * A class to store relevant data to use when writing a lightweight XDMF
 * file. The XDMF file in turn points to heavy data files (such as HDF5)
 * where the actual simulation data is stored.
 * This allows flexibility in arranging the data, and also
 * allows the mesh to be separated from the point data.
 */
class XDMFEntry
{
public:
  /**
   * Default constructor that creates an invalid object.
   */
  XDMFEntry();

  /**
   * Simplified constructor that calls the complete constructor for
   * cases where <code>solution_filename == mesh_filename</code>, and
   * <code>dim==spacedim</code>.
   */
  XDMFEntry(const std::string   &filename,
            const double         time,
            const std::uint64_t  nodes,
            const std::uint64_t  cells,
            const unsigned int   dim,
            const ReferenceCell &cell_type);

  /**
   * Simplified constructor that calls the complete constructor for
   * cases where <code>dim==spacedim</code>.
   */
  XDMFEntry(const std::string   &mesh_filename,
            const std::string   &solution_filename,
            const double         time,
            const std::uint64_t  nodes,
            const std::uint64_t  cells,
            const unsigned int   dim,
            const ReferenceCell &cell_type);

  /**
   * Constructor that sets all members to provided parameters.
   */
  XDMFEntry(const std::string   &mesh_filename,
            const std::string   &solution_filename,
            const double         time,
            const std::uint64_t  nodes,
            const std::uint64_t  cells,
            const unsigned int   dim,
            const unsigned int   spacedim,
            const ReferenceCell &cell_type);

  /**
   * Record an attribute and associated dimensionality.
   */
  void
  add_attribute(const std::string &attr_name, const unsigned int dimension);

  /**
   * Read or write the data of this object for serialization using the
   * [BOOST serialization
   * library](https://www.boost.org/doc/libs/1_74_0/libs/serialization/doc/index.html).
   */
  template <class Archive>
  void
  serialize(Archive &ar, const unsigned int /*version*/)
  {
    ar &valid &h5_sol_filename &h5_mesh_filename &entry_time &num_nodes
      &num_cells &dimension &space_dimension &cell_type &attribute_dims;
  }

  /**
   * Get the XDMF content associated with this entry.
   * If the entry is not valid, this returns an empty string.
   */
  std::string
  get_xdmf_content(const unsigned int indent_level) const;

private:
  /**
   * Whether this entry is valid and contains data to be written.
   */
  bool valid;

  /**
   * The name of the HDF5 heavy data solution file this entry references.
   */
  std::string h5_sol_filename;

  /**
   * The name of the HDF5 mesh file this entry references.
   */
  std::string h5_mesh_filename;

  /**
   * The simulation time associated with this entry.
   */
  double entry_time;

  /**
   * The number of data nodes.
   */
  std::uint64_t num_nodes;

  /**
   * The number of data cells.
   */
  std::uint64_t num_cells;

  /**
   * The dimension associated with the data.
   */
  unsigned int dimension;

  /**
   * The dimension of the space the data lives in.
   * Note that dimension <= space_dimension.
   */
  unsigned int space_dimension;

  /**
   * The type of cell in deal.II language. We currently only support
   * xdmf entries where all cells have the same type.
   */
  ReferenceCell cell_type;

  /**
   * The attributes associated with this entry and their dimension.
   */
  std::map<std::string, unsigned int> attribute_dims;
};



/* -------------------- inline functions ------------------- */

namespace DataOutBase
{
  inline bool
  EpsFlags::RgbValues::is_grey() const
  {
    return (red == green) && (red == blue);
  }


  /* -------------------- template functions ------------------- */

  /**
   * Output operator for an object of type <tt>DataOutBase::Patch</tt>. This
   * operator dumps the intermediate graphics format represented by the patch
   * data structure. It may later be converted into regular formats for a
   * number of graphics programs.
   */
  template <int dim, int spacedim>
  std::ostream &
  operator<<(std::ostream &out, const Patch<dim, spacedim> &patch);



  /**
   * Input operator for an object of type <tt>DataOutBase::Patch</tt>. This
   * operator reads the intermediate graphics format represented by the patch
   * data structure, using the format in which it was written using the
   * operator<<.
   */
  template <int dim, int spacedim>
  std::istream &
  operator>>(std::istream &in, Patch<dim, spacedim> &patch);
} // namespace DataOutBase


DEAL_II_NAMESPACE_CLOSE

#endif
