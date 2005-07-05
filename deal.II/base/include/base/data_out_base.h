//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__data_out_base_h
#define __deal2__data_out_base_h


#include <base/config.h>
#include <base/point.h>
#include <base/table.h>
#include <grid/geometry_info.h>

#include <vector>
#include <string>

// Only include the Tecplot API header if the appropriate files
// were detected by configure
#ifdef DEAL_II_HAVE_TECPLOT
#  include "TECIO.h"
#  include <string.h>
#endif

// we only need output streams, but older compilers did not provide
// them in a separate include file
#ifdef HAVE_STD_OSTREAM_HEADER
#  include <ostream>
#else
#  include <iostream>
#endif


class ParameterHandler;


/**
 * This is a base class for output of data on meshes of very general
 * form. It basically only provides a set of functions for several output
 * formats which take a list of patches and write them to an output
 * stream.
 *
 * By offering this interface to the different output formats, it is simple
 * to extend this class to new formats without depending on such things
 * as actual triangulations and handling of data vectors. These things shall
 * be provided by derived class which have a user callable interface then.
 *
 *
 * @section DataOutBaseInterface Interface
 * This class has an interface that is not usually called by a user directly;
 * also, it consists of <tt>static</tt> functions only. Usually, derived classes will
 * inherit this class <tt>protected</tt> to hide this interface to the users of thes
 * classes.
 *
 * The interface of this class basically consists of the declaration of a data
 * type describing a patch and a bunch of functions taking a list of patches
 * and writing them in one format or other to the stream. It is in the
 * responsibility of the derived classes to provide this list of patches.
 * In addition to the list of patches, a name for each data set may be given.
 *
 *
 * @section QueryingP Querying interface
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
 *
 * @subsection DataOutBasePatches Patches
 * Grids can be thought of as a collection of cells; if you want to write out
 * data on such a grid, you can do so by writing them one cell at a time.
 * The functions in this class therefore take a list of objects describing the
 * data on one cell each. This data for each cell usually consists of a list
 * of vertices for this cell, and a list of data values (for example solution
 * data, error information, etc) at each of these vertices.
 *
 * In some cases, this interface to a cell is too restricted, however. For
 * example, you may have higher order elements and printing the values at
 * the vertices only is not enough. For this reason, we not only provide
 * writing the data on the vertices only, but the data is organizes as a
 * tensor product grid on each cell. The parameter <tt>n_subdivision</tt>, which is
 * given for each patch separately, denotes how often the cell is to be
 * divided for output; for example, <tt>n_subdivisions==1</tt> yields no subdivision
 * of the cell, <tt>n_subdivisions==2</tt> will produce a grid of 3 times 3 points
 * in two spatial dimensions and 3 times 3 times 3 points in three dimensions,
 * <tt>n_subdivisions==2</tt> will yield 4 times 4 (times 4) points, etc. The actual
 * location of these points on the patch will be computed by a multilinear
 * transformation from the vertices given for this patch.
 *
 * Given these comments, the actual data to be printed on this patch of
 * points consists of several data sets each of which has a value at each
 * of the patch points. For example with <tt>n_subdivisions==2</tt> in two space
 * dimensions, each data set has to provide nine values, and since the
 * patch is to be printed as a tensor product (or its transformation to the
 * real space cell), its values are to be ordered like
 * <i>(x0,y0) (x0,y1) (x0,y2) (x1,y0) (x1,y1) (x1,y2) (x2,y0) (x2,y1) (x2,y2)</i>,
 * i.e. the z-coordinate runs fastest, then the y-coordinate, then x (if there
 * are that many space directions).
 *
 *
 * @subsection DataOutBaseGP Generalized patches
 *
 * In general, the patches as explained above might be too
 * restricted. For example, one might want to draw only the outer
 * faces of a domain in a three-dimensional computation, if one is not
 * interested in what happens inside. Then, the objects that should be
 * drawn are two-dimensional in a three-dimensional world. The
 * Patch class and associated output functions handle these
 * cases. The Patch class therefore takes two template parameters,
 * the first, named <tt>dim</tt> denoting the dimension of the object (in
 * the above example, this would be two), while the second, named
 * <tt>spacedim</tt>, denotes the dimension of the embedding space (this
 * would be three). The corner points of a patch have the dimension of
 * the space, while their number is determined by the dimension of the
 * patch. By default, the second template parameter has the same value
 * as the first, which would correspond to outputting a cell, rather
 * than a face or something else.
 *
 *
 * @section DataOutBaseFormats Supported output formats
 *
 * @subsection DataOutBaseOpenDX OpenDX (IBM Open Visualization Data Explorer}
 *
 * Since Data Explorer (DX) is distributed as OpenSource, there is a
 * well-accessible visualization tool for all (at least Unix-based)
 * platforms. Therefore, output in its natural file format is
 * included.
 *
 * 
 * @subsection DataOutBaseUCD AVS UCD format
 *
 * The UCD format is described in the AVS developer's guide. Due to
 * limitations in the present format, only node based data can be output,
 * which in one reason why we invented the patch concept. In order to
 * write higher order elements, you may split them up into several subdivisions
 * of each cell. These subcells will then, however, also appear as different
 * cells by programs which understand the UCD format.
 * 
 * No use is made of the possibility to give model data since these
 * are not supported by all UCD aware programs. You may give cell data
 * in derived classes by setting all values of a given data set on a
 * patch to the same value.
 *
 *
 * @subsection DataOutBaseGNUPLOT GNUPLOT format
 *
 * The GNUPLOT format is not able to handle data on unstructured grids
 * directly. Directly would mean that you only give the vertices and
 * the solution values thereon and the program constructs its own grid
 * to represent the data. This is only possible for a structured tensor
 * product grid in two dimensions. However, it is possible to give several
 * such patches within one file, which is exactly what the respective
 * function of this class does: writing each cell's data as a patch of
 * data, at least if the patches as passed from derived classes
 * represent cells. Note that the functions on patches need not be
 * continuous at interfaces between patches, so this method also works
 * for discontinuous elements. Note also, that GNUPLOT can do hidden
 * line removal for patched data.
 *
 * While this discussion applies to two spatial dimensions, it is more
 * complicated in 3d. The reason is that we could still use patches, but
 * it is difficult when trying to visualize them, since if we use a cut
 * through the data (by, for example, using x- and z-coordinates, a fixed
 * y-value and plot function values in z-direction, then the patched data
 * is not a patch in the sense GNUPLOT wants it any more. Therefore, we use
 * another approach, namely writing the data on the 3d grid as a sequence
 * of lines, i.e. two points each associated with one or more data sets.
 * There are therefore 12 lines for each subcells of a patch.
 *
 * Given the lines as described above, a cut through this data in Gnuplot
 * can then be achieved like this:
 * @verbatim
 *   set data style lines
 *   splot [:][:][0:] "T" using 1:2:($3==.5 ? $4 : -1)
 * @endverbatim
 * This command plots data in x- and y-direction unbounded, but in z-direction
 * only those data points which are above the x-y-plane (we assume here a
 * positive solution, if it has negative values, you might want to decrease the
 * lower bound). Furthermore, it only takes the data points with z-values (<tt>$3</tt>)
 * equal to 0.5, i.e. a cut through the domain at <tt>z=0.5</tt>. For the data points
 * on this plane, the data values of the first data set (<tt>$4</tt>) are raised in
 * z-direction above the x-y-plane; all other points are denoted the value
 * <tt>-1</tt> instead of the value of the data vector and are not plotted due to
 * the lower bound in z plotting direction, given in the third pair of brackets.
 *
 * Of course, more complex cuts are possible, including nonlinear
 * ones. Note however, that only those points which are actually on the
 * cut-surface are plotted.
 *
 *
 * @subsection DataOutBasePOVRAY POVRAY format
 *
 * Output in this format creates a povray source file, include standard
 * camera and light source definition for rendering with povray 3.1
 * At present, this format only supports output for two-dimensional data,
 * with values in the third direction taken from a data vector.
 *
 * The output uses two different povray-objects:
 *
 * <ul>
 * <li> <tt>BICUBIC_PATCH</tt>
 * A <tt>bicubic_patch</tt> is a 3-dimensional Bezier patch. It consists of 16 Points
 * describing the surface. The 4 corner points are touched by the object,
 * while the other 12 points pull and stretch the patch into shape.
 * One <tt>bicubic_patch</tt> is generated on each patch. Therefor the number of 
 * subdivisions has to be 3 to provide the patch with 16 points.
 * A bicubic patch is not exact but generates very smooth images.
 *
 * <li> <tt>MESH</tt>
 * The mesh object is used to store large number of triangles.
 * Every square of the patch data is split into one upper-left and one 
 * lower-right triangle. If the number of subdivisions is three, 32 triangle
 * are generated for every patch.
 * 
 * Using the smooth flag povray interpolates the normals on the triangles,
 * imitating a curved surface
 * </ul>
 *
 * All objects get one texture definition called Tex. This texture has to be
 * declared somewhere before the object data. This may be in an external 
 * data file or at the beginning of the output file.
 * Setting the <tt>external_data</tt> flag to false, an standard camera, light and
 * texture (scaled to fit the scene) is added to the outputfile. Set to true
 * an include file "data.inc" is included. This file is not generated by deal
 * and has to include camera, light and the texture definition Tex.
 *
 * You need povray (>=3.0) to render the scene. The minimum options for povray
 * are:
 * @verbatim
 *   povray +I<inputfile> +W<horiz. size> +H<ver. size> +L<include path>
 * @endverbatim
 * If the external file "data.inc" is used, the path to this file has to be
 * included in the povray options.
 *
 *
 * @subsection DataOutBaseEPS EPS (encapsulated Postscript<sup>TM</sup> format
 *
 * Output in this format circumvents the use of auxiliary graphic programs
 * converting some output format into a graphics format. This has the advantage
 * that output is easy and fast, and the disadvantage that you have to give a
 * whole bunch of parameters which determine the direction of sight, the mode of
 * colorization, the scaling of the height axis, etc. (Of course, all these
 * parameters have reasonable default values, which you may want to change from
 * time to time.) At present, this format only supports output for two-dimensional
 * data, with values in the third direction taken from a data vector.
 *
 * Basically, output consists of the mesh and the cells in between them. You can
 * draw either of these, or both, or none if you are really interested in an empty
 * picture. If written, the mesh uses black lines. The cells in between the mesh
 * are either not printed (this will result in a loss of hidden line removal, i.e.
 * you can "see through" the cells to lines behind), printed in white (which does
 * nothing apart from the hidden line removal), or colorized using one of the
 * data vectors (which need not be the same as the one used for computing the
 * height information) and a customizable color function. The default color
 * functions chooses the color between black, blue, green, red and white, with
 * growing values of the data field chosen for colorization. At present, cells
 * are displayed with one color per cell only, which is taken from the value of
 * the data field at the center of the cell; bilinear interpolation of the color
 * on a cell is not used.
 *
 * By default, the viewpoint is chosen like the default viewpoint in GNUPLOT, i.e.
 * with an angle of 60 degrees with respect to the positive z-axis and rotated
 * 30 degrees in positive sense (as seen from above) away from the negative y-axis.
 * Of course you can change these settings.
 *
 * EPS output is written without a border around the picture, i.e. the bounding box
 * is close to the output on all four sides. Coordinates are written using at most
 * five digits, to keep picture size at a reasonable size.
 *
 * All parameters along with their default values are listed in the documentation
 * of the <tt>EpsFlags</tt> member class of this class. See there for more and detailed
 * information.
 *
 * Please note that due to the various transformations each patch has to undergo
 * before actual outut, memory requirements may be rather large for large numbers
 * of patches.
 *
 *
 * @subsection DataOutBaseGMV GMV format
 *
 * The write_gmv() function writes the data in a format understood by
 * the GMV (general mesh viewer) program. This program is able to
 * generate 2d and 3d plots of almost arbitrarily many data sets,
 * along with shading, cuts through data sets and many other nifty
 * features.
 *
 * Data is written in the following format: nodes are considered the points
 * of the patches. In spatial dimensions less than three, zeroes are
 * inserted for the missing coordinates. The data vectors are written as
 * node or cell data, where for the first the data space is interpolated to
 * (bi-,tri-)linear elements.
 *
 *
 * @subsection DataOutBaseTecplot Tecplot format
 *
 * The write_tecplot() function writes the data in <a
 * href="http://www.amtec.com">Tecplot</a> FEBLOCK format. The program
 * supports 1, 2, and 3D data and has features such as contouring,
 * slicing, drawing streamlines, and animation. Patches are written as
 * a collection of quadrilaterals in 2D or bricks in 3D, with the
 * nodal values interpolated to (bi-,tri-) linear elements. These
 * functions will write Tecplot ASCII formatted files.
 * 
 * Additionally, Tecplot binary output is supported through
 * write_tecplot_binary().  For this to work properly
 * <tt>./configure</tt> checks for the Tecplot API at build time. To
 * write Tecplot binary files directly make sure that the TECHOME
 * environment variable points to the Tecplot installation directory,
 * and that the files @$TECHOME/include/TECIO.h and
 * @$TECHOME/lib/tecio.a are readable.  If these files are not
 * availabe (or in the case of 1D) <tt>write_tecplot_binary</tt> will
 * simply call <tt>write_tecplot</tt> and thus larger ASCII data files
 * will be produced rather than more efficient Tecplot binary files.
 * For more information consult the Tecplot Users and Reference
 * manuals.
 *
 *
 *
 * @subsection DataOutBaseVTK VTK format
 *
 * This is the file format used by the Visualization Toolkit <a
 * href="http://www.kitware.com/vtk.html">VTK</a>, as described in
 * their manual, section 14.3. It is similar to the GMV format, see
 * there for more information.
 *
 *
 * @subsection DataOutBaseD2 deal.II intermediate format
 *
 * In addition to all the other formats, this class can also write
 * data in the deal.II format. This is not a format that read by any
 * other graphics program, but is rather a direct dump of the
 * intermediate internal format used by deal.II. This internal format
 * is generated by the various classes that can generate output using
 * the DataOutBase class, for example from a finite element solution,
 * and is then converted in the present class to the final graphics
 * format. The reason why we offer to write out this intermediate
 * format is that it can be read back into a deal.II program using the
 * DataOutReader class, which is helpful in at least two contexts:
 * First, this can be used to later generate graphical output in any
 * other graphics format presently understood; this way, it is not
 * necessary to know at run-time which output format is requested, or
 * if multiple output files in different formats are needed. Secondly,
 * in contrast to almost all other graphics formats, it is possible to
 * merge several files that contain intermediate format data, and
 * generate a single output file from it, which may be again in
 * intermediate format or any of the final formats. This latter option
 * is most helpful for parallel programs: as demonstrated in the
 * step-17 example program, it is possible to let only one processor
 * generate the graphical output for the entire parallel program, but
 * this can become vastly inefficient if many processors are involved,
 * because the load is no longer balanced. The way out is to let each
 * processor generate intermediate graphical output for its chunk of
 * the domain, and the later merge the different files into one, which
 * is an operation that is much cheaper than the generation of the
 * intermediate data.
 *
 * Intermediate format deal.II data is usually stored in files with
 * the ending <tt>.d2</tt>.
 *
 *
 * @section DataOutBaseOP Output parameters
 *
 * All functions take a parameter which is a structure of type <tt>XFlags</tt>, where
 * <tt>X</tt> is the name of the output format. To find out what flags are presently
 * supported, read the documentation of the different structures.
 *
 * Note that usually the output formats used for scientific visualization
 * programs have no or very few parameters (apart from some compatibility flags)
 * because there the actual appearance of output is determined using the
 * visualization program and the files produced by this class store more or less
 * only raw data.
 *
 * The direct output formats, like Postscript or Povray need to be given a lot
 * more parameters, though, since there the output file has to contain all
 * details of the viewpoint, light source, etc.
 *
 *
 * @author Wolfgang Bangerth 1999, 2000, 2001; EPS output based on an earlier implementation by Stefan Nauber for the old DataOut class; Povray output by Thomas Richter 1999; OpenDX output by Guido Kanschat, 2001; Tecplot output by Benjamin Shelton Kirk, 2002.
 */
class DataOutBase 
{
  public:
				     /**
				      * Data structure describing a
				      * patch of data in <tt>dim</tt> space
				      * dimensions. See the general
				      * documentation of the
				      * <tt>DataOutBase</tt> class for more
				      * information on its contents
				      * and purposes.  In the case of
				      * two dimensions, the next
				      * picture ist an example of
				      * <tt>n_subdivision</tt> = 4 because
				      * the number of cells is
				      * equal to <tt>2^dim</tt>.
				      * @verbatim
				      *  __ __ __ __
				      * |  |  |  |  |
				      * |__|__|__|__| 
				      * |  |  |  |  |
				      * |__|__|__|__|
				      * |  |  |  |  |
				      * |__|__|__|__| 
				      * |  |  |  |  |
				      * |__|__|__|__|
				      * @endverbatim
				      * @author Wolfgang Bangerth
				      */
    template <int dim, int spacedim=dim>
    struct Patch
    {
					 /**
					  * Corner points of a patch.
					  * Inner points are computed by
					  * a multilinear transform of
					  * the unit cell to the cell
					  * specified by these corner
					  * points. The order of points
					  * is the same as for cells
					  * in the triangulation.
					  */
	Point<spacedim> vertices[GeometryInfo<dim>::vertices_per_cell];

					 /**
					  * Numbers of neighbors of a patch.
					  * OpenDX format requires
					  * neighbor information for
					  * advanced output. Here the
					  * neighborship relationship
					  * of patches is
					  * stored. During output,
					  * this must be transformed
					  * into neighborship of
					  * sub-grid cells.
					  */
	unsigned int neighbors[GeometryInfo<dim>::faces_per_cell];

					 /**
					  * Number of this
					  * patch. Since we are not
					  * sure patches are handled
					  * in the same order, always,
					  * we better store this.
					  */
	unsigned int patch_index;
	
					 /**
					  * Number of subdivisions with
					  * which this patch is to be
					  * written. <tt>1</tt> means no
					  * subdivision, <tt>2</tt> means
					  * bisection, <tt>3</tt> trisection,
					  * etc.
					  */
	unsigned int n_subdivisions;
	
					 /**
					  * Data vectors. The format is
					  * as follows:
					  * <tt>data(i,.)</tt> denotes the data
					  * belonging to the <tt>i</tt>th data
					  * vector. <tt>data.n()</tt>
					  * therefore equals the number
					  * of output points; this
					  * number is <tt>(subdivisions+1)^{dim</tt>}.
					  * <tt>data.m()</tt> equals the number of
					  * data vectors.
					  *
					  * Within each column,
					  * <tt>data(.,j)</tt> are the data
					  * values at the output point
					  * <tt>j</tt>, where <tt>j</tt> runs first
					  * over the last direction,
					  * then over the second last
					  * one etc, just as if it was
					  * organized as an array
					  * <tt>double[x][y][z]</tt>. This is
					  * also the order of points
					  * as provided by the
					  * <tt>QIterated</tt> class when
					  * used with the <tt>QTrapez</tt>
					  * class as subquadrature.
					  * Note that if
					  * <tt>subdivisions==1</tt>, the
					  * elements of <tt>data[i]</tt>
					  * correspond to vertices
					  * <tt>(0,1)</tt> in 1d, 
					  * <tt>(0, 3, 1,2)</tt> in 2d, and 
					  * <tt>(0, 4, 3, 7, 1, 5, 2, 6)</tt>
					  * in 3d as following:
					  * @verbatim
					  *  
					  *      7________6
					  *      /       /|
					  *     /       / |
					  *   3/______2/  |   
					  *   |   |   |   |
					  |   |   4___|___5
					  *   |  /    |  /
					  *   | /     | /
					  *  0|/______1/
					  * @endverbatim
					  * 
					  * For exemple in 2d: If
					  * <tt>subdivisions==2</tt> the
					  * elements of <tt>data[i]</tt> are
					  * given by the following
					  * numeration:
					  *
					  * @verbatim
					  *  2 ____5 ____8
					  *   |     |     |
					  *   |     |     | 
					  *   |     |     | 
					  *  1|____4|____7|
					  *   |     |     |
					  *   |     |     | 
					  *   |     |     | 
					  *  0|____3|____6|
					  * @endverbatim
					  *
					  * Since the number of data vectors
					  * is usually the same for all
					  * patches to be printed, <tt>data.size()</tt>
					  * should yield the same value for all
					  * patches provided.
					  */
	Table<2,float> data;
	
					 /**
					  * Default constructor. Sets
					  * <tt>n_subdivisions</tt> to one.
					  */
	Patch ();

					 /**
					  * Compare the present patch
					  * for equality with another
					  * one. This is used in a few
					  * of the automated tests in
					  * our testsuite.
					  */
	bool operator == (const Patch &patch) const;
	
					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
	
					 /**
					  * Value to be used if this
					  * patch has no neighbor on
					  * one side.
					  */
	static const unsigned int no_neighbor = deal_II_numbers::invalid_unsigned_int;
					 /** @addtogroup Exceptions
					  * @{ */
	
					 /**
					  * Exception
					  */
	DeclException2 (ExcInvalidCombinationOfDimensions,
			int, int,
			<< "It is not possible to have a structural dimension of " << arg1
			<< " to be larger than the space dimension of the surrounding"
			<< " space " << arg2);
					 //@}
    };

    				     /**
				      * Flags describing the details of
				      * output in OpenDX format. At
				      * present no flags are implemented.
				      */
    struct DXFlags 
    {
					 /**
					  * Write in multigrid format.
					  * This will be necessary to
					  * get streamlines right on
					  * locally refined grids.
					  */
	bool write_multigrid;

					 /**
					  * Write neighbor information.
					  */
	bool write_neighbors;

					 /**
					  * Constructor.
					  */
	DXFlags (const bool write_multigrid = false,
		 const bool write_neighbors = false);
        
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object.
					  */
	unsigned int memory_consumption () const;
    };

				     /**
				      * Flags describing the details of
				      * output in UCD format.
				      */
    struct UcdFlags 
    {
					 /**
					  * Write a comment at the beginning
					  * of the file stating the date of
					  * creation and some other data.
					  * While this is supported by the
					  * UCD format (and the AVS program),
					  * some other programs get confused
					  * by this, so you can switch it off
					  * this way.
					  *
					  * Default: <tt>true</tt>.
					  */
	bool write_preamble;
	
					 /**
					  * Constructor.
					  */
	UcdFlags (const bool write_preamble = true);

					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };

				     /**
				      * Flags describing the details of
				      * output in Gnuplot format. At
				      * present no flags are implemented.
				      */
    struct GnuplotFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * <tt>private</tt> as well).
					  */
	int dummy;

      public:
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };

    				     /**
				      * Flags describing the details
				      * of output in Povray
				      * format. Several flags are
				      * implemented, see their
				      * respective documentation.
				      */
    struct PovrayFlags 
    {
					 /**
					  * Normal vector interpolation,
					  * if set to true
					  *
					  * default = false
					  */
	bool smooth;
	
					 /**
					  * Use bicubic patches (b-splines)
					  * instead of triangles.
					  *
					  * default = false
					  */
	bool bicubic_patch;

					 /**
					  * include external "data.inc"
					  * with camera, light and
					  * texture definition for the
					  * scene.
					  *
					  * default = false
					  */
	bool external_data;
	
					 /**
					  * Constructor.
					  */
	PovrayFlags (const bool smooth = false,
		     const bool bicubic_patch = false,
		     const bool external_data = false);
	
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };


				     /**
				      * Flags describing the details of
				      * output in encapsulated postscript
				      * format.
				      */
    struct EpsFlags 
    {
					 /**
					  * This denotes the number of the
					  * data vector which shall be used
					  * for generating the height
					  * information. By default, the
					  * first data vector is taken,
					  * i.e. <tt>height_vector==0</tt>, if
					  * there is any data vector. If there
					  * is no data vector, no height
					  * information is generated.
					  */
	unsigned int height_vector;

					 /**
					  * Number of the vector which is
					  * to be taken to colorize cells.
					  * The same applies as for
					  * <tt>height_vector</tt>.
					  */
	unsigned int color_vector;
	
					 /**
					  * Enum denoting the possibilities
					  * whether the scaling should be done
					  * such that the given <tt>size</tt> equals
					  * the width or the height of
					  * the resulting picture.
					  */
	enum SizeType {
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
					  * Width or height of the output
					  * as given in postscript units
					  * This usually is given by the
					  * strange unit 1/72 inch. Whether
					  * this is height or width is
					  * specified by the flag
					  * <tt>size_type</tt>.
					  *
					  * Default is 300, which represents
					  * a size of roughly 10 cm.
					  */
	unsigned int size;

					 /**
					  * Width of a line in postscript
					  * units. Default is 0.5.
					  */
	double line_width;
	
					 /**
					  * Angle of the line origin-viewer
					  * against the z-axis in degrees.
					  *
					  * Default is the Gnuplot-default
					  * of 60.
					  */
	double azimut_angle;

					 /**
					  * Angle by which the viewers
					  * position projected onto the
					  * x-y-plane is rotated around
					  * the z-axis, in positive sense
					  * when viewed from above. The
					  * unit are degrees, and zero
					  * equals a position above or below
					  * the negative y-axis.
					  *
					  * Default is the
					  * Gnuplot-default of 30.
					  * An exemple of a
					  * Gnuplot-default of 0 is
					  * the following:
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
					  * Factor by which the z-axis is to
					  * be stretched as compared to the
					  * x- and y-axes. This is to compensate
					  * for the different sizes that
					  * coordinate and solution values may
					  * have and to prevent that the plot
					  * looks to much out-of-place (no
					  * elevation at all if solution values
					  * are much smaller than coordinate
					  * values, or the common "extremely
					  * mountainous area" in the opposite
					  * case.
					  *
					  * Default is <tt>1.0</tt>.
					  */
	double z_scaling;

					 /**
					  * Flag the determines whether the
					  * lines bounding the cells (or the
					  * parts of each patch) are to be
					  * plotted.
					  *
					  * Default: <tt>true</tt>.
					  */
	bool   draw_mesh;

					 /**
					  * Flag whether to fill the regions
					  * between the lines bounding the cells
					  * or not. If not, no hidden line removal
					  * is performed, which in this crude
					  * implementation is done through
					  * writing the cells in a back-to-front
					  * order, thereby hiding the cells in
					  * the background by cells in the
					  * foreground.
					  *
					  * If this flag is <tt>false</tt> and <tt>draw_mesh</tt>
					  * is <tt>false</tt> as well, nothing will be
					  * printed.
					  *
					  * If this flag is <tt>true</tt>, then the cells
					  * will be drawn either colored by one
					  * of the data sets (if <tt>shade_cells</tt> is
					  * <tt>true</tt>), or pure white (if
					  * <tt>shade_cells</tt> is false or if there are
					  * no data sets).
					  *
					  * Default is <tt>true</tt>.
					  */
	bool   draw_cells;

					 /**
					  * Flag to determine whether the cells
					  * shall be colorized by the data
					  * set denoted by <tt>color_vector</tt>, or
					  * simply be painted in white. This
					  * flag only makes sense if
					  * <tt>draw_cells==true</tt>. Colorization is
					  * done through the <tt>color_function</tt>.
					  *
					  * Default is <tt>true</tt>.
					  */
	bool   shade_cells;

					 /**
					  * Structure keeping the three color
					  * values in the RGB system.
					  */
	struct RgbValues
	{
	    float red;
	    float green;
	    float blue;

					     /**
					      * Return <tt>true</tt> if the
					      * color represented by
					      * the three color values
					      * is a grey scale,
					      * i.e. all components
					      * are equal.
					      */
	    bool is_grey () const;
	};

					 /**
					  * Definition of a function pointer
					  * type taking a value and returning
					  * a triple of color values in RGB
					  * values.
					  *
					  * Besides the actual value by which
					  * the color is to be computed, min
					  * and max values of the data to
					  * be colorized are given as well.
					  */
	typedef RgbValues (*ColorFunction) (const double value,
					    const double min_value,
					    const double max_value);

					 /**
					  * This is a pointer to the function
					  * which is used to colorize the cells.
					  * By default, it points to the
					  * static function <tt>default_color_function</tt>
					  * which is a member of this class.
					  */
	ColorFunction color_function;


					 /**
					  * Default colorization function. This
					  * one does what one usually wants:
					  * It shifts colors from black (lowest
					  * value) through blue, green and red
					  * to white (highest value). For the
					  * exact defition of the color scale
					  * refer to the implementation.
					  *
					  * This function was originally written
					  * by Stefan Nauber.
					  */
	static RgbValues
        default_color_function (const double value,
                                const double min_value,
                                const double max_value);

					 /**
					  * This is an alternative color
					  * function producing a grey scale
					  * between black (lowest values)
					  * and white (highest values). You
					  * may use it by setting the
					  * <tt>color_function</tt> variable to the
					  * address of this function.
					  */
	static RgbValues
        grey_scale_color_function (const double value,
                                   const double min_value,
                                   const double max_value);
	
					 /**
					  * This is one more
					  * alternative color function
					  * producing a grey scale
					  * between white (lowest
					  * values) and black (highest
					  * values), i.e. the scale is
					  * reversed to the previous
					  * one. You may use it by
					  * setting the
					  * <tt>color_function</tt>
					  * variable to the address of
					  * this function.
					  */
	static RgbValues
        reverse_grey_scale_color_function (const double value,
                                           const double min_value,
                                           const double max_value);
	
					 /**
					  * Constructor.
					  */
	EpsFlags (const unsigned int  height_vector = 0,
		  const unsigned int  color_vector  = 0,
		  const SizeType      size_type     = width,
		  const unsigned int  size          = 300,
		  const double        line_width    = 0.5,
		  const double        azimut_angle  = 60,
		  const double        turn_angle    = 30,
		  const double        z_scaling     = 1.0,
		  const bool          draw_mesh     = true,
		  const bool          draw_cells    = true,
		  const bool          shade_cells   = true,
		  const ColorFunction color_function= &default_color_function);

					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  *
					  * For coloring, only the color
					  * functions declared in this
					  * class are offered.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };

    				     /**
				      * Flags describing the details of
				      * output in gmv format. At
				      * present no flags are implemented.
				      */
    struct GmvFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * <tt>private</tt> as well).
					  */
	int dummy;

      public:
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };

    				     /**
				      * Flags describing the details of
				      * output in Tecplot format. At
				      * present no flags are implemented.
				      */
    struct TecplotFlags 
    {

      public:

					 /**
					  * This variable is needed to hold the
					  * output file name when using the
					  * Tecplot API to write binary files.
					  * If the user doesn't set the file
					  * name with this variable only
					  * ASCII Tecplot output will be
					  * produced.
					  */
        const char* tecplot_binary_file_name;
      
	                                 /**
	                                  * Constructor
	                                  **/
	TecplotFlags (const char* tecplot_binary_file_name = NULL);
      
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };

    				     /**
				      * Flags describing the details
				      * of output in VTK format. At
				      * present no flags are
				      * implemented.
				      */
    struct VtkFlags 
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * <tt>private</tt> as well).
					  */
	int dummy;

      public:
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };



    				     /**
				      * Flags describing the details
				      * of output in deal.II
				      * intermediate format. At
				      * present no flags are
				      * implemented.
				      */
    struct Deal_II_IntermediateFlags
    {
      private:
					 /**
					  * Dummy entry to suppress compiler
					  * warnings when copying an empty
					  * structure. Remove this member
					  * when adding the first flag to
					  * this structure (and remove the
					  * <tt>private</tt> as well).
					  */
	int dummy;

      public:
					 /**
					  * Declare all flags with name
					  * and type as offered by this
					  * class, for use in input files.
					  */
	static void declare_parameters (ParameterHandler &prm);

					 /**
					  * Read the parameters declared in
					  * <tt>declare_parameters</tt> and set the
					  * flags for this output format
					  * accordingly.
					  *
					  * The flags thus obtained overwrite
					  * all previous contents of this object.
					  */
	void parse_parameters (ParameterHandler &prm);

					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this
					  * object. Since sometimes
					  * the size of objects can
					  * not be determined exactly
					  * (for example: what is the
					  * memory consumption of an
					  * STL <tt>std::map</tt> type with a
					  * certain number of
					  * elements?), this is only
					  * an estimate. however often
					  * quite close to the true
					  * value.
					  */
	unsigned int memory_consumption () const;
    };
    
				     /**
				      * Provide a data type specifying
				      * the presently supported output
				      * formats.
				      */
    enum OutputFormat {
	  default_format,
					   /**
					    * Output for IBM OpenDX.
					    */
	  dx,
					   /**
					    * Output in AVS UCD format.
					    */
	  ucd,
					   /**
					    * Output for the gnuplot tool.
					    */
	  gnuplot,
					   /**
					    * Output for the povray raytracer.
					    */
	  povray,
					   /**
					    * Output in encapsulated
					    * PostScript.
					    */
	  eps,
					   /**
					    * Output for GMV.
					    */
	  gmv,
					   /**
					    * Output for tecplot in
					    * text format.
					    */
	  
	  tecplot,
					   /**
					    * Output for tecplot in
					    * binaryformat. Faster and
					    * smaller than text
					    * format.
					    */
	  tecplot_binary,
					   /**
					    * Output in VTK format.
					    */
	  vtk,
					   /**
					    * Output in deal.II
					    * intermediate format.
					    */
	  deal_II_intermediate
    };


    				     /**
				      * Write the given list of patches
				      * to the output stream in OpenDX
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_dx (const std::vector<Patch<dim,spacedim> > &patches,
			  const std::vector<std::string>          &data_names,
			  const DXFlags                          &flags,
			  std::ostream                            &out);
    
				     /**
				      * Write the given list of patches
				      * to the output stream in ucd
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_ucd (const std::vector<Patch<dim,spacedim> > &patches,
			   const std::vector<std::string>          &data_names,
			   const UcdFlags                          &flags,
			   std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in gnuplot
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_gnuplot (const std::vector<Patch<dim,spacedim> > &patches,
			       const std::vector<std::string>          &data_names,
			       const GnuplotFlags                      &flags,
			       std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in povray
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_povray (const std::vector<Patch<dim,spacedim> > &patches,
			      const std::vector<std::string>          &data_names,
			      const PovrayFlags                       &flags,
			      std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in eps
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_eps (const std::vector<Patch<dim,spacedim> > &patches,
			   const std::vector<std::string>          &data_names,
			   const EpsFlags                          &flags,
			   std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in gmv
				      * format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_gmv (const std::vector<Patch<dim,spacedim> > &patches,
			   const std::vector<std::string>          &data_names,
			   const GmvFlags                          &flags,
			   std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in Tecplot
				      * ASCII format. See the general
				      * documentation for more information
				      * on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_tecplot (const std::vector<Patch<dim,spacedim> > &patches,
			       const std::vector<std::string>          &data_names,
			       const TecplotFlags                      &flags,
			       std::ostream                            &out);

    				     /**
				      * Write the given list of patches
				      * to the output stream in Tecplot
				      * binary format. See the general
				      * documentation for more information
				      * on the parameters. <tt>tecplot_binary_file_name</tt>
				      * (specified through the TecplotFlags
				      * struct) indicates the name of the file
				      * to be written.  If the file name is not
				      * set ASCII output is produced.
				      *
				      * If the Tecplot API is not present this simply
				      * calls the standard write_tecplot file so that
				      * ASCII output is still produced.
				      */
    template <int dim, int spacedim>
    static void write_tecplot_binary (const std::vector<Patch<dim,spacedim> > &patches,
				      const std::vector<std::string>          &data_names,
				      const TecplotFlags                      &flags,
				      std::ostream                            &out);

    				     /**
				      * Write the given list of
				      * patches to the output stream
				      * in vtk format. See the general
				      * documentation for more
				      * information on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_vtk (const std::vector<Patch<dim,spacedim> > &patches,
			   const std::vector<std::string>          &data_names,
			   const VtkFlags                          &flags,
			   std::ostream                            &out);

    				     /**
				      * Write the given list of
				      * patches to the output stream
				      * in deal.II intermediate
				      * format. See the general
				      * documentation for more
				      * information on the parameters.
				      */
    template <int dim, int spacedim>
    static void write_deal_II_intermediate (const std::vector<Patch<dim,spacedim> > &patches,
					    const std::vector<std::string>          &data_names,
					    const Deal_II_IntermediateFlags         &flags,
					    std::ostream                            &out);


                                     /**
                                      * Given an input stream that contains
                                      * data written by
                                      * write_deal_II_intermediate, determine
                                      * the <tt>dim</tt> and <tt>spacedim</tt>
                                      * template parameters with which that
                                      * function was called, and return them
                                      * as a pair of values.
                                      *
                                      * Note that this function eats a number
                                      * of elements at the present position of
                                      * the stream, and therefore alters
                                      * it. In order to read from it using,
                                      * for example, the DataOutReader class,
                                      * you may wish to either reset the
                                      * stream to its previous position, or
                                      * close and reopen it.
                                      */
    static
    std::pair<unsigned int, unsigned int>
    determine_intermediate_format_dimensions (std::istream &input);
    
				     /**
				      * Return the <tt>OutputFormat</tt>
				      * value corresponding to the
				      * given string. If the string
				      * does not match any known
				      * format, an exception is
				      * thrown.
				      *
				      * Since this function does not
				      * need data from this object, it
				      * is static and can thus be
				      * called without creating an
				      * object of this class. Its main
				      * purpose is to allow a program
				      * to use any implemented output
				      * format without the need to
				      * extend the program's parser
				      * each time a new format is
				      * implemented.
				      *
				      * To get a list of presently
				      * available format names,
				      * e.g. to give it to the
				      * ParameterHandler class,
				      * use the function
				      * get_output_format_names().
				      */
    static OutputFormat parse_output_format (const std::string &format_name);

				     /**
				      * Return a list of implemented
				      * output formats. The different
				      * names are separated by
				      * vertical bar signs (<tt>`|'</tt>)
				      * as used by the
				      * ParameterHandler classes.
				      */
    static std::string get_output_format_names ();

				     /**
				      * Provide a function which tells us which
				      * suffix a file with a given output format
				      * usually has. At present the following
				      * formats are defined:
				      * <ul>
				      * <li> <tt>dx</tt>: <tt>.dx</tt>
				      * <li> <tt>ucd</tt>: <tt>.inp</tt>
				      * <li> <tt>gnuplot</tt>: <tt>.gnuplot</tt>
				      * <li> <tt>povray</tt>: <tt>.pov</tt>
				      * <li> <tt>eps</tt>: <tt>.eps</tt>
				      * <li> <tt>gmv</tt>: <tt>.gmv</tt>
				      * <li> <tt>vtk</tt>: <tt>.vtk</tt>.
				      * <li> <tt>deal_II_intermediate</tt>: <tt>.d2</tt>.
				      * </ul>
				      */
    static std::string default_suffix (const OutputFormat output_format);
    
				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    static unsigned int memory_consumption ();
    
				     /** @addtogroup Exceptions
				      * @{ */
    
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidDatasetSize,
		    int, int,
		    << "The number of points in this data set is " << arg1
		    << ", but we expected " << arg2 << " in each space direction.");
				     /**
				      * An output function did not
				      * receive any patches for
				      * writing.
				      */
    DeclException0 (ExcNoPatches);
				     /**
				       * Exception
				       */
    DeclException0 (ExcTecplotAPIError);
 				      /**
				       * Exception
				       */
    DeclException1 (ExcErrorOpeningTecplotFile,
		    char*,
		    << "There was an error opening Tecplot file " << arg1
		    << " for output");
    
				     //@}  
  private:
				     /**
				      * Class holding the data of one
				      * cell of a patch in two space
				      * dimensions for output. It is
				      * the projection of a cell in
				      * three-dimensional space (two
				      * coordinates, one height value)
				      * to the direction of sight.
				      */
    class EpsCell2d
    {
      public:
	
					 /**
					  * Vector of vertices of this cell.
					  */
	Point<2> vertices[4];
	
					 /**
					  * Data value from which the actual
					  * colors will be computed by
					  * the colorization function stated
					  * in the <tt>EpsFlags</tt> class.
					  */
	float color_value;
	
					 /**
					  * Depth into the picture, which
					  * is defined as the distance from
					  * an observer at an the origin in
					  * direction of the line of sight.
					  */
	float depth;
	
					 /**
					  * Comparison operator for
					  * sorting.
					  */
	bool operator < (const EpsCell2d &) const;
    };


				     /**
				      * This is a helper function for
				      * the <tt>write_gmv</tt>
				      * function. There, the data in
				      * the patches needs to be copied
				      * around as output is one
				      * variable globally at a time,
				      * rather than all data on each
				      * vertex at a time. This copying
				      * around can be done detached
				      * from the main thread, and is
				      * thus moved into this separate
				      * function.
				      *
				      * Note that because of the
				      * similarity of the formats,
				      * this function is also used by
				      * the Vtk and Tecplot output
				      * functions.
				      */
    template <int dim, int spacedim>
    static void
    write_gmv_reorder_data_vectors (const std::vector<Patch<dim,spacedim> > &patches,
				    Table<2,double>                         &data_vectors);

};




/**
 * This class is the interface to the <tt>DataOutBase</tt> class, as already its name
 * might suggest. It does not offer much functionality apart from a way to
 * access the implemented formats and a way to dynamically dispatch what output
 * format to chose.
 *
 * This class is thought as a base class to classes actually
 * generating data for output. It has two abstract virtual functions,
 * <tt>get_patches</tt> and <tt>get_dataset_names</tt> which are to produce the data
 * which is actually needed. These are the only functions that need to
 * be overloaded by a derived class.  In additional to that, it has a
 * function for each output format supported by the underlying base
 * class which gets the output data using these two virtual functions
 * and passes them to the raw output functions.
 *
 * The purpose of this class is mainly two-fold: to support storing flags
 * by which the output in the different output formats are controlled,
 * and means to work with output in a way where output format, flags and
 * other things are determined at run time. In addition to that it offers
 * the abstract interface to derived classes briefly discussed above.
 *
 *
 * @section DataOutInterfaceOF Output flags
 *
 * The way we treat flags in this class is very similar to that used in
 * the <tt>GridOut</tt> class. For detailed information on the why's and how's,
 * as well as an example of programming, we refer to the documentation
 * of that class.
 *
 * Basically, this class stores a set of flags for each output format
 * supported by the underlying <tt>DataOutBase</tt> class. These are used
 * whenever one of the <tt>write_*</tt> functions is used. By default, the
 * values of these flags are set to reasonable start-ups, but in case
 * you want to change them, you can create a structure holding the flags
 * for one of the output formats and set it using the <tt>set_flags</tt> functions
 * of this class to determine all future output the object might produce
 * by that output format.
 *
 * For information on what parameters are supported by different output
 * functions, please see the documentation of the <tt>DataOutBase</tt> class and
 * its member classes.
 *
 *
 * @section DataOutInterfaceSelectP Run time selection of output parameters
 *
 * In the output flags classes, described above, many flags are
 * defined for output in the different formats. In order to make them
 * available to the input file handler class <tt>ParameterHandler</tt>, each
 * of these has a function declaring these flags to the parameter
 * handler and to read them back from an actual input file. In order
 * to avoid that in user programs these functions have to be called
 * for each available output format and the respective flag class, the
 * present <tt>DataOutInterface</tt> class offers a function
 * <tt>declare_parameters</tt> which calls the respective function of all
 * known output format flags classes. The flags of each such format
 * are packed together in a subsection in the input file.
 * Likewise, there is a function <tt>parse_parameters</tt> which reads
 * these parameters and stores them in the flags associated with
 * this object (see above).
 * 
 * Using these functions, you do not have to track which formats are
 * presently implemented.
 *
 * Usage is as follows:
 * @code
 *                               // within function declaring parameters:
 *   ...
 *   prm.enter_subsection ("Output format options");
 *     DataOutInterface<dim>::declare_parameters (prm);
 *   prm.leave_subsection ();
 *   ...
 *
 *
 *                               // within function doing the output:
 *   ...
 *   DataOut<dim> out;
 *   prm.enter_subsection ("Output format options");
 *   out.parse_parameters (prm);
 *   prm.leave_subsection ();
 *   ...
 * @endcode
 * Note that in the present example, the class <tt>DataOut</tt> was used. However, any
 * other class derived from <tt>DataOutInterface</tt> would work alike.
 *
 *
 * @section DataOutInterfaceSelectF Run time selection of formats
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
 * The functions offering the different output format names are,
 * respectively, <tt>default_suffix</tt>, <tt>parse_output_format</tt>, and
 * <tt>get_output_format_names</tt>. They make the selection of ouput formats
 * in parameter files much easier, and especially independent of
 * the formats presently implemented. User programs need therefore not
 * be changed whenever a new format is implemented.
 *
 * Additionally, objects of this class have a default format, which
 * can be set by the parameter "Output format" of the parameter
 * file. Within a program, this can be changed by the member function
 * <tt>set_default_format</tt>. Using this default format, it is possible to leave
 * the format selection completely to the parameter file. A suitable
 * suffix for the output file name can be obtained by <tt>default_suffix</tt>
 * without arguments.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim, int spacedim=dim>
class DataOutInterface : private DataOutBase
{
  public:
                                     /**
                                      * Import a few names that were
                                      * previously in this class and have then
                                      * moved to the base class. Since the
                                      * base class is inherited from
                                      * privately, we need to re-import these
                                      * symbols to make sure that references
                                      * to DataOutInterface<dim,spacedim>::XXX
                                      * remain valid.
                                      */
    using DataOutBase::OutputFormat;
    using DataOutBase::default_format;
    using DataOutBase::dx;
    using DataOutBase::gnuplot;
    using DataOutBase::povray;
    using DataOutBase::eps;
    using DataOutBase::tecplot;
    using DataOutBase::tecplot_binary;
    using DataOutBase::vtk;
    using DataOutBase::deal_II_intermediate;
    using DataOutBase::parse_output_format;
    using DataOutBase::get_output_format_names;
    
                                     /**
                                      * Destructor. Does nothing, but is
                                      * declared virtual since this class has
                                      * virtual functions.
                                      */
    virtual ~DataOutInterface ();

				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in OpenDX
				      * format.
				      */
    void write_dx (std::ostream &out) const;

				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in UCD
				      * format.
				      */
    void write_ucd (std::ostream &out) const;

				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in GNUPLOT
				      * format.
				      */
    void write_gnuplot (std::ostream &out) const;

    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in POVRAY
				      * format.
				      */
    void write_povray (std::ostream &out) const;

				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in EPS
				      * format.
				      */
    void write_eps (std::ostream &out) const;

    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in GMV
				      * format.
				      */
    void write_gmv (std::ostream &out) const;


    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in Tecplot
				      * format.
				      */
    void write_tecplot (std::ostream &out) const;

    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it in the Tecplot binary
				      * output format. Note that the name
				      * of the output file must be specified
				      * through the TecplotFlags interface.
				      */
    void write_tecplot_binary (std::ostream &out) const;

    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function and
				      * write it to <tt>out</tt> in Vtk
				      * format.
				      */
    void write_vtk (std::ostream &out) const;

    				     /**
				      * Obtain data through the
				      * <tt>get_patches</tt> function
				      * and write it to <tt>out</tt>
				      * in deal.II intermediate
				      * format.
				      */
    void write_deal_II_intermediate (std::ostream &out) const;
    
				     /**
				      * Write data and grid to <tt>out</tt>
				      * according to the given data
				      * format. This function simply
				      * calls the appropriate
				      * <tt>write_*</tt> function. If no
				      * output format is requested,
				      * the <tt>default_format</tt> is
				      * written.
				      *
				      * An error occurs if no format
				      * is provided and the default
				      * format is <tt>default_format</tt>.
				      */
    void write (std::ostream       &out,
		const OutputFormat  output_format = default_format) const;

				     /**
				      * Set the default format. The
				      * value set here is used
				      * anytime, output for format
				      * <tt>default_format</tt> is
				      * requested.
				      */
    void set_default_format (const OutputFormat default_format);

				     /**
				      * Set the flags to be used for
				      * output in OpenDX format.
				      */
    void set_flags (const DXFlags &dx_flags);

				     /**
				      * Set the flags to be used for
				      * output in UCD format.
				      */
    void set_flags (const UcdFlags &ucd_flags);

    				     /**
				      * Set the flags to be used for
				      * output in GNUPLOT format.
				      */
    void set_flags (const GnuplotFlags &gnuplot_flags);

    				     /**
				      * Set the flags to be used for
				      * output in POVRAY format.
				      */
    void set_flags (const PovrayFlags &povray_flags);

    				     /**
				      * Set the flags to be used for
				      * output in EPS output.
				      */
    void set_flags (const EpsFlags &eps_flags);

    				     /**
				      * Set the flags to be used for
				      * output in GMV format.
				      */
    void set_flags (const GmvFlags &gmv_flags);

    				     /**
				      * Set the flags to be used for
				      * output in Tecplot format.
				      */
    void set_flags (const TecplotFlags &tecplot_flags);

    				     /**
				      * Set the flags to be used for
				      * output in VTK format.
				      */
    void set_flags (const VtkFlags &vtk_flags);

    				     /**
				      * Set the flags to be used for
				      * output in VTK format.
				      */
    void set_flags (const Deal_II_IntermediateFlags &deal_II_intermediate_flags);    

				     /**
				      * A function that returns the same
				      * string as the respective function in
				      * the base class does; the only
				      * exception being that if the parameter
				      * is omitted, then the value for the
				      * present default format is returned.
				      */
    std::string
    default_suffix (const OutputFormat output_format = default_format) const;

				     /**
				      * Declare parameters for all
				      * output formats by declaring
				      * subsections within the
				      * parameter file for each output
				      * format and call the respective
				      * <tt>declare_parameters</tt>
				      * functions of the flag classes
				      * for each output format.
				      *
				      * Some of the declared
				      * subsections may not contain
				      * entries, if the respective
				      * format does not export any
				      * flags.
				      *
				      * Note that the top-level
				      * parameters denoting the number
				      * of subdivisions per patch and
				      * the output format are not
				      * declared, since they are only
				      * passed to virtual functions
				      * and are not stored inside
				      * objects of this type. You have
				      * to declare them yourself.
				      */
    static void declare_parameters (ParameterHandler &prm);

				     /**
				      * Read the parameters declared
				      * in <tt>declare_parameters</tt> and
				      * set the flags for the output
				      * formats accordingly.
				      *
				      * The flags thus obtained
				      * overwrite all previous
				      * contents of the flag objects
				      * as default-constructed or set
				      * by the set_flags() function.
				      */
    void parse_parameters (ParameterHandler &prm);
    
				     /**
				      * Determine an estimate for
				      * the memory consumption (in
				      * bytes) of this
				      * object. Since sometimes
				      * the size of objects can
				      * not be determined exactly
				      * (for example: what is the
				      * memory consumption of an
				      * STL <tt>std::map</tt> type with a
				      * certain number of
				      * elements?), this is only
				      * an estimate. however often
				      * quite close to the true
				      * value.
				      */
    unsigned int memory_consumption () const;

  protected:
				     /**
				      * This is the abstract function
				      * through which derived classes
				      * propagate preprocessed data in
				      * the form of Patch
				      * structures (declared in the
				      * base class DataOutBase) to
				      * the actual output
				      * function. You need to overload
				      * this function to allow the
				      * output functions to know what
				      * they shall print.
				      */
    virtual const std::vector<typename DataOutBase::Patch<dim,spacedim> > &
    get_patches () const = 0;

				     /**
				      * Abstract virtual function
				      * through which the names of
				      * data sets are obtained by the
				      * output functions of the base
				      * class.
				      */
    virtual std::vector<std::string> get_dataset_names () const = 0;

  private:
				     /**
				      * Standard output format.  Use
				      * this format, if output format
				      * default_format is
				      * requested. It can be changed
				      * by the <tt>set_format</tt> function
				      * or in a parameter file.
				      */
    OutputFormat default_fmt;
    
				     /**
				      * Flags to be used upon output
				      * of OpenDX data. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    DXFlags     dx_flags;

				     /**
				      * Flags to be used upon output
				      * of UCD data. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    UcdFlags     ucd_flags;

				     /**
				      * Flags to be used upon output
				      * of GNUPLOT data. Can be
				      * changed by using the
				      * <tt>set_flags</tt> function.
				      */
    GnuplotFlags gnuplot_flags;

    				     /**
				      * Flags to be used upon output
				      * of POVRAY data. Can be changed
				      * by using the <tt>set_flags</tt>
				      * function.
				      */
    PovrayFlags povray_flags;

				     /**
				      * Flags to be used upon output
				      * of EPS data in one space
				      * dimension. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    EpsFlags     eps_flags;    

				     /**
				      * Flags to be used upon output
				      * of gmv data in one space
				      * dimension. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    GmvFlags     gmv_flags;

				     /**
				      * Flags to be used upon output
				      * of Tecplot data in one space
				      * dimension. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    TecplotFlags tecplot_flags;

				     /**
				      * Flags to be used upon output
				      * of vtk data in one space
				      * dimension. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    VtkFlags     vtk_flags;

				     /**
				      * Flags to be used upon output
				      * of vtk data in one space
				      * dimension. Can be changed by
				      * using the <tt>set_flags</tt>
				      * function.
				      */
    Deal_II_IntermediateFlags     deal_II_intermediate_flags;
};



/**
 * A class that is used to read data written in deal.II intermediate
 * format back in, so that it can be written out in any of the other
 * supported graphics formats. This class has two main purposes:
 *
 * The first use of this class is so that application programs can
 * defer the decision of which graphics format to use until after the
 * program has been run. The data is written in intermediate format
 * into a file, and later on it can then be converted into any
 * graphics format you wish. This may be useful, for example, if you
 * want to convert it to gnuplot format to get a quick glimpse and
 * later on want to convert it to OpenDX format as well to get a high
 * quality version of the data. The present class allows to read this
 * intermediate format back into the program, and allows it to be
 * written in any other supported format using the relevant functions
 * of the base class.
 *
 * The second use is mostly useful in parallel programs: rather than
 * having one central process generate the graphical output for the
 * entire program, one can let each process generate the graphical
 * data for the cells it owns, and write it into a separate file in
 * intermediate format. Later on, all these intermediate files can
 * then be read back in and merged together, a process that is fast
 * compared to generating the data in the first place. The use of the
 * intermediate format is mostly because it allows separate files to
 * be merged, while this is almost impossible once the data has been
 * written out in any of the supported established graphics formats.
 *
 * This second use scenario is explained in some detail in the step-18
 * example program.
 *
 * Both these applications are implemented in the step-19 example program.
 * There, a slight complication is also explained: in order to read data back
 * into this object, you have to know the template parameters for the space
 * dimension which were used when writing the data. If this knowledge is
 * available at compile time, then this is no problem. However, if it is not
 * (such as in a simple format converted), then it needs to be figured out at
 * run time, even though the compiler already needs it at compile time. A way
 * around using the DataOutBase::determine_intermediate_format_dimensions()
 * function is explained in step-19.
 * 
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, int spacedim=dim>
class DataOutReader : public DataOutInterface<dim,spacedim>
{
  public:
				     /**
				      * Read a sequence of patches as
				      * written previously by
				      * <tt>DataOutBase::write_deal_II_intermediate</tt>
				      * and store them in the present
				      * object. This overwrites any
				      * previous content.
				      */
    void read (std::istream &in);

                                     /**
                                      * This function can be used to
                                      * merge the patches read by the
                                      * other object into the patches
                                      * that this present object
                                      * stores. This is sometimes
                                      * handy if one has, for example,
                                      * a domain decomposition
                                      * algorithm where each block is
                                      * represented by a DoFHandler of
                                      * its own, but one wants to
                                      * output the solution on all the
                                      * blocks at the same
                                      * time. Alternatively, it may
                                      * also be used for parallel
                                      * programs, where each process
                                      * only generates output for its
                                      * share of the cells, even if
                                      * all processes can see all
                                      * cells.
                                      *
                                      * For this to work, the input
                                      * files for the present object
                                      * and the given argument need to
                                      * have the same number of output
                                      * vectors, and they need to use
                                      * the same number of
                                      * subdivisions per patch. The
                                      * output will probably look
                                      * rather funny if patches in
                                      * both objects overlap in space.
                                      *
                                      * If you call read() for this
                                      * object after merging in
                                      * patches, the previous state is
                                      * overwritten, and the merged-in
                                      * patches are lost.
				      *
                                      * This function will fail if
                                      * either this or the other
                                      * object did not yet set up any
                                      * patches.
				      *
				      * The use of this function is
				      * demonstrated in step-19.
                                      */
    void merge (const DataOutReader<dim,spacedim> &other);
    
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcNoPatches);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcIncompatibleDatasetNames);
                                     /**
                                      * Exception
                                      */
    DeclException0 (ExcIncompatiblePatchLists);

  protected:
				     /**
				      * This is the function
				      * through which this class
				      * propagates preprocessed data in
				      * the form of Patch
				      * structures (declared in the
				      * base class DataOutBase) to
				      * the actual output
				      * function.
				      *
				      * It returns the patches as read
				      * the last time a stream was
				      * given to the read() function.
				      */
    virtual const std::vector<typename ::DataOutBase::Patch<dim,spacedim> > &
    get_patches () const;

				     /**
				      * Abstract virtual function
				      * through which the names of
				      * data sets are obtained by the
				      * output functions of the base
				      * class.
				      *
				      * Return the names of the
				      * variables as read the last
				      * time we read a file.
				      */
    virtual std::vector<std::string> get_dataset_names () const;
    
  private:
    std::vector<typename ::DataOutBase::Patch<dim,spacedim> > patches;
    std::vector<std::string> dataset_names;
};




/* -------------------- inline functions ------------------- */

inline
bool
DataOutBase::EpsFlags::RgbValues::is_grey () const
{
  return (red == green) && (red == blue);
}


/* -------------------- template functions ------------------- */

/**
 * Output operator for an object of type
 * <tt>DataOutBase::Patch</tt>. This operator dumps the intermediate
 * graphics format represented by the patch data structure. It may
 * later be converted into regular formats for a number of graphics
 * programs.
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, int spacedim>
std::ostream &
operator << (std::ostream                           &out,
	     const DataOutBase::Patch<dim,spacedim> &patch);



/**
 * Input operator for an object of type
 * <tt>DataOutBase::Patch</tt>. This operator reads the intermediate
 * graphics format represented by the patch data structure, using the
 * format in which it was written using the operator<<.
 *
 * @author Wolfgang Bangerth, 2005
 */
template <int dim, int spacedim>
std::istream &
operator >> (std::istream                     &in,
	     DataOutBase::Patch<dim,spacedim> &patch);



#endif
