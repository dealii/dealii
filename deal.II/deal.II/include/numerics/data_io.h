/*----------------------------   data_io.h     ---------------------------*/
/*      <Id:>                 */
#ifndef __data_io_H
#define __data_io_H
/*----------------------------   data_io.h     ---------------------------*/

#include <base/exceptions.h>
#include <grid/geometry_info.h>
#include <vector>
#include <string>

template <int dim> class Triangulation;
template <int dim> class DoFHandler;

class dVector;






/**
 *  Structure which is passed to the #Triangulation<dim>::create_triangulation#
 *  function. It contains all data needed to construct a cell, namely the
 *  indices of the vertices and the material indicator.
 */
template <int dim>
struct CellData {
    int           vertices[GeometryInfo<dim>::vertices_per_cell];
    unsigned char material_id;
};





/**
 *  Structure to be passed to the #Triangulation<dim>::create_triangulation#
 *  function to describe boundary information.
 *
 *  This structure is the same for all dimensions, since we use an input
 *  function which is the same for all dimensions. The content of objects
 *  of this structure varies with the dimensions, however.
 *
 *  Since in one space dimension, there is no boundary information apart
 *  from the two end points of the interval, this structure does not contain
 *  anything and exists only for consistency, to allow a common interface
 *  for all space dimensions. All fields should always be empty.
 *
 *  Boundary data in 2D consists
 *  of a list of lines which belong to a given boundary component. A
 *  boundary component is a list of lines which are given a common
 *  number describing the boundary condition to hold on this part of the
 *  boundary. The triangulation creation function gives lines not in this
 *  list either the boundary indicator zero (if on the boundary) or 255
 *  (if in the interior). Explicitely giving a line the indicator 255
 *  will result in an error, as well as giving an interior line a boundary
 *  indicator.
 */
struct SubCellData {
				     /**
				      * Each record of this vector describes
				      * a line on the boundary and its boundary
				      * indicator.
				      */
    vector<CellData<1> > boundary_lines;

				     /**
				      * Each record of this vector describes
				      * a quad on the boundary and its boundary
				      * indicator.
				      */
    vector<CellData<2> > boundary_quads;

				     /**
				      * This function checks whether the vectors
				      * which may not be used in a given
				      * dimension are really empty. I.e.,
				      * whether the #boundary_*# arrays are
				      * empty when in one space dimension
				      * and whether the #boundary_quads#
				      * array is empty when in two dimensions.
				      *
				      * Since this structure is the same for all
				      * dimensions, the actual dimension has
				      * to be given as a parameter.
				      */
    bool check_consistency (const unsigned int dim) const;
};


    


/**
 * This class implements an input mechanism for grid data. It allows to
 * read a grid structure into a triangulation object. Future versions
 * will also allow to read data on this grid into vectors.
 *
 * At present, only UCD (unstructured cell data) is supported as input
 * format for grid data. Any numerical data after the block of topological
 * information is ignored.
 *
 * To read grid data, the triangulation to be fed with has to be empty.
 * When giving a file which does not contain the assumed information or
 * which does not keep to the right format, the state of the triangulation
 * will be undefined afterwards. Upon input, only lines in one dimension
 * and line and quads in two dimensions are accepted. All other cell types
 * (e.g. triangles in two dimensions, quads and hexes in 3d) are rejected.
 * The vertex and cell numbering in the UCD file, which
 * need not be consecutively, is lost upon transfer to the triangulation
 * object, since this one needs consecutively numbered elements.
 *
 * Material indicators are accepted to denote the material id of cells and
 * to denote boundary part indication for lines in 2D. Read the according
 * sections in the documentation of the \Ref{Triangulation} class for
 * further details.
 *
 *
 * \subsection{Structure of input grid data}
 * 
 * It is your duty to use a correct numbering of vertices in the cell list,
 * i.e. for lines, you have to first give the vertex with the lower coordinate
 * value, then that with the higher coordinate value. For quadrilaterals in
 * two dimensions, the vertex indices in the #quad# list have to be such that
 * the vertices are numbered in counter-clockwise sense.
 *
 * In two dimensions, another difficulty occurs, which has to do with the sense
 * of a quadrilateral. A quad consists of four lines which have a direction,
 * which is per definitionem as follows:
 * \begin{verbatim}
 *   3-->--2
 *   |     |
 *   ^     ^
 *   |     |
 *   0-->--1
 * \end{verbatim}
 * Now, two adjacent cells must have a vertex numbering such that the direction
 * of the common side is the same. For example, the following two quads
 * \begin{verbatim}
 *   3---4---5
 *   |   |   |
 *   0---1---2
 * \end{verbatim}
 * may be characterised by the vertex numbers (0 1 4 3) and (1 2 5 4), since
 * the middle line would get the direction #1->4# when viewed from both cells.
 * The numbering (0 1 4 3) and (5 4 1 2) would not be allowed, since the left
 * quad would give the common line the direction #1->4#, while the right one
 * would want to use #4->1#, leading to ambiguity. The #Triangulation# object
 * is capable of detecting this special case, which can be eliminated by
 * rotating the indices of the right quad by two. However, it would not
 * know what to do if you gave the vertex indices (4 1 2 5), since then it
 * would have to rotate by one element or three, the decision which to take is
 * not yet implemented.
 *
 * There are more ambiguous cases, where the triangulation may not know what
 * to do at all without the use of very sophisticated algorithms. On such example
 * is the following:
 * \begin{verbatim}
 *   9---10-----11
 *   |   |    / |
 *   6---7---8  |
 *   |   |   |  |
 *   3---4---5  |
 *   |   |    \ |
 *   0---1------2
 * \end{verbatim}
 * Assume that you had numbered the vertices in the cells at the left boundary
 * in a way, that the following line directions are induced:
 * \begin{verbatim}
 *   9->-10-----11
 *   ^   ^    / |
 *   6->-7---8  |
 *   ^   ^   |  |
 *   3->-4---5  |
 *   ^   ^    \ |
 *   0->-1------2
 * \end{verbatim}
 * (This could for example be done by using the indices (0 1 4 3), (3 4 7 6),
 * (6 7 10 9) for the three cells). Now, you will not find a way of giving
 * indices for the right cells, without introducing either ambiguity for
 * one line or other, or without violating that within each cells, there must be
 * one vertex from which both lines are directed away and the opposite one to
 * which both adjacent lines point to.
 *
 * The solution in this case is to renumber one of the three left cells, e.g.
 * by reverting the sense of the line between vertices 7 and 10 by numbering
 * the top left cell by (9 6 7 10).
 *
 * But this is a thing that the triangulation
 * object can't do for you, since it would involve backtracking to cells
 * already created when we find that we can't number the indices of one of
 * the rightmost cells consistently. It is neither clear how to do this
 * backtracking nor whether it can be done with a stopping algorithm, if
 * possible within polynomial time. This kind of numbering must be made
 * upon construction of the coarse grid, unfortunately.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>
class DataIn {
  public:
				     /**
				      * Constructor.
				      */
    DataIn ();
    
				     /**
				      * Attach this triangulation
				      * to be fed with the grid data.
				      */
    void attach_triangulation (Triangulation<dim> *tria);

				     /**
				      * Read grid data from an ucd file.
				      * Numerical data is ignored.
				      */
    void read_ucd (istream &);

				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoTriangulationSelected);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVertexIndex,
		    int, int,
		    << "Trying to access invalid vertex index " << arg2
		    << " while creating cell " << arg1);
				     /**
				      * Exception
				      */
    DeclException0 (ExcInternalError);
    
  private:
				     /**
				      * Store address of the triangulation to
				      * be fed with the data read in.
				      */
    Triangulation<dim> *tria;
};







/**
 * This class implements an output mechanism for grid and simulation data
 * in several formats.
 * At present it supports output in UCD (unstructured cell data) and
 * partly in GNUPLOT format.
 *
 * It allows the user to attach a degree of freedom handler object
 * (#DoFHandler#) which also gives access to the geometry data of the
 * underlying triangulation and to add data vectors of which the values
 * are to be written.
 *
 *
 * \subsection{Limitations}
 * 
 * At present, no grouping of components to vectors is implemented, i.e.
 * you can only write each component independent of the others. Also, it
 * is not possible to output calculations which were performed on elements
 * with more or less than one degree of freedom per vertex.
 *
 * 
 * \subsection{UCD format}
 *
 * The UCD format is described in the AVS developer's guide. Due to
 * limitations in the present format, only node based data can be output,
 * so higher order elements are only written with their node values, no
 * interior or line values are used. No use is made of the possibility
 * to give cell and model data since these are not supported by all
 * UCD aware programs.
 *
 * The ASCII UCD format is used. In future versions, a binary version may
 * follow up.
 *
 * Note that to enumerate the vertices, not the vertex index is used but
 * the index of the degree of freedom located on this vertex. This makes
 * the mapping between the vertices and the entries in the data vectors
 * much easier.
 *
 *
 * \subsection{GNUPLOT format}
 *
 * The GNUPLOT format is not able to handle data on unstructured grids
 * directly. Directly would mean that you only give the vertices and
 * the solution values thereon and the program constructs its own grid
 * to represent the data. This is only possible for a structured tensor
 * product grid in two dimensions.
 *
 * In one dimension, the format is obviously #x v1 v2 ...#, where #x#
 * is the coordinate value of a grid point, while the #vi# are the
 * vector elements referring to the present node. Within GNUPLOT,
 * call #plot "filename" using 1:x#. #x# denotes the number of the data set you
 * want to see plus one. For example #using 1:4# would mean to plot the
 * third data vector.
 *
 * For more than one dimension, the #DataOut<dim>::write_gnuplot()# somehow
 * duplicates the functionality of the #Triangulation<dim>::print_gnuplot()#
 * functions. These, however, offer more functionality in some respect.
 * The grid is represented as a sequence of lines, where each cell is
 * a sequence of five vertices (the first one is appended to close the
 * contour of the cell) with the data appended after each vertex. Each cell
 * is therefore a sequence of five lines #x y v1 v2 ...# forming together
 * the bounding line of this cell. After each cell, two newlines are inserted
 * to prevent GNUPLOT from joining the lines bounding two cells.
 *
 * To view the results in two dimensions, use #set data style lines#
 * within gnuplot and call #plot "filename"# to see the grid. Use
 * #set parametric# and #splot "filename" using 1:2:x# to get a 3d surface
 * plot of the (#x-2#)th data set. For example, using #x=4# would mean to
 * plot the second data set.
 *
 * @author Wolfgang Bangerth, 1998
 */
template <int dim>  
class DataOut {
  public:
				     /**
				      * Provide a data type specifying the
				      * presently supported output formats.
				      */
    enum OutputFormat { ucd, gnuplot };
    
				     /**
				      * Constructor
				      */
    DataOut ();
    
				     /**
				      * Designate a dof handler to be used
				      * to extract geometry data and the
				      * mapping between nodes and node values.
				      */
    void attach_dof_handler (DoFHandler<dim> &);

				     /**
				      * Add a data vector together with its
				      * name and the physical unit
				      * (for example meter, kelvin, etc). By
				      * default, "<dimensionless>" is assumed
				      * for the units.
				      *
				      * A pointer to the vector is stored, so
				      * you have to make sure the vector
				      * exists at that address at least as
				      * long as you call the
				      * #write_*# functions.
				      *
				      * It is assumed that the vector has the
				      * same number of components as there are
				      * degrees of freedom in the dof handler.
				      * Therefore, no block vectors are allowed
				      * at present.
				      */
    void add_data_vector (const dVector &data,
			  const string  &name,
			  const string  &units="<dimensionless>");

				     /**
				      * Release the pointers to the data
				      * vectors. You have to set all data
				      * entries again using the
				      * #add_data_vector# function. The pointer
				      * to the dof handler remains stored,
				      * however.
				      */
    void clear_data_vectors ();

				     /**
				      * Write the stored data to the given
				      * stream in UCD data. You may have
				      * written any comment to that stream
				      * before calling this function. Comments
				      * start with the \# character in the
				      * first column of a line and may only
				      * appear at the beginning of a file,
				      * without non-comment lines inbetween.
				      */
    void write_ucd (ostream &out) const;

				     /**
				      * Write data and grid in GNUPLOT format.
				      */
    void write_gnuplot (ostream &out) const;

				     /**
				      * Write data and grid to #out# according
				      * to the given data format. This function
				      * simply calles the appropriate
				      * #write_*# function.
				      */
    void write (ostream &out, const OutputFormat output_format) const;
    
				     /**
				      * Provide a function which tells us which
				      * suffix with a given output format
				      * usually has. At present the following
				      * formats are defined:
				      * \begin{itemize}
				      * \item UCD: #.inp#
				      * \item GNUPLOT: #.gnuplot#
				      * \end{itemize}
				      *
				      * Since this function does not need data
				      * from this object, it is static and can
				      * thus be called without creating an
				      * object of this class.
				      */
    static string default_suffix (const OutputFormat output_format);
    
				     /**
				      * Exception
				      */
    DeclException0 (ExcIncorrectDofsPerVertex);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNotImplemented);
				     /**
				      * Exception
				      */
    DeclException0 (ExcNoDoFHandlerSelected);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidVectorSize,
		    int, int,
		    << "The vector has size " << arg1
		    << " but the DoFHandler objects says there are " << arg2
		    << " degrees of freedom.");
    
  private:

				     /**
				      * Declare an entry in the list of
				      * data elements.
				      */
    struct DataEntry {
					 /**
					  * Default constructor, constructs
					  * invalid object, but is needed for
					  * the STL classes.
					  */
	DataEntry ();
			
					 /**
					  * Constructor
					  */
	DataEntry (const dVector *data, const string name, const string units);
	
					 /**
					  * Pointer to the data vector.
					  */
	const dVector *data;
					 /**
					  * Name of this component.
					  */
	string   name;
					 /**
					  * Physical unit name of this
					  * component.
					  */
	string   units;
    };

				     /**
				      * Pointer to the dof handler object.
				      */
    DoFHandler<dim>   *dofs;

				     /**
				      * List of data elements.
				      */
    vector<DataEntry>  data;


				     /**
				      * Return the number of faces in the
				      * triangulation which have a boundary
				      * indicator not equal to zero. Only
				      * these faces are explicitely printed
				      * in the #write_*# functions;
				      * all faces with indicator 255 are
				      * interior ones and an indicator with
				      * value zero for faces at the boundary
				      * are considered default.
				      *
				      * This function always returns an empty
				      * list in one dimension.
				      *
				      * The reason for this function is the
				      * same as for #write_ucd_faces#. See
				      * there for more information.
				      */
    unsigned int n_boundary_faces () const;

				     /**
				      * Write the grid information about
				      * faces to #out#. Only those faces
				      * are printed which are on the boundary
				      * and which have a boundary indicator
				      * not equal to zero, since the latter
				      * is the default for boundary faces.
				      *
				      * Since cells and faces are continuously
				      * numbered, the #starting_index# for
				      * the numbering of the faces is passed
				      * also.
				      *
				      * This function unfortunately can not
				      * be included in the regular #write_ucd#
				      * function, since it needs special
				      * treatment for the case #dim==1#, in
				      * which case the face iterators are
				      * #void*#'s and lack the member functions
				      * which are called. We would not actually
				      * call these functions, but the compiler
				      * would complain anyway when compiling
				      * the function for #dim==1#. Bad luck.
				      */
    void write_ucd_faces (ostream &out,
			  const unsigned int starting_index) const;
};




		


/*----------------------------   data_io.h     ---------------------------*/
/* end of #ifndef __data_io_H */
#endif
/*----------------------------   data_io.h     ---------------------------*/
