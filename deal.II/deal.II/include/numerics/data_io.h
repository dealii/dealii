/*----------------------------   io.h     ---------------------------*/
/*      <Id:>                 */
#ifndef __data_io_H
#define __data_io_H
/*----------------------------   io.h     ---------------------------*/

#include <base/exceptions.h>
#include <vector.h>

#ifdef DEBUG
#  define _G_NO_NRV     // don't use GNU's named return values in debug modes
#  include <String.h>   // because if we did we'd harvest tons of warnings.
#  undef _G_NO_NRV
#else
#  include <String.h>
#endif


template <int dim> class Triangulation;
template <int dim> class DoFHandler;

class dVector;

template <class Key, class T, class Compare> class map;
template <class T> struct less;


/**
  This class implements an input mechanism for grid data. It allows to
  read a grid structure into a triangulation object. Future versions
  will also allow to read data on this grid into vectors.

  At present, only UCD (unstructured cell data) is supported as input
  format for grid data. Any numerical data after the block of topological
  information is ignored.

  To read grid data, the triangulation to be fed with has to be empty.
  When giving a file which does not contain the assumed information or
  which does not keep to the right format, the state of the triangulation
  will be undefined afterwards. Upon input, only lines in one dimension
  and quads in two dimensions are accepted. All other cell types (e.g. lines
  or triangles in two dimensions, quads and hexes in 3d) are ignored. No
  warning is issued. The vertex and cell numbering in the UCD file, which
  need not be consecutively, is lost upon transfer to the triangulation
  object, since this one needs consecutively numbered elements.
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
    
  private:
				     /**
				      * Store address of the triangulation to
				      * be fed with the data read in.
				      */
    Triangulation<dim> *tria;
};







/**
  This class implements an output mechanism for grid and simulation data
  in several formats.
  At present it supports output in UCD (unstructured cell data) and
  partly in GNUPLOT format.

  It allows the user to attach a degree of freedom handler object
  (#DoFHandler#) which also gives access to the geometry data of the
  underlying triangulation and to add data vectors of which the values
  are to be written.


  {\bf Limitations}
  
  At present, no grouping of components to vectors is implemented, i.e.
  you can only write each component independent of the others. Also, it
  is not possible to output calculations which were performed on elements
  with more or less than one degree of freedom per vertex.

  
  {\bf UCD format}

  The UCD format is described in the AVS developer's guide. Due to
  limitations in the present format, only node based data can be output,
  so higher order elements are only written with their node values, no
  interior or line values are used. No use is made of the possibility
  to give cell and model data since these are not supported by all
  UCD aware programs.

  The ASCII UCD format is used. In future versions, a binary version may
  follow up.

  Note that to enumerate the vertices, not the vertex index is used but
  the index of the degree of freedom located on this vertex. This makes
  the mapping between the vertices and the entries in the data vectors
  much easier.


  {\bg GNUPLOT format}

  The gnuplot format is not able to handle data on unstructured grids, so
  the actual data is ignored in two and three space dimension. In one
  dimension, both mesh and data is written, in two and three dimensions
  only the grid is printed as a sequence of lines.

  For more than one dimension, the #DataOut<dim>::write_gnuplot()# somehow
  duplicates the functionality of the #Triangulation<dim>::print_gnuplot()#
  functions. These, however, offer more functionality.

  To view the results in two or three dimensions, use #set data style lines#
  within gnuplot and call #plot "filename"#. In one dimension call
  #plot "filename" using 1:x#. #x# denotes the number of the data set you
  want to see plus one. For example #using 1:4# would mean to plot the
  third data vector.
  */
template <int dim>  
class DataOut {
  public:
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
				      * (e.g. meter, kelvin, etc). By default,
				      * "<dimensionless>" is assumed for the
				      * units.
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
    void add_data_vector (dVector &data,
			  String  &name,
			  String  &units="<dimensionless>");

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
				      * Write data and grid in one dimension,
				      * only grid in two or three dimensions.
				      */
    void write_gnuplot (ostream &out) const;
    
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
    
  private:

				     /**
				      * Declare an entry in the list of
				      * data elements.
				      */
    struct DataEntry {
					 /**
					  * Pointer to the data vector.
					  */
	dVector *data;
					 /**
					  * Name of this component.
					  */
	String   name;
					 /**
					  * Physical unit name of this
					  * component.
					  */
	String   units;
    };

				     /**
				      * Pointer to the dof handler object.
				      */
    DoFHandler<dim>   *dofs;

				     /**
				      * List of data elements.
				      */
    vector<DataEntry>  data;
};





/*----------------------------   io.h     ---------------------------*/
/* end of #ifndef __data_io_H */
#endif
/*----------------------------   io.h     ---------------------------*/
