//----------------------------  data_out.h  ---------------------------
//    Version: $Name$
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  data_out.h  ---------------------------
#ifndef __deal2__data_out_h
#define __deal2__data_out_h


#include <base/data_out_base.h>
#include <base/smartpointer.h>
#include <base/multithread_info.h>


template <int dim> class DoFHandler;

/**
 * This is an abstract class which provides the functionality to generate
 * patches for output by base classes from data vectors on a grid. It alloes
 * to store a pointer to a @ref{DoFHandler} object and one or more pointers to
 * node and cell data denoting functions on the grid which shall later
 * be written in any of the implemented data formats.
 *
 *
 * @sect3{User visible interface}
 *
 * The user visible interface of this class consists of functions which allow
 * a user to make a @ref{DoFHandler} object known to this class and to add data
 * vectors which will later be written to a file in some format. Instead of
 * pondering about the different functions, an example is probably the best
 * way:
 * @begin{verbatim}
 *   ...
 *   ...   // compute solution, which is of type Vector<double>
 *   ...   // and contains nodal values
 *   ...
 *   ...   // compute error_estimator, which is of type Vector<double>
 *   ...   // and contains one value per cell
 *
 *   vector<string> solution_names;
 *   solution_names.push_back ("x-displacement");
 *   solution_names.push_back ("y-displacement");
 *
 *   DataOut<dim> data_out;
 *   data_out.attach_dof_handler (dof_handler);
 *   data_out.add_data_vector (solution, solution_names);
 *   data_out.add_data_vector (error_estimator, "estimated_error");
 *
 *   data_out.build_patches ();
 *
 *   ofstream output_file ("output");
 *   data_out.write_xxx (output_file);
 *
 *   data_out.clear();
 * @end{verbatim}
 *
 * @p{attach_dof_handler} tells this class that all future operations are to take
 * place with the @ref{DoFHandler} object and the triangulation it lives on. We then
 * add the solution vector and the error estimator; note that they have different
 * dimensions, because the solution is a nodal vector, here consisting of two
 * components ("x-displacement" and "y-displacement") while the error estimator
 * probably is a vector holding cell data. When attaching a data vector, you have
 * to give a name to each component of the vector, which is done through an object
 * of type @p{vector<string>} as second argument; if only one component is in the
 * vector, for example if we are adding cell data as in the second case, or if
 * the finite element used by the @ref{DoFHandler} has only one component, then you
 * can use the second @p{add_data_vector} function which takes a @p{string} instead
 * of the @p{vector<string>}.
 *
 * You should note that this class does not copy the vector given to it through
 * the @p{add_data_vector} functions, for memory consumption reasons. It only
 * stores a reference to it, so it is in your responsibility to make sure that
 * the data vectors exist long enough.
 *
 * After adding all data vectors, you need to call a function which generates
 * the patches for output from the stored data. This function is here called
 * @p{build_patches}, but the naming is up to the derived class that actually
 * implements this.
 *
 * Finally, you write the data in one format or other, indicated by @p{write_xxx},
 * to a file and may want to clear this object as soon as possible to reduce
 * memory requirements.
 *
 * Please note, that in the example above, an object of type @ref{DataOut} was
 * used, i.e. an object of a derived class. This is necessary since this
 * class does not provide means to actually generate the patches, only aids to
 * store and access data.
 *
 * Note that the base class of this class, @ref{DataOutInterface} offers several
 * functions to ease programming with run-time determinable output formats
 * (i.e. you need not use a fixed format by calling @p{write_xxx} in the above
 * example, but you can select it by a run-time parameter without having
 * to write the @p{if () ... else ...} clauses yourself), and also functions
 * and classes offering ways to control the appearance of the output by
 * setting flags for each output format.
 * 
 *
 * @sect3{Information for derived classes}
 *
 * What is actually missing this class is a way to produce the patches
 * for output itself, from the stored data and degree of freedom
 * information.  Since this task is often application dependent it is
 * left to derived classes. For example, in many applications, it
 * might be wanted to limit the depth of output to a certain number of
 * refinement levels and write data from finer cells only in a way
 * interpolated to coarser cells, to reduce the amount of
 * output. Also, it might be wanted to use different numbers of
 * subdivisions on different cells when forming a patch, for example
 * to accomplish for different polynomial degrees of the trial space
 * on different cells. Also, the output need not necessarily consist
 * of a patch for each cell, but might be made up of patches for
 * faces, of other things. Take a look at derived classes to what is
 * possible in this respect.
 *
 * For this reason, it is left to a derived class to provide a
 * function, named usually @p{build_patches} or the like, which fills
 * the @p{patches} array of this class.
 *
 * Regarding the templates of this class, it needs three values: first
 * the space dimension in which the triangulation and the DoF handler
 * operate, second the dimension of the objects which the patches
 * represent.  Although in most cases they are equal, there are also
 * classes for which this does not hold, for example if one outputs
 * the result of a computation exploiting rotational symmetry in the
 * original domain (in which the space dimension of the output would
 * be one higher than that of the DoF handler, see the
 * @ref{DataOut_Rotation} class), or one might conceive that one could
 * write a class that only outputs the solution on a cut through the
 * domain, in which case the space dimension of the output is less
 * than that of the DoF handler. The last template argument denotes
 * the dimension of the space into which the patches are embedded;
 * usually, this dimension is the same as the dimensio of the patches
 * themselves (which is also the default value of the template
 * parameter), but there might be cases where this is not so. For
 * example, in the @ref{DataOut_Faces} class, patches are generated
 * from faces of the triangulation. Thus, the dimension of the patch
 * is one less than the dimension of the embedding space, which is, in
 * this case, equal to the dimension of the triangulation and DoF
 * handler. However, for the cut through the domain mentioned above,
 * if the cut is a straight one, then the cut can be embedded into a
 * space of one dimension lower than the dimension of the
 * triangulation, so that the last template parameter has the same
 * value as the second one.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dof_handler_dim, int patch_dim, int patch_space_dim=patch_dim>
class DataOut_DoFData : public DataOutInterface<patch_dim,patch_space_dim>
{
  public:
				     /**
				      * Constructor
				      */
    DataOut_DoFData ();

				     /**
				      * Destructor.
				      */
    virtual ~DataOut_DoFData ();

    				     /**
				      * Designate a dof handler to be
				      * used to extract geometry data
				      * and the mapping between nodes
				      * and node values.
				      */
    void attach_dof_handler (const DoFHandler<dof_handler_dim> &);

				     /**
				      * Add a data vector together
				      * with its name and the physical
				      * unit (for example meter,
				      * kelvin, etc). By default,
				      * "<dimensionless>" is assumed
				      * for the units.
				      *
				      * A pointer to the vector is
				      * stored, so you have to make
				      * sure the vector exists at that
				      * address at least as long as
				      * you call the @p{write_*}
				      * functions.
				      *
				      * It is assumed that the vector
				      * has the same number of
				      * components as there are
				      * degrees of freedom in the dof
				      * handler, in which case it is
				      * assumed to be a vector storing
				      * nodal data; or the size may be
				      * the number of active cells on
				      * the present grid, in which
				      * case it is assumed to be a
				      * cell data vector.
				      *
				      * If it is a vector holding DoF
				      * data, the names given shall be
				      * one for each component, if the
				      * finite element in use is
				      * composed of several
				      * subelements.  If it is a
				      * finite element composed of
				      * only one subelement, then
				      * there is another function
				      * following which takes a single
				      * name instead of a vector of
				      * names.
				      *
				      * The names of a data vector
				      * shall only contain characters
				      * which are letters, underscore
				      * and a few other ones. Refer to
				      * the @p{ExcInvalidCharacter}
				      * exception declared in this
				      * class to see which characters
				      * are valid and which are not.
				      */
    void add_data_vector (const Vector<double>           &data,
			  const std::vector<std::string> &names);

				     /**
				      * This function is an
				      * abbreviation to the above one,
				      * intended for use with finite
				      * elements that are not composed
				      * of subelements. In this case,
				      * only one name per data vector
				      * needs to be given, which is
				      * what this function takes. It
				      * simply relays its arguments
				      * after a conversion of the
				      * @p{name} to a vector of
				      * strings, to the other
				      * @p{add_data_vector} function
				      * above.
				      *
				      * If @p{data} is a vector with
				      * multiple components this
				      * function will generate
				      * distinct names for all
				      * components by appending an
				      * underscore and the number of
				      * each component to @p{name}
				      */
    void add_data_vector (const Vector<double> &data,
			  const std::string    &name);

				     /**
				      * Release the pointers to the
				      * data vectors. This allows
				      * output of a new set of vectors
				      * without supplying the DoF
				      * handler again. Therefore, the
				      * @ref{DataOut} object can be used
				      * in an algebraic context. Note
				      * that besides the data vectors
				      * also the patches already
				      * computed are deleted.
				      */
    void clear_data_vectors ();

				     /**
				      * Release pointers to all input
				      * data elements, i.e. pointers
				      * to data vectors and to the DoF
				      * handler object. This function
				      * may be useful when you have
				      * called the @p{build_patches}
				      * function of derived class,
				      * since then the patches are
				      * built and the input data is no
				      * more needed, nor is there a
				      * need to reference it. You can
				      * then output the patches
				      * detached from the main thread
				      * and need not make sure anymore
				      * that the DoF handler object
				      * and vectors must not be
				      * deleted before the output
				      * thread is finished.
				      */
    void clear_input_data_references ();

				     /**
				      * Release the pointers to the
				      * data vectors and the DoF
				      * handler. You have to set all
				      * data entries again using the
				      * @p{add_data_vector}
				      * function. The pointer to the
				      * dof handler is cleared as
				      * well, along with all other
				      * data. In effect, this function
				      * resets everything to a virgin
				      * state.
				      */
    virtual void clear ();

    				     /**
				      * Determine an estimate for the
				      * memory consumption (in bytes)
				      * of this object.
				      */
    unsigned int memory_consumption () const;

				     /**
				      * Exception
				      */
    DeclException0 (ExcNoDoFHandlerSelected);
				     /**
				      * Exception
				      */
    DeclException3 (ExcInvalidVectorSize,
		    int, int, int,
		    << "The vector has size " << arg1
		    << " but the DoFHandler objects says there are " << arg2
		    << " degrees of freedom and there are " << arg3
		    << " active cells.");
				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidCharacter,
		    std::string,
		    << "Please use only the characters [a-zA-Z0-9_<>()] for" << std::endl
		    << "description strings since AVS will only accept these." << std::endl
		    << "The string you gave was <" << arg1 << ">.");
				     /**
				      * Exception
				      */
    DeclException0 (ExcOldDataStillPresent);
				     /**
				      * Exception
				      */
    DeclException2 (ExcInvalidNumberOfNames,
		    int, int,
		    << "You have to give one name per component in your "
		    << "data vector. The number you gave was " << arg1
		    << ", but the number of components is " << arg2);

  protected:
    				     /**
				      * Declare an entry in the list of
				      * data elements.
				      */
    struct DataEntry {
					 /**
					  * Constructor. If no arguments are
					  * given, an invalid object is
					  * constructed (we need a constructor
					  * with no explicit arguments for
					  * STL classes).
					  */
	DataEntry (const Vector<double>           *data = 0,
		   const std::vector<std::string> &names = std::vector<std::string>());

					 /**
					  * Determine an estimate for the
					  * memory consumption (in bytes)
					  * of this object.
					  */
	unsigned int memory_consumption () const;
	
					 /**
					  * Pointer to the data
					  * vector. Note that
					  * ownership of the vector
					  * pointed to remains with
					  * the caller of this class.
					  */
	const Vector<double> *data;
	
					 /**
					  * Names of the components of this
					  * data vector.
					  */
	std::vector<std::string> names;
	
					 /**
					  * Physical unit name of this
					  * component.
					  */
	std::string   units;
    };

				     /**
				      * Pointer to the dof handler object.
				      */
    SmartPointer<const DoFHandler<dof_handler_dim> > dofs;

				     /**
				      * List of data elements with vectors of
				      * values for each degree of freedom.
				      */
    typename std::vector<DataEntry>  dof_data;

				     /**
				      * List of data elements with vectors of
				      * values for each cell.
				      */
    typename std::vector<DataEntry>  cell_data;

				     /**
				      * This is a list of patches that is
				      * created each time @p{build_patches}
				      * is called. These patches are used
				      * in the output routines of the base
				      * classes.
				      */
    typename std::vector<DataOutBase::Patch<patch_dim,patch_space_dim> > patches;

				     /**
				      * Function by which the base
				      * class's functions get to know
				      * what patches they shall write
				      * to a file.
				      */
    virtual const typename std::vector<typename DataOutBase::Patch<patch_dim,patch_space_dim> > &
    get_patches () const;

				     /**
				      * Virtual function through
				      * which the names of data sets are
				      * obtained by the output functions
				      * of the base class.
				      */
    virtual std::vector<std::string> get_dataset_names () const;
};



/**
 * This class is an actual implementation of the functionality proposed by
 * the @ref{DataOut_DoFData} class. It offers a function @p{build_patches} that
 * generates the patches to be written in some graphics format from the data
 * stored in the base class. Most of the interface and an example of its
 * use is described in the documentation of the base class.
 *
 * The only thing this class offers is the function @p{build_patches} which
 * loops over all cells of the triangulation stored by the @p{attach_dof_handler}
 * function of the base class and convert the data on these to actual patches
 * which are the objects that are later output by the functions of the
 * base classes. You can give a parameter to the function which determines
 * how many subdivisions in each coordinate direction are to be performed,
 * i.e. of how many subcells each patch shall consist. Default is one, but
 * for quadratic elements you may want to choose two, for cubic elements three,
 * and so on.
 *
 * Note that after having called @p{build_patches} once, you can call one or
 * more of the @p{write_*} functions of the base classes. You can therefore
 * output the same data in more than one format without having to rebuild
 * the patches.
 *
 *
 * @sect3{User interface information}
 *
 * The base classes of this class, @ref{DataOutBase}, @ref{DataOutInterface} and
 * @ref{DataOut_DoFData} offer several interfaces of their own. Refer to the
 * @ref{DataOutBase} class's documentation for a discussion of the different
 * output formats presently supported, @ref{DataOutInterface} for ways of
 * selecting which format to use upon output at run-time and without
 * the need to adapt your program when new formats become available, as
 * well as for flags to determine aspects of output. The @ref{DataOut_DoFData}
 * class's documentation has an example of using nodal data to generate
 * output.
 *
 *
 * @sect3{Extensions}
 *
 * By default, this class produces patches for all active cells. Sometimes,
 * this is not what you want, maybe because they are simply too many (and too
 * small to be seen individually) or because you only want to see a certain
 * region of the domain, or for some other reason.
 *
 * For this, internally the @p{build_patches} function does not generate
 * the sequence of cells to be converted into patches itself, but relies
 * on the two functions @p{first_cell} and @p{next_cell}. By default, they
 * return the first active cell, and the next active cell, respectively.
 * Since they are @p{virtual} functions, you may overload them to select other
 * cells for output. If cells are not active, interpolated values are taken
 * for output instead of the exact values on active cells.
 *
 * The two functions are not constant, so you may store information within
 * your derived class about the last accessed cell. This is useful if the
 * information of the last cell which was accessed is not sufficient to
 * determine the next one.
 *
 * @author Wolfgang Bangerth, 1999
 */
template <int dim>
class DataOut : public DataOut_DoFData<dim,dim> 
{
  public:
    				     /**
				      * This is the central function
				      * of this class since it builds
				      * the list of patches to be
				      * written by the low-level
				      * functions of the base
				      * class. See the general
				      * documentation of this class
				      * for further information.
				      *
				      * The function supports
				      * multithreading, if deal.II is
				      * compiled in multithreading
				      * mode. The default number of
				      * threads to be used to build
				      * the patches is set to
				      * @p{multithread_info.n_default_threads}.
				      */
    virtual void build_patches (const unsigned int n_subdivisions = 1,
				const unsigned int n_threads      = multithread_info.n_default_threads);

				     /**
				      * Return the first cell which we
				      * want output for. The default
				      * implementation returns the
				      * first active cell, but you
				      * might want to return other
				      * cells in a derived class.
				      */
    virtual typename DoFHandler<dim>::cell_iterator
    first_cell ();
    
				     /**
				      * Return the next cell after
				      * @p{cell} which we want output
				      * for.  If there are no more
				      * cells, @p{dofs->end()} shall
				      * be returned.
				      *
				      * The default implementation
				      * returns the next active cell,
				      * but you might want to return
				      * other cells in a derived
				      * class. Note that the default
				      * implementation assumes that
				      * the given @p{cell} is active,
				      * which is guaranteed as long as
				      * @p{first_cell} is also used
				      * from the default
				      * implementation. Overloading
				      * only one of the two functions
				      * might not be a good idea.
				      */
    virtual typename DoFHandler<dim>::cell_iterator
    next_cell (const typename DoFHandler<dim>::cell_iterator &cell);

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumberOfSubdivisions,
		    int,
		    << "The number of subdivisions per patch, " << arg1
		    << ", is not valid.");
    
  private:
				     /**
				      * All data needed in one thread
				      * is gathered in the struct
				      * Data.
				      * The data is handled globally
				      * to avoid allocation of memory
				      * in the threads.
				      */
    struct Data 
    {
	unsigned int n_threads;
	unsigned int this_thread;
	unsigned int n_components;
	unsigned int n_datasets;
	unsigned int n_subdivisions;
	std::vector<double>          patch_values;
	std::vector<Vector<double> > patch_values_system;
	Data ()
	  {}
    };
				     /**
				      * Builds every @p{n_threads}'s
				      * patch. This function may be
				      * called in parallel.
				      * If multithreading is not
				      * used, the function is called
				      * once and generates all patches.
				      */
    void build_some_patches (Data data);
};


#endif
