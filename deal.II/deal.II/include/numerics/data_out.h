//---------------------------------------------------------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//---------------------------------------------------------------------------
#ifndef __deal2__data_out_h
#define __deal2__data_out_h



#include <base/config.h>
#include <base/smartpointer.h>
#include <base/data_out_base.h>
#include <dofs/dof_handler.h>
#include <fe/mapping.h>
#include <hp/q_collection.h>
#include <hp/fe_collection.h>
#include <hp/mapping_collection.h>
#include <hp/fe_values.h>
#include <numerics/data_postprocessor.h>
#include <numerics/data_component_interpretation.h>

#include <base/std_cxx1x/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

template <int, int> class FEValuesBase;

namespace internal
{
  namespace DataOut
  {
				     /**
				      * A data structure that holds
				      * all data needed in one thread
				      * when building patches in
				      * parallel.  These data
				      * structures are created
				      * globally rather than on each
				      * cell to avoid allocation of
				      * memory in the threads. This is
				      * a base class for the
				      * AdditionalData kind of data
				      * structure discussed in the
				      * documentation of the
				      * WorkStream class.
				      *
				      * The
				      * <code>cell_to_patch_index_map</code>
				      * is an array that stores for index
				      * <tt>[i][j]</tt> the number of the
				      * patch that associated with the cell
				      * with index @p j on level @p i. This
				      * information is set up prior to
				      * generation of the patches, and is
				      * needed to generate neighborship
				      * information.
				      *
				      * This structure is used by
				      * several of the DataOut*
				      * classes, which derived their
				      * own ParallelData classes from
				      * it for additional fields.
				      */
    template <int dim, int spacedim>
    struct ParallelDataBase
    {
	template <class FE>
	ParallelDataBase (const unsigned int n_components,
			  const unsigned int n_datasets,
			  const unsigned int n_subdivisions,
			  const unsigned int n_q_points,
			  const std::vector<unsigned int> &n_postprocessor_outputs,
			  const FE &finite_elements);
			  
	const unsigned int n_components;
	const unsigned int n_datasets;
	const unsigned int n_subdivisions;

	std::vector<double>                                patch_values;
	std::vector<dealii::Vector<double> >               patch_values_system;
	std::vector<Tensor<1,spacedim> >                   patch_gradients;
	std::vector<std::vector<Tensor<1,spacedim> > >     patch_gradients_system;
	std::vector<Tensor<2,spacedim> >                   patch_hessians;
	std::vector<std::vector<Tensor<2,spacedim> > >     patch_hessians_system;
	std::vector<std::vector<dealii::Vector<double> > > postprocessed_values;

	const dealii::hp::FECollection<dim,spacedim>      fe_collection;	
    };


				     /**
				      * A derived class for use in the
				      * DataOut class. This is
				      * a class for the
				      * AdditionalData kind of data
				      * structure discussed in the
				      * documentation of the
				      * WorkStream context.
				      */
    template <int dim, int spacedim>
    struct ParallelData : public ParallelDataBase<dim,spacedim>
    {
	template <class FE>
	ParallelData (const Quadrature<dim> &quadrature,
		      const unsigned int n_components,
		      const unsigned int n_datasets,
		      const unsigned int n_subdivisions,
		      const std::vector<unsigned int> &n_postprocessor_outputs,
		      const Mapping<dim,spacedim> &mapping,
		      const std::vector<std::vector<unsigned int> > &cell_to_patch_index_map,
		      const FE &finite_elements,
		      const UpdateFlags update_flags);

	const dealii::hp::QCollection<dim> q_collection;
	const dealii::hp::MappingCollection<dim,spacedim> mapping_collection;
	dealii::hp::FEValues<dim,spacedim> x_fe_values;

	const std::vector<std::vector<unsigned int> > *cell_to_patch_index_map;
    };
  }
}


//TODO: Most of the documentation of DataOut_DoFData applies to DataOut.

/**
 * This is an abstract class which provides the functionality to
 * generate patches for output by base classes from data vectors on a
 * grid. It allows to store a pointer to a DoFHandler object and one
 * or more pointers to node and cell data denoting functions on the
 * grid which shall later be written in any of the implemented data
 * formats.
 *
 *
 * <h3>User visible interface</h3>
 *
 * The user visible interface of this class consists of functions which allow
 * a user to make a DoFHandler object known to this class and to add data
 * vectors which will later be written to a file in some format. Instead of
 * pondering about the different functions, an example is probably the best
 * way:
 * @code
 *   ...
 *   ...   // compute solution, which contains nodal values
 *   ...
 *   ...   // compute error_estimator, which contains one value per cell
 *
 *   std::vector<std::string> solution_names;
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
 * @endcode
 *
 * attach_dof_handler() tells this class that all future operations
 * are to take place with the DoFHandler object and the triangulation
 * it lives on. We then add the solution vector and the error
 * estimator; note that they have different dimensions, because the
 * solution is a nodal vector, here consisting of two components
 * ("x-displacement" and "y-displacement") while the error estimator
 * probably is a vector holding cell data. When attaching a data
 * vector, you have to give a name to each component of the vector,
 * which is done through an object of type <tt>vector<string></tt> as
 * second argument; if only one component is in the vector, for
 * example if we are adding cell data as in the second case, or if the
 * finite element used by the DoFHandler has only one component, then
 * you can use the second add_data_vector() function which takes a @p
 * string instead of the <tt>vector<string></tt>.
 *
 * The add_data_vector() functions have additional arguments (with default
 * values) that can be used to specify certain transformations. In particular,
 * it allows to attach DataPostprocessor arguments to compute derived
 * information from a data vector at each quadrature point (for example, the
 * Mach number in hypersonic flow can be computed from density and velocities;
 * @ref step_29 "step-29" also shows an example); another piece of information
 * specified through arguments with default values is how certain output
 * components should be interpreted, i.e. whether each component of the data
 * is logically an independent scalar field, or whether some of them together
 * form logically a vector-field (see the
 * DataComponentInterpretation::DataComponentInterpretation enum, and the @ref
 * step_22 "step-22" tutorial program).
 *
 * It should be noted that this class does not copy the vector given to it through
 * the add_data_vector() functions, for memory consumption reasons. It only
 * stores a reference to it, so it is in your responsibility to make sure that
 * the data vectors exist long enough.
 *
 * After adding all data vectors, you need to call a function which generates
 * the patches for output from the stored data. Derived classes name this
 * function build_patches(). Finally, you write() the data in one format or other,
 * to a file.
 *
 * Please note, that in the example above, an object of type DataOut was
 * used, i.e. an object of a derived class. This is necessary since this
 * class does not provide means to actually generate the patches, only aids to
 * store and access data.
 *
 * Note that the base class of this class, DataOutInterface offers
 * several functions to ease programming with run-time determinable
 * output formats (i.e. you need not use a fixed format by calling
 * DataOutInterface::write_xxx in the above example, but you can
 * select it by a run-time parameter without having to write the
 * <tt>if () ... else ...</tt> clauses yourself), and also functions
 * and classes offering ways to control the appearance of the output
 * by setting flags for each output format.
 * 
 *
 * <h3>Information for derived classes</h3>
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
 * function, named usually build_patches() or the like, which fills
 * the #patches array of this class.
 *
 * Regarding the templates of this class, it needs three values: first
 * the space dimension in which the triangulation and the DoF handler
 * operate, second the dimension of the objects which the patches
 * represent.  Although in most cases they are equal, there are also
 * classes for which this does not hold, for example if one outputs
 * the result of a computation exploiting rotational symmetry in the
 * original domain (in which the space dimension of the output would
 * be one higher than that of the DoF handler, see the
 * DataOut_Rotation() class), or one might conceive that one could
 * write a class that only outputs the solution on a cut through the
 * domain, in which case the space dimension of the output is less
 * than that of the DoF handler. The last template argument denotes
 * the dimension of the space into which the patches are embedded;
 * usually, this dimension is the same as the dimensio of the patches
 * themselves (which is also the default value of the template
 * parameter), but there might be cases where this is not so. For
 * example, in the DataOut_Faces() class, patches are generated
 * from faces of the triangulation. Thus, the dimension of the patch
 * is one less than the dimension of the embedding space, which is, in
 * this case, equal to the dimension of the triangulation and DoF
 * handler. However, for the cut through the domain mentioned above,
 * if the cut is a straight one, then the cut can be embedded into a
 * space of one dimension lower than the dimension of the
 * triangulation, so that the last template parameter has the same
 * value as the second one.
 *
 * @ingroup output
 * @author Wolfgang Bangerth, 1999
 */
template <class DH, int patch_dim, int patch_space_dim=patch_dim>
class DataOut_DoFData : public DataOutInterface<patch_dim,patch_space_dim>
{
  public:

				     /**
				      * Typedef to the iterator type
				      * of the dof handler class under
				      * consideration.
				      */
    typedef typename DH::cell_iterator cell_iterator;
    typedef typename DH::active_cell_iterator active_cell_iterator;

  public:

				     /**
				      * Type describing what the
				      * vector given to
				      * add_data_vector() is: a
				      * vector that has one entry per
				      * degree of freedom in a
				      * DoFHandler object (such
				      * as solution vectors), or one
				      * entry per cell in the
				      * triangulation underlying the
				      * DoFHandler object (such
				      * as error per cell data). The
				      * value #type_automatic tells
				      * add_data_vector() to find
				      * out itself (see the
				      * documentation of
				      * add_data_vector() for the
				      * method used).
				      */
    enum DataVectorType
    {
					   /**
					    * Data vector entries are
					    * associated to degrees of freedom
					    */
	  type_dof_data,
	  
					   /**
					    * Data vector entries are one per
					    * grid cell
					    */
	  type_cell_data,
	  
					   /**
					    * Find out automatically
					    */
	  type_automatic
    };
    
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
    void attach_dof_handler (const DH &);

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
				      * you call the <tt>write_*</tt>
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
				      * cell data vector. As the
				      * number of degrees of freedom
				      * and of cells is usually not
				      * equal, the function can
				      * determine itself which type of
				      * vector it is given; however,
				      * there is one corner case,
				      * namely if you compute with
				      * piecewise constant elements
				      * and have one scalar quantity,
				      * then there are as many cells
				      * as there are degrees of
				      * freedom, but they may be
				      * ordered differently. In that
				      * case, you can change the last
				      * argument of the function from
				      * its default value
				      * #type_automatic to either
				      * #type_dof_data or #type_cell_data,
				      * depending on what the vector
				      * represents. Apart from this
				      * corner case, you can leave the
				      * argument at its default value
				      * and let the function determine
				      * the type of the vector itself.
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
				      * The data_component_interpretation
				      * argument contains information about
				      * how the individual components of
				      * output files that consist of more than
				      * one data set are to be interpreted.
				      *
				      * For example, if one has a finite
				      * element for the Stokes equations in
				      * 2d, representing components (u,v,p),
				      * one would like to indicate that the
				      * first two, u and v, represent a
				      * logical vector so that later on when
				      * we generate graphical output we can
				      * hand them off to a visualization
				      * program that will automatically know
				      * to render them as a vector field,
				      * rather than as two separate and
				      * independent scalar fields.
				      *
				      * The default value of this argument
				      * (i.e. an empty vector) corresponds is
				      * equivalent to a vector of values
				      * DataComponentInterpretation::component_is_scalar,
				      * indicating that all output components
				      * are independent scalar
				      * fields. However, if the given data
				      * vector represents logical vectors, you
				      * may pass a vector that contains values
				      * DataComponentInterpretation::component_is_part_of_vector. In
				      * the example above, one would pass in a
				      * vector with components
				      * (DataComponentInterpretation::component_is_part_of_vector,
				      * DataComponentInterpretation::component_is_part_of_vector,
				      * DataComponentInterpretation::component_is_scalar)
				      * for (u,v,p).
				      *
				      * The names of a data vector
				      * shall only contain characters
				      * which are letters, underscore
				      * and a few other ones. Refer to
				      * the ExcInvalidCharacter
				      * exception declared in this
				      * class to see which characters
				      * are valid and which are not.
				      *
				      * The actual type for the template
				      * argument may be any vector type from
				      * which FEValues can extract values
				      * on a cell using the
				      * FEValuesBase::get_function_values() function.
				      */
    template <class VECTOR>
    void add_data_vector (const VECTOR                   &data,
			  const std::vector<std::string> &names,
			  const DataVectorType            type = type_automatic,
			  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
			  = std::vector<DataComponentInterpretation::DataComponentInterpretation>());

				     /**
				      * This function is an abbreviation to
				      * the above one (see there for a
				      * discussion of the various arguments),
				      * intended for use with finite elements
				      * that are not composed of
				      * subelements. In this case, only one
				      * name per data vector needs to be
				      * given, which is what this function
				      * takes. It simply relays its arguments
				      * after a conversion of the @p name to a
				      * vector of strings, to the other
				      * add_data_vector() function above.
				      *
				      * If @p data is a vector with
				      * multiple components this
				      * function will generate
				      * distinct names for all
				      * components by appending an
				      * underscore and the number of
				      * each component to @p name
				      *
				      * The actual type for the template
				      * argument may be any vector type from
				      * which FEValues can extract values
				      * on a cell using the
				      * FEValuesBase::get_function_values() function.
				      */
    template <class VECTOR>
    void add_data_vector (const VECTOR         &data,
			  const std::string    &name,
			  const DataVectorType  type = type_automatic,
			  const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
			  = std::vector<DataComponentInterpretation::DataComponentInterpretation>());

				     /**
				      * This function is an alternative to the
				      * above ones, allowing the output of
				      * derived quantities instead of the given
				      * data. This converison has to be done in
				      * a class derived from DataPostprocessor.
				      *
				      * The names for these derived quantities
				      * are provided by the @p
				      * data_postprocessor argument. Likewise,
				      * the data_component_interpretation
				      * argument of the other
				      * add_data_vector() functions is
				      * provided by the data_postprocessor
				      * argument. As only data of type @p
				      * type_dof_data can be transformed, this
				      * type is also known implicitly and does
				      * not have to be given.
				      *
				      * The actual type for the template
				      * argument may be any vector type from
				      * which FEValues can extract values
				      * on a cell using the
				      * FEValuesBase::get_function_values() function.
				      */
    template <class VECTOR>
    void add_data_vector (const VECTOR                           &data,
			  const DataPostprocessor<DH::space_dimension> &data_postprocessor);

				     /**
				      * Release the pointers to the
				      * data vectors. This allows
				      * output of a new set of vectors
				      * without supplying the DoF
				      * handler again. Therefore, the
				      * DataOut object can be used
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
				      * called the @p build_patches
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
                                      * This function can be used to
                                      * merge the patches that were
                                      * created using the
                                      * @p build_patches function of
                                      * the object given as argument
                                      * into the list of patches
                                      * created by this object. This
                                      * is sometimes handy if one has,
                                      * for example, a domain
                                      * decomposition algorithm where
                                      * each block is represented by a
                                      * DoFHandler of its own,
                                      * but one wants to output the
                                      * solution on all the blocks at
                                      * the same time.
                                      *
                                      * For this to work, the given
                                      * argument and this object need
                                      * to have the same number of
                                      * output vectors, and they need
                                      * to use the same number of
                                      * subdivisions per patch. The
                                      * output will probably look
                                      * rather funny if patches in
                                      * both objects overlap in space.
                                      *
                                      * If you call
                                      * build_patches() for this
                                      * object after merging in
                                      * patches, the previous state is
                                      * overwritten, and the merged-in
                                      * patches are lost.
                                      *
				      * The second parameter allows to
				      * shift each node of the patches
				      * in the object passed in in the
				      * first parameter by a certain
				      * amount. This is sometimes
				      * useful to generate "exploded"
				      * views of a collection of
				      * blocks.
				      *
                                      * This function will fail if
                                      * either this or the other
                                      * object did not yet set up any
                                      * patches.
                                      */
    template <class DH2>
    void merge_patches (const DataOut_DoFData<DH2,patch_dim,patch_space_dim> &source,
			const Point<patch_space_dim> &shift = Point<patch_space_dim>());
    
				     /**
				      * Release the pointers to the
				      * data vectors and the DoF
				      * handler. You have to set all
				      * data entries again using the
				      * add_data_vector()
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
    DeclException0 (ExcDataPostprocessingIsNotPossibleForCellData);
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
    DeclException2 (ExcInvalidCharacter,
		    std::string, size_t,
		    << "Please use only the characters [a-zA-Z0-9_<>()] for" << std::endl
		    << "description strings since some graphics formats will only accept these."
		    << std::endl
		    << "The string you gave was <" << arg1
		    << ">, the invalid character is <" << arg1[arg2]
		    << ">." << std::endl);
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

    DeclException2 (ExcInvalidVectorDeclaration,
		    int, std::string,
		    << "When declaring that a number of components in a data\n"
		    << "set to be output logically form a vector instead of\n"
		    << "simply a set of scalar fields, you need to specify\n"
		    << "this for all relevant components. Furthermore,\n"
		    << "vectors must always consist of exactly <dim>\n"
		    << "components. However, the vector component at\n"
		    << "position " << arg1 << " with name <" << arg2
		    << "> does not satisfy these conditions.");
    
  protected:
    				     /**
                                      * For each vector that has been added
                                      * through the add_data_vector()
                                      * functions, we need to keep track of a
                                      * pointer to it, and allow data
                                      * extraction from it when we generate
                                      * patches. Unfortunately, we need to do
                                      * this for a number of different vector
                                      * types. Fortunately, they all have the
                                      * same interface. So the way we go is to
                                      * have a base class that provides the
                                      * functions to access the vector's
                                      * information, and to have a derived
                                      * template class that can be
                                      * instantiated for each vector
                                      * type. Since the vectors all have the
                                      * same interface, this is no big
                                      * problem, as they can all use the same
                                      * general templatized code.
                                      *
                                      * @author Wolfgang Bangerth, 2004
				      */
    class DataEntryBase
    {
      public:
					 /**
					  * Constructor. Give a list of names
					  * for the individual components of
					  * the vector and their
					  * interpretation as scalar or vector
					  * data. This constructor assumes
					  * that no postprocessor is going to
					  * be used.
					  */
	DataEntryBase (const std::vector<std::string> &names,
		       const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation);

					 /**
					  * Constructor when a data
					  * postprocessor is going to be
					  * used. In that case, the names and
					  * vector declarations are going to
					  * be aquired from the postprocessor.
					  */
	DataEntryBase (const DataPostprocessor<DH::space_dimension> *data_postprocessor);
	
                                         /**
                                          * Destructor made virtual.
                                          */
        virtual ~DataEntryBase ();
        
                                         /**
                                          * Assuming that the stored vector is
                                          * a cell vector, extract the given
                                          * element from it.
                                          */
        virtual
        double
        get_cell_data_value (const unsigned int cell_number) const = 0;

                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the values on the present
                                          * cell from the vector we actually
                                          * store.
                                          */
        virtual
        void
        get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<double>             &patch_values) const = 0;
        
                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the values on the present
                                          * cell from the vector we actually
                                          * store. This function does the same
                                          * as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<Vector<double> >    &patch_values_system) const = 0;

                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the gradients on the present
                                          * cell from the vector we actually
                                          * store.
                                          */
        virtual
        void
        get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
				std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const = 0;
        
                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the gradients on the present
                                          * cell from the vector we actually
                                          * store. This function does the same
                                          * as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
				std::vector<std::vector<Tensor<1,DH::space_dimension> > > &patch_gradients_system) const = 0;

                                         /**
                                          * Given a FEValuesBase object, extract
                                          * the second derivatives on the
                                          * present cell from the vector we
                                          * actually store.
                                          */
        virtual
        void
        get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			       std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const = 0;
        
                                         /**
                                          * Given a FEValuesBase object, extract
                                          * the second derivatives on the
                                          * present cell from the vector we
                                          * actually store. This function does
                                          * the same as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			       std::vector<std::vector< Tensor<2,DH::space_dimension> > > &patch_hessians_system) const = 0;

                                         /**
                                          * Clear all references to the
                                          * vectors.
                                          */
        virtual void clear () = 0;
        
					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this object.
					  */
	virtual unsigned int memory_consumption () const = 0;

					 /**
					  * Names of the components of this
					  * data vector.
					  */
	const std::vector<std::string> names;

					 /**
					  * A vector that for each of the
					  * n_output_variables variables of
					  * the current data set indicates
					  * whether they are scalar fields,
					  * parts of a vector-field, or any of
					  * the other supported kinds of data.
					  */
	const std::vector<DataComponentInterpretation::DataComponentInterpretation>
	data_component_interpretation;
	
					 /**
					  * Pointer to a DataPostprocessing
					  * object which shall be applied to
					  * this data vector.
					  */
	SmartPointer<const DataPostprocessor<DH::space_dimension> > postprocessor;

					 /**
					  * Number of output variables this
					  * dataset provides (either number of
					  * components in vector valued function
					  * / data vector or number of computed
					  * quantities, if DataPostprocessor is
					  * applied). This variable is
					  * determined via and thus equivalent
					  * to <tt>names.size()</tt>.
					  */
	unsigned int n_output_variables;
    };


                                     /**
                                      * Class that stores a pointer to a
                                      * vector of type equal to the template
                                      * argument, and provides the functions
                                      * to extract data from it.
                                      *
                                      * @author Wolfgang Bangerth, 2004
                                      */
    template <typename VectorType>
    class DataEntry : public DataEntryBase
    {
      public:
					 /**
					  * Constructor. Give a list of names
					  * for the individual components of
					  * the vector and their
					  * interpretation as scalar or vector
					  * data. This constructor assumes
					  * that no postprocessor is going to
					  * be used.
					  */
	DataEntry (const VectorType               *data,
		   const std::vector<std::string> &names,
		   const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation);

					 /**
					  * Constructor when a data
					  * postprocessor is going to be
					  * used. In that case, the names and
					  * vector declarations are going to
					  * be aquired from the postprocessor.
					  */
	DataEntry (const VectorType                       *data,
		   const DataPostprocessor<DH::space_dimension> *data_postprocessor);

                                         /**
                                          * Assuming that the stored vector is
                                          * a cell vector, extract the given
                                          * element from it.
                                          */
        virtual
        double
        get_cell_data_value (const unsigned int cell_number) const;

                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the values on the present
                                          * cell from the vector we actually
                                          * store.
                                          */
        virtual
        void
        get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<double>             &patch_values) const;
        
                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the values on the present
                                          * cell from the vector we actually
                                          * store. This function does the same
                                          * as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<Vector<double> >    &patch_values_system) const;

                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the gradients on the present
                                          * cell from the vector we actually
                                          * store.
                                          */
        virtual
        void
        get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
				std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const;
        
                                         /**
                                          * Given a FEValuesBase object,
                                          * extract the gradients on the present
                                          * cell from the vector we actually
                                          * store. This function does the same
                                          * as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
				std::vector<std::vector<Tensor<1,DH::space_dimension> > > &patch_gradients_system) const;

                                         /**
                                          * Given a FEValuesBase object, extract
                                          * the second derivatives on the
                                          * present cell from the vector we
                                          * actually store.
                                          */
        virtual
        void
        get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			       std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const;
        
                                         /**
                                          * Given a FEValuesBase object, extract
                                          * the second derivatives on the
                                          * present cell from the vector we
                                          * actually store. This function does
                                          * the same as the one above but for
                                          * vector-valued finite elements.
                                          */
        virtual
        void
        get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
			       std::vector<std::vector< Tensor<2,DH::space_dimension> > > &patch_hessians_system) const;

                                         /**
                                          * Clear all references to the
                                          * vectors.
                                          */
        virtual void clear ();
        
					 /**
					  * Determine an estimate for
					  * the memory consumption (in
					  * bytes) of this object.
					  */
	virtual unsigned int memory_consumption () const;

      private:
                                         /**
                                          * Pointer to the data
                                          * vector. Note that
                                          * ownership of the vector
                                          * pointed to remains with
                                          * the caller of this class.
                                          */
        const VectorType *vector;
    };


				     /**
				      * Abbreviate the somewhat lengthy
				      * name for the Patch class.
				      */
    typedef dealii::DataOutBase::Patch<patch_dim,patch_space_dim> Patch;

				     /**
				      * Pointer to the dof handler object.
				      */
    SmartPointer<const DH> dofs;

				     /**
				      * List of data elements with vectors of
				      * values for each degree of freedom.
				      */
    std::vector<std_cxx1x::shared_ptr<DataEntryBase> >  dof_data;

				     /**
				      * List of data elements with vectors of
				      * values for each cell.
				      */
    std::vector<std_cxx1x::shared_ptr<DataEntryBase> >  cell_data;

				     /**
				      * This is a list of patches that is
				      * created each time build_patches()
				      * is called. These patches are used
				      * in the output routines of the base
				      * classes.
				      */
    std::vector<Patch> patches;
    
				     /**
				      * Function by which the base
				      * class's functions get to know
				      * what patches they shall write
				      * to a file.
				      */
    virtual
    const std::vector<Patch> & get_patches () const;

				     /**
				      * Virtual function through
				      * which the names of data sets are
				      * obtained by the output functions
				      * of the base class.
				      */
    virtual
    std::vector<std::string> get_dataset_names () const;

				     /**
				      * Overload of the respective
				      * DataOutInterface::get_vector_data_ranges()
				      * function. See there for a more
				      * extensive documentation.
				      */
    virtual
    std::vector<std_cxx1x::tuple<unsigned int, unsigned int, std::string> >
    get_vector_data_ranges () const;
    
				     /**
				      * Make all template siblings
				      * friends. Needed for the
				      * merge_patches() function.
				      */
    template <class, int, int>
    friend class DataOut_DoFData;

#ifdef DEAL_II_NESTED_CLASS_FRIEND_BUG
                                     /**
                                      * Make DataEntry a friend. This should
                                      * not be strictly necessary, since
                                      * members are implicitly friends, but in
                                      * this case it seems as if icc needs
                                      * this. Otherwise, it complains that
                                      * DataEntry can't derive from
                                      * DataEntryBase since the latter is a
                                      * private member of DataOut_DoFData.
				      *
				      * For whatever weird reason, it is also
				      * not enough to make just DataEntry a
				      * friend, but we have to fully qualify
				      * it for icc, while gcc 2.95 insists on
				      * the non-qualified version...
                                      */
#  ifdef DEAL_II_NESTED_CLASS_TEMPL_FRIEND_BUG
    template <typename>
    friend class DataEntry;
#  else
    template <int N1, template <int> class DH1, int N2, int N3>
    template <typename>
    friend class DataOut_DoFData<DH1,N2,N3>::DataEntry;
#  endif
#endif
};



/**
 * This class is an actual implementation of the functionality proposed by
 * the DataOut_DoFData class. It offers a function build_patches() that
 * generates the patches to be written in some graphics format from the data
 * stored in the base class. Most of the interface and an example of its
 * use is described in the documentation of the base class.
 *
 * The only thing this class offers is the function build_patches() which
 * loops over all cells of the triangulation stored by the
 * attach_dof_handler() function of the base class and convert the data on
 * these to actual patches which are the objects that are later output by the
 * functions of the base classes. You can give a parameter to the function
 * which determines how many subdivisions in each coordinate direction are to
 * be performed, i.e. of how many subcells each patch shall consist. Default
 * is one, but you may want to choose a higher number for higher order
 * elements, so as two for quadratic elements, three for cubic elements three,
 * and so on. The purpose of this parameter is because most graphics programs
 * do not allow to specify higher order shape functions in the file formats:
 * only data at vertices can be plotted and is then shown as a bilinear
 * interpolation within the interior of cells. This may be insufficient if you
 * have higher order finite elements, and the only way to achieve better
 * output is to subdivide each cell of the mesh into several cells for
 * graphical output.
 *
 * Note that after having called build_patches() once, you can call one or
 * more of the write() functions of DataOutInterface. You can therefore
 * output the same data in more than one format without having to rebuild
 * the patches.
 *
 *
 * <h3>User interface information</h3>
 *
 * The base classes of this class, DataOutBase, DataOutInterface and
 * DataOut_DoFData() offer several interfaces of their own. Refer to the
 * DataOutBase class's documentation for a discussion of the different
 * output formats presently supported, DataOutInterface for ways of
 * selecting which format to use upon output at run-time and without
 * the need to adapt your program when new formats become available, as
 * well as for flags to determine aspects of output. The DataOut_DoFData()
 * class's documentation has an example of using nodal data to generate
 * output.
 *
 *
 * <h3>Extensions</h3>
 *
 * By default, this class produces patches for all active cells. Sometimes,
 * this is not what you want, maybe because they are simply too many (and too
 * small to be seen individually) or because you only want to see a certain
 * region of the domain (for example in parallel programs such as the @ref step_18 "step-18"
 * example program), or for some other reason.
 *
 * For this, internally build_patches() does not generate
 * the sequence of cells to be converted into patches itself, but relies
 * on the two functions first_cell() and next_cell(). By default, they
 * return the first active cell, and the next active cell, respectively.
 * Since they are @p virtual functions, you may overload them to select other
 * cells for output. If cells are not active, interpolated values are taken
 * for output instead of the exact values on active cells.
 *
 * The two functions are not constant, so you may store information within
 * your derived class about the last accessed cell. This is useful if the
 * information of the last cell which was accessed is not sufficient to
 * determine the next one.
 *
 * There is one caveat, however: if you have cell data (in contrast to nodal,
 * or dof, data) such as error indicators, then you must make sure that
 * first_cell() and next_cell() only walk over active cells, since cell data
 * cannot be interpolated to a coarser cell. If you do have cell data and use
 * this pair of functions and they return a non-active cell, then an exception
 * will be thrown.
 * 
 * @ingroup output
 * @author Wolfgang Bangerth, 1999
 */
template <int dim, class DH=DoFHandler<dim> >
class DataOut : public DataOut_DoFData<DH, DH::dimension, DH::space_dimension> 
{
  public:
				     /**
				      * Typedef to the iterator type
				      * of the dof handler class under
				      * consideration.
				      */
    typedef typename DataOut_DoFData<DH,DH::dimension,DH::space_dimension>::cell_iterator cell_iterator;
    typedef typename DataOut_DoFData<DH,DH::dimension,DH::space_dimension>::active_cell_iterator active_cell_iterator;

				     /**
				      * Enumeration describing the region of the
				      * domain in which curved cells shall be
				      * created.
				      */
    enum CurvedCellRegion
    {
	  no_curved_cells,
	  curved_boundary,
	  curved_inner_cells
    };
    
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
				      * The default value <tt>0</tt>
				      * of <tt>n_subdivisions</tt>
				      * indicates that the value
				      * stored in
				      * DataOutInterface::default_subdivisions
				      * is to be used.
				      */
    virtual void build_patches (const unsigned int n_subdivisions = 0);

				     /**
				      * Same as above, except that the
				      * additional first parameter
				      * defines a mapping that is to
				      * be used in the generation of
				      * output. If
				      * <tt>n_subdivisions>1</tt>, the
				      * points interior of subdivided
				      * patches which originate from
				      * cells at the boundary of the
				      * domain can be computed using the
				      * mapping, i.e. a higher order
				      * mapping leads to a
				      * representation of a curved
				      * boundary by using more
				      * subdivisions. Some mappings like
				      * MappingQEulerian result in curved cells
				      * in the interior of the domain. However,
				      * there is nor easy way to get this
				      * information from the Mapping. Thus the
				      * last argument @p curved_region take one
				      * of three values resulting in no curved
				      * cells at all, curved cells at the
				      * boundary (default) or curved cells in
				      * the whole domain.
				      *
				      * Even for non-curved cells the
				      * mapping argument can be used
				      * for the Eulerian mappings (see
				      * class MappingQ1Eulerian) where
				      * a mapping is used not only to
				      * determine the position of
				      * points interior to a cell, but
				      * also of the vertices.  It
				      * offers an opportunity to watch
				      * the solution on a deformed
				      * triangulation on which the
				      * computation was actually
				      * carried out, even if the mesh
				      * is internally stored in its
				      * undeformed configuration and
				      * the deformation is only
				      * tracked by an additional
				      * vector that holds the
				      * deformation of each vertex.
				      *
				      * @todo The @p mapping argument should be
				      * replaced by a hp::MappingCollection in
				      * case of a hp::DoFHandler.
				      */
    virtual void build_patches (const Mapping<DH::dimension,DH::space_dimension> &mapping,
				const unsigned int n_subdivisions = 0,
				const CurvedCellRegion curved_region = curved_boundary);
    
				     /**
				      * Return the first cell which we
				      * want output for. The default
				      * implementation returns the
				      * first active cell, but you
				      * might want to return other
				      * cells in a derived class.
				      */
    virtual cell_iterator first_cell ();
    
				     /**
				      * Return the next cell after
				      * @p cell which we want output
				      * for.  If there are no more
				      * cells, <tt>#dofs->end()</tt> shall
				      * be returned.
				      *
				      * The default implementation
				      * returns the next active cell,
				      * but you might want to return
				      * other cells in a derived
				      * class. Note that the default
				      * implementation assumes that
				      * the given @p cell is active,
				      * which is guaranteed as long as
				      * first_cell() is also used
				      * from the default
				      * implementation. Overloading
				      * only one of the two functions
				      * might not be a good idea.
				      */
    virtual cell_iterator next_cell (const cell_iterator &cell);

				     /**
				      * Exception
				      */
    DeclException1 (ExcInvalidNumberOfSubdivisions,
		    int,
		    << "The number of subdivisions per patch, " << arg1
		    << ", is not valid.");
    
  private:
				     /**
				      * Build one patch. This function
				      * is called in a WorkStream
				      * context.
				      */
    void build_one_patch (const std::pair<cell_iterator, unsigned int> *cell_and_index,
			  internal::DataOut::ParallelData<DH::dimension, DH::space_dimension> &data,
			  DataOutBase::Patch<DH::dimension, DH::space_dimension> &patch,
			  const CurvedCellRegion curved_cell_region);
};



// -------------------- template and inline functions ------------------------


namespace internal
{
  namespace DataOut
  {
    template <int dim, int spacedim>
    template <class FE>
    ParallelDataBase<dim,spacedim>::
    ParallelDataBase (const unsigned int n_components,
		      const unsigned int n_datasets,
		      const unsigned int n_subdivisions,
		      const unsigned int n_q_points,
		      const std::vector<unsigned int> &n_postprocessor_outputs,
		      const FE &finite_elements)
		    :
		    n_components (n_components),
		    n_datasets (n_datasets),
		    n_subdivisions (n_subdivisions),
		    patch_values (n_q_points),
		    patch_values_system (n_q_points),
		    patch_gradients (n_q_points),
		    patch_gradients_system (n_q_points),
		    patch_hessians (n_q_points),
		    patch_hessians_system (n_q_points),
		    postprocessed_values (n_postprocessor_outputs.size()),
		    fe_collection (finite_elements)
    {
      for (unsigned int k=0; k<n_q_points; ++k)
	{
	  patch_values_system[k].reinit(n_components);
	  patch_gradients_system[k].resize(n_components);
	  patch_hessians_system[k].resize(n_components);
	}

      for (unsigned int dataset=0; dataset<n_postprocessor_outputs.size(); ++dataset)
	if (n_postprocessor_outputs[dataset] != 0)
	  postprocessed_values[dataset]
	    .resize(n_q_points,
		    dealii::Vector<double>(n_postprocessor_outputs[dataset]));
    }
  }
}


template <class DH, int patch_dim, int patch_space_dim>
template <class DH2>
void
DataOut_DoFData<DH,patch_dim,patch_space_dim>::
merge_patches (const DataOut_DoFData<DH2,patch_dim,patch_space_dim> &source,
	       const Point<patch_space_dim> &shift)
{
  const std::vector<Patch> source_patches = source.get_patches ();
  Assert (patches.size () != 0,        ExcNoPatches ());
  Assert (source_patches.size () != 0, ExcNoPatches ());
                                   // check equality of component
                                   // names
  Assert (get_dataset_names() == source.get_dataset_names(),
          ExcIncompatibleDatasetNames());
                                   // make sure patches are compatible. we'll
                                   // assume that if the first respective
                                   // patches are ok that all the other ones
                                   // are ok as well
  Assert (patches[0].n_subdivisions == source_patches[0].n_subdivisions,
          ExcIncompatiblePatchLists());
  Assert (patches[0].data.n_rows() == source_patches[0].data.n_rows(),
          ExcIncompatiblePatchLists());
  Assert (patches[0].data.n_cols() == source_patches[0].data.n_cols(),
          ExcIncompatiblePatchLists());

				   // check equality of the vector data
				   // specifications
  Assert (get_vector_data_ranges().size() ==
	  source.get_vector_data_ranges().size(),
	  ExcMessage ("Both sources need to declare the same components "
		      "as vectors."));
  for (unsigned int i=0; i<get_vector_data_ranges().size(); ++i)
    {
      Assert (get_vector_data_ranges()[i].template get<0>() ==
	      source.get_vector_data_ranges()[i].template get<0>(),
	      ExcMessage ("Both sources need to declare the same components "
			  "as vectors."));
      Assert (get_vector_data_ranges()[i].template get<1>() ==
	      source.get_vector_data_ranges()[i].template get<1>(),
	      ExcMessage ("Both sources need to declare the same components "
			  "as vectors."));
      Assert (get_vector_data_ranges()[i].template get<2>() ==
	      source.get_vector_data_ranges()[i].template get<2>(),
	      ExcMessage ("Both sources need to declare the same components "
			  "as vectors."));
    }
  
                                   // merge patches. store old number
                                   // of elements, since we need to
                                   // adjust patch numbers, etc
                                   // afterwards
  const unsigned int old_n_patches = patches.size();
  patches.insert (patches.end(),
                  source_patches.begin(),
                  source_patches.end());

				   // perform shift, if so desired
  if (shift != Point<patch_space_dim>())
    for (unsigned int i=old_n_patches; i<patches.size(); ++i)
      for (unsigned int v=0; v<GeometryInfo<patch_dim>::vertices_per_cell; ++v)
	patches[i].vertices[v] += shift;
  
                                   // adjust patch numbers
  for (unsigned int i=old_n_patches; i<patches.size(); ++i)
    patches[i].patch_index += old_n_patches;
  
                                   // adjust patch neighbors
  for (unsigned int i=old_n_patches; i<patches.size(); ++i)
    for (unsigned int n=0; n<GeometryInfo<patch_dim>::faces_per_cell; ++n)
      if (patches[i].neighbors[n] != Patch::no_neighbor)
        patches[i].neighbors[n] += old_n_patches;
}


DEAL_II_NAMESPACE_CLOSE

#endif
