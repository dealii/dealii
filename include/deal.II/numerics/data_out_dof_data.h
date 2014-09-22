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

#ifndef __deal2__data_out_dof_data_h
#define __deal2__data_out_dof_data_h



#include <deal.II/base/config.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/data_out_base.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/grid/tria.h>
#include <deal.II/fe/mapping.h>
#include <deal.II/hp/q_collection.h>
#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/fe_values.h>
#include <deal.II/numerics/data_postprocessor.h>
#include <deal.II/numerics/data_component_interpretation.h>

#include <deal.II/base/std_cxx11/shared_ptr.h>

DEAL_II_NAMESPACE_OPEN

template <int, int> class FEValuesBase;

namespace internal
{
  namespace DataOut
  {
    /**
     * For each vector that has been added through the add_data_vector()
     * functions, we need to keep track of a pointer to it, and allow data
     * extraction from it when we generate patches. Unfortunately, we need to
     * do this for a number of different vector types. Fortunately, they all
     * have the same interface. So the way we go is to have a base class that
     * provides the functions to access the vector's information, and to have
     * a derived template class that can be instantiated for each vector
     * type. Since the vectors all have the same interface, this is no big
     * problem, as they can all use the same general templatized code.
     *
     * @author Wolfgang Bangerth, 2004
     */
    template <class DH>
    class DataEntryBase
    {
    public:
      /**
       * Constructor. Give a list of names for the individual components of
       * the vector and their interpretation as scalar or vector data. This
       * constructor assumes that no postprocessor is going to be used.
       */
      DataEntryBase (const DH                       *dofs,
                     const std::vector<std::string> &names,
                     const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation);

      /**
       * Constructor when a data postprocessor is going to be used. In that
       * case, the names and vector declarations are going to be acquired from
       * the postprocessor.
       */
      DataEntryBase (const DH                       *dofs,
                     const DataPostprocessor<DH::space_dimension> *data_postprocessor);

      /**
       * Destructor made virtual.
       */
      virtual ~DataEntryBase ();

      /**
       * Assuming that the stored vector is a cell vector, extract the given
       * element from it.
       */
      virtual
      double
      get_cell_data_value (const unsigned int cell_number) const = 0;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store.
       */
      virtual
      void
      get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<double>             &patch_values) const = 0;

      /**
       * Given a FEValuesBase object, extract the values on the present cell
       * from the vector we actually store. This function does the same as the
       * one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_values (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                           std::vector<dealii::Vector<double> >    &patch_values_system) const = 0;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store.
       */
      virtual
      void
      get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                              std::vector<Tensor<1,DH::space_dimension> >       &patch_gradients) const = 0;

      /**
       * Given a FEValuesBase object, extract the gradients on the present
       * cell from the vector we actually store. This function does the same
       * as the one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_gradients (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                              std::vector<std::vector<Tensor<1,DH::space_dimension> > > &patch_gradients_system) const = 0;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store.
       */
      virtual
      void
      get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<Tensor<2,DH::space_dimension> >       &patch_hessians) const = 0;

      /**
       * Given a FEValuesBase object, extract the second derivatives on the
       * present cell from the vector we actually store. This function does
       * the same as the one above but for vector-valued finite elements.
       */
      virtual
      void
      get_function_hessians (const FEValuesBase<DH::dimension,DH::space_dimension> &fe_patch_values,
                             std::vector<std::vector< Tensor<2,DH::space_dimension> > > &patch_hessians_system) const = 0;

      /**
       * Clear all references to the vectors.
       */
      virtual void clear () = 0;

      /**
       * Determine an estimate for the memory consumption (in bytes) of this
       * object.
       */
      virtual std::size_t memory_consumption () const = 0;

      /**
       * Pointer to the DoFHandler object that the vector is based on.
       */
      SmartPointer<const DH> dof_handler;

      /**
       * Names of the components of this data vector.
       */
      const std::vector<std::string> names;

      /**
       * A vector that for each of the n_output_variables variables of the
       * current data set indicates whether they are scalar fields, parts of a
       * vector-field, or any of the other supported kinds of data.
       */
      const std::vector<DataComponentInterpretation::DataComponentInterpretation>
      data_component_interpretation;

      /**
       * Pointer to a DataPostprocessing object which shall be applied to this
       * data vector.
       */
      SmartPointer<const dealii::DataPostprocessor<DH::space_dimension> > postprocessor;

      /**
       * Number of output variables this dataset provides (either number of
       * components in vector valued function / data vector or number of
       * computed quantities, if DataPostprocessor is applied). This variable
       * is determined via and thus equivalent to <tt>names.size()</tt>.
       */
      unsigned int n_output_variables;
    };


    /**
     * A data structure that holds all data needed in one thread when building
     * patches in parallel. These data structures are created globally rather
     * than on each cell to avoid allocation of memory in the threads. This is
     * a base class for the AdditionalData kind of data structure discussed in
     * the documentation of the WorkStream class.
     *
     * The <code>cell_to_patch_index_map</code> is an array that stores for
     * index <tt>[i][j]</tt> the number of the patch that associated with the
     * cell with index @p j on level @p i. This information is set up prior to
     * generation of the patches, and is needed to generate neighborship
     * information.
     *
     * This structure is used by several of the DataOut* classes, which
     * derived their own ParallelData classes from it for additional fields.
     */
    template <int dim, int spacedim>
    struct ParallelDataBase
    {
      ParallelDataBase (const unsigned int n_datasets,
                        const unsigned int n_subdivisions,
                        const std::vector<unsigned int> &n_postprocessor_outputs,
                        const Mapping<dim,spacedim> &mapping,
                        const std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > &finite_elements,
                        const UpdateFlags update_flags,
                        const bool        use_face_values);

      ParallelDataBase (const ParallelDataBase &data);

      template <typename DH>
      void reinit_all_fe_values(std::vector<std_cxx11::shared_ptr<DataEntryBase<DH> > > &dof_data,
                                const typename dealii::Triangulation<dim,spacedim>::cell_iterator &cell,
                                const unsigned int face = numbers::invalid_unsigned_int);

      const FEValuesBase<dim,spacedim> &
      get_present_fe_values(const unsigned int dataset) const;

      void resize_system_vectors(const unsigned int n_components);

      const unsigned int n_datasets;
      const unsigned int n_subdivisions;

      std::vector<double>                                patch_values;
      std::vector<dealii::Vector<double> >               patch_values_system;
      std::vector<Tensor<1,spacedim> >                   patch_gradients;
      std::vector<std::vector<Tensor<1,spacedim> > >     patch_gradients_system;
      std::vector<Tensor<2,spacedim> >                   patch_hessians;
      std::vector<std::vector<Tensor<2,spacedim> > >     patch_hessians_system;
      std::vector<std::vector<dealii::Vector<double> > > postprocessed_values;

      const dealii::hp::MappingCollection<dim,spacedim> mapping_collection;
      const std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<dim,spacedim> > > finite_elements;
      const UpdateFlags update_flags;

      std::vector<std_cxx11::shared_ptr<dealii::hp::FEValues<dim,spacedim> > > x_fe_values;
      std::vector<std_cxx11::shared_ptr<dealii::hp::FEFaceValues<dim,spacedim> > > x_fe_face_values;
    };
  }
}


//TODO: Most of the documentation of DataOut_DoFData applies to DataOut.

/**
 * This is an abstract class which provides the functionality to generate
 * patches for output by base classes from data vectors on a grid. It allows
 * to store one or more pointers to a DoFHandler and attached node and cell
 * data denoting functions on the grid which shall later be written in any of
 * the implemented data formats.
 *
 *
 * <h3>User visible interface</h3>
 *
 * The user visible interface of this class allows the user to specify data in
 * two different ways. One is to make a DoFHandler object known to this class
 * and to add data vectors that all correspond to this DoFHandler or the grid
 * cells which will later be written to a file in some format. The second
 * approach is to pass a DoFHandler object along with the vector. This allows
 * setting data from different DoFHandlers in a neat way (of course, they both
 * need to be based on the same triangulation). Instead of pondering about the
 * different functions, an example for the first kind is probably the best
 * explanation:
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
 * information from a data vector at each point at which the field will
 * be evaluated so that it can be written to a file (for example, the
 * Mach number in hypersonic flow can be computed from density and velocities;
 * step-29 also shows an example); another piece of information
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
 * Please note that in the example above, an object of type DataOut was
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
   * Typedef to the iterator type of the dof handler class under
   * consideration.
   */
  typedef typename Triangulation<DH::dimension,DH::space_dimension>::cell_iterator cell_iterator;
  typedef typename Triangulation<DH::dimension,DH::space_dimension>::active_cell_iterator active_cell_iterator;

public:

  /**
   * Type describing what the vector given to add_data_vector() is: a vector
   * that has one entry per degree of freedom in a DoFHandler object (such as
   * solution vectors), or one entry per cell in the triangulation underlying
   * the DoFHandler object (such as error per cell data). The value
   * #type_automatic tells add_data_vector() to find out itself (see the
   * documentation of add_data_vector() for the method used).
   */
  enum DataVectorType
  {
    /**
     * Data vector entries are associated to degrees of freedom
     */
    type_dof_data,

    /**
     * Data vector entries are one per grid cell
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
   * Designate a dof handler to be used to extract geometry data and the
   * mapping between nodes and node values. This call is not necessary if all
   * added data vectors are supplemented with a DoFHandler argument.
   *
   * This call is optional: If you add data vectors with specified DoFHandler
   * object, then that contains all information needed to generate the output.
   */
  void attach_dof_handler (const DH &);

  /**
   * Designate a triangulation to be used to extract geometry data and the
   * mapping between nodes and node values.
   *
   * This call is optional: If you add data vectors with specified DoFHandler
   * object, then that contains all information needed to generate the
   * output. This call is useful when you only output cell vectors and no
   * DoFHandler at all, in which case it provides the geometry.
   */
  void attach_triangulation (const Triangulation<DH::dimension,
                             DH::space_dimension> &);

  /**
   * Add a data vector together with its name.
   *
   * A pointer to the vector is stored, so you have to make sure the vector
   * exists at that address at least as long as you call the <tt>write_*</tt>
   * functions.
   *
   * It is assumed that the vector has the same number of components as there
   * are degrees of freedom in the dof handler, in which case it is assumed to
   * be a vector storing nodal data; or the size may be the number of active
   * cells on the present grid, in which case it is assumed to be a cell data
   * vector. As the number of degrees of freedom and of cells is usually not
   * equal, the function can determine itself which type of vector it is
   * given. However, there are corner cases where this automatic determination
   * does not work.  One example is if you compute with piecewise constant
   * elements and have a scalar solution, then there are as many cells as
   * there are degrees of freedom (though they may be numbered
   * differently). Another possibility is if you have a 1d mesh embedded in 2d
   * space and the mesh consists of a closed curve of cells; in this case,
   * there are as many nodes as there are cells, and when using a Q1 element
   * you will have as many degrees of freedom as there are cells.  In these
   * cases, you can change the last argument of the function from its default
   * value #type_automatic to either #type_dof_data or #type_cell_data,
   * depending on what the vector represents. Apart from such corner cases,
   * you can leave the argument at its default value and let the function
   * determine the type of the vector itself.
   *
   * If it is a vector holding DoF data, the names given shall be one for each
   * component of the underlying finite element.  If it is a finite element
   * composed of only one subelement, then there is another function following
   * which takes a single name instead of a vector of names.
   *
   * The data_component_interpretation argument contains information about how
   * the individual components of output files that consist of more than one
   * data set are to be interpreted.
   *
   * For example, if one has a finite element for the Stokes equations in 2d,
   * representing components (u,v,p), one would like to indicate that the
   * first two, u and v, represent a logical vector so that later on when we
   * generate graphical output we can hand them off to a visualization program
   * that will automatically know to render them as a vector field, rather
   * than as two separate and independent scalar fields.
   *
   * The default value of this argument (i.e. an empty vector) corresponds is
   * equivalent to a vector of values
   * DataComponentInterpretation::component_is_scalar, indicating that all
   * output components are independent scalar fields. However, if the given
   * data vector represents logical vectors, you may pass a vector that
   * contains values
   * DataComponentInterpretation::component_is_part_of_vector. In the example
   * above, one would pass in a vector with components
   * (DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_part_of_vector,
   * DataComponentInterpretation::component_is_scalar) for (u,v,p).
   *
   * The names of a data vector shall only contain characters which are
   * letters, underscore and a few other ones. Refer to the
   * ExcInvalidCharacter exception declared in this class to see which
   * characters are valid and which are not.
   *
   * @note The actual type for the vector argument may be any vector type from
   * which FEValues can extract values on a cell using the
   * FEValuesBase::get_function_values() function. In particular, this
   * includes all of the usual vector types, but also IndexSet (see step-41
   * for a use of this).
   */
  template <class VECTOR>
  void add_data_vector (const VECTOR                   &data,
                        const std::vector<std::string> &names,
                        const DataVectorType            type = type_automatic,
                        const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
                        = std::vector<DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * This function is an abbreviation to the above one (see there for a
   * discussion of the various arguments), intended for use with finite
   * elements that are not composed of subelements. In this case, only one
   * name per data vector needs to be given, which is what this function
   * takes. It simply relays its arguments after a conversion of the @p name
   * to a vector of strings, to the other add_data_vector() function above.
   *
   * If @p data is a vector with multiple components this function will
   * generate distinct names for all components by appending an underscore and
   * the number of each component to @p name
   *
   * The actual type for the template argument may be any vector type from
   * which FEValues can extract values on a cell using the
   * FEValuesBase::get_function_values() function.
   */
  template <class VECTOR>
  void add_data_vector (const VECTOR         &data,
                        const std::string    &name,
                        const DataVectorType  type = type_automatic,
                        const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
                        = std::vector<DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * This function is an extension of the above one (see there for a
   * discussion of the arguments except the first one) and allows to set a
   * vector with its own DoFHandler object. This DoFHandler needs to be
   * compatible with the other DoFHandler objects assigned with calls to @p
   * add_data_vector or @p attach_dof_handler, in the sense that the
   * DoFHandler objects need to be based on the same triangulation. This
   * function allows you to export data from multiple DoFHandlers that
   * describe different solution components.
   *
   * Since this function takes a DoFHandler object and hence naturally
   * represents dof data, the data vector type argument present in the other
   * methods above is skipped.
   */
  template <class VECTOR>
  void add_data_vector (const DH                       &dof_handler,
                        const VECTOR                   &data,
                        const std::vector<std::string> &names,
                        const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
                        = std::vector<DataComponentInterpretation::DataComponentInterpretation>());


  /**
   * This function is an abbreviation of the function above with only a scalar
   * @p dof_handler given and a single data name.
   */
  template <class VECTOR>
  void add_data_vector (const DH                       &dof_handler,
                        const VECTOR                   &data,
                        const std::string              &name,
                        const std::vector<DataComponentInterpretation::DataComponentInterpretation> &data_component_interpretation
                        = std::vector<DataComponentInterpretation::DataComponentInterpretation>());

  /**
   * This function is an alternative to the above ones, allowing the output of
   * derived quantities instead of the given data. This converison has to be
   * done in a class derived from DataPostprocessor.
   *
   * The names for these derived quantities are provided by the @p
   * data_postprocessor argument. Likewise, the data_component_interpretation
   * argument of the other add_data_vector() functions is provided by the
   * data_postprocessor argument. As only data of type @p type_dof_data can be
   * transformed, this type is also known implicitly and does not have to be
   * given.
   *
   * @note The actual type for the vector argument may be any vector type from
   * which FEValues can extract values on a cell using the
   * FEValuesBase::get_function_values() function. In particular, this
   * includes all of the usual vector types, but also IndexSet (see step-41
   * for a use of this).
   */
  template <class VECTOR>
  void add_data_vector (const VECTOR                           &data,
                        const DataPostprocessor<DH::space_dimension> &data_postprocessor);

  /**
   * Same function as above, but with a DoFHandler object that does not need
   * to coincide with the DoFHandler initially set. Note that the
   * postprocessor can only read data from the given DoFHandler and solution
   * vector, not other solution vectors or DoFHandlers.
   */
  template <class VECTOR>
  void add_data_vector (const DH                               &dof_handler,
                        const VECTOR                           &data,
                        const DataPostprocessor<DH::space_dimension> &data_postprocessor);

  /**
   * Release the pointers to the data vectors. This allows output of a new set
   * of vectors without supplying the DoF handler again. Therefore, the
   * DataOut object can be used in an algebraic context. Note that besides the
   * data vectors also the patches already computed are deleted.
   */
  void clear_data_vectors ();

  /**
   * Release pointers to all input data elements, i.e. pointers to data
   * vectors and to the DoF handler object. This function may be useful when
   * you have called the @p build_patches function of derived class, since
   * then the patches are built and the input data is no more needed, nor is
   * there a need to reference it. You can then output the patches detached
   * from the main thread and need not make sure anymore that the DoF handler
   * object and vectors must not be deleted before the output thread is
   * finished.
   */
  void clear_input_data_references ();

  /**
   * This function can be used to merge the patches that were created using
   * the @p build_patches function of the object given as argument into the
   * list of patches created by this object. This is sometimes handy if one
   * has, for example, a domain decomposition algorithm where each block is
   * represented by a DoFHandler of its own, but one wants to output the
   * solution on all the blocks at the same time.
   *
   * For this to work, the given argument and this object need to have the
   * same number of output vectors, and they need to use the same number of
   * subdivisions per patch. The output will probably look rather funny if
   * patches in both objects overlap in space.
   *
   * If you call build_patches() for this object after merging in patches, the
   * previous state is overwritten, and the merged-in patches are lost.
   *
   * The second parameter allows to shift each node of the patches in the
   * object passed in in the first parameter by a certain amount. This is
   * sometimes useful to generate "exploded" views of a collection of blocks.
   *
   * This function will fail if either this or the other object did not yet
   * set up any patches.
   */
  template <class DH2>
  void merge_patches (const DataOut_DoFData<DH2,patch_dim,patch_space_dim> &source,
                      const Point<patch_space_dim> &shift = Point<patch_space_dim>());

  /**
   * Release the pointers to the data vectors and the DoF handler. You have to
   * set all data entries again using the add_data_vector() function. The
   * pointer to the dof handler is cleared as well, along with all other
   * data. In effect, this function resets everything to a virgin state.
   */
  virtual void clear ();

  /**
   * Determine an estimate for the memory consumption (in bytes) of this
   * object.
   */
  std::size_t memory_consumption () const;

  /**
   * Exception
   */
  DeclException0 (ExcNoTriangulationSelected);

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
   * Abbreviate the somewhat lengthy name for the Patch class.
   */
  typedef dealii::DataOutBase::Patch<patch_dim,patch_space_dim> Patch;

  /**
   * Pointer to the triangulation object.
   */
  SmartPointer<const Triangulation<DH::dimension,DH::space_dimension> > triangulation;

  /**
   * Pointer to the optional handler object.
   */
  SmartPointer<const DH> dofs;

  /**
   * List of data elements with vectors of values for each degree of freedom.
   */
  std::vector<std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> > >  dof_data;

  /**
   * List of data elements with vectors of values for each cell.
   */
  std::vector<std_cxx11::shared_ptr<internal::DataOut::DataEntryBase<DH> > >  cell_data;

  /**
   * This is a list of patches that is created each time build_patches() is
   * called. These patches are used in the output routines of the base
   * classes.
   */
  std::vector<Patch> patches;

  /**
   * Function by which the base class's functions get to know what patches
   * they shall write to a file.
   */
  virtual
  const std::vector<Patch> &get_patches () const;

  /**
   * Virtual function through which the names of data sets are obtained by the
   * output functions of the base class.
   */
  virtual
  std::vector<std::string> get_dataset_names () const;

  /**
   * Extracts the finite elements stored in the dof_data object, including a
   * dummy object of FE_DGQ<dim>(0) in case only the triangulation is used.
   */
  std::vector<std_cxx11::shared_ptr<dealii::hp::FECollection<DH::dimension,DH::space_dimension> > >
  get_finite_elements() const;

  /**
   * Overload of the respective DataOutInterface::get_vector_data_ranges()
   * function. See there for a more extensive documentation.
   */
  virtual
  std::vector<std_cxx11::tuple<unsigned int, unsigned int, std::string> >
  get_vector_data_ranges () const;

  /**
   * Make all template siblings friends. Needed for the merge_patches()
   * function.
   */
  template <class, int, int>
  friend class DataOut_DoFData;
};



// -------------------- template and inline functions ------------------------


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
      Assert (std_cxx11::get<0>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<0>(source.get_vector_data_ranges()[i]),
              ExcMessage ("Both sources need to declare the same components "
                          "as vectors."));
      Assert (std_cxx11::get<1>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<1>(source.get_vector_data_ranges()[i]),
              ExcMessage ("Both sources need to declare the same components "
                          "as vectors."));
      Assert (std_cxx11::get<2>(get_vector_data_ranges()[i]) ==
              std_cxx11::get<2>(source.get_vector_data_ranges()[i]),
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
