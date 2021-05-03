// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

#ifndef dealii_meshworker_scratch_data_h
#define dealii_meshworker_scratch_data_h

#include <deal.II/base/config.h>

#include <deal.II/algorithms/general_data_storage.h>

#include <deal.II/differentiation/ad.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <boost/any.hpp>

#include <algorithm>
#include <map>
#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace MeshWorker
{
  /**
   * A helper class to simplify the parallel assembly of linear and non-linear
   * problems, and the evaluation of finite element fields.
   *
   * This class is a drop in ScratchData to use with the WorkStream::run()
   * function and with the MeshWorker::mesh_loop function().
   *
   * The ScratchData class has three main goals:
   * - to create FEValues, FEFaceValues, FESubfaceValues, and FEInterfaceValues
   *   for the current cell and for its neighbor cell *on demand* (only if they
   *   are necessary for the algorithm provided by the user), and to provide a
   *   uniform interface to access the FEValues objects when assembling cell,
   *   face, or subface contributions
   * - to store arbitrary data types (or references to arbitrary data types),
   *   that can be retrieved by name (for example, when the assembly needs to
   *   make reference to other objects such as previous time steps' solution
   *   vectors, previous nonlinear iteration vectors, geometry descriptions,
   *   etc.) in an object of type
   * - to provide a reasonable interface for those use cases where the user may
   *   need temporary vectors of data at quadrature points, allowing the
   *   construction of these temporaries *on demand*, and easy access to
   *   values, gradients, etc., of already computed solution vectors.
   *
   * The methods in the section "Methods to work on current cell"
   * initialize on demand internal FEValues, FEFaceValues, FESubfaceValues, and
   * FEInterfaceValues objects on the current cell, allowing the use of this
   * class as a single substitute for four different objects used to integrate
   * and query finite element values on cells, faces, and subfaces.
   *
   * Similarly, the methods in the section "Methods to work on neighbor cell"
   * initialize on demand
   * (different) internal FEValues, FEFaceValues, and FESubfaceValues, on
   * neighbor cells, allowing the use of this class also as a single substitute
   * for the additional three objects you would typically need to integrate on
   * the neighbor cell, and on its faces and subfaces (for example, in
   * discontinuous Galerkin methods).
   *
   * If you need to retrieve values or gradients of finite element solution
   * vectors, on the cell, face, or subface that has just been initialized
   * with one of the functions in the section "Methods to work on current cell",
   * you can use the methods in the section "Evaluation of finite element fields
   * and their derivatives on the current cell".
   *
   * An example usage for this class is given by the following snippet of code:
   *
   * @code
   * ScratchData<dim,spacedim> scratch(fe,
   *                                   cell_quadrature,
   *                                   update_values | update_gradients);
   *
   * FEValuesExtractors::Vector velocity(0);
   * FEValuesExtractors::Scalar pressure(dim);
   *
   * ...
   *
   * for(const auto &cell : dof_handler.active_cell_iterators())
   * {
   *    const auto &fe_values = scratch.reinit(cell);
   *    scratch.extract_local_dof_values("current_solution",
   *                                     current_solution);
   *
   *    scratch.extract_local_dof_values("previous_solution",
   *                                     previous_solution);
   *
   *    const auto &local_solution_values =
   *          scratch.get_local_dof_values("current_solution",
   *                                       current_solution);
   *
   *    const auto &local_previous_solution_values =
   *          scratch.get_local_dof_values("previous_solution",
   *                                       previous_solution);
   *
   *    const auto &previous_solution_values_velocity =
   *          scratch.get_values("previous_solution", velocity);
   *    const auto &previous_solution_values_pressure =
   *          scratch.get_values("previous_solution", pressure);
   *
   *    const auto &solution_values_velocity =
   *          scratch.get_values("solution", velocity);
   *    const auto &solution_values_pressure =
   *          scratch.get_values("solution", pressure);
   *
   *    const auto &solution_symmetric_gradients_velocity =
   *          scratch.get_symmetric_gradients("solution", velocity);
   *
   *    // Do something with the above.
   * }
   * @endcode
   *
   * The order in which you call functions of this class matters: if you call
   * the ScratchData::reinit() function that takes an active cell iterator,
   * then subsequent calls to methods that internally need an FEValuesBase
   * object will use the internal FEValues object initialized with the given
   * cell to perform their calculations. On the other hand, if you have called
   * the ScratchData::reinit() method that also takes a face index, all
   * subsequent calls to methods that need an FEValuesBase object, will use an
   * internally stored FEFaceValues object, initialized with the cell and face
   * index passed to the ScratchData::reinit() function. The same applies for
   * the ScratchData::reinit() method that takes three arguments: the cell, the
   * face index, and the subface index.
   *
   * The user code should be structured without interleaving work on cells and
   * work on faces.
   *
   * Consider, for example, the following snippet of code:
   *
   * @code
   * ScratchData<dim,spacedim> scratch(...);
   * FEValuesExtractors::Scalar temperature(0);
   *
   * for(const auto &cell : dof_handler.active_cell_iterators())
   * {
   *   const auto &fe_values = scratch.reinit(cell);
   *   const auto &local_dof_values =
   *         scratch.extract_local_dof_values("solution", solution);
   *
   *   // This will contain all values of the temperature on the cell
   *   // quadrature points
   *   const auto &solution_values_cell =
   *         scratch.get_values("solution", temperature);
   *
   *   // Do something with values on the cell
   *   ...
   *
   *   // Now start working on the faces
   *   for(unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
   *   {
   *      const auto &fe_face_values = scratch.reinit(cell, f);
   *
   *      // Now we'll have access to the values of the temperature on the faces
   *      const auto &solution_values_face =
   *            scratch.get_values("solution", temperature);
   *
   *      // Notice how the function call is the same, but the result depends
   *      // on what was the last reinit() function that was called. In this
   *      // case, we called reinit(cell, f), triggering internal work on
   *      // faces.
   *   }
   *
   *   // At this point, we would like to go back and work on cells,
   *   // for example querying the values of the gradients:
   *   const auto &solution_gradients =
   *         scratch.get_gradients("solution", temperature);
   *
   *   // This assertion would be triggered in debug mode
   *   AssertDimension(solution_gradients.size(), quadrature_cell.size());
   *
   *   // However, with the above call the content of solution_gradients is the
   *   // gradient of the temperature on the quadrature points of the last face
   *   // visited in the previous loop.
   *   // If you really want to have the values of the gradients on the cell,
   *   // at this point your only option is to call
   *   scratch.reinit(cell);
   *
   *   // (again) before querying for the gradients:
   *   const auto &solution_gradients_cell =
   *         scratch.get_gradients("solution", temperature);
   *
   *   // The call to reinit(cell), however, is an expensive one. You
   *   // should make sure you only call it once per cell by grouping together
   *   // all queries that concern cells in the same place (right after you
   *   // call the reinit(cell) method).
   *   // A similar argument holds for all calls on each face and subface.
   * }
   * @endcode
   *
   * When using this class, please cite
   * @code{.bib}
   * @article{SartoriGiulianiBardelloni-2018-a,
   * 	Author  = {Sartori, Alberto and Giuliani, Nicola and
   *            Bardelloni, Mauro and Heltai, Luca},
   * 	Journal = {SoftwareX},
   * 	Pages   = {318--327},
   * 	Title   = {{deal2lkit: A toolkit library for high performance
   *              programming in deal.II}},
   *  Doi     = {10.1016/j.softx.2018.09.004},
   * 	Volume  = {7},
   * 	Year    = {2018}}
   * @endcode
   */
  template <int dim, int spacedim = dim>
  class ScratchData
  {
  public:
    /**
     * Create an empty ScratchData object. A SmartPointer pointing to
     * @p mapping and @p fe is stored internally. Make sure they live longer
     * than this class instance.
     *
     * The constructor does not initialize any of the internal FEValues objects.
     * These are initialized the first time one of the reinit() functions is
     * called, using the arguments passed here.
     *
     * @param mapping The mapping to use in the internal FEValues objects
     * @param fe The finite element
     * @param quadrature The cell quadrature
     * @param update_flags UpdateFlags for the current cell FEValues and
     * neighbor cell FEValues
     * @param face_quadrature Face quadrature, used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     * @param face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     */
    ScratchData(
      const Mapping<dim, spacedim> &      mapping,
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default);

    /**
     * Similar to the other constructor, but this one allows to specify
     * different flags for neighbor cells and faces.
     *
     * @param mapping The mapping to use in the internal FEValues objects
     * @param fe The finite element
     * @param quadrature The cell quadrature
     * @param update_flags UpdateFlags for the current cell FEValues
     * @param neighbor_update_flags UpdateFlags for the neighbor cell FEValues
     * @param face_quadrature Face quadrature, used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     * @param face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for the current cell
     * @param neighbor_face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for the neighbor cell
     */
    ScratchData(
      const Mapping<dim, spacedim> &      mapping,
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const UpdateFlags &                 neighbor_update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default,
      const UpdateFlags &        neighbor_face_update_flags = update_default);

    /**
     * Same as the other constructor, using the default MappingQ1.
     *
     * @param fe The finite element
     * @param quadrature The cell quadrature
     * @param update_flags UpdateFlags for the current cell FEValues and
     * neighbor cell FEValues
     * @param face_quadrature Face quadrature, used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     * @param face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     */
    ScratchData(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default);

    /**
     * Same as the other constructor, using the default MappingQ1.
     *
     * @param fe The finite element
     * @param quadrature The cell quadrature
     * @param update_flags UpdateFlags for the current cell FEValues
     * @param neighbor_update_flags UpdateFlags for the neighbor cell FEValues
     * @param face_quadrature Face quadrature, used for FEFaceValues and
     * FESubfaceValues for both the current cell and the neighbor cell
     * @param face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for the current cell
     * @param neighbor_face_update_flags UpdateFlags used for FEFaceValues and
     * FESubfaceValues for the neighbor cell
     */
    ScratchData(
      const FiniteElement<dim, spacedim> &fe,
      const Quadrature<dim> &             quadrature,
      const UpdateFlags &                 update_flags,
      const UpdateFlags &                 neighbor_update_flags,
      const Quadrature<dim - 1> &face_quadrature   = Quadrature<dim - 1>(),
      const UpdateFlags &        face_update_flags = update_default,
      const UpdateFlags &        neighbor_face_update_flags = update_default);

    /**
     * Deep copy constructor. FEValues objects are not copied.
     */
    ScratchData(const ScratchData<dim, spacedim> &scratch);

    /**
     * @name Methods to work on current cell
     */
    /**@{*/ // CurrentCellMethods

    /**
     * Initialize the internal FEValues with the given @p cell, and return
     * a reference to it.
     *
     * After calling this function, get_current_fe_values() will return the
     * same object of this method, as an FEValuesBase reference.
     */
    const FEValues<dim, spacedim> &
    reinit(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell);

    /**
     * Initialize the internal FEFaceValues to use the given @p face_no on the given
     * @p cell, and return a reference to it.
     *
     * After calling this function, get_current_fe_values() will return the
     * same object of this method, as an FEValuesBase reference.
     */
    const FEFaceValues<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no);

    /**
     * Initialize the internal FESubfaceValues to use the given @p subface_no,
     * on @p face_no, on the given @p cell, and return a reference to it.
     *
     * After calling this function, get_current_fe_values() will return the
     * same object of this method, as an FEValuesBase reference.
     *
     * If @p subface_no is numbers::invalid_unsigned_int, the reinit() function
     * that takes only the @p cell and the @p face_no is called.
     */
    const FEFaceValuesBase<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no,
           const unsigned int subface_no);

    /**
     * Initialize the internal FEInterfaceValues with the given arguments, and
     * return a reference to it.
     *
     * After calling this function, get_local_dof_indices(),
     * get_quadrature_points(), get_normal_vectors(), and get_JxW_values() will
     * be forwarded to the local FEInterfaceValues object. The methods
     * get_current_fe_values() will return the FEValuesBase associated to the
     * current cell, while get_neighbor_fe_values() will be associated with the
     * neighbor cell. The method get_local_dof_indices() will return the
     * same result of FEInterfaceValues::get_interface_dof_to_dof_indices(),
     * while the get_neighbor_dof_indices() will return the local dof indices
     * of the neighbor cell.
     */
    const FEInterfaceValues<dim, spacedim> &
    reinit(const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
           const unsigned int face_no,
           const unsigned int sub_face_no,
           const typename DoFHandler<dim, spacedim>::active_cell_iterator
             &                cell_neighbor,
           const unsigned int face_no_neighbor,
           const unsigned int sub_face_no_neighbor);

    /**
     * Get the currently initialized FEValues.
     *
     * This function will return the internal FEValues if the
     * reinit(cell) function was called last. If the reinit(cell, face_no)
     * function was called, then this function returns the internal
     * FEFaceValues, and if the reinit(cell, face_no, subface_no) function was
     * called (with a valid @p subface_no argument), it returns the internal
     * FESubfaceValues object.
     */
    const FEValuesBase<dim, spacedim> &
    get_current_fe_values() const;

    /**
     * Return the quadrature points of the internal FEValues object.
     */
    const std::vector<Point<spacedim>> &
    get_quadrature_points() const;

    /**
     * Return the JxW values of the internal FEValues object.
     */
    const std::vector<double> &
    get_JxW_values() const;

    /**
     * Return the last computed normal vectors.
     */
    const std::vector<Tensor<1, spacedim>> &
    get_normal_vectors() const;

    /**
     * Return the local dof indices for the cell passed the last time the
     * reinit() function was called.
     */
    const std::vector<types::global_dof_index> &
    get_local_dof_indices() const;

    /** @} */ // CurrentCellMethods

    /**
     * @name Methods to work on neighbor cell
     */
    /** @{ */ // NeighborCellMethods

    /**
     * Initialize the internal neighbor FEValues to use the given @p cell, and
     * return a reference to it.
     *
     * After calling this function, get_current_neighbor_fe_values() will return
     * the same object of this method, as an FEValuesBase reference.
     */
    const FEValues<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell);

    /**
     * Initialize the internal FEFaceValues to use the given @p face_no on the
     * given @p cell, and return a reference to it.
     *
     * After calling this function, get_current_neighbor_fe_values() will return
     * the same object of this method, as an FEValuesBase reference.
     */
    const FEFaceValues<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int                                              face_no);

    /**
     * Initialize the internal FESubfaceValues to use the given @p subface_no,
     * on @p face_no, on the given @p cell, and return a reference to it.
     *
     * After calling this function, get_current_neighbor_fe_values() will return
     * the same object of this method, as an FEValuesBase reference.
     *
     * If @p subface_no is numbers::invalid_unsigned_int, the reinit() function
     * that takes only the @p cell and the @p face_no is called.
     */
    const FEFaceValuesBase<dim, spacedim> &
    reinit_neighbor(
      const typename DoFHandler<dim, spacedim>::active_cell_iterator &cell,
      const unsigned int                                              face_no,
      const unsigned int subface_no);

    /**
     * Get the currently initialized neighbor FEValues.
     *
     * This function will return the neighbor FEValues if the
     * reinit_neighbor(cell) function was called last. If the
     * reinit_neighbor(cell, face_no) function was called, then this function
     * returns the internal neighbor FEFaceValues, and if the
     * reinit_neighbor(cell, face_no, subface_no) function was
     * called (with a valid @p subface_no argument), it returns the internal neighbor
     * FESubfaceValues object.
     */
    const FEValuesBase<dim, spacedim> &
    get_current_neighbor_fe_values() const;

    /**
     * Return the JxW values of the neighbor FEValues object.
     */
    const std::vector<double> &
    get_neighbor_JxW_values() const;

    /**
     * Return the last computed normal vectors on the neighbor.
     */
    const std::vector<Tensor<1, spacedim>> &
    get_neighbor_normal_vectors();

    /**
     * Return the local dof indices of the neighbor passed the last time the
     * reinit_neighbor() function was called.
     */
    const std::vector<types::global_dof_index> &
    get_neighbor_dof_indices() const;

    /** @} */ // NeighborCellMethods

    /**
     * Return a GeneralDataStorage object that can be used to store any amount
     * of data, of any type, which is then made accessible by an identifier
     * string.
     */
    GeneralDataStorage &
    get_general_data_storage();

    /**
     * Return a GeneralDataStorage object that can be used to store any amount
     * of data, of any type, which is then made accessible by an identifier
     * string.
     */
    const GeneralDataStorage &
    get_general_data_storage() const;

    /**
     * @name Evaluation of finite element fields and their derivatives on the current cell
     */
    /** @{ */ // CurrentCellEvaluation

    /**
     * Extract the local dof values associated with the internally initialized
     * cell.
     *
     * Before you call this function, you have to make sure you have previously
     * called one of the reinit() functions.
     *
     * At every call of this function, a new vector of dof values is generated
     * and stored internally, unless a previous vector with the same name is
     * found. If this is the case, the content of that vector is overwritten.
     *
     * If you give a unique @p global_vector_name, then for each cell you are
     * guaranteed to get a unique vector of independent dofs of the same type
     * as the dummy variable. If you use an automatic differentiation number
     * type (like Sacado::Fad::DFad<double>,
     * Sacado::Fad::DFad<Sacado::Fad::DFad<double>>, etc.) this method will
     * also initialize the independent variables internally, allowing you
     * to perform automatic differentiation.
     *
     * You can access the extracted local dof values by calling the method
     * get_local_dof_values() with the same @p global_vector_name argument
     * you passed here.
     *
     * Notice that using this initialization strategy renders the use of this
     * ScratchData object incompatible with the AD helper classes (since they
     * do their own data management). In particular, it is necessary for the
     * user to manage all of the AD data themselves (both before and after this
     * call).
     */
    template <typename VectorType, typename Number = double>
    void
    extract_local_dof_values(const std::string &global_vector_name,
                             const VectorType & input_vector,
                             const Number       dummy = Number(0));

    /**
     * After calling extract_local_dof_values(), you can retrieve the stored
     * information through this method.
     *
     * Both the argument @p global_vector_name and the type of the @p dummy
     * variable should match the ones you passed to the
     * extract_local_dof_values() function.
     */
    template <typename Number = double>
    const std::vector<Number> &
    get_local_dof_values(const std::string &global_vector_name,
                         Number             dummy = Number(0)) const;

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the values of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_values flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_value_type<Number>> &
    get_values(const std::string &global_vector_name,
               const Extractor &  variable,
               const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the gradients of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_gradients flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_gradient_type<Number>> &
    get_gradients(const std::string &global_vector_name,
                  const Extractor &  variable,
                  const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the symmetric gradients of the function at the quadrature points, and
     * return a vector with the correct type deduced by the Extractor you passed
     * as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_gradients flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_symmetric_gradient_type<Number>> &
    get_symmetric_gradients(const std::string &global_vector_name,
                            const Extractor &  variable,
                            const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the divergences of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_gradients flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_divergence_type<Number>> &
    get_divergences(const std::string &global_vector_name,
                    const Extractor &  variable,
                    const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the curls of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_gradients flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_curl_type<Number>> &
    get_curls(const std::string &global_vector_name,
              const Extractor &  variable,
              const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the hessians of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_hessians flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_hessian_type<Number>> &
    get_hessians(const std::string &global_vector_name,
                 const Extractor &  variable,
                 const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the Laplacians of the function at the quadrature points, and return a
     * vector with the correct type deduced by the Extractor you passed as the
     * @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_hessians flag must be an element of the
     * list of UpdateFlags that you passed to the constructor of this object.
     * See "The interplay of UpdateFlags, Mapping, and FiniteElement" in
     * the documentation of the FEValues class for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_laplacian_type<Number>> &
    get_laplacians(const std::string &global_vector_name,
                   const Extractor &  variable,
                   const Number       dummy = Number(0));

    /**
     * For the solution vector identified by @p global_vector_name, compute
     * the third_derivatives of the function at the quadrature points, and
     * return a vector with the correct type deduced by the Extractor you passed
     * as the @p variable argument.
     *
     * Before you can call this method, you need to call the
     * extract_local_dof_values() method at least once, passing the same
     * @p global_vector_name string, and the same type for the variable @p dummy.
     *
     * If you have not previously called the extract_local_dof_values() method,
     * this function will throw an exception.
     *
     * For this function to work properly, the underlying FEValues,
     * FEFaceValues, or FESubfaceValues object for which you called one of the
     * reinit() functions, must have computed the information you are
     * requesting. To do so, the update_3rd_derivatives flag must be an
     * element of the list of UpdateFlags that you passed to the constructor of
     * this object. See "The interplay of UpdateFlags, Mapping, and
     * FiniteElement" in the documentation of the FEValues for more information.
     */
    template <typename Extractor, typename Number = double>
    const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                        template solution_third_derivative_type<Number>> &
    get_third_derivatives(const std::string &global_vector_name,
                          const Extractor &  variable,
                          const Number       dummy = Number(0));

    /** @} */ // CurrentCellEvaluation

    /**
     * Return a reference to the used mapping.
     */
    const Mapping<dim, spacedim> &
    get_mapping() const;

  private:
    /**
     * Construct a unique name to store vectors of values, gradients,
     * divergences, etc., in the internal GeneralDataStorage object.
     */
    template <typename Extractor, typename Number = double>
    std::string
    get_unique_name(const std::string &global_vector_name,
                    const Extractor &  variable,
                    const std::string &object_type,
                    const unsigned int size,
                    const Number &     exemplar_number) const;

    /**
     * Construct a unique name to store local dof values.
     */
    template <typename Number = double>
    std::string
    get_unique_dofs_name(const std::string &global_vector_name,
                         const unsigned int size,
                         const Number &     exemplar_number) const;

    /**
     * The mapping used by the internal FEValues. Make sure it lives
     * longer than this class.
     */
    SmartPointer<const Mapping<dim, spacedim>> mapping;

    /**
     * The finite element used by the internal FEValues. Make sure it lives
     * longer than this class.
     */
    SmartPointer<const FiniteElement<dim, spacedim>> fe;

    /**
     * Quadrature formula used to integrate on the current cell, and on its
     * neighbor.
     */
    Quadrature<dim> cell_quadrature;

    /**
     * Quadrature formula used to integrate on faces, subfaces, and neighbor
     * faces and subfaces.
     */
    Quadrature<dim - 1> face_quadrature;

    /**
     * UpdateFlags to use when initializing the cell FEValues object.
     */
    UpdateFlags cell_update_flags;

    /**
     * UpdateFlags to use when initializing the neighbor cell FEValues objects.
     */
    UpdateFlags neighbor_cell_update_flags;

    /**
     * UpdateFlags to use when initializing FEFaceValues and FESubfaceValues
     * objects.
     */
    UpdateFlags face_update_flags;

    /**
     * UpdateFlags to use when initializing neighbor FEFaceValues and
     * FESubfaceValues objects.
     */
    UpdateFlags neighbor_face_update_flags;

    /**
     * Finite element values on the current cell.
     */
    std::unique_ptr<FEValues<dim, spacedim>> fe_values;

    /**
     * Finite element values on the current face.
     */
    std::unique_ptr<FEFaceValues<dim, spacedim>> fe_face_values;

    /**
     * Finite element values on the current subface.
     */
    std::unique_ptr<FESubfaceValues<dim, spacedim>> fe_subface_values;

    /**
     * Finite element values on the neighbor cell.
     */
    std::unique_ptr<FEValues<dim, spacedim>> neighbor_fe_values;

    /**
     * Finite element values on the neighbor face.
     */
    std::unique_ptr<FEFaceValues<dim, spacedim>> neighbor_fe_face_values;

    /**
     * Finite element values on the neighbor subface.
     */
    std::unique_ptr<FESubfaceValues<dim, spacedim>> neighbor_fe_subface_values;

    /**
     * Interface values on facets.
     */
    std::unique_ptr<FEInterfaceValues<dim, spacedim>> interface_fe_values;

    /**
     * Dof indices on the current cell.
     */
    std::vector<types::global_dof_index> local_dof_indices;

    /**
     * Dof indices on the neighbor cell.
     */
    std::vector<types::global_dof_index> neighbor_dof_indices;

    /**
     * User data storage.
     */
    GeneralDataStorage user_data_storage;

    /**
     * Internal data storage.
     */
    GeneralDataStorage internal_data_storage;

    /**
     * A pointer to the last used FEValues/FEFaceValues, or FESubfaceValues
     * object on this cell.
     */
    SmartPointer<const FEValuesBase<dim, spacedim>> current_fe_values;

    /**
     * A pointer to the last used FEValues/FEFaceValues, or FESubfaceValues
     * object on the neighbor cell.
     */
    SmartPointer<const FEValuesBase<dim, spacedim>> current_neighbor_fe_values;
  };

#ifndef DOXYGEN
  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  std::string
  ScratchData<dim, spacedim>::get_unique_name(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const std::string &object_type,
    const unsigned int size,
    const Number &     exemplar_number) const
  {
    return global_vector_name + "_" + variable.get_name() + "_" + object_type +
           "_" + Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <int dim, int spacedim>
  template <typename Number>
  std::string
  ScratchData<dim, spacedim>::get_unique_dofs_name(
    const std::string &global_vector_name,
    const unsigned int size,
    const Number &     exemplar_number) const
  {
    return global_vector_name + "_independent_local_dofs_" +
           Utilities::int_to_string(size) + "_" +
           Utilities::type_to_string(exemplar_number);
  }



  template <int dim, int spacedim>
  template <typename VectorType, typename Number>
  void
  ScratchData<dim, spacedim>::extract_local_dof_values(
    const std::string &global_vector_name,
    const VectorType & input_vector,
    const Number       dummy)
  {
    const unsigned int n_dofs = local_dof_indices.size();

    const std::string name =
      get_unique_dofs_name(global_vector_name, n_dofs, dummy);

    auto &independent_local_dofs =
      internal_data_storage
        .template get_or_add_object_with_name<std::vector<Number>>(name,
                                                                   n_dofs);

    AssertDimension(independent_local_dofs.size(), n_dofs);

    if (Differentiation::AD::is_tapeless_ad_number<Number>::value == true)
      for (unsigned int i = 0; i < n_dofs; ++i)
        Differentiation::AD::internal::Marking<Number>::independent_variable(
          input_vector(local_dof_indices[i]),
          i,
          n_dofs,
          independent_local_dofs[i]);
    else
      for (unsigned int i = 0; i < n_dofs; ++i)
        independent_local_dofs[i] = input_vector(local_dof_indices[i]);
  }



  template <int dim, int spacedim>
  template <typename Number>
  const std::vector<Number> &
  ScratchData<dim, spacedim>::get_local_dof_values(
    const std::string &global_vector_name,
    Number             dummy) const
  {
    const unsigned int n_dofs =
      get_current_fe_values().get_fe().n_dofs_per_cell();

    const std::string dofs_name =
      get_unique_dofs_name(global_vector_name, n_dofs, dummy);

    Assert(
      internal_data_storage.stores_object_with_name(dofs_name),
      ExcMessage(
        "You did not call yet extract_local_dof_values with the right types!"));

    return internal_data_storage
      .template get_object_with_name<std::vector<Number>>(dofs_name);
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_value_type<Number>> &
  ScratchData<dim, spacedim>::get_values(const std::string &global_vector_name,
                                         const Extractor &  variable,
                                         const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_values_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_value_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_values_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_gradient_type<Number>> &
  ScratchData<dim, spacedim>::get_gradients(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_gradients_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_gradient_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_gradients_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_hessian_type<Number>> &
  ScratchData<dim, spacedim>::get_hessians(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_hessians_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_hessian_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_hessians_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_laplacian_type<Number>> &
  ScratchData<dim, spacedim>::get_laplacians(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_laplacians_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_laplacian_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_laplacians_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_third_derivative_type<Number>> &
  ScratchData<dim, spacedim>::get_third_derivatives(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_third_derivatives_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_third_derivative_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_third_derivatives_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_symmetric_gradient_type<Number>> &
  ScratchData<dim, spacedim>::get_symmetric_gradients(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_symmetric_gradient_q", n_q_points, dummy);


    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_symmetric_gradient_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_symmetric_gradients_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }


  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_divergence_type<Number>> &
  ScratchData<dim, spacedim>::get_divergences(
    const std::string &global_vector_name,
    const Extractor &  variable,
    const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_divergence_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_divergence_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);


    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_divergences_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }



  template <int dim, int spacedim>
  template <typename Extractor, typename Number>
  const std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                      template solution_curl_type<Number>> &
  ScratchData<dim, spacedim>::get_curls(const std::string &global_vector_name,
                                        const Extractor &  variable,
                                        const Number       dummy)
  {
    const std::vector<Number> &independent_local_dofs =
      get_local_dof_values(global_vector_name, dummy);

    const FEValuesBase<dim, spacedim> &fev = get_current_fe_values();

    const unsigned int n_q_points = fev.n_quadrature_points;

    const std::string name = get_unique_name(
      global_vector_name, variable, "_curl_q", n_q_points, dummy);

    // Now build the return type
    using RetType =
      std::vector<typename FEValuesViews::View<dim, spacedim, Extractor>::
                    template solution_curl_type<Number>>;

    RetType &ret =
      internal_data_storage.template get_or_add_object_with_name<RetType>(
        name, n_q_points);

    AssertDimension(ret.size(), n_q_points);

    fev[variable].get_function_curls_from_local_dof_values(
      independent_local_dofs, ret);
    return ret;
  }

#endif

} // namespace MeshWorker

DEAL_II_NAMESPACE_CLOSE

#endif
