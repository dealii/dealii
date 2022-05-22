// ---------------------------------------------------------------------
//
// Copyright (C) 2021 - 2022 by the deal.II authors
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

#ifndef dealii_non_matching_fe_values
#define dealii_non_matching_fe_values

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/smartpointer.h>
#include <deal.II/base/std_cxx17/optional.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/fe_update_flags.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>
#include <deal.II/hp/mapping_collection.h>
#include <deal.II/hp/q_collection.h>

#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/non_matching/quadrature_generator.h>

#include <deque>

DEAL_II_NAMESPACE_OPEN

namespace NonMatching
{
  /**
   * Struct storing UpdateFlags for the 3 regions of a cell, $K$, that is
   * defined by the sign of a level set function, $\psi$:
   * @f[
   * N = \{x \in K : \psi(x) < 0 \}, \\
   * P = \{x \in K : \psi(x) > 0 \}, \\
   * S = \{x \in K : \psi(x) = 0 \}.
   * @f]
   * As in the QuadratureGenerator class, we refer to $N$, $P$ and $S$ as the
   * inside, outside, and surface region. RegionUpdateFlags is used to describe
   * how the FEValues objects, which are created by NonMatching::FEValues,
   * should be updated.
   */
  struct RegionUpdateFlags
  {
    /**
     * Constructor, sets the UpdateFlags for each region to update_default.
     */
    RegionUpdateFlags();

    /**
     * Flags for the region $\{x \in K : \psi(x) < 0 \}$
     */
    UpdateFlags inside;

    /**
     * Flags for the region $\{x \in K : \psi(x) > 0 \}$
     */
    UpdateFlags outside;

    /**
     * Flags for the region $\{x \in K : \psi(x) = 0 \}$
     */
    UpdateFlags surface;
  };


  /**
   * This class is intended to facilitate assembling in immersed (in the sense
   * of cut) finite element methods when the domain is described by a level set
   * function, $\psi : \mathbb{R}^{dim} \to \mathbb{R}$. In this type of
   * method, we typically need to integrate over 3 different regions of each
   * cell, $K$:
   * @f[
   * N = \{x \in K : \psi(x) < 0 \}, \\
   * P = \{x \in K : \psi(x) > 0 \}, \\
   * S = \{x \in K : \psi(x) = 0 \}.
   * @f]
   * Thus we need quadrature rules for these 3 regions:
   * @image html immersed_quadratures.svg
   * As in the QuadratureGenerator class, we refer to $N$, $P$, and $S$ as the
   * inside, outside, and surface regions. The constructor of this class takes a
   * discrete level set function described by a DoFHandler and a Vector. When
   * the reinit() function is called, the QuadratureGenerator will be called in
   * the background to create these immersed quadrature rules. This class then
   * creates dealii::FEValues objects for the inside/outside regions and an
   * FEImmersedSurfaceValues object for the surface region. These objects can
   * then be accessed through one of the functions: get_inside_fe_values(),
   * get_outside_fe_values(), or
   * get_surface_fe_values().
   * Since a cut between a cell and the domain can be arbitrarily small, the
   * underlying algorithm may generate a quadrature rule with 0 points. This
   * can, for example, happen if the relative size of the cut is similar to the
   * floating-point accuracy. Since the FEValues-like objects are not allowed to
   * contain 0 points, the object that get_inside/outside/surface_fe_values()
   * returns is wrapped in a std_cxx17::optional. This requires us to check if
   * the returned FEValues-like object contains a value before we use it:
   * @code
   * NonMatching::FEValues<dim> fe_values(...);
   *
   * for (const auto &cell : dof_handler.active_cell_iterators())
   *  {
   *    fe_values.reinit(cell);
   *
   *    const std_cxx17::optional<FEValues<dim>> &fe_values_inside =
   *      fe_values.get_inside_fe_values();
   *
   *    if (fe_values_inside)
   *      {
   *        // Assemble locally
   *        for (const unsigned int q_index :
   *             fe_values_inside->quadrature_point_indices())
   *          {
   *            // ...
   *          }
   *      }
   *  }
   * @endcode
   *
   * Of course, it is somewhat expensive to generate the immersed quadrature
   * rules and create FEValues objects with the generated quadratures. To reduce
   * the amount of work, the reinit() function of this class uses the
   * MeshClassifier passed to the constructor to check how the incoming cell
   * relates to the level set function. It only generates the immersed
   * quadrature rules if the cell is intersected. If the cell is completely
   * inside or outside, it returns a cached FEValues object created with a
   * quadrature over the reference cell: $[0, 1]^{dim}$.
   */
  template <int dim>
  class FEValues
  {
  public:
    using AdditionalData = typename QuadratureGenerator<dim>::AdditionalData;

    /**
     * Constructor.
     *
     * @param fe_collection Collection of FiniteElements to be used.
     * @param quadrature 1-dimensional quadrature rule used to generate the
     * immersed quadrature rules. See the QuadratureGenerator class. On the non
     * intersected elements a tensor product of this quadrature will be used.
     * @param mesh_classifier Object used to determine when the immersed
     * quadrature rules need to be generated.
     * @param region_update_flags Struct storing UpdateFlags for the
     * inside/outside/surface region of the cell.
     * @param dof_handler The DoFHandler associated with the discrete level set
     * function.
     * @param level_set The degrees of freedom of the discrete level set function.
     * @param additional_data Additional data passed to the QuadratureGenerator
     * class.
     *
     * @note  Pointers to @p fe_collection, @p mesh_classifier, @p dof_handler,
     * and @p level_set are stored internally, so these need to have a longer life
     * span than the instance of this class.
     */
    template <class VectorType>
    FEValues(const hp::FECollection<dim> &fe_collection,
             const Quadrature<1> &        quadrature,
             const RegionUpdateFlags      region_update_flags,
             const MeshClassifier<dim> &  mesh_classifier,
             const DoFHandler<dim> &      dof_handler,
             const VectorType &           level_set,
             const AdditionalData &       additional_data = AdditionalData());

    /**
     * Constructor.
     *
     * @param mapping_collection Collection of Mappings to be used.
     * @param fe_collection Collection of FiniteElements to be used.
     * @param q_collection Collection of Quadrature rules over $[0, 1]^{dim}$
     * that should be used when a cell is not intersected and we do not need to
     * generate immersed quadrature rules.
     * @param q_collection_1D Collection of 1-dimensional quadrature rules used
     * to generate the immersed quadrature rules. See the QuadratureGenerator
     * class.
     * @param mesh_classifier Object used to determine when the immersed
     * quadrature rules need to be generated.
     * @param region_update_flags Struct storing UpdateFlags for the
     * inside/outside/surface region of the cell.
     * @param dof_handler The DoFHandler associated with the discrete level set
     * function.
     * @param level_set The degrees of freedom of the discrete level set function.
     * @param additional_data Additional data passed to the QuadratureGenerator
     * class.
     *
     * @note Pointers to @p mapping_collection, @p fe_collection,
     * @p mesh_classifier, @p dof_handler, and @p level_set are stored
     * internally, so these need to have a longer life span than the instance of
     * this class.
     */
    template <class VectorType>
    FEValues(const hp::MappingCollection<dim> &mapping_collection,
             const hp::FECollection<dim> &     fe_collection,
             const hp::QCollection<dim> &      q_collection,
             const hp::QCollection<1> &        q_collection_1D,
             const RegionUpdateFlags           region_update_flags,
             const MeshClassifier<dim> &       mesh_classifier,
             const DoFHandler<dim> &           dof_handler,
             const VectorType &                level_set,
             const AdditionalData &additional_data = AdditionalData());

    /**
     * Reinitialize the various FEValues-like objects for the 3 different
     * regions of the cell. After calling this function an FEValues-like object
     * can be retrieved for each region using the functions
     * get_inside_fe_values(),
     * get_outside_fe_values(), or
     * get_surface_fe_values().
     */
    template <bool level_dof_access>
    void
    reinit(
      const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell);

    /**
     * Return an dealii::FEValues object reinitialized with a quadrature for the
     * inside region of the cell: $\{x \in K : \psi(x) < 0 \}$.
     *
     * @note If the quadrature rule over the region is empty, e.g. because the
     * cell is completely located in the outside domain, the returned
     * optional will not contain a value.
     */
    const std_cxx17::optional<dealii::FEValues<dim>> &
    get_inside_fe_values() const;

    /**
     * Return an dealii::FEValues object reinitialized with a quadrature for the
     * outside region of the cell: $\{x \in K : \psi(x) > 0 \}$.
     *
     * @note If the quadrature rule over the region is empty, e.g. because the
     * cell is completely located in the inside domain, the returned
     * optional will not contain a value.
     */
    const std_cxx17::optional<dealii::FEValues<dim>> &
    get_outside_fe_values() const;

    /**
     * Return an dealii::FEValues object reinitialized with a quadrature for the
     * surface region of the cell: $\{x \in K : \psi(x) = 0 \}$.
     *
     * @note If the quadrature rule over the region is empty, e.g. because the
     * cell is not intersected, the returned optional will not contain a value.
     */
    const std_cxx17::optional<FEImmersedSurfaceValues<dim>> &
    get_surface_fe_values() const;

  private:
    /**
     * Do work common to the constructors. The incoming QCollection should be
     * quadratures integrating over $[0, 1]^{dim}$. These will be used on the
     * non-intersected cells.
     */
    void
    initialize(const hp::QCollection<dim> &q_collection);

    /**
     * A pointer to the collection of mappings to be used.
     */
    const SmartPointer<const hp::MappingCollection<dim>> mapping_collection;

    /**
     * A pointer to the collection of finite elements to be used.
     */
    const SmartPointer<const hp::FECollection<dim>> fe_collection;

    /**
     * Collection of 1-dimensional quadrature rules that are used by
     * QuadratureGenerator as base for generating the immersed quadrature rules.
     */
    const hp::QCollection<1> q_collection_1D;

    /**
     * Location of the last cell that reinit was called with.
     */
    LocationToLevelSet current_cell_location;

    /**
     * Active fe index of the last cell that reinit was called with.
     */
    unsigned int active_fe_index;

    /**
     * The update flags passed to the constructor.
     */
    const RegionUpdateFlags region_update_flags;

    /**
     * Pointer to the MeshClassifier passed to the constructor.
     */
    const SmartPointer<const MeshClassifier<dim>> mesh_classifier;

    /**
     * For each element in the FECollection passed to the constructor,
     * this object contains an dealii::FEValues object created with a quadrature
     * rule over the full reference cell: $[0, 1]^{dim}$ and UpdateFlags for the
     * inside region. Thus, these optionals should always contain a value.
     *
     * When LocationToLevelSet of the cell is INSIDE (and we do not need
     * to generate an immersed quadrature), we return the dealii::FEValues
     * object in this container corresponding to the cell's active_fe_index.
     *
     * This container is a std::deque, which is compatible with the `FEValues`
     * class that does not have a copy-constructor.
     */
    std::deque<std_cxx17::optional<dealii::FEValues<dim>>>
      fe_values_inside_full_quadrature;

    /**
     * For each element in the FECollection passed to the constructor,
     * this object contains an dealii::FEValues object created with a quadrature
     * rule over the full reference cell: $[0, 1]^{dim}$ and UpdateFlags for the
     * outside region. Thus, these optionals should always contain a value.
     *
     * When LocationToLevelSet of the cell is OUTSIDE (and we do not need
     * to generate an immersed quadrature), we return the dealii::FEValues
     * object in this container corresponding to the cell's active_fe_index.
     *
     * This container is a std::deque, which is compatible with the `FEValues`
     * class that does not have a copy-constructor.
     */
    std::deque<std_cxx17::optional<dealii::FEValues<dim>>>
      fe_values_outside_full_quadrature;

    /**
     * FEValues object created with a quadrature rule integrating over the
     * inside region, $\{x \in B : \psi(x) < 0 \}$, that was generated in the
     * last call to reinit(..). If the cell in the last call was not intersected
     * or if 0 quadrature points were generated, this optional will not contain
     * a value.
     */
    std_cxx17::optional<dealii::FEValues<dim>> fe_values_inside;

    /**
     * FEValues object created with a quadrature rule integrating over the
     * outside region, $\{x \in B : \psi(x) > 0 \}$, that was generated in the
     * last call to reinit(..). If the cell in the last call was not intersected
     * or if 0 quadrature points were generated, this optional will not contain
     * a value.
     */
    std_cxx17::optional<dealii::FEValues<dim>> fe_values_outside;

    /**
     * FEImmersedSurfaceValues object created with a quadrature rule integrating
     * over the surface region, $\{x \in B : \psi(x) = 0 \}$, that was generated
     * in the last call to reinit(..). If the cell in the last call was not
     * intersected or if 0 quadrature points were generated, this optional will
     * not contain a value.
     */
    std_cxx17::optional<NonMatching::FEImmersedSurfaceValues<dim>>
      fe_values_surface;

    /**
     * Object that generates the immersed quadrature rules.
     */
    DiscreteQuadratureGenerator<dim> quadrature_generator;
  };



  /**
   * This class is intended to facilitate assembling interface terms on faces in
   * immersed (in the sense of cut) finite element methods. These types of terms
   * occur mainly in cut discontinuous Galerkin methods. This class works
   * analogously to NonMatching::FEValues. The domain is assumed to be described
   * by a level set function,
   * $\psi : \mathbb{R}^{dim} \to \mathbb{R}$,
   * and this class assumes that we want to integrate over two different regions
   * of each face, $F$:
   * @f[
   * N = \{x \in F : \psi(x) < 0 \}, \\
   * P = \{x \in F : \psi(x) > 0 \},
   * @f]
   * which we as before refer to as the "inside" and "outside" regions of the
   * face.
   * @image html non_matching_fe_interface_values.svg
   * When calling one of the reinit functions of this class, this class
   * will create quadrature rules over these regions in the background and set
   * up dealii::FEInterfaceValues objects with these quadratures. These objects
   * can then be accessed through one of the functions: get_inside_fe_values()
   * or get_outside_fe_values(). For the same reason as described in
   * NonMatching::FEValues, the returned the FEInterfaceValues objects that
   * get_inside/outside_fe_values() returns is wrapped in a
   * std_cxx17::optional. Assembling using this class would then typically be
   * done in the following way:
   *
   * @code
   *  NonMatching::FEInterfaceValues<dim> fe_interface_values(...)
   *
   *  for (const auto &cell : dof_handler.active_cell_iterators())
   *    {
   *      for (unsigned int face_index : cell->face_indices())
   *        {
   *          if (cell->at_boundary(face_index))
   *            fe_interface_values.reinit(cell, face_index);
   *          else
   *            fe_interface_values.reinit(cell,
   *                                       face_index,
   *                                       subface_index,
   *                                       cell->neighbor(face_index),
   *                                       cell->neighbor_of_neighbor(face_index),
   *                                       neighbor_subface_index);
   *
   *          const std_cxx17::optional<dealii::FEInterfaceValues<dim>>
   *            &inside_fe_values = fe_interface_values.get_inside_fe_values();
   *
   *          // Check if the returned optional contains a value
   *          if (inside_fe_values)
   *            {
   *              // Assemble locally
   *            }
   *        }
   *    }
   * @endcode
   *
   * To reduce the amount of work, the reinit() function of this class uses the
   * MeshClassifier passed to the constructor to check how the incoming cell
   * relates to the level set function. The immersed quadrature rules are only
   * generated if the cell is intersected. If the cell is completely inside or
   * outside, it returns a cached FEInterfaceValues object created with a
   * quadrature over the reference cell: $[0, 1]^{dim-1}$.
   */
  template <int dim>
  class FEInterfaceValues
  {
  public:
    using AdditionalData =
      typename FaceQuadratureGenerator<dim>::AdditionalData;

    /**
     * Constructor.
     *
     * @param fe_collection Collection of FiniteElements to be used.
     * @param quadrature 1-dimensional quadrature rule used to generate the
     * immersed quadrature rules. See the QuadratureGenerator class. On the non
     * intersected faces a tensor product of this quadrature will be used.
     * @param mesh_classifier Object used to determine when the immersed
     * quadrature rules need to be generated.
     * @param region_update_flags Struct storing UpdateFlags for the
     * inside/outside region of the cell.
     * @param dof_handler The DoFHandler associated with the discrete level set
     * function.
     * @param level_set The degrees of freedom of the discrete level set function.
     * @param additional_data Additional data passed to the QuadratureGenerator
     * class.
     *
     * @note  Pointers to @p fe_collection, @p mesh_classifier, @p dof_handler,
     * and @p level_set are stored internally, so these need to have a longer life
     * span than the instance of this class.
     */
    template <class VectorType>
    FEInterfaceValues(const hp::FECollection<dim> &fe_collection,
                      const Quadrature<1> &        quadrature,
                      const RegionUpdateFlags      region_update_flags,
                      const MeshClassifier<dim> &  mesh_classifier,
                      const DoFHandler<dim> &      dof_handler,
                      const VectorType &           level_set,
                      const AdditionalData &additional_data = AdditionalData());

    /**
     * Constructor.
     *
     * @param mapping_collection Collection of Mappings to be used.
     * @param fe_collection Collection of FiniteElements to be used.
     * @param q_collection Collection of Quadrature rules over $[0, 1]^{dim-1}$
     * that should be used when a face is not intersected and we do not need to
     * generate immersed quadrature rules.
     * @param q_collection_1D Collection of 1-dimensional quadrature rules used
     * to generate the immersed quadrature rules. See the QuadratureGenerator
     * class.
     * @param mesh_classifier Object used to determine when the immersed
     * quadrature rules need to be generated.
     * @param region_update_flags Struct storing UpdateFlags for the
     * inside/outside region of the cell.
     * @param dof_handler The DoFHandler associated with the discrete level set
     * function.
     * @param level_set The degrees of freedom of the discrete level set function.
     * @param additional_data Additional data passed to the QuadratureGenerator
     * class.
     *
     * @note Pointers to @p mapping_collection, @p fe_collection,
     * @p mesh_classifier, @p dof_handler, and @p level_set are stored
     * internally, so these need to have a longer life span than the instance of
     * this class.
     */
    template <class VectorType>
    FEInterfaceValues(const hp::MappingCollection<dim> &mapping_collection,
                      const hp::FECollection<dim> &     fe_collection,
                      const hp::QCollection<dim - 1> &  q_collection,
                      const hp::QCollection<1> &        q_collection_1D,
                      const RegionUpdateFlags           region_update_flags,
                      const MeshClassifier<dim> &       mesh_classifier,
                      const DoFHandler<dim> &           dof_handler,
                      const VectorType &                level_set,
                      const AdditionalData &additional_data = AdditionalData());

    /**
     * Reinitialize on the shared face between two neighboring cells.
     * After calling this function, an FEInterfaceValues object can be retrieved
     * for each region using the functions get_inside_fe_values() or
     * get_outside_fe_values().
     *
     * @note Currently @p cell and @p cell_neighbor need to have the same
     * finite element associated with them.
     *
     * @note Specifying sub_face_no/sub_face_no_neighbor is currently not
     * implemented, the passed values of these must equal
     * numbers::invalid_unsigned_int.
     */
    template <class CellIteratorType, class CellNeighborIteratorType>
    void
    reinit(const CellIteratorType &        cell,
           const unsigned int              face_no,
           const unsigned int              sub_face_no,
           const CellNeighborIteratorType &cell_neighbor,
           const unsigned int              face_no_neighbor,
           const unsigned int              sub_face_no_neighbor);


    /**
     * Reinitialize on a face on the boundary of the domain,
     * when there is no neighbor and a single face @p face_no of
     * the @p cell. After calling this function, an
     * FEInterfaceValues object can be retrieved for each region using the
     * functions get_inside_fe_values() or get_outside_fe_values().
     */
    template <class CellIteratorType>
    void
    reinit(const CellIteratorType &cell, const unsigned int face_no);

    /**
     * Return an dealii::FEInterfaceValues object reinitialized with a
     * quadrature for the inside region of the cell:
     * $\{x \in K : \psi(x) < 0 \}$.
     *
     * @note If the quadrature rule over the region is empty, e.g. because the
     * cell is completely located in the outside domain, the returned
     * optional will not contain a value.
     */
    const std_cxx17::optional<dealii::FEInterfaceValues<dim>> &
    get_inside_fe_values() const;

    /**
     * Return an dealii::FEInterfaceValues object reinitialized with a
     * quadrature for the outside region of the cell:
     * $\{x \in K : \psi(x) > 0 \}$.
     *
     * @note If the quadrature rule over the region is empty, e.g. because the
     * cell is completely located in the inside domain, the returned
     * optional will not contain a value.
     */
    const std_cxx17::optional<dealii::FEInterfaceValues<dim>> &
    get_outside_fe_values() const;

  private:
    /**
     * Do work common to the constructors. The incoming QCollection should be
     * quadratures integrating over $[0, 1]^{dim-1}$. These will be used on the
     * non-intersected cells.
     */
    void
    initialize(const hp::QCollection<dim - 1> &q_collection);

    /**
     * Do work that is common two the two reinit functions of the class.
     * @p call_reinit is a std::function that describes how we should call
     * reinit on a single dealii::FEInterfaceValues object, which is what
     * differs between the two reinit functions.
     */
    template <bool level_dof_access>
    void
    do_reinit(
      const TriaIterator<DoFCellAccessor<dim, dim, level_dof_access>> &cell,
      const unsigned int                                               face_no,
      const std::function<void(dealii::FEInterfaceValues<dim> &)> &call_reinit);

    /**
     * A pointer to the collection of mappings to be used.
     */
    const SmartPointer<const hp::MappingCollection<dim>> mapping_collection;

    /**
     * A pointer to the collection of finite elements to be used.
     */
    const SmartPointer<const hp::FECollection<dim>> fe_collection;

    /**
     * Collection of 1-dimensional quadrature rules that are used by
     * QuadratureGenerator as base for generating the immersed quadrature rules.
     */
    const hp::QCollection<1> q_collection_1D;

    /**
     * Location of the last face that reinit was called with.
     */
    LocationToLevelSet current_face_location;

    /**
     * Active fe index of the last cell that reinit was called with.
     */
    unsigned int active_fe_index;

    /**
     * The update flags passed to the constructor.
     */
    const RegionUpdateFlags region_update_flags;

    /**
     * Pointer to the MeshClassifier passed to the constructor.
     */
    const SmartPointer<const MeshClassifier<dim>> mesh_classifier;

    /**
     * For each element in the FECollection passed to the constructor,
     * this object contains an dealii::FEInterfaceValues object created with a
     * quadrature rule over the full reference cell: $[0, 1]^{dim-1}$ and
     * UpdateFlags for the inside region. Thus, these optionals should always
     * contain a value.
     *
     * When LocationToLevelSet of the cell is INSIDE (and we do not need
     * to generate an immersed quadrature), we return the
     * dealii::FEInterfaceValues object in this container corresponding to the
     * cell's active_fe_index.
     *
     * This container is a std::deque, which is compatible with the
     * `FEInterfaceValues` class that does not have a copy-constructor.
     */
    std::deque<std_cxx17::optional<dealii::FEInterfaceValues<dim>>>
      fe_values_inside_full_quadrature;

    /**
     * For each element in the FECollection passed to the constructor,
     * this object contains an dealii::FEInterfaceValues object created with a
     * quadrature rule over the full reference cell: $[0, 1]^{dim-1}$ and
     * UpdateFlags for the outside region. Thus, these optionals should always
     * contain a value.
     *
     * When LocationToLevelSet of the cell is OUTSIDE (and we do not need
     * to generate an immersed quadrature), we return the dealii::FEValues
     * object in this container corresponding to the cell's active_fe_index.
     *
     * This container is a std::deque, which is compatible with the
     * `FEInterfaceValues` class that does not have a copy-constructor.
     */
    std::deque<std_cxx17::optional<dealii::FEInterfaceValues<dim>>>
      fe_values_outside_full_quadrature;

    /**
     * FEInterfaceValues object created with a quadrature rule integrating over
     * the inside region, $\{x \in B : \psi(x) < 0 \}$, that was generated in
     * the last call to reinit(..). If the cell in the last call was not
     * intersected or if 0 quadrature points were generated, this optional will
     * not contain a value.
     */
    std_cxx17::optional<dealii::FEInterfaceValues<dim>> fe_values_inside;

    /**
     * FEInterfaceValues object created with a quadrature rule integrating over
     * the outside region, $\{x \in B : \psi(x) > 0 \}$, that was generated in
     * the last call to reinit(..). If the cell in the last call was not
     * intersected or if 0 quadrature points were generated, this optional will
     * not contain a value.
     */
    std_cxx17::optional<dealii::FEInterfaceValues<dim>> fe_values_outside;

    /**
     * Object that generates the immersed quadrature rules.
     */
    DiscreteFaceQuadratureGenerator<dim> face_quadrature_generator;
  };


#ifndef DOXYGEN

  /*---------------------- Inline functions ---------------------*/

  template <int dim>
  template <class CellIteratorType>
  inline void
  FEInterfaceValues<dim>::reinit(const CellIteratorType &cell,
                                 const unsigned int      face_no)
  {
    // Lambda describing how we should call reinit on a single
    // dealii::FEInterfaceValues object.
    const auto reinit_operation =
      [&cell, face_no](dealii::FEInterfaceValues<dim> &fe_interface_values) {
        fe_interface_values.reinit(cell, face_no);
      };

    do_reinit(cell, face_no, reinit_operation);
  }



  template <int dim>
  template <class CellIteratorType, class CellNeighborIteratorType>
  inline void
  FEInterfaceValues<dim>::reinit(const CellIteratorType &        cell,
                                 const unsigned int              face_no,
                                 const unsigned int              sub_face_no,
                                 const CellNeighborIteratorType &cell_neighbor,
                                 const unsigned int face_no_neighbor,
                                 const unsigned int sub_face_no_neighbor)
  {
    Assert(sub_face_no == numbers::invalid_unsigned_int, ExcNotImplemented());
    Assert(sub_face_no_neighbor == numbers::invalid_unsigned_int,
           ExcNotImplemented());

    // Lambda describing how we should call reinit on a single
    // dealii::FEInterfaceValues object.
    const auto reinit_operation =
      [&cell,
       face_no,
       sub_face_no,
       &cell_neighbor,
       face_no_neighbor,
       sub_face_no_neighbor](
        dealii::FEInterfaceValues<dim> &fe_interface_values) {
        fe_interface_values.reinit(cell,
                                   face_no,
                                   sub_face_no,
                                   cell_neighbor,
                                   face_no_neighbor,
                                   sub_face_no_neighbor);
      };

    do_reinit(cell, face_no, reinit_operation);
  }

#endif

} // namespace NonMatching
DEAL_II_NAMESPACE_CLOSE

#endif
