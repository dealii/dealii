// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2016 - 2024 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_quadrature_point_data_h
#define dealii_quadrature_point_data_h

#include <deal.II/base/config.h>

#include <deal.II/base/enable_observer_pointer.h>
#include <deal.II/base/quadrature.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/fe.h>
#include <deal.II/fe/fe_tools.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>

#include <deal.II/lac/vector.h>

#include <map>
#include <optional>
#include <type_traits>
#include <vector>

DEAL_II_NAMESPACE_OPEN

/**
 * @addtogroup Quadrature
 * @{
 */

/**
 * A class for storing at each cell represented by iterators of type @p CellIteratorType
 * a vector of data @p DataType .
 * The underlying structure and the initialize() method of this class are
 * designed in such a way that one could use different child classes derived
 * from the base DataType class to store data on a given cell. This implies the
 * usage of pointers, in our case -- std::shared_ptr.
 *
 * The type @p DataType is arbitrary, but when using a class derived from
 * TransferableQuadraturePointData one can use the facilities of
 * parallel::distributed::ContinuousQuadratureDataTransfer.
 *
 * @note The data type stored on each cell can be different.
 * However, within the cell this class stores a vector of objects of a single
 * data type. For this reason, this class may not be sufficiently flexible when,
 * for example, adopting a level-set approach to describe material behavior.
 */
template <typename CellIteratorType, typename DataType>
class CellDataStorage : public EnableObserverPointer
{
public:
  /**
   * Default constructor.
   */
  CellDataStorage() = default;

  /**
   * Default destructor.
   */
  ~CellDataStorage() override = default;

  /**
   * Initialize data on the @p cell to store @p number_of_data_points_per_cell of objects of type @p T .
   * The typename @p T is possibly another class which is derived from the
   * base @p DataType class. In order to initialize the vector of objects
   * we have to assume that the class @p T has a default constructor.
   * This function has to be called on every cell where data is to be stored.
   *
   * After the data is initialized, it can be modified using the get_data()
   * function.
   *
   * @note Subsequent calls of this function with the same @p cell will not
   * alter the objects associated with it. In order to remove the stored data,
   * use the erase() function.
   *
   * @note It is possible to use different types @p T for different cells which
   * may reflect, for example, different constitutive models of continuum
   * mechanics in different parts of the domain.
   *
   * @note The first time this method is called, it stores a ObserverPointer to the
   * Triangulation object that owns the cell. The future invocations of this
   * method expects the cell to be from the same stored triangulation.
   *
   * @pre The type @p T needs to either equal @p DataType, or be a class derived
   * from @p DataType. @p T needs to be default constructible.
   */
  template <typename T = DataType>
  void
  initialize(const CellIteratorType &cell,
             const unsigned int      number_of_data_points_per_cell);

  /**
   * Same as above but for a range of iterators starting at @p cell_start
   * until, but not including, @p cell_end for all locally owned cells, i.e.
   * for which `cell->is_locally_owned()==true` .
   */
  template <typename T = DataType>
  void
  initialize(
    const CellIteratorType                                          &cell_start,
    const typename std_cxx20::type_identity<CellIteratorType>::type &cell_end,
    const unsigned int number_of_data_points_per_cell);

  /**
   * Removes data stored at the @p cell. Returns true if the data was removed.
   * If no data is attached to the @p cell, this function will not do anything
   * and returns false.
   *
   * @note This function will also check that there are no
   * outstanding references to the data stored on this cell. That is to say,
   * that the only references to the stored data are that made by this class.
   */
  bool
  erase(const CellIteratorType &cell);

  /**
   * Clear all the data stored in this object.
   */
  void
  clear();

  /**
   * Get a vector of the data located at @p cell .
   * A possible additional typename @p T is the class to which the base class
   * DataType could be cast. Since @p DataType is stored as shared pointers,
   * there is minimal overhead in returning a vector by value instead of by
   * reference.
   * This allows flexibility if class @p T is not the same as @p DataType on a
   * cell-by-cell basis.
   *
   * @pre The type @p T needs to match the class provided to initialize() .
   *
   * @pre @p cell must be from the same Triangulation that is used to
   * initialize() the cell data.
   */
  template <typename T = DataType>
  std::vector<std::shared_ptr<T>>
  get_data(const CellIteratorType &cell);

  /**
   * Get a vector of constant pointers to data located at @p cell .
   * A possible additional typename @p T is the class to which the base class
   * DataType could be cast. Since @p DataType is stored as shared pointers,
   * there is minimal overhead in returning a vector by value instead of by
   * reference.
   * This allows flexibility if class @p T is not the same as @p DataType on a
   * cell-by-cell basis.
   *
   * @pre The type @p T needs to match the class provided to initialize() .
   *
   * @pre @p cell must be from the same Triangulation that is used to
   * initialize() the cell data.
   */
  template <typename T = DataType>
  std::vector<std::shared_ptr<const T>>
  get_data(const CellIteratorType &cell) const;

  /**
   * Returns a std::optional indicating whether @p cell contains an
   * associated data or not. If data is available, dereferencing the
   * std::optional reveals a vector of pointers to the underlying data
   * at the quadrature points.
   * A possible additional typename @p T is the class to which the base class
   * DataType could be cast. Since @p DataType is stored as shared pointers,
   * there is minimal overhead in returning a vector by value instead of by
   * reference.
   * This allows flexibility if class @p T is not the same as @p DataType on a
   * cell-by-cell basis.
   *
   * @pre The type @p T needs to match the class provided to initialize().
   * @pre @p cell must be from the same Triangulation that is used to
   * initialize() the cell data.
   */
  template <typename T = DataType>
  std::optional<std::vector<std::shared_ptr<T>>>
  try_get_data(const CellIteratorType &cell);

  /**
   * Returns a std::optional indicating whether @p cell contains an
   * associated data or not. If data is available, dereferencing the
   * std::optional reveals a vector of constant pointers to the
   * underlying data at the quadrature points.
   * A possible additional typename @p T is the class to which the base class
   * DataType could be cast. Since @p DataType is stored as shared pointers,
   * there is minimal overhead in returning a vector by value instead of by
   * reference.
   * This allows flexibility if class @p T is not the same as @p DataType on a
   * cell-by-cell basis.
   *
   * @pre The type @p T needs to match the class provided to initialize().
   * @pre @p cell must be from the same Triangulation that is used to
   * initialize() the cell data.
   */
  template <typename T = DataType>
  std::optional<std::vector<std::shared_ptr<const T>>>
  try_get_data(const CellIteratorType &cell) const;

private:
  /**
   * Number of dimensions
   */
  static constexpr unsigned int dimension =
    CellIteratorType::AccessorType::dimension;

  /**
   * Number of space dimensions
   */
  static constexpr unsigned int space_dimension =
    CellIteratorType::AccessorType::space_dimension;

  /**
   * To ensure that all the cells in the CellDataStorage come from the same
   * Triangulation, we need to store a reference to that Triangulation within
   * the class.
   */
  ObserverPointer<const Triangulation<dimension, space_dimension>,
                  CellDataStorage<CellIteratorType, DataType>>
    tria;

  /**
   * A map to store a vector of data on each cell.
   * We need to use CellId as the key because it remains unique during
   * adaptive refinement.
   */
  std::map<CellId, std::vector<std::shared_ptr<DataType>>> map;

  /**
   * @addtogroup Exceptions
   */
  DeclExceptionMsg(
    ExcCellDataTypeMismatch,
    "Cell data is being retrieved with a type which is different than the type used to initialize it");

  /**
   * @addtogroup Exceptions
   */
  DeclExceptionMsg(
    ExcTriangulationMismatch,
    "The provided cell iterator does not belong to the triangulation that corresponds to the CellDataStorage object.");
};


/**
 * An abstract class which specifies requirements for data on
 * a single quadrature point to be transferable during refinement or
 * repartitioning.
 *
 * This class provides a framework by which derived classes representing data at
 * quadrature points can declare how many scalar values they store, and then
 * implement functions that pack and unpack these scalars into arrays.
 * These arrays are used to transfer data from quadrature points of one cell to
 * quadrature points of another cell as well as between processors
 * upon mesh refinement and repartitioning.
 * The transfer of quadrature point data between parent and child cells requires
 * some kind of projection and/or interpolation.
 * One possible implementation is via the L2 projection and prolongation
 * matrices as implemented in
 * parallel::distributed::ContinuousQuadratureDataTransfer class.
 *
 * To store and access instances of classes derived from this class, see the
 * CellDataStorage class.
 */
class TransferableQuadraturePointData
{
public:
  /**
   * Default constructor.
   */
  TransferableQuadraturePointData() = default;

  /**
   * Default virtual destructor.
   */
  virtual ~TransferableQuadraturePointData() = default;

  /**
   * Return the total number of values which will be
   * packed/unpacked from the user's DataType class. Consequently it is also
   * the size of the vectors in pack_values() and unpack_values() .
   */
  virtual unsigned int
  number_of_values() const = 0;

  /**
   * A virtual function that have to be implemented in derived classes to
   * pack all data stored in the derived class into a vector @p values .
   * This vector will contain all scalar and/or Tensorial data local to each
   * quadrature point.
   *
   * @note  The function will be called with @p values of size number_of_values().
   * The implementation may still have an assert to check that it is indeed the
   * case.
   */
  virtual void
  pack_values(std::vector<double> &values) const = 0;

  /**
   * The opposite of the above, namely
   * unpack a vector @p values into the data stored in this class.
   *
   * @note  The function will be called with @p values of size number_of_values().
   * The implementation may still have an assert to check that it is indeed the
   * case.
   */
  virtual void
  unpack_values(const std::vector<double> &values) = 0;
};


#ifdef DEAL_II_WITH_P4EST
namespace parallel
{
  namespace distributed
  {
    /**
     * A class for the transfer of continuous data stored at quadrature points
     * when performing h-adaptive refinement of
     * parallel::distributed::Triangulation .
     *
     * <h3>Implementation details</h3>
     *
     * This class implements the transfer of the quadrature point data between
     * cells in case of adaptive refinement using L2 projection. That also
     * includes automatic shipping of information between different processors.
     *
     * To that end, the constructor of the class is provided with three main
     * objects:
     * scalar FiniteElement @p projection_fe, @p mass_quadrature and @p data_quadrature
     * Quadrature rules.
     * First, the data located at @p data_quadrature of each cell is L2-projected
     * to the continuous space defined by a single FiniteElement @p projection_fe .
     * This is achieved using
     * FETools::compute_projection_from_quadrature_points_matrix(). In  doing so
     * the @ref GlossMassMatrix "mass matrix" of this element is required, which will be calculated
     * with the @p mass_quadrature rule . Should the cell now belong to another processor,
     * the data is then sent to this processor. The class makes use of a feature
     * of p4est (and parallel::distributed::Triangulation) that allows one to
     * attach information to cells during mesh refinement and rebalancing. On
     * receiving information on the target cell, the data is projected back to
     * the quadrature points using the matrix calculated by
     * FETools::compute_interpolation_to_quadrature_points_matrix() .
     * In the case that local refinement is performed, this class first
     * project local DoF values of the parent element to each child.
     *
     *
     * This class is templated by @p DataType type, however the user's @p DataType class
     * has to be derived from the TransferableQuadraturePointData class. In
     * practice that amounts to implementing the following three functions shown
     * below for a quadrature point data with 2 scalars:
     * @code
     * class MyQData : public TransferableQuadraturePointData
     * {
     * public:
     *   double elasticity_parameter_lambda;
     *   double elasticity_parameter_mu;
     *   unsigned int number_of_values() const
     *   {
     *      return 2;
     *   }
     *
     *   // a function to pack scalars into a vector
     *   void pack_values(std::vector<double> &values) const
     *   {
     *     Assert (values.size()==2, ExcInternalError());
     *     values[0] = elasticity_parameter_lambda;
     *     values[1] = elasticity_parameter_mu;
     *   }
     *
     *   void unpack_values(const std::vector<double> &values)
     *   {
     *     Assert (values.size() ==2, ExcInternalError());
     *     elasticity_parameter_lambda = values[0];
     *     elasticity_parameter_mu     = values[1];
     *   }
     * };
     * @endcode
     * Note that the order of packing and unpacking has to be the same.
     *
     * This class can then be use with CellDataStorage in the following way:
     * @code
     * CellDataStorage<typename Triangulation<dim,dim>::cell_iterator,MyQData>
     *   data_storage;
     * parallel::distributed::ContinuousQuadratureDataTransfer<dim,MyQData>
     * data_transfer(FE_Q<dim>(2),QGauss<dim>(3),QGauss<dim>(4));
     * //...populate data for all active cells in data_storage
     * //...mark cells for refinement...
     * data_transfer.prepare_for_coarsening_and_refinement(triangulation,data_storage);
     * triangulation.execute_coarsening_and_refinement();
     * //...initialize quadrature point data on new cells by calling
     * // CellDataStorage::initialize()
     * data_transfer.interpolate();
     * @endcode
     * This approach can be extended to quadrature point data with Tensor
     * objects of arbitrary order, although with a little bit more work in
     * packing and unpacking of data inside MyQData class.
     *
     * @note Currently coarsening is not supported.
     *
     * @note The functionality provided by this class can alternatively be achieved
     * using parallel::distributed::SolutionTransfer. However, that would
     * require the following steps: (i) create an auxiliary DoFHandler with a
     * (discontinuous Galerkin) FiniteElement which has enough components to
     * represent all data stored at the quadrature points; (ii) project the data
     * to the FiniteElement space and thereby store results in global vectors;
     * (iii) use parallel::distributed::SolutionTransfer to project FE vectors
     * to the new mesh; and (iv) finally project the data back to the quadrature
     * points on the new mesh via FEValues class. The
     * ContinuousQuadratureDataTransfer class aims at simplifying the whole
     * process by only requiring that the quadrature point data class is derived
     * from the TransferableQuadraturePointData. Everything else will be done
     * automatically.
     *
     * @note This class is not well suited to situations where the values stored
     * at quadrature points represent samples from a discontinuous field. An
     * example for such a situation would be where the data stored at the
     * quadrature points represents the elastic or plastic state of a material,
     * i.e., a property that varies discontinuously within a solid. In such
     * cases, trying to transfer data from the quadrature points to a finite
     * element field that is continuous (at least within one cell) will likely
     * yield over and undershoots that, once evaluated at a different set of
     * quadrature points (on child or parent cells) results in values that will
     * not make much sense.
     */
    template <int dim, typename DataType>
    class ContinuousQuadratureDataTransfer
    {
    public:
      static_assert(
        std::is_base_of_v<TransferableQuadraturePointData, DataType>,
        "User's DataType class should be derived from TransferableQuadraturePointData");

      /**
       * An alias for a cell.
       */
      using CellIteratorType =
        typename parallel::distributed::Triangulation<dim>::cell_iterator;

      /**
       * Constructor which takes the FiniteElement @p projection_fe , the quadrature
       * rule @p mass_quadrature used to integrate its local @ref GlossMassMatrix "mass matrix" and
       * finally the quadrature rule @p data_quadrature which is used to store @p DataType.
       *
       * @pre @p projection_fe has to be scalar-valued.
       *
       * @note Since this class does projection on cell-by-cell basis,
       * @p projection_fe is only required to be continuous within the cell.
       */
      ContinuousQuadratureDataTransfer(const FiniteElement<dim> &projection_fe,
                                       const Quadrature<dim> &mass_quadrature,
                                       const Quadrature<dim> &data_quadrature);

      /**
       * Prepare for coarsening and refinement of a triangulation @p tria .
       * @p data_storage represents the cell data which should be transferred
       * and it should be initialized for each locally owned active cell.
       *
       * @note Although CellDataStorage class allows storing on different cells
       * different objects derived from the base class, here we assume that
       * @p data_storage contains objects of the same type, more specifically
       * they pack/unpack the same data.
       */
      void
      prepare_for_coarsening_and_refinement(
        parallel::distributed::Triangulation<dim>   &tria,
        CellDataStorage<CellIteratorType, DataType> &data_storage);

      /**
       * Interpolate the data previously stored in this object before the mesh
       * was refined or coarsened onto the quadrature points of the currently
       * active set of cells.
       *
       * @note Before calling this function the user is expected to populate the
       * data stored in the @p data_storage object provided to prepare_for_coarsening_and_refinement()
       * at new cells using CellDataStorage::initialize(). If that is not the
       * case, an exception will be thrown in debug mode.
       */
      void
      interpolate();

    private:
      /**
       * A callback function used to pack the data on the current mesh into
       * objects that can later be retrieved after refinement, coarsening and
       * repartitioning.
       */
      std::vector<char>
      pack_function(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
                        &cell,
        const CellStatus status);

      /**
       * A callback function used to unpack the data on the current mesh that
       * has been packed up previously on the mesh before refinement,
       * coarsening and repartitioning.
       */
      void
      unpack_function(
        const typename parallel::distributed::Triangulation<dim>::cell_iterator
                        &cell,
        const CellStatus status,
        const boost::iterator_range<std::vector<char>::const_iterator>
          &data_range);

      /**
       * FiniteElement used to project data from and to quadrature points.
       */
      const std::unique_ptr<const FiniteElement<dim>> projection_fe;

      /**
       * The size of the data that will be sent, which depends on the DataType
       * class.
       */
      std::size_t data_size_in_bytes;

      /**
       * Number of quadrature points at which DataType is stored.
       */
      const unsigned int n_q_points;

      /**
       * Projection matrix from the quadrature points to local DoFs for a single
       * scalar.
       */
      FullMatrix<double> project_to_fe_matrix;

      /**
       * Projection matrix from the local DoFs to quadrature points for a single
       * scalar.
       */
      FullMatrix<double> project_to_qp_matrix;

      /**
       * Auxiliary matrix which represents projection of each internal value
       * stored
       * at the quadrature point (second index) to the local DoFs of the @p projection_fe
       * (first index).
       */
      FullMatrix<double> matrix_dofs;

      /**
       * Projection of @p matrix_dofs to each child cell in case of adaptive refinement.
       */
      FullMatrix<double> matrix_dofs_child;

      /**
       * Auxiliary matrix which represents data (second index) stored at each
       * quadrature point (first index).
       */
      FullMatrix<double> matrix_quadrature;

      /**
       * The handle that the parallel::distributed::Triangulation has assigned
       * to this object while registering the pack_callback function.
       */
      unsigned int handle;

      /**
       * A pointer to the CellDataStorage class whose data will be transferred.
       */
      CellDataStorage<CellIteratorType, DataType> *data_storage;

      /**
       * A pointer to the distributed triangulation to which cell data is
       * attached.
       */
      parallel::distributed::Triangulation<dim> *triangulation;
    };

  } // namespace distributed

} // namespace parallel

#endif

/** @} */

#ifndef DOXYGEN

// -------------------  inline and template functions ----------------

//--------------------------------------------------------------------
//                         CellDataStorage
//--------------------------------------------------------------------

template <typename CellIteratorType, typename DataType>
template <typename T>
inline void
CellDataStorage<CellIteratorType, DataType>::initialize(
  const CellIteratorType &cell,
  const unsigned int      n_q_points)
{
  static_assert(std::is_base_of_v<DataType, T>,
                "User's T class should be derived from user's DataType class");
  // The first time this method is called, it has to initialize the reference
  // to the triangulation object
  if (!tria)
    tria = &cell->get_triangulation();
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto key = cell->id();
  if (map.find(key) == map.end())
    {
      map[key] = std::vector<std::shared_ptr<DataType>>(n_q_points);
      // we need to initialize one-by-one as the std::vector<>(q, T())
      // will end with a single same T object stored in each element of the
      // vector:
      const auto it = map.find(key);
      for (unsigned int q = 0; q < n_q_points; ++q)
        it->second[q] = std::make_shared<T>();
    }
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline void
CellDataStorage<CellIteratorType, DataType>::initialize(
  const CellIteratorType                                          &cell_start,
  const typename std_cxx20::type_identity<CellIteratorType>::type &cell_end,
  const unsigned int                                               number)
{
  for (CellIteratorType it = cell_start; it != cell_end; ++it)
    if (it->is_locally_owned())
      initialize<T>(it, number);
}



template <typename CellIteratorType, typename DataType>
inline bool
CellDataStorage<CellIteratorType, DataType>::erase(const CellIteratorType &cell)
{
  const auto key = cell->id();
  const auto it  = map.find(key);
  if (it == map.end())
    return false;
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());
  for (unsigned int i = 0; i < it->second.size(); ++i)
    {
      Assert(
        it->second[i].use_count() == 1,
        ExcMessage(
          "Can not erase the cell data multiple objects reference its data."));
    }

  return (map.erase(key) == 1);
}



template <typename CellIteratorType, typename DataType>
inline void
CellDataStorage<CellIteratorType, DataType>::clear()
{
  // Do not call
  // map.clear();
  // as we want to be sure no one uses the stored objects. Loop manually:
  auto it = map.begin();
  while (it != map.end())
    {
      // loop over all objects and see if no one is using them
      for (unsigned int i = 0; i < it->second.size(); ++i)
        {
          Assert(
            it->second[i].use_count() == 1,
            ExcMessage(
              "Can not erase the cell data, multiple objects reference it."));
        }
      it = map.erase(it);
    }
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::vector<std::shared_ptr<T>>
CellDataStorage<CellIteratorType, DataType>::get_data(
  const CellIteratorType &cell)
{
  static_assert(std::is_base_of_v<DataType, T>,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  Assert(it != map.end(), ExcMessage("Could not find data for the cell"));

  // It would be nice to have a specialized version of this function for
  // T==DataType. However explicit (i.e full) specialization of a member
  // template is only allowed when the enclosing class is also explicitly (i.e
  // fully) specialized. Thus, stick with copying of shared pointers even when
  // the T==DataType:
  std::vector<std::shared_ptr<T>> res(it->second.size());
  for (unsigned int q = 0; q < res.size(); ++q)
    {
      res[q] = std::dynamic_pointer_cast<T>(it->second[q]);
      Assert(res[q], ExcCellDataTypeMismatch());
    }
  return res;
}



template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::vector<std::shared_ptr<const T>>
CellDataStorage<CellIteratorType, DataType>::get_data(
  const CellIteratorType &cell) const
{
  static_assert(std::is_base_of_v<DataType, T>,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  Assert(it != map.end(), ExcMessage("Could not find QP data for the cell"));

  // Cast base class to the desired class. This has to be done irrespectively of
  // T==DataType as we need to return shared_ptr<const T> to make sure the user
  // does not modify the content of QP objects
  std::vector<std::shared_ptr<const T>> res(it->second.size());
  for (unsigned int q = 0; q < res.size(); ++q)
    {
      res[q] = std::dynamic_pointer_cast<const T>(it->second[q]);
      Assert(res[q], ExcCellDataTypeMismatch());
    }
  return res;
}

template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::optional<std::vector<std::shared_ptr<T>>>
CellDataStorage<CellIteratorType, DataType>::try_get_data(
  const CellIteratorType &cell)
{
  static_assert(std::is_base_of_v<DataType, T>,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  if (it != map.end())
    {
      // Cast base class to the desired class. This has to be done
      // irrespectively of T==DataType as we need to return
      // shared_ptr<const T> to make sure the user
      // does not modify the content of QP objects
      std::vector<std::shared_ptr<T>> result(it->second.size());
      for (unsigned int q = 0; q < result.size(); ++q)
        {
          result[q] = std::dynamic_pointer_cast<T>(it->second[q]);
          Assert(result[q], ExcCellDataTypeMismatch());
        }
      return {result};
    }
  else
    {
      return {};
    }
}

template <typename CellIteratorType, typename DataType>
template <typename T>
inline std::optional<std::vector<std::shared_ptr<const T>>>
CellDataStorage<CellIteratorType, DataType>::try_get_data(
  const CellIteratorType &cell) const
{
  static_assert(std::is_base_of_v<DataType, T>,
                "User's T class should be derived from user's DataType class");
  Assert(&cell->get_triangulation() == tria, ExcTriangulationMismatch());

  const auto it = map.find(cell->id());
  if (it != map.end())
    {
      // Cast base class to the desired class. This has to be done
      // irrespectively of T==DataType as we need to return
      // shared_ptr<const T> to make sure the user
      // does not modify the content of QP objects
      std::vector<std::shared_ptr<const T>> result(it->second.size());
      for (unsigned int q = 0; q < result.size(); ++q)
        {
          result[q] = std::dynamic_pointer_cast<const T>(it->second[q]);
          Assert(result[q], ExcCellDataTypeMismatch());
        }
      return {result};
    }
  else
    {
      return {};
    }
}

//--------------------------------------------------------------------
//                    ContinuousQuadratureDataTransfer
//--------------------------------------------------------------------


/*
 * Pack cell data of type @p DataType stored using @p data_storage in @p cell
 * at each quadrature point to @p matrix_data. Here @p matrix_data is a matrix
 * whose first index corresponds to different quadrature points on the cell
 * whereas the second index represents different values stored at each
 * quadrature point in the DataType class.
 */
template <typename CellIteratorType, typename DataType>
inline void
pack_cell_data(const CellIteratorType                            &cell,
               const CellDataStorage<CellIteratorType, DataType> *data_storage,
               FullMatrix<double>                                &matrix_data)
{
  static_assert(std::is_base_of_v<TransferableQuadraturePointData, DataType>,
                "User's DataType class should be derived from QPData");

  if (const auto qpd = data_storage->try_get_data(cell))
    {
      const unsigned int m = qpd->size();
      Assert(m > 0, ExcInternalError());
      const unsigned int n = (*qpd)[0]->number_of_values();
      matrix_data.reinit(m, n);

      std::vector<double> single_qp_data(n);
      for (unsigned int q = 0; q < m; ++q)
        {
          (*qpd)[q]->pack_values(single_qp_data);
          AssertDimension(single_qp_data.size(), n);

          for (unsigned int i = 0; i < n; ++i)
            matrix_data(q, i) = single_qp_data[i];
        }
    }
  else
    {
      matrix_data.reinit({0, 0});
    }
}



/*
 * the opposite of the pack function above.
 */
template <typename CellIteratorType, typename DataType>
inline void
unpack_to_cell_data(const CellIteratorType                      &cell,
                    const FullMatrix<double>                    &values_at_qp,
                    CellDataStorage<CellIteratorType, DataType> *data_storage)
{
  static_assert(std::is_base_of_v<TransferableQuadraturePointData, DataType>,
                "User's DataType class should be derived from QPData");

  if (const auto qpd = data_storage->try_get_data(cell))
    {
      const unsigned int n = values_at_qp.n();
      AssertDimension((*qpd)[0]->number_of_values(), n);

      std::vector<double> single_qp_data(n);
      AssertDimension(qpd->size(), values_at_qp.m());

      for (unsigned int q = 0; q < qpd->size(); ++q)
        {
          for (unsigned int i = 0; i < n; ++i)
            single_qp_data[i] = values_at_qp(q, i);
          (*qpd)[q]->unpack_values(single_qp_data);
        }
    }
}


#  ifdef DEAL_II_WITH_P4EST

namespace parallel
{
  namespace distributed
  {
    template <int dim, typename DataType>
    inline ContinuousQuadratureDataTransfer<dim, DataType>::
      ContinuousQuadratureDataTransfer(const FiniteElement<dim> &projection_fe_,
                                       const Quadrature<dim>    &lhs_quadrature,
                                       const Quadrature<dim>    &rhs_quadrature)
      : projection_fe(
          std::unique_ptr<const FiniteElement<dim>>(projection_fe_.clone()))
      , data_size_in_bytes(0)
      , n_q_points(rhs_quadrature.size())
      , project_to_fe_matrix(projection_fe->n_dofs_per_cell(), n_q_points)
      , project_to_qp_matrix(n_q_points, projection_fe->n_dofs_per_cell())
      , handle(numbers::invalid_unsigned_int)
      , data_storage(nullptr)
      , triangulation(nullptr)
    {
      Assert(
        projection_fe->n_components() == 1,
        ExcMessage(
          "ContinuousQuadratureDataTransfer requires scalar FiniteElement"));

      FETools::compute_projection_from_quadrature_points_matrix(
        *projection_fe.get(),
        lhs_quadrature,
        rhs_quadrature,
        project_to_fe_matrix);

      FETools::compute_interpolation_to_quadrature_points_matrix(
        *projection_fe.get(), rhs_quadrature, project_to_qp_matrix);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::
      prepare_for_coarsening_and_refinement(
        parallel::distributed::Triangulation<dim>   &tr_,
        CellDataStorage<CellIteratorType, DataType> &data_storage_)
    {
      Assert(data_storage == nullptr,
             ExcMessage("This function can be called only once"));
      triangulation = &tr_;
      data_storage  = &data_storage_;

      handle = triangulation->register_data_attach(
        [this](const typename parallel::distributed::Triangulation<
                 dim>::cell_iterator &cell,
               const CellStatus       status) {
          return this->pack_function(cell, status);
        },
        /*returns_variable_size_data=*/true);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::interpolate()
    {
      triangulation->notify_ready_to_unpack(
        handle,
        [this](const typename parallel::distributed::Triangulation<
                 dim>::cell_iterator &cell,
               const CellStatus       status,
               const boost::iterator_range<std::vector<char>::const_iterator>
                 &data_range) {
          this->unpack_function(cell, status, data_range);
        });

      // invalidate the pointers
      data_storage  = nullptr;
      triangulation = nullptr;
    }



    template <int dim, typename DataType>
    inline std::vector<char>
    ContinuousQuadratureDataTransfer<dim, DataType>::pack_function(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
        &cell,
      const CellStatus /*status*/)
    {
      pack_cell_data(cell, data_storage, matrix_quadrature);

      // project to FE
      const unsigned int number_of_values = matrix_quadrature.n();
      matrix_dofs.reinit(project_to_fe_matrix.m(), number_of_values);
      if (number_of_values > 0)
        project_to_fe_matrix.mmult(matrix_dofs, matrix_quadrature);

      return Utilities::pack(matrix_dofs, /*allow_compression=*/false);
    }



    template <int dim, typename DataType>
    inline void
    ContinuousQuadratureDataTransfer<dim, DataType>::unpack_function(
      const typename parallel::distributed::Triangulation<dim>::cell_iterator
                      &cell,
      const CellStatus status,
      const boost::iterator_range<std::vector<char>::const_iterator>
        &data_range)
    {
      Assert((status != CellStatus::children_will_be_coarsened),
             ExcNotImplemented());
      (void)status;

      matrix_dofs =
        Utilities::unpack<FullMatrix<double>>(data_range.begin(),
                                              data_range.end(),
                                              /*allow_compression=*/false);
      const unsigned int number_of_values = matrix_dofs.n();
      if (number_of_values == 0)
        return;

      matrix_quadrature.reinit(n_q_points, number_of_values);

      if (cell->has_children())
        {
          // we need to first use prolongation matrix to get dofvalues on child
          // cells based on dofvalues stored in the parent's data_store
          matrix_dofs_child.reinit(projection_fe->n_dofs_per_cell(),
                                   number_of_values);
          for (unsigned int child = 0; child < cell->n_children(); ++child)
            if (cell->child(child)->is_locally_owned())
              {
                projection_fe
                  ->get_prolongation_matrix(child, cell->refinement_case())
                  .mmult(matrix_dofs_child, matrix_dofs);

                // now we do the usual business of evaluating FE on quadrature
                // points:
                project_to_qp_matrix.mmult(matrix_quadrature,
                                           matrix_dofs_child);

                // finally, put back into the map:
                unpack_to_cell_data(cell->child(child),
                                    matrix_quadrature,
                                    data_storage);
              }
        }
      else
        {
          // if there are no children, evaluate FE field at
          // rhs_quadrature points.
          project_to_qp_matrix.mmult(matrix_quadrature, matrix_dofs);

          // finally, put back into the map:
          unpack_to_cell_data(cell, matrix_quadrature, data_storage);
        }
    }

  } // namespace distributed

} // namespace parallel

#  endif // DEAL_II_WITH_P4EST

#endif // DOXYGEN
DEAL_II_NAMESPACE_CLOSE

#endif
