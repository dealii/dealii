// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2021 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_mpi_remote_point_evaluation_h
#define dealii_mpi_mpi_remote_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_tags.h>

#include <deal.II/dofs/dof_handler.h>

#include <boost/signals2/connection.hpp>


DEAL_II_NAMESPACE_OPEN

// Forward declarations
namespace GridTools
{
  template <int dim, int spacedim>
  class Cache;

  namespace internal
  {
    template <int dim, int spacedim>
    struct DistributedComputePointLocationsInternal;
  }
} // namespace GridTools

namespace Utilities
{
  namespace MPI
  {
    /**
     * @brief Communicate values between a mesh and arbitrary points
     *
     * This class enables evaluation of finite element functions defined
     * on a Triangulation at arbitrary points inside the mesh. The underlying
     * finite element mesh can be distributed among processes, which makes
     * the operations more involved due to communication.
     *
     * Once initialized, the data exchange is handled by the functions
     * process_and_evaluate(), that transfers data from the Triangulation
     * to the set of points, and evaluate_and_process(), that performs
     * the inverse operation (taking data defined at each point and making
     * it available in the corresponding cell of the Triangulation).
     *
     * The helper functions VectorTools::point_values() and
     * VectorTools::point_gradients (the overloads taking an instance
     * of RemotePointEvaluation) can be used in place of
     * process_and_evaluate() and perform the evaluation of a finite
     * element function.
     *
     * The template argument `DataType` of the process_and_evaluate() and
     * evaluate_and_process() functions can be a scalar (`float`, `double`),
     * `Tensor<rank,dim>`, or `Point<dim>`. Additionally, one can specify
     * the number of components (in case `DataType` is a scalar).
     *
     * @note Usage and implementation details of this class and the
     *   helper functions above are explained in detail in step-87.
     *
     * @note The name of the fields are chosen with the method
     *   evaluate_and_process() in mind. Here, quantities are
     *   computed at specified arbitrary positioned points (and even on remote
     *   processes in the MPI universe) cell by cell and these values are sent
     *   to requesting processes, which receive the result and resort the
     *   result according to the points.
     */
    template <int dim, int spacedim = dim>
    class RemotePointEvaluation
    {
    public:
      /**
       * AdditionalData structure that can be used to tweak parameters
       * of RemotePointEvaluation.
       */
      struct AdditionalData
      {
      public:
        /**
         *  Constructor.
         */
        AdditionalData(
          const double       tolerance                              = 1e-6,
          const bool         enforce_unique_mapping                 = false,
          const unsigned int rtree_level                            = 0,
          const std::function<std::vector<bool>()> &marked_vertices = {});

        /**
         * Tolerance in terms of unit cell coordinates for determining all cells
         * around a point passed to RemotePointEvaluation during reinit().
         * Depending on the problem, it might be necessary to adjust the
         * tolerance in order to be able to identify a cell. Floating point
         * arithmetic implies that a point will, in general, not lie exactly on
         * a vertex, edge, or face.
         */
        double tolerance;

        /**
         * Enforce unique mapping, i.e., (one-to-one) relation of points and
         * cells.
         */
        bool enforce_unique_mapping;

        /**
         * RTree level to be used during the construction of the bounding boxes.
         */
        unsigned int rtree_level;

        /**
         * Function that marks relevant vertices to make search of active cells
         * around point more efficient.
         */
        std::function<std::vector<bool>()> marked_vertices;
      };

      /**
       * Constructor.
       *
       * @param additional_data Configure options for RemotePointEvaluation.
       */
      RemotePointEvaluation(
        const AdditionalData &additional_data = AdditionalData());

      /**
       * Constructor. This constructor is deprecated. Use the other constructor
       * taking AdditionalData instead.
       *
       * @param tolerance Tolerance in terms of unit cell coordinates for
       *   determining all cells around a point passed to the class during
       *   reinit(). Depending on the problem, it might be necessary to adjust
       *   the tolerance in order to be able to identify a cell.
       *   Floating point arithmetic implies that a point will, in general, not
       *   lie exactly on a vertex, edge, or face.
       * @param enforce_unique_mapping Enforce unique mapping, i.e.,
       *   (one-to-one) relation of points and cells.
       * @param rtree_level RTree level to be used during the construction of the bounding boxes.
       * @param marked_vertices Function that marks relevant vertices to make search
       *   of active cells around point more efficient.
       *
       * @deprecated
       */
      DEAL_II_DEPRECATED_WITH_COMMENT(
        "Use the constructor with AdditionalData struct.")
      RemotePointEvaluation(
        const double       tolerance,
        const bool         enforce_unique_mapping                 = false,
        const unsigned int rtree_level                            = 0,
        const std::function<std::vector<bool>()> &marked_vertices = {});

      /**
       * Destructor.
       */
      ~RemotePointEvaluation();

      /**
       * Set up internal data structures and communication pattern based on
       * a list of points @p points and mesh description (@p tria and @p
       * mapping).
       *
       * @warning This is a collective call that needs to be executed by all
       *   processors in the communicator.
       *
       * @note If you want to be sure that all points have been found, call
       *   all_points_found() after calling this function.
       */
      void
      reinit(const std::vector<Point<spacedim>> &points,
             const Triangulation<dim, spacedim> &tria,
             const Mapping<dim, spacedim>       &mapping);

      /**
       * Set up internal data structures and communication pattern based on
       * an existing Cache and a list of points @p points.
       *
       * @warning This is a collective call that needs to be executed by all
       *   processors in the communicator.
       *
       * @note If you want to be sure that all points have been found, call
       *   all_points_found() after calling this function.
       */
      void
      reinit(const GridTools::Cache<dim, spacedim> &cache,
             const std::vector<Point<spacedim>>    &points);

      /**
       * Set up internal data structures and communication pattern based on
       * GridTools::internal::DistributedComputePointLocationsInternal.
       *
       * This function is called internally by the reinit() function above.
       * Having it as a separate function makes it possible to set up the class
       * if it is known in which cells corresponding reference points are
       * located (e.g. if intersections of cells are known).
       */
      void
      reinit(const GridTools::internal::
               DistributedComputePointLocationsInternal<dim, spacedim> &data,
             const Triangulation<dim, spacedim>                        &tria,
             const Mapping<dim, spacedim> &mapping);

      /**
       * Helper class to store and to access data of points positioned in
       * processed cells.
       */
      class CellData
      {
      public:
        /**
         * Constructor.
         *
         * @param triangulation Triangulation of the domain.
         */
        CellData(const Triangulation<dim, spacedim> &triangulation);

        /**
         * Return an object that can be thought of as an array containing all
         * indices from zero (inclusive) to the total number of cells
         * where points reside given by `cells.size()` (exclusive). This allows
         * one to write code using range-based `for` loops of the following
         * kind:
         *
         * @code
         * for (const auto cell : cell_data.cell_indices())
         *   {
         *     const auto cell_dofs =
         *       cell_data.get_active_cell_iterator(cell)->as_dof_handler_iterator(
         *         dof_handler);
         *
         *     const auto unit_points = cell_data.get_unit_points(cell);
         *     const auto local_value = cell_data.get_data_view(cell, values);
         *
         *     // user code: not shown
         *   }
         * @endcode
         *
         * Here, we are looping over all cells where points reside and
         * use the index to call the functions get_active_cell_iterator(),
         * get_unit_points(), and get_data_view().
         */
        std_cxx20::ranges::iota_view<unsigned int, unsigned int>
        cell_indices() const;

        /**
         * Return active cell iterator of the processed cell @p cell.
         */
        typename Triangulation<dim, spacedim>::active_cell_iterator
        get_active_cell_iterator(const unsigned int cell) const;

        /**
         * Return unit points of the processed cell @p cell.
         */
        ArrayView<const Point<dim>>
        get_unit_points(const unsigned int cell) const;

        /**
         * Return local view of the processed cell @p cell for the vector @p values.
         */
        template <typename DataType>
        ArrayView<DataType>
        get_data_view(const unsigned int         cell,
                      const ArrayView<DataType> &values) const;

        /**
         * Level and index of processed cells.
         */
        std::vector<std::pair<int, int>> cells;

        /**
         * Pointers to the start and end of the (reference) points
         * associated to the cell.
         */
        std::vector<unsigned int> reference_point_ptrs;

        /**
         * Reference points in the interval [0,1]^dim.
         */
        std::vector<Point<dim>> reference_point_values;

      private:
        /**
         * Reference to the underlying triangulation needed for
         * get_active_cell_iterator().
         */
        const Triangulation<dim, spacedim> &triangulation;
      };

      /**
       * Return internal data structure storing the data of points
       * positioned in cells.
       */
      const CellData &
      get_cell_data() const;

      /**
       * Evaluate function @p evaluation_function in the given points and
       * triangulation. The result is stored in @p output.
       *
       * @p buffer is a temporary buffer than can be used to avoid extra
       * allocations.
       *
       * @note If the map of points to cells is not a
       *   one-to-one relation (is_map_unique()==false), the result needs to be
       *   processed with the help of get_point_ptrs(). This
       *   might be the case if a point coincides with a geometric entity (e.g.,
       *   vertex) that is shared by multiple cells or a point is outside of the
       *   computational domain.
       *
       * @warning This is a collective call that needs to be executed by all
       *   processors in the communicator.
       *
       * @note This function can handle arbitrary types including scalar and
       *   Tensor objects. In the case that scalar types are used, one can
       *   specify the number of components @p n_components. This allows to
       *   provide unrolled tensors, which is useful, e.g., if its dimension
       *   and its rank is not known at compile time.
       */
      template <typename DataType, unsigned int n_components = 1>
      void
      evaluate_and_process(
        std::vector<DataType> &output,
        std::vector<DataType> &buffer,
        const std::function<void(const ArrayView<DataType> &, const CellData &)>
                  &evaluation_function,
        const bool sort_data = true) const;

      /**
       * Same as above but with the result provided as return value and
       * without external allocation of a user-provided buffer.
       */
      template <typename DataType, unsigned int n_components = 1>
      std::vector<DataType>
      evaluate_and_process(
        const std::function<void(const ArrayView<DataType> &, const CellData &)>
                  &evaluation_function,
        const bool sort_data = true) const;

      /**
       * This method is the inverse of the method evaluate_and_process(). It
       * makes the data at the points, provided by @p input, available in the
       * function @p evaluation_function.
       *
       * @warning This is a collective call that needs to be executed by all
       *   processors in the communicator.
       *
       * @note This function can handle arbitrary types including scalar and
       *   Tensor objects. In the case that scalar types are used, one can
       *   specify the number of components @p n_components. This allows to
       *   provide unrolled tensors, which is useful, e.g., if its dimension
       *   and its rank is not known at compile time.
       */
      template <typename DataType, unsigned int n_components = 1>
      void
      process_and_evaluate(
        const std::vector<DataType>                 &input,
        std::vector<DataType>                       &buffer,
        const std::function<void(const ArrayView<const DataType> &,
                                 const CellData &)> &evaluation_function,
        const bool                                   sort_data = true) const;

      /**
       * Same as above but without external allocation of a user-provided
       * buffer.
       */
      template <typename DataType, unsigned int n_components = 1>
      void
      process_and_evaluate(
        const std::vector<DataType>                 &input,
        const std::function<void(const ArrayView<const DataType> &,
                                 const CellData &)> &evaluation_function,
        const bool                                   sort_data = true) const;

      /**
       * Return a CRS-like data structure to determine the position of the
       * result corresponding a point and the amount.
       */
      const std::vector<unsigned int> &
      get_point_ptrs() const;

      /**
       * Return if points and cells have a one-to-one relation. This is not the
       * case if a point is not owned by any cell (the point is outside of the
       * domain) or if multiple cells own the point (the point is positioned
       * on a geometric entity shared by neighboring cells).
       */
      bool
      is_map_unique() const;

      /**
       * Return if all points could be found in the domain.
       */
      bool
      all_points_found() const;

      /**
       * Return if point @p i could be found in the domain.
       */
      bool
      point_found(const unsigned int i) const;

      /**
       * Return the Triangulation object used during reinit().
       */
      const Triangulation<dim, spacedim> &
      get_triangulation() const;

      /**
       * Return the Mapping object used during reinit().
       */
      const Mapping<dim, spacedim> &
      get_mapping() const;

      /**
       * Return if the internal data structures have been set up and if yes
       * whether they are still valid (and not invalidated due to changes of the
       * Triangulation).
       */
      bool
      is_ready() const;

      /**
       * Return permutation needed for sending. If this is applied, evaluation
       * results do not have to be sorted according to the receiving ranks
       * during evaluate_and_process().
       */
      const std::vector<unsigned int> &
      get_send_permutation() const;

      /**
       * Data is contiguous rank by rank directly after MPI communication.
       * If no sorting is requested during evaluate_and_process(),
       * this permutation allows to access all values corresponding to an
       * evaluation point (sorted into different cells).
       */
      const std::vector<unsigned int> &
      get_inverse_recv_permutation() const;


    private:
      /**
       * Additional data with basic settings.
       */
      const AdditionalData additional_data;

      /**
       * Storage for the status of the triangulation signal.
       */
      boost::signals2::connection tria_signal;

      /**
       * Flag indicating if the reinit() function has been called and if yes
       * the triangulation has not been modified since then (potentially
       * invalidating the communication pattern).
       */
      bool ready_flag;

      /**
       * Reference to the Triangulation object used during reinit().
       */
      ObserverPointer<const Triangulation<dim, spacedim>> tria;

      /**
       * Reference to the Mapping object used during reinit().
       */
      ObserverPointer<const Mapping<dim, spacedim>> mapping;

      /**
       * (One-to-one) relation of points and cells.
       */
      bool unique_mapping;

      /**
       * Cache if all points passed in during reinit() have been found.
       */
      bool all_points_found_flag;

      /**
       * Since for each point multiple or no results can be available, the
       * pointers in this vector indicate the first and last entry associated
       * with a point in a CRS-like fashion.
       */
      std::vector<unsigned int> point_ptrs;

      /**
       * Permutation index within a recv buffer.
       */
      std::vector<unsigned int> recv_permutation;

      /**
       * Inverse of permutation index within a recv buffer.
       */
      std::vector<unsigned int> recv_permutation_inv;

      /**
       * Pointers of ranges within a receive buffer that are filled by ranks
       * specified by recv_ranks.
       */
      std::vector<unsigned int> recv_ptrs;

      /**
       * Ranks from where data is received.
       */
      std::vector<unsigned int> recv_ranks;

      /**
       * Point data sorted according to cells so that evaluation (incl. reading
       * of degrees of freedoms) needs to be performed only once per cell.
       */
      std::unique_ptr<CellData> cell_data;

      /**
       * Permutation index within a send buffer.
       */
      std::vector<unsigned int> send_permutation;

      /**
       * Inverse of permutation index within a send buffer.
       */
      std::vector<unsigned int> send_permutation_inv;

      /**
       * Ranks to send to.
       */
      std::vector<unsigned int> send_ranks;

      /**
       * Pointers of ranges within a send buffer to be sent to the ranks
       * specified by send_ranks.
       */
      std::vector<unsigned int> send_ptrs;

      /**
       * Buffer size (if internal sorting is requested) determined as max of:
       *
       * needed during evaluate_and_process()
       * - point data (sorted according to cells) -> input from cell loop
       * - point data (sorted according to ranks) -> to be sent
       * - memory for receiving data (by one process at a time)
       *
       * needed during process_and_evaluate()
       * - point data (sorted according to cells) -> output to cell loop
       * - point data (sorted according to ranks) -> to be sent
       * - memory for receiving data (by one process at a time)
       */
      unsigned int buffer_size_with_sorting;

      /**
       * Buffer size (if sorting is not requested); corresponds to the
       * number of evaluation points.
       */
      unsigned int buffer_size_without_sorting;
    };

    namespace internal
    {
#ifdef DEAL_II_WITH_MPI
      /**
       * Pack @p data and send it via MPI_Isend.
       */
      template <typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == false, void>
      pack_and_isend(const ArrayView<const T>       &data,
                     const unsigned int              rank,
                     const unsigned int              tag,
                     const MPI_Comm                  comm,
                     std::vector<std::vector<char>> &buffers,
                     std::vector<MPI_Request>       &requests)
      {
        requests.emplace_back(MPI_Request());

        buffers.emplace_back(Utilities::pack(
          std::vector<T>(data.data(), data.data() + data.size()), false));

        const int ierr = MPI_Isend(buffers.back().data(),
                                   buffers.back().size(),
                                   MPI_CHAR,
                                   rank,
                                   tag,
                                   comm,
                                   &requests.back());
        AssertThrowMPI(ierr);
      }



      /**
       * Above function specialized for data types supported by MPI
       * so that one can skip packing.
       */
      template <typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == true, void>
      pack_and_isend(const ArrayView<const T> &data,
                     const unsigned int        rank,
                     const unsigned int        tag,
                     const MPI_Comm            comm,
                     std::vector<std::vector<char>> & /*buffers*/,
                     std::vector<MPI_Request> &requests)
      {
        requests.emplace_back(MPI_Request());

        const int ierr = MPI_Isend(data.data(),
                                   data.size(),
                                   Utilities::MPI::mpi_type_id_for_type<T>,
                                   rank,
                                   tag,
                                   comm,
                                   &requests.back());
        AssertThrowMPI(ierr);
      }



      /**
       * Above function specialized for Tenors objects. The underlying data type
       * might be supported by MPI so that one can skip packing.
       */
      template <int rank_, int dim, typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == true, void>
      pack_and_isend(const ArrayView<const Tensor<rank_, dim, T>> &data,
                     const unsigned int                            rank,
                     const unsigned int                            tag,
                     const MPI_Comm                                comm,
                     std::vector<std::vector<char>>               &buffers,
                     std::vector<MPI_Request>                     &requests)
      {
        ArrayView<const T> data_(reinterpret_cast<const T *>(data.data()),
                                 data.size() * Utilities::pow(dim, rank_));

        pack_and_isend(data_, rank, tag, comm, buffers, requests);
      }


      /**
       * Receive message, unpack it, and store the result in @p data.
       */
      template <typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == false, void>
      recv_and_unpack(const ArrayView<T> &data,
                      const MPI_Comm      comm,
                      const MPI_Status   &status,
                      std::vector<char>  &buffer)
      {
        int message_length;
        int ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
        AssertThrowMPI(ierr);

        buffer.resize(message_length);

        ierr = MPI_Recv(buffer.data(),
                        buffer.size(),
                        MPI_CHAR,
                        status.MPI_SOURCE,
                        internal::Tags::remote_point_evaluation,
                        comm,
                        MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);

        // unpack data
        const auto temp = Utilities::unpack<std::vector<T>>(buffer, false);

        for (unsigned int i = 0; i < data.size(); ++i)
          data[i] = temp[i];
      }



      /**
       * Above function specialized for data types supported by MPI
       * so that one can skip unpacking.
       */
      template <typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == true, void>
      recv_and_unpack(const ArrayView<T> &data,
                      const MPI_Comm      comm,
                      const MPI_Status   &status,
                      std::vector<char> & /*buffer*/)
      {
        const auto ierr = MPI_Recv(data.data(),
                                   data.size(),
                                   Utilities::MPI::mpi_type_id_for_type<T>,
                                   status.MPI_SOURCE,
                                   internal::Tags::remote_point_evaluation,
                                   comm,
                                   MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
      }



      /**
       * Above function specialized for Tensor objects. The underlying data type
       * might be supported by MPI so that one can skip unpacking.
       */
      template <int rank_, int dim, typename T>
      std::enable_if_t<Utilities::MPI::is_mpi_type<T> == true, void>
      recv_and_unpack(const ArrayView<Tensor<rank_, dim, T>> &data,
                      const MPI_Comm                          comm,
                      const MPI_Status                       &status,
                      std::vector<char>                      &buffer)
      {
        const ArrayView<T> data_(reinterpret_cast<T *>(data.data()),
                                 data.size() * Utilities::pow(dim, rank_));

        recv_and_unpack(data_, comm, status, buffer);
      }
#endif
    } // namespace internal



    template <int dim, int spacedim>
    template <typename DataType>
    ArrayView<DataType>
    RemotePointEvaluation<dim, spacedim>::CellData::get_data_view(
      const unsigned int         cell,
      const ArrayView<DataType> &values) const
    {
      AssertIndexRange(cell, cells.size());
      return {values.data() + reference_point_ptrs[cell],
              reference_point_ptrs[cell + 1] - reference_point_ptrs[cell]};
    }



    template <int dim, int spacedim>
    template <typename DataType, unsigned int n_components>
    void
    RemotePointEvaluation<dim, spacedim>::evaluate_and_process(
      std::vector<DataType> &output,
      std::vector<DataType> &buffer,
      const std::function<void(const ArrayView<DataType> &, const CellData &)>
                &evaluation_function,
      const bool sort_data) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)output;
      (void)buffer;
      (void)evaluation_function;
      (void)sort_data;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_mpi_communicator());

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_mpi_communicator());

      // allocate memory for output and buffer
      output.resize(point_ptrs.back() * n_components);

      buffer.resize(
        (sort_data ? buffer_size_with_sorting : buffer_size_without_sorting) *
        n_components);

      // ... for evaluation
      ArrayView<DataType> buffer_eval(buffer.data(),
                                      send_permutation.size() * n_components);

      // ... for communication (send)
      ArrayView<DataType> buffer_send(
        sort_data ? (buffer.data() + send_permutation.size() * n_components) :
                    buffer.data(),
        send_permutation.size() * n_components);

      // more arrays
      std::vector<MPI_Request>       send_requests;
      std::vector<std::vector<char>> send_buffers_packed;
      std::vector<char>              recv_buffer_packed;

      // evaluate functions at points
      evaluation_function(buffer_eval, *cell_data);

      // sort for communication (optional)
      if (sort_data)
        {
          const auto my_rank_local_recv_ptr =
            std::find(recv_ranks.begin(), recv_ranks.end(), my_rank);

          if (my_rank_local_recv_ptr != recv_ranks.end())
            {
              const unsigned int my_rank_local_recv =
                std::distance(recv_ranks.begin(), my_rank_local_recv_ptr);
              const unsigned int my_rank_local_send = std::distance(
                send_ranks.begin(),
                std::find(send_ranks.begin(), send_ranks.end(), my_rank));
              const unsigned int  start = send_ptrs[my_rank_local_send];
              const unsigned int  end   = send_ptrs[my_rank_local_send + 1];
              const unsigned int *recv_ptr =
                recv_permutation.data() + recv_ptrs[my_rank_local_recv];
              for (unsigned int i = 0; i < send_permutation.size(); ++i)
                {
                  const unsigned int send_index = send_permutation[i];

                  if (start <= send_index && send_index < end)
                    // local data -> can be copied to output directly
                    for (unsigned int c = 0; c < n_components; ++c)
                      output[recv_ptr[send_index - start] * n_components + c] =
                        buffer_eval[i * n_components + c];
                  else
                    // data to be sent
                    for (unsigned int c = 0; c < n_components; ++c)
                      buffer_send[send_index * n_components + c] =
                        buffer_eval[i * n_components + c];
                }
            }
          else
            {
              for (unsigned int i = 0; i < send_permutation.size(); ++i)
                for (unsigned int c = 0; c < n_components; ++c)
                  buffer_send[send_permutation[i] * n_components + c] =
                    buffer_eval[i * n_components + c];
            }
        }

      // send data
      send_buffers_packed.reserve(send_ranks.size());
      send_requests.reserve(send_ranks.size());

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            continue;

          internal::pack_and_isend(
            ArrayView<const DataType>(
              buffer_send.begin() + send_ptrs[i] * n_components,
              (send_ptrs[i + 1] - send_ptrs[i]) * n_components),
            send_ranks[i],
            internal::Tags::remote_point_evaluation,
            tria->get_mpi_communicator(),
            send_buffers_packed,
            send_requests);
        }

      // copy local data directly to output if no sorting is requested
      if (!sort_data)
        {
          const auto my_rank_local_recv_ptr =
            std::find(recv_ranks.begin(), recv_ranks.end(), my_rank);

          if (my_rank_local_recv_ptr != recv_ranks.end())
            {
              const unsigned int my_rank_local_recv =
                std::distance(recv_ranks.begin(), my_rank_local_recv_ptr);
              const unsigned int my_rank_local_send = std::distance(
                send_ranks.begin(),
                std::find(send_ranks.begin(), send_ranks.end(), my_rank));

              for (unsigned int j = recv_ptrs[my_rank_local_recv],
                                k = send_ptrs[my_rank_local_send];
                   j < recv_ptrs[my_rank_local_recv + 1];
                   ++j, ++k)
                for (unsigned int c = 0; c < n_components; ++c)
                  output[j * n_components + c] =
                    buffer_eval[k * n_components + c];
            }
        }

      // receive data
      for (unsigned int i = 0; i < recv_ranks.size(); ++i)
        {
          if (recv_ranks[i] == my_rank)
            continue;

          // receive remote data
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_mpi_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          const auto ptr =
            std::find(recv_ranks.begin(), recv_ranks.end(), status.MPI_SOURCE);

          Assert(ptr != recv_ranks.end(), ExcNotImplemented());

          const unsigned int j = std::distance(recv_ranks.begin(), ptr);

          // ... for communication (recv)
          ArrayView<DataType> recv_buffer(
            sort_data ?
              (buffer.data() + send_permutation.size() * 2 * n_components) :
              (output.data() + recv_ptrs[j] * n_components),
            (recv_ptrs[j + 1] - recv_ptrs[j]) * n_components);

          internal::recv_and_unpack(recv_buffer,
                                    tria->get_mpi_communicator(),
                                    status,
                                    recv_buffer_packed);

          if (sort_data)
            {
              // sort data into output vector (optional)
              for (unsigned int i = recv_ptrs[j], k = 0; i < recv_ptrs[j + 1];
                   ++i, ++k)
                for (unsigned int c = 0; c < n_components; ++c)
                  output[recv_permutation[i] * n_components + c] =
                    recv_buffer[k * n_components + c];
            }
        }

      // make sure all messages have been sent
      if (!send_requests.empty())
        {
          const int ierr = MPI_Waitall(send_requests.size(),
                                       send_requests.data(),
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }
#endif
    }


    template <int dim, int spacedim>
    template <typename DataType, unsigned int n_components>
    std::vector<DataType>
    RemotePointEvaluation<dim, spacedim>::evaluate_and_process(
      const std::function<void(const ArrayView<DataType> &, const CellData &)>
                &evaluation_function,
      const bool sort_data) const
    {
      std::vector<DataType> output;
      std::vector<DataType> buffer;

      this->evaluate_and_process<DataType, n_components>(output,
                                                         buffer,
                                                         evaluation_function,
                                                         sort_data);

      return output;
    }



    template <int dim, int spacedim>
    template <typename DataType, unsigned int n_components>
    void
    RemotePointEvaluation<dim, spacedim>::process_and_evaluate(
      const std::vector<DataType>                 &input,
      std::vector<DataType>                       &buffer,
      const std::function<void(const ArrayView<const DataType> &,
                               const CellData &)> &evaluation_function,
      const bool                                   sort_data) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)input;
      (void)buffer;
      (void)evaluation_function;
      (void)sort_data;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_mpi_communicator());

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_mpi_communicator());

      // allocate memory for buffer
      const auto &point_ptrs = this->get_point_ptrs();

      AssertDimension(input.size(),
                      sort_data ? ((point_ptrs.size() - 1) * n_components) :
                                  (point_ptrs.back() * n_components));
      buffer.resize(
        (sort_data ? buffer_size_with_sorting : buffer_size_without_sorting) *
        n_components);

      // ... for evaluation
      ArrayView<DataType> buffer_eval(buffer.data(),
                                      send_permutation.size() * n_components);

      // ... for communication (send)
      ArrayView<DataType> buffer_send(
        sort_data ? (buffer.data() + send_permutation.size() * n_components) :
                    const_cast<DataType *>(input.data()),
        point_ptrs.back() * n_components);

      // more arrays
      std::vector<MPI_Request>       send_requests;
      std::vector<std::vector<char>> send_buffers_packed;
      std::vector<char>              recv_buffer_packed;

      // sort for communication (and duplicate data if necessary)

      if (sort_data)
        {
          const auto my_rank_local_recv_ptr =
            std::find(recv_ranks.begin(), recv_ranks.end(), my_rank);

          if (my_rank_local_recv_ptr != recv_ranks.end())
            {
              // optimize the case if we have our own rank among the possible
              // list
              const unsigned int my_rank_local_recv =
                std::distance(recv_ranks.begin(), my_rank_local_recv_ptr);
              const unsigned int my_rank_local_send = std::distance(
                send_ranks.begin(),
                std::find(send_ranks.begin(), send_ranks.end(), my_rank));

              const unsigned int  start = recv_ptrs[my_rank_local_recv];
              const unsigned int  end   = recv_ptrs[my_rank_local_recv + 1];
              const unsigned int *send_ptr =
                send_permutation_inv.data() + send_ptrs[my_rank_local_send];
              for (unsigned int i = 0, k = 0; i < point_ptrs.size() - 1; ++i)
                {
                  const unsigned int next = point_ptrs[i + 1];
                  for (unsigned int j = point_ptrs[i]; j < next; ++j, ++k)
                    {
                      const unsigned int recv_index = recv_permutation_inv[k];

                      // local data -> can be copied to final buffer directly
                      if (start <= recv_index && recv_index < end)
                        for (unsigned int c = 0; c < n_components; ++c)
                          buffer_eval[send_ptr[recv_index - start] *
                                        n_components +
                                      c] = input[i * n_components + c];
                      else // data to be sent
                        for (unsigned int c = 0; c < n_components; ++c)
                          buffer_send[recv_index * n_components + c] =
                            input[i * n_components + c];
                    }
                }
            }
          else
            {
              for (unsigned int i = 0, k = 0; i < point_ptrs.size() - 1; ++i)
                for (unsigned int j = point_ptrs[i]; j < point_ptrs[i + 1];
                     ++j, ++k)
                  for (unsigned int c = 0; c < n_components; ++c)
                    buffer_send[recv_permutation_inv[k] * n_components + c] =
                      input[i * n_components + c];
            }
        }

      // send data
      send_buffers_packed.reserve(recv_ranks.size());
      send_requests.reserve(recv_ranks.size());

      for (unsigned int i = 0; i < recv_ranks.size(); ++i)
        {
          if (recv_ranks[i] == my_rank)
            continue;

          internal::pack_and_isend(
            ArrayView<const DataType>(
              buffer_send.begin() + recv_ptrs[i] * n_components,
              (recv_ptrs[i + 1] - recv_ptrs[i]) * n_components),
            recv_ranks[i],
            internal::Tags::remote_point_evaluation,
            tria->get_mpi_communicator(),
            send_buffers_packed,
            send_requests);
        }

      // copy local data directly to output if no sorting is requested
      if (!sort_data)
        {
          const auto my_rank_local_recv_ptr =
            std::find(recv_ranks.begin(), recv_ranks.end(), my_rank);

          if (my_rank_local_recv_ptr != recv_ranks.end())
            {
              const unsigned int my_rank_local_recv =
                std::distance(recv_ranks.begin(), my_rank_local_recv_ptr);
              const unsigned int my_rank_local_send = std::distance(
                send_ranks.begin(),
                std::find(send_ranks.begin(), send_ranks.end(), my_rank));

              for (unsigned int j = recv_ptrs[my_rank_local_recv],
                                k = send_ptrs[my_rank_local_send];
                   j < recv_ptrs[my_rank_local_recv + 1];
                   ++j, ++k)
                for (unsigned int c = 0; c < n_components; ++c)
                  buffer_eval[k * n_components + c] =
                    input[j * n_components + c];
            }
        }

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            continue;

          // receive remote data
          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_mpi_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          // write data into buffer vector
          const auto ptr =
            std::find(send_ranks.begin(), send_ranks.end(), status.MPI_SOURCE);

          Assert(ptr != send_ranks.end(), ExcNotImplemented());

          const unsigned int j = std::distance(send_ranks.begin(), ptr);

          ArrayView<DataType> recv_buffer(
            sort_data ?
              (buffer.data() +
               (point_ptrs.back() + send_permutation.size()) * n_components) :
              (buffer.data() + send_ptrs[j] * n_components),
            (send_ptrs[j + 1] - send_ptrs[j]) * n_components);

          internal::recv_and_unpack(recv_buffer,
                                    tria->get_mpi_communicator(),
                                    status,
                                    recv_buffer_packed);

          if (sort_data)
            {
              for (unsigned int i = send_ptrs[j], k = 0; i < send_ptrs[j + 1];
                   ++i, ++k)
                for (unsigned int c = 0; c < n_components; ++c)
                  buffer_eval[send_permutation_inv[i] * n_components + c] =
                    recv_buffer[k * n_components + c];
            }
        }

      if (!send_requests.empty())
        {
          const int ierr = MPI_Waitall(send_requests.size(),
                                       send_requests.data(),
                                       MPI_STATUSES_IGNORE);
          AssertThrowMPI(ierr);
        }

      // evaluate function at points
      evaluation_function(buffer_eval, *cell_data);
#endif
    }



    template <int dim, int spacedim>
    template <typename DataType, unsigned int n_components>
    void
    RemotePointEvaluation<dim, spacedim>::process_and_evaluate(
      const std::vector<DataType>                 &input,
      const std::function<void(const ArrayView<const DataType> &,
                               const CellData &)> &evaluation_function,
      const bool                                   sort_data) const
    {
      std::vector<DataType> buffer;
      this->process_and_evaluate<DataType, n_components>(input,
                                                         buffer,
                                                         evaluation_function,
                                                         sort_data);
    }

  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
