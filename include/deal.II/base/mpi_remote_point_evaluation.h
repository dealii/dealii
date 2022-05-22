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

#ifndef dealii_mpi_mpi_remote_point_evaluation_h
#define dealii_mpi_mpi_remote_point_evaluation_h

#include <deal.II/base/config.h>

#include <deal.II/dofs/dof_handler.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    /**
     * Helper class to access values on non-matching grids.
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
       * Constructor.
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
       */
      RemotePointEvaluation(
        const double       tolerance                              = 1e-6,
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
             const Mapping<dim, spacedim> &      mapping);

      /**
       * Data of points positioned in a cell.
       */
      struct CellData
      {
        /**
         * Level and index of cells.
         */
        std::vector<std::pair<int, int>> cells;

        /**
         * Pointers to beginning and ending of the (reference) points
         * associated to the cell.
         */
        std::vector<unsigned int> reference_point_ptrs;

        /**
         * Reference points in the interval [0,1]^dim.
         */
        std::vector<Point<dim>> reference_point_values;
      };

      /**
       * Evaluate function @p evaluation_function in the given  points and
       * triangulation. The result is stored in @p output.
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
       */
      template <typename T>
      void
      evaluate_and_process(
        std::vector<T> &output,
        std::vector<T> &buffer,
        const std::function<void(const ArrayView<T> &, const CellData &)>
          &evaluation_function) const;

      /**
       * This method is the inverse of the method evaluate_and_process(). It
       * makes the data at the points, provided by @p input, available in the
       * function @p evaluation_function.
       *
       * @warning This is a collective call that needs to be executed by all
       *   processors in the communicator.
       */
      template <typename T>
      void
      process_and_evaluate(
        const std::vector<T> &input,
        std::vector<T> &      buffer,
        const std::function<void(const ArrayView<const T> &, const CellData &)>
          &evaluation_function) const;

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

    private:
      /**
       * Tolerance to be used while determining the surrounding cells of a
       * point.
       */
      const double tolerance;

      /**
       * Enforce unique mapping, i.e., (one-to-one) relation of points and
       * cells.
       */
      const bool enforce_unique_mapping;

      /**
       * RTree level to be used during the construction of the bounding boxes.
       */
      const unsigned int rtree_level;

      /**
       * Function that marks relevant vertices to make search of active cells
       * around point more efficient.
       */
      const std::function<std::vector<bool>()> marked_vertices;

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
      SmartPointer<const Triangulation<dim, spacedim>> tria;

      /**
       * Reference to the Mapping object used during reinit().
       */
      SmartPointer<const Mapping<dim, spacedim>> mapping;

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
       * of degrees of freedoms) needs to performed only once per cell.
       */
      CellData cell_data;

      /**
       * Permutation index within a send buffer.
       */
      std::vector<unsigned int> send_permutation;

      /**
       * Ranks to send to.
       */
      std::vector<unsigned int> send_ranks;

      /**
       * Pointers of ranges within a send buffer to be sent to the ranks
       * specified by send_ranks.
       */
      std::vector<unsigned int> send_ptrs;
    };


    template <int dim, int spacedim>
    template <typename T>
    void
    RemotePointEvaluation<dim, spacedim>::evaluate_and_process(
      std::vector<T> &output,
      std::vector<T> &buffer,
      const std::function<void(const ArrayView<T> &, const CellData &)>
        &evaluation_function) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)output;
      (void)buffer;
      (void)evaluation_function;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_communicator());

      output.resize(point_ptrs.back());
      buffer.resize(send_permutation.size() * 2);
      ArrayView<T> buffer_1(buffer.data(), buffer.size() / 2);
      ArrayView<T> buffer_2(buffer.data() + buffer.size() / 2,
                            buffer.size() / 2);

      // evaluate functions at points
      evaluation_function(buffer_1, cell_data);

      // sort for communication
      for (unsigned int i = 0; i < send_permutation.size(); ++i)
        buffer_2[send_permutation[i]] = buffer_1[i];

      // process remote quadrature points and send them away
      std::map<unsigned int, std::vector<char>> temp_map;

      std::vector<MPI_Request> requests;
      requests.reserve(send_ranks.size());

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_communicator());

      std::map<unsigned int, std::vector<T>> temp_recv_map;

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            {
              // process locally-owned values
              temp_recv_map[my_rank] =
                std::vector<T>(buffer_2.begin() + send_ptrs[i],
                               buffer_2.begin() + send_ptrs[i + 1]);
              continue;
            }

          temp_map[send_ranks[i]] =
            Utilities::pack(std::vector<T>(buffer_2.begin() + send_ptrs[i],
                                           buffer_2.begin() + send_ptrs[i + 1]),
                            false);

          auto &buffer = temp_map[send_ranks[i]];

          requests.push_back(MPI_Request());

          const int ierr = MPI_Isend(buffer.data(),
                                     buffer.size(),
                                     MPI_CHAR,
                                     send_ranks[i],
                                     internal::Tags::remote_point_evaluation,
                                     tria->get_communicator(),
                                     &requests.back());
          AssertThrowMPI(ierr);
        }

      for (const auto recv_rank : recv_ranks)
        {
          if (recv_rank == my_rank)
            continue;

          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int message_length;
          ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
          AssertThrowMPI(ierr);

          std::vector<char> buffer(message_length);

          ierr = MPI_Recv(buffer.data(),
                          buffer.size(),
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          internal::Tags::remote_point_evaluation,
                          tria->get_communicator(),
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          temp_recv_map[status.MPI_SOURCE] =
            Utilities::unpack<std::vector<T>>(buffer, false);
        }

      // make sure all messages have been sent
      const int ierr =
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      // copy received data into output vector
      auto it = recv_permutation.begin();
      for (const auto &j : temp_recv_map)
        for (const auto &i : j.second)
          {
            output[*it] = i;
            it++;
          }
#endif
    }


    template <int dim, int spacedim>
    template <typename T>
    void
    RemotePointEvaluation<dim, spacedim>::process_and_evaluate(
      const std::vector<T> &input,
      std::vector<T> &      buffer,
      const std::function<void(const ArrayView<const T> &, const CellData &)>
        &evaluation_function) const
    {
#ifndef DEAL_II_WITH_MPI
      Assert(false, ExcNeedsMPI());
      (void)input;
      (void)buffer;
      (void)evaluation_function;
#else
      static CollectiveMutex      mutex;
      CollectiveMutex::ScopedLock lock(mutex, tria->get_communicator());

      const auto &ptr = this->get_point_ptrs();

      std::map<unsigned int, std::vector<T>> temp_recv_map;

      for (unsigned int i = 0; i < recv_ranks.size(); ++i)
        temp_recv_map[recv_ranks[i]].resize(recv_ptrs[i + 1] - recv_ptrs[i]);

      const unsigned int my_rank =
        Utilities::MPI::this_mpi_process(tria->get_communicator());

#  ifdef DEBUG
      {
        unsigned int i = 0;

        for (auto &j : temp_recv_map)
          i += j.second.size();

        AssertDimension(recv_permutation.size(), i);
      }
#  endif

      {
        // duplicate data to be able to sort it more easily in the next step
        std::vector<T> buffer_(ptr.back());
        for (unsigned int i = 0, c = 0; i < ptr.size() - 1; ++i)
          {
            const auto n_entries = ptr[i + 1] - ptr[i];

            for (unsigned int j = 0; j < n_entries; ++j, ++c)
              buffer_[c] = input[i];
          }

        // sort data according to the ranks
        auto it = recv_permutation.begin();
        for (auto &j : temp_recv_map)
          for (auto &i : j.second)
            {
              i = buffer_[*it];
              it++;
            }
      }

      // buffer.resize(point_ptrs.back());
      buffer.resize(send_permutation.size() * 2);
      ArrayView<T> buffer_1(buffer.data(), buffer.size() / 2);
      ArrayView<T> buffer_2(buffer.data() + buffer.size() / 2,
                            buffer.size() / 2);

      // process remote quadrature points and send them away
      std::map<unsigned int, std::vector<char>> temp_map;

      std::vector<MPI_Request> requests;
      requests.reserve(recv_ranks.size());

      for (const auto recv_rank : recv_ranks)
        {
          if (recv_rank == my_rank)
            continue;

          temp_map[recv_rank] =
            Utilities::pack(temp_recv_map[recv_rank], false);

          auto &buffer_send = temp_map[recv_rank];

          requests.push_back(MPI_Request());

          const int ierr = MPI_Isend(buffer_send.data(),
                                     buffer_send.size(),
                                     MPI_CHAR,
                                     recv_rank,
                                     internal::Tags::remote_point_evaluation,
                                     tria->get_communicator(),
                                     &requests.back());
          AssertThrowMPI(ierr);
        }

      for (unsigned int i = 0; i < send_ranks.size(); ++i)
        {
          if (send_ranks[i] == my_rank)
            {
              const auto &buffer_send = temp_recv_map[send_ranks[i]];
              // process locally-owned values
              const unsigned int j = std::distance(send_ranks.begin(),
                                                   std::find(send_ranks.begin(),
                                                             send_ranks.end(),
                                                             my_rank));

              AssertDimension(buffer_send.size(),
                              send_ptrs[j + 1] - send_ptrs[j]);

              for (unsigned int i = send_ptrs[j], c = 0; i < send_ptrs[j + 1];
                   ++i, ++c)
                buffer_1[i] = buffer_send[c];

              continue;
            }

          MPI_Status status;
          int        ierr = MPI_Probe(MPI_ANY_SOURCE,
                               internal::Tags::remote_point_evaluation,
                               tria->get_communicator(),
                               &status);
          AssertThrowMPI(ierr);

          int message_length;
          ierr = MPI_Get_count(&status, MPI_CHAR, &message_length);
          AssertThrowMPI(ierr);

          std::vector<char> recv_buffer(message_length);

          ierr = MPI_Recv(recv_buffer.data(),
                          recv_buffer.size(),
                          MPI_CHAR,
                          status.MPI_SOURCE,
                          internal::Tags::remote_point_evaluation,
                          tria->get_communicator(),
                          MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          const auto recv_buffer_unpacked =
            Utilities::unpack<std::vector<T>>(recv_buffer, false);

          auto ptr =
            std::find(send_ranks.begin(), send_ranks.end(), status.MPI_SOURCE);

          Assert(ptr != send_ranks.end(), ExcNotImplemented());

          const unsigned int j = std::distance(send_ranks.begin(), ptr);

          AssertDimension(recv_buffer_unpacked.size(),
                          send_ptrs[j + 1] - send_ptrs[j]);

          for (unsigned int i = send_ptrs[j], c = 0; i < send_ptrs[j + 1];
               ++i, ++c)
            {
              AssertIndexRange(i, buffer_1.size());
              AssertIndexRange(c, recv_buffer_unpacked.size());
              buffer_1[i] = recv_buffer_unpacked[c];
            }
        }

      const int ierr =
        MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
      AssertThrowMPI(ierr);

      // sort for easy access during function call
      for (unsigned int i = 0; i < send_permutation.size(); ++i)
        buffer_2[i] = buffer_1[send_permutation[i]];

      // evaluate function at points
      evaluation_function(buffer_2, cell_data);
#endif
    }

  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
