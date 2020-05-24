// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2020 by the deal.II authors
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

#ifndef dealii_mpi_consensus_algorithm_templates_h
#define dealii_mpi_consensus_algorithm_templates_h

#include <deal.II/base/config.h>

#include <deal.II/base/exceptions.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi_consensus_algorithms.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/tensor.h>

#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/lapack_full_matrix.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/vector.h>

#include <vector>

DEAL_II_NAMESPACE_OPEN

namespace Utilities
{
  namespace MPI
  {
    namespace ConsensusAlgorithms
    {
      template <typename T1, typename T2>
      void
      Process<T1, T2>::answer_request(const unsigned int,
                                      const std::vector<T1> &,
                                      std::vector<T2> &)
      {
        // nothing to do
      }



      template <typename T1, typename T2>
      void
      Process<T1, T2>::create_request(const unsigned int, std::vector<T1> &)
      {
        // nothing to do
      }



      template <typename T1, typename T2>
      void
      Process<T1, T2>::prepare_buffer_for_answer(const unsigned int,
                                                 std::vector<T2> &)
      {
        // nothing to do
      }



      template <typename T1, typename T2>
      void
      Process<T1, T2>::read_answer(const unsigned int, const std::vector<T2> &)
      {
        // nothing to do
      }



      template <typename T1, typename T2>
      Interface<T1, T2>::Interface(Process<T1, T2> &process,
                                   const MPI_Comm & comm)
        : process(process)
        , comm(comm)
        , my_rank(this_mpi_process(comm))
        , n_procs(n_mpi_processes(comm))
      {}



      template <typename T1, typename T2>
      NBX<T1, T2>::NBX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::run()
      {
        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, this->comm);

        // 1) send requests and start receiving the answers
        start_communication();

        // 2) answer requests and check if all requests of this process have
        // been
        //    answered
        while (!check_own_state())
          answer_requests();

        // 3) signal to all other processes that all requests of this process
        // have
        //    been answered
        signal_finish();

        // 4) nevertheless, this process has to keep on answering (potential)
        //    incoming requests until all processes have received the
        //    answer to all requests
        while (!check_global_state())
          answer_requests();

        // 5) process the answer to all requests
        clean_up_and_end_communication();
      }



      template <typename T1, typename T2>
      bool
      NBX<T1, T2>::check_own_state()
      {
#ifdef DEAL_II_WITH_MPI
        int        all_receive_requests_are_done;
        const auto ierr = MPI_Testall(recv_requests.size(),
                                      recv_requests.data(),
                                      &all_receive_requests_are_done,
                                      MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);

        return all_receive_requests_are_done;
#else
        return true;
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::signal_finish()
      {
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
        const auto ierr = MPI_Ibarrier(this->comm, &barrier_request);
        AssertThrowMPI(ierr);
#  else
        AssertThrow(
          false,
          ExcMessage(
            "ConsensusAlgorithms::NBX uses MPI 3.0 features. You should compile with at least MPI 3.0."));
#  endif
#endif
      }



      template <typename T1, typename T2>
      bool
      NBX<T1, T2>::check_global_state()
      {
#ifdef DEAL_II_WITH_MPI
        int        all_ranks_reached_barrier;
        const auto ierr = MPI_Test(&barrier_request,
                                   &all_ranks_reached_barrier,
                                   MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);
        return all_ranks_reached_barrier;
#else
        return true;
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::answer_requests()
      {
#ifdef DEAL_II_WITH_MPI

        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_process_deliver;

        // check if there is a request pending
        MPI_Status status;
        int        request_is_pending;
        const auto ierr = MPI_Iprobe(MPI_ANY_SOURCE,
                                     tag_request,
                                     this->comm,
                                     &request_is_pending,
                                     &status);
        AssertThrowMPI(ierr);

        if (request_is_pending) // request is pending
          {
            // get rank of requesting process
            const auto other_rank = status.MPI_SOURCE;

#  ifdef DEBUG
            Assert(requesting_processes.find(other_rank) ==
                     requesting_processes.end(),
                   ExcMessage("Process is requesting a second time!"));
            requesting_processes.insert(other_rank);
#  endif

            std::vector<T1> buffer_recv;
            // get size of incoming message
            int  number_amount;
            auto ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
            AssertThrowMPI(ierr);

            // allocate memory for incoming message
            Assert(number_amount % sizeof(T1) == 0, ExcInternalError());
            buffer_recv.resize(number_amount / sizeof(T1));
            ierr = MPI_Recv(buffer_recv.data(),
                            number_amount,
                            MPI_BYTE,
                            other_rank,
                            tag_request,
                            this->comm,
                            &status);
            AssertThrowMPI(ierr);

            // allocate memory for answer message
            request_buffers.emplace_back(std::make_unique<std::vector<T2>>());
            request_requests.emplace_back(std::make_unique<MPI_Request>());

            // process request
            auto &request_buffer = *request_buffers.back();
            this->process.answer_request(other_rank,
                                         buffer_recv,
                                         request_buffer);

            // start to send answer back
            ierr = MPI_Isend(request_buffer.data(),
                             request_buffer.size() * sizeof(T2),
                             MPI_BYTE,
                             other_rank,
                             tag_deliver,
                             this->comm,
                             request_requests.back().get());
            AssertThrowMPI(ierr);
          }
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::start_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // 1)
        targets              = this->process.compute_targets();
        const auto n_targets = targets.size();

        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_process_deliver;

        // 2) allocate memory
        recv_buffers.resize(n_targets);
        recv_requests.resize(n_targets);
        send_requests.resize(n_targets);
        send_buffers.resize(n_targets);

        {
          // 4) send and receive
          for (unsigned int i = 0; i < n_targets; i++)
            {
              const unsigned int rank  = targets[i];
              const unsigned int index = i;

              // translate index set to a list of pairs
              auto &send_buffer = send_buffers[index];
              this->process.create_request(rank, send_buffer);

              // start to send data
              auto ierr = MPI_Isend(send_buffer.data(),
                                    send_buffer.size() * sizeof(T1),
                                    MPI_BYTE,
                                    rank,
                                    tag_request,
                                    this->comm,
                                    &send_requests[index]);
              AssertThrowMPI(ierr);

              // start to receive data
              auto &recv_buffer = recv_buffers[index];
              this->process.prepare_buffer_for_answer(rank, recv_buffer);
              ierr = MPI_Irecv(recv_buffer.data(),
                               recv_buffer.size() * sizeof(T2),
                               MPI_BYTE,
                               rank,
                               tag_deliver,
                               this->comm,
                               &recv_requests[index]);
              AssertThrowMPI(ierr);
            }
        }
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::clean_up_and_end_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // clean up
        {
          if (send_requests.size() > 0)
            {
              const int ierr = MPI_Waitall(send_requests.size(),
                                           send_requests.data(),
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }

          if (recv_requests.size() > 0)
            {
              const int ierr = MPI_Waitall(recv_requests.size(),
                                           recv_requests.data(),
                                           MPI_STATUSES_IGNORE);
              AssertThrowMPI(ierr);
            }


          const int ierr = MPI_Wait(&barrier_request, MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          for (auto &i : request_requests)
            {
              const auto ierr = MPI_Wait(i.get(), MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr);
            }

#  ifdef DEBUG
          // note: IBarrier seems to make problem during testing, this
          // additional Barrier seems to help
          MPI_Barrier(this->comm);
#  endif
        }

        // unpack data
        {
          for (unsigned int i = 0; i < targets.size(); i++)
            this->process.read_answer(targets[i], recv_buffers[i]);
        }
#endif
      }



      template <typename T1, typename T2>
      PEX<T1, T2>::PEX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::run()
      {
        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, this->comm);

        // 1) send requests and start receiving the answers
        //    especially determine how many requests are expected
        const unsigned int n_requests = start_communication();

        // 2) answer requests
        for (unsigned int request = 0; request < n_requests; request++)
          answer_requests(request);

        // 3) process answers
        clean_up_and_end_communication();
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::answer_requests(int index)
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;

        MPI_Status status;
        auto ierr = MPI_Probe(MPI_ANY_SOURCE, tag_request, this->comm, &status);
        AssertThrowMPI(ierr);

        // get rank of incoming message
        const auto other_rank = status.MPI_SOURCE;

        std::vector<T1> buffer_recv;

        // get size of incoming message
        int number_amount;
        ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
        AssertThrowMPI(ierr);

        // allocate memory for incoming message
        Assert(number_amount % sizeof(T1) == 0, ExcInternalError());
        buffer_recv.resize(number_amount / sizeof(T1));
        ierr = MPI_Recv(buffer_recv.data(),
                        number_amount,
                        MPI_BYTE,
                        other_rank,
                        tag_request,
                        this->comm,
                        &status);
        AssertThrowMPI(ierr);

        // process request
        auto &request_buffer = requests_buffers[index];
        this->process.answer_request(other_rank, buffer_recv, request_buffer);

        // start to send answer back
        ierr = MPI_Isend(request_buffer.data(),
                         request_buffer.size() * sizeof(T2),
                         MPI_BYTE,
                         other_rank,
                         tag_deliver,
                         this->comm,
                         &requests_answers[index]);
        AssertThrowMPI(ierr);
#else
        (void)index;
#endif
      }



      template <typename T1, typename T2>
      unsigned int
      PEX<T1, T2>::start_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // 1) determine with which processes this process wants to communicate
        targets = this->process.compute_targets();

        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;

        // 2) determine who wants to communicate with this process
        sources =
          compute_point_to_point_communication_pattern(this->comm, targets);

        const auto n_targets = targets.size();
        const auto n_sources = sources.size();

        // 2) allocate memory
        recv_buffers.resize(n_targets);
        send_buffers.resize(n_targets);
        send_and_recv_buffers.resize(2 * n_targets);

        requests_answers.resize(n_sources);
        requests_buffers.resize(n_sources);

        // 4) send and receive
        for (unsigned int i = 0; i < n_targets; i++)
          {
            const unsigned int rank = targets[i];

            // pack data which should be sent
            auto &send_buffer = send_buffers[i];
            this->process.create_request(rank, send_buffer);

            // start to send data
            auto ierr = MPI_Isend(send_buffer.data(),
                                  send_buffer.size() * sizeof(T1),
                                  MPI_BYTE,
                                  rank,
                                  tag_request,
                                  this->comm,
                                  &send_and_recv_buffers[n_targets + i]);
            AssertThrowMPI(ierr);

            // start to receive data
            auto &recv_buffer = recv_buffers[i];
            this->process.prepare_buffer_for_answer(rank, recv_buffer);
            ierr = MPI_Irecv(recv_buffer.data(),
                             recv_buffer.size() * sizeof(T2),
                             MPI_BYTE,
                             rank,
                             tag_deliver,
                             this->comm,
                             &send_and_recv_buffers[i]);
            AssertThrowMPI(ierr);
          }

        return sources.size();
#else
        return 0;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::clean_up_and_end_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // finalize all MPI_Requests
        if (send_and_recv_buffers.size() > 0)
          {
            auto ierr = MPI_Waitall(send_and_recv_buffers.size(),
                                    send_and_recv_buffers.data(),
                                    MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
          }

        if (requests_answers.size() > 0)
          {
            auto ierr = MPI_Waitall(requests_answers.size(),
                                    requests_answers.data(),
                                    MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
          }

        // unpack received data
        for (unsigned int i = 0; i < targets.size(); i++)
          this->process.read_answer(targets[i], recv_buffers[i]);
#endif
      }



      template <typename T1, typename T2>
      Selector<T1, T2>::Selector(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {
        // Depending on the number of processes we switch between
        // implementations. We reduce the threshold for debug mode to be able to
        // test also the non-blocking implementation. This feature is tested by:
        // tests/multigrid/transfer_matrix_free_06.with_mpi=true.with_p4est=true.with_trilinos=true.mpirun=15.output
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
#    ifdef DEBUG
        if (Utilities::MPI::n_mpi_processes(comm) > 14)
#    else
        if (Utilities::MPI::n_mpi_processes(comm) > 99)
#    endif
          consensus_algo.reset(new NBX<T1, T2>(process, comm));
        else
#  endif
#endif
          consensus_algo.reset(new PEX<T1, T2>(process, comm));
      }



      template <typename T1, typename T2>
      void
      Selector<T1, T2>::run()
      {
        consensus_algo->run();
      }


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
