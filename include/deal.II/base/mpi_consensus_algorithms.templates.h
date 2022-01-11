// ---------------------------------------------------------------------
//
// Copyright (C) 2011 - 2021 by the deal.II authors
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
      namespace
      {
        /**
         * Return whether a vector of targets (MPI ranks) has only unique
         * elements.
         *
         * This function is only used within assertions, which causes GCC
         * to issue a warning in release mode that due to -Werror then causes an
         * error. We suppress this by using the [[gnu::unused]] error (because
         * the
         * [[maybe_unused]] attribute is only supported from C++17 forward).
         *
         * Unfortunately, in contrast to what the standard says, the Microsoft
         * compiler does not ignore the gnu::unused attribute as it should,
         * and then produces an error of its own. So we disable the attribute
         * for that compiler.
         */
#ifndef DEAL_II_MSVC
        [[gnu::unused]]
#endif
        bool
        has_unique_elements(const std::vector<unsigned int> &targets)
        {
          std::vector<unsigned int> my_destinations = targets;
          std::sort(my_destinations.begin(), my_destinations.end());
          return (std::adjacent_find(my_destinations.begin(),
                                     my_destinations.end()) ==
                  my_destinations.end());
        }
      } // namespace



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
        , job_supports_mpi(Utilities::MPI::job_supports_mpi())
        , my_rank(job_supports_mpi ? this_mpi_process(comm) : 0)
        , n_procs(job_supports_mpi ? n_mpi_processes(comm) : 1)
      {}



      template <typename T1, typename T2>
      NBX<T1, T2>::NBX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      NBX<T1, T2>::run()
      {
        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, this->comm);

        // 1) Send data to identified targets and start receiving
        //    the answers from these very same processes.
        start_communication();

        // 2) Until all posted receive operations are known to have completed,
        //    answer requests and keep checking whether all requests of
        //    this process have been answered.
        //
        //    The requests that we catch in the answer_requests() function
        //    originate elsewhere, that is, they are not in response
        //    to our own messages
        //
        //    Note also that we may not catch all incoming requests in
        //    the following two lines: our own requests may have been
        //    satisfied before we've dealt with all incoming requests.
        //    That's ok: We will get around to dealing with all remaining
        //    message later. We just want to move on to the next step
        //    as early as possible.
        while (all_locally_originated_receives_are_completed() == false)
          maybe_answer_one_request();

        // 3) Signal to all other processes that all requests of this process
        //    have been answered
        signal_finish();

        // 4) Nevertheless, this process has to keep on answering (potential)
        //    incoming requests until all processes have received the
        //    answer to all requests
        while (all_remotely_originated_receives_are_completed() == false)
          maybe_answer_one_request();

        // 5) process the answer to all requests
        clean_up_and_end_communication();

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::start_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // 1)
        targets = this->process.compute_targets();
        Assert(has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        const auto n_targets = targets.size();

        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_answer_request;

        // 2) allocate memory
        send_requests.resize(n_targets);
        send_buffers.resize(n_targets);

        {
          // 4) send and receive
          for (unsigned int index = 0; index < n_targets; ++index)
            {
              const unsigned int rank = targets[index];
              AssertIndexRange(rank,
                               Utilities::MPI::n_mpi_processes(this->comm));

              auto &send_buffer = send_buffers[index];
              this->process.create_request(rank, send_buffer);

              // Post a request to send data
              auto ierr = MPI_Isend(send_buffer.data(),
                                    send_buffer.size() * sizeof(T1),
                                    MPI_BYTE,
                                    rank,
                                    tag_request,
                                    this->comm,
                                    &send_requests[index]);
              AssertThrowMPI(ierr);
            }
        }
#endif
      }



      template <typename T1, typename T2>
      bool
      NBX<T1, T2>::all_locally_originated_receives_are_completed()
      {
#ifdef DEAL_II_WITH_MPI
        // We know that all requests have come in when we have pending
        // messages from all targets with the right tag. We can check
        // for pending messages with MPI_IProbe, which returns
        // immediately with a return code that indicates whether
        // it has found a message from a given process with a given
        // tag
        for (const unsigned int target : targets)
          {
            const int tag_deliver = Utilities::MPI::internal::Tags::
              consensus_algorithm_nbx_process_deliver;

            int        request_is_pending;
            const auto ierr = MPI_Iprobe(target,
                                         tag_deliver,
                                         this->comm,
                                         &request_is_pending,
                                         MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            // If there is no pending message from that process,
            // then we are clearly not done receiving everything
            // yet -- so return false:
            if (request_is_pending == 0)
              return false;
          }

        // If we have made it here, then we have received an answer
        // from everyone and can return true:
        return true;

#else
        return true;
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::maybe_answer_one_request()
      {
#ifdef DEAL_II_WITH_MPI

        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_nbx_process_deliver;

        // Check if there is a request pending. By selecting the
        // tag_request tag, these are other processes asking for
        // our own replies, not these other processes' replies
        // to our own requests.
        //
        // There may be multiple such pending messages. We
        // only answer one.
        MPI_Status status;
        int        request_is_pending;
        const auto ierr = MPI_Iprobe(MPI_ANY_SOURCE,
                                     tag_request,
                                     this->comm,
                                     &request_is_pending,
                                     &status);
        AssertThrowMPI(ierr);

        if (request_is_pending != 0)
          {
            // Get the rank of the requesting process and add it to the
            // list of requesting processes (which may contain duplicates).
            const auto other_rank = status.MPI_SOURCE;

            Assert(requesting_processes.find(other_rank) ==
                     requesting_processes.end(),
                   ExcMessage("Process is requesting a second time!"));
            requesting_processes.insert(other_rank);

            // get size of incoming message
            int  number_amount;
            auto ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
            AssertThrowMPI(ierr);

            // allocate memory for incoming message
            Assert(number_amount % sizeof(T1) == 0, ExcInternalError());
            std::vector<T1> buffer_recv(number_amount / sizeof(T1));
            ierr = MPI_Recv(buffer_recv.data(),
                            number_amount,
                            MPI_BYTE,
                            other_rank,
                            tag_request,
                            this->comm,
                            MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            // Allocate memory for an answer message to the current request,
            // and ask the 'process' object to produce an answer:
            request_buffers.emplace_back(std::make_unique<std::vector<T2>>());
            auto &request_buffer = *request_buffers.back();
            this->process.answer_request(other_rank,
                                         buffer_recv,
                                         request_buffer);

            // Then initiate sending the answer back to the requester.
            request_requests.emplace_back(std::make_unique<MPI_Request>());
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
      NBX<T1, T2>::all_remotely_originated_receives_are_completed()
      {
#ifdef DEAL_II_WITH_MPI
        int        all_ranks_reached_barrier;
        const auto ierr = MPI_Test(&barrier_request,
                                   &all_ranks_reached_barrier,
                                   MPI_STATUSES_IGNORE);
        AssertThrowMPI(ierr);
        return all_ranks_reached_barrier != 0;
#else
        return true;
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

          int ierr = MPI_Wait(&barrier_request, MPI_STATUS_IGNORE);
          AssertThrowMPI(ierr);

          for (auto &i : request_requests)
            {
              ierr = MPI_Wait(i.get(), MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr);
            }

#  ifdef DEBUG
          // note: IBarrier seems to make problem during testing, this
          // additional Barrier seems to help
          ierr = MPI_Barrier(this->comm);
          AssertThrowMPI(ierr);
#  endif
        }

        // We know from the various calls to MPI_Iprobe that all of our
        // requests have received answers, but we have not actually
        // gotten the data from MPI. Do so and unpack the data.
        {
          for (const unsigned int target : targets)
            {
              const int tag_deliver = Utilities::MPI::internal::Tags::
                consensus_algorithm_nbx_process_deliver;

              // Just to be sure, double check that the message really
              // is already here -- in other words, that the following
              // MPI_Recv is going to return immediately. But as part of
              // this, we also use the status object to query the size
              // of the message so that we can resize the receive buffer.
              int        request_is_pending;
              MPI_Status status;
              {
                const int ierr = MPI_Iprobe(target,
                                            tag_deliver,
                                            this->comm,
                                            &request_is_pending,
                                            &status);
                AssertThrowMPI(ierr);
              }

              (void)request_is_pending;
              Assert(request_is_pending, ExcInternalError());

              // OK, so yes, a message is here. Receive it.
              int message_size;
              {
                const int ierr =
                  MPI_Get_count(&status, MPI_BYTE, &message_size);
                AssertThrowMPI(ierr);
              }
              Assert(message_size % sizeof(T2) == 0, ExcInternalError());
              std::vector<T2> recv_buffer(message_size / sizeof(T2));

              {
                const int ierr = MPI_Recv(recv_buffer.data(),
                                          recv_buffer.size() * sizeof(T2),
                                          MPI_BYTE,
                                          target,
                                          tag_deliver,
                                          this->comm,
                                          MPI_STATUS_IGNORE);
                AssertThrowMPI(ierr);
              }

              this->process.read_answer(target, recv_buffer);
            }
        }
#endif
      }



      template <typename T1, typename T2>
      PEX<T1, T2>::PEX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      PEX<T1, T2>::run()
      {
        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, this->comm);

        // 1) Send requests and start receiving the answers.
        //    In particular, determine how many requests we should expect
        //    on the current process.
        const unsigned int n_requests = start_communication();

        // 2) Answer requests:
        for (unsigned int request = 0; request < n_requests; ++request)
          answer_one_request(request);

        // 3) Process answers:
        clean_up_and_end_communication();

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename T1, typename T2>
      unsigned int
      PEX<T1, T2>::start_communication()
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;


        // 1) determine with which processes this process wants to communicate
        // with
        targets = this->process.compute_targets();
        Assert(has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        // 2) determine who wants to communicate with this process
        sources =
          compute_point_to_point_communication_pattern(this->comm, targets);

        const unsigned int n_targets = targets.size();
        const unsigned int n_sources = sources.size();

        // 2) allocate memory
        recv_buffers.resize(n_targets);
        send_buffers.resize(n_targets);
        send_request_and_recv_answer_requests.resize(2 * n_targets);

        send_answer_requests.resize(n_sources);
        requests_buffers.resize(n_sources);

        // 4) send and receive
        for (unsigned int i = 0; i < n_targets; ++i)
          {
            const unsigned int rank = targets[i];
            AssertIndexRange(rank, Utilities::MPI::n_mpi_processes(this->comm));

            // pack data which should be sent
            auto &send_buffer = send_buffers[i];
            this->process.create_request(rank, send_buffer);

            // start to send data
            auto ierr =
              MPI_Isend(send_buffer.data(),
                        send_buffer.size() * sizeof(T1),
                        MPI_BYTE,
                        rank,
                        tag_request,
                        this->comm,
                        &send_request_and_recv_answer_requests[n_targets + i]);
            AssertThrowMPI(ierr);

            // Post the operation that receives the answers
            auto &recv_buffer = recv_buffers[i];
            this->process.prepare_buffer_for_answer(rank, recv_buffer);
            ierr = MPI_Irecv(recv_buffer.data(),
                             recv_buffer.size() * sizeof(T2),
                             MPI_BYTE,
                             rank,
                             tag_deliver,
                             this->comm,
                             &send_request_and_recv_answer_requests[i]);
            AssertThrowMPI(ierr);
          }

        return sources.size();
#else
        return 0;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::answer_one_request(const unsigned int index)
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;

        // Wait until we have a message ready for retrieval, though we don't
        // care which process it is from. We know that the source must be
        // listed in the 'sources' array, though.
        MPI_Status status;
        int ierr = MPI_Probe(MPI_ANY_SOURCE, tag_request, this->comm, &status);
        AssertThrowMPI(ierr);

        // Get rank of incoming message and verify that it makes sense
        const unsigned int other_rank = status.MPI_SOURCE;
        Assert(std::find(sources.begin(), sources.end(), other_rank) !=
                 sources.end(),
               ExcInternalError());

        Assert(requesting_processes.find(other_rank) ==
                 requesting_processes.end(),
               ExcMessage(
                 "A process is sending a request after a request from "
                 "the same process has previously already been "
                 "received. This algorithm does not expect this to happen."));
        requesting_processes.insert(other_rank);

        std::vector<T1> buffer_recv;

        // Actually get the incoming message:
        int number_amount;
        ierr = MPI_Get_count(&status, MPI_BYTE, &number_amount);
        AssertThrowMPI(ierr);
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

        // Process request by asking the user-provided function for
        // the answer and post a send for it.
        auto &request_buffer = requests_buffers[index];
        this->process.answer_request(other_rank, buffer_recv, request_buffer);

        ierr = MPI_Isend(request_buffer.data(),
                         request_buffer.size() * sizeof(T2),
                         MPI_BYTE,
                         other_rank,
                         tag_deliver,
                         this->comm,
                         &send_answer_requests[index]);
        AssertThrowMPI(ierr);
#else
        (void)index;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::clean_up_and_end_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // Finalize all MPI_Request objects for both the
        // send-request and receive-answer operations.
        if (send_request_and_recv_answer_requests.size() > 0)
          {
            const int ierr =
              MPI_Waitall(send_request_and_recv_answer_requests.size(),
                          send_request_and_recv_answer_requests.data(),
                          MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
          }

        // Then also check the send-answer requests.
        if (send_answer_requests.size() > 0)
          {
            const int ierr = MPI_Waitall(send_answer_requests.size(),
                                         send_answer_requests.data(),
                                         MPI_STATUSES_IGNORE);
            AssertThrowMPI(ierr);
          }

        // We now know that all answers to the requests we have sent
        // have been received and put in their respective buffers.
        // Pass them on to the user-provided functions:
        for (unsigned int i = 0; i < targets.size(); ++i)
          this->process.read_answer(targets[i], recv_buffers[i]);
#endif
      }



      template <typename T1, typename T2>
      Serial<T1, T2>::Serial(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Serial<T1, T2>::run()
      {
        const auto targets = this->process.compute_targets();

        // The only valid target for a serial program is itself.
        if (targets.size() != 0)
          {
            Assert(targets.size() == 1,
                   ExcMessage(
                     "On a single process, the only valid target "
                     "is process zero (the process itself), which can only be "
                     "listed once."));
            AssertDimension(targets[0], 0);

            // Since the caller indicates that there is a target, and since we
            // know that it is the current process, let the process send
            // something to itself.
            std::vector<T1> send_buffer;
            std::vector<T2> recv_buffer;
            std::vector<T2> request_buffer;

            this->process.create_request(0, send_buffer);
            this->process.prepare_buffer_for_answer(0, recv_buffer);
            this->process.answer_request(0, send_buffer, request_buffer);
            recv_buffer = request_buffer;
            this->process.read_answer(0, recv_buffer);
          }

        return targets; // nothing to do
      }



      template <typename T1, typename T2>
      Selector<T1, T2>::Selector(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {
        // Depending on the number of processes we switch between
        // implementations. We reduce the threshold for debug mode to be
        // able to test also the non-blocking implementation. This feature
        // is tested by:
        // tests/multigrid/transfer_matrix_free_06.with_mpi=true.with_p4est=true.with_trilinos=true.mpirun=10.output
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
#    ifdef DEBUG
        if (this->n_procs > 10)
#    else
        if (this->n_procs > 99)
#    endif
          consensus_algo.reset(new NBX<T1, T2>(process, comm));
        else
#  endif
#endif
          if (this->n_procs > 1)
          consensus_algo.reset(new PEX<T1, T2>(process, comm));
        else
          consensus_algo.reset(new Serial<T1, T2>(process, comm));
      }



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Selector<T1, T2>::run()
      {
        return consensus_algo->run();
      }


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
