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
      Process<T1, T2>::read_answer(const unsigned int, const std::vector<T2> &)
      {
        // nothing to do
      }



      template <typename T1, typename T2>
      Interface<T1, T2>::Interface(Process<T1, T2> &process,
                                   const MPI_Comm & comm)
        : process(&process)
        , comm(comm)
      {}



      template <typename T1, typename T2>
      Interface<T1, T2>::Interface()
        : process(nullptr)
        , comm(MPI_COMM_NULL)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Interface<T1, T2>::run()
      {
        Assert(process != nullptr,
               ExcMessage("This function can only be called if the "
                          "deprecated non-default constructor of this class "
                          "has previously been called to set the Process "
                          "object and a communicator."));
        return run(*process, comm);
      }



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Interface<T1, T2>::run(Process<T1, T2> &process, const MPI_Comm &comm)
      {
        // Unpack the 'process' object and call the function that takes
        // function objects for all operations.
        return run(
          process.compute_targets(),
          /* create_request: */
          [&process](const unsigned int target) {
            std::vector<T1> request;
            process.create_request(target, request);
            return request;
          },
          /* answer_request: */
          [&process](const unsigned int     source,
                     const std::vector<T1> &request) {
            std::vector<T2> answer;
            process.answer_request(source, request, answer);
            return answer;
          },
          /* process_answer: */
          [&process](const unsigned int target, const std::vector<T2> &answer) {
            process.read_answer(target, answer);
          },
          comm);
      }



      template <typename T1, typename T2>
      NBX<T1, T2>::NBX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      NBX<T1, T2>::run(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &create_request,
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &answer_request,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
        Assert(has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, comm);

        // 1) Send data to identified targets and start receiving
        //    the answers from these very same processes.
        start_communication(targets, create_request, comm);

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
        while (all_locally_originated_receives_are_completed(process_answer,
                                                             comm) == false)
          maybe_answer_one_request(answer_request, comm);

        // 3) Signal to all other processes that all requests of this process
        //    have been answered
        signal_finish(comm);

        // 4) Nevertheless, this process has to keep on answering (potential)
        //    incoming requests until all processes have received the
        //    answer to all requests
        while (all_remotely_originated_receives_are_completed() == false)
          maybe_answer_one_request(answer_request, comm);

        // 5) process the answer to all requests
        clean_up_and_end_communication(comm);

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::start_communication(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &             create_request,
        const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
        // 1)
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
              AssertIndexRange(rank, Utilities::MPI::n_mpi_processes(comm));

              auto &send_buffer = send_buffers[index];
              send_buffer       = create_request(rank);

              // Post a request to send data
              auto ierr = MPI_Isend(send_buffer.data(),
                                    send_buffer.size() * sizeof(T1),
                                    MPI_BYTE,
                                    rank,
                                    tag_request,
                                    comm,
                                    &send_requests[index]);
              AssertThrowMPI(ierr);
            }

          // Also record that we expect an answer from each target we sent
          // a request to:
          n_outstanding_answers = n_targets;
        }
#else
        (void)targets;
        (void)create_request;
        (void)comm;
#endif
      }



      template <typename T1, typename T2>
      bool
      NBX<T1, T2>::all_locally_originated_receives_are_completed(
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
        // We know that all requests have come in when we have pending
        // messages from all targets with the right tag (some of which we may
        // have already taken care of below, after discovering their existence).
        // We can check for pending messages with MPI_IProbe, which returns
        // immediately with a return code that indicates whether
        // it has found a message from any process with a given
        // tag.
        if (n_outstanding_answers == 0)
          return true;
        else
          {
            const int tag_deliver = Utilities::MPI::internal::Tags::
              consensus_algorithm_nbx_process_deliver;

            int        request_is_pending;
            MPI_Status status;
            const auto ierr = MPI_Iprobe(
              MPI_ANY_SOURCE, tag_deliver, comm, &request_is_pending, &status);
            AssertThrowMPI(ierr);

            // If there is no pending message with this tag,
            // then we are clearly not done receiving everything
            // yet -- so return false.
            if (request_is_pending == 0)
              return false;
            else
              {
                // OK, so we have gotten a reply to our answer from
                // one rank. Let us process it, after double checking
                // that it is indeed one we were still expecting:
                const auto target = status.MPI_SOURCE;

                // Then query the size of the message, allocate enough memory,
                // receive the data, and process it.
                int message_size;
                {
                  const int ierr =
                    MPI_Get_count(&status, MPI_BYTE, &message_size);
                  AssertThrowMPI(ierr);
                }
                Assert(message_size % sizeof(T2) == 0, ExcInternalError());
                std::vector<T2> recv_buffer(message_size / sizeof(T2));

                {
                  const int tag_deliver = Utilities::MPI::internal::Tags::
                    consensus_algorithm_nbx_process_deliver;

                  const int ierr = MPI_Recv(recv_buffer.data(),
                                            recv_buffer.size() * sizeof(T2),
                                            MPI_BYTE,
                                            target,
                                            tag_deliver,
                                            comm,
                                            MPI_STATUS_IGNORE);
                  AssertThrowMPI(ierr);
                }

                process_answer(target, recv_buffer);

                // Finally, remove this rank from the list of outstanding
                // targets:
                --n_outstanding_answers;

                // We could do another go-around from the top of this
                // else-branch to see whether there are actually other messages
                // that are currently pending. But that would mean spending
                // substantial time in receiving answers while we should also be
                // sending answers to requests we have received from other
                // places. So let it be enough for now. If there are outstanding
                // answers, we will get back to this function before long and
                // can take care of them then.
                return (n_outstanding_answers == 0);
              }
          }

#else
        (void)process_answer;
        (void)comm;

        return true;
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::maybe_answer_one_request(
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &             answer_request,
        const MPI_Comm &comm)
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
        const auto ierr = MPI_Iprobe(
          MPI_ANY_SOURCE, tag_request, comm, &request_is_pending, &status);
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
                            comm,
                            MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            // Allocate memory for an answer message to the current request,
            // and ask the 'process' object to produce an answer:
            request_buffers.emplace_back(std::make_unique<std::vector<T2>>());
            auto &request_buffer = *request_buffers.back();
            request_buffer       = answer_request(other_rank, buffer_recv);

            // Then initiate sending the answer back to the requester.
            request_requests.emplace_back(std::make_unique<MPI_Request>());
            ierr = MPI_Isend(request_buffer.data(),
                             request_buffer.size() * sizeof(T2),
                             MPI_BYTE,
                             other_rank,
                             tag_deliver,
                             comm,
                             request_requests.back().get());
            AssertThrowMPI(ierr);
          }
#else
        (void)answer_request;
        (void)comm;
#endif
      }



      template <typename T1, typename T2>
      void
      NBX<T1, T2>::signal_finish(const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
        const auto ierr = MPI_Ibarrier(comm, &barrier_request);
        AssertThrowMPI(ierr);
#  else
        AssertThrow(false,
                    ExcMessage(
                      "ConsensusAlgorithms::NBX uses MPI 3.0 features. "
                      "You should compile with at least MPI 3.0."));
#  endif
#else
        (void)comm;
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
      NBX<T1, T2>::clean_up_and_end_communication(const MPI_Comm &comm)
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
          ierr = MPI_Barrier(comm);
          AssertThrowMPI(ierr);
#  endif
        }
#else
        (void)comm;
#endif
      }



      template <typename T1, typename T2>
      PEX<T1, T2>::PEX(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      PEX<T1, T2>::run(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &create_request,
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &answer_request,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
        Assert(has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, comm);

        // 1) Send requests and start receiving the answers.
        //    In particular, determine how many requests we should expect
        //    on the current process.
        const unsigned int n_requests =
          start_communication(targets, create_request, comm);

        // 2) Answer requests:
        for (unsigned int request = 0; request < n_requests; ++request)
          answer_one_request(request, answer_request, comm);

        // 3) Process answers:
        process_incoming_answers(targets.size(), process_answer, comm);

        // 4) Make sure all sends have successfully terminated:
        clean_up_and_end_communication();

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename T1, typename T2>
      unsigned int
      PEX<T1, T2>::start_communication(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &             create_request,
        const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;

        // 1) determine with which processes this process wants to communicate
        // with
        const unsigned int n_targets = targets.size();

        // 2) determine who wants to communicate with this process
        const unsigned int n_sources =
          compute_n_point_to_point_communications(comm, targets);

        // 2) allocate memory
        recv_buffers.resize(n_targets);
        send_buffers.resize(n_targets);
        send_request_requests.resize(n_targets);

        send_answer_requests.resize(n_sources);
        requests_buffers.resize(n_sources);

        // 4) send and receive
        for (unsigned int i = 0; i < n_targets; ++i)
          {
            const unsigned int rank = targets[i];
            AssertIndexRange(rank, Utilities::MPI::n_mpi_processes(comm));

            // pack data which should be sent
            auto &send_buffer = send_buffers[i];
            send_buffer       = create_request(rank);

            // start to send data
            auto ierr = MPI_Isend(send_buffer.data(),
                                  send_buffer.size() * sizeof(T1),
                                  MPI_BYTE,
                                  rank,
                                  tag_request,
                                  comm,
                                  &send_request_requests[i]);
            AssertThrowMPI(ierr);
          }

        return n_sources;
#else
        (void)targets;
        (void)create_request;
        (void)comm;
        return 0;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::answer_one_request(
        const unsigned int index,
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &             answer_request,
        const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_request = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_answer_request;
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;

        // Wait until we have a message ready for retrieval, though we don't
        // care which process it is from.
        MPI_Status status;
        int        ierr = MPI_Probe(MPI_ANY_SOURCE, tag_request, comm, &status);
        AssertThrowMPI(ierr);

        // Get rank of incoming message and verify that it makes sense
        const unsigned int other_rank = status.MPI_SOURCE;

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
                        comm,
                        &status);
        AssertThrowMPI(ierr);

        // Process request by asking the user-provided function for
        // the answer and post a send for it.
        auto &request_buffer = requests_buffers[index];
        request_buffer       = answer_request(other_rank, buffer_recv);

        ierr = MPI_Isend(request_buffer.data(),
                         request_buffer.size() * sizeof(T2),
                         MPI_BYTE,
                         other_rank,
                         tag_deliver,
                         comm,
                         &send_answer_requests[index]);
        AssertThrowMPI(ierr);
#else
        (void)answer_request;
        (void)comm;
        (void)index;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::process_incoming_answers(
        const unsigned int n_targets,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
#ifdef DEAL_II_WITH_MPI
        const int tag_deliver = Utilities::MPI::internal::Tags::
          consensus_algorithm_pex_process_deliver;

        // We know how many targets we have sent requests to. These
        // targets will all eventually send us their responses, but
        // we need not process them in order -- rather, just see what
        // comes in and then look at message originators' ranks and
        // message sizes
        for (unsigned int i = 0; i < n_targets; ++i)
          {
            MPI_Status status;
            {
              const int ierr =
                MPI_Probe(MPI_ANY_SOURCE, tag_deliver, comm, &status);
              AssertThrowMPI(ierr);
            }

            const auto other_rank = status.MPI_SOURCE;
            int        message_size;
            {
              const int ierr = MPI_Get_count(&status, MPI_BYTE, &message_size);
              AssertThrowMPI(ierr);
            }
            Assert(message_size % sizeof(T2) == 0, ExcInternalError());
            std::vector<T2> recv_buffer(message_size / sizeof(T2));

            // Now actually receive the answer. Because the MPI_Probe
            // above blocks until we have a message, we know that the
            // following MPI_Recv call will immediately succeed.
            {
              const int ierr = MPI_Recv(recv_buffer.data(),
                                        recv_buffer.size() * sizeof(T2),
                                        MPI_BYTE,
                                        other_rank,
                                        tag_deliver,
                                        comm,
                                        MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr);
            }

            process_answer(other_rank, recv_buffer);
          }
#else
        (void)n_targets;
        (void)process_answer;
        (void)comm;
#endif
      }



      template <typename T1, typename T2>
      void
      PEX<T1, T2>::clean_up_and_end_communication()
      {
#ifdef DEAL_II_WITH_MPI
        // Finalize all MPI_Request objects for both the
        // send-request and receive-answer operations.
        if (send_request_requests.size() > 0)
          {
            const int ierr = MPI_Waitall(send_request_requests.size(),
                                         send_request_requests.data(),
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
#endif
      }



      template <typename T1, typename T2>
      Serial<T1, T2>::Serial(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Serial<T1, T2>::run(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &create_request,
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &answer_request,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
        (void)comm;
        Assert((Utilities::MPI::job_supports_mpi() == false) ||
                 (Utilities::MPI::n_mpi_processes(comm) == 1),
               ExcMessage("You shouldn't use the 'Serial' class on "
                          "communicators that have more than one process "
                          "associated with it."));

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
            const std::vector<T1> request = create_request(0);
            const std::vector<T2> answer  = answer_request(0, request);
            process_answer(0, answer);
          }

        return targets; // nothing to do
      }



      template <typename T1, typename T2>
      Selector<T1, T2>::Selector(Process<T1, T2> &process, const MPI_Comm &comm)
        : Interface<T1, T2>(process, comm)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      Selector<T1, T2>::run(
        const std::vector<unsigned int> &targets,
        const std::function<std::vector<T1>(const unsigned int)>
          &create_request,
        const std::function<std::vector<T2>(const unsigned int,
                                            const std::vector<T1> &)>
          &answer_request,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &             process_answer,
        const MPI_Comm &comm)
      {
        // Depending on the number of processes we switch between
        // implementations. We reduce the threshold for debug mode to be
        // able to test also the non-blocking implementation. This feature
        // is tested by:
        // tests/multigrid/transfer_matrix_free_06.with_mpi=true.with_p4est=true.with_trilinos=true.mpirun=10.output

        const unsigned int n_procs = (Utilities::MPI::job_supports_mpi() ?
                                        Utilities::MPI::n_mpi_processes(comm) :
                                        1);
#ifdef DEAL_II_WITH_MPI
#  if DEAL_II_MPI_VERSION_GTE(3, 0)
#    ifdef DEBUG
        if (n_procs > 10)
#    else
        if (n_procs > 99)
#    endif
          consensus_algo.reset(new NBX<T1, T2>());
        else
#  endif
#endif
          if (n_procs > 1)
          consensus_algo.reset(new PEX<T1, T2>());
        else
          consensus_algo.reset(new Serial<T1, T2>());

        return consensus_algo->run(
          targets, create_request, answer_request, process_answer, comm);
      }


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
