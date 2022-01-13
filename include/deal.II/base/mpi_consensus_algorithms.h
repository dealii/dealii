// ---------------------------------------------------------------------
//
// Copyright (C) 2020 - 2021 by the deal.II authors
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

#ifndef dealii_mpi_consensus_algorithm_h
#define dealii_mpi_consensus_algorithm_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    /**
     * A namespace for consensus algorithms designed for dynamic-sparse
     * communication patterns.
     *
     * @ingroup MPI
     */
    namespace ConsensusAlgorithms
    {
      /**
       * A base class for concrete implementations of classes that
       * provide the information that the algorithms derived from
       * the ConsensusAlgorithms::Interface base class require. The main
       * functionality of this class is to return a list of process
       * ranks this process wants data from and to deal with the optional
       * payload of the messages sent/received by the ConsensusAlgorithm
       * classes.
       *
       * There are two kinds of messages:
       * - send/request message: A message consisting of a data request
       *   which should be answered by another process. This message is
       *   considered as a request message by the receiving rank.
       * - receive message: The answer to a send/request message.
       *
       * @tparam T1 The type of the elements of the vector to sent.
       * @tparam T2 The type of the elements of the vector to received.
       *
       * @note Since the payloads of the messages are optional, users have
       *    to deal with buffers themselves. The ConsensusAlgorithm classes
       *    (1) deliver only references to empty vectors (of size 0) the data
       *    to be sent can be inserted to or read from, and (2) communicate
       *    these vectors blindly.
       */
      template <typename T1, typename T2>
      class Process
      {
      public:
        /**
         * Destructor. Made `virtual` to ensure that one can work with
         * derived classes.
         */
        virtual ~Process() = default;

        /**
         * @return A vector of ranks this process wants to send a request to.
         *
         * @note This is the only method which has to be implemented since the
         *       payloads of the messages are optional.
         */
        virtual std::vector<unsigned int>
        compute_targets() = 0;

        /**
         * Add a payload to the request to the process with the specified rank.
         *
         * @param[in]  other_rank Rank of the process.
         * @param[out] send_buffer data to be sent part of the request
         * (optional).
         *
         * @note The buffer is empty. Before using it, you have to set its size.
         */
        virtual void
        create_request(const unsigned int other_rank,
                       std::vector<T1> &  send_buffer);

        /**
         * Prepare the buffer where the payload of the answer of the request to
         * the process with the specified rank is saved in. The most obvious
         * task is to resize the buffer, since it is empty when the function is
         * called.
         *
         * @param[in]  other_rank Rank of the process.
         * @param[out] recv_buffer Data to be sent as part of the request
         * (optional).
         */
        virtual void
        prepare_buffer_for_answer(const unsigned int other_rank,
                                  std::vector<T2> &  recv_buffer);

        /**
         * Prepare the buffer where the payload of the answer of the request to
         * the process with the specified rank is saved in.
         *
         * @param[in]  other_rank Rank of the process.
         * @param[in]  buffer_recv Received payload (optional).
         * @param[out] request_buffer Payload to be sent as part of the request
         *             (optional).
         *
         * @note The request_buffer is empty. Before using it, you have to set
         *       its size.
         */
        virtual void
        answer_request(const unsigned int     other_rank,
                       const std::vector<T1> &buffer_recv,
                       std::vector<T2> &      request_buffer);

        /**
         * Process the payload of the answer of the request to the process with
         * the specified rank.
         *
         * @param[in] other_rank rank of the process
         * @param[in] recv_buffer data to be sent part of the request (optional)
         */
        virtual void
        read_answer(const unsigned int     other_rank,
                    const std::vector<T2> &recv_buffer);
      };



      /**
       * A base class for algorithms that implement the task of coming up with
       * communication patterns to retrieve data from other processes in a
       * dynamic-sparse way. In computer science, this is often called a
       * <a href="https://en.wikipedia.org/wiki/Consensus_algorithm">consensus
       * problem</a>.
       *
       * The problem consensus algorithms are trying to solve is this: Let's
       * say you have $P$ processes that work together via MPI. Each (or at
       * least some) of these want to send information to some of the other
       * processes, or request information from other processes. No process
       * knows which other process wants to communicate with them. The challenge
       * is to determine who needs to talk to whom and what information needs to
       * be sent, and to come up with an algorithm that ensures that this
       * communication happens.
       *
       * That this is not a trivial problem can be seen by an analogy of the
       * postal service. There, some senders may request information from some
       * other participants in the postal service. So they send a letter that
       * requests the information, but the recipients do not know how many such
       * letters they need to expect (or that they should expect any at all).
       * They also do not know how long they need to keep checking their mailbox
       * for incoming requests. The recipients can be considered reliable,
       * however: We can assume that everyone who is sent a request puts a
       * letter with the answer in the mail. This time at least the recipients
       * of these answers know that they are waiting for these answers because
       * they have previously sent a request. They do not know in advance,
       * however, when the answer will arrive and how long to wait. The goal of
       * a consensus algorithm is then to come up with a strategy in which every
       * participant can say who they want to send requests to, what that
       * request is, and is then guaranteed an answer. The algorithm will only
       * return when all requests by all participants have been answered and the
       * answer delivered to the requesters.
       *
       * The problem is generally posed in terms of *requests* and *answers*.
       * In practice, either of these two may be empty messages. For example,
       * processes may simply want to send information to others that they know
       * these others need; in this case, the "answer" message may be empty
       * and its meaning is simply an affirmation that the information was
       * received. Similarly, in some cases processes simply need to inform
       * others that they want information, but the destination process knows
       * what information is being requested (based on where in the program
       * the request happens) and can send that information without there be
       * any identifying information in the request; in that case, the
       * request message may be empty and simply serve to identify the
       * requester. (Each message can be queried for its sender.)
       *
       * As mentioned in the first paragraph, the algorithms we are interested
       * in are what are call "dynamic-sparse". Dynamic-sparse in this context
       * means the following:
       * - Dynamic: By the time tha algorithm is called, the other processes do
       *   not know yet that they have to answer requests.
       * - Sparse: Each process only has to communicate with a small subset of
       *   processes of the MPI communicator.
       *
       * In order to run the communication algorithms, users of this class have
       * to provide a number of pieces of information:
       * - An MPI communicator.
       * - On each process, a list of ranks of processes to communicate with.
       * - Functionality to pack/unpack data to send as either the original
       *   request or as part of the answer.
       * This information is encoded through the `virtual` functions of classes
       * derived from ConsensusAlgorithm::Process, and the constructor of the
       * current class receives an object of a type derived from that class as
       * an argument.
       *
       * This base class only introduces a basic interface to achieve
       * these goals, while derived classes implement different algorithms
       * to actually compute such communication patterns and perform the
       * communication.
       *
       * @tparam T1 The type of the elements of the vector to be sent.
       * @tparam T2 The type of the elements of the vector to be received.
       */
      template <typename T1, typename T2>
      class Interface
      {
      public:
        /**
         * Constructor. @p process is an object that provides information
         * about what processes the current process wants to communicate with,
         * and the data to be sent/received. @p comm is the communicator on
         * which this communication is to happen.
         */
        Interface(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * Destructor. Made `virtual` to ensure that one can work with
         * derived classes.
         */
        virtual ~Interface() = default;

        /**
         * Run the consensus algorithm and return the requesting processes.
         */
        virtual std::vector<unsigned int>
        run() = 0;

      protected:
        /**
         * Reference to the process provided by the user.
         */
        Process<T1, T2> &process;

        /**
         * MPI communicator.
         */
        const MPI_Comm &comm;

        /**
         * Cache if job supports MPI.
         */
        const bool job_supports_mpi;

        /**
         * Rank of this process.
         */
        const unsigned int my_rank;

        /**
         * Number of processes in the communicator.
         */
        const unsigned int n_procs;
      };


      /**
       * This class implements a concrete algorithm for the
       * ConsensusAlgorithms::Interface base class, using only point-to-point
       * communications and a single IBarrier.
       *
       * @note This class closely follows @cite hoefler2010scalable. Since the
       *       algorithm shown there is not considering payloads, the algorithm
       *       has been modified here in such a way that synchronous sends
       *       (Issend) have been replaced by equivalent Isend/Irecv, where
       *       Irecv receives the answer to a request (with payload).
       *
       * @tparam T1 The type of the elements of the vector to be sent.
       * @tparam T2 The type of the elements of the vector to be received.
       */
      template <typename T1, typename T2>
      class NBX : public Interface<T1, T2>
      {
      public:
        /**
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator
         */
        NBX(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * Destructor.
         */
        virtual ~NBX() = default;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * List of processes this process wants to send requests to.
         */
        std::vector<unsigned int> targets;

        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<T1>> send_buffers;

        /**
         * Requests for sending requests.
         */
        std::vector<MPI_Request> send_requests;

        /**
         * Buffers for receiving answers to requests.
         */
        std::vector<std::vector<T2>> recv_buffers;


        /**
         * Requests for receiving answers to requests.
         */
        std::vector<MPI_Request> recv_requests;

        /**
         * Buffers for sending answers to requests.
         */
        std::vector<std::unique_ptr<std::vector<T2>>> request_buffers;

        /**
         * Requests for sending answers to requests.
         */
        std::vector<std::unique_ptr<MPI_Request>> request_requests;

        // request for barrier
        MPI_Request barrier_request;
#endif

        /**
         * List of processes who have made a request to this process.
         */
        std::set<unsigned int> requesting_processes;

        /**
         * Check whether all of the requests for answers that were created by
         * the communication posted from the current process to other ranks
         * have been satisfied.
         */
        bool
        all_locally_originated_receives_are_completed();

        /**
         * Signal to all other ranks that this rank has received all request
         * answers via entering IBarrier.
         */
        void
        signal_finish();

        /**
         * Check whether all of the requests for answers that were created by
         * communication posted from other processes to the current rank
         * have been satisfied.
         */
        bool
        all_remotely_originated_receives_are_completed();

        /**
         * Check whether a request message from another rank has been received,
         * and if so, process the request by storing the data and sending an
         * answer.
         */
        void
        maybe_answer_one_request();

        /**
         * Start to send all requests via ISend and post IRecvs for the incoming
         * answer messages.
         */
        void
        start_communication();

        /**
         * After all rank has received all answers, the MPI data structures can
         * be freed and the received answers can be processed.
         */
        void
        clean_up_and_end_communication();
      };

      /**
       * This class implements a concrete algorithm for the
       * ConsensusAlgorithms::Interface base class, using a two step approach.
       * In the first step the source ranks are determined and in the second
       * step a static sparse data exchange is performed.
       *
       * @note In contrast to NBX, this class splits the same
       *   task into two distinct steps. In the first step, all processes
       *   are identified who want to send a request to this process. In the
       *   second step, the data is exchanged. However, since - in the
       *   second step - now it is clear how many requests have to be answered,
       *   i.e. when this process can stop waiting for requests, no IBarrier is
       *   needed.
       *
       * @note Under the hood, this function uses
       *   Utilities::MPI::compute_point_to_point_communication_pattern()
       *   to determine the source processes, which itself is based on the
       *   NBX-algorithm from @cite hoefler2010scalable that is implemented
       *   in the ConsensusAlgorithms::NBX class (a sister class to the
       *   current one).
       *
       * @tparam T1 The type of the elements of the vector to be sent.
       * @tparam T2 The type of the elements of the vector to be received.
       */
      template <typename T1, typename T2>
      class PEX : public Interface<T1, T2>
      {
      public:
        /**
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator
         */
        PEX(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * Destructor.
         */
        virtual ~PEX() = default;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * List of ranks of processes this processes wants to send a request to.
         */
        std::vector<unsigned int> targets;

        /**
         * List of ranks of processes wanting to send a request to this process.
         */
        std::vector<unsigned int> sources;

        // data structures to send and receive requests

        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<T1>> send_buffers;

        /**
         * Buffers for receiving answers to requests.
         */
        std::vector<std::vector<T2>> recv_buffers;

        /**
         * Requests for sending requests and receiving answers to requests.
         */
        std::vector<MPI_Request> send_and_recv_buffers;

        /**
         * Buffers for sending answers to requests.
         */
        std::vector<std::vector<T2>> requests_buffers;

        /**
         * Requests for sending answers to requests.
         */
        std::vector<MPI_Request> requests_answers;
#endif
        /**
         * List of processes who have made a request to this process.
         */
        std::set<unsigned int> requesting_processes;

        /**
         * Start to send all requests via ISend and post IRecvs for the incoming
         * answer messages.
         */
        unsigned int
        start_communication();

        /**
         * The ith request message from another rank has been received: process
         * the request and send an answer.
         */
        void
        answer_requests(int index);

        /**
         * After all answers have been exchanged, the MPI data structures can be
         * freed and the received answers can be processed.
         */
        void
        clean_up_and_end_communication();
      };



      /**
       * A serial fall back for the above classes to allow programming
       * independently of whether MPI is used or not.
       */
      template <typename T1, typename T2>
      class Serial : public Interface<T1, T2>
      {
      public:
        /**
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator (ignored)
         */
        Serial(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run() override;
      };

      /**
       * A class which delegates its task to other
       * ConsensusAlgorithms::Interface implementations depending on the number
       * of processes in the MPI communicator. For a small number of processes
       * it uses PEX and for a large number of processes NBX. The threshold
       * depends if the program is compiled in debug or release mode.
       *
       * @tparam T1 The type of the elements of the vector to be sent.
       * @tparam T2 The type of the elements of the vector to be received.
       */
      template <typename T1, typename T2>
      class Selector : public Interface<T1, T2>
      {
      public:
        /**
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator.
         */
        Selector(Process<T1, T2> &process, const MPI_Comm &comm);

        /**
         * Destructor.
         */
        virtual ~Selector() = default;

        /**
         * @copydoc Interface::run()
         *
         * @note The function call is delegated to another ConsensusAlgorithms::Interface implementation.
         */
        virtual std::vector<unsigned int>
        run() override;

      private:
        // Pointer to the actual ConsensusAlgorithms::Interface implementation.
        std::shared_ptr<Interface<T1, T2>> consensus_algo;
      };

      /**
       * This class implements Utilities::MPI::ConsensusAlgorithms::Process,
       * using user-provided function wrappers.
       * The advantage of this class is that users do not have to write their
       * own implementation but can register lambda functions directly.
       */
      template <typename T1, typename T2>
      class AnonymousProcess : public Process<T1, T2>
      {
      public:
        /**
         * Register functions that should be called for implementing the
         * interface of Process.
         *
         * @param function_compute_targets called during `compute_targets`.
         * @param function_create_request called during `create_request`.
         * @param function_answer_request called during `answer_request`.
         * @param function_prepare_buffer_for_answer called during
         *   `prepare_buffer_for_answer`.
         * @param function_read_answer called during `read_answer`.
         */
        AnonymousProcess(
          const std::function<std::vector<unsigned int>()>
            &function_compute_targets,
          const std::function<void(const unsigned int, std::vector<T1> &)>
            &function_create_request = {},
          const std::function<void(const unsigned int,
                                   const std::vector<T1> &,
                                   std::vector<T2> &)>
            &function_answer_request = {},
          const std::function<void(const unsigned int, std::vector<T2> &)>
            &function_prepare_buffer_for_answer = {},
          const std::function<void(const unsigned int, const std::vector<T2> &)>
            &function_read_answer = {});

        /**
         * @copydoc Process::compute_targets()
         */
        std::vector<unsigned int>
        compute_targets() override;

        /**
         * @copydoc Process::create_request()
         */
        void
        create_request(const unsigned int other_rank,
                       std::vector<T1> &  send_buffer) override;

        /**
         * @copydoc Process::answer_request()
         */
        void
        answer_request(const unsigned int     other_rank,
                       const std::vector<T1> &buffer_recv,
                       std::vector<T2> &      request_buffer) override;

        /**
         * @copydoc Process::prepare_buffer_for_answer()
         */
        void
        prepare_buffer_for_answer(const unsigned int other_rank,
                                  std::vector<T2> &  recv_buffer) override;

        /**
         * @copydoc Process::read_answer()
         */
        void
        read_answer(const unsigned int     other_rank,
                    const std::vector<T2> &recv_buffer) override;

      private:
        const std::function<std::vector<unsigned int>()>
          function_compute_targets;
        const std::function<void(const int, std::vector<T1> &)>
          function_create_request;
        const std::function<
          void(const unsigned int, const std::vector<T1> &, std::vector<T2> &)>
          function_answer_request;
        const std::function<void(const int, std::vector<T2> &)>
          function_prepare_buffer_for_answer;
        const std::function<void(const int, const std::vector<T2> &)>
          function_read_answer;
      };



      template <typename T1, typename T2>
      AnonymousProcess<T1, T2>::AnonymousProcess(
        const std::function<std::vector<unsigned int>()>
          &function_compute_targets,
        const std::function<void(const unsigned int, std::vector<T1> &)>
          &                                           function_create_request,
        const std::function<void(const unsigned int,
                                 const std::vector<T1> &,
                                 std::vector<T2> &)> &function_answer_request,
        const std::function<void(const unsigned int, std::vector<T2> &)>
          &function_prepare_buffer_for_answer,
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &function_read_answer)
        : function_compute_targets(function_compute_targets)
        , function_create_request(function_create_request)
        , function_answer_request(function_answer_request)
        , function_prepare_buffer_for_answer(function_prepare_buffer_for_answer)
        , function_read_answer(function_read_answer)
      {}



      template <typename T1, typename T2>
      std::vector<unsigned int>
      AnonymousProcess<T1, T2>::compute_targets()
      {
        return function_compute_targets();
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::create_request(const unsigned int other_rank,
                                               std::vector<T1> &  send_buffer)
      {
        if (function_create_request)
          function_create_request(other_rank, send_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::answer_request(
        const unsigned int     other_rank,
        const std::vector<T1> &buffer_recv,
        std::vector<T2> &      request_buffer)
      {
        if (function_answer_request)
          function_answer_request(other_rank, buffer_recv, request_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::prepare_buffer_for_answer(
        const unsigned int other_rank,
        std::vector<T2> &  recv_buffer)
      {
        if (function_prepare_buffer_for_answer)
          function_prepare_buffer_for_answer(other_rank, recv_buffer);
      }



      template <typename T1, typename T2>
      void
      AnonymousProcess<T1, T2>::read_answer(const unsigned int     other_rank,
                                            const std::vector<T2> &recv_buffer)
      {
        if (function_read_answer)
          function_read_answer(other_rank, recv_buffer);
      }



    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
