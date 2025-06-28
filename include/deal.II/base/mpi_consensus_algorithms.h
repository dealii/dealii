// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_mpi_consensus_algorithm_h
#define dealii_mpi_consensus_algorithm_h

#include <deal.II/base/config.h>

#include <deal.II/base/mpi.h>
#include <deal.II/base/mpi.templates.h>
#include <deal.II/base/mpi_tags.h>

#include <set>
#include <vector>

DEAL_II_NAMESPACE_OPEN


namespace Utilities
{
  namespace MPI
  {
    /**
     * A namespace for algorithms that implement the task of communicating
     * in a dynamic-sparse way. In computer science, this is often called a
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
     * in are "dynamic-sparse":
     * - Dynamic: By the time the algorithm is called, the other processes do
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
     * This information is typically either provided as direct objects (for
     * the first two of the points above), or as function objects (for the
     * third point above). In the latter case, the function objects are often
     * simply lambda functions declared right in the context where one wants
     * to run a consensus algorithm; these lambda functions may then reference
     * variables that are active at the point of declaration of the lambda
     * function, such as variables local to the surrounding function.
     *
     *
     * <h3>Available implementations</h3>
     *
     * There are many ways to implement the general functionality required
     * for these "consensus algorithms". This namespace provides several
     * implementations of consensus algorithms, specifically the
     * NBX and PEX algorithms, along with a serial one for the case where
     * one wants to run such an algorithm on a single process. The key
     * entry points to these algorithms are the
     * nbx(), pex(), serial(), and selector() functions that take a
     * communicator, a list of targets, and a number of functions
     * as argument. The selector() function redirects to the other
     * implementations based on the number of processes that participate
     * in an MPI universe, since some implementations are better or worse
     * suited for large or small parallel computations.
     *
     * This namespace also implements specializations of each of the
     * functions for the specific case where a calling process is not
     * actually interested in receiving and processing answers -- that
     * is, the goal is simply to *send* messages to a number of targets,
     * but no answer is required; all we want to know is that by the end of
     * the call, all targets have been sent their respective data. This,
     * strictly speaking, does not fall under the umbrella of "consensus
     * algorithms", but is really just a "some-to-some" communication.
     * (This operation is also provided by the Utilities::MPI::some_to_some()
     * function, though with a different interface.) For this special
     * case, the functions in this namespace only need to receive
     * an MPI communicator, a list of targets, and function objects that
     * encode and decode the messages to be sent, but no functions for
     * encoding a reply, or processing a reply.
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
       * @tparam RequestType The type of the elements of the vector to sent.
       * @tparam AnswerType The type of the elements of the vector to received.
       *
       * @note Since the payloads of the messages are optional, users have
       *    to deal with buffers themselves. The ConsensusAlgorithm classes
       *    (1) deliver only references to empty vectors (of size 0) the data
       *    to be sent can be inserted to or read from, and (2) communicate
       *    these vectors blindly.
       *
       * @deprecated Instead of deriving a class from this base class and
       *   providing a corresponding object to one of the run() functions,
       *   use the free functions in this namespace that take function
       *   objects as arguments.
       */
      template <typename RequestType, typename AnswerType>
      class DEAL_II_DEPRECATED_EARLY Process
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
        create_request(const unsigned int other_rank, RequestType &send_buffer);

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
        answer_request(const unsigned int other_rank,
                       const RequestType &buffer_recv,
                       AnswerType        &request_buffer);

        /**
         * Process the payload of the answer of the request to the process with
         * the specified rank.
         *
         * @param[in] other_rank rank of the process
         * @param[in] recv_buffer data to be sent part of the request (optional)
         */
        virtual void
        read_answer(const unsigned int other_rank,
                    const AnswerType  &recv_buffer);
      };



      /**
       * A base class for algorithms that implement consensus algorithms,
       * see the documentation of the surrounding namespace for more
       * information.
       *
       * This base class only introduces a basic interface to achieve
       * these goals, while derived classes implement different algorithms
       * to actually compute such communication patterns and perform the
       * communication.
       *
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      class Interface
      {
      public:
        /**
         * Default constructor.
         */
        Interface() = default;

        /**
         * Destructor. Made `virtual` to ensure that one can work with
         * derived classes.
         */
        virtual ~Interface() = default;

        /**
         * Run the consensus algorithm and return a vector of process ranks
         * that have requested answers from the current process.
         *
         * This version of the run() function simply unpacks the functions
         * packaged in `process` and calls the version of the run() function
         * that takes a number of `std::function` arguments.
         *
         * @deprecated Instead of deriving a class from the Process base class and
         *   providing a corresponding object to this function,
         *   use the other run() function in this class that takes function
         *   objects as arguments.
         */
        DEAL_II_DEPRECATED_EARLY
        std::vector<unsigned int>
        run(Process<RequestType, AnswerType> &process, const MPI_Comm comm);

        /**
         * Run the consensus algorithm and return a vector of process ranks
         * that have requested answers from the current process.
         *
         * @param[in] targets A vector that contains the ranks of processes
         *   to which requests should be sent and from which answers need
         *   to be received.
         * @param[in] create_request A function object that takes the rank
         *   of a target process as argument and returns the message that
         *   forms the request to this target.
         * @param[in] answer_request A function that takes as arguments the
         *   rank of the process that has sent a request to us, along with
         *   the message of the request, and returns the message that forms
         *   the answer that should be sent back to the requesting process.
         * @param[in] process_answer A function object that takes as argument
         *   the rank of a process from which we have received an answer
         *   to a previously sent request, along with the message that
         *   forms this answer. This function is used to describe what
         *   the caller of the consensus algorithm wants to do with the
         *   received answer.
         * @param[in] comm The MPI communicator on which the whole algorithm
         *   is to be performed.
         */
        virtual std::vector<unsigned int>
        run(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm) = 0;
      };


      /**
       * This class implements a concrete algorithm for the
       * ConsensusAlgorithms::Interface base class, using only point-to-point
       * communications and a single IBarrier. This algorithm is suitable
       * for very large process counts because it does not require the
       * allocation of arrays with size proportional to the number of processes.
       *
       * @note This class closely follows @cite hoefler2010scalable, but our
       *   implementation also deals with payloads.
       *
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      class NBX : public Interface<RequestType, AnswerType>
      {
      public:
        /**
         * Default constructor.
         */
        NBX() = default;

        /**
         * Destructor.
         */
        virtual ~NBX() = default;

        // Import the declarations from the base class.
        using Interface<RequestType, AnswerType>::run;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm) override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<char>> send_buffers;

        /**
         * Requests for sending requests.
         */
        std::vector<MPI_Request> send_requests;

        /**
         * Buffers for sending answers to requests. We use a vector of
         * pointers because that guarantees that the buffers themselves
         * are newer moved around in memory, even if the vector is
         * resized and consequently its elements (the pointers) are moved
         * around.
         */
        std::vector<std::unique_ptr<std::vector<char>>> request_buffers;

        /**
         * Requests for sending answers to requests.
         */
        std::vector<std::unique_ptr<MPI_Request>> request_requests;

        /**
         * The number of processes from which we are still expecting answers.
         */
        unsigned int n_outstanding_answers;

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
        all_locally_originated_receives_are_completed(
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm);

        /**
         * Signal to all other ranks that this rank has received all request
         * answers via entering IBarrier.
         */
        void
        signal_finish(const MPI_Comm comm);

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
        maybe_answer_one_request(
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const MPI_Comm                                        comm);

        /**
         * Start to send all requests via ISend and post IRecvs for the incoming
         * answer messages.
         */
        void
        start_communication(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const MPI_Comm                                        comm);

        /**
         * After all rank has received all answers, the MPI data structures can
         * be freed and the received answers can be processed.
         */
        void
        clean_up_and_end_communication(const MPI_Comm comm);
      };


      /**
       * This function implements a concrete algorithm for the
       * consensus algorithms problem (see the documentation of the
       * surrounding namespace), using only point-to-point
       * communications and a single IBarrier. This algorithm is suitable
       * for very large process counts because it does not require the
       * allocation of arrays with size proportional to the number of processes.
       *
       * @note This class closely follows @cite hoefler2010scalable, but our
       *   implementation also deals with payloads.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] answer_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and returns the message that forms
       *   the answer that should be sent back to the requesting process.
       * @param[in] process_answer A function object that takes as argument
       *   the rank of a process from which we have received an answer
       *   to a previously sent request, along with the message that
       *   forms this answer. This function is used to describe what
       *   the caller of the consensus algorithm wants to do with the
       *   received answer.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       * @tparam AnswerType The type of the object to be received.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()`,
       *   `answer_request()`, and `process_answer()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm);

      /**
       * This function provides a specialization of the one above for
       * the case where a sending process does not require an answer.
       * Strictly speaking, the name "request" is then incorrect, as it
       * is simply one process sending a message to another, but we keep
       * the name for symmetry with the function above that processes both
       * requests and answers.
       *
       * Since the function does not deal with answers, the algorithm
       * implemented is really just a "some-to-some algorithm", as
       * also provided by the Utilities::MPI::some_to_some() function.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] process_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and processes that message.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()` and
       *   `process_request()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<void(const unsigned int, const RequestType &)>
                        &process_request,
          const MPI_Comm comm);

      /**
       * This class implements a concrete algorithm for the
       * ConsensusAlgorithms::Interface base class, using a two step approach.
       * In the first step the source ranks are determined and in the second
       * step a static sparse data exchange is performed. This algorithm is most
       * suitable for relatively small process counts -- say, less than 100.
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
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      class PEX : public Interface<RequestType, AnswerType>
      {
      public:
        /**
         * Default constructor.
         */
        PEX() = default;

        /**
         * Destructor.
         */
        virtual ~PEX() = default;

        // Import the declarations from the base class.
        using Interface<RequestType, AnswerType>::run;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm) override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<char>> send_buffers;

        /**
         * Buffers for receiving answers to requests.
         */
        std::vector<std::vector<char>> recv_buffers;

        /**
         * MPI request objects for sending request messages.
         */
        std::vector<MPI_Request> send_request_requests;

        /**
         * Buffers for sending answers to requests.
         */
        std::vector<std::vector<char>> requests_buffers;

        /**
         * Requests for sending answers to requests.
         */
        std::vector<MPI_Request> send_answer_requests;
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
        start_communication(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const MPI_Comm                                        comm);

        /**
         * The `index`th request message from another rank has been received:
         * process the request and send an answer.
         */
        void
        answer_one_request(
          const unsigned int                                    index,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const MPI_Comm                                        comm);

        /**
         * Receive and process all of the incoming responses to the
         * requests we sent.
         */
        void
        process_incoming_answers(
          const unsigned int n_targets,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm);

        /**
         * After all answers have been exchanged, the MPI data structures can be
         * freed and the received answers can be processed.
         */
        void
        clean_up_and_end_communication();
      };



      /**
       * This function implements a concrete algorithm for the
       * consensus algorithms problem (see the documentation of the
       * surrounding namespace), using a two step approach.
       * In the first step the source ranks are determined and in the second
       * step a static sparse data exchange is performed. This algorithm is most
       * suitable for relatively small process counts -- say, less than 100.
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
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] answer_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and returns the message that forms
       *   the answer that should be sent back to the requesting process.
       * @param[in] process_answer A function object that takes as argument
       *   the rank of a process from which we have received an answer
       *   to a previously sent request, along with the message that
       *   forms this answer. This function is used to describe what
       *   the caller of the consensus algorithm wants to do with the
       *   received answer.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       * @tparam AnswerType The type of the object to be received.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()`,
       *   `answer_request()`, and `process_answer()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm);

      /**
       * This function provides a specialization of the one above for
       * the case where a sending process does not require an answer.
       * Strictly speaking, the name "request" is then incorrect, as it
       * is simply one process sending a message to another, but we keep
       * the name for symmetry with the function above that processes both
       * requests and answers.
       *
       * Since the function does not deal with answers, the algorithm
       * implemented is really just a "some-to-some algorithm", as
       * also provided by the Utilities::MPI::some_to_some() function.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] process_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and processes that message.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()` and
       *   `process_request()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<void(const unsigned int, const RequestType &)>
                        &process_request,
          const MPI_Comm comm);


      /**
       * A serial fall back for the above classes to allow programming
       * independently of whether MPI is used or not.
       */
      template <typename RequestType, typename AnswerType>
      class Serial : public Interface<RequestType, AnswerType>
      {
      public:
        /**
         * Default constructor.
         */
        Serial() = default;

        // Import the declarations from the base class.
        using Interface<RequestType, AnswerType>::run;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm) override;
      };



      /**
       * This function implements a concrete algorithm for the
       * consensus algorithms problem (see the documentation of the
       * surrounding namespace), as a fall-back option for the case
       * where the communicator provided has only one rank (or when
       * MPI is simply not used at all).
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] answer_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and returns the message that forms
       *   the answer that should be sent back to the requesting process.
       * @param[in] process_answer A function object that takes as argument
       *   the rank of a process from which we have received an answer
       *   to a previously sent request, along with the message that
       *   forms this answer. This function is used to describe what
       *   the caller of the consensus algorithm wants to do with the
       *   received answer.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed. Since this function is supposed to be run
       *   only for serial cases, the function throws an exception if
       *   the provided communicator denotes an MPI universe with more
       *   than one process.
       *
       * @tparam RequestType The type of the object to be sent.
       * @tparam AnswerType The type of the object to be received.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      serial(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm);

      /**
       * This function provides a specialization of the one above for
       * the case where a sending process does not require an answer.
       * Strictly speaking, the name "request" is then incorrect, as it
       * is simply one process sending a message to another, but we keep
       * the name for symmetry with the function above that processes both
       * requests and answers.
       *
       * Since the function does not deal with answers, the algorithm
       * implemented is really just a "some-to-some algorithm", as
       * also provided by the Utilities::MPI::some_to_some() function.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] process_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and processes that message.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed. Since this function is supposed to be run
       *   only for serial cases, the function throws an exception if
       *   the provided communicator denotes an MPI universe with more
       *   than one process.
       *
       * @tparam RequestType The type of the object to be sent.
       */
      template <typename RequestType>
      std::vector<unsigned int>
      serial(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<void(const unsigned int, const RequestType &)>
                      &process_request,
        const MPI_Comm comm);



      /**
       * A class which delegates its task to other
       * ConsensusAlgorithms::Interface implementations depending on the number
       * of processes in the MPI communicator. For a small number of processes
       * it uses PEX and for a large number of processes NBX. The threshold
       * depends if the program is compiled in debug or release mode, but the
       * goal is to always use the most efficient algorithm for however many
       * processes participate in the communication.
       *
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      class Selector : public Interface<RequestType, AnswerType>
      {
      public:
        /**
         * Default constructor.
         */
        Selector() = default;

        /**
         * Destructor.
         */
        virtual ~Selector() = default;

        // Import the declarations from the base class.
        using Interface<RequestType, AnswerType>::run;

        /**
         * @copydoc Interface::run()
         *
         * @note The function call is delegated to another ConsensusAlgorithms::Interface implementation.
         */
        virtual std::vector<unsigned int>
        run(
          const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm) override;

      private:
        // Pointer to the actual ConsensusAlgorithms::Interface implementation.
        std::shared_ptr<Interface<RequestType, AnswerType>> consensus_algo;
      };



      /**
       * This function implements a concrete algorithm for the
       * consensus algorithms problem (see the documentation of the
       * surrounding namespace). In particular, it delegates its work
       * to one of the other functions in this namespace depending on the number
       * of processes in the MPI communicator. For a small number of processes
       * it uses pex() and for a large number of processes nbx(). The threshold
       * depends if the program is compiled in debug or release mode, but the
       * goal is to always use the most efficient algorithm for however many
       * processes participate in the communication.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] answer_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and returns the message that forms
       *   the answer that should be sent back to the requesting process.
       * @param[in] process_answer A function object that takes as argument
       *   the rank of a process from which we have received an answer
       *   to a previously sent request, along with the message that
       *   forms this answer. This function is used to describe what
       *   the caller of the consensus algorithm wants to do with the
       *   received answer.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       * @tparam AnswerType The type of the object to be received.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()`,
       *   `answer_request()`, and `process_answer()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      selector(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm);

      /**
       * This function provides a specialization of the one above for
       * the case where a sending process does not require an answer.
       * Strictly speaking, the name "request" is then incorrect, as it
       * is simply one process sending a message to another, but we keep
       * the name for symmetry with the function above that processes both
       * requests and answers.
       *
       * Since the function does not deal with answers, the algorithm
       * implemented is really just a "some-to-some algorithm", as
       * also provided by the Utilities::MPI::some_to_some() function.
       *
       * @param[in] targets A vector that contains the ranks of processes
       *   to which requests should be sent and from which answers need
       *   to be received.
       * @param[in] create_request A function object that takes the rank
       *   of a target process as argument and returns the message that
       *   forms the request to this target.
       * @param[in] process_request A function that takes as arguments the
       *   rank of the process that has sent a request to us, along with
       *   the message of the request, and processes that message.
       * @param[in] comm The MPI communicator on which the whole algorithm
       *   is to be performed.
       *
       * @tparam RequestType The type of the object to be sent.
       *
       * @note Nothing good will generally happen if any of the function
       *   objects passed as arguments throws an exception when called.
       *   This is because when that happens, one MPI process will stop
       *   participating in MPI communications, deadlocking the other
       *   processes. As a consequence, the `create_request()`
       *   and `process_request()` functions should
       *   not throw any exceptions; if they encounter error
       *   conditions, they should instead call `MPI_Abort()` or use
       *   another way to reliably print an error message and then
       *   bring the MPI universe down.
       */
      template <typename RequestType>
      std::vector<unsigned int>
      selector(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<void(const unsigned int, const RequestType &)>
                      &process_request,
        const MPI_Comm comm);



#ifndef DOXYGEN
      // Implementation of the functions in this namespace.

      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm)
      {
        return NBX<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<void(const unsigned int, const RequestType &)>
                        &process_request,
          const MPI_Comm comm)
      {
        // TODO: For the moment, simply implement this special case by
        // forwarding to the other function with rewritten function
        // objects and using an empty type as answer type. This way,
        // we have the interface in place and can provide a more
        // efficient implementation later on.
        using EmptyType = std::tuple<>;

        return nbx<RequestType, EmptyType>(
          targets,
          create_request,
          // answer_request:
          [&process_request](const unsigned int source_rank,
                             const RequestType &request) -> EmptyType {
            process_request(source_rank, request);
            // Return something. What it is is arbitrary here, except that
            // we want it to be as small an object as possible. Using
            // std::tuple<> is interpreted as an empty object that is packed
            // down to a zero-length char array.
            return {};
          },
          // process_answer:
          [](const unsigned int /*target_rank */,
             const EmptyType & /*answer*/) {},
          comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<AnswerType(const unsigned int,
                                         const RequestType &)> &answer_request,
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm)
      {
        return PEX<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int>                      &targets,
          const std::function<RequestType(const unsigned int)> &create_request,
          const std::function<void(const unsigned int, const RequestType &)>
                        &process_request,
          const MPI_Comm comm)
      {
        // TODO: For the moment, simply implement this special case by
        // forwarding to the other function with rewritten function
        // objects and using an empty type as answer type. This way,
        // we have the interface in place and can provide a more
        // efficient implementation later on.
        using EmptyType = std::tuple<>;

        return pex<RequestType, EmptyType>(
          targets,
          create_request,
          // answer_request:
          [&process_request](const unsigned int source_rank,
                             const RequestType &request) -> EmptyType {
            process_request(source_rank, request);
            // Return something. What it is is arbitrary here, except that
            // we want it to be as small an object as possible. Using
            // std::tuple<> is interpreted as an empty object that is packed
            // down to a zero-length char array.
            return {};
          },
          // process_answer:
          [](const unsigned int /*target_rank */,
             const EmptyType & /*answer*/) {},
          comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      serial(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        return Serial<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType>
      std::vector<unsigned int>
      serial(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<void(const unsigned int, const RequestType &)>
                      &process_request,
        const MPI_Comm comm)
      {
        // TODO: For the moment, simply implement this special case by
        // forwarding to the other function with rewritten function
        // objects and using an empty type as answer type. This way,
        // we have the interface in place and can provide a more
        // efficient implementation later on.
        using EmptyType = std::tuple<>;

        return serial<RequestType, EmptyType>(
          targets,
          create_request,
          // answer_request:
          [&process_request](const unsigned int source_rank,
                             const RequestType &request) -> EmptyType {
            process_request(source_rank, request);
            // Return something. What it is is arbitrary here, except that
            // we want it to be as small an object as possible. Using
            // std::tuple<> is interpreted as an empty object that is packed
            // down to a zero-length char array.
            return {};
          },
          // process_answer:
          [](const unsigned int /*target_rank */,
             const EmptyType & /*answer*/) {},
          comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      selector(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        return Selector<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType>
      std::vector<unsigned int>
      selector(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<void(const unsigned int, const RequestType &)>
                      &process_request,
        const MPI_Comm comm)
      {
        // TODO: For the moment, simply implement this special case by
        // forwarding to the other function with rewritten function
        // objects and using an empty type as answer type. This way,
        // we have the interface in place and can provide a more
        // efficient implementation later on.
        using EmptyType = std::tuple<>;

        return selector<RequestType, EmptyType>(
          targets,
          create_request,
          // answer_request:
          [&process_request](const unsigned int source_rank,
                             const RequestType &request) -> EmptyType {
            process_request(source_rank, request);
            // Return something. What it is is arbitrary here, except that
            // we want it to be as small an object as possible. Using
            // std::tuple<> is interpreted as an empty object that is packed
            // down to a zero-length char array.
            return {};
          },
          // process_answer:
          [](const unsigned int /*target_rank */,
             const EmptyType & /*answer*/) {},
          comm);
      }

#endif


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities



#ifndef DOXYGEN

// ----------------- Implementation of template functions

namespace Utilities
{
  namespace MPI
  {
    namespace ConsensusAlgorithms
    {
      namespace internal
      {
        /**
         * Return whether a vector of targets (MPI ranks) has only unique
         * elements.
         */
        inline bool
        has_unique_elements(const std::vector<unsigned int> &targets)
        {
          std::vector<unsigned int> my_destinations = targets;
          std::sort(my_destinations.begin(), my_destinations.end());
          return (std::adjacent_find(my_destinations.begin(),
                                     my_destinations.end()) ==
                  my_destinations.end());
        }



        /**
         * Handle exceptions inside the ConsensusAlgorithm::run() functions.
         */
        inline void
        handle_exception(std::exception_ptr &&exception, const MPI_Comm comm)
        {
#  ifdef DEAL_II_WITH_MPI
          // an exception within a ConsensusAlgorithm likely causes an
          // MPI deadlock. Abort with a reasonable error message instead.
          try
            {
              std::rethrow_exception(exception);
            }
          catch (ExceptionBase &exc)
            {
              // report name of the deal.II exception:
              std::cerr
                << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
              std::cerr
                << "Exception '" << exc.get_exc_name() << "'"
                << " on rank " << Utilities::MPI::this_mpi_process(comm)
                << " on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

              // Then bring down the whole MPI world
              MPI_Abort(comm, 255);
            }
          catch (std::exception &exc)
            {
              std::cerr
                << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
              std::cerr
                << "Exception within ConsensusAlgorithm"
                << " on rank " << Utilities::MPI::this_mpi_process(comm)
                << " on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

              // Then bring down the whole MPI world
              MPI_Abort(comm, 255);
            }
          catch (...)
            {
              std::cerr
                << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
              std::cerr
                << "Unknown exception within ConsensusAlgorithm!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

              // Then bring down the whole MPI world
              MPI_Abort(comm, 255);
            }
#  else
          (void)comm;

          // No need to be concerned about deadlocks without MPI.
          // Defer to exception handling further up the callstack.
          std::rethrow_exception(exception);
#  endif
        }
      } // namespace internal



      template <typename RequestType, typename AnswerType>
      void
      Process<RequestType, AnswerType>::answer_request(const unsigned int,
                                                       const RequestType &,
                                                       AnswerType &)
      {
        // nothing to do
      }



      template <typename RequestType, typename AnswerType>
      void
      Process<RequestType, AnswerType>::create_request(const unsigned int,
                                                       RequestType &)
      {
        // nothing to do
      }



      template <typename RequestType, typename AnswerType>
      void
      Process<RequestType, AnswerType>::read_answer(const unsigned int,
                                                    const AnswerType &)
      {
        // nothing to do
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      Interface<RequestType, AnswerType>::run(
        Process<RequestType, AnswerType> &process,
        const MPI_Comm                    comm)
      {
        // Unpack the 'process' object and call the function that takes
        // function objects for all operations.
        return run(
          process.compute_targets(),
          /* create_request: */
          [&process](const unsigned int target) {
            RequestType request;
            process.create_request(target, request);
            return request;
          },
          /* answer_request: */
          [&process](const unsigned int source, const RequestType &request) {
            AnswerType answer;
            process.answer_request(source, request, answer);
            return answer;
          },
          /* process_answer: */
          [&process](const unsigned int target, const AnswerType &answer) {
            process.read_answer(target, answer);
          },
          comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      NBX<RequestType, AnswerType>::run(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        Assert(internal::has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, comm);

        try
          {
            // 1) Send data to identified targets and start receiving
            //    the answers from these very same processes.
            start_communication(targets, create_request, comm);

            // 2) Until all posted receive operations are known to have
            //    completed, answer requests and keep checking whether all
            //    requests of this process have been answered.
            //
            //    The requests that we catch in the answer_requests()
            //    function originate elsewhere, that is, they are not in
            //    response to our own messages
            //
            //    Note also that we may not catch all incoming requests in
            //    the following two lines: our own requests may have been
            //    satisfied before we've dealt with all incoming requests.
            //    That's ok: We will get around to dealing with all
            //    remaining message later. We just want to move on to the
            //    next step as early as possible.
            while (all_locally_originated_receives_are_completed(process_answer,
                                                                 comm) == false)
              maybe_answer_one_request(answer_request, comm);

            // 3) Signal to all other processes that all requests of this
            //    process have been answered
            signal_finish(comm);

            // 4) Nevertheless, this process has to keep on answering
            //    (potential) incoming requests until all processes have
            //    received the answer to all requests
            while (all_remotely_originated_receives_are_completed() == false)
              maybe_answer_one_request(answer_request, comm);

            // 5) process the answer to all requests
            clean_up_and_end_communication(comm);
          }
        catch (...)
          {
            internal::handle_exception(std::current_exception(), comm);
          }

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename RequestType, typename AnswerType>
      void
      NBX<RequestType, AnswerType>::start_communication(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const MPI_Comm                                        comm)
      {
#  ifdef DEAL_II_WITH_MPI
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
              send_buffer =
                (create_request ? Utilities::pack(create_request(rank), false) :
                                  std::vector<char>());

              // Post a request to send data
              auto ierr = MPI_Isend(send_buffer.data(),
                                    send_buffer.size(),
                                    MPI_CHAR,
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
#  else
        (void)targets;
        (void)create_request;
        (void)comm;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      bool
      NBX<RequestType, AnswerType>::
        all_locally_originated_receives_are_completed(
          const std::function<void(const unsigned int, const AnswerType &)>
                        &process_answer,
          const MPI_Comm comm)
      {
#  ifdef DEAL_II_WITH_MPI
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
                // OK, so we have gotten a reply to our request from
                // one rank. Let us process it.
                const auto target = status.MPI_SOURCE;

                // Then query the size of the message, allocate enough memory,
                // receive the data, and process it.
                int message_size;
                {
                  const int ierr =
                    MPI_Get_count(&status, MPI_CHAR, &message_size);
                  AssertThrowMPI(ierr);
                }
                std::vector<char> recv_buffer(message_size);

                {
                  const int tag_deliver = Utilities::MPI::internal::Tags::
                    consensus_algorithm_nbx_process_deliver;

                  const int ierr = MPI_Recv(recv_buffer.data(),
                                            recv_buffer.size(),
                                            MPI_CHAR,
                                            target,
                                            tag_deliver,
                                            comm,
                                            MPI_STATUS_IGNORE);
                  AssertThrowMPI(ierr);
                }

                if (process_answer)
                  process_answer(target,
                                 Utilities::unpack<AnswerType>(recv_buffer,
                                                               false));

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

#  else
        (void)process_answer;
        (void)comm;

        return true;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      NBX<RequestType, AnswerType>::maybe_answer_one_request(
        const std::function<AnswerType(const unsigned int, const RequestType &)>
                      &answer_request,
        const MPI_Comm comm)
      {
#  ifdef DEAL_II_WITH_MPI

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
            auto ierr = MPI_Get_count(&status, MPI_CHAR, &number_amount);
            AssertThrowMPI(ierr);

            // allocate memory for incoming message
            std::vector<char> buffer_recv(number_amount);
            ierr = MPI_Recv(buffer_recv.data(),
                            number_amount,
                            MPI_CHAR,
                            other_rank,
                            tag_request,
                            comm,
                            MPI_STATUS_IGNORE);
            AssertThrowMPI(ierr);

            // Allocate memory for an answer message to the current request,
            // and ask the 'process' object to produce an answer:
            request_buffers.emplace_back(std::make_unique<std::vector<char>>());
            auto &request_buffer = *request_buffers.back();
            if (answer_request)
              request_buffer =
                Utilities::pack(answer_request(other_rank,
                                               Utilities::unpack<RequestType>(
                                                 buffer_recv, false)),
                                false);

            // Then initiate sending the answer back to the requester.
            request_requests.emplace_back(std::make_unique<MPI_Request>());
            ierr = MPI_Isend(request_buffer.data(),
                             request_buffer.size(),
                             MPI_CHAR,
                             other_rank,
                             tag_deliver,
                             comm,
                             request_requests.back().get());
            AssertThrowMPI(ierr);
          }
#  else
        (void)answer_request;
        (void)comm;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      NBX<RequestType, AnswerType>::signal_finish(const MPI_Comm comm)
      {
#  ifdef DEAL_II_WITH_MPI
        const auto ierr = MPI_Ibarrier(comm, &barrier_request);
        AssertThrowMPI(ierr);
#  else
        (void)comm;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      bool
      NBX<RequestType,
          AnswerType>::all_remotely_originated_receives_are_completed()
      {
#  ifdef DEAL_II_WITH_MPI
        int        all_ranks_reached_barrier;
        const auto ierr = MPI_Test(&barrier_request,
                                   &all_ranks_reached_barrier,
                                   MPI_STATUS_IGNORE);
        AssertThrowMPI(ierr);
        return all_ranks_reached_barrier != 0;
#  else
        return true;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      NBX<RequestType, AnswerType>::clean_up_and_end_communication(
        const MPI_Comm comm)
      {
        (void)comm;
#  ifdef DEAL_II_WITH_MPI
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

          if constexpr (running_in_debug_mode())
            {
              // note: IBarrier seems to make problem during testing, this
              // additional Barrier seems to help
              ierr = MPI_Barrier(comm);
              AssertThrowMPI(ierr);
            }
        }
#  endif
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      PEX<RequestType, AnswerType>::run(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        Assert(internal::has_unique_elements(targets),
               ExcMessage("The consensus algorithms expect that each process "
                          "only sends a single message to another process, "
                          "but the targets provided include duplicates."));

        static CollectiveMutex      mutex;
        CollectiveMutex::ScopedLock lock(mutex, comm);

        try
          {
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
          }
        catch (...)
          {
            internal::handle_exception(std::current_exception(), comm);
          }

        return std::vector<unsigned int>(requesting_processes.begin(),
                                         requesting_processes.end());
      }



      template <typename RequestType, typename AnswerType>
      unsigned int
      PEX<RequestType, AnswerType>::start_communication(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const MPI_Comm                                        comm)
      {
#  ifdef DEAL_II_WITH_MPI
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
            if (create_request)
              send_buffer = Utilities::pack(create_request(rank), false);

            // start to send data
            auto ierr = MPI_Isend(send_buffer.data(),
                                  send_buffer.size(),
                                  MPI_CHAR,
                                  rank,
                                  tag_request,
                                  comm,
                                  &send_request_requests[i]);
            AssertThrowMPI(ierr);
          }

        return n_sources;
#  else
        (void)targets;
        (void)create_request;
        (void)comm;
        return 0;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      PEX<RequestType, AnswerType>::answer_one_request(
        const unsigned int index,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
                      &answer_request,
        const MPI_Comm comm)
      {
#  ifdef DEAL_II_WITH_MPI
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

        // Actually get the incoming message:
        int number_amount;
        ierr = MPI_Get_count(&status, MPI_CHAR, &number_amount);
        AssertThrowMPI(ierr);

        std::vector<char> buffer_recv(number_amount);
        ierr = MPI_Recv(buffer_recv.data(),
                        number_amount,
                        MPI_CHAR,
                        other_rank,
                        tag_request,
                        comm,
                        &status);
        AssertThrowMPI(ierr);

        // Process request by asking the user-provided function for
        // the answer and post a send for it.
        auto &request_buffer = requests_buffers[index];
        request_buffer =
          (answer_request ?
             Utilities::pack(answer_request(other_rank,
                                            Utilities::unpack<RequestType>(
                                              buffer_recv, false)),
                             false) :
             std::vector<char>());

        ierr = MPI_Isend(request_buffer.data(),
                         request_buffer.size(),
                         MPI_CHAR,
                         other_rank,
                         tag_deliver,
                         comm,
                         &send_answer_requests[index]);
        AssertThrowMPI(ierr);
#  else
        (void)answer_request;
        (void)comm;
        (void)index;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      PEX<RequestType, AnswerType>::process_incoming_answers(
        const unsigned int n_targets,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
#  ifdef DEAL_II_WITH_MPI
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
              const int ierr = MPI_Get_count(&status, MPI_CHAR, &message_size);
              AssertThrowMPI(ierr);
            }
            std::vector<char> recv_buffer(message_size);

            // Now actually receive the answer. Because the MPI_Probe
            // above blocks until we have a message, we know that the
            // following MPI_Recv call will immediately succeed.
            {
              const int ierr = MPI_Recv(recv_buffer.data(),
                                        recv_buffer.size(),
                                        MPI_CHAR,
                                        other_rank,
                                        tag_deliver,
                                        comm,
                                        MPI_STATUS_IGNORE);
              AssertThrowMPI(ierr);
            }

            if (process_answer)
              process_answer(other_rank,
                             Utilities::unpack<AnswerType>(recv_buffer, false));
          }
#  else
        (void)n_targets;
        (void)process_answer;
        (void)comm;
#  endif
      }



      template <typename RequestType, typename AnswerType>
      void
      PEX<RequestType, AnswerType>::clean_up_and_end_communication()
      {
#  ifdef DEAL_II_WITH_MPI
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
#  endif
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      Serial<RequestType, AnswerType>::run(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        Assert(Utilities::MPI::n_mpi_processes(comm) == 1,
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
            const RequestType request =
              (create_request ? create_request(0) : RequestType());
            const AnswerType answer =
              (answer_request ? answer_request(0, request) : AnswerType());

            if (process_answer)
              process_answer(0, answer);
          }

        return targets; // nothing to do
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      Selector<RequestType, AnswerType>::run(
        const std::vector<unsigned int>                      &targets,
        const std::function<RequestType(const unsigned int)> &create_request,
        const std::function<AnswerType(const unsigned int, const RequestType &)>
          &answer_request,
        const std::function<void(const unsigned int, const AnswerType &)>
                      &process_answer,
        const MPI_Comm comm)
      {
        // Depending on the number of processes we switch between
        // implementations. We reduce the threshold for debug mode to be
        // able to test also the non-blocking implementation. This feature
        // is tested by:
        // tests/multigrid/transfer_matrix_free_06.with_mpi=true.with_p4est=true.with_trilinos=true.mpirun=10.output

        const unsigned int n_procs = (Utilities::MPI::job_supports_mpi() ?
                                        Utilities::MPI::n_mpi_processes(comm) :
                                        1);
#  ifdef DEAL_II_WITH_MPI
#    ifdef DEBUG
        if (n_procs > 10)
#    else
        if (n_procs > 99)
#    endif
          consensus_algo.reset(new NBX<RequestType, AnswerType>());
        else
#  endif
          if (n_procs > 1)
          consensus_algo.reset(new PEX<RequestType, AnswerType>());
        else
          consensus_algo.reset(new Serial<RequestType, AnswerType>());

        return consensus_algo->run(
          targets, create_request, answer_request, process_answer, comm);
      }


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities

#endif // DOXYGEN


DEAL_II_NAMESPACE_CLOSE

#endif
