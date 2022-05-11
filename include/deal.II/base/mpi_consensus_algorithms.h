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
     * This namespace provides several implementations of consensus algorithms,
     * such as the nbx(), pex(), serial(), and selector() functions.
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
       * @tparam RequestType The type of the elements of the vector to sent.
       * @tparam AnswerType The type of the elements of the vector to received.
       *
       * @note Since the payloads of the messages are optional, users have
       *    to deal with buffers themselves. The ConsensusAlgorithm classes
       *    (1) deliver only references to empty vectors (of size 0) the data
       *    to be sent can be inserted to or read from, and (2) communicate
       *    these vectors blindly.
       */
      template <typename RequestType, typename AnswerType>
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
        create_request(const unsigned int        other_rank,
                       std::vector<RequestType> &send_buffer);

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
        answer_request(const unsigned int              other_rank,
                       const std::vector<RequestType> &buffer_recv,
                       std::vector<AnswerType> &       request_buffer);

        /**
         * Process the payload of the answer of the request to the process with
         * the specified rank.
         *
         * @param[in] other_rank rank of the process
         * @param[in] recv_buffer data to be sent part of the request (optional)
         */
        virtual void
        read_answer(const unsigned int             other_rank,
                    const std::vector<AnswerType> &recv_buffer);
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
        Interface();

        /**
         * Constructor. @p process is an object that provides information
         * about what processes the current process wants to communicate with,
         * and the data to be sent/received. @p comm is the communicator on
         * which this communication is to happen.
         *
         * @deprecated This constructor stores the Process object and the
         *   communicator so that one can later call the run() function
         *   without arguments. This approach is deprecated. Instead, use
         *   the default constructor of this class along with the run()
         *   function that takes an argument.
         */
        DEAL_II_DEPRECATED
        Interface(Process<RequestType, AnswerType> &process,
                  const MPI_Comm &                  comm);

        /**
         * Destructor. Made `virtual` to ensure that one can work with
         * derived classes.
         */
        virtual ~Interface() = default;

        /**
         * Run the consensus algorithm and return a vector of process ranks
         * that have requested answers from the current process.
         *
         * @deprecated This function is deprecated. It can be called
         *   if the Process object and communicator to be used have previously
         *   been provided to the non-default constructor. Use the run()
         *   functions taking arguments instead.
         */
        DEAL_II_DEPRECATED
        std::vector<unsigned int>
        run();

        /**
         * Run the consensus algorithm and return a vector of process ranks
         * that have requested answers from the current process.
         *
         * This version of the run() function simply unpacks the functions
         * packaged in `process` and calls the version of the run() function
         * that takes a number of `std::function` arguments.
         */
        std::vector<unsigned int>
        run(Process<RequestType, AnswerType> &process, const MPI_Comm &comm);

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
        run(const std::vector<unsigned int> &targets,
            const std::function<std::vector<RequestType>(const unsigned int)>
              &                                   create_request,
            const std::function<std::vector<AnswerType>(
              const unsigned int,
              const std::vector<RequestType> &)> &answer_request,
            const std::function<void(const unsigned int,
                                     const std::vector<AnswerType> &)>
              &             process_answer,
            const MPI_Comm &comm) = 0;

      private:
        /**
         * Reference to the process provided by the user.
         *
         * This member variable is only used in the deprecated constructor
         * and the run() function without argument. It is a `nullptr`
         * otherwise
         */
        DEAL_II_DEPRECATED
        Process<RequestType, AnswerType> *process;

        /**
         * MPI communicator.
         *
         * This member variable is only used in the deprecated constructor
         * and the run() function without argument.
         */
        DEAL_II_DEPRECATED
        MPI_Comm comm;
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
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator
         *
         * @deprecated This constructor stores the Process object and the
         *   communicator so that one can later call the run() function
         *   without arguments. This approach is deprecated. Instead, use
         *   the default constructor of this class along with the run()
         *   function that takes an argument.
         */
        DEAL_II_DEPRECATED
        NBX(Process<RequestType, AnswerType> &process, const MPI_Comm &comm);

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
        run(const std::vector<unsigned int> &targets,
            const std::function<std::vector<RequestType>(const unsigned int)>
              &                                   create_request,
            const std::function<std::vector<AnswerType>(
              const unsigned int,
              const std::vector<RequestType> &)> &answer_request,
            const std::function<void(const unsigned int,
                                     const std::vector<AnswerType> &)>
              &             process_answer,
            const MPI_Comm &comm) override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<RequestType>> send_buffers;

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
        std::vector<std::unique_ptr<std::vector<AnswerType>>> request_buffers;

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
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm);

        /**
         * Signal to all other ranks that this rank has received all request
         * answers via entering IBarrier.
         */
        void
        signal_finish(const MPI_Comm &comm);

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
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const MPI_Comm &                      comm);

        /**
         * Start to send all requests via ISend and post IRecvs for the incoming
         * answer messages.
         */
        void
        start_communication(
          const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &             create_request,
          const MPI_Comm &comm);

        /**
         * After all rank has received all answers, the MPI data structures can
         * be freed and the received answers can be processed.
         */
        void
        clean_up_and_end_communication(const MPI_Comm &comm);
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
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &                                   create_request,
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm);


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
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator
         *
         * @deprecated This constructor stores the Process object and the
         *   communicator so that one can later call the run() function
         *   without arguments. This approach is deprecated. Instead, use
         *   the default constructor of this class along with the run()
         *   function that takes an argument.
         */
        DEAL_II_DEPRECATED
        PEX(Process<RequestType, AnswerType> &process, const MPI_Comm &comm);

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
        run(const std::vector<unsigned int> &targets,
            const std::function<std::vector<RequestType>(const unsigned int)>
              &                                   create_request,
            const std::function<std::vector<AnswerType>(
              const unsigned int,
              const std::vector<RequestType> &)> &answer_request,
            const std::function<void(const unsigned int,
                                     const std::vector<AnswerType> &)>
              &             process_answer,
            const MPI_Comm &comm) override;

      private:
#ifdef DEAL_II_WITH_MPI
        /**
         * Buffers for sending requests.
         */
        std::vector<std::vector<RequestType>> send_buffers;

        /**
         * Buffers for receiving answers to requests.
         */
        std::vector<std::vector<AnswerType>> recv_buffers;

        /**
         * MPI request objects for sending request messages.
         */
        std::vector<MPI_Request> send_request_requests;

        /**
         * Buffers for sending answers to requests.
         */
        std::vector<std::vector<AnswerType>> requests_buffers;

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
          const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &             create_request,
          const MPI_Comm &comm);

        /**
         * The `index`th request message from another rank has been received:
         * process the request and send an answer.
         */
        void
        answer_one_request(const unsigned int                    index,
                           const std::function<std::vector<AnswerType>(
                             const unsigned int,
                             const std::vector<RequestType> &)> &answer_request,
                           const MPI_Comm &                      comm);

        /**
         * Receive and process all of the incoming responses to the
         * requests we sent.
         */
        void
        process_incoming_answers(
          const unsigned int n_targets,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm);

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
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &                                   create_request,
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm);


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

        /**
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator (ignored)
         *
         * @deprecated This constructor stores the Process object and the
         *   communicator so that one can later call the run() function
         *   without arguments. This approach is deprecated. Instead, use
         *   the default constructor of this class along with the run()
         *   function that takes an argument.
         */
        DEAL_II_DEPRECATED
        Serial(Process<RequestType, AnswerType> &process, const MPI_Comm &comm);

        // Import the declarations from the base class.
        using Interface<RequestType, AnswerType>::run;

        /**
         * @copydoc Interface::run()
         */
        virtual std::vector<unsigned int>
        run(const std::vector<unsigned int> &targets,
            const std::function<std::vector<RequestType>(const unsigned int)>
              &                                   create_request,
            const std::function<std::vector<AnswerType>(
              const unsigned int,
              const std::vector<RequestType> &)> &answer_request,
            const std::function<void(const unsigned int,
                                     const std::vector<AnswerType> &)>
              &             process_answer,
            const MPI_Comm &comm) override;
      };



      /**
       * This function implements a concrete algorithm for the
       * consensus algorithms problem (see the documentation of the
       * surrounding namespace), as a fall-back option for the case
       * where the communicator provided has only one rank (or when
       * MPI is simply not used at all).
       *
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &                                   create_request,
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm);



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
         * Constructor.
         *
         * @param process Process to be run during consensus algorithm.
         * @param comm MPI Communicator.
         *
         * @deprecated This constructor stores the Process object and the
         *   communicator so that one can later call the run() function
         *   without arguments. This approach is deprecated. Instead, use
         *   the default constructor of this class along with the run()
         *   function that takes an argument.
         */
        DEAL_II_DEPRECATED
        Selector(Process<RequestType, AnswerType> &process,
                 const MPI_Comm &                  comm);

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
        run(const std::vector<unsigned int> &targets,
            const std::function<std::vector<RequestType>(const unsigned int)>
              &                                   create_request,
            const std::function<std::vector<AnswerType>(
              const unsigned int,
              const std::vector<RequestType> &)> &answer_request,
            const std::function<void(const unsigned int,
                                     const std::vector<AnswerType> &)>
              &             process_answer,
            const MPI_Comm &comm) override;

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
       * @tparam RequestType The type of the elements of the vector to be sent.
       * @tparam AnswerType The type of the elements of the vector to be received.
       */
      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      selector(const std::vector<unsigned int> &targets,
               const std::function<std::vector<RequestType>(const unsigned int)>
                 &                                   create_request,
               const std::function<std::vector<AnswerType>(
                 const unsigned int,
                 const std::vector<RequestType> &)> &answer_request,
               const std::function<void(const unsigned int,
                                        const std::vector<AnswerType> &)>
                 &             process_answer,
               const MPI_Comm &comm);



      /**
       * This class implements Utilities::MPI::ConsensusAlgorithms::Process,
       * using user-provided function wrappers.
       * The advantage of this class is that users do not have to write their
       * own implementation but can register lambda functions directly.
       */
      template <typename RequestType, typename AnswerType>
      class AnonymousProcess : public Process<RequestType, AnswerType>
      {
      public:
        /**
         * Register functions that should be called for implementing the
         * interface of Process.
         *
         * @param function_compute_targets called during `compute_targets`.
         * @param function_create_request called during `create_request`.
         * @param function_answer_request called during `answer_request`.
         * @param function_read_answer called during `read_answer`.
         */
        AnonymousProcess(
          const std::function<std::vector<unsigned int>()>
            &function_compute_targets,
          const std::function<void(const unsigned int,
                                   std::vector<RequestType> &)>
            &function_create_request = {},
          const std::function<void(const unsigned int,
                                   const std::vector<RequestType> &,
                                   std::vector<AnswerType> &)>
            &function_answer_request = {},
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
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
        create_request(const unsigned int        other_rank,
                       std::vector<RequestType> &send_buffer) override;

        /**
         * @copydoc Process::answer_request()
         */
        void
        answer_request(const unsigned int              other_rank,
                       const std::vector<RequestType> &buffer_recv,
                       std::vector<AnswerType> &       request_buffer) override;

        /**
         * @copydoc Process::read_answer()
         */
        void
        read_answer(const unsigned int             other_rank,
                    const std::vector<AnswerType> &recv_buffer) override;

      private:
        const std::function<std::vector<unsigned int>()>
          function_compute_targets;
        const std::function<void(const int, std::vector<RequestType> &)>
          function_create_request;
        const std::function<void(const unsigned int,
                                 const std::vector<RequestType> &,
                                 std::vector<AnswerType> &)>
          function_answer_request;
        const std::function<void(const int, const std::vector<AnswerType> &)>
          function_read_answer;
      };


#ifndef DOXYGEN
      // Implementation of the functions in this namespace.

      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      nbx(const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &                                   create_request,
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm)
      {
        return NBX<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      pex(const std::vector<unsigned int> &targets,
          const std::function<std::vector<RequestType>(const unsigned int)>
            &                                   create_request,
          const std::function<std::vector<AnswerType>(
            const unsigned int,
            const std::vector<RequestType> &)> &answer_request,
          const std::function<void(const unsigned int,
                                   const std::vector<AnswerType> &)>
            &             process_answer,
          const MPI_Comm &comm)
      {
        return PEX<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      serial(const std::vector<unsigned int> &targets,
             const std::function<std::vector<RequestType>(const unsigned int)>
               &                                   create_request,
             const std::function<std::vector<AnswerType>(
               const unsigned int,
               const std::vector<RequestType> &)> &answer_request,
             const std::function<void(const unsigned int,
                                      const std::vector<AnswerType> &)>
               &             process_answer,
             const MPI_Comm &comm)
      {
        return Serial<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      selector(const std::vector<unsigned int> &targets,
               const std::function<std::vector<RequestType>(const unsigned int)>
                 &                                   create_request,
               const std::function<std::vector<AnswerType>(
                 const unsigned int,
                 const std::vector<RequestType> &)> &answer_request,
               const std::function<void(const unsigned int,
                                        const std::vector<AnswerType> &)>
                 &             process_answer,
               const MPI_Comm &comm)
      {
        return Selector<RequestType, AnswerType>().run(
          targets, create_request, answer_request, process_answer, comm);
      }



      template <typename RequestType, typename AnswerType>
      AnonymousProcess<RequestType, AnswerType>::AnonymousProcess(
        const std::function<std::vector<unsigned int>()>
          &function_compute_targets,
        const std::function<void(const unsigned int,
                                 std::vector<RequestType> &)>
          &function_create_request,
        const std::function<void(const unsigned int,
                                 const std::vector<RequestType> &,
                                 std::vector<AnswerType> &)>
          &function_answer_request,
        const std::function<void(const unsigned int,
                                 const std::vector<AnswerType> &)>
          &function_read_answer)
        : function_compute_targets(function_compute_targets)
        , function_create_request(function_create_request)
        , function_answer_request(function_answer_request)
        , function_read_answer(function_read_answer)
      {}



      template <typename RequestType, typename AnswerType>
      std::vector<unsigned int>
      AnonymousProcess<RequestType, AnswerType>::compute_targets()
      {
        return function_compute_targets();
      }



      template <typename RequestType, typename AnswerType>
      void
      AnonymousProcess<RequestType, AnswerType>::create_request(
        const unsigned int        other_rank,
        std::vector<RequestType> &send_buffer)
      {
        if (function_create_request)
          function_create_request(other_rank, send_buffer);
      }



      template <typename RequestType, typename AnswerType>
      void
      AnonymousProcess<RequestType, AnswerType>::answer_request(
        const unsigned int              other_rank,
        const std::vector<RequestType> &buffer_recv,
        std::vector<AnswerType> &       request_buffer)
      {
        if (function_answer_request)
          function_answer_request(other_rank, buffer_recv, request_buffer);
      }



      template <typename RequestType, typename AnswerType>
      void
      AnonymousProcess<RequestType, AnswerType>::read_answer(
        const unsigned int             other_rank,
        const std::vector<AnswerType> &recv_buffer)
      {
        if (function_read_answer)
          function_read_answer(other_rank, recv_buffer);
      }

#endif


    } // namespace ConsensusAlgorithms
  }   // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
