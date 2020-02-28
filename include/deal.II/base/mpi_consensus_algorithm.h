// ---------------------------------------------------------------------
//
// Copyright (C) 2020 by the deal.II authors
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
     * This class implements Utilities::MPI::ConsensusAlgorithmProcess,
     * using user-provided function wrappers.
     * The advantage of this class is that users do not have to write their
     * own implementation but can register lambda functions directly.
     */
    template <typename T1, typename T2>
    class AnonymousConsensusAlgorithmProcess
      : public ConsensusAlgorithmProcess<T1, T2>
    {
    public:
      /**
       * Register functions that should be called for implementing the interface
       * of ConsensusAlgorithmProcess.
       *
       * @param function_compute_targets called during `compute_targets`.
       * @param function_create_request called during `create_request`.
       * @param function_answer_request called during `answer_request`.
       * @param function_prepare_buffer_for_answer called during
       *   `prepare_buffer_for_answer`.
       * @param function_read_answer called during `read_answer`.
       */
      AnonymousConsensusAlgorithmProcess(
        const std::function<std::vector<unsigned int>()>
          &function_compute_targets,
        const std::function<void(const unsigned int, std::vector<T1> &)>
          &function_create_request =
            [](const unsigned int, std::vector<T1> &) {},
        const std::function<void(const unsigned int,
                                 const std::vector<T1> &,
                                 std::vector<T2> &)> &function_answer_request =
          [](const unsigned int, const std::vector<T1> &, std::vector<T2> &) {},
        const std::function<void(const unsigned int, std::vector<T2> &)>
          &function_prepare_buffer_for_answer =
            [](const unsigned int, std::vector<T2> &) {},
        const std::function<void(const unsigned int, const std::vector<T2> &)>
          &function_read_answer =
            [](const unsigned int, const std::vector<T2> &) {});

      /**
       * @copydoc ConsensusAlgorithmProcess::compute_targets()
       */
      std::vector<unsigned int>
      compute_targets() override;

      /**
       * @copydoc ConsensusAlgorithmProcess::create_request()
       */
      void
      create_request(const unsigned int other_rank,
                     std::vector<T1> &  send_buffer) override;

      /**
       * @copydoc ConsensusAlgorithmProcess::answer_request()
       */
      void
      answer_request(const unsigned int     other_rank,
                     const std::vector<T1> &buffer_recv,
                     std::vector<T2> &      request_buffer) override;

      /**
       * @copydoc ConsensusAlgorithmProcess::prepare_buffer_for_answer()
       */
      void
      prepare_buffer_for_answer(const unsigned int other_rank,
                                std::vector<T2> &  recv_buffer) override;

      /**
       * @copydoc ConsensusAlgorithmProcess::read_answer()
       */
      void
      read_answer(const unsigned int     other_rank,
                  const std::vector<T2> &recv_buffer) override;

    private:
      const std::function<std::vector<unsigned int>()> function_compute_targets;
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
    AnonymousConsensusAlgorithmProcess<T1, T2>::
      AnonymousConsensusAlgorithmProcess(
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
    AnonymousConsensusAlgorithmProcess<T1, T2>::compute_targets()
    {
      return function_compute_targets();
    }



    template <typename T1, typename T2>
    void
    AnonymousConsensusAlgorithmProcess<T1, T2>::create_request(
      const unsigned int other_rank,
      std::vector<T1> &  send_buffer)
    {
      function_create_request(other_rank, send_buffer);
    }



    template <typename T1, typename T2>
    void
    AnonymousConsensusAlgorithmProcess<T1, T2>::answer_request(
      const unsigned int     other_rank,
      const std::vector<T1> &buffer_recv,
      std::vector<T2> &      request_buffer)
    {
      function_answer_request(other_rank, buffer_recv, request_buffer);
    }



    template <typename T1, typename T2>
    void
    AnonymousConsensusAlgorithmProcess<T1, T2>::prepare_buffer_for_answer(
      const unsigned int other_rank,
      std::vector<T2> &  recv_buffer)
    {
      function_prepare_buffer_for_answer(other_rank, recv_buffer);
    }



    template <typename T1, typename T2>
    void
    AnonymousConsensusAlgorithmProcess<T1, T2>::read_answer(
      const unsigned int     other_rank,
      const std::vector<T2> &recv_buffer)
    {
      function_read_answer(other_rank, recv_buffer);
    }


  } // end of namespace MPI
} // end of namespace Utilities


DEAL_II_NAMESPACE_CLOSE

#endif
