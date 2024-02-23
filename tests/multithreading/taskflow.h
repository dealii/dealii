// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2020 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


#include <taskflow/taskflow.hpp>

namespace taskflow_v1
{
  template <typename Worker,
            typename Copier,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator                             &begin,
      const std_cxx20::type_identity_t<Iterator> &end,
      Worker                                      worker,
      Copier                                      copier,
      const ScratchData                          &sample_scratch_data,
      const CopyData                             &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    if (MultithreadInfo::n_threads() == 1)
      {
        // need to copy the sample since it is marked const
        ScratchData scratch_data = sample_scratch_data;
        CopyData    copy_data    = sample_copy_data; // NOLINT

        for (Iterator i = begin; i != end; ++i)
          {
            // need to check if the function is not the zero function. To
            // check zero-ness, create a C++ function out of it and check that
            if (static_cast<const std::function<
                  void(const Iterator &, ScratchData &, CopyData &)> &>(worker))
              worker(i, scratch_data, copy_data);
            if (static_cast<const std::function<void(const CopyData &)> &>(
                  copier))
              copier(copy_data);
          }

        return;
      }

    tf::Executor &executor = MultithreadInfo::get_taskflow_executor();
    tf::Taskflow  taskflow;

    ScratchData scratch_data = sample_scratch_data;
    CopyData    copy_data    = sample_copy_data; // NOLINT

    tf::Task last_copier;

    std::vector<std::unique_ptr<CopyData>> copy_datas;

    unsigned int idx = 0;
    for (Iterator i = begin; i != end; ++i, ++idx)
      {
        copy_datas.emplace_back();

        auto worker_task = taskflow
                             .emplace([it = i,
                                       idx,
                                       &sample_scratch_data,
                                       &copy_datas,
                                       &sample_copy_data,
                                       &worker]() {
                               // std::cout << "worker " << idx << std::endl;
                               ScratchData scratch = sample_scratch_data;
                               auto       &copy    = copy_datas[idx];
                               copy =
                                 std::make_unique<CopyData>(sample_copy_data);

                               worker(it, scratch, *copy.get());
                             })
                             .name("worker");

        tf::Task copier_task = taskflow
                                 .emplace([idx, &copy_datas, &copier]() {
                                   copier(*copy_datas[idx].get());
                                   copy_datas[idx].reset();
                                 })
                                 .name("copy");

        worker_task.precede(copier_task);

        if (!last_copier.empty())
          last_copier.precede(copier_task);
        last_copier = copier_task;
      }

    executor.run(taskflow).wait();
    if (false)
      {
        std::ofstream f("graph.dia");
        taskflow.dump(f);
        f.close();
      }
  }


  template <typename MainClass,
            typename Iterator,
            typename ScratchData,
            typename CopyData>
  void
  run(const Iterator                             &begin,
      const std_cxx20::type_identity_t<Iterator> &end,
      MainClass                                  &main_object,
      void (MainClass::*worker)(const Iterator &, ScratchData &, CopyData &),
      void (MainClass::*copier)(const CopyData &),
      const ScratchData &sample_scratch_data,
      const CopyData    &sample_copy_data,
      const unsigned int queue_length = 2 * MultithreadInfo::n_threads(),
      const unsigned int chunk_size   = 8)
  {
    // forward to the other function
    run(
      begin,
      end,
      [&main_object, worker](const Iterator &iterator,
                             ScratchData    &scratch_data,
                             CopyData       &copy_data) {
        (main_object.*worker)(iterator, scratch_data, copy_data);
      },
      [&main_object, copier](const CopyData &copy_data) {
        (main_object.*copier)(copy_data);
      },
      sample_scratch_data,
      sample_copy_data,
      queue_length,
      chunk_size);
  }
} // namespace taskflow_v1
