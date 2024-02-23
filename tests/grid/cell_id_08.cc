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


// benchmark for creation and serialization of CellId objects

#include <deal.II/base/timer.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/cell_id.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/iostreams/device/array.hpp>
#include <boost/iostreams/device/back_inserter.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <iostream>
#include <vector>

#include "../tests.h"


template <typename Container>
std::size_t
size_in_bytes(const Container &container)
{
  return sizeof(Container) +
         (sizeof(typename Container::value_type) * container.size());
}



template <int dim>
void
benchmark(const unsigned int n_global_refinements)
{
  // be generous with the estimate
  const std::size_t size_estimate_per_cellid_in_bytes = 96;

  TimerOutput computing_timer(std::cout,
                              TimerOutput::summary,
                              TimerOutput::wall_times);

  Triangulation<dim> tria;
  GridGenerator::hyper_cube(tria);
  tria.refine_global(n_global_refinements);

  const unsigned int n_cells = tria.n_active_cells();


  //
  // create
  //
  std::vector<CellId> cell_ids;
  cell_ids.reserve(n_cells);
  {
    TimerOutput::Scope t(computing_timer, "cell->id");

    for (const auto &cell : tria.active_cell_iterators())
      cell_ids.push_back(cell->id());
  }
  deallog << "Creation successful." << std::endl;
  std::cout << "Number of CellId objects: " << cell_ids.size() << std::endl;


  //
  // boost binary archive
  //
  {
    std::string buffer;
    buffer.reserve(n_cells * size_estimate_per_cellid_in_bytes);
    std::vector<CellId> deserialized;
    deserialized.reserve(n_cells);
    {
      TimerOutput::Scope t(computing_timer, "boost binary_archive save");

      boost::iostreams::filtering_ostreambuf fosb;
      fosb.push(boost::iostreams::back_inserter(buffer));
      boost::archive::binary_oarchive oa(fosb);
      oa                             &cell_ids;
    }
    {
      TimerOutput::Scope t(computing_timer, "boost binary_archive load");

      boost::iostreams::filtering_istreambuf fisb;
      fisb.push(boost::iostreams::array_source(buffer.data(), buffer.size()));
      boost::archive::binary_iarchive ia(fisb);
      ia                             &deserialized;
    }
    AssertThrow(cell_ids == deserialized,
                ExcMessage("Serialization with boost binary archives failed."));
    deallog << "Serialization with boost binary archives successful."
            << std::endl;
    std::cout << "Buffer size in bytes with ..."
              << "  boost binary_archive : " << size_in_bytes(buffer)
              << std::endl;
  }


  //
  // Utilities uncompressed
  //
  {
    std::vector<char> buffer;
    buffer.reserve(n_cells * size_estimate_per_cellid_in_bytes);
    std::vector<std::size_t> sizes;
    sizes.reserve(n_cells);
    std::vector<CellId> deserialized;
    deserialized.reserve(n_cells);
    {
      TimerOutput::Scope t(computing_timer, "Utilities::pack uncompressed");

      for (const auto &cell_id : cell_ids)
        sizes.push_back(
          Utilities::pack(cell_id, buffer, /*allow_compression*/ false));
    }
    {
      TimerOutput::Scope t(computing_timer, "Utilities::unpack uncompressed");

      std::vector<char>::const_iterator cbegin = buffer.cbegin(), cend;
      for (const auto &size : sizes)
        {
          cend = cbegin + size;
          deserialized.push_back(Utilities::unpack<CellId>(
            cbegin, cend, /*allow_compression*/ false));
          cbegin = cend;
        }
      Assert(cend == buffer.cend(), ExcInternalError());
    }
    AssertThrow(cell_ids == deserialized,
                ExcMessage("Serialization with Utilities failed."));
    deallog << "Serialization with Utilities successful." << std::endl;
    std::cout << "  Utilities            : " << size_in_bytes(buffer)
              << std::endl;
  }


  //
  // binary representation
  //
  {
    std::vector<CellId::binary_type> buffer;
    buffer.reserve(n_cells);
    std::vector<CellId> deserialized;
    deserialized.reserve(n_cells);
    {
      TimerOutput::Scope t(computing_timer, "CellId::to_binary");

      for (const auto &cell_id : cell_ids)
        buffer.push_back(cell_id.to_binary<dim>());
    }
    {
      TimerOutput::Scope t(computing_timer, "CellId from binary");

      for (const auto &binary : buffer)
        deserialized.emplace_back(binary);
    }
    AssertThrow(cell_ids == deserialized,
                ExcMessage("Serialization with binary representation failed."));
    deallog << "Serialization with binary representation successful."
            << std::endl;
    std::cout << "  binary representation: " << size_in_bytes(buffer)
              << std::endl;
  }

  deallog << "OK" << std::endl;
}



int
main()
{
  initlog();

  benchmark<2>(4);
  // benchmark<3>(8);
}
