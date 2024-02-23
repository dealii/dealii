// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


/*
 * check GraphColoring::color_sparsity_pattern with 5 nodes with
 * no connections with each other. All nodes will be colored with
 * one color.
 */


#include <deal.II/base/graph_coloring.h>

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"


void
fill_graph(DynamicSparsityPattern &graph)
{
  // Edges in only one direction
  graph.add(0, 0);
  graph.add(1, 1);
  graph.add(2, 2);
  graph.add(3, 3);
  graph.add(4, 4);
}


int
main(int argc, char **argv)
{
  // Initialize MPI and Zoltan
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    all;

  // Number of nodes
  unsigned int num_indices = 5;

  // Create temporary object to hold graph connection info.
  DynamicSparsityPattern dynamic_sparse_graph;
  dynamic_sparse_graph.reinit(num_indices, num_indices);

  // fill dynamic sparsity pattern
  fill_graph(dynamic_sparse_graph);

  // Create sparity pattern to hold graph connection info.
  SparsityPattern sp_graph;
  sp_graph.copy_from(dynamic_sparse_graph);

  Assert(num_indices == sp_graph.n_rows(), ExcInternalError());

  std::vector<unsigned int> color_indices;
  unsigned int              num_colors;

  color_indices.resize(num_indices);
  num_colors = GraphColoring::color_sparsity_pattern(sp_graph, color_indices);

  // color
  deallog << "Coloring" << std::endl;
  deallog << "Number of colors used: " << num_colors << std::endl;
  for (unsigned int i = 0; i < num_indices; ++i)
    {
      deallog << i << ' ' << color_indices[i] << std::endl;
    }
}
