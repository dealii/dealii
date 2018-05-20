// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

// check wrapper function GraphColoring::color_sparsity_pattern

#include "../tests.h"
#include <deal.II/base/graph_coloring.h>

void
fill_graph(DynamicSparsityPattern& graph)
{
  //Edges in only one direction
  graph.add(0, 1);
  graph.add(0, 2);
  graph.add(1, 2);
  graph.add(1, 4);
  graph.add(4, 3);
  graph.add(2, 4);
  graph.add(3, 2);

  //Edges in opposite direction
  graph.add(1, 0);
  graph.add(2, 0);
  graph.add(2, 1);
  graph.add(4, 1);
  graph.add(3, 4);
  graph.add(4, 2);
  graph.add(2, 3);
}

int
main(int argc, char** argv)
{
  //Initialize MPI and Zoltan
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv);
  MPILogInitAll                    all;

  //Number of nodes
  unsigned int num_indices = 5;

  //Create temporary object to hold graph connection info.
  DynamicSparsityPattern dynamic_sparse_graph;
  dynamic_sparse_graph.reinit(num_indices, num_indices);

  //fill dynamic sparsity pattern
  fill_graph(dynamic_sparse_graph);

  //Create sparity pattern to hold graph connection info.
  SparsityPattern sp_graph;
  sp_graph.copy_from(dynamic_sparse_graph);

  Assert(num_indices == sp_graph.n_rows(), ExcInternalError());

  std::vector<unsigned int> color_indices;
  unsigned int              num_colors;

  color_indices.resize(num_indices);
  num_colors = GraphColoring::color_sparsity_pattern(sp_graph, color_indices);

  //color
  deallog << "Coloring" << std::endl;
  deallog << "Number of colors used: " << num_colors << std::endl;
  for(unsigned int i = 0; i < num_indices; i++)
    {
      deallog << i << " " << color_indices[i] << std::endl;
    }
}
