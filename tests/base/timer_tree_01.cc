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

// tests basic functionality of TimerTree

#include <deal.II/base/mpi.h>
#include <deal.II/base/timer_tree.h>

#include <fstream>
#include <iostream>
#include <sstream>

using namespace dealii;

void
test1()
{
  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // clang-format off
  pcout << std::endl
        << std::endl
        << "Timer: test 1 (basic)"
        << std::endl
        << std::endl;
  // clang-format on

  Timer timer;
  timer.restart();

  TimerTree tree;

  for (unsigned int i = 0; i < 1000; ++i)
    {
      tree.insert({"General"}, 20.0);

      tree.insert({"General", "Part 1"}, 2.0);

      tree.insert({"General", "Part 2"}, 3.0);

      tree.insert({"General", "Part 2", "Sub a"}, 0.75);

      tree.insert({"General", "Part 2", "Sub b"}, 0.9);

      tree.insert({"General", "Part 3"}, 4.0);

      tree.insert({"General", "Part 3", "Sub a"}, 0.5);

      tree.insert({"General", "Part 3", "Sub a", "sub-sub a"}, 0.04);

      tree.insert({"General", "Part 3", "Sub b"}, 0.98765);

      tree.insert({"General"}, 2.0);
    }

  pcout << std::endl << "timings for level = 0:" << std::endl;
  tree.print_level(pcout, 0);
  pcout << std::endl << "timings for level = 1:" << std::endl;
  tree.print_level(pcout, 1);
  pcout << std::endl << "timings for level = 2:" << std::endl;
  tree.print_level(pcout, 2);
  pcout << std::endl << "timings all:" << std::endl;
  tree.print_plain(pcout);
}

void
test2()
{
  ConditionalOStream pcout(std::cout,
                           Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) ==
                             0);

  // clang-format off
  pcout << std::endl
        << std::endl
        << "TimerTree: test 2 (modular coupling)"
        << std::endl
        << std::endl;
  // clang-format on

  Timer timer;
  timer.restart();

  TimerTree tree;
  tree.insert({"FSI"}, 100.);

  std::shared_ptr<TimerTree> tree_fluid;
  tree_fluid.reset(new TimerTree());
  // overall time can be inserted first ...
  tree_fluid->insert({"Fluid"}, 70.);
  tree_fluid->insert({"Fluid", "Pressure Poisson"}, 40.);
  tree_fluid->insert({"Fluid", "Postprocessing"}, 10.);
  tree_fluid->insert({"Fluid", "ALE update"}, 15.);


  std::shared_ptr<TimerTree> tree_structure;
  tree_structure.reset(new TimerTree());
  tree_structure->insert({"Elasticity", "Right-hand side"}, 2.);
  tree_structure->insert({"Elasticity", "Assemble"}, 9.);
  tree_structure->insert({"Elasticity", "Solve"}, 14.);
  // ... but can also be inserted last
  tree_structure->insert({"Elasticity"}, 25.);

  // Inserting sub-trees

  // inserting a sub-tree multiple times does not change the results
  tree.insert({"FSI"}, tree_fluid);
  tree.insert({"FSI"}, tree_fluid);

  // a new name can be provided for a sub-tree
  tree.insert({"FSI"}, tree_structure, "Structure");

  pcout << std::endl << "timings for level = 0:" << std::endl;
  tree.print_level(pcout, 0);
  pcout << std::endl << "timings for level = 1:" << std::endl;
  tree.print_level(pcout, 1);
  pcout << std::endl << "timings for level = 2:" << std::endl;
  tree.print_level(pcout, 2);
  pcout << std::endl << "timings all:" << std::endl;
  tree.print_plain(pcout);

  tree.clear();

  // should be empty after clear()
  pcout << std::endl << "timings for level = 0:" << std::endl;
  tree.print_level(pcout, 0);
  pcout << std::endl << "timings for level = 1:" << std::endl;
  tree.print_level(pcout, 1);
  pcout << std::endl << "timings for level = 2:" << std::endl;
  tree.print_level(pcout, 2);
  pcout << std::endl << "timings all:" << std::endl;
  tree.print_plain(pcout);

  // clear() must no touch sub-trees
  pcout << std::endl << "timings for level = 0:" << std::endl;
  tree_structure->print_level(pcout, 0);
  pcout << std::endl << "timings for level = 1:" << std::endl;
  tree_structure->print_level(pcout, 1);
  pcout << std::endl << "timings for level = 2:" << std::endl;
  tree_structure->print_level(pcout, 2);
  pcout << std::endl << "timings all:" << std::endl;
  tree_structure->print_plain(pcout);
}

int
main(int argc, char **argv)
{
  try
    {
      Utilities::MPI::MPI_InitFinalize mpi(argc, argv, 1);

      deallog.depth_console(0);

      test1();

      test2();
    }
  catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

  return 0;
}
