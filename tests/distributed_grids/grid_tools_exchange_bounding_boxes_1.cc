// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// test GridTools::exchange_local_bounding_boxes

#include <deal.II/base/bounding_box.h>
#include <deal.II/base/mpi.h>
#include <deal.II/base/point.h>

#include <deal.II/grid/grid_tools.h>

#include "../tests.h"

template <int spacedim>
void
test_exchange_bbox()
{
  // For process i the number of boxes n_bboxes[i%7] is created
  std::vector<unsigned int> n_bboxes = {2, 4, 3, 5, 1, 3, 8};

  const MPI_Comm mpi_communicator = MPI_COMM_WORLD;
  unsigned int   n_procs = Utilities::MPI::n_mpi_processes(mpi_communicator);
  unsigned int   proc    = Utilities::MPI::this_mpi_process(mpi_communicator);

  deallog << "Test for: dimension " << spacedim << std::endl;
  deallog << n_procs << "  mpi processes" << std::endl;

  // Creating the local bboxes
  // The first bbox first coordinate is proc: then use
  // the natural numbers increasingly
  std::vector<BoundingBox<spacedim>> loc_bboxes(n_bboxes[proc % 7]);
  for (unsigned int b = 0; b < loc_bboxes.size(); ++b)
    {
      Point<spacedim> pt1, pt2;
      for (unsigned int d = 0; d < spacedim; ++d)
        {
          pt1[d] = (proc + 2 * b * spacedim + d) * (proc + 1);
          pt2[d] = (proc + (2 * b + 1) * spacedim + d) * (proc + 1);
        }
      BoundingBox<spacedim> bbox(std::make_pair(pt1, pt2));
      loc_bboxes[b] = bbox;
    }

  std::vector<std::vector<BoundingBox<spacedim>>> global_boxes =
    GridTools::exchange_local_bounding_boxes(loc_bboxes, mpi_communicator);


  bool passed = true;
  // First check if the dimensions are correct:
  for (unsigned int i = 0; i < n_procs; ++i)
    {
      if (global_boxes[i].size() != n_bboxes[i % 7])
        {
          deallog << "Test FAILED: dimension check failed for process " << i
                  << std::endl;
          deallog << "Expected number of bboxes: " << n_bboxes[i % 7]
                  << std::endl;
          deallog << "Received number of bboxes: " << global_boxes[i].size()
                  << std::endl;
          passed = false;
        }
    }

  // Checking if all received data is correct
  if (passed)
    {
      for (unsigned int p = 0; p < n_procs; ++p)
        {
          for (unsigned int b = 0; b < n_bboxes[p % 7]; ++b)
            {
              std::vector<Point<spacedim>> boundary_pt(2);
              boundary_pt[0] = global_boxes[p][b].get_boundary_points().first;
              boundary_pt[1] = global_boxes[p][b].get_boundary_points().second;
              for (unsigned int d = 0; d < spacedim; ++d)
                for (unsigned int pt = 0; pt < 2; ++pt)
                  if (std::abs(boundary_pt[pt][d] -
                               (p + (2 * b + pt) * spacedim + d) * (p + 1)) >
                      1e-10)
                    {
                      passed = false;
                      deallog << "Test FAILED, value not corresponding for:"
                              << std::endl;
                      deallog << "Current proc: " << proc << std::endl;
                      deallog << "Value for proc: " << p << std::endl;
                      deallog << "BBox " << b << " Point " << pt
                              << " Dimension " << d << std::endl;
                      deallog << "Value found: " << (int)boundary_pt[pt][d]
                              << std::endl;
                      deallog << "Value expected: "
                              << (p + (2 * b + pt) * spacedim + d) * (p + 1)
                              << std::endl;
                    }
            }
        }
    }
  if (passed)
    deallog << "Test passed" << std::endl;
  else
    {
      deallog << "Test failed" << std::endl;
      deallog << "Current proc values" << std::endl;
      for (auto b : loc_bboxes)
        {
          deallog << b.get_boundary_points().first << " and " << std::endl;
          deallog << b.get_boundary_points().second << std::endl;
        }
      deallog << "Received values" << std::endl;
      for (auto vec_box : global_boxes)
        {
          for (auto b : vec_box)
            {
              deallog << b.get_boundary_points().first << " and " << std::endl;
              deallog << b.get_boundary_points().second << std::endl;
            }
        }
    }
}

int
main(int argc, char *argv[])
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  deallog << "Test: GridTools::exchange_local_bounding_boxes " << std::endl;

  test_exchange_bbox<1>();
  test_exchange_bbox<2>();
  test_exchange_bbox<3>();
}
